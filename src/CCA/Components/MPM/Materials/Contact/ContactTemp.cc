#include "LRContact_CoulombAdhesive_Ancient.cc"

void LRContact_CoulombAdhesive::exMomInterpolated(const ProcessorGroup*,
                                        const PatchSubset* patches,
                                        const MaterialSubset* matls,
                                        DataWarehouse* old_dw,
                                        DataWarehouse* new_dw)
{
  if(d_oneOrTwoStep==2) {

   int numMatls = d_materialManager->getNumMatls( "MPM" );
   ASSERTEQ(numMatls, matls->size());

   // Need access to all velocity fields at once
   std::vector<constNCVariable<double> >  gmass(numMatls);
   std::vector<constNCVariable<double> >  gvolume(numMatls);
   std::vector<constNCVariable<double> >  gmatlprominence(numMatls);
   std::vector<NCVariable<Vector> >       gvelocity(numMatls);

   Ghost::GhostType  gnone = Ghost::None;

   for(int p=0;p<patches->size();p++)  {
	 const Patch* patch = patches->get(p);
	 Vector dx = patch->dCell();

	 bool equalGrid = false;
	 if (dx.x() == dx.y() == dx.z()) equalGrid = true;
	 double hPerp = dx.x();
	 double cellVolume = dx.x()*dx.y()*dx.z();
	 double invCellVolume = 1.0/cellVolume;


     constNCVariable<double> NC_CCweight;
     constNCVariable<int> alphaMaterial;
     constNCVariable<Vector> normAlphaToBeta;

     old_dw->get(NC_CCweight,      lb->NC_CCweightLabel,     0,patch, gnone, 0);
     new_dw->get(alphaMaterial,    lb->gAlphaMaterialLabel,  0,patch, gnone, 0);
     new_dw->get(normAlphaToBeta,  lb->gNormAlphaToBetaLabel,0,patch, gnone, 0);

     delt_vartype delT;
     old_dw->get(delT, lb->delTLabel, getLevel(patches));

     // Load per material arrays.
     for(int m=0;m<numMatls;m++){
       int dwi = matls->get(m);
       new_dw->get(gmass[m],          lb->gMassLabel,     			dwi, patch, gnone, 0);
       new_dw->get(gvolume[m],        lb->gVolumeLabel,   			dwi, patch, gnone, 0);
       new_dw->get(gmatlprominence[m],lb->gMatlProminenceLabel,		dwi, patch, gnone, 0);
       new_dw->getModifiable(gvelocity[m],   lb->gVelocityLabel,	dwi, patch );
     } // Load material arrays.

     const Vector projDummy(1.0, 1.0, 1.0);  // Dummy vector to decompose normals and tangents.
     for(NodeIterator iter = patch->getNodeIterator(); !iter.done();iter++)  {
       IntVector nodeIndex = *iter;
       int alpha = alphaMaterial[nodeIndex];
       if( alpha >= 0 )  {  // Only work on nodes where alpha!=-99 (i.e. multipler materials)

    	 // Calculate nodal volume for axisymmetric problems
      	 if(flag->d_axisymmetric)	{  // Nodal volume isn't constant for axisymmetry
      	   // volume = r*dr*dtheta*dy  (dtheta = 1 radian)
      	   double r = min((patch->getNodePosition(nodeIndex)).x(),.5*dx.x());
      	   cellVolume =  r*dx.x()*dx.y();
     	 }

      	 double nodalWeighting = 8.0 * NC_CCweight[nodeIndex];
    	 // Calculate nodal CoM quantities
      	 Vector p_CoM(0.0, 0.0, 0.0); // center of mass momentum
      	 double m_CoM=0.0;            // nodal mass for center of mass calcs
      	 double nodalVolume = 0.0;    // total nodal volume
      	 for (int matlIndex = 0; matlIndex < numMatls; ++matlIndex) {
           if (d_matls.requested(matlIndex)) {
      		 p_CoM += gvelocity[matlIndex][nodeIndex] * gmass[matlIndex][nodeIndex];
      		 m_CoM += gmass[matlIndex][nodeIndex];
      		 nodalVolume += gvolume[matlIndex][nodeIndex] * nodalWeighting;
           } // Material in contact definition
      	 } // Iterate over materials

      	// Vector v_CoM = p_CoM/m_CoM;

      	 // Only apply contact if the node is full relative to a constraint
      	 if ( (nodalVolume*invCellVolume) > d_vol_const) {
      	   // Only directions for normal are +/- the alpha normal, so do the vector mechanics once.
      	   Vector alphaNormal = normAlphaToBeta[nodeIndex];
      	   double alphaProminence = gmatlprominence[alpha][nodeIndex];
      	   double alphaMass = gmass[alpha][nodeIndex];
      	   alphaNormal.safe_normalize();  // Ensure the alpha normal is unit

      	   // Determine the perpendicular distance metric only if grid is not uniform.
      	   if (!equalGrid) { // Need to account for non-uniform grid spacings.
          	 Vector alphaTangent = (projDummy - Dot(projDummy,alphaNormal)*alphaNormal);
          	 alphaTangent.safe_normalize();
      		 Vector a = alphaTangent / dx;
      		 a *= a;
      		 Vector b = Cross(alphaNormal, alphaTangent)/dx;
      		 b *= b;
      		 hPerp = cellVolume * sqrt(Dot(a,a)*Dot(b,b));
      	   } // Grid is equal in all directions

      	   double separationOffset = d_AdhesiveOffset*hPerp;


      	   // Material normal and alpha normal are in opposite directions
      	   // Tangents are in the same direction.
      	   Vector matlNormal = -alphaNormal;
      	   for (int matlIndex = 0; matlIndex < numMatls; ++matlIndex) {
      		 double matlMass = gmass[matlIndex][nodeIndex];
      		 double remainderMass = alphaMass - matlMass;
//      		 Vector v_matl = gvelocity[matlIndex][nodeIndex];
      		 Vector p_matl = gvelocity[matlIndex][nodeIndex]*matlMass;

      		 if (d_matls.requested(matlIndex) && (matlIndex != alpha) && (matlMass > 1.0e-16)) {
      		   // Calculate surface separation
      		   double d_LR = gmatlprominence[matlIndex][nodeIndex] - alphaProminence;
      		   bool contact = (d_LR < 0);
      		   bool nearContact = (!contact && d_LR < separationOffset);
	      	   // Find the contact area for this material.
			   double matlVolume = gvolume[matlIndex][nodeIndex] * nodalWeighting;
			   double remainderVolume = nodalVolume - matlVolume;
			   double minVolume = Min(matlVolume,remainderVolume);
			   double contactArea = sqrt(2.0*nodalVolume*minVolume)/hPerp;
	      	   Vector p_correct(0.0, 0.0, 0.0);
			   if (contact or nearContact) { // Close enough to care about contact calculations
	      	     Vector delta_p = p_CoM - p_matl;
	      	     double N_Ac_dt = -Dot(delta_p, matlNormal);
	      	     Vector matlTangent = delta_p + N_Ac_dt*matlNormal;
	      	     double S_stick_Ac_dt = matlTangent.safe_normalize(1.0e-20);
	      	     Vector normalCorrect = -N_Ac_dt * matlNormal;
	      	     if (contact) { // Surfaces already in contact or trying to penetrate one another
	      	       double S_adhere_Ac_dt = d_shearAdhesion * contactArea * delT;
	      	       if (N_Ac_dt <= 0.0) { // Surfaces in compression
	      	    	 double S_slide_Ac_dt = S_adhere_Ac_dt + d_mu * N_Ac_dt;
	      	    	 if (S_stick_Ac_dt < S_slide_Ac_dt) {
	      	    	   p_correct = -N_Ac_dt * matlNormal + S_stick_Ac_dt * matlTangent; // Should be equal to delta_p
	      	    	 } // stick < slide
	      	    	 else {
	      	    	   p_correct = -N_Ac_dt * matlNormal + S_slide_Ac_dt * matlTangent; // Use frictional momentum instead
	      	    	 } // stick >= slide
	      	       } else { // Surface in tension
	      	    	   double shear_term = S_stick_Ac_dt / S_adhere_Ac_dt;
	      	    	   double normal_term = -N_Ac_dt/(d_normalAdhesion*contactArea*delT);
	      	    	   bool adhesionBroken = ((shear_term*shear_term + normal_term*normal_term) > 1.0);
	      	    	   if (!adhesionBroken) p_correct = -N_Ac_dt * matlNormal + S_stick_Ac_dt * matlTangent; // If still adhering, then stick
	      	       } // tension
	      	     } // contact
	      	     else { // Near contact, only check adhesion
	      	       double S_adhere_Ac_dt = d_shearAdhesion * contactArea * delT;
	      		   double shear_term = S_stick_Ac_dt / S_adhere_Ac_dt;
	      		   double normal_term = -N_Ac_dt/(d_normalAdhesion*contactArea*delT);
	      		   bool adhesionBroken = ((shear_term*shear_term + normal_term*normal_term) > 1.0);
	      		   if (!adhesionBroken) p_correct = -N_Ac_dt * matlNormal + S_stick_Ac_dt * matlTangent; // If still adhering, then stick
	      	     } // else (near Contact)
			   } // Contact or near contact
			   double fudgeFactor = max(1.0,(separationOffset - d_LR)/separationOffset); // Why?
			   p_correct = p_correct * fudgeFactor;
			   gvelocity[matlIndex][nodeIndex] += p_correct/matlMass;
			   gvelocity[alpha][nodeIndex]     -= p_correct/alphaMass;
             }  // d_matls.requested && matl != alpha && mass > 1e-16
      	   }  // Loop over materials
         }  // Volume fraction above d_vol_const
       }  // multimaterial node (alpha > 0)
     }  // Loop over nodes
   } // Loop over patches
  } // d_oneOrTwoSteps == 2
} // LRContact_CoulombAdhesive::exMomInterpolated
