#include <CCA/Components/MPM/Materials/Contact/LRContact_CoulombAdhesive.cc>

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
     }

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

      	 Vector v_CoM = p_CoM/m_CoM;

      	 // Only apply contact if the node is full relative to a constraint
      	 if ( (nodalVolume*invCellVolume) > d_vol_const) {
      	   // Only directions for normal are +/- the alpha normal, so do the vector mechanics once.
      	   Vector alphaNormal = normAlphaToBeta[nodeIndex];
      	   double alphaProminence = gmatlprominence[alpha][nodeIndex];
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
      	   }

      	   double separationOffset = 0.01*hPerp;


      	   // Material normal and alpha normal are in opposite directions
      	   // Tangents are in the same direction.
      	   Vector matlNormal = -alphaNormal;
      	   for (int matlIndex = 0; matlIndex < numMatls; ++matlIndex) {
      		 double matlMass = gmass[matlIndex][nodeIndex];
      		 Vector v_matl = gvelocity[matlIndex][nodeIndex];

      		 if (d_matls.requested(matlIndex) && (matlIndex != alpha) && (matlMass > 1.0e-16)) {
      		   // Calculate surface separation
      		   double separation = gmatlprominence[matlIndex][nodeIndex] - alphaProminence;
      		   bool contact = (separation <= separationOffset);
	      	   // Find the contact area for this material.
			   double matlVolume = gvolume[matlIndex][nodeIndex] * nodalWeighting;
			   double remainderVolume = nodalVolume - matlVolume;
			   double minVolume = Min(matlVolume,remainderVolume);
			   double contactArea = sqrt(2.0*nodalVolume*minVolume)/hPerp;
      		   if (contact) {
      			   Vector delta_v = v_matl - v_CoM;
      			   double dvDotn = Dot(delta_v,matlNormal);
      			   Vector v_adjusted(0.0, 0.0, 0.0);
      			   Vector v_normal_correct = matlNormal * dvDotn;
      			   Vector matlTangent = delta_v - v_normal_correct;
      			   double dvDott = matlTangent.safe_normalize(1.0e-20);
//      			   double mu_prime;
//      			   double mag_dvDotn = fabs(dvDotn);
      			   double v_Stick_Ac_dt = dvDott;
      			   if (dvDotn > 0.0) { // material in compression
      				 // Push normal to enforce contact
//      				 mu_prime = Min(d_mu,dvDott/mag_dvDotn);
      				 // dvDotn = -N * contactArea * delta_time
      				 // d_ShearAdhesion units = pressure = F/Area
      				 // F/Area * Area * delta_time = m*acceleration * delta_time = m*v = momentum
      				 // mass * dvDott = Stick * Area * delta_time
      				 // Stick = area * delta_time / (mass * velocity)
      				 double inv_contactArea_times_time = 1.0/(contactArea * delT);
      				 double v_Slide_Ac_dt = (d_ShearAdhesion)/(matlMass) - d_mu * dvDotn*inv_contactArea_times_time;
      				 double mu_prime_Ac_dt = Min(v_Stick_Ac_dt,v_Slide_Ac_dt);
      			   } else { // dvDotn < 0.0 means material is in tension

      			   }
      		   }

      		   Vector delta_p = matlMass * (v_CoM - gvelocity[matlIndex][nodeIndex]);
      		   double dpDotn = Dot(delta_p,matlNormal);  // d_n = -N*Ac*dt
      		   Vector t_dpDott = delta_p - dpDotn * matlNormal;
      		   double dpDott = t_dpDott.length(); // d_t = S_stick * Ac * dt
      		   Vector matlTangent = t_dpDott/(dpDott + 1.0e-100);

      		   // Determine if interface is under compression and/or in contact
      		   bool compression = (dpDotn > 0.0);
      		   bool stuck = false;
      		   Vector dp_corrected(0.0, 0.0, 0.0);
	      	   // Find the contact area for this material.
			   double matlVolume = gvolume[matlIndex][nodeIndex] * nodalWeighting;
			   double remainderVolume = nodalVolume - matlVolume;
			   double minVolume = Min(matlVolume,remainderVolume);
			   double A_c = sqrt(2.0*nodalVolume*minVolume)/hPerp;
      		   if (contact) {
      			   if (compression) {
      				 double S_stick_Ac_dt = dpDott;
      				 double Slide_Ac_dt = d_ShearAdhesion * A_c * delT - d_mu * m * dpDotn;
      				 double N_Ac_dt = dpDotn;
      				 Vector v_stick =
      				 double S_slide = d_ShearAdhesion * A_c * delT - d_mu * dpDotn;
      				 double S_coeff = Min(dpDott,S_slide);

      			   } else { // In contact, tensile stress

      				 double shear = dpDott / (d_ShearAdhesion * A_c * delT);
      				 double tensile = dpDotn / (d_TensileAdhesion * A_c * delT);

      				 stuck = (tensile*tensile + shear*shear) < 1.0;
      			   }

      		   }
      		   // Calculate momentum correction


      		   // Find the various frictional and adhesive quantities
      		   double stress_denom = 1.0/(A_c * delT);
      		   double S_slide_ac_dt = - d_mu * dpDotn;

      		   // Materials in compressive contact if separation < saparationOffset && dvDotN > 0.0
      		   // delta_v = Dot(dV,n)n + Dot(dV,t)t
      		   // delta_v - Dot(delta_v,n)*n = Dot(delta_v,t)t
      		   // Dot(delta_v,t)t = delta_v - Dot(delta_v,n)*n = delta_v - dvDotN*n = dvDotT*t
      		   // dP_matl = m_matl*(-dV_matl)
      		   // S_stick * Ac * del_T = Dot(m_matl*(-dV_matl),tangent)
      		   // dP_matl = m_matl*v_CoM - m_matl*v_matl = m_matl * delta_v (Nairn, Bard, Smith Eq. 12
      		   // dP_matl = (dp_matl.normal)normal + (dp_matl.tangent)tangent = (-N*A_c*del_t)normal + (Sstick*A_c*del_t)tangent
      		   // S_stick * Ac * del_T = Dot(dP_matl,tangent)
      		   // -N * A_c * del_T = Dot(dP_matl,normal)

      		   // Materials not adhering if (S_stick/S_a)^2 + (N/N_a)^2 > 1
      		   // N_a = adhesive strength
      		   // Materials in -adhesive- contact if

      		 }
      	   }

    	 // Loop over materials.  Only proceed if velocity field mass
    	 // is nonzero (not numerical noise) and the difference from
   	     // the centerOfMassVelocity is nonzero (More than one velocity
    	 // field is contributing to grid vertex).
    	   for(int matlIndex = 0; matlIndex < numMatls; matlIndex++) {
    		 double mass = gmass[matlIndex][nodeIndex];
    	     if (d_matls.requested(matlIndex) &&  (matlIndex != alpha) && (mass > 1.e-16)) {
    		 // Find separation of this material to the alpha material
    		   double separation = gmatlprominence[matlIndex][nodeIndex] - gmatlprominence[alpha][nodeIndex];
    		   Vector deltaVelocity = gvelocity[matlIndex][nodeIndex] - velocity_CoM;
    		   Vector surfaceNormal = -1.0*normAlphaToBeta[nodeIndex];
    		   double dvDotN = Dot(deltaVelocity,surfaceNormal);
    		   Vector dvProjN = surfaceNormal*dvDotN;
    		   Vector dvProjT = deltaVelocity - dvProjN;
    		   double S_stick = dvProjT.length(); // Magnitude of frictional force to maintain stick conditions.
    		   Vector surfaceTangent = dvProjT/(S_stick+1.e-100);
    		   double hPerp = dx.x();

    		   if (!equalGrid) { // Calculate our distance metric based on unequal grid sizes
    			   Vector a = surfaceTangent / dx;
    			   a*= a;
    			   Vector b = Cross(surfaceNormal,surfaceTangent)/dx;
    			   b*= b;
    			   hPerp = cellVolume * sqrt(Dot(a,a)*Dot(b,b));
    		   }

    		   double overlapOffsetDistance = hPerp*d_overlapOffsetMult;
    		   // If separation is less than threshold, materials are overlapped
    		   // Contact conditions:  separation < overlapOffsetDistance and dvDotN > 0.0
    		   // Adhesive conditions:  even if separation > overlapOffsetDistance or dvDotN < 0.0 adhesion
    		   //   may cause stick.
    		   if (separation <= overlapOffsetDistance) {
    			 double minVol=Min(gvolume[matlIndex][c],totalNodalVol-gvolume[matlIndex][c]);
    			 double contactArea = sqrt(2.0*totalNodalVol*minVol)/hPerp;

    			       Vector dV_friction(0.0, 0.0, 0.0);
    			       if (dvDotN > 0.0) { // Surface is moving into contact
    			    	   double muEffective = min(d_mu, S_stick);
    			       }

    									   // Velocity difference dictates displacement vector:
    									   //   U = u_n + u_t = UNCORRECTED displacement vector
    									   //   u_n = delta_velocity * delta_t * normal_vector
            			Vector deltaV_friction(0.0, 0.0, 0.0);
            			if (diffDotN > 0.0) { // Surface is moving toward other surface; probably not right for adhesion.
            				Vector diffProjN = normal_direction * diffDotN;
            				Vector diffProjT = velocityDifference - diffProjN; // Remainder is tangential component of velocity
            				Vector tangent_direction = diffProjT/(diffProjT.length()+1.e-100); // Instead diffProjT.normal()
            				double diffDotT = Dot(velocityDifference,tangent_direction);
            				// Friction coefficient is the lesser of the friction coefficient or a virtual
            				//   friction coefficient needed to fully express the tangent velocity change needed.
            				double frictionCoefficient = Min(d_mu, diffDotT/fabs(diffDotN));

            				// Re-implement velocity difference but subject to friction.  Not right for adhesion yet.
            				deltaV_friction = -diffProjN + tangent_direction*frictionCoefficient*fabs(diffDotN);

            				// Reduce friction if not in direct contact; linear reduction.
            				double reductionFactor = max(1.0,(0.01*dx.x() - separation)/0.01*dx.x());
            				deltaV_friction *= reductionFactor;

            				Vector deltaV_alphaMaterial = -deltaV_friction * gmass[matlIndex][nodeIndex]/gmass[alpha][nodeIndex];
            				gvelocity[matlIndex][nodeIndex] 	+= deltaV_friction;
            				gvelocity[alpha][c] += deltaV_alphaMaterial;
            			}
            		}

            	}

           if(!d_matls.requested(matlIndex)) continue;
           if(matlIndex==alpha) continue;
            double mass=gmass[matlIndex][nodeIndex];
            if(mass>1.e-16){ // There is mass of material beta at this node
              // Check relative separation of the material prominence
              double separation = gmatlprominence[matlIndex][nodeIndex] -
                                  gmatlprominence[alpha][nodeIndex];
              // If that separation is negative, the matls have overlapped
//              if(separation <= 0.0){
              if(separation <= d_overlapOffsetMult*dx.x()){
               Vector deltaVelocity=gvelocity[matlIndex][nodeIndex] - centerOfMassVelocity;
               Vector normal = -1.0*normAlphaToBeta[nodeIndex];
               double normalDeltaVel=Dot(deltaVelocity,normal);
               Vector Dv(0.,0.,0.);
               if(normalDeltaVel > 0.0){
                 Vector normal_normaldV = normal*normalDeltaVel;
                 Vector dV_normalDV = deltaVelocity - normal_normaldV;
                 Vector surfaceTangent = dV_normalDV/(dV_normalDV.length()+1.e-100);
                 double tangentDeltaVelocity=Dot(deltaVelocity,surfaceTangent);
                 double frictionCoefficient=
                         Min(d_mu,tangentDeltaVelocity/fabs(normalDeltaVel));

                 // Calculate velocity change needed to enforce contact
                 Dv = -normal_normaldV
                   -surfaceTangent*frictionCoefficient*fabs(normalDeltaVel);

#if 0
                // Define contact algorithm imposed strain, find maximum
                Vector epsilon=(Dv/dx)*delT;
                double epsilon_max=
                  Max(fabs(epsilon.x()),fabs(epsilon.y()),fabs(epsilon.z()));
                if(!compare(epsilon_max,0.0)){
                   epsilon_max *= Max(1.0, mass/(centerOfMassMass-mass));

                   // Scale velocity change if contact algorithm
                   // imposed strain is too large.
                   double ff=Min(epsilon_max,.5)/epsilon_max;
                   Dv=Dv*ff;
                }
#endif
                double ff = max(1.0,(.01*dx.x() - separation)/.01*dx.x());
                Dv=Dv*ff;
                Vector DvAlpha = -Dv*gmass[matlIndex][c]/gmass[alpha][nodeIndex];
                gvelocity[matlIndex][c]    +=Dv;
                gvelocity[alpha][nodeIndex]+=DvAlpha;
              } // if (relative velocity) * normal < 0
             }  // if separation
            }   // if !compare && !compare
          }     // matls
        }       // if (volume constraint)
      }         // if(alpha > 0)
    }           // NodeIterator
  }             // patches
 }              // if d_oneOrTwoStep
}
