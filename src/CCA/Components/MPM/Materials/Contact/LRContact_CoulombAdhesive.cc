/*
 * The MIT License
 *
 * Copyright (c) 1997-2019 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <Core/Math/Matrix3.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Components/MPM/Materials/MPMMaterial.h>
#include <CCA/Components/MPM/Materials/Contact/LRContact_CoulombAdhesive.h>
#include <vector>
#include <iostream>

using namespace Uintah;
using std::vector;
using std::string;

using namespace std;


LRContact_CoulombAdhesive::LRContact_CoulombAdhesive(const ProcessorGroup* myworld,
                                 ProblemSpecP& ps,MaterialManagerP& d_sS,
                                 MPMLabel* Mlb,MPMFlags* MFlag)
  : Contact(myworld, Mlb, MFlag, ps)
{
  ps->require("mu",d_mu);
  ps->getWithDefault("adhesive_layer_fraction",	d_adhesiveThickness, 	0.01);
  ps->getWithDefault("shear_adhesion", 			d_shearAdhesion, 		0.0);
  ps->getWithDefault("normal_adhesion", 		d_normalAdhesion, 		0.0);
  ps->getWithDefault("volume_constraint",		d_vol_const, 			0.0);
  ps->getWithDefault("scale_adhesion",          d_scaleAdhesion,        false);
  ps->getWithDefault("OneOrTwoStep",     		d_oneOrTwoStep, 		2);

  d_materialManager = d_sS;

  if(flag->d_8or27==8){
    NGP=1;
    NGN=1;
  } else{
    NGP=2;
    NGN=2;
  }
}

LRContact_CoulombAdhesive::~LRContact_CoulombAdhesive()
{
  // Destructor
}

void LRContact_CoulombAdhesive::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type", "LRCoulombAdhesive");
  contact_ps->appendElement("mu",                		d_mu);
  contact_ps->appendElement("adhesive_layer_fraction", 	d_adhesiveThickness);
  contact_ps->appendElement("shear_adhesion",			d_shearAdhesion);
  contact_ps->appendElement("normal_adhesion",			d_normalAdhesion);
  contact_ps->appendElement("volume_constraint", 		d_vol_const);
  contact_ps->appendElement("OneOrTwoStep",      		d_oneOrTwoStep);
  d_matls.outputProblemSpec(contact_ps);
}

void LRContact_CoulombAdhesive::exMomInterpolated(const ProcessorGroup	*			,
                                        		  const PatchSubset		* patches	,
												  const MaterialSubset	* matls		,
												  	    DataWarehouse	* old_dw	,
														DataWarehouse	* new_dw	)
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
	 if ((dx.x() == dx.y()) && (dx.y() == dx.z()) && (dx.x() == dx.z())) equalGrid = true;
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

     const Vector projDummy(Vector(1.0, 1.0, 1.0).safe_normalize(1e-50));  // Dummy vector to decompose normals and tangents.

     for(NodeIterator iter = patch->getNodeIterator(); !iter.done(); ++iter)  {
       IntVector nodeIndex = *iter;
       int alpha = alphaMaterial[nodeIndex];
       if( alpha >= 0 )  {  // Only work on nodes where alpha!=-99 (i.e. multiple materials)

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
      	   double alphaMass = gmass[alpha][nodeIndex];
      	   //alphaNormal.safe_normalize();  // Ensure the alpha normal is unit

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

      	   double separationOffset = d_adhesiveThickness*hPerp;


      	   // Material normal and alpha normal are in opposite directions
      	   // Tangents are in the same direction.
      	   Vector matlNormal = -alphaNormal; // matlNormal = -normAlphaToBeta[nodeIndex]
      	   for (int matlIndex = 0; matlIndex < numMatls; ++matlIndex) {
      		 double matlMass = gmass[matlIndex][nodeIndex];
      		 Vector p_matl = gvelocity[matlIndex][nodeIndex]*matlMass;

      		 if (d_matls.requested(matlIndex) && (matlIndex != alpha) && (matlMass > 1.0e-16)) {
      		   bool adhesive = (d_shearAdhesion > 0) && (d_normalAdhesion > 0);
      		   // Calculate surface separation
      		   double d_LR = gmatlprominence[matlIndex][nodeIndex] - alphaProminence;
      		   bool contact = (d_LR < 0);
//      		   bool contact = (d_LR < separationOffset);
      		   bool nearContact = (!contact && d_LR < separationOffset);
	      	   // Find the contact area for this material.
			   double matlVolume = gvolume[matlIndex][nodeIndex] * nodalWeighting;
			   double remainderVolume = nodalVolume - matlVolume;
			   double minVolume = Min(matlVolume,remainderVolume);
			   double contactArea = sqrt(2.0*nodalVolume*minVolume)/hPerp;
			   adhesive = (adhesive && contactArea > 1.0e-16);

	      	   Vector p_correct(0.0, 0.0, 0.0);
			   if (contact or nearContact) { // Close enough to care about contact calculations
	      	     Vector delta_p = matlMass * v_CoM - p_matl; // m * (v_CoM - v) = m * (-deltaVelocity)
	      	     double N_Ac_dt = -Dot(delta_p, matlNormal); // -(n1*(P1CoM-p1M)+n2*(P2CoM-p2M)+n3*(P3CoM-p3M))

	      	     Vector matlTangent = delta_p + N_Ac_dt*matlNormal;
	      	     double S_stick_Ac_dt = matlTangent.safe_normalize(1.0e-20);
	      	     Vector normalCorrect = -N_Ac_dt * matlNormal;
	      	     if (contact) { // Surfaces already in contact or trying to penetrate one another
	      	       double S_adhere_Ac_dt = d_shearAdhesion * contactArea * delT;
	      	       if (N_Ac_dt > 0.0) { // Surfaces in compression
	      	    	 double S_slide_Ac_dt = S_adhere_Ac_dt + d_mu * N_Ac_dt;
	      	    	 if (S_stick_Ac_dt < S_slide_Ac_dt) {
	      	    	   p_correct = normalCorrect + S_stick_Ac_dt * matlTangent; // Should be equal to delta_p
	      	    	 } // stick < slide
	      	    	 else {
	      	    	   p_correct = normalCorrect + S_slide_Ac_dt * matlTangent; // Use frictional momentum instead
	      	    	 } // stick >= slide
	      	       }
	      	       else if (adhesive){ // Surface in tension
	      	    	   double shear_term = S_stick_Ac_dt / S_adhere_Ac_dt;
	      	    	   double normal_term = -N_Ac_dt/(d_normalAdhesion*contactArea*delT);
	      	    	   bool adhesionBroken = ((shear_term*shear_term + normal_term*normal_term) > 1.0);
	      	    	   if (!adhesionBroken) p_correct = normalCorrect + S_stick_Ac_dt * matlTangent; // If still adhering, then stick
	      	       } // tension
	      	     } // contact
	      	     else if (adhesive) { // Near contact, only check adhesion
		      	   double scaleFactor = 1.0;
		      	   if (d_scaleAdhesion) scaleFactor = 1.01 - (d_LR/separationOffset); // Linearly scales adhesion as pulloff occurs
	      	       double S_adhere_Ac_dt = scaleFactor * d_shearAdhesion * contactArea * delT;
	      		   double shear_term = S_stick_Ac_dt / S_adhere_Ac_dt;
	      		   double normal_term = -N_Ac_dt/(scaleFactor * d_normalAdhesion*contactArea*delT);
	      		   bool adhesionBroken = ((shear_term*shear_term + normal_term*normal_term) > 1.0);
	      		   if (!adhesionBroken) p_correct = normalCorrect + S_stick_Ac_dt * matlTangent; // If still adhering, then stick
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

void LRContact_CoulombAdhesive::exMomIntegrated(const ProcessorGroup	*			,
                                        		const PatchSubset		* patches	,
												const MaterialSubset	* matls		,
										  	    	  DataWarehouse		* old_dw	,
													  DataWarehouse		* new_dw	)
{
  if(d_oneOrTwoStep==2) {

   int numMatls = d_materialManager->getNumMatls( "MPM" );
   ASSERTEQ(numMatls, matls->size());

   // Need access to all velocity fields at once
   std::vector<constNCVariable<double> >  gmass(numMatls);
   std::vector<constNCVariable<double> >  gvolume(numMatls);
   std::vector<constNCVariable<double> >  gmatlprominence(numMatls);
   std::vector<NCVariable<Vector> >       gvelocity_star(numMatls);

   Ghost::GhostType  gnone = Ghost::None;

   for(int p=0;p<patches->size();p++)  {
	 const Patch* patch = patches->get(p);
	 Vector dx = patch->dCell();

	 bool equalGrid = false;
	 if ((dx.x() == dx.y()) && (dx.y() == dx.z()) && (dx.x() == dx.z())) equalGrid = true;
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
       new_dw->get(gmass[m],          			lb->gMassLabel,     			dwi, patch, gnone, 0);
       new_dw->get(gvolume[m],        			lb->gVolumeLabel,   			dwi, patch, gnone, 0);
       new_dw->get(gmatlprominence[m],			lb->gMatlProminenceLabel,		dwi, patch, gnone, 0);
       new_dw->getModifiable(gvelocity_star[m],	lb->gVelocityStarLabel,			dwi, patch 			);
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
      		 p_CoM += gvelocity_star[matlIndex][nodeIndex] * gmass[matlIndex][nodeIndex];
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

      	   double separationOffset = d_adhesiveThickness*hPerp;


      	   // Material normal and alpha normal are in opposite directions
      	   // Tangents are in the same direction.
      	   Vector matlNormal = -alphaNormal;
      	   for (int matlIndex = 0; matlIndex < numMatls; ++matlIndex) {
      		 double matlMass = gmass[matlIndex][nodeIndex];
      		 Vector p_matl = gvelocity_star[matlIndex][nodeIndex]*matlMass;

      		 if (d_matls.requested(matlIndex) && (matlIndex != alpha) && (matlMass > 1.0e-16)) {
      		   bool adhesive = (d_shearAdhesion > 0) && (d_normalAdhesion > 0);
      		   // Calculate surface separation
      		   double d_LR = gmatlprominence[matlIndex][nodeIndex] - alphaProminence;
      		   bool contact = (d_LR < 0);
      		   bool nearContact = (!contact && d_LR < separationOffset);
	      	   // Find the contact area for this material.
			   double matlVolume = gvolume[matlIndex][nodeIndex] * nodalWeighting;
			   double remainderVolume = nodalVolume - matlVolume;
			   double minVolume = Min(matlVolume,remainderVolume);
			   double contactArea = sqrt(2.0*nodalVolume*minVolume)/hPerp;
			   adhesive = adhesive && (contactArea > 1.0e-16);
	      	   Vector p_correct(0.0, 0.0, 0.0);
			   if (contact or nearContact) { // Close enough to care about contact calculations
	      	     Vector delta_p = matlMass*v_CoM - p_matl;
	      	     double N_Ac_dt = -Dot(delta_p, matlNormal);
	      	     Vector matlTangent = delta_p + N_Ac_dt*matlNormal;
	      	     double S_stick_Ac_dt = matlTangent.safe_normalize(1.0e-20);
	      	     Vector normalCorrect = -N_Ac_dt * matlNormal;
	      	     if (contact) { // Surfaces already in contact or trying to penetrate one another
	      	       double S_adhere_Ac_dt = d_shearAdhesion * contactArea * delT;
	      	       if (N_Ac_dt <= 0.0) { // Surfaces in compression
	      	    	 double S_slide_Ac_dt = S_adhere_Ac_dt + d_mu * N_Ac_dt;
	      	    	 if (S_stick_Ac_dt < S_slide_Ac_dt) {
	      	    	   p_correct = normalCorrect + S_stick_Ac_dt * matlTangent; // Should be equal to delta_p
	      	    	 } // stick < slide
	      	    	 else {
	      	    	   p_correct = normalCorrect + S_slide_Ac_dt * matlTangent; // Use frictional momentum instead
	      	    	 } // stick >= slide
	      	       } else if (adhesive) { // Surface in tension
	      	    	   double shear_term = S_stick_Ac_dt / S_adhere_Ac_dt;
	      	    	   double normal_term = -N_Ac_dt/(d_normalAdhesion*contactArea*delT);
	      	    	   bool adhesionBroken = ((shear_term*shear_term + normal_term*normal_term) > 1.0);
	      	    	   if (!adhesionBroken) p_correct = normalCorrect + S_stick_Ac_dt * matlTangent; // If still adhering, then stick
	      	       } // tension
	      	     } // contact
	      	     else if (adhesive) { // Near contact, only check adhesion
	      	       double scaleFactor = 1.0;
	      	       if (d_scaleAdhesion) scaleFactor = 1.01 - (d_LR/separationOffset); // Linearly scales adhesion as pulloff occurs
	      	       double S_adhere_Ac_dt = scaleFactor * d_shearAdhesion * contactArea * delT;
	      		   double shear_term = S_stick_Ac_dt / S_adhere_Ac_dt;
	      		   double normal_term = -N_Ac_dt/(scaleFactor * d_normalAdhesion*contactArea*delT);
	      		   bool adhesionBroken = ((shear_term*shear_term + normal_term*normal_term) > 1.0);
	      		   if (!adhesionBroken) p_correct = normalCorrect + S_stick_Ac_dt * matlTangent; // If still adhering, then stick
	      	     } // else (near Contact)
			   } // Contact or near contact
			   double fudgeFactor = max(1.0,(separationOffset - d_LR)/separationOffset); // Why?
			   p_correct = p_correct * fudgeFactor;
			   gvelocity_star[matlIndex][nodeIndex] += p_correct/matlMass;
			   gvelocity_star[alpha][nodeIndex]     -= p_correct/alphaMass;
             }  // d_matls.requested && matl != alpha && mass > 1e-16
      	   }  // Loop over materials
         }  // Volume fraction above d_vol_const
       }  // multimaterial node (alpha > 0)
     }  // Loop over nodes
   } // Loop over patches
  } // d_oneOrTwoSteps == 2
} // LRContact_CoulombAdhesive::exMomIntegrated

void LRContact_CoulombAdhesive::addComputesAndRequiresInterpolated(SchedulerP & sched,
                                                        const PatchSet* patches,
                                                        const MaterialSet* ms)
{
  Task * t = scinew Task("Friction::exMomInterpolated", 
                      this, &LRContact_CoulombAdhesive::exMomInterpolated);

  MaterialSubset* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();
  
  const MaterialSubset* mss = ms->getUnion();
  t->requires(Task::OldDW, lb->delTLabel);
  t->requires(Task::NewDW, lb->gMassLabel,                  		Ghost::None);
  t->requires(Task::NewDW, lb->gVolumeLabel,                		Ghost::None);
  t->requires(Task::NewDW, lb->gMatlProminenceLabel,        		Ghost::None);
  t->requires(Task::NewDW, lb->gAlphaMaterialLabel,         		Ghost::None);
  t->requires(Task::NewDW, lb->gNormAlphaToBetaLabel,	z_matl,		Ghost::None);
  t->requires(Task::OldDW, lb->NC_CCweightLabel,		z_matl,     Ghost::None);
  t->modifies(			   lb->gVelocityLabel,      	mss					   );

  sched->addTask(t, patches, ms);

  if (z_matl->removeReference())
    delete z_matl; // shouln't happen, but...
}

void LRContact_CoulombAdhesive::addComputesAndRequiresIntegrated(SchedulerP & sched,
                                                       const PatchSet* patches,
                                                       const MaterialSet* ms) 
{
  Task * t = scinew Task("Friction::exMomIntegrated", 
                      this, &LRContact_CoulombAdhesive::exMomIntegrated);

  MaterialSubset* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();
  
  const MaterialSubset* mss = ms->getUnion();
  t->requires(Task::OldDW, lb->delTLabel);
  t->requires(Task::NewDW, lb->gMassLabel,                  		Ghost::None);
  t->requires(Task::NewDW, lb->gVolumeLabel,                		Ghost::None);
  t->requires(Task::NewDW, lb->gMatlProminenceLabel,        		Ghost::None);
  t->requires(Task::NewDW, lb->gAlphaMaterialLabel,         		Ghost::None);
  t->requires(Task::OldDW, lb->NC_CCweightLabel,		z_matl,     Ghost::None);
  t->requires(Task::NewDW, lb->gNormAlphaToBetaLabel,	z_matl,		Ghost::None);
  t->modifies(             lb->gVelocityStarLabel,  	mss					   );

  sched->addTask(t, patches, ms);

  if (z_matl->removeReference())
    delete z_matl; // shouln't happen, but...
}
