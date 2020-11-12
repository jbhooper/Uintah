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
  ps->getWithDefault("adhesive_layer_fraction",	d_adhesiveThickness, 	0.05);
  ps->getWithDefault("shear_adhesion", 			d_shearAdhesion, 		1e-100);
  ps->getWithDefault("normal_adhesion", 		d_normalAdhesion, 		1e-100);
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
  contact_ps->appendElement("type", "LR_coulomb_adhesive");
  contact_ps->appendElement("mu",                		d_mu);
  contact_ps->appendElement("adhesive_layer_fraction", 	d_adhesiveThickness);
  contact_ps->appendElement("shear_adhesion",			d_shearAdhesion);
  contact_ps->appendElement("normal_adhesion",			d_normalAdhesion);
  contact_ps->appendElement("volume_constraint", 		d_vol_const);
  contact_ps->appendElement("OneOrTwoStep",      		d_oneOrTwoStep);
  d_matls.outputProblemSpec(contact_ps);
}

void LRContact_CoulombAdhesive::enforceContact(const	PatchSubset		*	patches
											  ,const	MaterialSubset	*	materials
											  ,			DataWarehouse	*	old_dw
											  ,			DataWarehouse	*	new_dw
											  ,const    VarLabel        *   label_velocity_to_use
											  )
{
  // Referenced papers:
  // [CPM]   Nairn et. al., Computational Particle Mechanics, 5, 285-296 (2018)
  // [CMES]  Nairn,         Comput. Modeling in Engrg. & Sci., 92, 271-299 (2013)
  // [CMAME[ Nairn et. al., Comput. Methods Appl. Mech. Engrg. 362, 112859 (2020)

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

	  double hPerp = dx.x();
	  bool equalGrid = ((dx.x() == dx.y()) && (dx.y() == dx.z()) && (dx.x() == dx.z()));
	  double cellVolume = dx.x()*dx.y()*dx.z();
	  double invCellVolume = 1.0/(cellVolume);

	  delt_vartype delT;
	  old_dw->get(delT, lb->delTLabel, getLevel(patches));

	  constNCVariable<double> NC_CCweight;
	  constNCVariable<int> alphaMaterial;
	  constNCVariable<Vector> normAlphaToBeta;

	  old_dw->get(NC_CCweight,      lb->NC_CCweightLabel,     0,patch, gnone, 0);
	  new_dw->get(alphaMaterial,    lb->gAlphaMaterialLabel,  0,patch, gnone, 0);
	  new_dw->get(normAlphaToBeta,  lb->gNormAlphaToBetaLabel,0,patch, gnone, 0);

	  // Load per material arrays.
	  for(int m=0;m<numMatls;m++){
		  int dwi = materials->get(m);
		  new_dw->get(gmass[m],          lb->gMassLabel,     		dwi, patch, gnone, 0);
		  new_dw->get(gvolume[m],        lb->gVolumeLabel,   		dwi, patch, gnone, 0);
		  new_dw->get(gmatlprominence[m],lb->gMatlProminenceLabel,	dwi, patch, gnone, 0);

		  new_dw->getModifiable(gvelocity[m],   label_velocity_to_use,	dwi, patch );
	  } // Load material arrays.

	  for(NodeIterator iter = patch->getNodeIterator(); !iter.done(); ++iter)  {
		  IntVector nodeIndex = *iter;
		  const int alpha = alphaMaterial[nodeIndex];

		  // Calculate nodal volume for axisymmetric problems
		  if(flag->d_axisymmetric)	{  // Nodal volume isn't constant for axisymmetry
			  // volume = r*dr*dtheta*dy  (dtheta = 1 radian)
			  double r = min((patch->getNodePosition(nodeIndex)).x(),.5*dx.x());
			  cellVolume =  r*dx.x()*dx.y();
			  invCellVolume = 1.0/cellVolume;
		  }

		  Vector p_CoM(0.0, 0.0, 0.0);
		  double m_CoM = 0.0;
		  double nodalVolume = 0.0;
		  for (int matlIndex = 0; matlIndex < numMatls; ++matlIndex) {
			  if (d_matls.requested(matlIndex)) {
				  p_CoM += gvelocity[matlIndex][nodeIndex] * gmass[matlIndex][nodeIndex];
				  m_CoM += gmass[matlIndex][nodeIndex];
				  nodalVolume += gvolume[matlIndex][nodeIndex] * 8.0 * NC_CCweight[nodeIndex];
			  }
		  }

		  // Only work on nodes where alpha != -99 (i.e. multiple materials)
		  if( alpha >= 0 )  {
			  Vector v_CoM = p_CoM/m_CoM; // Calculate the center of mass velocity
//			  double nodalWeighting = 8.0 * NC_CCweight[nodeIndex];
//			  Vector p_CoM(0.0, 0.0, 0.0); // center of mass momentum
//			  double m_CoM=0.0;            // nodal mass for center of mass calcs
//			  double nodalVolume = 0.0;    // total nodal volume
//			  // Calculate nodal CoM quantities
//			  for (int matlIndex = 0; matlIndex < numMatls; ++matlIndex) {
//				  if (d_matls.requested(matlIndex)) { // Work only on requested materials
//					  p_CoM += gvelocity[matlIndex][nodeIndex] * gmass[matlIndex][nodeIndex];
//					  m_CoM += gmass[matlIndex][nodeIndex];
//					  nodalVolume += gvolume[matlIndex][nodeIndex] * nodalWeighting;
//				  } // Material is requested
//			  } // Loop over all materials
			  // Mommentum change necessary to enforce perfect contact:
			  // dP^a_i = m^a_i * (v_CoM - v^a_i) = d_n n^ + d_t t^ [CPM] eq. 12 & eq. 15
			  // d_n = -N * Ac * delT ; d_t = S_stick * Ac * delT   [CPM] eq. 15
			  double massAlpha = gmass[alpha][nodeIndex];
			  Vector deltaPAlpha = massAlpha*(v_CoM - gvelocity[alpha][nodeIndex]);
			  //  Contact cases:
			  //    Compression:  Materials in contact range and Dot(v^b-v^a,n^) < 0 [CMAME] eq. 14
			  //    Adhesive:     Materials in contact range, and (S_Stick/S_a)^2 + (N/N_a)^2 > 1 [CPM] eq. 25

			  // Only check contact if the node is full relative to a constraint
			  if ( (nodalVolume*invCellVolume) > d_vol_const) {
				  // Only directions for normal are +/- the alpha normal per LR description [CMAME]
				  Vector alphaNormal = normAlphaToBeta[nodeIndex];
				  double d_n = Dot(deltaPAlpha, alphaNormal); // [CPM] eq. 15
				  const Vector dP_alpha_normal = d_n * alphaNormal;
				  Vector alphaTangent = (deltaPAlpha - dP_alpha_normal).safe_normal();
				  double d_t = Dot(deltaPAlpha, alphaTangent); // [CPM] eq. 15
				  const Vector dP_alpha_tangent = d_t * alphaTangent;

				  // Calculate cell "width" for LR plane angle if grid not equal [CMES] eq. 20 [CPM] eq. 11
				  if (!equalGrid) { // JBH FIXME TODO This is untested! -- 8/2020
					  const Vector projDummy(Vector(1.0, 1.0, 1.0).safe_normalize(1e-100));
					  Vector alphaTangent = (projDummy - Dot(projDummy,alphaNormal)*alphaNormal);
					  alphaTangent.safe_normalize(1e-100);
					  Vector b = Cross(alphaNormal, alphaTangent)/dx; // b should be normal by definition
					  Vector a = alphaTangent / dx;
					  a *= a;
					  b *= b;
					  hPerp = cellVolume * sqrt(Dot(a,a)*Dot(b,b));
				  } // Calculate thickness on non-equal grids

				  // Calculate effective thickness of the contacting volume
				  //   perpendicular to the contacting surface. (Nairn, CMES 92, 271, 2013)

				  double alphaVol = gvolume[alpha][nodeIndex];  // Volume of alpha material
				  double betaVol = nodalVolume - alphaVol; 		// Volume of non-alpha materials
				  double Ac = sqrt(2.0*nodalVolume*min(alphaVol,betaVol))/hPerp; // Contact area

				  double momentum_to_stress = 1.0/(Ac * delT);
				  double normal_stress = -d_n * momentum_to_stress; // N
				  double tangent_stress = d_t * momentum_to_stress; // S_stick

				  // It's easiest to use the alpha correction as a reference
				  double alphaMass = gmass[alpha][nodeIndex];
//				  Vector v_alpha_normal = dP_alpha_normal / alphaMass;
//				  Vector v_alpha_tangent = dP_alpha_tangent / alphaMass;

				  // separationOffset: Maximum separation for contact interactions.
				  double separationOffset = d_adhesiveThickness*hPerp;
				  double zeroOffset = 0.01*hPerp;
				  if (hPerp != dx.x()) {
					  std::cerr << "ERROR:  hPerp != dx.x(): hPerp: " << hPerp << " dx.x(): " << dx.x() << " Grid is equal: " << equalGrid << " \n";
				  }

				  double betaMass = m_CoM-alphaMass; // Lumped mass of all non-alpha materials
				  for (int matlIndex = 0; matlIndex < numMatls; ++matlIndex) {
					  double matlMass = gmass[matlIndex][nodeIndex]; // Not the same as betaMass
					  double betaFraction = matlMass/betaMass;  // Fraction of beta material this material comprises
					  double massFraction = matlMass/m_CoM;     // Fraction of total mass made up of this specific material
					  if (d_matls.requested(matlIndex) && (matlIndex != alpha) && (betaFraction > 1e-5) && (massFraction > 1e-5)) {
						  // Material in contact model, not alpha, and mass fractions substantial enough to consider
						  double matlSeparation = gmatlprominence[matlIndex][nodeIndex] - gmatlprominence[alpha][nodeIndex]; //[CMAME] eq. 24
						  //  	   matlSeparation <= 0 				 => In contact. [CMAME] eq. 14; 24
						  //   0 < matlSeparation < separationOffset => Possible adhesion
						  if (matlSeparation <= separationOffset) { // In contact
							  Vector delta_v_material = gvelocity[matlIndex][nodeIndex] - v_CoM;
							  Vector deltaVelocity = gvelocity[matlIndex][nodeIndex] - v_CoM; //!! Jim TODO
							  Vector normal = -1.0*normAlphaToBeta[nodeIndex]; //!! Jim TODO
							  double normalDeltaVel = Dot(deltaVelocity,normal); //!! Jim TODO
							  // All non-alpha materials have a normal opposite the alpha normal
							  double matl_normal_dot_delta_v = Dot(-alphaNormal,delta_v_material);
							  if (fabs(matl_normal_dot_delta_v - normalDeltaVel) > 1e-16) {
							    std::cerr << "normalDeltaVel: " << normalDeltaVel << " matl_normal_dot_delta_v: " << matl_normal_dot_delta_v << " Difference: " << matl_normal_dot_delta_v - normalDeltaVel << std::endl;
							  }
							  Vector dv_beta(0.0, 0.0, 0.0);
							  if ((matlSeparation <= zeroOffset) && (matl_normal_dot_delta_v > 0.0)) { // Potentially compressive
								  Vector normal_normalDV = normal*normalDeltaVel; //!! Jim TODO
								  if ((normal_normalDV - dP_alpha_normal/betaMass).length() > 1.0e-16) {
								    std::cerr << "normal_normalDV: " << normal_normalDV << " dP_alpha_normal: " << dP_alpha_normal/betaMass << std::endl;
								  }
								  Vector dV_normalDV = deltaVelocity - normal_normalDV; //!! Jim TODO
								  Vector surfaceTangent = dV_normalDV/(dV_normalDV.length()+1.e-100);//!! Jim TODO
								  if ((surfaceTangent - alphaTangent).length() > 1.0e-16) {
									  std::cerr << "SurfaceTangent: " << surfaceTangent << " alphaTangent: " << alphaTangent << " differential length: " << (surfaceTangent-alphaTangent).length() << std::endl;
								  }
								  std::cerr << "Alpha_delta_V:" << normal_normalDV + dV_normalDV << " delta_v_from_dP: " << deltaPAlpha/betaMass << std::endl;
								  double tangentDeltaVelocity = Dot(deltaVelocity, surfaceTangent);//!! Jim TODO
								  std::cerr << "TangentDeltaVelocity: " << tangentDeltaVelocity << " tangentStress/betaMass: " << tangent_stress/betaMass << std::endl;
								  double frictionCoefficient = Min(d_mu,tangentDeltaVelocity/fabs(normalDeltaVel));//!! Jim TODO
								  // If materials in contact, then friction can apply:
								  // S_slide is the stress threshold necessary to transition from stick to slide.
								  double S_slide = d_shearAdhesion - d_mu * normal_stress; // [CPM Eq. 24]
								  if (fabs( fabs(S_slide*Ac*delT/betaMass) - fabs( d_mu*fabs(normalDeltaVel) ) ) > 1.0e-16) {
								    std::cerr << " frictionCoeff*fabs(normalDeltaVel): " << frictionCoefficient*fabs(normalDeltaVel) << " S_slide * Ac * delT: " << S_slide * Ac * delT/betaMass << " mu*N: " << d_mu*fabs(normalDeltaVel) << std::endl;
								  }
								  if (tangent_stress > S_slide) { // Required stress > that required to slide
									  // Tangential momentum of sliding:  S_slide * Ac * delta_t
									  dv_beta = (-dP_alpha_normal - (S_slide * Ac * delT)*alphaTangent)/betaMass;
									  std::cerr << "Slide:  Jim Dv: " << (-normal_normalDV - surfaceTangent*frictionCoefficient*fabs(normalDeltaVel)) << " My Dv: " << dv_beta << std::endl;
								  } else {
									  dv_beta = (dP_alpha_normal - dP_alpha_tangent)/betaMass; // Stick condition
									  std::cerr << "Stick:  Jim Dv: " << (-normal_normalDV - surfaceTangent*frictionCoefficient*fabs(normalDeltaVel)) << " My Dv: " << dv_beta << std::endl;
								  }
							  } else if (fabs(d_normalAdhesion) > 0 && fabs(d_shearAdhesion)> 0) { // Check for adhesion
								  double normal_stress_ratio = normal_stress / d_normalAdhesion;  // S_stick/S_a
								  double tangent_stress_ratio = tangent_stress / d_shearAdhesion; // N/N_a
								  bool is_adhered = (normal_stress_ratio*normal_stress_ratio + tangent_stress_ratio*tangent_stress_ratio) <= 1;
								  if (is_adhered) {
									  // std::cerr << "ADHESIVE! Normal Stress: " << normal_stress << " Tangent Stress: " << tangent_stress << "\n";
									  dv_beta = (dP_alpha_normal - dP_alpha_tangent)/betaMass; // Adhered, so stick conditions
								  } else {
									  // Do nothing; we have no constraints on the motion of the material.
									  // FIXME TODO Should we instead subtract the necessary velocity to overcome the adhesion?
								  }
							  }
//							   I have no idea what thi sis supposed to do.  It seems to always be 1? FIXME TODO JBH
//							   double fudgeFactor = max(1.0,(separationOffset - d_LR)/separationOffset); // Why?
//							   p_correct = p_correct * fudgeFactor;
//							   gvelocity[matlIndex][nodeIndex] += p_correct/matlMass;
//							   gvelocity[alpha][nodeIndex]     -= p_correct/alphaMass;

							  // Update this material's velocity and the alpha velocity.
							  // We update alpha velocity incrementally in case some portions of the ostensible beta mass
							  //   don't end up contributing.
							  gvelocity[alpha][nodeIndex] -= (dv_beta)*(gmass[matlIndex][nodeIndex]/alphaMass);
							  gvelocity[matlIndex][nodeIndex] += dv_beta;
						  } // Surface within contact/adhesion distance
					  } // Material is requested, not alpha, and has fractional sufficient fractional representation as part of overall and
					  //   response mass fraction.
				  }  // Loop over materials
			  }  // Volume fraction above d_vol_const
		  }  // multimaterial node (alpha > 0)
	  }  // Loop over nodes
  } // Loop over patches
}

void LRContact_CoulombAdhesive::exMomInterpolated(const ProcessorGroup	*			,
                                        		  const PatchSubset		* patches	,
												  const MaterialSubset	* matls		,
												  	    DataWarehouse	* old_dw	,
														DataWarehouse	* new_dw	)
{
  if(d_oneOrTwoStep==2) {
	  this->enforceContact(patches, matls, old_dw, new_dw, lb->gVelocityLabel);
  } // d_oneOrTwoSteps == 2
} // LRContact_CoulombAdhesive::exMomInterpolated

void LRContact_CoulombAdhesive::exMomIntegrated(const ProcessorGroup	*			,
                                        		const PatchSubset		* patches	,
												const MaterialSubset	* matls		,
										  	    	  DataWarehouse		* old_dw	,
													  DataWarehouse		* new_dw	)
{
  if(d_oneOrTwoStep==2) {
	  this->enforceContact(patches, matls, old_dw, new_dw, lb->gVelocityStarLabel);
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
