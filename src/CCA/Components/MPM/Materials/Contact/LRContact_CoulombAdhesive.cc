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

//#define JIMCOMP

using namespace Uintah;
using std::vector;
using std::string;

using namespace std;


LRContact_CoulombAdhesive::LRContact_CoulombAdhesive(const 	ProcessorGroup		* myworld
                                 	 	 	 	 	,	 	ProblemSpecP		& ps
												    ,		MaterialManagerP	& d_sS
													,		MPMLabel			* Mlb
													,		MPMFlags			* MFlag
  	  	  	  	  	  	  	  	  	  	  	  	  	): Contact(myworld, Mlb, MFlag, ps)
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
	// [CMAME] Nairn et. al., Comput. Methods Appl. Mech. Engrg. 362, 112859 (2020)

	Ghost::GhostType  gnone = Ghost::None;
	delt_vartype delT;
	old_dw->get(delT, lb->delTLabel, getLevel(patches));
	const double minimumMassToCalculate = 1.0e-16;

	int numMatls = d_materialManager->getNumMatls( "MPM" );
	ASSERTEQ(numMatls, matls->size());

	for(int p=0;p<patches->size();p++)  {
		const Patch* patch = patches->get(p);
		Vector dx = patch->dCell();

		double cellVolume = dx.x() * dx.y() * dx.z();
		double invCellVolume = 1.0/cellVolume;
		bool equalGrid = (dx.x() == dx.y() && dx.y() == dx.z());
		double hPerp = dx.x();

		constNCVariable<double> NC_CCweight;  	old_dw->get(NC_CCweight,      lb->NC_CCweightLabel,     0,patch, gnone, 0);
		constNCVariable<int> alphaMaterial;		new_dw->get(alphaMaterial,    lb->gAlphaMaterialLabel,  0,patch, gnone, 0);
		constNCVariable<Vector> normAlphaToBeta;	new_dw->get(normAlphaToBeta,  lb->gNormAlphaToBetaLabel,0,patch, gnone, 0);

		// Need access to all velocity fields at once
		std::vector<constNCVariable<double> >  gmass(numMatls);
		std::vector<constNCVariable<double> >  gvolume(numMatls);
		std::vector<constNCVariable<double> >  gmatlprominence(numMatls);
		std::vector<NCVariable<Vector> >       gvelocity(numMatls);

		for(int m=0;m<numMatls;m++){
			new_dw->get(gmass[m],          		lb->gMassLabel,     		materials->get(m), patch, gnone, 0);
			new_dw->get(gvolume[m],        		lb->gVolumeLabel,   		materials->get(m), patch, gnone, 0);
			new_dw->get(gmatlprominence[m],		lb->gMatlProminenceLabel,	materials->get(m), patch, gnone, 0);
			new_dw->getModifiable(gvelocity[m],   label_velocity_to_use,		materials->get(m), patch );
		}

		for(NodeIterator iter = patch->getNodeIterator(); !iter.done(); ++iter)  {
			IntVector nodeIndex = *iter;
			const double nodalWeight = 8.0 * NC_CCweight[nodeIndex];
			const int alpha = alphaMaterial[nodeIndex];
			// Calculate nodal volume for axisymmetric problems
			if(flag->d_axisymmetric)	{  // Nodal volume isn't constant for axisymmetry  // volume = r*dr*dtheta*dy  (dtheta = 1 radian)
				double r = min((patch->getNodePosition(nodeIndex)).x(),.5*dx.x());
				cellVolume =  r*dx.x()*dx.y();
				invCellVolume = 1.0/cellVolume;
			}
			if (alpha >= 0) { // Alpha != -99 means node has an interface and normal to interface
				Vector v_CoM(0.0, 0.0, 0.0);
				Vector p_CoM(0.0, 0.0, 0.0);
				double m_CoM = 0.0;
				double nodalVolume = 0.0;
				for (int matlIndex = 0; matlIndex < numMatls; ++matlIndex) {
//					if (d_matls.requested(matlIndex)) {
						p_CoM += gvelocity[matlIndex][nodeIndex] * gmass[matlIndex][nodeIndex];
						m_CoM += gmass[matlIndex][nodeIndex];
						nodalVolume += gvolume[matlIndex][nodeIndex] * nodalWeight;
//					}
				}
				double alphaMass = gmass[alpha][nodeIndex];

				bool debug_interface = false;
				v_CoM = p_CoM/m_CoM;
				if ( (nodalVolume * invCellVolume) > d_vol_const) { // Ensure node is populated enough to calculate on
					const Vector betaNormal = -1.0*normAlphaToBeta[nodeIndex];
					const double betaMass = m_CoM - alphaMass;
					for (int matlIndex = 0; matlIndex < numMatls; ++matlIndex) {
						if (d_matls.requested(matlIndex) && matlIndex != alpha) {
							double matlMass = gmass[matlIndex][nodeIndex];
							if (matlMass > minimumMassToCalculate) {
								if (debug_interface) {
								  std::cerr << "Patch: " << p << " Num Matls: " << numMatls;
								  std::cerr << " Node: " << nodeIndex << " alpha Material: " << alpha;
								  std::cerr << " Other material interacting: " << matlIndex;
								  std::cerr << " m_CoM: " << m_CoM << " Alpha Mass: " << alphaMass;
								  std::cerr << " Beta mass: " << betaMass << "  Matl mass: " << matlMass << std::endl;
								  if (betaMass <= 0) {
									  std::cerr << "      ";
									  for (int matlIndex = 0; matlIndex < numMatls; ++matlIndex) {
										  std::cerr << " mass[" << matlIndex << "]: " << gmass[matlIndex][nodeIndex];
									  }
									  std::cerr << std::endl;
								  }

								}
								double materialSeparation = gmatlprominence[matlIndex][nodeIndex] - gmatlprominence[alpha][nodeIndex];
								const Vector deltaVelocity = gvelocity[matlIndex][nodeIndex] - v_CoM;
								const double deltaV_normal = Dot(deltaVelocity, betaNormal);
								Vector normBeta_deltaV = deltaV_normal*betaNormal;
								Vector tangBeta_deltaV = deltaVelocity - normBeta_deltaV;
								double deltaV_tangent = tangBeta_deltaV.length();
								Vector betaTangent = tangBeta_deltaV/(deltaV_tangent+1.0e-100);
								double deltaV_tangent_reproject = Dot(deltaVelocity, betaTangent);
								betaTangent.safe_normalize(1.0e-100);
								if (!equalGrid) { // Untested -- TODO Test FIXME JBH 8/2020
									Vector a = betaTangent / dx;
									Vector b = Cross(betaNormal, betaTangent)/dx;
									a *= a;
									b *= b;
									hPerp = cellVolume*sqrt(Dot(a,a)*Dot(b,b));  // [CMES] eq. 22
								}
								double alphaVol = gvolume[alpha][nodeIndex];
								// According to [CMES] eq. 21:  V_cell = dx*dy*dz = hPerp * Ac
								double ac_vol = cellVolume/hPerp;
								// According to [CPM] eq. 11 contact area:
								double Ac = sqrt(2.0*nodalVolume*min(alphaVol,nodalVolume - alphaVol))/hPerp;
								double velocity_to_stress = betaMass / (Ac * delT);
								double stress_to_velocity = (Ac * delT)/betaMass;
								// Sstick Ac dt = betaMass * deltaV_tangent
								double S_stick = deltaV_tangent * velocity_to_stress;
								S_stick = deltaV_tangent_reproject * velocity_to_stress;
								// -N Ac dt = betaMass * deltaV_normal
								double N = - deltaV_normal * velocity_to_stress;
								// S_slide Ac dt = Sa Ac dt - mu (-N Ac dt) // [CPM] Eqs. 24, 15
								double S = d_shearAdhesion - d_mu * N; // [CPM] Eqn. 24
								double shear_ratio = S_stick*S_stick/(d_shearAdhesion * d_shearAdhesion);
								double normal_ratio = N * N / (d_normalAdhesion * d_normalAdhesion);

								// Set distances to check against
								const double adhesiveWidth = d_adhesiveThickness * hPerp; // Entire cell width = hPerp
								const double contactThreshold = 0.01*hPerp;

								Vector deltaV_corr(0.0, 0.0, 0.0);
								if (materialSeparation <= adhesiveWidth) {
									double v_shear_adhesion = d_shearAdhesion / velocity_to_stress;
#ifdef JIMCOMP
									Vector deltaVelocity_JIM(0.0, 0.0, 0.0);
									Vector normal_JIM(0.0, 0.0, 0.0);
									double normalDeltaVel_JIM = 0.0;
									Vector Dv_JIM(0.0, 0.0, 0.0);
									Vector normal_normaldV_JIM(0.0, 0.0, 0.0);
									Vector dV_normalDV_JIM(0.0, 0.0, 0.0);
									Vector surfaceTangent_JIM(0.0, 0.0, 0.0);
									double tangentDeltaVelocity_JIM = 0.0;
									double frictionCoefficient_JIM = 0.0;
									if (materialSeparation <= 0.01*dx.x()) {
										deltaVelocity_JIM = gvelocity[matlIndex][nodeIndex] - v_CoM;
										normal_JIM = -1.0 * normAlphaToBeta[nodeIndex];
										normalDeltaVel_JIM = Dot(deltaVelocity_JIM,normal_JIM);
										if (normalDeltaVel_JIM > 0.0) {
											normal_normaldV_JIM = normal_JIM*normalDeltaVel_JIM;
											dV_normalDV_JIM = deltaVelocity_JIM - normal_normaldV_JIM;
											surfaceTangent_JIM = dV_normalDV_JIM/(dV_normalDV_JIM.length()+1.e-100);
											tangentDeltaVelocity_JIM = Dot(deltaVelocity_JIM, surfaceTangent_JIM);
											frictionCoefficient_JIM = Min(d_mu, tangentDeltaVelocity_JIM/fabs(normalDeltaVel_JIM));

											Dv_JIM = -normal_normaldV_JIM - surfaceTangent_JIM * frictionCoefficient_JIM*fabs(normalDeltaVel_JIM);
										}
									}
#endif

									bool friction = false;
									bool adhesive = false;
									double v_slide = 0.0;
									double min_tang_magnitude = 0.0;
									double minSlideMagnitude = 0.0;
									Vector tangent_deltaV;
									Vector normal_corr_from_stress = N * stress_to_velocity * betaNormal;
									Vector tangent_corr_from_stress = Min(S, S_stick) * stress_to_velocity * betaTangent;
									if ((materialSeparation <= contactThreshold) && (deltaV_normal > 0.0)) { // Materials in direct compressive contact
				                        v_slide = (v_shear_adhesion + d_mu*deltaV_normal);
				                        min_tang_magnitude = Min(v_slide, deltaV_tangent_reproject);
										tangent_deltaV = min_tang_magnitude*betaTangent;
										minSlideMagnitude = Min(S,S_stick);
										deltaV_corr = normal_corr_from_stress - tangent_corr_from_stress;
#ifdef JIMCOMP
//										std::cerr << " FRICTION:" << std::endl;
										friction = true;
#endif
									} else if ((shear_ratio + normal_ratio) <= 1.0) { // Adhesion maintains contact
										deltaV_corr = -normBeta_deltaV - tangBeta_deltaV;
#ifdef JIMCOMP
//										std::cerr << "ADHESIVE:" << std::endl;
										adhesive = true;
#endif
									} // Correct based on separation of surfaces
#ifdef JIMCOMP
									if ((deltaV_corr - Dv_JIM).length() > 1.0e-16) {
										std::ios_base::fmtflags oldFlags = std::cerr.flags();
										int oldPrecision = std::cerr.precision();
										std::cerr.precision(8);
										std::cerr.flags(std::ios::scientific | std::ios::showpoint);
										std::cerr << scientific << setprecision(8);
										if (friction) {
											std::cerr << "FRICTION:" << std::endl;
											std::cerr << " Adhesion Strength: " << d_shearAdhesion << " hPerp: " << hPerp << " Ac: " << Ac << " dt: " << delT << " cell_vol: " << cellVolume << " alphaVol: " << alphaVol << " nodalVol: " << nodalVolume << " matlVol: " << gvolume[matlIndex][nodeIndex] << " matlMass: " << matlMass << " matlIndex: " << matlIndex << " alphaIndex: " << alpha << std::endl;
											std::cerr << "Slide velocity: " << v_slide << " |--| " << " Stick velocity: " << deltaV_tangent_reproject << std::endl;
										}
										if (adhesive) std::cerr << "ADHESION:" << std::endl;
										std::cerr << "\tComparison	    :            JIM             --            Mine        :  MagnitudeOfDifference\n";
										std::cerr << "\tNormal            : " << normal_JIM 				<< " -- " << betaNormal 				<< " : " << (normal_JIM-betaNormal).length() << std::endl;
										std::cerr << "\tdV                : " << deltaVelocity_JIM 			<< " -- " << deltaVelocity 				<< " : " << (deltaVelocity_JIM-deltaVelocity).length() << std::endl;
										std::cerr << "\tdV_along_normal   : " << normalDeltaVel_JIM 		<< " -- " << deltaV_normal 				<< " : " << fabs(normalDeltaVel_JIM - deltaV_normal) << std::endl;
										std::cerr << "\tnormal dV vector  : " << normal_normaldV_JIM 		<< " -- " << normBeta_deltaV 			<< " : " << (normal_normaldV_JIM-normBeta_deltaV).length() << std::endl;
										std::cerr << "\ttangent dV vector : " << dV_normalDV_JIM 			<< " -- " << tangBeta_deltaV 			<< " : " << (dV_normalDV_JIM-tangBeta_deltaV).length() << std::endl;
										std::cerr << "\ttangent vector    : " << surfaceTangent_JIM 		<< " -- " << betaTangent 				<< " : " << (surfaceTangent_JIM-betaTangent).length() << std::endl;
										std::cerr << "\tdV_along_tangent  : " << tangentDeltaVelocity_JIM 	<< " -- " << deltaV_tangent_reproject 	<< " : " << fabs(tangentDeltaVelocity_JIM - deltaV_tangent_reproject)<< std::endl;
										std::cerr << "\tdV_correction     : " << Dv_JIM 					<< " -- " << deltaV_corr 				<< " : " << (Dv_JIM - deltaV_corr).length() << std::endl;
										std::cerr << "\tnormal_correction : " << normal_normaldV_JIM        << " -- " << -normal_corr_from_stress   << " : " << (normal_normaldV_JIM + normal_corr_from_stress).length() << std::endl;
										std::cerr << "\ttangent_correction: " << surfaceTangent_JIM * frictionCoefficient_JIM * fabs(normalDeltaVel_JIM) << " -- " << tangent_corr_from_stress << " : " << ((surfaceTangent_JIM * frictionCoefficient_JIM * fabs(normalDeltaVel_JIM))-tangent_corr_from_stress).length() << std::endl;
										std::cerr << "\tslide vs. stick   : " << v_slide                    << " -- " << deltaV_tangent_reproject   << " : " << min_tang_magnitude << std::endl;
										std::cerr << "\tS/S_stick         : " << S                          << " -- " << S_stick                    << " : " << minSlideMagnitude << std::endl;
										std::cerr << "\tS/S_stick as vel  : " << S * stress_to_velocity     << " -- " << S_stick * stress_to_velocity  << " : " << minSlideMagnitude * stress_to_velocity << std::endl;
										std::cerr.flags(oldFlags);
										std::cerr.precision(oldPrecision);
									}
#endif
								} // Material is in surface interaction range
								// We change the velocity only, so we have to properly account for momentum balance by proper mass weighting.
								gvelocity[matlIndex][nodeIndex] += deltaV_corr;  					//  dp = (deltaV_corr) 						* matlMass
								gvelocity[alpha][nodeIndex] -= deltaV_corr * matlMass/alphaMass;	// -dp =-(deltaV_corr * matlMass/alphaMass) * alphaMass
							} // Matl mass high enough for calculation
						} // Material is in contact specification and is not the alpha material
					} // Loop over all materials
				} // Nodal volume occupied enough to justify calculation
			} // Node is multimaterial node.
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

void LRContact_CoulombAdhesive::addComputesAndRequiresInterpolated(		 SchedulerP 	& sched
                                                        		  ,const PatchSet		* patches
																  ,const MaterialSet	* ms
																  )
{
  Task * t = scinew Task("Friction::exMomInterpolated", this, &LRContact_CoulombAdhesive::exMomInterpolated);

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

void LRContact_CoulombAdhesive::addComputesAndRequiresIntegrated(		SchedulerP 	& sched
                                                       	   	    ,const 	PatchSet	* patches
																,const 	MaterialSet	* ms
																)
{
  Task * t = scinew Task("Friction::exMomIntegrated", this, &LRContact_CoulombAdhesive::exMomIntegrated);

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
