//
// $Id$
//

// NullContact.cc
//
// One of the derived Contact classes.  This particular
// class is used when no contact is desired.  This would
// be used for example when a single velocity field is
// present in the problem, so doing contact wouldn't make
// sense.

#include "NullContact.h"
#include <SCICore/Geometry/Vector.h>
#include <SCICore/Geometry/IntVector.h>
#include <Uintah/Grid/Array3Index.h>
#include <Uintah/Grid/Grid.h>
#include <Uintah/Grid/Level.h>
#include <Uintah/Grid/NCVariable.h>
#include <Uintah/Grid/Patch.h>
#include <Uintah/Grid/NodeIterator.h>
#include <Uintah/Grid/ReductionVariable.h>
#include <Uintah/Grid/SimulationState.h>
#include <Uintah/Grid/SimulationStateP.h>
#include <Uintah/Interface/DataWarehouse.h>
#include <Uintah/Grid/Task.h>
#include <Uintah/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <Uintah/Components/MPM/MPMLabel.h>
using namespace Uintah::MPM;

NullContact::NullContact(ProblemSpecP& ps, SimulationStateP& d_sS)
{
  // Constructor
 
  IntVector v_f;
  ps->require("vel_fields",v_f);
  std::cout << "vel_fields = " << v_f << endl;

  d_sharedState = d_sS;

}

NullContact::~NullContact()
{
  // Destructor

}

void NullContact::initializeContact(const Patch* /*patch*/,
                                    int /*vfindex*/,
                                    DataWarehouseP& /*new_dw*/)
{

}

void NullContact::exMomInterpolated(const ProcessorGroup*,
				    const Patch* patch,
				    DataWarehouseP& /*old_dw*/,
				    DataWarehouseP& new_dw)
{

  //  All this does is carry forward the array from gVelocityLabel
  //  to gMomExedVelocityLabel

  int numMatls = d_sharedState->getNumMatls();
  int NVFs = d_sharedState->getNumVelFields();

  //  const MPMLabel* lb = MPMLabel::getLabels();

  // Retrieve necessary data from DataWarehouse
  vector<NCVariable<double> > gmass(NVFs);
  vector<NCVariable<Vector> > gvelocity(NVFs);
  for(int m = 0; m < numMatls; m++){
    Material* matl = d_sharedState->getMaterial( m );
    MPMMaterial* mpm_matl = dynamic_cast<MPMMaterial*>(matl);
    if(mpm_matl){
      int vfindex = matl->getVFIndex();
      new_dw->get(gvelocity[vfindex], lb->gVelocityLabel, vfindex, patch,
                  Ghost::None, 0);

      new_dw->put(gvelocity[vfindex], lb->gMomExedVelocityLabel, vfindex, patch);
    }
  }

  

}

void NullContact::exMomIntegrated(const ProcessorGroup*,
				  const Patch* patch,
                                  DataWarehouseP& /*old_dw*/,
                                  DataWarehouseP& new_dw)
{

  //  All this does is carry forward the array from gVelocityStarLabel
  //  and gAccelerationLabel to gMomExedVelocityStarLabel and 
  //  gMomExedAccelerationLabel respectively

  int numMatls = d_sharedState->getNumMatls();
  int NVFs = d_sharedState->getNumVelFields();
  //  const MPMLabel* lb = MPMLabel::getLabels();

  vector<NCVariable<double> > gmass(NVFs);
  vector<NCVariable<Vector> > gvelocity_star(NVFs);
  vector<NCVariable<Vector> > gacceleration(NVFs);
  for(int m = 0; m < numMatls; m++){
    Material* matl = d_sharedState->getMaterial( m );
    MPMMaterial* mpm_matl = dynamic_cast<MPMMaterial*>(matl);
    if(mpm_matl){
      int vfindex = matl->getVFIndex();
      new_dw->get(gvelocity_star[vfindex], lb->gVelocityStarLabel, vfindex,
                  patch, Ghost::None, 0);
      new_dw->get(gacceleration[vfindex], lb->gAccelerationLabel, vfindex,
                  patch, Ghost::None, 0);

    new_dw->put(gvelocity_star[vfindex], lb->gMomExedVelocityStarLabel,
							 vfindex, patch);
    new_dw->put(gacceleration[vfindex], lb->gMomExedAccelerationLabel,
							 vfindex, patch);
    }
  }

  
}

void NullContact::addComputesAndRequiresInterpolated( Task* t,
                                             const MPMMaterial* matl,
                                             const Patch* patch,
                                             DataWarehouseP& /*old_dw*/,
                                             DataWarehouseP& new_dw) const
{
  //  const MPMLabel* lb = MPMLabel::getLabels();
  int idx = matl->getDWIndex();
  t->requires( new_dw, lb->gVelocityLabel, idx, patch, Ghost::None);

  t->computes( new_dw, lb->gMomExedVelocityLabel, idx, patch );


}

void NullContact::addComputesAndRequiresIntegrated( Task* t,
                                             const MPMMaterial* matl,
                                             const Patch* patch,
                                             DataWarehouseP& /*old_dw*/,
                                             DataWarehouseP& new_dw) const
{

  //  const MPMLabel* lb = MPMLabel::getLabels();
  int idx = matl->getDWIndex();
  t->requires(new_dw, lb->gVelocityStarLabel, idx, patch, Ghost::None);
  t->requires(new_dw, lb->gAccelerationLabel, idx, patch, Ghost::None);

  t->computes( new_dw, lb->gMomExedVelocityStarLabel, idx, patch);
  t->computes( new_dw, lb->gMomExedAccelerationLabel, idx, patch);


}


// $Log$
// Revision 1.16  2000/09/25 20:23:20  sparker
// Quiet g++ warnings
//
// Revision 1.15  2000/07/05 23:43:36  jas
// Changed the way MPMLabel is used.  No longer a Singleton class.  Added
// MPMLabel* lb to various classes to retain the original calling
// convention.  Still need to actually fill the d_particleState with
// the various VarLabels that are used.
//
// Revision 1.14  2000/06/17 07:06:38  sparker
// Changed ProcessorContext to ProcessorGroup
//
// Revision 1.13  2000/06/03 05:25:46  sparker
// Added a new for pSurfLabel (was uninitialized)
// Uncommented pleaseSaveIntegrated
// Minor cleanups of reduction variable use
// Removed a few warnings
//
// Revision 1.12  2000/05/30 20:19:09  sparker
// Changed new to scinew to help track down memory leaks
// Changed region to patch
//
// Revision 1.11  2000/05/26 21:37:35  jas
// Labels are now created and accessed using Singleton class MPMLabel.
//
// Revision 1.10  2000/05/25 23:05:09  guilkey
// Created addComputesAndRequiresInterpolated and addComputesAndRequiresIntegrated
// for each of the three derived Contact classes.  Also, got the NullContact
// class working.  It doesn't do anything besides carry forward the data
// into the "MomExed" variable labels.
//
// Revision 1.9  2000/05/11 20:10:17  dav
// adding MPI stuff.  The biggest change is that old_dws cannot be const and so a large number of declarations had to change.
//
// Revision 1.8  2000/05/10 20:02:48  sparker
// Added support for ghost cells on node variables and particle variables
//  (work for 1 patch but not debugged for multiple)
// Do not schedule fracture tasks if fracture not enabled
// Added fracture directory to MPM sub.mk
// Be more uniform about using IntVector
// Made patches have a single uniform index space - still needs work
//
// Revision 1.7  2000/05/08 18:42:46  guilkey
// Added an initializeContact function to all contact classes.  This is
// a null function for all but the FrictionContact.
//
// Revision 1.6  2000/05/02 06:07:14  sparker
// Implemented more of DataWarehouse and SerialMPM
//
// Revision 1.5  2000/04/27 21:28:58  jas
// Contact is now created using a factory.
//
// Revision 1.4  2000/04/26 06:48:20  sparker
// Streamlined namespaces
//
// Revision 1.3  2000/03/20 23:50:44  dav
// renames SingleVel to SingleVelContact
//
// Revision 1.2  2000/03/20 17:17:12  sparker
// Made it compile.  There are now several #idef WONT_COMPILE_YET statements.
//
// Revision 1.1  2000/03/16 01:05:13  guilkey
// Initial commit for Contact base class, as well as a NullContact
// class and SingleVel, a class which reclaims the single velocity
// field result from a multiple velocity field problem.
//
