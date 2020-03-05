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

// LRContact_CoulombAdhesive.h

#ifndef __LRCONTACT_COULOMBADHESIVE_H__
#define __LRCONTACT_COULOMBADHESIVE_H__

#include <CCA/Components/MPM/Materials/Contact/Contact.h>
#include <CCA/Components/MPM/Materials/Contact/ContactMaterialSpec.h> 
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManagerP.h>


namespace Uintah {
/**************************************

CLASS
   LRContact_CoulombAdhesiveADH
   
   This version of contact is based on John Nairn and Chad
   Hammerquist's 2019 manuscript that describes the use of logistic
   regression to find a common normal between objects, and uses
   particle geometry to find the most prominent portion of a particle
   at each node, and applies contact when contacting materials'
   prominences overlap.

   This model ads adhesion and is adapted from code constructed by Jim
   Guilkey.

GENERAL INFORMATION

   LRContact_CoulombAdhesive.h

   Justin B Hooper
   Department of Materials Science & Engineering
   University of Utah

  

KEYWORDS
   Contact_Model_Friction_Adhesion

  
****************************************/

      class LRContact_CoulombAdhesive : public Contact {
      private:
         
         // Prevent copying of this class
         // copy constructor
         LRContact_CoulombAdhesive(const LRContact_CoulombAdhesive &con);
         LRContact_CoulombAdhesive& operator=(const LRContact_CoulombAdhesive &con);
         
         MaterialManagerP    d_materialManager;
         
         // Coefficient of friction
         double d_mu;
         // Adhesive yield strength
         double d_adhesiveYield;
         // Nodal volume fraction that must occur before contact is applied
         double d_vol_const;
         // Separation Factor to consider materials not in contact.  (Not currently used)
         // double d_sepFac;
         int NGP;
         int NGN;

      public:
         // Constructor
         LRContact_CoulombAdhesive(const ProcessorGroup* myworld,
                         ProblemSpecP& ps, MaterialManagerP& d_sS,MPMLabel* lb,
                         MPMFlags* MFlag);
         
         // Destructor
         virtual ~LRContact_CoulombAdhesive();

         virtual void outputProblemSpec(ProblemSpecP& ps);

         // Basic contact methods
         virtual void exMomInterpolated(const ProcessorGroup*,
                                        const PatchSubset* patches,
                                        const MaterialSubset* matls,
                                        DataWarehouse* old_dw,
                                        DataWarehouse* new_dw);
         
         virtual void exMomIntegrated(const ProcessorGroup*,
                                      const PatchSubset* patches,
                                      const MaterialSubset* matls,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw);
         
         virtual void addComputesAndRequiresInterpolated(SchedulerP & sched,
                                             const PatchSet* patches,
                                             const MaterialSet* matls);

         virtual void addComputesAndRequiresIntegrated(SchedulerP & sched,
                                             const PatchSet* patches,
                                             const MaterialSet* matls);
      };
} // End namespace Uintah

#endif /* __LRCONTACT_COULOMBADHESIVE_H__ */
