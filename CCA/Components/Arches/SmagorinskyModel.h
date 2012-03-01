/*

The MIT License

Copyright (c) 1997-2011 Center for the Simulation of Accidental Fires and 
Explosions (CSAFE), and  Scientific Computing and Imaging Institute (SCI), 
University of Utah.

License for the specific language governing rights and limitations under
Permission is hereby granted, free of charge, to any person obtaining a 
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation 
the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.

*/


//----- SmagorinskyModel.h --------------------------------------------------

#ifndef Uintah_Component_Arches_SmagorinskyModel_h
#define Uintah_Component_Arches_SmagorinskyModel_h

/**************************************
CLASS
   SmagorinskyModel
   
   Class SmagorinskyModel is an LES model for
   computing sub-grid scale turbulent viscosity.


GENERAL INFORMATION
   SmagorinskyModel.h - declaration of the class
   
   Author: Rajesh Rawat (rawat@crsim.utah.edu)

   All major modifications since 01.01.2004 done by:
   Stanislav Borodai(borodai@crsim.utah.edu)
      
   Creation Date:   Mar 1, 2000
   
   C-SAFE 
   
   Copyright U of U 2000

KEYWORDS


DESCRIPTION
   Class SmagorinskyModel is an LES model for
   computing sub-grid scale turbulent viscosity.


WARNING
   none
****************************************/

#include <CCA/Components/Arches/Arches.h>
#include <CCA/Components/Arches/TurbulenceModel.h>
#include <Core/Grid/Variables/SFCXVariable.h>
#include <Core/Grid/Variables/SFCYVariable.h>
#include <Core/Grid/Variables/SFCZVariable.h>

namespace Uintah {
class PhysicalConstants;
class BoundaryCondition;


class SmagorinskyModel: public TurbulenceModel {

public:
  
  // GROUP: Constructors:
  ////////////////////////////////////////////////////////////////////////
  // Blank constructor for SmagorinskyModel.
  SmagorinskyModel(const ArchesLabel* label, 
                   const MPMArchesLabel* MAlb,
                   PhysicalConstants* phyConsts,
                   BoundaryCondition* bndryCondition);

  // GROUP: Destructors:
  ////////////////////////////////////////////////////////////////////////
  // Virtual destructor for SmagorinskyModel.
  virtual ~SmagorinskyModel();

  // GROUP: Problem Setup :
  ///////////////////////////////////////////////////////////////////////
  // Set up the problem specification database
  virtual void problemSetup(const ProblemSpecP& db);

  // GROUP: Schedule Action :
  ///////////////////////////////////////////////////////////////////////
  // Schedule the recomputation of Turbulence Model data
  //    [in] 
  //        data User data needed for solve 
  virtual void sched_reComputeTurbSubmodel(SchedulerP&,
                                           const PatchSet* patches,
                                           const MaterialSet* matls,
                                         const TimeIntegratorLabel* timelabels);



  // GROUP: Schedule Action :
  ///////////////////////////////////////////////////////////////////////
  // Schedule the computation of Turbulence Model data
  //    [in] 
  //        data User data needed for solve 
  virtual void sched_computeScalarVariance(SchedulerP&,
                                           const PatchSet* patches,
                                           const MaterialSet* matls,
                                           const TimeIntegratorLabel* timelabels);
                                           
  virtual void sched_computeScalarDissipation(SchedulerP&,
                                              const PatchSet* patches,
                                              const MaterialSet* matls,
                                              const TimeIntegratorLabel* timelabels);
  // GROUP: Access Methods :
  ///////////////////////////////////////////////////////////////////////
  // Get the molecular viscosity
  double getMolecularViscosity() const; 

  ////////////////////////////////////////////////////////////////////////
  // Get the Smagorinsky model constant
  double getSmagorinskyConst() const {
    return d_CF;
  }
  inline void set3dPeriodic(bool periodic) {}
  inline double getTurbulentPrandtlNumber() const {
    return d_turbPrNo;
  }
  inline void setTurbulentPrandtlNumber(double turbPrNo) {
    d_turbPrNo = turbPrNo;
  }
  inline bool getDynScalarModel() const {
    return false;
  }

protected:
      PhysicalConstants* d_physicalConsts;
      BoundaryCondition* d_boundaryCondition;

private:

  // GROUP: Constructors (not instantiated):
  ////////////////////////////////////////////////////////////////////////
  // Blank constructor for SmagorinskyModel.
  SmagorinskyModel();

  // GROUP: Action Methods (private)  :
  ///////////////////////////////////////////////////////////////////////
  // Actually reCalculate the Turbulence sub model
  //    [in] 
  //        documentation here
  void reComputeTurbSubmodel(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset* matls,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw,
                             const TimeIntegratorLabel* timelabels);

  ///////////////////////////////////////////////////////////////////////
  // Actually Calculate the subgrid scale variance
  //    [in] 
  //        documentation here
  void computeScalarVariance(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset* matls,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw,
                             const TimeIntegratorLabel* timelabels);
                             
  void computeScalarDissipation(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw,
                                const TimeIntegratorLabel* timelabels);

 protected:
      double d_CF; //model constant
      double d_factorMesh; // lengthscale = fac_mesh*meshsize
      double d_filterl; // prescribed filter length scale
      double d_CFVar; // model constant for mixture fraction variance
      double d_turbPrNo; // turbulent prandtl number

 private:

      // const VarLabel* variables 

      inline double compute_smag_viscos ( constSFCXVariable<double>& u, constSFCYVariable<double> v, constSFCZVariable<double> w, 
          constCCVariable<Vector>& UCC,
          constCCVariable<double>& den, double pmixl, Vector Dx, IntVector c ){ 

        IntVector cxp = c + IntVector(1,0,0);
        IntVector cxm = c - IntVector(1,0,0);
        IntVector cyp = c + IntVector(0,1,0);
        IntVector cym = c - IntVector(0,1,0);
        IntVector czp = c + IntVector(0,0,1);
        IntVector czm = c - IntVector(0,0,1);

        double uep = u[cxp];
        double uwp = u[c];
        double unp = 0.50 * UCC[cyp].x();
        double usp = 0.50 * UCC[cym].x();
        double utp = 0.50 * UCC[czp].x();
        double ubp = 0.50 * UCC[czm].x();

        double vep = 0.50 * UCC[cxp].y();
        double vwp = 0.50 * UCC[cxm].y();
        double vnp = v[cyp];
        double vsp = v[c];
        double vtp = 0.50 * UCC[czp].y();
        double vbp = 0.50 * UCC[czm].y();

        double wep = 0.50 * UCC[cxp].z();
        double wwp = 0.50 * UCC[cxm].z();
        double wnp = 0.50 * UCC[cyp].z();
        double wsp = 0.50 * UCC[cym].z();
        double wtp = w[czp];
        double wbp = w[c];

        double s11 = (uep-uwp)/Dx.x();
        double s22 = (vnp-vsp)/Dx.y();
        double s33 = (wtp-wbp)/Dx.z();
        double s12 = 0.50 * ((unp-usp)/Dx.y() + (vep-vwp)/Dx.x());
        double s13 = 0.50 * ((utp-ubp)/Dx.z() + (wep-wwp)/Dx.x());
        double s23 = 0.50 * ((vtp-vbp)/Dx.z() + (wnp-wsp)/Dx.y());

        double IsI = 2.0 * ( pow(s11,2.0) + pow(s22,2.0) + pow(s33,2.0)
            + 2.0 * ( pow(s12,2) + pow(s13,2) + pow(s23,2) ) ); 

        IsI = std::sqrt( IsI ); 

        double turb_viscos = pow(pmixl, 2.0) * den[c] * IsI; 

        return turb_viscos; 

      };

}; // End class SmagorinskyModel
} // End namespace Uintah
  
  

#endif

// $Log : $

