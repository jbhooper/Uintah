/*
 * The MIT License
 *
 * Copyright (c) 1997-2013 The University of Utah
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

#ifndef UINTAH_MD_ELECTROSTATICS_SPME_H
#define UINTAH_MD_ELECTROSTATICS_SPME_H

#include <CCA/Components/MD/Electrostatics.h>
#include <CCA/Components/MD/ShiftedCardinalBSpline.h>
#include <CCA/Components/MD/SimpleGrid.h>
#include <CCA/Components/Schedulers/OnDemandDataWarehouse.h>
#include <CCA/Components/MD/SPMEPatch.h>
#include <CCA/Components/MD/PatchMaterialKey.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Thread/ConditionVariable.h>

#include <vector>

#include <sci_defs/fftw_defs.h>

namespace Uintah {

  using namespace SCIRun;

  typedef std::complex<double> dblcomplex;
  typedef int particleIndex;
  typedef int particleId;

  class MDSystem;
  class SPMEMapPoint;
  class ParticleSubset;
  class MDLabel;

  /**
   *  @class SPME
   *  @ingroup MD
   *  @author Alan Humphrey and Justin Hooper
   *  @date   January, 2013
   *
   *  @brief
   *
   *  @param
   */
  class SPME : public Electrostatics {

    public:

      /**
       * @brief
       * @param
       */
      SPME();

      /**
       * @brief
       * @param
       * @param
       * @param
       * @param
       */
      SPME(MDSystem* system,
           const double ewaldBeta,
           const bool polarizable,
           const double polarizationTolerance,
           const IntVector& kLimits,
           const int splineOrder,
           const int maxPolarizableIterations = 0);

      /**
       * @brief
       * @param
       */
      ~SPME();

      /**
       * @brief
       * @param
       * @param
       * @param
       * @return
       */
      void initialize(const ProcessorGroup* pg,
                      const PatchSubset* patches,
                      const MaterialSubset* materials,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw);

      /**
       * @brief
       * @param None
       * @return None
       */
      void setup(const ProcessorGroup* pg,
                 const PatchSubset* patches,
                 const MaterialSubset* materials,
                 DataWarehouse* old_dw,
                 DataWarehouse* new_dw);

      /**
       * @brief
       * @param
       * @return None
       */
      void calculate(const ProcessorGroup* pg,
                     const PatchSubset* perProcPatches,
                     const MaterialSubset* materials,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw,
                     SchedulerP& subscheduler,
                     const LevelP& level,
                     SimulationStateP& sharedState);

      /**
       * @brief
       * @param None
       * @return None
       */
      void finalize(const ProcessorGroup* pg,
                    const PatchSubset* patches,
                    const MaterialSubset* materials,
                    DataWarehouse* old_dw,
                    DataWarehouse* new_dw);

      /**
       * @brief
       * @param None
       * @return
       */
      inline ElectrostaticsType getType() const
      {
        return d_electrostaticMethod;
      }

      /**
       * @brief
       * @param
       * @return
       */
      inline void setMDLabel(MDLabel* lb)
      {
        d_lb = lb;
      }

    private:

      /**
       * @brief
       * @param
       * @return
       */
      void scheduleCalculatePreTransform(SchedulerP& sched,
                                         const ProcessorGroup* pg,
                                         const PatchSet* patches,
                                         const MaterialSet* materials,
                                         DataWarehouse* subOldDW,
                                         DataWarehouse* subNewDW);

      /**
       * @brief
       * @param
       * @return
       */
      void scheduleReduceNodeLocalQ(SchedulerP& sched,
                                    const ProcessorGroup* pg,
                                    const PatchSet* patches,
                                    const MaterialSet* materials,
                                    DataWarehouse* subOldDW,
                                    DataWarehouse* subNewDW);

      /**
       * @brief
       * @param
       * @return
       */
      void scheduleTransformRealToFourier(SchedulerP& sched,
                                          const ProcessorGroup* pg,
                                          const PatchSet* perProcPatches,
                                          const MaterialSet* materials,
                                          DataWarehouse* subOldDW,
                                          DataWarehouse* subNewDW,
                                          const LevelP& level);

      /**
       * @brief
       * @param
       * @return
       */
      void scheduleCalculateInFourierSpace(SchedulerP& sched,
                                           const ProcessorGroup* pg,
                                           const PatchSet* patches,
                                           const MaterialSet* materials,
                                           DataWarehouse* subOldDW,
                                           DataWarehouse* subNewDW);

      /**
       * @brief
       * @param
       * @return
       */
      void scheduleTransformFourierToReal(SchedulerP& sched,
                                          const ProcessorGroup* pg,
                                          const PatchSet* perProcPatches,
                                          const MaterialSet* materials,
                                          DataWarehouse* subOldDW,
                                          DataWarehouse* subNewDW,
                                          const LevelP& level);

      /**
       * @brief
       * @param
       * @return
       */
      void scheduleDistributeNodeLocalQ(SchedulerP& sched,
                                        const ProcessorGroup* pg,
                                        const PatchSet* patches,
                                        const MaterialSet* materials,
                                        DataWarehouse* subOldDW,
                                        DataWarehouse* subNewDW);

      /**
       * @brief
       * @param
       * @param
       * @return
       */
      void generateChargeMap(std::vector<SPMEMapPoint>* chargeMap,
                             ParticleSubset* pset,
                             constParticleVariable<Point>& particlePositions,
                             constParticleVariable<long64>& particleIDs);

      /**
       * @brief Map points (charges) onto the underlying grid.
       * @param
       * @param
       * @param
       * @return
       */
      void mapChargeToGrid(SPMEPatch* spmePatch,
                           const std::vector<SPMEMapPoint>* gridMap,
                           ParticleSubset* pset,
                           constParticleVariable<double>& charges);

      /**
       * @brief Map forces from grid back to points.
       * @param
       * @return
       */
      void mapForceFromGrid(SPMEPatch* spmePatch,
                            const std::vector<SPMEMapPoint>* gridMap,
                            ParticleSubset* pset,
                            constParticleVariable<double>& charges,
                            ParticleVariable<Vector>& pforcenew);

      /**
       * @brief
       * @param
       * @param
       * @return
       */
      vector<double> calculateOrdinalSpline(const int orderMinusOne,
                                            const int splineOrder);

      /**
       * @brief Generates the local portion of the B grid (see. Essmann et. al., J. Phys. Chem. 103, p 8577, 1995)
       *          Equation 4.8
       * @param
       * @param Extents - The number of internal grid points on the current processor
       * @param Offsets - The global mapping of the local (0,0,0) coordinate into the global grid index
       *
       * @return None
       */
      void calculateBGrid(SimpleGrid<double>& BGrid,
                                        const IntVector& localExtents,
                                        const IntVector& globalOffset) const;

      /**
       * @brief
       * @param
       * @param
       * @param
       * @param
       * @return
       */
      void generateBVector(std::vector<dblcomplex>& bVector,
                           const std::vector<double>& mFractional,
                           const int initialIndex,
                           const int localGridExtent) const;

      /**
       * @brief
       * @param
       * @param
       * @param
       * @param
       * @return
       */
      void generateBVectorChunk(std::vector<dblcomplex>& bVector,
                                const int m_initial,
                                const int localGridExtent,
                                const int K) const;

      /**
       * @brief Generates the local portion of the C grid (see. Essmann et. al., J. Phys. Chem. 103, p 8577, 1995)
       *          Equation 3.9
       * @param Extents - The number of internal grid points on the current processor
       * @param Offsets - The global mapping of the local (0,0,0) coordinate into the global grid index
       *
       * @return None
       */
      void calculateCGrid(SimpleGrid<double>& CGrid,
                          const IntVector& extents,
                          const IntVector& offsets) const;

      /**
       * @brief
       * @param stressPrefactor Pointer to the SimpleGrid of stress mat3
       * @param extents - The number of internal grid points on the current processor
       * @param offsets - The global mapping of the local (0,0,0) coordinate into the global grid index
       *
       * @return None.
       */
      void calculateStressPrefactor(SimpleGrid<Matrix3>* stressPrefactor,
                                    const IntVector& extents,
                                    const IntVector& offset);

      /**
       * @brief Generates split grid vector.
       *        Generates the vector of points from 0..K/2 in the first half of the array, followed by -K/2..-1
       * @param mPrime - A reference to the vector to populate
       * @param kMax - Maximum number of grid points for direction
       * @param spline - CenteredCardinalBSpline that determines the number of wrapping points necessary
       * @return std::vector<double> of (0..[m=K/2],[K/2-K]..-1);
       */
      inline void generateMPrimeVector(std::vector<double>& mPrime, unsigned int kMax) const
      {
        size_t halfMax = kMax / 2;

        for (size_t idx = 0; idx <= halfMax; ++idx) {
          mPrime[idx] = static_cast<double>(idx);
        }

        for (size_t Index = halfMax + 1; Index < kMax; ++Index) {
          mPrime[Index] = static_cast<double>(static_cast<double>(Index) - static_cast<int>(kMax));
        }
      }

      /**
       * @brief Generates reduced Fourier grid vector. Generates the vector of values i/K_i for i = 0...K-1
       * @param KMax - Maximum number of grid points for direction
       * @param InterpolatingSpline - CenteredCardinalBSpline that determines the number of wrapping points necessary
       * @return std::vector<double> of the reduced coordinates for the local grid along the input lattice direction
       */
      inline void generateMFractionalVector(std::vector<double>& mFractional,
                                            size_t kMax) const
      {
        double kMaxInv = 1.0 / static_cast<double>(kMax);
        for (size_t idx = 0; idx < kMax; ++idx) {
          mFractional[idx] = static_cast<double>(idx) * kMaxInv;
        }
      }

      /**
       * @brief Perform all calculations preceding the FFT transform of the charge grid to Fourier space.
       * @param const ProcessorGroup* pg -- All processors processing SPME patches
       * @param const PatchSubset* patches -- Patches to be processed by this thread
       * @param const MaterialSubset* materials -- Material subset belonging to this patch
       * @param DataWarehouse* old_dw -- Last time step's data warehouse
       * @param DataWarehouse* new_dw -- This time step's data warehouse
       * @return None
       */
      void calculatePreTransform(const ProcessorGroup* pg,
                                 const PatchSubset* patches,
                                 const MaterialSubset* materials,
                                 DataWarehouse* old_dw,
                                 DataWarehouse* new_dw);

      /**
       * @brief Perform all Fourier space calculations.
       * @param const ProcessorGroup* pg -- All processors processing SPME patches
       * @param const PatchSubset* patches -- Patches to be processed by this thread
       * @param const MaterialSubset* materials -- Material subset belonging to this patch
       * @param DataWarehouse* old_dw -- Last time step's data warehouse
       * @param DataWarehouse* new_dw -- This time step's data warehouse
       * @return None
       */
      void calculateInFourierSpace(const ProcessorGroup* pg,
                                   const PatchSubset* patches,
                                   const MaterialSubset* materials,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw);

      /**
       * @brief Perform calculations proceeding the FFT transform from Fourier to real space.
       * @param const ProcessorGroup* pg -- All processors processing SPME patches
       * @param const PatchSubset* patches -- Patches to be processed by this thread
       * @param const MaterialSubset* materials -- Material subset belonging to this patch
       * @param DataWarehouse* old_dw -- Last time step's data warehouse
       * @param DataWarehouse* new_dw -- This time step's data warehouse
       * @return None
       */
      void calculatePostTransform(const ProcessorGroup* pg,
                                  const PatchSubset* patches,
                                  const MaterialSubset* materials,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw);
      /**
       * @brief Perform necessary operation to transform Q grid to fourier space
       * @param const ProcessorGroup* pg -- All processors processing SPME patches
       * @param const PatchSubset* patches -- Patches to be processed by this thread
       * @param const MaterialSubset* materials -- Material subset belonging to this patch
       * @param DataWarehouse* old_dw -- Last time step's data warehouse
       * @param DataWarehouse* new_dw -- This time step's data warehouse
       * @return None
       *
       */
      void transformRealToFourier(const ProcessorGroup* pg,
                                  const PatchSubset* patches,
                                  const MaterialSubset* materials,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw);
      /**
       * @brief Perform necessary operation to transform Q grid to fourier space
       * @param const ProcessorGroup* pg -- All processors processing SPME patches
       * @param const PatchSubset* patches -- Patches to be processed by this thread
       * @param const MaterialSubset* materials -- Material subset belonging to this patch
       * @param DataWarehouse* old_dw -- Last time step's data warehouse
       * @param DataWarehouse* new_dw -- This time step's data warehouse
       * @return None
       *
       */
      void transformFourierToReal(const ProcessorGroup*pg,
                                  const PatchSubset* patches,
                                  const MaterialSubset* materials,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw);

      /**
       * @brief Reduce node-local Q data for the global reduction and FFT
       * @param None
       * @return None
       */
      void reduceNodeLocalQ(const ProcessorGroup* pg,
                            const PatchSubset* patches,
                            const MaterialSubset* materials,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw);
      /*
       * @brief Copy internal only patch local Q data to the node-local Q copy for global reduction
       * @param None
       * @return None
       */
      void copyToNodeLocalQ(const ProcessorGroup* pg,
                            const PatchSubset* patches,
                            const MaterialSubset* materials,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw);

      /**
       * @brief redistribute node-local Q data (force)
       * @param None
       * @return None
       */
      void distributeNodeLocalQ(const ProcessorGroup* pg,
                                const PatchSubset* patches,
                                const MaterialSubset* materials,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw);

      /**
       * @brief Checks for convergence of polarizability calculation
       * @param None
       * @return Bool - true if converged, false if not
       */
      bool checkConvergence() const;

      /**
       * @brief
       * @param
       * @return
       */
      inline bool getPolarizableCalculation() const
      {
        return d_polarizable;
      }

      // Values fixed on instantiation
      ElectrostaticsType d_electrostaticMethod;         //!< Implementation type for long range electrostatics
      MDSystem* d_system;                               //!< A handle to the MD simulation system object
      MDLabel* d_lb;                                    //!< A handle on the set of MD specific labels
      double d_ewaldBeta;						                    //!< The Ewald calculation damping coefficient
      bool d_polarizable;				                    	  //!< Use polarizable Ewald formulation
      double d_polarizationTolerance;                   //!< Tolerance threshold for polarizable system
      IntVector d_kLimits;                              //!< Number of grid divisions in each direction
      int d_maxPolarizableIterations;                   //!< Max number of polarization iterations to do
      ShiftedCardinalBSpline d_interpolatingSpline;     //!< Spline object to hold info for spline calculation
      SimpleGrid<dblcomplex>* d_Q_nodeLocal;            //!< The local version of the global Q grid
      SimpleGrid<dblcomplex>* d_Q_nodeLocalScratch;     //!< Scratch version of the global Q grid - for Allreduce

      struct LocalFFTData {
          fftw_complex* complexData;
          ptrdiff_t numElements;
          ptrdiff_t startAddress;
      };

      fftw_plan d_forwardPlan;                          //!< Forward FFTW MPI plan
      fftw_plan d_backwardPlan;                         //!< Reverse FFTW MPI plan
      LocalFFTData d_localFFTData;                      //!< The local portion of the global 3D FFT data

      // Variables we'll get from the MDSystem instance to make life easier
      Matrix3 d_unitCell;           //!< Unit cell lattice parameters
      Matrix3 d_inverseUnitCell;    //!< Inverse lattice parameters
      double d_systemVolume;        //!< Volume of the unit cell

      std::map<int, SPMEPatch*> d_spmePatchMap;  //!< These are the pieces of the K-space grid, map to Uintah patches

      Mutex d_Qlock;               //!< for local reductions on d_Q_nodeLocal (contention on overlapping ghost cells)
      mutable CrowdMonitor d_spmeLock;

  };

}  // End namespace Uintah

#endif
