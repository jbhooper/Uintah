/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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

#ifndef UINTAH_CCA_COMPONENTS_PARENT_MULTISCALESWITCHER_H
#define UINTAH_CCA_COMPONENTS_PARENT_MULTISCALESWITCHER_H

#include <CCA/Ports/SimulationInterface.h>

#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Parallel/UintahParallelComponent.h>

#include <map>
#include <set>

namespace Uintah {

  class MultiScaleSwitcher : public UintahParallelComponent, public SimulationInterface {

  public:

             MultiScaleSwitcher( const ProcessorGroup* myworld, ProblemSpecP& ups, bool doAMR, const std::string & uda );

    virtual ~MultiScaleSwitcher();

    virtual void preGridProblemSetup(const ProblemSpecP&        params,
                                           GridP&               grid,
                                           SimulationStateP&    state);

    virtual void problemSetup( const ProblemSpecP     & params,
                               const ProblemSpecP     & restart_prob_spec,
                                     GridP            & grid,
                                     SimulationStateP & );

    virtual void outputProblemSpec(ProblemSpecP& ps);
    virtual void outputPS(Dir& dir);
    virtual void scheduleInitialize(            const LevelP& level, SchedulerP& sched );
    virtual void scheduleRestartInitialize(     const LevelP& level, SchedulerP& sched );
    virtual void scheduleComputeStableTimestep( const LevelP& level, SchedulerP& sched );
    virtual void scheduleTimeAdvance(           const LevelP& level, SchedulerP& sched );

    virtual void scheduleSwitchTest(            const LevelP& level, SchedulerP& sched );
    virtual void scheduleInitNewVars(           const LevelP& level, SchedulerP& sched );
    virtual void scheduleCarryOverVars(         const LevelP& level, SchedulerP& sched );
    virtual void scheduleSwitchInitialization(  const LevelP& level, SchedulerP& sched );
    virtual void scheduleFinalizeTimestep(      const LevelP& level, SchedulerP& sched );

    virtual bool   needRecompile(double time, double delt, const GridP& grid);
    virtual void   restartInitialize();
    virtual bool   restartableTimesteps();
    virtual double recomputeTimestep(double);

    enum switchState {
      idle,
      switching
    };


  private:

    SimulationInterface* matchComponentToLevelset(const LevelP& level);

    bool matchIndependentStatusToLevelset(const LevelP& level);

    bool matchAMRStatusToLevelset(const LevelP& level);

    void switchTest( const ProcessorGroup * /*pg*/,
                     const PatchSubset    * patches,
                     const MaterialSubset * matls,
                           DataWarehouse  * old_dw,
                           DataWarehouse  * new_dw );

    void initNewVars( const ProcessorGroup * /*pg*/,
                      const PatchSubset    *    patches,
                      const MaterialSubset * matls,
                            DataWarehouse  *  old_dw,
                            DataWarehouse  *  new_dw );

    void carryOverVars( const ProcessorGroup * /*pg*/,
                        const PatchSubset    * patches,
                        const MaterialSubset * matls,
                              DataWarehouse  * old_dw,
                              DataWarehouse  * new_dw );
                    
    void readSwitcherState( const ProblemSpecP     & /*ps*/,
                                  SimulationStateP & state );

    ProblemSpecP d_master_ups;
    switchState  d_switchState;
    SchedulerP   d_subScheduler;

    // since tasks are scheduled per-level, we can't turn the switch flag off
    // until they all are done, and since we need to turn it off during compilation,
    // we need to keep track of which levels we've switched
    std::vector<bool> d_doSwitching;

    bool d_restarting;

    SimulationInterface* d_sim;
    SimulationStateP     d_sharedState;
    unsigned int         d_numComponents;
    unsigned int         d_componentIndex;
    std::vector<bool>    d_componentIsAMR;
    
    struct initVars {
      std::vector<std::string>        varNames;
      std::vector<std::string>        matlSetNames;
      std::vector<const MaterialSet*> matls;
      std::vector<int>                levels;
      std::vector<VarLabel*>          varLabels;
    };
    
    std::map<int, initVars*>                     d_initVars;
    std::set<const VarLabel*, VarLabel::Compare> d_computedVars;
    

    std::vector<std::string>          d_in_file;                  // contains the name of all the subcomponent inputfiles
    std::vector<std::string>          d_carryOverVars;
    std::vector<VarLabel*>            d_carryOverVarLabels;
    std::vector<MaterialSubset*>      d_carryOverVarMatls;
    std::vector<bool>                 d_carryOverFinestLevelOnly; // either all levels or finest only
    std::vector<std::vector<bool> >   d_doCarryOverVarPerLevel;   // size to numlevels

    std::map<std::string, int>        d_componentNameIndexMap;    // Maps a component name string ("mpm","md",etc)
                                                                  // to the local index from which it is referenced
    // disable copy and assignment
    MultiScaleSwitcher(const MultiScaleSwitcher&);
    MultiScaleSwitcher& operator=(const MultiScaleSwitcher&);
	 
  };

}

#endif  // UINTAH_CCA_COMPONENTS_PARENT_MULTISCALESWITCHER_H
