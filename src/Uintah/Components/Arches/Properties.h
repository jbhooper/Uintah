
#ifndef Uintah_Component_Arches_Properties_h
#define Uintah_Component_Arches_Properties_h

/*
 * Placeholder - nothing here yet
 */

#include <Uintah/Parallel/UintahParallelComponent.h>
#include <Uintah/Interface/CFDInterface.h>
#include <Uintah/Grid/Region.h>
#include <Uintah/Parallel/ProcessorContext.h>
namespace Uintah {
namespace ArchesSpace {

class Properties {
public:
    Properties();
    ~Properties();

    void problemSetup(const ProblemSpecP& params);
    void sched_computeProps(const LevelP& level,
			    SchedulerP&, DataWarehouseP& old_dw,
			    DataWarehouseP& new_dw);
    int getNumMixVars() const{
      return d_numMixingVars;
    }
private:
    static const int MAX_MIXSTREAMS = 10;
    void computeProps(const ProcessorContext*,
		      const Region* region,
		      const DataWarehouseP& old_dw,
		      DataWarehouseP& new_dw);
    Properties(const Properties&);
    Properties& operator=(const Properties&);
    int d_numMixingVars;
    double d_denUnderrelax;
//temp soln...limits num of mixing streams
    struct Stream {
      double d_density;
      double d_temperature;
      Stream();
      void problemSetup(ProblemSpecP&);
    };
    Stream d_streams[MAX_MIXSTREAMS]; 
};

} // end namespace ArchesSpace
} // end namespace Uintah


#endif

