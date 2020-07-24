/*
 * MPMBase.h
 *
 *  Created on: Jul 13, 2020
 *      Author: jbhooper
 */

#ifndef SRC_CCA_COMPONENTS_MPM_MPMBASE_H_
#define SRC_CCA_COMPONENTS_MPM_MPMBASE_H_

#include <CCA/Components/Application/ApplicationCommon.h>

#include <CCA/Ports/DataWarehouseP.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Util/DebugStream.h>

#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>

namespace Uintah {

	class ProcessorGroup;

	class MPMFlags;
	class MPMLabel;

	class MPMBase:  public ApplicationCommon {
	public:
		MPMBase(const 	ProcessorGroup 		* 	myWorld
			   ,		MaterialManagerP		materialManager
			   );
		virtual ~MPMBase();

		// Interfaces inherited form ApplicationCommon
		virtual void problemSetup(const	ProblemSpecP	&	probSpec
								 ,const	ProblemSpecP	&	restart_probSpec
								 ,		GridP			&	grid);

		virtual	void scheduleInitialize(const	LevelP		&	level
									   ,		SchedulerP	&	scheduler);

		virtual void scheduleRestartInitialize(const	LevelP		&	level
											  ,			SchedulerP	&	scheduler);

		virtual	void scheduleComputeStableTimeStep(const	LevelP		&	level
												  ,			SchedulerP	&	scheduler);



	private:
		MPMFlags	*	m_flags		=	nullptr;

		static const std::string m_className = "MPMBase";

	protected:
		Ghost::GhostType 	particle_ghost_type{Ghost::None};
		int					particle_ghost_layer{0};

	public:
		// Particle state
		std::vector<std::vector<const VarLabel* > > d_particleState;
		std::vector<std::vector<const VarLabel* > > d_particleState_preReloc;

		std::vector<std::vector<const VarLabel* > > d_cohesiveZoneState;
		std::vector<std::vector<const VarLabel* > > d_cohesiveZoneState_preReloc;

		MPMLabel*	m_labels {nullptr};


	};
}


#endif /* SRC_CCA_COMPONENTS_MPM_MPMBASE_H_ */
