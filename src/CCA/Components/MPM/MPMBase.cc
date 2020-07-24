/*
 * MPMBase.cc
 *
 *  Created on: Jul 13, 2020
 *      Author: jbhooper
 */

#include <MPMBase.h>



namespace Uintah {

Dout mpm_base_doing("BaseMPMRoutines", "MPMBase", "Track the action being performed by the MPMBase class instance.", false);
Dout mpm_base_schedule("BaseMPMScheduling", "MPMBase", "Track the scheduling of MPMBase tasks.", false);

	MPMBase::MPMBase(const	ProcessorGroup 		* 	myWorld
					,		MaterialManagerP		materialManager
					):ApplicationCommon(myWorld, materialManager)
	{

	}


	void MPMBase::problemSetup(const 	ProblemSpecP	&	probSpec
						  	  ,const	ProblemSpecP	&	restart_probSpec
							  ,			GridP			&	grid) {
	// Set up debug output file.

	std::string this_routine = m_className + "::problemSetup";
	std::stringstream messageOut;
	messageOut << "Doing " << this_routine << std::endl;
	DOUT(mpm_base_doing, messageOut.str());

	m_scheduler->setPositionVar(m_labels->pXLabel);

	ProblemSpecP restart_matSpec = nullptr;
	ProblemSpecP matSpec = probSpec->findBlockWithOutAttribute("MaterialProperties");

	bool isRestarting = false;
	if (matSpec) {
		restart_matSpec = matSpec;
	} else if (restart_probSpec) {
		isRestarting = true;
		restart_matSpec = restart_probSpec;
	} else {
		restart_matSpec = probSpec;
	}

	ProblemSpecP mpm_matSpec = restart_matSpec->findBlock("MPM");
	if (!mpm_matSpec) {
		std::ostringstream error_msg;
		error_msg << "ERROR:  " << this_routine
				  << " Missing MPM section in the input file." << std::endl;
		throw ProblemSetupException(error_msg.str(), __FILE__, __LINE__);
	}

	// Read all MPM flags (check MPMFlags.cc)
	m_flags->readMPMFlags(restart_matSpec, m_output);
	if (m_flags->d_integrator_type == "implicit") {
		std::ostringstream error_msg;
		error_msg << "ERROR:  Implicit MPM is not currently supported via the MPMBase class."
				  << std::endl;
		throw ProblemSetupException(error_msg.str(), __FILE__, __LINE__);
	}

	// Read in boundary tractions (Stress; generic read in boundary loads)
	std::vector<std::string>::const_iterator face_it;
	std::vector<std::string> & face_list = m_flags->d_bndy_face_txt_list;
	for(face_it = face_list.begin(); face_it != face_list.end(); ++face_it) {
		Patch::FaceType face = Patch::invalidFace;
		for (Patch::FaceType ft = Patch::startFace; ft <= Patch::endFace; ft=Patch::nextFace(ft)) {
			if (Patch::getFaceName(ft) == *face_it) face = ft;
		}
		if (face != Patch::invalidFace) {
			d_bndy_traction_faces.push_back(face);
		}
	}

	// Parse AMR related section from the input file
	ProblemSpecP amr_probSpec = probSpec->findBlock("AMR");
	if (amr_probSpec) {
		ProblemSpecP mpm_amr_probSpec = amr_probSpec->findBlock("MPM");
		if (!mpm_amr_probSpec) {
			std::ostringstream error_msg;
			error_msg << "ERROR -->  BaseMPM:\n"
					  <<"\tMPM section is missing from the AMR section of the input file." << std::endl;
			throw ProblemSetupException(error_msg.str(), __FILE__, __LINE__);
		}
		mpm_amr_probSpec->getWithDefault("min_grid_level", m_flags->d_minGridLevel, 0);
		mpm_amr_probSpec->getWithDefault("max_grid_level", m_flags->d_maxGridLevel, 1000);
		ProblemSpecP refine_spec = mpm_amr_probSpec->findBlock("Refinement_Criteria_Thresholds");

		// Pull ou tthe refinment threshold criteria
		if (refine_spec != nullptr) {
			for (ProblemSpecP var_spec = refine_spec->findBlock("Variable") ;
				 var_spec != nullptr ;
				 var_spec = var_spec->findNextBlock("Variable")) {
				thresholdVar data;
				std::string name, value, matl;

			}
		}
	}


	}
}
