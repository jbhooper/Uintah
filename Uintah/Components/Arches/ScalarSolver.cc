//----- ScalarSolver.cc ----------------------------------------------

/* REFERENCED */
static char *id="@(#) $Id$";

#include <Uintah/Components/Arches/ScalarSolver.h>
#include <Uintah/Components/Arches/RBGSSolver.h>
#include <Uintah/Components/Arches/Discretization.h>
#include <Uintah/Components/Arches/Source.h>
#include <Uintah/Components/Arches/BoundaryCondition.h>
#include <Uintah/Components/Arches/TurbulenceModel.h>
#include <Uintah/Components/Arches/PhysicalConstants.h>
#include <Uintah/Exceptions/InvalidValue.h>
#include <Uintah/Interface/Scheduler.h>
#include <Uintah/Interface/ProblemSpec.h>
#include <Uintah/Grid/Level.h>
#include <Uintah/Grid/Patch.h>
#include <Uintah/Grid/Task.h>
#include <SCICore/Util/NotFinished.h>
#include <Uintah/Components/Arches/Arches.h>

using namespace Uintah::ArchesSpace;
using namespace std;

//****************************************************************************
// Private constructor for PressureSolver
//****************************************************************************
ScalarSolver::ScalarSolver()
{
}

//****************************************************************************
// Default constructor for PressureSolver
//****************************************************************************
ScalarSolver::ScalarSolver(TurbulenceModel* turb_model,
			   BoundaryCondition* bndry_cond,
			   PhysicalConstants* physConst) : 
                                 d_turbModel(turb_model), 
                                 d_boundaryCondition(bndry_cond),
				 d_physicalConsts(physConst),
				 d_generation(0)
{
  d_scalarLabel = scinew VarLabel("scalar",
				    CCVariable<Vector>::getTypeDescription() );
  // BB : (tmp) velocity is set as CCVariable (should be FCVariable)
  d_velocityLabel = scinew VarLabel("velocity",
				    CCVariable<Vector>::getTypeDescription() );
  d_densityLabel = scinew VarLabel("density",
				   CCVariable<double>::getTypeDescription() );
  d_viscosityLabel = scinew VarLabel("viscosity",
				     CCVariable<double>::getTypeDescription() );
  d_scalarCoefLabel = scinew VarLabel("scalarCoeff",
				   CCVariable<Vector>::getTypeDescription() );
  d_scalarLinSrcLabel = scinew VarLabel("scalarLinearSource",
				     CCVariable<Vector>::getTypeDescription() );
  d_scalarNonLinSrcLabel = scinew VarLabel("scalarNonlinearSource",
					CCVariable<Vector>::getTypeDescription() );
}

//****************************************************************************
// Destructor
//****************************************************************************
ScalarSolver::~ScalarSolver()
{
}

//****************************************************************************
// Problem Setup
//****************************************************************************
void 
ScalarSolver::problemSetup(const ProblemSpecP& params)
{
  ProblemSpecP db = params->findBlock("MixtureFractionSolver");
  string finite_diff;
  db->require("finite_difference", finite_diff);
  if (finite_diff == "second") 
    d_discretize = scinew Discretization();
  else {
    throw InvalidValue("Finite Differencing scheme "
		       "not supported: " + finite_diff);
    //throw InvalidValue("Finite Differencing scheme "
	//	       "not supported: " + finite_diff, db);
  }
  // make source and boundary_condition objects
  d_source = scinew Source(d_turbModel, d_physicalConsts);
  string linear_sol;
  db->require("linear_solver", linear_sol);
  if (linear_sol == "linegs")
    d_linearSolver = scinew RBGSSolver();
  else {
    throw InvalidValue("linear solver option"
		       " not supported" + linear_sol);
    //throw InvalidValue("linear solver option"
	//	       " not supported" + linear_sol, db);
  }
  d_linearSolver->problemSetup(db);
}

//****************************************************************************
// Schedule solve of linearized scalar equation
//****************************************************************************
void 
ScalarSolver::solve(const LevelP& level,
		    SchedulerP& sched,
		    DataWarehouseP& old_dw,
		    DataWarehouseP& new_dw,
		    double time, double delta_t, int index)
{
  //create a new data warehouse to store matrix coeff
  // and source terms. It gets reinitialized after every 
  // pressure solve.
  DataWarehouseP matrix_dw = sched->createDataWarehouse(d_generation);
  ++d_generation;

  //computes stencil coefficients and source terms
  sched_buildLinearMatrix(level, sched, new_dw, matrix_dw, delta_t, index);
    
  d_linearSolver->sched_scalarSolve(level, sched, new_dw, matrix_dw, index);
    
}

//****************************************************************************
// Schedule build of linear matrix
//****************************************************************************
void 
ScalarSolver::sched_buildLinearMatrix(const LevelP& level,
				      SchedulerP& sched,
				      DataWarehouseP& old_dw,
				      DataWarehouseP& new_dw,
				      double delta_t, int index)
{
  for(Level::const_patchIterator iter=level->patchesBegin();
      iter != level->patchesEnd(); iter++){
    const Patch* patch=*iter;
    {
      // steve: requires two arguments
      //Task* tsk = scinew Task("ScalarSolver::BuildCoeff",
	//		      patch, old_dw, new_dw, this,
	//		      Discretization::buildLinearMatrix,
	//		      delta_t, index);
      Task* tsk = scinew Task("ScalarSolver::BuildCoeff",
			      patch, old_dw, new_dw, this,
			      &ScalarSolver::buildLinearMatrix,
			      delta_t, index);

      int numGhostCells = 0;
      int matlIndex = 0;
      tsk->requires(old_dw, d_scalarLabel, matlIndex, patch, Ghost::None,
		    numGhostCells);
      tsk->requires(old_dw, d_velocityLabel, matlIndex, patch, Ghost::None,
		    numGhostCells);
      tsk->requires(old_dw, d_densityLabel, matlIndex, patch, Ghost::None,
		    numGhostCells);
      tsk->requires(old_dw, d_viscosityLabel, matlIndex, patch, Ghost::None,
		    numGhostCells);

      /// requires convection coeff because of the nodal
      // differencing
      // computes all the components of velocity
      // added one more argument of index to specify scalar component
      tsk->computes(new_dw, d_scalarCoefLabel, matlIndex, patch);
      tsk->computes(new_dw, d_scalarLinSrcLabel, matlIndex, patch);
      tsk->computes(new_dw, d_scalarNonLinSrcLabel, matlIndex, patch);

      sched->addTask(tsk);
    }

  }
}

//****************************************************************************
// Actually build of linear matrix
//****************************************************************************
void ScalarSolver::buildLinearMatrix(const ProcessorContext* pc,
				     const Patch* patch,
				     DataWarehouseP& old_dw,
				     DataWarehouseP& new_dw,
				     double delta_t, int index)
{
  // compute ith componenet of velocity stencil coefficients
  d_discretize->calculateScalarCoeff(pc, patch, old_dw,
				       new_dw, delta_t, index);
  d_source->calculateScalarSource(pc, patch, old_dw,
				    new_dw, delta_t, index);
#ifdef WONT_COMPILE_YET
  d_boundaryCondition->scalarBC(pc, patch, old_dw,
				  new_dw, delta_t, index);
#endif
  // similar to mascal
  d_source->modifyScalarMassSource(pc, patch, old_dw,
				   new_dw, delta_t, index);
  d_discretize->calculateScalarDiagonal(pc, patch, old_dw,
				     new_dw, index);
}

//
// $Log$
// Revision 1.3  2000/06/04 22:40:15  bbanerje
// Added Cocoon stuff, changed task, require, compute, get, put arguments
// to reflect new declarations. Changed sub.mk to include all the new files.
//
//

