
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include <Uintah/Components/ICE/ICE.h>
#include <Uintah/Interface/DataWarehouse.h>
#include <Uintah/Grid/Grid.h>
#include <Uintah/Grid/Task.h>
#include <Uintah/Grid/Level.h>
#include <Uintah/Interface/Scheduler.h>
#include <Uintah/Grid/CCVariable.h>
#include <Uintah/Grid/NCVariable.h>
#include <Uintah/Grid/FCVariable.h>
#include <Uintah/Interface/ProblemSpec.h>
#include <Uintah/Grid/Patch.h>
#include <Uintah/Grid/CellIterator.h>
#include <Uintah/Grid/SoleVariable.h>
#include <SCICore/Geometry/Vector.h>
#include <SCICore/Geometry/IntVector.h>
#include <Uintah/Grid/VarLabel.h>
#include <Uintah/Grid/VarTypes.h>
#include <SCICore/Malloc/Allocator.h>

#include "ICE_switches.i"
#include "nrutil+.h"
#include "functionDeclare.h"
#include "parameters.h"
#include "switches.h"
#include "macros.h"
#include "cpgplot.h" /*must have this for plotting to work   */

using SCICore::Geometry::Vector;
using SCICore::Geometry::IntVector;
using std::cerr;
using std::endl;

/* ---------------------------------------------------------------------
GENERAL INFORMATION

 FILE NAME:  ICE.cc
 Purpose:    This is the main component for the Uintah ICE cfd code. 
.
History: 
Version   Programmer         Date       Description                      
     -------   ----------         ----       -----------                 
        1.0     Todd Harman       02/22/99                               
.                                                                    
    Programming Conventions
        i, j, k         Loop indices for the x, y, z directions respectively
        f               is a loop index for face-centered values.
        m               Loop index for the different materials
.
                                 ________ 
                                /  1    /|
                               /_______/ |
                              |       | ______(3)
                       (4)____| I,J,K |  |     
                              |       | /      
                              |_______|/
                                  |               (6) = back face
                                 (2)              (5) = front face
.
 STEPS:
    - Set some eviromnental variables required for PGPLOT
    - Initialize some variables that are mainly used in testing
    - MEMORY SECTION: Allocate the memory needed for all of the arrays
      For all of the face-centered arrays set equate the common face addresses
      [i][j][k][RIGHT][m] = [i-1][j][k][LEFT][m]
    - PROBLEM INITIALIZATION SECTION: Read in the input file, test the inputs,
      set the boundary condtions, generate the grid
    - MAIN LOOP
        to be filled in
    
 ---------------------------------------------------------------------*/   
extern "C" void audit();


using Uintah::ICESpace::ICE;

ICE::ICE()
{
    /*__________________________________
    *   Cell-centered variables
    *___________________________________*/
    delTLabel         = scinew VarLabel("delT",      delt_vartype::getTypeDescription() );
    
    press_CCLabel     = scinew VarLabel("press_CC",  CCVariable<double>::getTypeDescription() );
    press_CCLabel_0   = scinew VarLabel("press_CC_0",  CCVariable<double>::getTypeDescription() );
    press_CCLabel_1   = scinew VarLabel("press_CC_1",CCVariable<double>::getTypeDescription() );
    press_CCLabel_2   = scinew VarLabel("press_CC_2",CCVariable<double>::getTypeDescription() );
    press_CCLabel_3   = scinew VarLabel("press_CC_3",CCVariable<double>::getTypeDescription() );
    press_CCLabel_4   = scinew VarLabel("press_CC_4",CCVariable<double>::getTypeDescription() );
    press_CCLabel_5   = scinew VarLabel("press_CC_5",CCVariable<double>::getTypeDescription() );
    press_CCLabel_6_7 = scinew VarLabel("press_CC_6_7",CCVariable<double>::getTypeDescription() );

    rho_CCLabel       = scinew VarLabel("rho_CC",    CCVariable<double>::getTypeDescription() );
    rho_CCLabel_0     = scinew VarLabel("rho_CC_0",  CCVariable<double>::getTypeDescription() );
    rho_CCLabel_1     = scinew VarLabel("rho_CC_1",  CCVariable<double>::getTypeDescription() );
    rho_CCLabel_2     = scinew VarLabel("rho_CC_2",  CCVariable<double>::getTypeDescription() );
    rho_CCLabel_3     = scinew VarLabel("rho_CC_3",  CCVariable<double>::getTypeDescription() );
    rho_CCLabel_4     = scinew VarLabel("rho_CC_4",  CCVariable<double>::getTypeDescription() );
    rho_CCLabel_5     = scinew VarLabel("rho_CC_5",  CCVariable<double>::getTypeDescription() );
    rho_CCLabel_6_7   = scinew VarLabel("rho_CC_6_7",  CCVariable<double>::getTypeDescription() );
 
    temp_CCLabel      = scinew VarLabel("temp_CC",   CCVariable<double>::getTypeDescription() );
    temp_CCLabel_0    = scinew VarLabel("temp_CC_0", CCVariable<double>::getTypeDescription() );
    temp_CCLabel_1    = scinew VarLabel("temp_CC_1", CCVariable<double>::getTypeDescription() );
    temp_CCLabel_2    = scinew VarLabel("temp_CC_2", CCVariable<double>::getTypeDescription() );
    temp_CCLabel_3    = scinew VarLabel("temp_CC_3", CCVariable<double>::getTypeDescription() );
    temp_CCLabel_4    = scinew VarLabel("temp_CC_4", CCVariable<double>::getTypeDescription() );
    temp_CCLabel_5    = scinew VarLabel("temp_CC_5", CCVariable<double>::getTypeDescription() );
    temp_CCLabel_6_7  = scinew VarLabel("temp_CC_6_7", CCVariable<double>::getTypeDescription() );
 
    vel_CCLabel       = scinew VarLabel("vel_CC",    CCVariable<Vector>::getTypeDescription() );
    vel_CCLabel_0     = scinew VarLabel("vel_CC_0",  CCVariable<Vector>::getTypeDescription() );
    vel_CCLabel_1     = scinew VarLabel("vel_CC_1",  CCVariable<Vector>::getTypeDescription() );
    vel_CCLabel_2     = scinew VarLabel("vel_CC_2",  CCVariable<Vector>::getTypeDescription() );
    vel_CCLabel_3     = scinew VarLabel("vel_CC_3",  CCVariable<Vector>::getTypeDescription() );
    vel_CCLabel_4     = scinew VarLabel("vel_CC_4",  CCVariable<Vector>::getTypeDescription() );
    vel_CCLabel_5     = scinew VarLabel("vel_CC_5",  CCVariable<Vector>::getTypeDescription() );
    vel_CCLabel_6_7   = scinew VarLabel("vel_CC_6_7",  CCVariable<Vector>::getTypeDescription() );
 
    cv_CCLabel        = scinew VarLabel("cv_CC",     CCVariable<double>::getTypeDescription() );
    div_velfc_CCLabel = scinew VarLabel("div_velfc_CC", CCVariable<double>::getTypeDescription() );

  // Face centered variables
    vel_FCLabel       = scinew VarLabel("vel_FC",    FCVariable<Vector>::getTypeDescription() );
    press_FCLabel     = scinew VarLabel("press_FC",  FCVariable<double>::getTypeDescription() );
    tau_FCLabel       = scinew VarLabel("tau_FC",    FCVariable<Vector>::getTypeDescription() );

    /*__________________________________
    *   Plotting variables
    *___________________________________*/
    stat = putenv("PGPLOT_DIR=" PGPLOT_DIR);
    stat = putenv("PGPLOT_I_AM_HERE=0");              
    stat = putenv("PGPLOT_PLOTTING_ON_OFF=1");
    stat = putenv("PGPLOT_OPEN_NEW_WINDOWS=1");  
    
    /*__________________________________
     *   Allocate memory for the arrays
     *___________________________________*/
#include "allocate_memory.i"
    
}





ICE::~ICE()
{
  /*__________________________________
   *   Now deallocate memory
   *___________________________________*/
  fprintf(stderr,"Now deallocating memory");
#include "free_memory.i"
}




/* ---------------------------------------------------------------------
   GENERAL INFORMATION
   Function:  ICE::problemSetup--
   Filename:  ICE_actual.cc
   Purpose:   
   
   History: 
   Version   Programmer         Date       Description                      
   -------   ----------         ----       -----------                 
   1.0     John Schmidt      06/23/00                              
   _____________________________________________________________________*/
void ICE::problemSetup(const ProblemSpecP& prob_spec, GridP&,
		       SimulationStateP&)
{
  
  double  viscosity, thermal_conductivity, specific_heat, speed_of_sound;
  double ideal_gas_constant, d_gamma;
    
  printSwitch = 1;    
  t           = 0.0;  
  m           = 1;
  fileNum     = 1;
  /*__________________________________
   *   Read in from the spec file
   *___________________________________*/
  ProblemSpecP mat_ps = prob_spec->findBlock("MaterialProperties");
  
  ProblemSpecP ice_mat_ps = mat_ps->findBlock("ICE");
  
  for (ProblemSpecP ps = ice_mat_ps->findBlock("material"); ps != 0;
       ps = ps->findNextBlock("material") ) 
 {
    ps->require(    "viscosity",            viscosity);
    ps->require(    "thermal_conductivity", thermal_conductivity);
    ps->require(    "specific_heat",        specific_heat);
    ps->require(    "speed_of_sound",       speed_of_sound);
    ps->require(    "ideal_gas_constant",   ideal_gas_constant);
    ps->require(    "gamma",                d_gamma);
  }
  
  cerr << "viscosity " << viscosity << endl;
  cerr << "thermal_conductivity " << thermal_conductivity << endl;
  cerr << "specific_heat " << specific_heat << endl;
  cerr << "speed_of_sound " << speed_of_sound << endl;
  cerr << "ideal_gas_constant " << ideal_gas_constant << endl;
  cerr << "gamma " << d_gamma << endl;
  
  audit();
  
  /*______________________________________________________________________
   *
   *  P  R  O  B  L  E  M     I  N  I  T  I  A  L  I  Z  A  T  I  O  N  
   *  - read input file
   *   - test the input variables
   *   - Equate the address of the face centered variables
   *   - Generate a grid
   *   - zero all of the face-centered arrays
   *   
   *                  
   * -----------------------------------------------------------------------  */
  
        readInputFile(   &xLoLimit,      &yLoLimit,      &zLoLimit,     
                        &xHiLimit,      &yHiLimit,      &zHiLimit,
                        &delX,          &delY,          &delZ,
                        uvel_CC,        vvel_CC,        wvel_CC, 
                        Temp_CC,        press_CC,       rho_CC,
                        scalar1_CC,     scalar2_CC,     scalar3_CC,
                        viscosity_CC,   thermalCond_CC, cv_CC,
                        R,              gamma,
                        &t_final,       t_output_vars,  delt_limits,
                        output_file_basename,           output_file_desc,       
                        grav,           speedSound,
                        BC_inputs,      BC_Values,      &CFL,
                        &nMaterials);      
    
    testInputFile(      xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        delX,           delY,           delZ,
                        Temp_CC,        press_CC,       rho_CC,
                        viscosity_CC,   thermalCond_CC, cv_CC,
                        speedSound,      
                        t_final,        t_output_vars,  delt_limits,
                        BC_inputs,      printSwitch,    CFL,
                        nMaterials); 
                   
    definition_of_different_physical_boundary_conditions(              
                        BC_inputs,      BC_types,       BC_float_or_fixed,
                        BC_Values,      nMaterials  );  
                        
/*__________________________________
* Now make sure that the face centered
* values know about each other.
* for example 
* [i][j][k][RIGHT][m] = [i-1][j][k][LEFT][m]
*___________________________________*/  

    equate_ptr_addresses_adjacent_cell_faces(              
                        x_FC,           y_FC,           z_FC,
                        uvel_FC,        vvel_FC,        wvel_FC,
                        press_FC,
                        tau_X_FC,       tau_Y_FC,       tau_Z_FC,
                        nMaterials);   

    /*__________________________________
    * Generate a grid
    *___________________________________*/ 
    generateGrid(       xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        delX,           delY,           delZ,
                        x_CC,           y_CC,           z_CC,   Vol_CC,  
                        x_FC,           y_FC,           z_FC );
    /*__________________________________
    *   zero the face-centered arrays
    *___________________________________*/
    zero_arrays_6d(
                        xLoLimit,       yLoLimit,       zLoLimit,             
                        xHiLimit,       yHiLimit,       zHiLimit,
                        1,              N_CELL_FACES,
                        1,              nMaterials,     
                        7,             
                        uvel_FC,        vvel_FC,        wvel_FC,
                        press_FC,
                        tau_X_FC,       tau_Y_FC,       tau_Z_FC);                         
    stat = putenv("PGPLOT_PLOTTING_ON_OFF=1");
  
  
  /*__________________________________
   *   overide the initial conditions
   *___________________________________*/
#if switchOveride_Initial_Conditions                               
#include "overide_initial_conds.i"
#endif 
  
  /*__________________________________
   *  If desired plot the inputs
   *___________________________________*/
#if switchDebug_main_input
#define switchInclude_main_1 1
#include "debugcode.i"
#undef switchInclude_main_1
#endif 
  
  /*__________________________________
   *   For the first time through
   *   set some variables
   *___________________________________*/
  delt    = delt_limits[3];              
  t       = delt;
  fprintf(stderr,"\nInitial time %f, timestep is %f\n",t,delt); 
  
  
#if 0
  CCVariable<Vector> vel_CC;
  ds.put("vel_CC", vel_CCLabel,0,patch);
  
#else
  cerr << "put vel_CC not done\n";
#endif
  
}





/* --------------------------------------------------------------------- 
GENERAL INFORMATION
 Function:  ICE::actuallyInitialize--
 Filename:  ICE_actual.cc
 Purpose:   -allocate variables in the new DW
            - Convert the NR array data into the UCF format
            - Put the data into the DW
            

History: 
Version   Programmer         Date       Description                      
-------   ----------         ----       -----------                 
  1.0     John Schmidt   06/23/00                              
_____________________________________________________________________*/ 
void ICE::actuallyInitialize(
    const ProcessorGroup*,
    const Patch* patch,
    DataWarehouseP& /* old_dw */,
    DataWarehouseP& to_dw)
{
    int m;
     cerr <<"Doing actuallyInitialize . . ." << endl;

/*______________________________________________________________________
*     T  O     D  W     W  R  A  P  P  E  R  
*--------------------------------------*/
#if switch_UCF_initialize
for (m = 1; m<= nMaterials; m++)
{
    CCVariable<double>  press_cc;
    CCVariable<double>  rho_cc;
    CCVariable<double>  temp_cc;
    CCVariable<Vector>  vel_cc;
    
    to_dw->allocate(    press_cc, press_CCLabel,  m,patch);
    to_dw->allocate(    rho_cc,   rho_CCLabel,    m,patch);
    to_dw->allocate(    temp_cc,  temp_CCLabel,   m,patch);
    to_dw->allocate(    vel_cc,   vel_CCLabel,    m,patch);

    /*__________________________________
    *   UCF NR
    *___________________________________*/
    ICE::after_each_step_wrapper(   patch,
                        m,
                        press_cc,  press_CC,
                        rho_cc,    rho_CC,
                        temp_cc,   Temp_CC,
                        vel_cc,    uvel_CC,     vvel_CC,        wvel_CC);
 
    to_dw->put(         press_cc,  press_CCLabel,  m,patch);
    to_dw->put(         rho_cc,    rho_CCLabel,    m,patch);
    to_dw->put(         temp_cc,   temp_CCLabel,   m,patch);
    to_dw->put(         vel_cc,    vel_CCLabel,    m,patch);
}
#endif   
/*_______________________________________________________________________*/
}




/* --------------------------------------------------------------------- 
GENERAL INFORMATION
 Function:  ICE::actuallyComputeStableTimestep--
 Filename:  ICE_actual.cc
 Purpose:   Compute the stable time step based on the courant condition
            
History: 
Version   Programmer         Date       Description                      
-------   ----------         ----       -----------                 
  1.0     John Schmidt   06/23/00                              
_____________________________________________________________________*/ 
void ICE::actuallyComputeStableTimestep(const ProcessorGroup*,
    const Patch* patch,
    DataWarehouseP& fromDW,
    DataWarehouseP& toDW)
{
  bool include_ghost_cells = NO;
  cerr << "Doing actuallyComputeStableTimestep . . ." << endl;

    /*__________________________________
    * convert UCF data into NR arrays 
    *___________________________________*/
  CCVariable<Vector> vel_cc;    
  for (m = 1; m<= nMaterials; m++) {
    toDW->get(vel_cc, vel_CCLabel,m, patch,Ghost::None,0);

    ICE::convertUCFToNR_4d(patch,       vel_cc,
                        uvel_CC,        vvel_CC,        wvel_CC,
                        include_ghost_cells,
                        xLoLimit,       xHiLimit,
		          yLoLimit,       yHiLimit,
                        zLoLimit,       zHiLimit,
                        m);
    }

    //    int numMatls = d_sharedState->getNumMatls();
    int numMatls = nMaterials;
    
    double A,B, delt_CFL;
    
    for (int m = 1; m < numMatls; m++) {
#if 0
       Material* matl = d_sharedState->getMaterial(m);
       ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);
      if (ice_matl) {
	int matlindex = matl->getDWIndex();
	inf vfindex = matl->getVFIndex();
#endif
	//
	int matlindex = m;
	toDW->get(vel_cc, vel_CCLabel,matlindex, patch,Ghost::None,0);


	// Once we get the data on the grid, in this case the velocities,
	// now we can get an iterator (which will be the equivalent to
	// writing a triply nested loop (i,j,k) to loop over all of the
	// velocity values of the patch.  

	// We need to pass the getBox to the getCellIterator, since it
	// indicates the extents of the patch in question.

	// The iter.done() is just a test condition to see if we have
	// visited all of the cells.

	// Note:  This is the preferred way of traversing over the grid
	// (or cells).

	for (CellIterator iter = patch->getCellIterator(patch->getBox());
	     !iter.done();
	     iter++) {

	  // Get the patch spacing.
	  double delx = patch->dCell().x();
	  double dely = patch->dCell().y();

	  double fudge_factor = 1.0;
#if 1
	  A = fudge_factor*CFL*delx/fabs(vel_cc[*iter].x() + SMALL_NUM);
	  B = fudge_factor*CFL*dely/fabs(vel_cc[*iter].y() + SMALL_NUM);
	  cout << "A = " << A << " B = " << B << " delt_CFL = " << delt_CFL << endl;

	  delt_CFL = DMIN(A,delt_CFL);
	  delt_CFL = DMIN(B,delt_CFL);
#endif
	  // Do other steps

	  //	}
	}

      }
      cout << "delt = " << delt << endl;
      delt = delt_CFL;
      cout << "delt = " << delt << endl;

  
    /*__________________________________
    *   Find the new time step based on the
    *   Courant condition
    *___________________________________*/        
    find_delta_time_based_on_CC_vel(
                        xLoLimit,        yLoLimit,      zLoLimit,
                        xHiLimit,        yHiLimit,      zHiLimit,
                        &delt,           delt_limits,
                        delX,            delY,          delZ,
                        uvel_CC,         vvel_CC,       wvel_CC,
                        speedSound,      CFL,           nMaterials );
    cout << " delt computed = " << delt << endl;
    delt_vartype dt(delt);
    toDW->put(dt, delTLabel);
}



/* ---------------------------------------------------------------------
GENERAL INFORMATION
 Function:  ICE::actually_Top_of_main_loop--
 Filename:  ICE_actual.cc
 Purpose:   - Include the pgplot variables and set environmental vars
            - Update the face and cell centered variables in the ghost cells
            - Before you enter the main loop find the time step
            

History: 
Version   Programmer         Date       Description                      
-------   ----------         ----       -----------                 
  1.0     John Schmidt   06/23/00                              
_____________________________________________________________________*/ 
void ICE::actually_Top_of_main_loop(const ProcessorGroup*,
			const Patch* patch,
			DataWarehouseP& from_dw,
			DataWarehouseP& to_dw)
{
    int should_I_write_output;
    int m;    

/*______________________________________________________________________
*    F  R  O  M     D  W     W  R  A  P  P  E  R 
*--------------------------------------*/
#if switch_UCF_stepTop_of_main_loopOnOff 
    for (m = 1; m<= nMaterials; m++)
    {
        CCVariable<double>  press_cc;
        CCVariable<double>  rho_cc;
        CCVariable<double>  temp_cc;
        CCVariable<Vector>  vel_cc;
 
        from_dw->get(       press_cc,       press_CCLabel,  m, patch, Ghost::None, 0);
        from_dw->get(       rho_cc,         rho_CCLabel,    m, patch, Ghost::None, 0);
        from_dw->get(       temp_cc,        temp_CCLabel,   m, patch, Ghost::None, 0);
        from_dw->get(       vel_cc,         vel_CCLabel,    m, patch, Ghost::None, 0);
        /*__________________________________
        *   UCF NR
        *___________________________________*/
        ICE::before_each_step_wrapper(  patch,
                            m,
                            press_cc,       press_CC,
                            rho_cc,         rho_CC,
                            temp_cc,        Temp_CC,
                            vel_cc,
                            uvel_CC,        vvel_CC,        wvel_CC);
    }
#endif
/*______________________________________________________________________*/  

    fprintf(stderr,"\n\n________________________________\n");
    cerr << "Actually at top of main loop" << endl;
  
    should_I_write_output = Is_it_time_to_write_output( t, t_output_vars  ); 
  /*__________________________________
   * update the physical boundary conditions
   * and initialize some arrays
   *___________________________________*/                        
    update_CC_FC_physical_boundary_conditions( 
                        xLoLimit,       yLoLimit,       zLoLimit,             
                        xHiLimit,       yHiLimit,       zHiLimit,             
                        delX,           delY,           delZ,
                        BC_types,       BC_float_or_fixed,
                        BC_Values, 
                        nMaterials,     3,                 
                        uvel_CC,        UVEL,           uvel_FC,
                        vvel_CC,        VVEL,           vvel_FC,
                        wvel_CC,        WVEL,           wvel_FC);
  
    update_CC_physical_boundary_conditions( 
                        xLoLimit,       yLoLimit,       zLoLimit,             
                        xHiLimit,       yHiLimit,       zHiLimit,             
                        delX,           delY,           delZ,
                        BC_types,       BC_float_or_fixed,
                        BC_Values, 
                        nMaterials,     3,                 
                        Temp_CC,TEMP,   rho_CC,DENSITY, press_CC,PRESS);
                        
    zero_arrays_4d(
                        xLoLimit,       yLoLimit,       zLoLimit,             
                        xHiLimit,       yHiLimit,       zHiLimit,
                        1,              nMaterials,     8,             
                        mass_source,    delPress_CC,    int_eng_source,  
                        xmom_source,    ymom_source,    zmom_source,
                        Vol_L_CC,       mass_CC);


/*`==========TESTING==========*/ 
    /*__________________________________
    *   Find the new time step based on the
    *   Courant condition
    *___________________________________*/        
 /*    find_delta_time_based_on_CC_vel(
                        xLoLimit,        yLoLimit,      zLoLimit,
                        xHiLimit,        yHiLimit,      zHiLimit,
                        &delt,           delt_limits,
                        delX,            delY,          delZ,
                        uvel_CC,         vvel_CC,       wvel_CC,
                        speedSound,      CFL,           nMaterials ); */
 /*==========TESTING==========`*/ 
    
/*______________________________________________________________________
*     T  O     D  W     W  R  A  P  P  E  R  
*--------------------------------------*/
#if switch_UCF_stepTop_of_main_loopOnOff 
    for (m = 1; m<= nMaterials; m++)
    {
        CCVariable<double>  press_cc_0;
        CCVariable<double>  rho_cc_0;
        CCVariable<double>  temp_cc_0;
        CCVariable<Vector>  vel_cc_0;

        to_dw->allocate(    press_cc_0, press_CCLabel_0,  m,patch);
        to_dw->allocate(    rho_cc_0,   rho_CCLabel_0,    m,patch);
        to_dw->allocate(    temp_cc_0,  temp_CCLabel_0,   m,patch);
        to_dw->allocate(    vel_cc_0,   vel_CCLabel_0,    m,patch);

        /*__________________________________
        *   UCF NR
        *___________________________________*/
        ICE::after_each_step_wrapper(   patch,
                            m,
                            press_cc_0,  press_CC,
                            rho_cc_0,    rho_CC,
                            temp_cc_0,   Temp_CC,
                            vel_cc_0,    uvel_CC,     vvel_CC,        wvel_CC);

        to_dw->put(         press_cc_0,  press_CCLabel_0,  m,patch);
        to_dw->put(         rho_cc_0,    rho_CCLabel_0,    m,patch);
        to_dw->put(         temp_cc_0,   temp_CCLabel_0,   m,patch);
        to_dw->put(         vel_cc_0,    vel_CCLabel_0,    m,patch);
    }
#endif   
/*_______________________________________________________________________*/
cerr << "Actually at top of main loop 2" << endl;
  /*__________________________________
    *   Quite full warn remarks
    *___________________________________*/
    should_I_write_output = should_I_write_output;                   
}





/* ---------------------------------------------------------------------
                                                        S  T  E  P     1 
GENERAL INFORMATION
 Function:  ICE::actuallyStep1--
 Filename:  ICE_actual.cc
 Purpose:   STEP 1 
            compute the cell-centered pressure using the equation of state
            and the speed of sound 

History: 
Version   Programmer         Date       Description                      
-------   ----------         ----       -----------                 
  1.0     John Schmidt   06/23/00                              
_____________________________________________________________________*/ 
void ICE::actuallyStep1(const ProcessorGroup*,
			const Patch* patch,
			DataWarehouseP& from_dw,
			DataWarehouseP& to_dw)
{
  int m;
 
 
/*______________________________________________________________________
*    F  R  O  M     D  W     W  R  A  P  P  E  R 
*--------------------------------------*/
#if switch_UCF_step1OnOff  
    for (m = 1; m<= nMaterials; m++)
    {
        CCVariable<double>  press_cc;
        CCVariable<double>  rho_cc;
        CCVariable<double>  temp_cc;
        CCVariable<Vector>  vel_cc;
 
        from_dw->get(       press_cc,     press_CCLabel_0,m, patch, Ghost::None, 0);
        from_dw->get(       rho_cc,       rho_CCLabel_0,  m, patch, Ghost::None, 0);
        from_dw->get(       temp_cc,      temp_CCLabel_0, m, patch, Ghost::None, 0);
        from_dw->get(       vel_cc,       vel_CCLabel_0,  m, patch, Ghost::None, 0);
        /*__________________________________
        *   UCF NR
        *___________________________________*/
        ICE::before_each_step_wrapper(  patch,
                            m,
                            press_cc,     press_CC,
                            rho_cc,       rho_CC,
                            temp_cc,      Temp_CC,
                            vel_cc,
                            uvel_CC,        vvel_CC,        wvel_CC);
    }
#endif    
/*______________________________________________________________________*/
  /*__________________________________
   *  Use the equation of state to get
   *  P at the cell center
   *___________________________________*/
#if switch_step1_OnOff
    cerr << "Actually doing step 1 " << endl;

    // INPUTS:  R, rho_CC, Temp_CC
    // OUTPUTS: press_CC
    equation_of_state(
                    xLoLimit,       yLoLimit,       zLoLimit,
                    xHiLimit,       yHiLimit,       zHiLimit,
                    R,
                    press_CC,       rho_CC,         Temp_CC,
                    cv_CC,          nMaterials   );  
    // INPUTS: gamma, R, Temp_CC
    // OUTPUS: speedSound
    speed_of_sound(
                    xLoLimit,       yLoLimit,       zLoLimit,       
                    xHiLimit,       yHiLimit,       zHiLimit,       
                    gamma,          R,              Temp_CC,     
                    speedSound,     nMaterials   );
#endif

/*______________________________________________________________________
*     T  O     D  W     W  R  A  P  P  E  R  
*--------------------------------------*/
#if switch_UCF_step1OnOff
    for (m = 1; m<= nMaterials; m++)
    {
        CCVariable<double>  press_cc_1;
        CCVariable<double>  rho_cc_1;
        CCVariable<double>  temp_cc_1;
        CCVariable<Vector>  vel_cc_1;
 
        to_dw->allocate(    press_cc_1, press_CCLabel_1,  m,patch);
        to_dw->allocate(    rho_cc_1,   rho_CCLabel_1,    m,patch);
        to_dw->allocate(    temp_cc_1,  temp_CCLabel_1,   m,patch);
        to_dw->allocate(    vel_cc_1,   vel_CCLabel_1,    m,patch);

        /*__________________________________
        *   UCF NR
        *___________________________________*/
        ICE::after_each_step_wrapper(   patch,
                           m,
                           press_cc_1,  press_CC,
                           rho_cc_1,    rho_CC,
                           temp_cc_1,   Temp_CC,
                           vel_cc_1,    uvel_CC,     vvel_CC,        wvel_CC);
 
        to_dw->put(        press_cc_1,  press_CCLabel_1,  m,patch);
        to_dw->put(        rho_cc_1,    rho_CCLabel_1,    m,patch);
        to_dw->put(        temp_cc_1,   temp_CCLabel_1,   m,patch);
        to_dw->put(        vel_cc_1,    vel_CCLabel_1,    m,patch);
    }
#endif   

/*_______________________________________________________________________*/
}






/* ---------------------------------------------------------------------
                                                        S  T  E  P     2 
GENERAL INFORMATION
 Function:  ICE::actuallyStep2--
 Filename:  ICE_actual.cc
 Purpose:   - Compute the face-centered velocities using time n data
            - Compute the divergence of the face centered velocities
            - Compute the change in the cell centered pressure (explicit)   
History: 
Version   Programmer         Date       Description                      
-------   ----------         ----       -----------                 
  1.0     John Schmidt   06/23/00                              
_____________________________________________________________________*/
void ICE::actuallyStep2(const ProcessorGroup*,
			const Patch* patch,
			DataWarehouseP& from_dw,
			DataWarehouseP& to_dw)
{
  int m;
/*______________________________________________________________________
*    F  R  O  M     D  W     W  R  A  P  P  E  R 
*--------------------------------------*/
#if switch_UCF_step2OnOff  
    for (m = 1; m<= nMaterials; m++)
    {
        CCVariable<double>  press_cc_1;
        CCVariable<double>  rho_cc_1;
        CCVariable<double>  temp_cc_1;
        CCVariable<Vector>  vel_cc_1;
 
        from_dw->get(       press_cc_1,     press_CCLabel_1,m, patch, Ghost::None, 0);
        from_dw->get(       rho_cc_1,       rho_CCLabel_1,  m, patch, Ghost::None, 0);
        from_dw->get(       temp_cc_1,      temp_CCLabel_1, m, patch, Ghost::None, 0);
        from_dw->get(       vel_cc_1,       vel_CCLabel_1,  m, patch, Ghost::None, 0);
        /*__________________________________
        *   UCF NR
        *___________________________________*/
        ICE::before_each_step_wrapper(  patch,
                            m,
                            press_cc_1,     press_CC,
                            rho_cc_1,       rho_CC,
                            temp_cc_1,      Temp_CC,
                            vel_cc_1,
                            uvel_CC,        vvel_CC,        wvel_CC);
    }
#endif    
/*______________________________________________________________________*/
 
    /*__________________________________
    *   Now do the work
    *___________________________________*/    
    cerr << "Actually doing step 2" << endl;


    stat = putenv("PGPLOT_PLOTTING_ON_OFF=1"); 
    // INPUTS: uvel_CC, vvel_CC, wvel_CC
    // OUTPUTS: uvel_FC, vvel_FC, wvel_FC
    compute_face_centered_velocities( 
                    xLoLimit,       yLoLimit,       zLoLimit,
                    xHiLimit,       yHiLimit,       zHiLimit,
                    delX,           delY,           delZ,
                    delt,           
                    BC_types,       BC_float_or_fixed,
                    BC_Values,
                    rho_CC,         grav,           press_CC,
                    uvel_CC,        vvel_CC,        wvel_CC,
                    uvel_FC,        vvel_FC,        wvel_FC,
                    nMaterials ); 



    // INPUTS: delX, delY, delZ, uvel_FC, vvel_FC, wvel_FC
    // OUTPUTS: div_velFC_CC
    divergence_of_face_centered_velocity(  
                    xLoLimit,       yLoLimit,       zLoLimit,
                    xHiLimit,       yHiLimit,       zHiLimit,
                    delX,           delY,           delZ,
                    uvel_FC,        vvel_FC,        wvel_FC,
                    div_velFC_CC,   nMaterials); 
    stat = putenv("PGPLOT_PLOTTING_ON_OFF=1");


#if switch_step2_OnOff                        
    // INPUTS: delX, delY, delZ, div_velFC_CC, delPress_CC, press_CC, rho_CC,
    //         delt, speedSound               
    // OUTPUTS: ?? delPress_CC    Is this right?
    explicit_delPress
             (  
                    xLoLimit,       yLoLimit,       zLoLimit,
                    xHiLimit,       yHiLimit,       zHiLimit,
                    delX,           delY,           delZ,
                    div_velFC_CC,
                    delPress_CC,    press_CC,
                    rho_CC,         delt,           speedSound,
                    nMaterials );
    
    // INPUTS:
    // OUTPUTS:  
    update_CC_physical_boundary_conditions( 
                    xLoLimit,       yLoLimit,       zLoLimit,             
                    xHiLimit,       yHiLimit,       zHiLimit,             
                    delX,           delY,           delZ,
                    BC_types,       BC_float_or_fixed,
                    BC_Values, 
                    nMaterials,     1,                 
                    delPress_CC,    DELPRESS);
                                            
    #endif
/*______________________________________________________________________
*     T  O     D  W     W  R  A  P  P  E  R  
*--------------------------------------*/
#if switch_UCF_step2OnOff
    for (m = 1; m<= nMaterials; m++)
    {
        CCVariable<double>  press_cc_2;
        CCVariable<double>  rho_cc_2;
        CCVariable<double>  temp_cc_2;
        CCVariable<Vector>  vel_cc_2;
 
        to_dw->allocate(    press_cc_2, press_CCLabel_2,  m,patch);
        to_dw->allocate(    rho_cc_2,   rho_CCLabel_2,    m,patch);
        to_dw->allocate(    temp_cc_2,  temp_CCLabel_2,   m,patch);
        to_dw->allocate(    vel_cc_2,   vel_CCLabel_2,    m,patch);

        /*__________________________________
        *   UCF NR
        *___________________________________*/
        ICE::after_each_step_wrapper(   patch,
                           m,
                           press_cc_2,  press_CC,
                           rho_cc_2,    rho_CC,
                           temp_cc_2,   Temp_CC,
                           vel_cc_2,    uvel_CC,     vvel_CC,        wvel_CC);
 
        to_dw->put(        press_cc_2,  press_CCLabel_2,  m,patch);
        to_dw->put(        rho_cc_2,    rho_CCLabel_2,    m,patch);
        to_dw->put(        temp_cc_2,   temp_CCLabel_2,   m,patch);
        to_dw->put(        vel_cc_2,    vel_CCLabel_2,    m,patch);
    }
#endif   
/*_______________________________________________________________________*/

}


/* ---------------------------------------------------------------------
                                                        S  T  E  P     3 
GENERAL INFORMATION
 Function:  ICE::actuallyStep3--
 Filename:  ICE_actual.cc
 Purpose:   - Compute the face centered pressure
             
History: 
Version   Programmer         Date       Description                      
-------   ----------         ----       -----------                 
  1.0     John Schmidt   06/23/00                              
_____________________________________________________________________*/
void ICE::actuallyStep3(const ProcessorGroup*,
			const Patch* patch,
			DataWarehouseP& from_dw,
			DataWarehouseP& to_dw)
{
    int m;
/*______________________________________________________________________
*    F  R  O  M     D  W     W  R  A  P  P  E  R 
*--------------------------------------*/
#if switch_UCF_step3OnOff  
    for (m = 1; m<= nMaterials; m++)
    {
        CCVariable<double>  press_cc_2;
        CCVariable<double>  rho_cc_2;
        CCVariable<double>  temp_cc_2;
        CCVariable<Vector>  vel_cc_2;
 
        from_dw->get(       press_cc_2,     press_CCLabel_2,m, patch, Ghost::None, 0);
        from_dw->get(       rho_cc_2,       rho_CCLabel_2,  m, patch, Ghost::None, 0);
        from_dw->get(       temp_cc_2,      temp_CCLabel_2, m, patch, Ghost::None, 0);
        from_dw->get(       vel_cc_2,       vel_CCLabel_2,  m, patch, Ghost::None, 0);
        /*__________________________________
        *   UCF NR
        *___________________________________*/
        ICE::before_each_step_wrapper(  patch,
                            m,
                            press_cc_2,     press_CC,
                            rho_cc_2,       rho_CC,
                            temp_cc_2,      Temp_CC,
                            vel_cc_2,
                            uvel_CC,        vvel_CC,        wvel_CC);
    }
#endif    
/*______________________________________________________________________*/


#if switch_step3_OnOff                                  
    cerr << "Actually doing step 3" << endl;
    press_face(         
                    xLoLimit,       yLoLimit,       zLoLimit,
                    xHiLimit,       yHiLimit,       zHiLimit,
                    delX,           delY,           delZ,
                    BC_types,       BC_float_or_fixed, BC_Values,
                    press_CC,       press_FC,       rho_CC, 
                    nMaterials );
#endif

/*______________________________________________________________________
*     T  O     D  W     W  R  A  P  P  E  R  
*--------------------------------------*/
#if switch_UCF_step3OnOff
    for (m = 1; m<= nMaterials; m++)
    {
        CCVariable<double>  press_cc_3;
        CCVariable<double>  rho_cc_3;
        CCVariable<double>  temp_cc_3;
        CCVariable<Vector>  vel_cc_3;
 
        to_dw->allocate(    press_cc_3, press_CCLabel_3,  m,patch);
        to_dw->allocate(    rho_cc_3,   rho_CCLabel_3,    m,patch);
        to_dw->allocate(    temp_cc_3,  temp_CCLabel_3,   m,patch);
        to_dw->allocate(    vel_cc_3,   vel_CCLabel_3,    m,patch);

        /*__________________________________
        *   UCF NR
        *___________________________________*/
        ICE::after_each_step_wrapper(   patch,
                            m,
                            press_cc_3,  press_CC,
                            rho_cc_3,    rho_CC,
                            temp_cc_3,   Temp_CC,
                            vel_cc_3,    uvel_CC,     vvel_CC,        wvel_CC);
 
        to_dw->put(         press_cc_3,  press_CCLabel_3,  m,patch);
        to_dw->put(         rho_cc_3,    rho_CCLabel_3,    m,patch);
        to_dw->put(         temp_cc_3,   temp_CCLabel_3,   m,patch);
        to_dw->put(         vel_cc_3,    vel_CCLabel_3,    m,patch);
    }
#endif   
/*_______________________________________________________________________*/
}





/* ---------------------------------------------------------------------
                                                        S  T  E  P     4 
GENERAL INFORMATION
 Function:  ICE::actuallyStep4--
 Filename:  ICE_actual.cc
 Purpose:   Compute sources of mass, momentum and energy
            Sources due to mass conversion, gravity
            pressure, divergence of the stressv and momentum exchange
            
             
History: 
Version   Programmer         Date       Description                      
-------   ----------         ----       -----------                 
  1.0     John Schmidt   06/23/00                              
_____________________________________________________________________*/
void ICE::actuallyStep4(const ProcessorGroup*,
			const Patch* patch,
			DataWarehouseP& from_dw,
			DataWarehouseP& to_dw)
{
    int m;
/*______________________________________________________________________
*    F  R  O  M     D  W     W  R  A  P  P  E  R 
*--------------------------------------*/
#if switch_UCF_step4OnOff  
    for (m = 1; m<= nMaterials; m++)
    {
        CCVariable<double>  press_cc_3;
        CCVariable<double>  rho_cc_3;
        CCVariable<double>  temp_cc_3;
        CCVariable<Vector>  vel_cc_3;
 
        from_dw->get(       press_cc_3,     press_CCLabel_3,m, patch, Ghost::None, 0);
        from_dw->get(       rho_cc_3,       rho_CCLabel_3,  m, patch, Ghost::None, 0);
        from_dw->get(       temp_cc_3,      temp_CCLabel_3, m, patch, Ghost::None, 0);
        from_dw->get(       vel_cc_3,       vel_CCLabel_3,  m, patch, Ghost::None, 0);
        /*__________________________________
        *   UCF NR
        *___________________________________*/
        ICE::before_each_step_wrapper(  patch,
                            m,
                            press_cc_3,     press_CC,
                            rho_cc_3,       rho_CC,
                            temp_cc_3,      Temp_CC,
                            vel_cc_3,
                            uvel_CC,        vvel_CC,        wvel_CC);
    }
#endif    
/*______________________________________________________________________*/
 
#if (switch_step4_OnOff == 1 && switch_Compute_burgers_eq == 0) 
    cerr << "Actually doing step 4" << endl;
    accumulate_momentum_source_sinks(
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        delt,
                        delX,           delY,           delZ,
                        grav,
                        mass_CC,        rho_CC,         press_FC,
                        Temp_CC,        cv_CC,
                        uvel_CC,        vvel_CC,        wvel_CC,
                        tau_X_FC,       tau_Y_FC,       tau_Z_FC,
                        viscosity_CC,
                        xmom_source,    ymom_source,    zmom_source,
                        nMaterials   );

     accumulate_energy_source_sinks(
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        delt,
                        delX,           delY,           delZ,
                        grav,           mass_CC,        rho_CC,
                        press_CC,       delPress_CC,    Temp_CC,
                        cv_CC,          speedSound,
                        uvel_CC,        vvel_CC,        wvel_CC,
                        div_velFC_CC,
                        int_eng_source,
                        nMaterials   );

    #endif
    
/*______________________________________________________________________
*     T  O     D  W     W  R  A  P  P  E  R  
*--------------------------------------*/
#if switch_UCF_step4OnOff
    for (m = 1; m<= nMaterials; m++)
    {
        CCVariable<double>  press_cc_4;
        CCVariable<double>  rho_cc_4;
        CCVariable<double>  temp_cc_4;
        CCVariable<Vector>  vel_cc_4;
 
        to_dw->allocate(    press_cc_4, press_CCLabel_4,  m,patch);
        to_dw->allocate(    rho_cc_4,   rho_CCLabel_4,    m,patch);
        to_dw->allocate(    temp_cc_4,  temp_CCLabel_4,   m,patch);
        to_dw->allocate(    vel_cc_4,   vel_CCLabel_4,    m,patch);

        /*__________________________________
        *   UCF NR
        *___________________________________*/
        ICE::after_each_step_wrapper(   patch,
                            m,
                            press_cc_4,  press_CC,
                            rho_cc_4,    rho_CC,
                            temp_cc_4,   Temp_CC,
                            vel_cc_4,    uvel_CC,     vvel_CC,        wvel_CC);
 
        to_dw->put(         press_cc_4,  press_CCLabel_4,  m,patch);
        to_dw->put(         rho_cc_4,    rho_CCLabel_4,    m,patch);
        to_dw->put(         temp_cc_4,   temp_CCLabel_4,   m,patch);
        to_dw->put(         vel_cc_4,    vel_CCLabel_4,    m,patch);
    }
#endif   
/*_______________________________________________________________________*/

}

/* ---------------------------------------------------------------------
                                                        S  T  E  P     5 
GENERAL INFORMATION
 Function:  ICE::actuallyStep5--
 Filename:  ICE_actual.cc
 Computes:  Lagrangian volume (currently not used in algoritm)
            - Converts primative variables into flux form
            - computes lagrangian mass, momentum and energy   
            
History: 
Version   Programmer         Date       Description                      
-------   ----------         ----       -----------                 
  1.0     John Schmidt   06/23/00                              
_____________________________________________________________________*/
void ICE::actuallyStep5(const ProcessorGroup*,
			const Patch* patch,
			DataWarehouseP& from_dw,
			DataWarehouseP& to_dw)
{
    int m;
/*______________________________________________________________________
*    F  R  O  M     D  W     W  R  A  P  P  E  R 
*--------------------------------------*/
#if switch_UCF_step5OnOff  
    for (m = 1; m<= nMaterials; m++)
    {
        CCVariable<double>  press_cc_4;
        CCVariable<double>  rho_cc_4;
        CCVariable<double>  temp_cc_4;
        CCVariable<Vector>  vel_cc_4;
 
        from_dw->get(       press_cc_4,     press_CCLabel_4,m, patch, Ghost::None, 0);
        from_dw->get(       rho_cc_4,       rho_CCLabel_4,  m, patch, Ghost::None, 0);
        from_dw->get(       temp_cc_4,      temp_CCLabel_4, m, patch, Ghost::None, 0);
        from_dw->get(       vel_cc_4,       vel_CCLabel_4,  m, patch, Ghost::None, 0);
        /*__________________________________
        *   UCF NR
        *___________________________________*/
        ICE::before_each_step_wrapper(  patch,
                            m,
                            press_cc_4,     press_CC,
                            rho_cc_4,       rho_CC,
                            temp_cc_4,      Temp_CC,
                            vel_cc_4,
                            uvel_CC,        vvel_CC,        wvel_CC);
    }
#endif    
/*______________________________________________________________________*/
     /*__________________________________
    *    S  T  E  P     5                        
    *   Compute Lagrangian values for the volume 
    *   mass, momentum and energy.
    *   Lagrangian values are the sum of the time n
    *   values and the sources computed in 4
    *___________________________________*/
 cerr << "Actually doing step 5" << endl;
#if switch_step5_OnOff 
    lagrangian_vol(     xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        delX,           delY,           delZ,
                        delt,           
                        Vol_L_CC,       Vol_CC,
                        uvel_FC,        vvel_FC,        wvel_FC,
                        nMaterials);
                        
    calc_flux_or_primitive_vars(    -1,           
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        rho_CC,         Vol_CC,         
                        uvel_CC,        vvel_CC,        wvel_CC,        
                        xmom_CC,        ymom_CC,        zmom_CC,
                        cv_CC,          int_eng_CC,     Temp_CC,
                        nMaterials );                       
                        
    lagrangian_values(  
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        Vol_L_CC,       Vol_CC,         rho_CC,
                        rho_L_CC,
                        xmom_CC,        ymom_CC,        zmom_CC,
                        uvel_CC,        vvel_CC,        wvel_CC,
                        xmom_L_CC,      ymom_L_CC,      zmom_L_CC,
                        mass_L_CC,      mass_source,    
                        xmom_source,    ymom_source,    zmom_source,
                        int_eng_CC,     int_eng_L_CC,   int_eng_source,
                        nMaterials);
    #endif  
/*______________________________________________________________________
*     T  O     D  W     W  R  A  P  P  E  R  
*--------------------------------------*/
#if switch_UCF_step5OnOff
    for (m = 1; m<= nMaterials; m++)
    {
        CCVariable<double>  press_cc_5;
        CCVariable<double>  rho_cc_5;
        CCVariable<double>  temp_cc_5;
        CCVariable<Vector>  vel_cc_5;
 
        to_dw->allocate(    press_cc_5, press_CCLabel_5,  m,patch);
        to_dw->allocate(    rho_cc_5,   rho_CCLabel_5,    m,patch);
        to_dw->allocate(    temp_cc_5,  temp_CCLabel_5,   m,patch);
        to_dw->allocate(    vel_cc_5,   vel_CCLabel_5,    m,patch);

        /*__________________________________
        *   UCF NR
        *___________________________________*/
        ICE::after_each_step_wrapper(   patch,
                            m,
                            press_cc_5,  press_CC,
                            rho_cc_5,    rho_CC,
                            temp_cc_5,   Temp_CC,
                            vel_cc_5,    uvel_CC,     vvel_CC,        wvel_CC);
 
        to_dw->put(         press_cc_5,  press_CCLabel_5,  m,patch);
        to_dw->put(         rho_cc_5,    rho_CCLabel_5,    m,patch);
        to_dw->put(         temp_cc_5,   temp_CCLabel_5,   m,patch);
        to_dw->put(         vel_cc_5,    vel_CCLabel_5,    m,patch);
    }
#endif   
/*_______________________________________________________________________*/

}





/* ---------------------------------------------------------------------
                                            S  T  E  P     6     &     7  
GENERAL INFORMATION
 Function:  ICE::actuallyStep6--
 Filename:  ICE_actual.cc
 Computes:  Advects the mass, momentum and energy 
            Computes time advances mass, momentum and energy
            coverts flux mass, momentum and energy into the 
            primative varibles   
            
History: 
Version   Programmer         Date       Description                      
-------   ----------         ----       -----------                 
  1.0     John Schmidt   06/23/00                              
_____________________________________________________________________*/
void ICE::actuallyStep6and7(const ProcessorGroup*,
			const Patch* patch,
			DataWarehouseP& from_dw,
			DataWarehouseP& to_dw)
{
    int m;
/*______________________________________________________________________
*    F  R  O  M     D  W     W  R  A  P  P  E  R 
*--------------------------------------*/
#if switch_UCF_step6_7OnOff  
    for (m = 1; m<= nMaterials; m++)
    {
        CCVariable<double>  press_cc_5;
        CCVariable<double>  rho_cc_5;
        CCVariable<double>  temp_cc_5;
        CCVariable<Vector>  vel_cc_5;
 
        from_dw->get(       press_cc_5,     press_CCLabel_5,m, patch, Ghost::None, 0);
        from_dw->get(       rho_cc_5,       rho_CCLabel_5,  m, patch, Ghost::None, 0);
        from_dw->get(       temp_cc_5,      temp_CCLabel_5, m, patch, Ghost::None, 0);
        from_dw->get(       vel_cc_5,       vel_CCLabel_5,  m, patch, Ghost::None, 0);
        /*__________________________________
        *   UCF NR
        *___________________________________*/
        ICE::before_each_step_wrapper(  patch,
                            m,
                            press_cc_5,     press_CC,
                            rho_cc_5,       rho_CC,
                            temp_cc_5,      Temp_CC,
                            vel_cc_5,
                            uvel_CC,        vvel_CC,        wvel_CC);
    }
#endif    
/*______________________________________________________________________*/
    /*_________________________________   
    *    S  T  E  P     6                            
    *   Compute the advection of mass,
    *   momentum and energy.  These
    *   quantities are advected using the face
    *   centered velocities velocities from 2
    *                  
    *    S  T  E  P     7 
    *   Compute the time advanced values for
    *   mass, momentum and energy.  "Time advanced"
    *   means the sum of the "Lagrangian" values,
    *   found in 5 and the advection contribution
    *   from 6                      
    *______________________________ */  
    #if (switch_step7_OnOff== 1 || switch_step6_OnOff == 1)
    cerr << "Actually doing step 6" << endl;
     advect_and_advance_in_time(   
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        delX,           delY,           delZ,
                        Vol_CC,         rho_CC,
                        xmom_CC,        ymom_CC,        zmom_CC,
                        Vol_L_CC,       rho_L_CC,       mass_L_CC,
                        xmom_L_CC,      ymom_L_CC,      zmom_L_CC,
                        int_eng_CC,     int_eng_L_CC,
                        uvel_FC,        vvel_FC,        wvel_FC,
                        delt,           nMaterials);
    /*__________________________________
    *   Backout the velocities from the 
    *   the momentum
    *___________________________________*/                        
    calc_flux_or_primitive_vars(    1,           
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        rho_CC,         Vol_CC,         
                        uvel_CC,        vvel_CC,        wvel_CC,        
                        xmom_CC,        ymom_CC,        zmom_CC,
                        cv_CC,          int_eng_CC,     Temp_CC,
                        nMaterials ); 
    #endif
/*______________________________________________________________________
*     T  O     D  W     W  R  A  P  P  E  R  
*--------------------------------------*/
#if switch_UCF_step6_7OnOff
    for (m = 1; m<= nMaterials; m++)
    {
        CCVariable<double>  press_cc_6_7;
        CCVariable<double>  rho_cc_6_7;
        CCVariable<double>  temp_cc_6_7;
        CCVariable<Vector>  vel_cc_6_7;
 
        to_dw->allocate(    press_cc_6_7, press_CCLabel_6_7,  m,patch);
        to_dw->allocate(    rho_cc_6_7,   rho_CCLabel_6_7,    m,patch);
        to_dw->allocate(    temp_cc_6_7,  temp_CCLabel_6_7,   m,patch);
        to_dw->allocate(    vel_cc_6_7,   vel_CCLabel_6_7,    m,patch);

        /*__________________________________
        *   UCF NR
        *___________________________________*/
        ICE::after_each_step_wrapper(   patch,
                            m,
                            press_cc_6_7,  press_CC,
                            rho_cc_6_7,    rho_CC,
                            temp_cc_6_7,   Temp_CC,
                            vel_cc_6_7,    uvel_CC,     vvel_CC,        wvel_CC);
 
        to_dw->put(         press_cc_6_7,  press_CCLabel_6_7,  m,patch);
        to_dw->put(         rho_cc_6_7,    rho_CCLabel_6_7,    m,patch);
        to_dw->put(         temp_cc_6_7,   temp_CCLabel_6_7,   m,patch);
        to_dw->put(         vel_cc_6_7,    vel_CCLabel_6_7,    m,patch);
    }
#endif   
/*_______________________________________________________________________*/
}



/* ---------------------------------------------------------------------
                                              
GENERAL INFORMATION
 Function:  ICE::actually_Bottom_of_main_loop--
 Filename:  ICE_actual.cc
 Purpose:   - output tecplot files
            - PGPLOT output
            - increment time
History: 
Version   Programmer         Date       Description                      
-------   ----------         ----       -----------                 
  1.0     John Schmidt   06/23/00                              
_____________________________________________________________________*/
void ICE::actually_Bottom_of_main_loop(const ProcessorGroup*,
                    const Patch* patch,
		      DataWarehouseP& from_dw,
		      DataWarehouseP& to_dw)
{
    double t      = this->cheat_t;
    double delt   = this->cheat_delt;
    int   should_I_write_output;
    int m;  

  /*__________________________________
   *   Plotting variables
   *___________________________________*/
#if (switchDebug_main == 1|| switchDebug_main == 2 || switchDebug_main_input == 1)
    #include "plot_declare_vars.h"   
#endif
    cerr << "Actually doing at bottom of main loop" << endl;
/*______________________________________________________________________
*    F  R  O  M     D  W     W  R  A  P  P  E  R 
*--------------------------------------*/
#if switch_UCF_Bottom_of_main_loopOnOff 
    for (m = 1; m<= nMaterials; m++)
    { 
        CCVariable<double>  press_cc_6_7;
        CCVariable<double>  rho_cc_6_7;
        CCVariable<double>  temp_cc_6_7;
        CCVariable<Vector>  vel_cc_6_7;

        from_dw->get(       press_cc_6_7,   press_CCLabel_6_7,m, patch, Ghost::None, 0);
        from_dw->get(       rho_cc_6_7,     rho_CCLabel_6_7,  m, patch, Ghost::None, 0);
        from_dw->get(       temp_cc_6_7,    temp_CCLabel_6_7, m, patch, Ghost::None, 0);
        from_dw->get(       vel_cc_6_7,     vel_CCLabel_6_7,  m, patch, Ghost::None, 0);
        /*__________________________________
        *   UCF NR
        *___________________________________*/
        ICE::before_each_step_wrapper(  patch,
                            m,
                            press_cc_6_7,   press_CC,
                            rho_cc_6_7,     rho_CC,
                            temp_cc_6_7,    Temp_CC,
                            vel_cc_6_7,
                            uvel_CC,        vvel_CC,        wvel_CC);

    }
#endif
/*______________________________________________________________________*/

    should_I_write_output = Is_it_time_to_write_output( t, t_output_vars  );
    /*__________________________________
    *    T  E  C  P  L  O  T  
    *___________________________________*/     
    #if tecplot    
    if ( should_I_write_output == YES)
    {                     
        tecplot_CC(         
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        x_CC,           y_CC,           z_CC,
                        uvel_CC,        vvel_CC,        wvel_CC,
                        press_CC,       Temp_CC,        rho_CC,
                        scalar1_CC,     scalar2_CC,     scalar3_CC,
                        fileNum,        output_file_basename,       output_file_desc,
                        nMaterials);

        tecplot_FC(         
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        x_FC,           y_FC,           z_FC,
                        uvel_FC,        vvel_FC,        wvel_FC,
                        fileNum,        output_file_basename,       output_file_desc,
                        nMaterials );
                            
        fileNum ++;
    } 
    #endif 
    /*__________________________________
     *  Plotting section
     *___________________________________*/
    #if switchDebug_main
    if ( should_I_write_output == YES)
    {
      #define switchInclude_main 1
      #include "debugcode.i"
      #undef switchInclude_main 
    }
    #endif
    /*__________________________________
     *  Clean up the plotting windows 
     *___________________________________*/
     putenv("PGPLOT_I_AM_HERE=1");              
     /* tell the plotting routine that   */
     /* you're at the bottom of main     */
     putenv("PGPLOT_OPEN_NEW_WINDOWS=1"); 
    
    
    /*__________________________________
     *   Advance time
     *___________________________________*/
    t = t + delt;
    fprintf(stderr,"\nTime is %f, timestep is %f\n",t,delt);  
    
    fprintf(stderr, "press return to continue \n");
    getchar();
 




/*______________________________________________________________________
*     T  O     D  W     W  R  A  P  P  E  R  
*--------------------------------------*/
#if switch_UCF_Bottom_of_main_loopOnOff
   for (m = 1; m<= nMaterials; m++)
   {
       CCVariable<Vector>  new_vel_cc;
       CCVariable<double>  new_temp_cc;
       CCVariable<double>  new_rho_cc;
       CCVariable<double>  new_press_cc;
 
       to_dw->allocate(    new_vel_cc,     vel_CCLabel,    m,patch);
       to_dw->allocate(    new_temp_cc,    temp_CCLabel,   m,patch);
       to_dw->allocate(    new_rho_cc,     rho_CCLabel,    m,patch);
       to_dw->allocate(    new_press_cc,   press_CCLabel,  m,patch);
       /*__________________________________
       *   UCF NR
       *___________________________________*/
       ICE::after_each_step_wrapper(  patch,
                           m,
                           new_press_cc,   press_CC,
                           new_rho_cc,     rho_CC,
                           new_temp_cc,    Temp_CC,
                           new_vel_cc,
                           uvel_CC,         vvel_CC,        wvel_CC);
 
       to_dw->put(         new_press_cc,   press_CCLabel,  m,patch);
       to_dw->put(         new_rho_cc,     rho_CCLabel,    m,patch);
       to_dw->put(         new_temp_cc,    temp_CCLabel,   m,patch);
       to_dw->put(         new_vel_cc,     vel_CCLabel,    m,patch);
   }
#endif  
/*______________________________________________________________________*/

}










