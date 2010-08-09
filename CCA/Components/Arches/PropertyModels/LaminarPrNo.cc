#include <CCA/Components/Arches/PropertyModels/LaminarPrNo.h>

using namespace Uintah; 

//---------------------------------------------------------------------------
//Method: Constructor
//---------------------------------------------------------------------------
LaminarPrNo::LaminarPrNo( std::string prop_name, SimulationStateP& shared_state ) : PropertyModelBase( prop_name, shared_state )
{
  std::string varlabel_name = "laminar_pr"; 
  // the prop is the pr number.  Along with this, we will also give access to D and mu
  _prop_label = VarLabel::create( varlabel_name, CCVariable<double>::getTypeDescription() ); 

  // additional local labels as needed by this class (delete this if it isn't used): 
  std::string name = "laminar_viscosity";
  _mu_label = VarLabel::create( name, CCVariable<double>::getTypeDescription() ); // Note: you need to add the label to the .h file
  _extra_local_labels.push_back( _mu_label ); 

  name = "laminar_diffusion_coef";
  _D_label = VarLabel::create( name, CCVariable<double>::getTypeDescription() ); // Note: you need to add the label to the .h file
  _extra_local_labels.push_back( _D_label ); 

}

//---------------------------------------------------------------------------
//Method: Destructor
//---------------------------------------------------------------------------
LaminarPrNo::~LaminarPrNo( )
{
  // Destroying all local VarLabels stored in _extra_local_labels: 
  for (std::vector<const VarLabel*>::iterator iter = _extra_local_labels.begin(); iter != _extra_local_labels.end(); iter++){

    VarLabel::destroy( *iter ); 

  }
  // Clean up anything else here ... 
}


//---------------------------------------------------------------------------
//Method: Problem Setup
//---------------------------------------------------------------------------
void LaminarPrNo::problemSetup( const ProblemSpecP& inputdb )
{
  ProblemSpecP db = inputdb; 

  _binary_mixture = true; //set this to true for now.

  db->require( "mix_frac_label", _mix_frac_label_name ); 
  db->require( "atm_pressure",   _pressure ); 

  //Fuel
  db->findBlock("fuel")->require( "molar_mass", _molar_mass_a ); 
  db->findBlock("fuel")->require( "critical_temperature", _crit_pressure_a ); 
  db->findBlock("fuel")->require( "dipole_moment", _dipole_moment_a ); 
  db->findBlock("fuel")->require( "lennard_jones_length", _lj_sigma_a ); 
  db->findBlock("fuel")->require( "lennard_jones_energy", _lj_ek_a ); 
  db->findBlock("fuel")->require( "viscosity", _viscosity_a ); 

  //Oxidizer
  db->findBlock("oxidizer")->require( "molar_mass", _molar_mass_b ); 
  db->findBlock("oxidizer")->require( "critical_temperature", _crit_pressure_b ); 
  db->findBlock("oxidizer")->require( "dipole_moment", _dipole_moment_b ); 
  db->findBlock("oxidizer")->require( "lennard_jones_length", _lj_sigma_b ); 
  db->findBlock("oxidizer")->require( "lennard_jones_energy", _lj_ek_b ); 
  db->findBlock("oxidizer")->require( "viscosity", _viscosity_b ); 

}

//---------------------------------------------------------------------------
//Method: Schedule Compute Property
//---------------------------------------------------------------------------
void LaminarPrNo::sched_computeProp( const LevelP& level, SchedulerP& sched, int time_substep )
{
  std::string taskname = "LaminarPrNo::computeProp"; 
  Task* tsk = scinew Task( taskname, this, &LaminarPrNo::computeProp, time_substep ); 

  if ( !(_has_been_computed) ) {

    if ( time_substep == 0 ) {
      
      tsk->computes( _prop_label ); 

      for (std::vector<const VarLabel*>::iterator iter = _extra_local_labels.begin(); iter != _extra_local_labels.end(); iter++){
        tsk->computes( *iter ); 
      }

    } else {

      tsk->modifies( _prop_label ); 

      for (std::vector<const VarLabel*>::iterator iter = _extra_local_labels.begin(); iter != _extra_local_labels.end(); iter++){
        tsk->modifies( *iter ); 
      }

    }

    const VarLabel* den_label = VarLabel::find("densityCP"); 
    tsk->requires( Task::NewDW, den_label,     Ghost::None, 0 );  
    //const VarLabel* T_label   = VarLabel::find("tempIN"); 
    //tsk->requires( Task::NewDW, T_label,     Ghost::None, 0 );  
    const VarLabel* f_label   = VarLabel::find("scalarSP"); 
    tsk->requires( Task::NewDW, f_label,     Ghost::None, 0 );  

    sched->addTask( tsk, level->eachPatch(), _shared_state->allArchesMaterials() ); 
    
    _has_been_computed = true; 

  }
}

//---------------------------------------------------------------------------
//Method: Actually Compute Property
//---------------------------------------------------------------------------
void LaminarPrNo::computeProp(const ProcessorGroup* pc, 
                                    const PatchSubset* patches, 
                                    const MaterialSubset* matls, 
                                    DataWarehouse* old_dw, 
                                    DataWarehouse* new_dw, 
                                    int time_substep )
{
  //patch loop
  for (int p=0; p < patches->size(); p++){

    const Patch* patch = patches->get(p);
    int archIndex = 0;
    int matlIndex = _shared_state->getArchesMaterial(archIndex)->getDWIndex(); 

    CCVariable<double> Pr;  // Prandlt number 
    CCVariable<double> D;   // Diffusion coef
    CCVariable<double> mu;  // viscosity
    constCCVariable<double> f;   // mixture fraction 
    constCCVariable<double> T;   // temperature 
    constCCVariable<double> rho; // density

    if ( new_dw->exists( _prop_label, matlIndex, patch ) ){
      new_dw->getModifiable( Pr, _prop_label, matlIndex, patch ); 
      new_dw->getModifiable( D,    _D_label, matlIndex, patch ); 
      new_dw->getModifiable( mu,   _mu_label, matlIndex, patch ); 
    } else {
      new_dw->allocateAndPut( Pr, _prop_label, matlIndex, patch ); 
      new_dw->allocateAndPut( D,    _D_label, matlIndex, patch ); 
      new_dw->allocateAndPut( mu,   _mu_label, matlIndex, patch ); 
      Pr.initialize(0.0); 
      D.initialize(0.0); 
      mu.initialize(0.0); 
    }

    const VarLabel* f_label = VarLabel::find( "scalarSP" ); 
    new_dw->get( f, f_label, matlIndex, patch, Ghost::None, 0 ); 
    //const VarLabel* T_label = VarLabel::find( "tempIN" ); 
    //new_dw->get( T, T_label, matlIndex, patch, Ghost::None, 0 ); 
    const VarLabel* den_label = VarLabel::find( "densityCP" ); 
    new_dw->get( rho, den_label, matlIndex, patch, Ghost::None, 0 ); 

    CellIterator iter = patch->getCellIterator(); 

    for (iter.begin(); !iter.done(); iter++){

      IntVector c = *iter; //i,j,k location

      double T = 298.0; 

      // viscosity 
      mu[c] = getVisc( f[c], T );
      // diffusion coefficient 
      D[c]  = getDiffCoef( T );
      // prandlt number
      Pr[c] = mu[c] / ( rho[c] * D[c] ); 

    }
  }
}

//---------------------------------------------------------------------------
//Method: Scheduler for Dummy Initialization
//---------------------------------------------------------------------------
void LaminarPrNo::sched_dummyInit( const LevelP& level, SchedulerP& sched )
{

  std::string taskname = "LaminarPrNo::dummyInit"; 

  Task* tsk = scinew Task(taskname, this, &LaminarPrNo::dummyInit);
  tsk->computes(_prop_label); 
  tsk->requires( Task::OldDW, _prop_label, Ghost::None, 0 ); 

  sched->addTask(tsk, level->eachPatch(), _shared_state->allArchesMaterials());

}

//---------------------------------------------------------------------------
//Method: Actually do the Dummy Initialization
//---------------------------------------------------------------------------
void LaminarPrNo::dummyInit( const ProcessorGroup* pc, 
                                            const PatchSubset* patches, 
                                            const MaterialSubset* matls, 
                                            DataWarehouse* old_dw, 
                                            DataWarehouse* new_dw )
{
  //patch loop
  for (int p=0; p < patches->size(); p++){

    const Patch* patch = patches->get(p);
    int archIndex = 0;
    int matlIndex = _shared_state->getArchesMaterial(archIndex)->getDWIndex(); 

    CCVariable<double> prop; 
    constCCVariable<double> old_prop; 

    new_dw->allocateAndPut( prop, _prop_label, matlIndex, patch ); 
    old_dw->get( old_prop, _prop_label, matlIndex, patch, Ghost::None, 0); 

    //prop.initialize(0.0); <--- Careful, don't reinitialize if you don't want to 
    prop.copyData( old_prop );

  }
}

//---------------------------------------------------------------------------
//Method: Scheduler for Initializing the Property
//---------------------------------------------------------------------------
void LaminarPrNo::sched_initialize( const LevelP& level, SchedulerP& sched )
{
  std::string taskname = "LaminarPrNo::initialize"; 

  Task* tsk = scinew Task(taskname, this, &LaminarPrNo::initialize);
  tsk->computes(_prop_label); 

  sched->addTask(tsk, level->eachPatch(), _shared_state->allArchesMaterials());
}

//---------------------------------------------------------------------------
//Method: Actually Initialize the Property
//---------------------------------------------------------------------------
void LaminarPrNo::initialize( const ProcessorGroup* pc, 
                                    const PatchSubset* patches, 
                                    const MaterialSubset* matls, 
                                    DataWarehouse* old_dw, 
                                    DataWarehouse* new_dw )
{
  //patch loop
  for (int p=0; p < patches->size(); p++){

    const Patch* patch = patches->get(p);
    int archIndex = 0;
    int matlIndex = _shared_state->getArchesMaterial(archIndex)->getDWIndex(); 

    CCVariable<double> prop; 

    new_dw->allocateAndPut( prop, _prop_label, matlIndex, patch ); 
    prop.initialize(0.0); 

    PropertyModelBase::base_initialize( patch, prop ); // generic initialization functionality 

  }
}
