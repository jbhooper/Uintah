#ifndef Uintah_Component_Arches_KScalarRHS_h
#define Uintah_Component_Arches_KScalarRHS_h

/*
 * The MIT License
 *
 * Copyright (c) 1997-2018 The University of Utah
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

#include <CCA/Components/Arches/Task/TaskInterface.h>
#include <CCA/Components/Arches/Transport/TransportHelper.h>
#include <CCA/Components/Arches/GridTools.h>
#include <CCA/Components/Arches/ConvectionHelper.h>
#include <CCA/Components/Arches/Directives.h>
#include <CCA/Components/Arches/BoundaryConditions/BoundaryFunctors.h>
#include <CCA/Components/Arches/UPSHelper.h>

namespace Uintah{

  template<typename T, typename PT>
  class KScalarRHS : public TaskInterface {

public:

    KScalarRHS<T, PT>( std::string task_name, int matl_index );
    ~KScalarRHS<T, PT>();

    typedef std::vector<ArchesFieldContainer::VariableInformation> ArchesVIVector;

    void problemSetup( ProblemSpecP& db );

    void register_initialize( ArchesVIVector& variable_registry , const bool pack_tasks);

    void register_timestep_init( ArchesVIVector& variable_registry , const bool packed_tasks);

    void register_timestep_eval( ArchesVIVector& variable_registry,
                                 const int time_substep , const bool packed_tasks);

    void register_compute_bcs( ArchesVIVector& variable_registry,
                               const int time_substep , const bool packed_tasks);

    void compute_bcs( const Patch* patch, ArchesTaskInfoManager* tsk_info );

    void initialize( const Patch* patch, ArchesTaskInfoManager* tsk_info );

    void timestep_init( const Patch* patch, ArchesTaskInfoManager* tsk_info );

    void eval( const Patch* patch, ArchesTaskInfoManager* tsk_info );

    void create_local_labels();

    //Build instructions for this (KScalarRHS) class.
    class Builder : public TaskInterface::TaskBuilder {

      public:

      Builder( std::string task_name, int matl_index )
      : m_task_name(task_name), m_matl_index(matl_index){}
      ~Builder(){}

      KScalarRHS* build()
      { return scinew KScalarRHS<T, PT>( m_task_name, m_matl_index ); }

      private:

      std::string m_task_name;
      int m_matl_index;

    };

private:

    typedef typename ArchesCore::VariableHelper<T>::ConstType CT;
    typedef typename ArchesCore::VariableHelper<T>::XFaceType FXT;
    typedef typename ArchesCore::VariableHelper<T>::YFaceType FYT;
    typedef typename ArchesCore::VariableHelper<T>::ZFaceType FZT;
    typedef typename ArchesCore::VariableHelper<CT>::XFaceType CFXT;
    typedef typename ArchesCore::VariableHelper<CT>::YFaceType CFYT;
    typedef typename ArchesCore::VariableHelper<CT>::ZFaceType CFZT;

    typedef typename ArchesCore::VariableHelper<PT>::XFaceType FluxXT;
    typedef typename ArchesCore::VariableHelper<PT>::YFaceType FluxYT;
    typedef typename ArchesCore::VariableHelper<PT>::ZFaceType FluxZT;

    //typedef typename ArchesCore::VariableHelper<T>::ConstXFaceType CFXT;
    //typedef typename ArchesCore::VariableHelper<T>::ConstYFaceType CFYT;
    //typedef typename ArchesCore::VariableHelper<T>::ConstZFaceType CFZT;

    std::string m_D_name;
    std::string m_premultiplier_name;
    std::string m_x_velocity_name;
    std::string m_y_velocity_name;
    std::string m_z_velocity_name;
    std::string m_eps_name;

    ArchesCore::EQUATION_CLASS m_eqn_class;

    std::vector<std::string> m_eqn_names;
    std::vector<bool> m_do_diff;
    std::vector<bool> m_do_clip;
    std::vector<double> m_low_clip;
    std::vector<double> m_high_clip;
    std::vector<double> m_init_value;

    bool m_has_D;

    int m_total_eqns;

    struct SourceInfo{
      std::string name;
      double weight;
    };

    std::vector<std::vector<SourceInfo> > m_source_info;
    std::map<std::string, double> m_scaling_info;

    std::vector<LIMITER> m_conv_scheme;

    ArchesCore::BCFunctors<T>* m_boundary_functors;

  };

  //------------------------------------------------------------------------------------------------
  template <typename T, typename PT>
  KScalarRHS<T, PT>::KScalarRHS( std::string task_name, int matl_index ) :
  TaskInterface( task_name, matl_index ) {

    m_boundary_functors = scinew ArchesCore::BCFunctors<T>();

  }

  //------------------------------------------------------------------------------------------------
  template <typename T, typename PT>
  KScalarRHS<T, PT>::~KScalarRHS(){

    delete m_boundary_functors;

  }

  //------------------------------------------------------------------------------------------------
  template <typename T, typename PT> void
  KScalarRHS<T, PT>::problemSetup( ProblemSpecP& input_db ){

    m_total_eqns = 0;

    ConvectionHelper* conv_helper = scinew ConvectionHelper();

    std::string eqn_class = "density_weighted";
    if ( input_db->findAttribute("class") ){
      input_db->getAttribute("class", eqn_class);
    }

    m_eqn_class = ArchesCore::assign_eqn_class_enum( eqn_class );

    for( ProblemSpecP db = input_db->findBlock("eqn"); db != nullptr; db = db->findNextBlock("eqn") ) {

      //Equation name
      std::string eqn_name;
      db->getAttribute("label", eqn_name);
      m_eqn_names.push_back(eqn_name);

      //Check for something other than density weighted:

      std::string eqn_class_str = "density_weighted";
      if ( db->findBlock("dqmom")){
        eqn_class_str = "dqmom";
      } else if ( db->findBlock("volumetric")){
        eqn_class_str = "volumetric";
      }
      m_eqn_class = ArchesCore::assign_eqn_class_enum( eqn_class_str );

      //Convection
      if ( db->findBlock("convection")){
        std::string conv_scheme;
        db->findBlock("convection")->getAttribute("scheme", conv_scheme);
        m_conv_scheme.push_back(conv_helper->get_limiter_from_string(conv_scheme));
      } else {
        m_conv_scheme.push_back(NOCONV);
      }

      //Diffusion
      m_has_D = false;
      if ( db->findBlock("diffusion")){
        m_do_diff.push_back(true);
        m_has_D = true;
      } else {
        m_do_diff.push_back(false);
      }

      //Clipping
      if ( db->findBlock("clip")){
        m_do_clip.push_back(true);
        double low; double high;
        db->findBlock("clip")->getAttribute("low", low);
        db->findBlock("clip")->getAttribute("high", high);
        m_low_clip.push_back(low);
        m_high_clip.push_back(high);
      } else {
        m_do_clip.push_back(false);
        m_low_clip.push_back(-999.9);
        m_high_clip.push_back(999.9);
      }

      //Scaling Constant
      if ( db->findBlock("scaling")){

        double scaling_constant;
        db->findBlock("scaling")->getAttribute("value", scaling_constant);
        m_scaling_info.insert(std::make_pair(eqn_name, scaling_constant));

      }

      //Initial Value
      if ( db->findBlock("initialize") ){
        double value;
        db->findBlock("initialize")->getAttribute("value",value);
        m_init_value.push_back(value);
      }
      else {
        m_init_value.push_back(0.0);
      }

      std::vector<SourceInfo> eqn_srcs;
      for ( ProblemSpecP src_db = db->findBlock("src"); src_db != nullptr; src_db = src_db->findNextBlock("src") ) {

        std::string src_label;
        double weight = 1.0;

        src_db->getAttribute("label",src_label);

        if ( src_db->findBlock("weight")){
          src_db->findBlock("weight")->getAttribute("value",weight);
        }

        SourceInfo info;
        info.name = src_label;
        info.weight = weight;

        eqn_srcs.push_back(info);

      }

      m_source_info.push_back(eqn_srcs);

    }

    // setup the boundary conditions for this eqn set
    m_boundary_functors->create_bcs( input_db, m_eqn_names );

    delete conv_helper;

    using namespace ArchesCore;

    ArchesCore::GridVarMap<T> var_map;
    var_map.problemSetup( input_db );
    m_eps_name = var_map.vol_frac_name;
    m_x_velocity_name = var_map.uvel_name;
    m_y_velocity_name = var_map.vvel_name;
    m_z_velocity_name = var_map.wvel_name;

    m_premultiplier_name = get_premultiplier_name(m_eqn_class);

    if ( input_db->findBlock("velocity") ){
      // can overide the global velocity space with this:
      input_db->findBlock("velocity")->getAttribute("xlabel",m_x_velocity_name);
      input_db->findBlock("velocity")->getAttribute("ylabel",m_y_velocity_name);
      input_db->findBlock("velocity")->getAttribute("zlabel",m_z_velocity_name);
    }

    // Diffusion coeff -- assuming the same one across all eqns.
    m_D_name = "NA";

    if ( m_has_D ){
      if ( input_db->findBlock("diffusion_coef") ) {
        input_db->findBlock("diffusion_coef")->getAttribute("label",m_D_name);
      } else {
        std::stringstream msg;
        msg << "Error: Diffusion specified for task " << _task_name << std::endl
        << "but no diffusion coefficient label specified." << std::endl;
        throw ProblemSetupException(msg.str(),__FILE__, __LINE__);
      }
    }

  }

  //------------------------------------------------------------------------------------------------
  template <typename T, typename PT>
  void
  KScalarRHS<T, PT>::create_local_labels(){

    const int istart = 0;
    int iend = m_eqn_names.size();
    for (int i = istart; i < iend; i++ ){
      register_new_variable<T>( m_eqn_names[i] );
      if ( m_premultiplier_name != "none" )
        register_new_variable<T>( m_premultiplier_name+m_eqn_names[i] );
      register_new_variable<T>( m_eqn_names[i]+"_rhs" );
      register_new_variable<FXT>( m_eqn_names[i]+"_x_flux" );
      register_new_variable<FYT>( m_eqn_names[i]+"_y_flux" );
      register_new_variable<FZT>( m_eqn_names[i]+"_z_flux" );
    }
  }

  template <typename T, typename PT> void
  KScalarRHS<T, PT>::register_initialize(
    std::vector<ArchesFieldContainer::VariableInformation>& variable_registry,
    const bool packed_tasks ){

    const int istart = 0;
    const int iend = m_eqn_names.size();
    for (int ieqn = istart; ieqn < iend; ieqn++ ){
      register_variable(  m_eqn_names[ieqn], ArchesFieldContainer::COMPUTES , variable_registry );
      if ( m_premultiplier_name != "none" )
        register_variable( m_premultiplier_name+m_eqn_names[ieqn], ArchesFieldContainer::COMPUTES , variable_registry );
      register_variable(  m_eqn_names[ieqn]+"_rhs", ArchesFieldContainer::COMPUTES , variable_registry );
      register_variable(  m_eqn_names[ieqn]+"_x_flux", ArchesFieldContainer::COMPUTES , variable_registry, _task_name );
      register_variable(  m_eqn_names[ieqn]+"_y_flux", ArchesFieldContainer::COMPUTES , variable_registry, _task_name );
      register_variable(  m_eqn_names[ieqn]+"_z_flux", ArchesFieldContainer::COMPUTES , variable_registry, _task_name );
    }

    for ( auto i = m_scaling_info.begin(); i != m_scaling_info.end(); i++ ){
      register_variable( i->first+"_unscaled", ArchesFieldContainer::COMPUTES, variable_registry, _task_name );
    }
  }

  //------------------------------------------------------------------------------------------------
  template <typename T, typename PT> void
  KScalarRHS<T, PT>::initialize( const Patch* patch, ArchesTaskInfoManager* tsk_info ){

    const int istart = 0;
    const int iend = m_eqn_names.size();
    for (int ieqn = istart; ieqn < iend; ieqn++ ){

      double scalar_init_value = m_init_value[ieqn];

      T& rhs = tsk_info->get_uintah_field_add<T>(m_eqn_names[ieqn]+"_rhs");
      T& phi = tsk_info->get_uintah_field_add<T>(m_eqn_names[ieqn]);

      rhs.initialize(0.0);
      phi.initialize(scalar_init_value);

      FXT& x_flux = tsk_info->get_uintah_field_add<FXT>(m_eqn_names[ieqn]+"_x_flux");
      FYT& y_flux = tsk_info->get_uintah_field_add<FYT>(m_eqn_names[ieqn]+"_y_flux");
      FZT& z_flux = tsk_info->get_uintah_field_add<FZT>(m_eqn_names[ieqn]+"_z_flux");

      x_flux.initialize(0.0);
      y_flux.initialize(0.0);
      z_flux.initialize(0.0);

      if ( m_premultiplier_name != "none" ){
        T& rho_phi  = tsk_info->get_uintah_field_add<T>(m_premultiplier_name+m_eqn_names[ieqn]);
        rho_phi.initialize(0.0);
      }

    } //eqn loop

    for ( auto i = m_scaling_info.begin(); i != m_scaling_info.end(); i++ ){
      T& phi_unscaled = tsk_info->get_uintah_field_add<T>( i->first+"_unscaled");
      phi_unscaled.initialize(0.0);
    }

  }

  //------------------------------------------------------------------------------------------------
  template <typename T, typename PT> void
  KScalarRHS<T, PT>::register_timestep_init( std::vector<ArchesFieldContainer::VariableInformation>& variable_registry , const bool packed_tasks){
    const int istart = 0;
    const int iend = m_eqn_names.size();
    for (int ieqn = istart; ieqn < iend; ieqn++ ){
      if ( m_premultiplier_name != "none" ){
        register_variable( m_premultiplier_name+m_eqn_names[ieqn], ArchesFieldContainer::COMPUTES , variable_registry, _task_name  );
        register_variable( m_premultiplier_name+m_eqn_names[ieqn], ArchesFieldContainer::REQUIRES , 0 , ArchesFieldContainer::OLDDW , variable_registry, _task_name );
      }
      register_variable( m_eqn_names[ieqn], ArchesFieldContainer::COMPUTES , variable_registry, _task_name  );
      register_variable( m_eqn_names[ieqn]+"_rhs", ArchesFieldContainer::COMPUTES , variable_registry, _task_name  );
      register_variable( m_eqn_names[ieqn], ArchesFieldContainer::REQUIRES , 0 , ArchesFieldContainer::OLDDW , variable_registry, _task_name  );
      register_variable( m_eqn_names[ieqn]+"_x_flux", ArchesFieldContainer::COMPUTES, variable_registry, _task_name );
      register_variable( m_eqn_names[ieqn]+"_y_flux", ArchesFieldContainer::COMPUTES, variable_registry, _task_name );
      register_variable( m_eqn_names[ieqn]+"_z_flux", ArchesFieldContainer::COMPUTES, variable_registry, _task_name );
    }
  }

  //------------------------------------------------------------------------------------------------
  template <typename T, typename PT> void
  KScalarRHS<T, PT>::timestep_init( const Patch* patch, ArchesTaskInfoManager* tsk_info ){

    typedef typename ArchesCore::VariableHelper<T>::ConstType CT;
    const int istart = 0;
    const int iend = m_eqn_names.size();

    for (int ieqn = istart; ieqn < iend; ieqn++ ){

      T& phi = tsk_info->get_uintah_field_add<T>( m_eqn_names[ieqn] );
      T& rhs = tsk_info->get_uintah_field_add<T>( m_eqn_names[ieqn]+"_rhs" );
      CT& old_phi = tsk_info->get_const_uintah_field_add<CT>( m_eqn_names[ieqn] );

      if ( m_premultiplier_name != "none" ){
        CT& old_rho_phi = tsk_info->get_const_uintah_field_add<CT>(m_premultiplier_name+m_eqn_names[ieqn] );
        T& rho_phi      = tsk_info->get_uintah_field_add<T>( m_premultiplier_name+m_eqn_names[ieqn] );
        rho_phi.copyData(old_rho_phi);
      }

      FXT& x_flux = tsk_info->get_uintah_field_add<FXT>(m_eqn_names[ieqn]+"_x_flux");
      FYT& y_flux = tsk_info->get_uintah_field_add<FYT>(m_eqn_names[ieqn]+"_y_flux");
      FZT& z_flux = tsk_info->get_uintah_field_add<FZT>(m_eqn_names[ieqn]+"_z_flux");

      phi.copyData(old_phi);
      rhs.initialize(0.0);
      x_flux.initialize(0.0);
      y_flux.initialize(0.0);
      z_flux.initialize(0.0);

    } //equation loop
  }

  //------------------------------------------------------------------------------------------------
  template <typename T, typename PT> void
  KScalarRHS<T, PT>::
  register_timestep_eval( std::vector<ArchesFieldContainer::VariableInformation>& variable_registry,
                          const int time_substep , const bool packed_tasks ){

    const int istart = 0;
    const int iend = m_eqn_names.size();
    for (int ieqn = istart; ieqn < iend; ieqn++ ){

      register_variable( m_eqn_names[ieqn], ArchesFieldContainer::REQUIRES, 1, ArchesFieldContainer::NEWDW, variable_registry, time_substep, _task_name );
      if ( m_premultiplier_name != "none" )
        register_variable( m_premultiplier_name+m_eqn_names[ieqn], ArchesFieldContainer::REQUIRES, 1, ArchesFieldContainer::NEWDW, variable_registry, time_substep, _task_name );
      register_variable( m_eqn_names[ieqn]+"_rhs", ArchesFieldContainer::MODIFIES, variable_registry, time_substep, _task_name );
      register_variable( m_eqn_names[ieqn]+"_x_flux", ArchesFieldContainer::MODIFIES, variable_registry, time_substep, _task_name );
      register_variable( m_eqn_names[ieqn]+"_y_flux", ArchesFieldContainer::MODIFIES, variable_registry, time_substep, _task_name );
      register_variable( m_eqn_names[ieqn]+"_z_flux", ArchesFieldContainer::MODIFIES, variable_registry, time_substep, _task_name );
      if ( m_do_diff[ieqn] ){
        register_variable( m_eqn_names[ieqn]+"_x_dflux", ArchesFieldContainer::REQUIRES, 1, ArchesFieldContainer::NEWDW, variable_registry, time_substep, _task_name );
        register_variable( m_eqn_names[ieqn]+"_y_dflux", ArchesFieldContainer::REQUIRES, 1, ArchesFieldContainer::NEWDW, variable_registry, time_substep, _task_name );
        register_variable( m_eqn_names[ieqn]+"_z_dflux", ArchesFieldContainer::REQUIRES, 1, ArchesFieldContainer::NEWDW, variable_registry, time_substep, _task_name );
      }
      if ( m_conv_scheme[ieqn] != NOCONV ){
        register_variable( m_eqn_names[ieqn]+"_x_psi", ArchesFieldContainer::REQUIRES, 1, ArchesFieldContainer::NEWDW, variable_registry, time_substep, _task_name, packed_tasks );
        register_variable( m_eqn_names[ieqn]+"_y_psi", ArchesFieldContainer::REQUIRES, 1, ArchesFieldContainer::NEWDW, variable_registry, time_substep, _task_name, packed_tasks );
        register_variable( m_eqn_names[ieqn]+"_z_psi", ArchesFieldContainer::REQUIRES, 1, ArchesFieldContainer::NEWDW, variable_registry, time_substep, _task_name, packed_tasks );
      }

      typedef std::vector<SourceInfo> VS;
      for (typename VS::iterator i = m_source_info[ieqn].begin(); i != m_source_info[ieqn].end(); i++){
        register_variable( i->name, ArchesFieldContainer::REQUIRES, 0, ArchesFieldContainer::NEWDW, variable_registry, time_substep, _task_name );
      }

    }

    //globally common variables
    register_variable( m_x_velocity_name, ArchesFieldContainer::REQUIRES, 1 , ArchesFieldContainer::LATEST, variable_registry, time_substep, _task_name );
    register_variable( m_y_velocity_name, ArchesFieldContainer::REQUIRES, 1 , ArchesFieldContainer::LATEST, variable_registry, time_substep, _task_name );
    register_variable( m_z_velocity_name, ArchesFieldContainer::REQUIRES, 1 , ArchesFieldContainer::LATEST, variable_registry, time_substep, _task_name );
    register_variable( m_eps_name, ArchesFieldContainer::REQUIRES, 1 , ArchesFieldContainer::OLDDW, variable_registry, time_substep, _task_name );

  }

  //------------------------------------------------------------------------------------------------
  template <typename T, typename PT> void
  KScalarRHS<T, PT>::eval( const Patch* patch, ArchesTaskInfoManager* tsk_info ){

    Vector Dx = patch->dCell();
    double ax = Dx.y() * Dx.z();
    double ay = Dx.z() * Dx.x();
    double az = Dx.x() * Dx.y();
    double V = Dx.x()*Dx.y()*Dx.z();

    CFXT& u = tsk_info->get_const_uintah_field_add<CFXT>(m_x_velocity_name);
    CFYT& v = tsk_info->get_const_uintah_field_add<CFYT>(m_y_velocity_name);
    CFZT& w = tsk_info->get_const_uintah_field_add<CFZT>(m_z_velocity_name);
    CT& eps = tsk_info->get_const_uintah_field_add<CT>(m_eps_name);
    Uintah::BlockRange range_cl_to_ech( patch->getCellLowIndex(),
                                        patch->getExtraCellHighIndex() );

    const int istart = 0;
    const int iend = m_eqn_names.size();
    for (int ieqn = istart; ieqn < iend; ieqn++ ){

      CT* rho_phi_ptr = nullptr;
      if ( m_premultiplier_name != "none" )
        rho_phi_ptr = tsk_info->get_const_uintah_field<CT>(m_premultiplier_name+m_eqn_names[ieqn]);

      CT& rho_phi = *rho_phi_ptr;
      CT& phi     = tsk_info->get_const_uintah_field_add<CT>(m_eqn_names[ieqn]);
      T& rhs      = tsk_info->get_uintah_field_add<T>(m_eqn_names[ieqn]+"_rhs");

      rhs.initialize(0.0);

      //Convection:
      FXT& x_flux = tsk_info->get_uintah_field_add<FXT>(m_eqn_names[ieqn]+"_x_flux");
      FYT& y_flux = tsk_info->get_uintah_field_add<FYT>(m_eqn_names[ieqn]+"_y_flux");
      FZT& z_flux = tsk_info->get_uintah_field_add<FZT>(m_eqn_names[ieqn]+"_z_flux");

      if ( m_conv_scheme[ieqn] != NOCONV ){

        FluxXT* x_psi_ptr;
        FluxYT* y_psi_ptr;
        FluxZT* z_psi_ptr;

        FieldTool<FluxXT> x_field_tool(tsk_info);
        x_psi_ptr = x_field_tool.get(m_eqn_names[ieqn]+"_x_psi");

        FieldTool<FluxYT> y_field_tool(tsk_info);
        y_psi_ptr = y_field_tool.get(m_eqn_names[ieqn]+"_y_psi");

        FieldTool<FluxZT> z_field_tool(tsk_info);
        z_psi_ptr = z_field_tool.get(m_eqn_names[ieqn]+"_z_psi");

        if ( m_premultiplier_name != "none" ){
          Uintah::ComputeConvectiveFlux<FluxXT, FluxYT, FluxZT >
            get_flux( rho_phi, u, v, w, (*x_psi_ptr), (*y_psi_ptr), (*z_psi_ptr),
                      x_flux, y_flux, z_flux, eps );

          Uintah::parallel_for( range_cl_to_ech, get_flux );
        } else {
          Uintah::ComputeConvectiveFlux<FluxXT, FluxYT, FluxZT >
            get_flux( phi, u, v, w, (*x_psi_ptr), (*y_psi_ptr), (*z_psi_ptr),
                      x_flux, y_flux, z_flux, eps );

          Uintah::parallel_for( range_cl_to_ech, get_flux );
        }

      }

      //Diffusion:
      if ( m_do_diff[ieqn] ) {

        CFXT& x_dflux = tsk_info->get_const_uintah_field_add<CFXT>(m_eqn_names[ieqn]+"_x_dflux");
        CFYT& y_dflux = tsk_info->get_const_uintah_field_add<CFYT>(m_eqn_names[ieqn]+"_y_dflux");
        CFZT& z_dflux = tsk_info->get_const_uintah_field_add<CFZT>(m_eqn_names[ieqn]+"_z_dflux");

        GET_EXTRACELL_BUFFERED_PATCH_RANGE(0,0);

        Uintah::BlockRange range_diff(low_patch_range, high_patch_range);

        Uintah::parallel_for( range_diff, [&](int i, int j, int k){

          rhs(i,j,k) += ax * ( x_dflux(i+1,j,k) - x_dflux(i,j,k) ) +
                        ay * ( y_dflux(i,j+1,k) - y_dflux(i,j,k) ) +
                        az * ( z_dflux(i,j,k+1) - z_dflux(i,j,k) );

        });
      }

      //Sources:
      typedef std::vector<SourceInfo> VS;
      for (typename VS::iterator isrc = m_source_info[ieqn].begin();
        isrc != m_source_info[ieqn].end(); isrc++){

        CT& src = *(tsk_info->get_const_uintah_field<CT>((*isrc).name));
        double weight = (*isrc).weight;
        Uintah::BlockRange src_range(patch->getCellLowIndex(), patch->getCellHighIndex());

        Uintah::parallel_for( src_range, [&](int i, int j, int k){

          rhs(i,j,k) += weight * src(i,j,k) * V;

        });
      }
    } // equation loop
  }

//--------------------------------------------------------------------------------------------------
  template <typename T, typename PT> void
  KScalarRHS<T, PT>::register_compute_bcs(
    std::vector<ArchesFieldContainer::VariableInformation>& variable_registry,
    const int time_substep , const bool packed_tasks){

    for ( auto i = m_eqn_names.begin(); i != m_eqn_names.end(); i++ ){
      register_variable( *i, ArchesFieldContainer::MODIFIES, variable_registry );
    }

    std::vector<std::string> bc_dep;
    m_boundary_functors->get_bc_dependencies( m_eqn_names, m_bcHelper, bc_dep );
    for ( auto i = bc_dep.begin(); i != bc_dep.end(); i++ ){
      register_variable( *i, ArchesFieldContainer::REQUIRES, 0 , ArchesFieldContainer::NEWDW,
                         variable_registry );
    }

    std::vector<std::string> bc_mod;
    m_boundary_functors->get_bc_modifies( m_eqn_names, m_bcHelper, bc_mod );
    for ( auto i = bc_mod.begin(); i != bc_mod.end(); i++ ){
      register_variable( *i, ArchesFieldContainer::MODIFIES, variable_registry );
    }

  }

//--------------------------------------------------------------------------------------------------
  template <typename T, typename PT> void
  KScalarRHS<T, PT>::compute_bcs( const Patch* patch, ArchesTaskInfoManager* tsk_info ){
    m_boundary_functors->apply_bc( m_eqn_names, m_bcHelper, tsk_info, patch );
  }
}
#endif
