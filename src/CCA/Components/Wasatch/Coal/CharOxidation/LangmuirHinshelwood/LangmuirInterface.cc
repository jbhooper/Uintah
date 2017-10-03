#include "LangmuirInterface.h"

#include <stdexcept>
#include <sstream>

#include <CCA/Components/Wasatch/Coal/CharOxidation/CO2andCO_RHS.h>
#include <CCA/Components/Wasatch/Coal/CharOxidation/O2RHS.h>
#include <CCA/Components/Wasatch/Coal/CharOxidation/CharRHS.h>
#include <CCA/Components/Wasatch/Coal/CharOxidation/H2andH2ORHS.h>
#include <CCA/Components/Wasatch/Coal/Devolatilization/CPD/c0_fun.h>
#include <CCA/Components/Wasatch/Coal/CoalData.h>

#include "CharOxidation.h"
#include "GasificationCO2.h"
#include "GasificationH2O.h"

using std::ostringstream;
using std::endl;

namespace LH{

  using CHAR::CharModel;
  using CHAR::CharGasSpecies;
  using CHAR::CharOxidationData;
  using Expr::ExpressionID;

  //------------------------------------------------------------------

  template< typename FieldT >
  LangmuirInterface<FieldT>::
  LangmuirInterface( WasatchCore::GraphCategories& gc,
                     const Tag& pDiamTag,
                     const Tag& pTempTag,
                     const Tag& gTempTag,
                     const Tag& co2MassFracTag,
                     const Tag& o2MassFracTag,
                     const Tag& h2oMassFracTag,
                     const Tag& mixMWTag,
                     const Tag& pDensTag,
                     const Tag& gPressTag,
                     const Tag& pMassTag,
                     const Tag& pMass0Tag,
                     const Coal::CoalType coalType,
                     const DEV::DevModel dvmodel,
                     const CharModel chmodel )
    : pTempTag_      ( pTempTag          ),
      gTempTag_      ( gTempTag          ),
      mixMWTag_      ( mixMWTag          ),
      pDensTag_      ( pDensTag          ),
      gPressTag_     ( gPressTag         ),
      pDiamTag_      ( pDiamTag          ),
      pMassTag_      ( pMassTag          ),
      pMass0Tag_     ( pMass0Tag         ),
      o2MassFracTag_ ( o2MassFracTag     ),
      h2oMassFracTag_( h2oMassFracTag    ),
      co2MassFracTag_( co2MassFracTag    ),
      o2_rhsTag_     ( Coal::StringNames::self().char_o2_rhs,  Expr::STATE_NONE ),
      h2o_rhsTag_    ( Coal::StringNames::self().char_h2o_rhs, Expr::STATE_NONE ),
      h2_rhsTag_     ( Coal::StringNames::self().char_h2_rhs,  Expr::STATE_NONE ),
      co2_rhsTag_    ( Coal::StringNames::self().char_co2_rhs, Expr::STATE_NONE ),
      co_rhsTag_     ( Coal::StringNames::self().char_co_rhs,  Expr::STATE_NONE ),
      initDevChar_   ( dvmodel == DEV::CPDM ),
      charModel_     ( chmodel  ),
      charData_      ( coalType ),
      gc_            ( gc                ),
      sNames_        ( Coal::StringNames::self() )
  {
    std::cout<< "Setting up char model: "
             << CHAR::char_model_name(charModel_)
             << std::endl;

    parse_equations();
    set_tags();
    register_expressions();
  }

  //------------------------------------------------------------------

  template< typename FieldT >
  void
  LangmuirInterface<FieldT>::
  parse_equations()
  {
    std::cout<<"Parsing equations...\n";

        // Get initial mass fraction of char within coal volatiles
        double c0 = 0.0;
        if( initDevChar_ ){
          c0 = CPD::c0_fun(charData_.get_C(), charData_.get_O());
        }

         // Calculate initial mass fraction of char within coal
         double char0 = charData_.get_fixed_C()+ charData_.get_vm()*c0;

        std::cout << std::endl
                  << "Initial char mass fraction in coal volatiles is : "
                  << charData_.get_vm()*c0
                  << std::endl
                  << "Initial mass fraction of char in coal is        : "
                  << char0 << std::endl;

    charEqn_ = new Coal::CoalEquation( sNames_.char_mass, pMassTag_, char0, gc_ );

    eqns_.clear();
    eqns_.push_back( charEqn_ );
  }

  //------------------------------------------------------------------

  template< typename FieldT >
  void
  LangmuirInterface<FieldT>::
  set_tags()
  {
    charMassTag_      = charEqn_->solution_variable_tag();
    charMass_rhsTag_  = charEqn_->rhs_tag();
    oxidation_rhsTag_ = Tag( sNames_.char_oxid_rhs,  STATE_NONE );
    heteroCo2Tag_     = Tag( sNames_.char_gasifco2,  STATE_NONE );
    heteroH2oTag_     = Tag( sNames_.char_gasifh2o,  STATE_NONE );
    co2CoRatioTag_    = Tag( sNames_.char_coco2ratio,STATE_NONE );

    speciesSrcTags_.clear();
    speciesSrcTags_.push_back( co2_rhsTag_ );
    speciesSrcTags_.push_back( co_rhsTag_  );
    speciesSrcTags_.push_back( o2_rhsTag_  );
    speciesSrcTags_.push_back( h2_rhsTag_  );
    speciesSrcTags_.push_back( h2o_rhsTag_ );

    co2CoTags_.clear();
    co2CoTags_.push_back( co2_rhsTag_ );
    co2CoTags_.push_back( co_rhsTag_  );

    char_co2coTags_.clear();
    char_co2coTags_.push_back( oxidation_rhsTag_ );
    char_co2coTags_.push_back( co2CoRatioTag_    );

    h2andh2o_rhsTags_.clear();
    h2andh2o_rhsTags_.push_back( h2_rhsTag_  );
    h2andh2o_rhsTags_.push_back( h2o_rhsTag_ );
  }

  //------------------------------------------------------------------

  template< typename FieldT >
  void
  LangmuirInterface<FieldT>::
  register_expressions()
  {
    std::cout<<"Registering expressions...\n";

    Expr::ExpressionFactory& factory = *(gc_[WasatchCore::ADVANCE_SOLUTION]->exprFactory);

    oxidationRHSID_ = factory.register_expression( new typename CHAR::CharOxidation  <FieldT>::
                                                   Builder( char_co2coTags_, pDiamTag_, pTempTag_,
                                                            gTempTag_, o2MassFracTag_, mixMWTag_,
                                                            pDensTag_, charMassTag_, gPressTag_,
                                                            pMass0Tag_, charData_, charModel_) );

    gasifco2ID_     = factory.register_expression( new typename CHAR::GasificationCO2<FieldT>::
                                                   Builder( heteroCo2Tag_,charMassTag_, pDiamTag_,
                                                            pTempTag_,co2MassFracTag_, mixMWTag_,
                                                            gPressTag_) );

    gasifh2oID_     = factory.register_expression( new typename CHAR::GasificationH2O<FieldT>::
                                                   Builder( heteroH2oTag_, charMassTag_, pDiamTag_,
                                                            pTempTag_, h2oMassFracTag_, mixMWTag_,
                                                            gPressTag_) );

    h2Andh2oRHSID_  = factory.register_expression( new typename CHAR::H2andH2ORHS    <FieldT>::
                                                   Builder( h2andh2o_rhsTags_, heteroH2oTag_) );

    co2coRHSID_     = factory.register_expression( new typename CHAR::CO2andCO_RHS   <FieldT>::
                                                   Builder( co2CoTags_, oxidation_rhsTag_, co2CoRatioTag_,
                                                            heteroCo2Tag_, heteroH2oTag_) );

    o2RHSID_        = factory.register_expression( new typename O2RHS          <FieldT>::
                                                   Builder( o2_rhsTag_, oxidation_rhsTag_, co2CoRatioTag_) );

    charRHSID_      = factory.register_expression( new typename CHAR::CharRHS        <FieldT>::
                                                   Builder( charMass_rhsTag_, oxidation_rhsTag_,
                                                            heteroCo2Tag_, heteroH2oTag_) );
  }

  //------------------------------------------------------------------

  template< typename FieldT >
  const Tag
  LangmuirInterface<FieldT>::
  gas_species_src_tag( const CharGasSpecies spec ) const
  {
    if( spec == CHAR::O2 ) return o2_rhsTag_;
    if( spec == CHAR::CO2) return co2_rhsTag_;
    if( spec == CHAR::CO ) return co_rhsTag_;
    if( spec == CHAR::H2 ) return h2_rhsTag_;
    if( spec == CHAR::H2O) return h2o_rhsTag_;
    return Tag();
  }

//==========================================================================
// Explicit template instantiation for supported versions of this expression
template class LangmuirInterface< SpatialOps::Particle::ParticleField >;
//==========================================================================


} // namespace LH
