#ifndef TurbulentAggregationCoefficient_Expr_h
#define TurbulentAggregationCoefficient_Expr_h
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVStaggeredOperatorTypes.h>

#include <expression/Expression.h>

/**
 *  \ingroup WasatchExpressions
 *  \class TurbulentAggregationCoefficient
 *  \author Alex Abboud
 *  \date June 2012
 *
 *  \tparam FieldT the type of field.
 *
 *  \brief Calculates the coefficent used for Turbulent aggregation
 *  \f$ (4/3)*(3 \pi /10)^{1/2} (\epsilon / \nu ) ^{1/2} \f$
 *  \f$ \epsilon \f$ is the energy dissipation, \f$ \nu \f$ is the kinematic viscosity
 */
template< typename FieldT >
class TurbulentAggregationCoefficient
: public Expr::Expression<FieldT>
{
  const Expr::Tag kinViscTag_, dissipationTag_;  //kinematic viscoisty and energy dissipation tags
  const double coefVal_;
  const FieldT* kinVisc_;
  const FieldT* dissipation_;
  
  TurbulentAggregationCoefficient( const Expr::Tag& kinViscTag,
                                   const Expr::Tag& dissipationTag,
                                   const double coefVal );
  
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    Builder( const Expr::Tag& result,
             const Expr::Tag& kinViscTag,
             const Expr::Tag& dissipationTag,
             const double coefVal )
    : ExpressionBuilder(result),
    kinvisct_(kinViscTag),
    dissipationt_(dissipationTag),
    coefval_(coefVal)
    {}
    
    ~Builder(){}
    
    Expr::ExpressionBase* build() const
    {
      return new TurbulentAggregationCoefficient<FieldT>( kinvisct_, dissipationt_, coefval_);
    }
    
  private:
    const Expr::Tag kinvisct_ ;
    const Expr::Tag dissipationt_;
    const double coefval_;
  };
  
  ~TurbulentAggregationCoefficient();
  
  void advertise_dependents( Expr::ExprDeps& exprDeps );
  void bind_fields( const Expr::FieldManagerList& fml );
  void bind_operators( const SpatialOps::OperatorDatabase& opDB );
  void evaluate();
  
};



// ###################################################################
//
//                          Implementation
//
// ###################################################################



template< typename FieldT >
TurbulentAggregationCoefficient<FieldT>::
TurbulentAggregationCoefficient( const Expr::Tag& kinViscTag,
                                 const Expr::Tag& dissipationTag,
                                 const double coefVal )
: Expr::Expression<FieldT>(),
kinViscTag_(kinViscTag),
dissipationTag_(dissipationTag),
coefVal_(coefVal)
{}

//--------------------------------------------------------------------

template< typename FieldT >
TurbulentAggregationCoefficient<FieldT>::
~TurbulentAggregationCoefficient()
{}

//--------------------------------------------------------------------

template< typename FieldT >
void
TurbulentAggregationCoefficient<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( kinViscTag_ );
  exprDeps.requires_expression( dissipationTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
TurbulentAggregationCoefficient<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  const typename Expr::FieldMgrSelector<FieldT>::type& fm = fml.template field_manager<FieldT>();
  kinVisc_ = &fm.field_ref( kinViscTag_ );
  dissipation_ = &fm.field_ref( dissipationTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
TurbulentAggregationCoefficient<FieldT>::
bind_operators( const SpatialOps::OperatorDatabase& opDB )
{}

//--------------------------------------------------------------------

template< typename FieldT >
void
TurbulentAggregationCoefficient<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& result = this->value();
  result <<= coefVal_ * sqrt( *dissipation_ / *kinVisc_ );
}

//--------------------------------------------------------------------

#endif // TurbulentAggregationCoefficient_Expr_h
