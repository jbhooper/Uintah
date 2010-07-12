#ifndef Wasatch_GraphHelperTools_h
#define Wasatch_GraphHelperTools_h

#include <list>
#include <map>
#include <set>

#include <expression/Expr_ExpressionID.h>

/** \file */

namespace Expr{
  class ExpressionBuilder;
  class ExpressionFactory;
  class TransportEquation;
}

namespace Wasatch{

  /**
   *  \enum Category
   *  \brief defines the broad categories for various kinds of tasks.
   *
   *  Tasks associated with a particular category are registered in
   *  the associated Wasatch method that Uintah calls to register
   *  tasks.  They are generally combined into one or more Expression
   *  trees that are wrapped using the Wasatch::TaskInterface.
   */
  enum Category{
    INITIALIZATION,
    TIMESTEP_SELECTION,
    ADVANCE_SOLUTION
  };

  /**
   *  \brief a list of transport equations to be solved.
   */
  typedef std::list<Expr::TransportEquation*> TransEqns;

  /**
   *  \brief a set of ExpressionID generally to be used to store the
   *  "root" nodes of an ExpressionTree.
   */
  typedef std::set< Expr::ExpressionID > IDSet;

  /**
   *  \struct GraphHelper
   *  \author James C. Sutherland
   *  \date   June, 2010
   *
   *  \brief Contains information helpful in constructing graphs from ExprLib
   */
  struct GraphHelper
  {
    Expr::ExpressionFactory* const exprFactory;  ///< The factory used to build expressions
    IDSet rootIDs;                               ///< The root IDs used to create the graph
    GraphHelper( Expr::ExpressionFactory* ef );
  };

  /**
   *  Defines a map that provides GraphHelper objects given the \ref
   *  Wasatch::Category "Category" that they belong to.
   */
  typedef std::map< Category, GraphHelper* > GraphCategories;

} // namespace Wasatch


#endif // Wasatch_GraphHelperTools_h
