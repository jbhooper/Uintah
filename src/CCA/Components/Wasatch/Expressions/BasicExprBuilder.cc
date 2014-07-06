/*
 * The MIT License
 *
 * Copyright (c) 2012 The University of Utah
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

//-- Uintah includes --//
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>


//-- Wasatch includes --//
#include <CCA/Components/Wasatch/Expressions/BasicExprBuilder.h>
#include <CCA/Components/Wasatch/FieldAdaptor.h>
#include <CCA/Components/Wasatch/FieldTypes.h>
#include <CCA/Components/Wasatch/ParseTools.h>
#include <CCA/Components/Wasatch/Expressions/MMS/TaylorVortex.h>
#include <CCA/Components/Wasatch/Expressions/MMS/Varden2DMMS.h>
#include <CCA/Components/Wasatch/Expressions/MMS/Functions.h>
#include <CCA/Components/Wasatch/Expressions/MMS/VardenMMS.h>
#include <CCA/Components/Wasatch/Expressions/ExprAlgebra.h>
#include <CCA/Components/Wasatch/Expressions/TimeDerivative.h>
#include <CCA/Components/Wasatch/Expressions/GeometryBased.h>
#include <CCA/Components/Wasatch/Expressions/Turbulence/WallDistance.h>
#include <CCA/Components/Wasatch/OldVariable.h>

#include <CCA/Components/Wasatch/TagNames.h>

#include <CCA/Components/Wasatch/Expressions/PBE/BrownianAggregationCoefficient.h>
#include <CCA/Components/Wasatch/Expressions/PBE/TurbulentAggregationCoefficient.h>
#include <CCA/Components/Wasatch/Expressions/PBE/MultiEnvMixingModel.h>
#include <CCA/Components/Wasatch/Expressions/PBE/Precipitation/PrecipitationBulkDiffusionCoefficient.h>
#include <CCA/Components/Wasatch/Expressions/PBE/Precipitation/PrecipitationMonosurfaceCoefficient.h>
#include <CCA/Components/Wasatch/Expressions/PBE/Precipitation/PrecipitationClassicNucleationCoefficient.h>
#include <CCA/Components/Wasatch/Expressions/PBE/Precipitation/PrecipitationRCritical.h>
#include <CCA/Components/Wasatch/Expressions/PBE/Precipitation/PrecipitationSource.h>
#include <CCA/Components/Wasatch/Expressions/PBE/Precipitation/ParticleVolumeFraction.h>
#include <CCA/Components/Wasatch/Expressions/PBE/Precipitation/PrecipitateEffectiveViscosity.h>
#include <CCA/Components/Wasatch/Expressions/PBE/Precipitation/CylindricalDiffusionCoefficient.h>
#include <CCA/Components/Wasatch/Expressions/PBE/Precipitation/KineticGrowthCoefficient.h>
#include <CCA/Components/Wasatch/Expressions/PBE/Precipitation/HomogeneousNucleationCoefficient.h>
#include <CCA/Components/Wasatch/Expressions/PBE/Precipitation/CriticalSurfaceEnergy.h>

#include <CCA/Components/Wasatch/Expressions/PostProcessing/Vorticity.h>
#include <CCA/Components/Wasatch/Expressions/PostProcessing/KineticEnergy.h>
#include <CCA/Components/Wasatch/Expressions/PostProcessing/VelocityMagnitude.h>
#include <CCA/Components/Wasatch/Expressions/PostProcessing/InterpolateExpression.h>

#include <CCA/Components/Wasatch/Expressions/Particles/ParticleInitialization.h>


// BC Expressions Includes
#include <CCA/Components/Wasatch/Expressions/BoundaryConditions/BoundaryConditions.h>
#include <CCA/Components/Wasatch/Expressions/BoundaryConditions/TurbulentInletBC.h>
#include <CCA/Components/Wasatch/Expressions/BoundaryConditions/BoundaryConditionBase.h>
#include <CCA/Components/Wasatch/Expressions/BoundaryConditions/VardenMMSBCs.h>

//-- ExprLib includes --//
#include <expression/ExprLib.h>


#include <string>
#ifndef PI
#  define PI 3.1415926535897932384626433832795
#endif

using std::endl;

namespace Wasatch{
  
  //------------------------------------------------------------------
  // Special parser for particle expressions
  
  template<typename FieldT>
  Expr::ExpressionBuilder*
  build_basic_particle_expr( Uintah::ProblemSpecP params )
  {
    const Expr::Tag tag = parse_nametag( params->findBlock("NameTag") );    
    Expr::ExpressionBuilder* builder = NULL;

    if( params->findBlock("Constant") ){
      double val;  params->get("Constant",val);
      typedef typename Expr::ConstantExpr<FieldT>::Builder Builder;
      builder = scinew Builder( tag, val );
    } else if (params->findBlock("ParticlePositionIC")) {
      Uintah::ProblemSpecP valParams = params->findBlock("ParticlePositionIC");
      // parse coordinate
      std::string coord;
      valParams->getAttribute("coordinate",coord);
      // check what type of bounds we are using: specified or patch based?
      std::string boundsType;
      valParams->getAttribute("bounds",boundsType);
      const bool usePatchBounds = (boundsType == "PATCHBASED");
      double lo = 0.0, hi = 1.0;

      if (valParams->findBlock("Bounds")) {
        valParams->findBlock("Bounds")->getAttribute("low", lo);
        valParams->findBlock("Bounds")->getAttribute("high", hi);
      }
      
      if (valParams->findBlock("Uniform")) {
        bool transverse = false;
        valParams->findBlock("Uniform")->getAttribute("transversedir", transverse);
        builder = scinew ParticleUniformIC::Builder( tag, lo, hi, transverse, coord, usePatchBounds );
      } else if (valParams->findBlock("Random")) {
        double seed = 0.0;
        valParams->findBlock("Random")->getAttribute("seed",seed);
        builder = scinew ParticleRandomIC::Builder( tag, coord, lo, hi, seed, usePatchBounds );
      }
      
    } else if ( params->findBlock("RandomField") ) {
      Uintah::ProblemSpecP valParams = params->findBlock("RandomField");
      double lo, hi, seed;
      valParams->getAttribute("low",lo);
      valParams->getAttribute("high",hi);
      valParams->getAttribute("seed",seed);
      const std::string coord="";
      const bool usePatchBounds = false;
      builder = scinew ParticleRandomIC::Builder( tag, coord, lo, hi, seed, usePatchBounds );
    } else if( params->findBlock("LinearFunction") ){
      double slope, intercept;
      Uintah::ProblemSpecP valParams = params->findBlock("LinearFunction");
      valParams->getAttribute("slope",slope);
      valParams->getAttribute("intercept",intercept);
      const Expr::Tag indepVarTag = parse_nametag( valParams->findBlock("NameTag") );
      typedef typename Expr::LinearFunction<FieldT>::Builder Builder;
      builder = scinew Builder( tag, indepVarTag, slope, intercept );
    } else if ( params->findBlock("SineFunction") ) {
      double amplitude, frequency, offset;
      Uintah::ProblemSpecP valParams = params->findBlock("SineFunction");
      valParams->getAttribute("amplitude",amplitude);
      valParams->getAttribute("frequency",frequency);
      valParams->getAttribute("offset",offset);
      const Expr::Tag indepVarTag = parse_nametag( valParams->findBlock("NameTag") );
      typedef typename Expr::SinFunction<FieldT>::Builder Builder;
      builder = scinew Builder( tag, indepVarTag, amplitude, frequency, offset);
    } else {
      std::ostringstream msg;
      msg << "ERROR: unsupported BasicExpression for Particles. Note that not all BasicExpressions are supported by particles. Please revise your input file." << std::endl;
      throw Uintah::ProblemSetupException( msg.str(), __FILE__, __LINE__ );
    }
    return builder;
  }

  //------------------------------------------------------------------
  
  template<typename FieldT>
  Expr::ExpressionBuilder*
  build_basic_expr( Uintah::ProblemSpecP params )
  {
    const Expr::Tag tag = parse_nametag( params->findBlock("NameTag") );
    
    const TagNames& tagNames = TagNames::self();

    Expr::ExpressionBuilder* builder = NULL;
    
    //    std::string exprType;
    //    Uintah::ProblemSpecP valParams = params->get( "value", exprType );
    if( params->findBlock("Constant") ){
      double val;  params->get("Constant",val);
      typedef typename Expr::ConstantExpr<FieldT>::Builder Builder;
      builder = scinew Builder( tag, val );
    }
    else if( params->findBlock("LinearFunction") ){
      double slope, intercept;
      Uintah::ProblemSpecP valParams = params->findBlock("LinearFunction");
      valParams->getAttribute("slope",slope);
      valParams->getAttribute("intercept",intercept);
      const Expr::Tag indepVarTag = parse_nametag( valParams->findBlock("NameTag") );
      typedef typename Expr::LinearFunction<FieldT>::Builder Builder;
      builder = scinew Builder( tag, indepVarTag, slope, intercept );
    }
    
    else if ( params->findBlock("SineFunction") ) {
      double amplitude, frequency, offset;
      Uintah::ProblemSpecP valParams = params->findBlock("SineFunction");
      valParams->getAttribute("amplitude",amplitude);
      valParams->getAttribute("frequency",frequency);
      valParams->getAttribute("offset",offset);
      const Expr::Tag indepVarTag = parse_nametag( valParams->findBlock("NameTag") );
      typedef typename Expr::SinFunction<FieldT>::Builder Builder;
      builder = scinew Builder( tag, indepVarTag, amplitude, frequency, offset);
    }
    
    else if ( params->findBlock("ParabolicFunction") ) {
      double a, b, c;
      Uintah::ProblemSpecP valParams = params->findBlock("ParabolicFunction");
      valParams->getAttribute("a",a);
      valParams->getAttribute("b",b);
      valParams->getAttribute("c",c);
      const Expr::Tag indepVarTag = parse_nametag( valParams->findBlock("NameTag") );
      typedef typename Expr::ParabolicFunction<FieldT>::Builder Builder;
      builder = scinew Builder( tag, indepVarTag, a, b, c);
    }
    
    else if ( params->findBlock("GaussianFunction") ) {
      double amplitude, deviation, mean, baseline;
      Uintah::ProblemSpecP valParams = params->findBlock("GaussianFunction");
      valParams->getAttribute("amplitude",amplitude);
      valParams->getAttribute("deviation",deviation);
      valParams->getAttribute("mean",mean);
      valParams->getAttribute("baseline",baseline);
      const Expr::Tag indepVarTag = parse_nametag( valParams->findBlock("NameTag") );
      typedef typename Expr::GaussianFunction<FieldT>::Builder Builder;
      builder = scinew Builder( tag, indepVarTag, amplitude, deviation, mean, baseline);
    }
    
    else if ( params->findBlock("DoubleTanhFunction") ) {
      double amplitude, width, midpointUp, midpointDown;
      Uintah::ProblemSpecP valParams = params->findBlock("DoubleTanhFunction");
      valParams->getAttribute("amplitude",amplitude);
      valParams->getAttribute("width",width);
      valParams->getAttribute("midpointUp",midpointUp);
      valParams->getAttribute("midpointDown",midpointDown);
      const Expr::Tag indepVarTag = parse_nametag( valParams->findBlock("NameTag") );
      typedef typename Expr::DoubleTanhFunction<FieldT>::Builder Builder;
      builder = scinew Builder( tag, indepVarTag, midpointUp, midpointDown, width, amplitude);
    }
    
    else if ( params->findBlock("SineTime") ) {
      const Expr::Tag timeVarTag( "time", Expr::STATE_NONE );
      typedef typename SineTime<FieldT>::Builder Builder;
      builder = scinew Builder( tag, timeVarTag );
    }
    
    else if ( params->findBlock("ExprAlgebra") ) {
      std::string algebraicOperation;
      Uintah::ProblemSpecP valParams = params->findBlock("ExprAlgebra");
      
      
      Expr::TagList srcFieldTagList;
      
      for( Uintah::ProblemSpecP exprParams = valParams->findBlock("NameTag");
          exprParams != 0;
          exprParams = exprParams->findNextBlock("NameTag") )
      {
        srcFieldTagList.push_back( parse_nametag( exprParams ) );
      }
      
      valParams->getAttribute("algebraicOperation",algebraicOperation);
      
      // for now, only support parsing for fields of same type.  In the future,
      // we could extend parsing support for differing source field types.
      typedef ExprAlgebra<FieldT> AlgExpr;
      typename AlgExpr::OperationType optype;
      if      (algebraicOperation == "SUM"       ) optype = AlgExpr::SUM;
      else if (algebraicOperation == "DIFFERENCE") optype = AlgExpr::DIFFERENCE;
      else if (algebraicOperation == "PRODUCT"   ) optype = AlgExpr::PRODUCT;
      else {
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__ << std::endl
        << "ERROR: The operator " << algebraicOperation
        << " is not supported in ExprAlgebra." << std::endl;
        throw std::invalid_argument( msg.str() );
      }
      builder = scinew typename AlgExpr::Builder( tag, srcFieldTagList, optype );
    }    
    
    else if( params->findBlock("WallDistanceFunction") ){
      Uintah::ProblemSpecP valParams = params->findBlock("WallDistanceFunction");
      const Expr::Tag indepVarTag = parse_nametag( valParams->findBlock("NameTag") );
      typedef typename WallDistance::Builder Builder;
      builder = scinew Builder( tag, indepVarTag );
    }
    
    else if( params->findBlock("ReadFromFile") ){
      std::string fieldType;
      params->getAttribute("type",fieldType);
      
      Uintah::ProblemSpecP valParams = params->findBlock("ReadFromFile");
      std::string fileName;
      valParams->get("FileName",fileName);
      
      const Expr::Tag xTag("X" + fieldType, Expr::STATE_NONE);
      const Expr::Tag yTag("Y" + fieldType, Expr::STATE_NONE);
      const Expr::Tag zTag("Z" + fieldType, Expr::STATE_NONE);
      
      typedef typename ReadFromFileExpression<FieldT>::Builder Builder;
      builder = scinew Builder( tag, xTag, yTag, zTag, fileName );
    }
    
    else if ( params->findBlock("StepFunction") ) {
      Uintah::ProblemSpecP valParams = params->findBlock("StepFunction");
      double transitionPoint, lowValue, highValue;
      valParams->getAttribute("transitionPoint",transitionPoint);
      valParams->getAttribute("lowValue",lowValue);
      valParams->getAttribute("highValue",highValue);
      const Expr::Tag indepVarTag = parse_nametag( valParams->findBlock("NameTag") );
      typedef typename StepFunction<FieldT>::Builder Builder;
      builder = scinew Builder( tag, indepVarTag, transitionPoint, lowValue, highValue );
    }
    
    else if ( params->findBlock("RayleighTaylor") ) {
      Uintah::ProblemSpecP valParams = params->findBlock("RayleighTaylor");
      double transitionPoint, lowValue, highValue, frequency, amplitude;
      std::string x1, x2;
      valParams->getAttribute("transitionPoint",transitionPoint);
      valParams->getAttribute("lowValue",lowValue);
      valParams->getAttribute("highValue",highValue);
      valParams->getAttribute("frequency",frequency);
      valParams->getAttribute("amplitude",amplitude);
      
      valParams->getAttribute("x1",x1);
      valParams->getAttribute("x2",x2);
      const Expr::Tag x1Tag(x1,Expr::STATE_NONE);
      const Expr::Tag x2Tag(x2,Expr::STATE_NONE);
      const Expr::Tag indepVarTag = parse_nametag( valParams->findBlock("NameTag") );
      typedef typename RayleighTaylor<FieldT>::Builder Builder;
      builder = scinew Builder( tag, indepVarTag, x1Tag, x2Tag, transitionPoint, lowValue, highValue, frequency, amplitude );
    }

    else if ( params->findBlock("VarDen1DMMSMixFracSrc") ) {
      Uintah::ProblemSpecP valParams = params->findBlock("VarDen1DMMSMixFracSrc");
      double D, rho0, rho1;
      valParams->getAttribute("D",    D);
      valParams->getAttribute("rho0", rho0);
      valParams->getAttribute("rho1", rho1);
      const Expr::Tag xTag = parse_nametag( valParams->findBlock("Coordinate")->findBlock("NameTag") );
      typedef typename VarDen1DMMSMixFracSrc<FieldT>::Builder Builder;
      builder = scinew Builder( tag, xTag, tagNames.time, D, rho0, rho1 );
    }
    
    else if ( params->findBlock("ExponentialVortex") ) {
      Uintah::ProblemSpecP valParams = params->findBlock("ExponentialVortex");
      double x0, y0, G, R, U, V;
      std::string velocityComponent;
      
      valParams->getAttribute("x0",x0);
      valParams->getAttribute("y0",y0);
      valParams->getAttribute("G",G);
      valParams->getAttribute("R",R);
      valParams->getAttribute("U",U);
      valParams->getAttribute("V",V);
      valParams->getAttribute("velocityComponent",velocityComponent);
      
      typedef ExponentialVortex<FieldT> ExpVortex;
      typename ExpVortex::VelocityComponent velComponent;
      if      (velocityComponent == "X1") velComponent = ExpVortex::X1;
      else if (velocityComponent == "X2") velComponent = ExpVortex::X2;
      else {
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__ << std::endl
        << "ERROR: The velocity component " << velocityComponent
        << " is not supported in the ExponentialVortex expression." << std::endl;
        throw std::invalid_argument( msg.str() );
      }
      
      const Expr::Tag xTag = parse_nametag( valParams->findBlock("Coordinate1")->findBlock("NameTag") );
      const Expr::Tag yTag = parse_nametag( valParams->findBlock("Coordinate2")->findBlock("NameTag") );
      
      builder = scinew typename ExpVortex::Builder( tag, xTag, yTag, x0, y0, G, R, U, V, velComponent );
    }
    
    else if ( params->findBlock("LambsDipole") ) {
      Uintah::ProblemSpecP valParams = params->findBlock("LambsDipole");
      double x0, y0, G, R, U;
      std::string velocityComponent;
      
      valParams->getAttribute("x0",x0);
      valParams->getAttribute("y0",y0);
      valParams->getAttribute("G",G);
      valParams->getAttribute("R",R);
      valParams->getAttribute("U",U);
      valParams->getAttribute("velocityComponent",velocityComponent);
      
      typedef LambsDipole<FieldT> Dipole;
      typename Dipole::VelocityComponent velComponent;
      if      (velocityComponent == "X1") velComponent = Dipole::X1;
      else if (velocityComponent == "X2") velComponent = Dipole::X2;
      else {
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__ << std::endl
        << "ERROR: The velocity component " << velocityComponent
        << " is not supported in the ExponentialVortex expression." << std::endl;
        throw std::invalid_argument( msg.str() );
      }
      
      const Expr::Tag xTag = parse_nametag( valParams->findBlock("Coordinate1")->findBlock("NameTag") );
      const Expr::Tag yTag = parse_nametag( valParams->findBlock("Coordinate2")->findBlock("NameTag") );
      builder = scinew typename Dipole::Builder( tag, xTag, yTag, x0, y0, G, R, U, velComponent );
    }
    
    else if( params->findBlock("RandomField") ){
      Uintah::ProblemSpecP valParams = params->findBlock("RandomField");
      double low, high, seed;
      valParams->getAttribute("low",low);
      valParams->getAttribute("high",high);
      valParams->getAttribute("seed",seed);
      typedef typename RandomField<FieldT>::Builder Builder;
      builder = scinew Builder( tag, low, high, seed );
    }
    
    else if( params->findBlock("TimeDerivative") ){
      Uintah::ProblemSpecP valParams = params->findBlock("TimeDerivative");
      const Expr::Tag srcTag = parse_nametag( valParams->findBlock("NameTag") );
      // create an old variable
      OldVariable& oldVar = OldVariable::self();
      oldVar.add_variable<FieldT>( ADVANCE_SOLUTION, srcTag);
      
      Expr::Tag srcOldTag = Expr::Tag( srcTag.name() + "_old", Expr::STATE_NONE );
      const TagNames& tagNames = TagNames::self();
      typedef typename TimeDerivative<FieldT>::Builder Builder;
      builder = scinew Builder( tag, srcTag, srcOldTag, tagNames.dt );
    }
    
    else if ( params->findBlock("BurnsChristonAbskg") ){
      typedef BurnsChristonAbskg<FieldT> BurnsChristonAbskgExpr;
      std::string fieldType;
      params->getAttribute("type",fieldType);
      const Expr::Tag xTag("X" + fieldType, Expr::STATE_NONE);
      const Expr::Tag yTag("Y" + fieldType, Expr::STATE_NONE);
      const Expr::Tag zTag("Z" + fieldType, Expr::STATE_NONE);
      
      builder = scinew typename BurnsChristonAbskgExpr::Builder( tag, xTag, yTag, zTag  );
    }

    else if ( params->findBlock("GeometryBased") ) {
      std::multimap <Uintah::GeometryPieceP, double > geomObjectsMap;
      double outsideValue = 1.0;
      Uintah::ProblemSpecP geomBasedSpec = params->findBlock("GeometryBased");
      geomBasedSpec->getAttribute("value", outsideValue);
      // parse all intrusions
      std::vector<Uintah::GeometryPieceP> geomObjects;
      
      for( Uintah::ProblemSpecP intrusionParams = geomBasedSpec->findBlock("Intrusion");
          intrusionParams != 0;
          intrusionParams = intrusionParams->findNextBlock("Intrusion") )
      {
        std::cout << "intrusion " << std::endl;
        Uintah::GeometryPieceFactory::create(intrusionParams->findBlock("geom_object"),geomObjects);
        double insideValue = 0.0;
        intrusionParams->getAttribute("value", insideValue);
        geomObjectsMap.insert(std::pair<Uintah::GeometryPieceP, double>(geomObjects.back(), insideValue)); // set a value inside the geometry object
      }
      builder = scinew typename GeometryBased<FieldT>::Builder(tag, geomObjectsMap, outsideValue);
    }

    return builder;
  }
  
  //------------------------------------------------------------------
  
  template<typename FieldT>
  Expr::ExpressionBuilder*
  build_taylor_vortex_mms_expr( Uintah::ProblemSpecP params )
  {
    const Expr::Tag tag = parse_nametag( params->findBlock("NameTag") );
    
    const TagNames& tagNames = TagNames::self();

    Expr::ExpressionBuilder* builder = NULL;
    
    //std::string exprType;
    //bool valParams = params->get("value",exprType);
    
    if( params->findBlock("VelocityX") ){
      double amplitude,viscosity;
      Uintah::ProblemSpecP valParams = params->findBlock("VelocityX");
      valParams->getAttribute("amplitude",amplitude);
      valParams->getAttribute("viscosity",viscosity);
      const Expr::Tag indepVarTag1 = parse_nametag( valParams->findBlock("XCoordinate")->findBlock("NameTag") );
      const Expr::Tag indepVarTag2 = parse_nametag( valParams->findBlock("YCoordinate")->findBlock("NameTag") );
      typedef typename VelocityX<FieldT>::Builder Builder;
      builder = scinew Builder( tag, indepVarTag1, indepVarTag2, tagNames.time, amplitude, viscosity );
    }
    
    else if( params->findBlock("VelocityY") ){
      double amplitude, viscosity;
      Uintah::ProblemSpecP valParams = params->findBlock("VelocityY");
      valParams->getAttribute("amplitude",amplitude);
      valParams->getAttribute("viscosity",viscosity);
      const Expr::Tag indepVarTag1 = parse_nametag( valParams->findBlock("XCoordinate")->findBlock("NameTag") );
      const Expr::Tag indepVarTag2 = parse_nametag( valParams->findBlock("YCoordinate")->findBlock("NameTag") );
      typedef typename VelocityY<FieldT>::Builder Builder;
      builder = scinew Builder( tag, indepVarTag1, indepVarTag2, tagNames.time, amplitude, viscosity );
    }
    
    else if( params->findBlock("GradPX") ){
      double amplitude, viscosity;
      Uintah::ProblemSpecP valParams = params->findBlock("GradPX");
      valParams->getAttribute("amplitude",amplitude);
      valParams->getAttribute("viscosity",viscosity);
      const Expr::Tag indepVarTag1 = parse_nametag( valParams->findBlock("XCoordinate")->findBlock("NameTag") );
      const Expr::Tag indepVarTag2 = parse_nametag( valParams->findBlock("YCoordinate")->findBlock("NameTag") );
      typedef typename GradPX<FieldT>::Builder Builder;
      builder = scinew Builder( tag, indepVarTag1, indepVarTag2, tagNames.time, amplitude, viscosity );
    }
    
    else if( params->findBlock("GradPY") ){
      double amplitude, viscosity;
      Uintah::ProblemSpecP valParams = params->findBlock("GradPY");
      valParams->getAttribute("amplitude",amplitude);
      valParams->getAttribute("viscosity",viscosity);
      const Expr::Tag indepVarTag1 = parse_nametag( valParams->findBlock("XCoordinate")->findBlock("NameTag") );
      const Expr::Tag indepVarTag2 = parse_nametag( valParams->findBlock("YCoordinate")->findBlock("NameTag") );
      typedef typename GradPY<FieldT>::Builder Builder;
      builder = scinew Builder( tag, indepVarTag1, indepVarTag2, tagNames.time, amplitude, viscosity );
    }
    
    else if( params->findBlock("TGVel3D") ){
      double angle;
      std::string velComponent;
      Uintah::ProblemSpecP valParams = params->findBlock("TGVel3D");
      valParams->getAttribute("angle",angle);
      valParams->get("VelComponent", velComponent);
      const Expr::Tag XCoordinate = parse_nametag( valParams->findBlock("XCoordinate")->findBlock("NameTag") );
      const Expr::Tag YCoordinate = parse_nametag( valParams->findBlock("YCoordinate")->findBlock("NameTag") );
      const Expr::Tag ZCoordinate = parse_nametag( valParams->findBlock("ZCoordinate")->findBlock("NameTag") );
      typedef typename TaylorGreenVel3D<FieldT>::Builder Builder;
      // shuffle the x, y, and z coordinates based on the velocity component
      if (velComponent=="X") {
        angle += 2*PI/3.0;
        builder = scinew Builder( tag, XCoordinate, YCoordinate, ZCoordinate, angle );
      } else if (velComponent=="Y") {
        angle -= 2*PI/3.0;
        builder = scinew Builder( tag, YCoordinate, XCoordinate, ZCoordinate, angle );
      } else if (velComponent=="Z") {
        builder = scinew Builder( tag, ZCoordinate, XCoordinate, YCoordinate, angle );
      }
    }
    
    return builder;
  }

  //------------------------------------------------------------------
  
  template<typename FieldT>
  Expr::ExpressionBuilder*
  build_precipitation_expr( Uintah::ProblemSpecP params , Uintah::ProblemSpecP wasatchParams )
  {
    const Expr::Tag tag = parse_nametag( params->findBlock("NameTag") );
    Expr::ExpressionBuilder* builder = NULL;

    const double kB = 1.3806488e-23;
    const double nA = 6.023e23;
    const double R = 8.314;

    if (params->findBlock("PrecipitationBulkDiffusionCoefficient") ) {
      double coef, molecularVolume, diffusionCoefficient, sMin;
      Uintah::ProblemSpecP coefParams = params->findBlock("PrecipitationBulkDiffusionCoefficient");
      coefParams -> getAttribute("Molec_Vol",molecularVolume);
      coefParams -> getAttribute("Diff_Coef",diffusionCoefficient);
      sMin = 0.0;
      if (coefParams->getAttribute("S_Min", sMin) )
        coefParams->getAttribute("S_Min", sMin);
      coef = molecularVolume*diffusionCoefficient;
      Expr::Tag sBarTag;
      if (coefParams->findBlock("SBar") ) 
        sBarTag = parse_nametag( coefParams->findBlock("SBar")->findBlock("NameTag") );
      const Expr::Tag saturationTag = parse_nametag( coefParams->findBlock("Supersaturation")->findBlock("NameTag") );
      const Expr::Tag eqTag  = parse_nametag( coefParams->findBlock("EquilibriumConcentration")->findBlock("NameTag") );
      typedef typename PrecipitationBulkDiffusionCoefficient<FieldT>::Builder Builder;
      builder = scinew Builder(tag, saturationTag, eqTag, sBarTag, coef, sMin);
    }
    
    else if (params->findBlock("CylindricalDiffusionCoefficient") ) {
      double coef, molecularVolume, diffusionCoefficient, sMin;
      Uintah::ProblemSpecP coefParams = params->findBlock("CylindricalDiffusionCoefficient");
      coefParams -> getAttribute("Molec_Vol",molecularVolume);
      coefParams -> getAttribute("Diff_Coef",diffusionCoefficient);
      sMin = 0.0;
      if (coefParams->getAttribute("S_Min", sMin) )
        coefParams->getAttribute("S_Min", sMin);
      coef = molecularVolume*diffusionCoefficient* 7.0/6.0/log(0.5);
      const Expr::Tag saturationTag = parse_nametag( coefParams->findBlock("Supersaturation")->findBlock("NameTag") );
      const Expr::Tag eqTag  = parse_nametag( coefParams->findBlock("EquilibriumConcentration")->findBlock("NameTag") );
      Expr::Tag sBarTag;
      if (coefParams->findBlock("SBar") ) 
        sBarTag = parse_nametag( coefParams->findBlock("SBar")->findBlock("NameTag") );
      typedef typename CylindricalDiffusionCoefficient<FieldT>::Builder Builder;
      builder = scinew Builder( tag, saturationTag, eqTag, sBarTag, coef, sMin);
    }
    
    else if (params->findBlock("KineticGrowthCoefficient") ) {
      double coef, sMax, sMin;
      Uintah::ProblemSpecP coefParams = params->findBlock("KineticGrowthCoefficient");
      coefParams -> getAttribute("K_A",coef);
      sMin = 0.0;
      sMax = 1e10;
      if( coefParams->getAttribute("S_Max",sMax) )
        coefParams->getAttribute("S_Max",sMax);
      if( coefParams->getAttribute("S_Min",sMin) )
        coefParams->getAttribute("S_Min",sMin);
      const Expr::Tag saturationTag = parse_nametag( coefParams->findBlock("Supersaturation")->findBlock("NameTag") );
      Expr::Tag sBarTag;
      if (coefParams->findBlock("SBar") ) 
        sBarTag = parse_nametag( coefParams->findBlock("SBar")->findBlock("NameTag") );
      typedef typename KineticGrowthCoefficient<FieldT>::Builder Builder;
      builder = scinew Builder( tag, saturationTag, sBarTag, coef, sMax, sMin);
    }
    
    else if (params->findBlock("PrecipitationMonosurfaceCoefficient") ) {
      double coef, expcoef, molecularDiameter, diffusionCoefficient, surfaceEnergy , T;
      Uintah::ProblemSpecP coefParams = params->findBlock("PrecipitationMonosurfaceCoefficient");
      coefParams -> getAttribute("Molec_D",molecularDiameter);
      coefParams -> getAttribute("Diff_Coef",diffusionCoefficient);
      coefParams -> getAttribute("Surf_Eng", surfaceEnergy);
      coefParams -> getAttribute("Temperature",T);
      coef = diffusionCoefficient * PI / molecularDiameter / molecularDiameter /molecularDiameter;
      expcoef = - surfaceEnergy * surfaceEnergy * molecularDiameter * molecularDiameter * PI / kB / kB / T / T;
      const Expr::Tag saturationTag = parse_nametag( coefParams->findBlock("Supersaturation")->findBlock("NameTag") );
      typedef typename PrecipitationMonosurfaceCoefficient<FieldT>::Builder Builder;
      builder = scinew Builder(tag, saturationTag, coef, expcoef);
    }
    
    else if (params->findBlock("PrecipitationClassicNucleationCoefficient") ) {
      double expcoef, SurfaceEnergy, MolecularVolume, T;
      Uintah::ProblemSpecP coefParams = params->findBlock("PrecipitationClassicNucleationCoefficient");
      coefParams -> getAttribute("Molec_Vol",MolecularVolume);
      coefParams -> getAttribute("Surf_Eng",SurfaceEnergy);
      coefParams -> getAttribute("Temperature",T);
      expcoef = -16 * PI / 3 * SurfaceEnergy * SurfaceEnergy * SurfaceEnergy / kB / kB / kB / T / T / T * MolecularVolume * MolecularVolume / nA / nA;
      const Expr::Tag saturationTag = parse_nametag( coefParams->findBlock("Supersaturation")->findBlock("NameTag") );
      typedef typename PrecipitationClassicNucleationCoefficient<FieldT>::Builder Builder;
      builder = scinew Builder(tag, saturationTag, expcoef);
    }
    
    else if (params->findBlock("HomogeneousNucleationCoefficient") ) {
      double molecularVolume, T, D, sRatio;
      double surfaceEnergy = 1.0;
      Uintah::ProblemSpecP coefParams = params->findBlock("HomogeneousNucleationCoefficient");
      coefParams -> getAttribute("Molar_Vol",molecularVolume);
      if (coefParams->getAttribute("Surf_Eng",surfaceEnergy) )
        coefParams -> getAttribute("Surf_Eng",surfaceEnergy);
      coefParams -> getAttribute("Temperature",T);
      coefParams -> getAttribute("Diff_Coef",D);
      coefParams -> getAttribute("S_Ratio", sRatio);
      molecularVolume = molecularVolume/6.02214129e23; //convert molar to molecular volume in this term
      Expr::Tag surfaceEngTag;
      if ( coefParams->findBlock("SurfaceEnergy") )
        surfaceEngTag = parse_nametag( coefParams->findBlock("SurfaceEnergy")->findBlock("NameTag") );
      const Expr::Tag saturationTag = parse_nametag( coefParams->findBlock("Supersaturation")->findBlock("NameTag") );
      const Expr::Tag eqConcTag = parse_nametag( coefParams->findBlock("EquilibriumConcentration")->findBlock("NameTag") );
      typedef typename HomogeneousNucleationCoefficient<FieldT>::Builder Builder;
      builder = scinew Builder(tag, saturationTag, eqConcTag,  surfaceEngTag, molecularVolume, surfaceEnergy, T, D, sRatio);
    }
    
    else if (params->findBlock("PrecipitationSimpleRStarValue") ) {
      double rknot, coef, cfCoef;
      Uintah::ProblemSpecP coefParams = params->findBlock("PrecipitationSimpleRStarValue");
      coefParams -> getAttribute("R0", rknot);
      coefParams -> getAttribute("Conversion_Fac", cfCoef);
      coef = rknot*cfCoef;
      const Expr::Tag saturationTag = parse_nametag( coefParams->findBlock("Supersaturation")->findBlock("NameTag") );
      const Expr::Tag surfaceEngTag; //dummy tag since this uses same function as classic rStar
      typedef typename PrecipitationRCritical<FieldT>::Builder Builder;
      builder = scinew Builder(tag, saturationTag, surfaceEngTag, coef);
    }
    
    else if (params->findBlock("PrecipitationClassicRStarValue") ) {
      double surfaceEnergy, molecularVolume, T, cfCoef, coef;
      Uintah::ProblemSpecP coefParams = params->findBlock("PrecipitationClassicRStarValue");
      coefParams -> getAttribute("Surf_Eng", surfaceEnergy);
      coefParams -> getAttribute("Conversion_Fac", cfCoef);
      coefParams -> getAttribute("Temperature", T);
      coefParams -> getAttribute("Molec_Vol",molecularVolume);
      coef = 2.0*surfaceEnergy*molecularVolume/R/T*cfCoef;
      Expr::Tag surfaceEngTag;
      if (coefParams->findBlock("SurfaceEnergy") ) {
        surfaceEngTag = parse_nametag( coefParams->findBlock("SurfaceEnergy")->findBlock("NameTag") ) ;
        coef = 2.0*molecularVolume/R/T*cfCoef;  //calcualte coefficient without the surface energy
      }
      
      const Expr::Tag saturationTag = parse_nametag( coefParams->findBlock("Supersaturation")->findBlock("NameTag") );
      typedef typename PrecipitationRCritical<FieldT>::Builder Builder;
      builder = scinew Builder(tag, saturationTag, surfaceEngTag, coef);
      //Note: both RStars are same basic form, same builder, but different coefficient parse
    }
    
    else if(params->findBlock("CriticalSurfaceEnergy") ) {
      double bulkSurfaceEnergy, T, molarVolume, coef, tolmanL;
      Uintah::ProblemSpecP coefParams = params->findBlock("CriticalSurfaceEnergy");
      coefParams -> getAttribute("Temperature",T);
      coefParams -> getAttribute("Bulk_Surf_Eng",bulkSurfaceEnergy);
      coefParams -> getAttribute("Molar_Vol",molarVolume);
      coefParams -> getWithDefault("TolmanLength",tolmanL,0.2);
      const Expr::Tag saturationTag = parse_nametag( coefParams->findBlock("Supersaturation")->findBlock("NameTag") );
      double r1 = pow(3.0*molarVolume/nA/4.0/PI,1.0/3.0); //convert molar vol to molec radius
      coef = 4.0 * tolmanL * R * T*bulkSurfaceEnergy* r1/molarVolume;
      typedef typename CriticalSurfaceEnergy<FieldT>::Builder Builder;
      builder = scinew Builder(tag, saturationTag, bulkSurfaceEnergy, coef);
    }
    
    else if (params->findBlock("BrownianAggregationCoefficient") ) {
      double T, coef;
      Uintah::ProblemSpecP coefParams = params->findBlock("BrownianAggregationCoefficient");
      coefParams -> getAttribute("Temperature", T);
      double ConvFac = 1.0;
      if (coefParams->getAttribute("Conversion_Fac", ConvFac) )
        coefParams->getAttribute("Conversion_Fac", ConvFac);
      coef = 2.0 * kB * T / 3.0 * ConvFac ;
      const Expr::Tag densityTag = parse_nametag( coefParams->findBlock("Density")->findBlock("NameTag") );
      typedef typename BrownianAggregationCoefficient<FieldT>::Builder Builder;
      builder = scinew Builder(tag, densityTag, coef);
    }
    
    else if (params->findBlock("TurbulentAggregationCoefficient") ) {
      Uintah::ProblemSpecP coefParams = params->findBlock("TurbulentAggregationCoefficient");
      const Expr::Tag kinematicViscosityTag = parse_nametag( coefParams->findBlock("KinematicViscosity")->findBlock("NameTag") );
      const Expr::Tag energyDissipationTag = parse_nametag( coefParams->findBlock("EnergyDissipation")->findBlock("NameTag") );
      double coef;
      double convFac = 1.0;
      if (coefParams->getAttribute("Conversion_Fac", convFac) )
        coefParams->getAttribute("Conversion_Fac", convFac);
      coef = (4.0 / 3.0) * sqrt(3.0 * PI / 10.0) * convFac;
      typedef typename TurbulentAggregationCoefficient<FieldT>::Builder Builder;
      builder = scinew Builder(tag, kinematicViscosityTag, energyDissipationTag, coef);
    }
    
    else if (params->findBlock("PrecipitateEffectiveViscosity") ) {
      Uintah::ProblemSpecP coefParams = params->findBlock("PrecipitateEffectiveViscosity");
      double corrFac, power, baseViscos, minStrain;
      const Expr::Tag volFracTag = parse_nametag( coefParams->findBlock("VolumeFraction")->findBlock("NameTag") );
      const Expr::Tag strainMagTag = parse_nametag( coefParams->findBlock("StrainMagnitude")->findBlock("NameTag") );
      coefParams -> getAttribute("Correction_Fac", corrFac);
      coefParams -> getAttribute("Power", power);
      coefParams -> getAttribute("BaseViscosity", baseViscos);
      coefParams -> getAttribute("MinStrain", minStrain);
      typedef typename PrecipitateEffectiveViscosity<FieldT>::Builder Builder;
      builder= scinew Builder(tag, volFracTag, strainMagTag, corrFac, baseViscos, power, minStrain);
    }
    
    else if (params->findBlock("ParticleVolumeFraction") ) {
      Uintah::ProblemSpecP coefParams = params->findBlock("ParticleVolumeFraction");
      double convFac;
      coefParams -> get("Conversion_Fac", convFac);
      std::string basePhiName;
      Expr::Tag momentTag;
      Expr::TagList zerothMomentTags;
      Expr::TagList firstMomentTags;
      //assumes all moment transport eqns are for particles
      for ( Uintah::ProblemSpecP momentParams=wasatchParams->findBlock("MomentTransportEquation");
           momentParams != 0;
           momentParams = momentParams->findNextBlock("MomentTransportEquation") ) {
        momentParams ->get("PopulationName",basePhiName);
        if (momentParams->findBlock("MultiEnvMixingModel") ) {
          momentTag = Expr::Tag("m_" + basePhiName + "_0_ave", Expr::STATE_NONE);
          zerothMomentTags.push_back(momentTag);
          momentTag = Expr::Tag("m_" + basePhiName + "_1_ave", Expr::STATE_NONE);
          firstMomentTags.push_back(momentTag);
        } else {
          momentTag = Expr::Tag("m_" + basePhiName + "_0", Expr::STATE_N);
          zerothMomentTags.push_back(momentTag);
          momentTag = Expr::Tag("m_" + basePhiName + "_1", Expr::STATE_N);
          firstMomentTags.push_back(momentTag);
        }
      }
      typedef typename ParticleVolumeFraction<FieldT>::Builder Builder;
      builder = scinew Builder(tag, zerothMomentTags, firstMomentTags, convFac);
    }
    
    else if (params->findBlock("MultiEnvMixingModel") ) {
      std::stringstream wID;
      std::string baseName, stateType;
      double maxDt;
      Uintah::ProblemSpecP multiEnvParams = params->findBlock("MultiEnvMixingModel");
      multiEnvParams -> getAttribute("MaxDt",maxDt);
      Expr::TagList multiEnvWeightsTags;
      params->findBlock("NameTag")->getAttribute("name",baseName);
      params->findBlock("NameTag")->getAttribute("state",stateType);
      const int numEnv = 3;
      //create the expression for weights and derivatives if block found
      for (int i=0; i<numEnv; i++) {
        wID.str(std::string());
        wID << i;
        if (stateType == "STATE_N") {
          multiEnvWeightsTags.push_back(Expr::Tag("w_" + baseName + "_" + wID.str(), Expr::STATE_N) );
          multiEnvWeightsTags.push_back(Expr::Tag("dwdt_" + baseName + "_" + wID.str(), Expr::STATE_N) );
        } else if (stateType == "STATE_NONE") {
          multiEnvWeightsTags.push_back(Expr::Tag("w_" + baseName + "_" + wID.str(), Expr::STATE_NONE) );
          multiEnvWeightsTags.push_back(Expr::Tag("dwdt_" + baseName + "_" + wID.str(), Expr::STATE_NONE) );
        }
      }
      
      const Expr::Tag mixFracTag = parse_nametag( multiEnvParams->findBlock("MixtureFraction")->findBlock("NameTag") );
      const Expr::Tag scalarVarTag = parse_nametag( multiEnvParams->findBlock("ScalarVariance")->findBlock("NameTag") );
      const Expr::Tag scalarDissTag = parse_nametag( multiEnvParams->findBlock("ScalarDissipation")->findBlock("NameTag") );
      builder = scinew typename MultiEnvMixingModel<FieldT>::Builder(multiEnvWeightsTags, mixFracTag, scalarVarTag, scalarDissTag, maxDt);
    }
    
    else if (params->findBlock("PrecipitationSource") ) {
      //this loops over all possible non-convective/non-diffusive rhs terms and creates a taglist
      Uintah::ProblemSpecP coefParams = params->findBlock("PrecipitationSource");
      std::vector<double> molecVolumes;
      Expr::TagList sourceTagList;
      Expr::Tag sourceTag;
      Expr::Tag midEnvWeightTag; //tag for central weight
      double molecVol;
      std::string modelType, basePhiName;
      
      const Expr::Tag etaScaleTag = parse_nametag( coefParams->findBlock("EtaScale")->findBlock("NameTag") );
      const Expr::Tag densityTag = parse_nametag( coefParams->findBlock("Density")->findBlock("NameTag") );
      
      if (coefParams->findBlock("MultiEnvWeight") ) {
        midEnvWeightTag = parse_nametag( coefParams->findBlock("MultiEnvWeight")->findBlock("NameTag") );
      }
      
      for( Uintah::ProblemSpecP momentParams=wasatchParams->findBlock("MomentTransportEquation");
           momentParams != 0;
           momentParams = momentParams->findNextBlock("MomentTransportEquation") ){
        momentParams->get("MolecVol", molecVol);
        momentParams->get("PopulationName", basePhiName);
        
        for( Uintah::ProblemSpecP growthParams=momentParams->findBlock("GrowthExpression");
             growthParams != 0;
             growthParams = growthParams->findNextBlock("GrowthExpression") ){
          molecVolumes.push_back(molecVol);
          growthParams->get("GrowthModel", modelType);
          sourceTag = Expr::Tag( "m_" + basePhiName + "_3_growth_" + modelType, Expr::STATE_NONE);
          sourceTagList.push_back(sourceTag);
        }
        for( Uintah::ProblemSpecP birthParams=momentParams->findBlock("BirthExpression");
             birthParams != 0;
             birthParams = birthParams->findNextBlock("BirthExpression") ){
          molecVolumes.push_back(molecVol);
          birthParams->get("BirthModel", modelType);
          sourceTag = Expr::Tag("m_" + basePhiName + "_3_birth_" + modelType, Expr::STATE_NONE);
          sourceTagList.push_back(sourceTag);
        }
        for( Uintah::ProblemSpecP deathParams=momentParams->findBlock("Dissolution");
            deathParams != 0;
            deathParams = deathParams->findNextBlock("Dissolution") ){
          molecVolumes.push_back(molecVol);
          sourceTag = Expr::Tag( "m_" + basePhiName + "_3_death", Expr::STATE_NONE);
          sourceTagList.push_back(sourceTag);
        }
      }
      typedef typename PrecipitationSource<FieldT>::Builder Builder;
      builder = scinew Builder(tag, sourceTagList, etaScaleTag, densityTag, midEnvWeightTag, molecVolumes);
    }
    return builder;
  }
  
  //------------------------------------------------------------------
  
  template<typename FieldT>
  Expr::ExpressionBuilder*
  build_post_processing_expr( Uintah::ProblemSpecP params )
  {
    const Expr::Tag tag = parse_nametag( params->findBlock("NameTag") );
    
    Expr::ExpressionBuilder* builder = NULL;
    if( params->findBlock("VelocityMagnitude") ){
      Uintah::ProblemSpecP valParams = params->findBlock("VelocityMagnitude");
      
      Expr::Tag xVelTag = Expr::Tag();
      if (valParams->findBlock("XVelocity"))
        xVelTag = parse_nametag( valParams->findBlock("XVelocity")->findBlock("NameTag") );
      
      Expr::Tag yVelTag = Expr::Tag();
      if (valParams->findBlock("YVelocity"))
        yVelTag = parse_nametag( valParams->findBlock("YVelocity")->findBlock("NameTag") );
      
      Expr::Tag zVelTag = Expr::Tag();
      if (valParams->findBlock("ZVelocity"))
        zVelTag = parse_nametag( valParams->findBlock("ZVelocity")->findBlock("NameTag") );
      
      typedef typename VelocityMagnitude<SVolField, XVolField, YVolField, ZVolField>::Builder Builder;
      builder = scinew Builder(tag, xVelTag, yVelTag, zVelTag);
    }
    
    else if( params->findBlock("Vorticity") ){
      Uintah::ProblemSpecP valParams = params->findBlock("Vorticity");
      std::string vorticityComponent;
      valParams->require("Component",vorticityComponent);
      
      Expr::Tag vel1Tag = Expr::Tag();
      if (valParams->findBlock("Vel1"))
        vel1Tag = parse_nametag( valParams->findBlock("Vel1")->findBlock("NameTag") );
      Expr::Tag vel2Tag = Expr::Tag();
      if (valParams->findBlock("Vel2"))
        vel2Tag = parse_nametag( valParams->findBlock("Vel2")->findBlock("NameTag") );
      if (vorticityComponent == "X") {
        typedef typename Vorticity<SVolField, ZVolField, YVolField>::Builder Builder;
        builder = scinew Builder(tag, vel1Tag, vel2Tag);
      } else if (vorticityComponent == "Y") {
        typedef typename Vorticity<SVolField, XVolField, ZVolField>::Builder Builder;
        builder = scinew Builder(tag, vel1Tag, vel2Tag);
      } else if (vorticityComponent == "Z") {
        typedef typename Vorticity<SVolField, YVolField, XVolField>::Builder Builder;
        builder = scinew Builder(tag, vel1Tag, vel2Tag);
      }
    }
    
    else if( params->findBlock("InterpolateExpression") ){
      Uintah::ProblemSpecP valParams = params->findBlock("InterpolateExpression");
      std::string srcFieldType;
      std::string destFieldType;
      valParams->getAttribute("type",srcFieldType);
      Expr::Tag srcTag = Expr::Tag();
      srcTag = parse_nametag( valParams->findBlock("NameTag") );
      
      switch( get_field_type(srcFieldType) ){
        case SVOL : {
          typedef typename InterpolateExpression<SVolField, FieldT>::Builder Builder;
          builder = scinew Builder(tag, srcTag);
          break;
        }
        case XVOL : {
          typedef typename InterpolateExpression<XVolField, FieldT>::Builder Builder;
          builder = scinew Builder(tag, srcTag);
          break;
        }
        case YVOL : {
          typedef typename InterpolateExpression<YVolField, FieldT>::Builder Builder;
          builder = scinew Builder(tag, srcTag);
          break;
        }
        case ZVOL : {
          typedef typename InterpolateExpression<ZVolField, FieldT>::Builder Builder;
          builder = scinew Builder(tag, srcTag);
          break;
        }
        default:
          std::ostringstream msg;
          msg << "ERROR: unsupported field type '" << srcFieldType << "'" << "while parsing an InterpolateExpression." << std::endl;
          throw Uintah::ProblemSetupException( msg.str(), __FILE__, __LINE__ );
      }
    } else if( params->findBlock("KineticEnergy") ) {
      Uintah::ProblemSpecP keSpec = params->findBlock("KineticEnergy");
      
      Expr::Tag xVelTag = Expr::Tag();
      if (keSpec->findBlock("XVelocity"))
        xVelTag = parse_nametag( keSpec->findBlock("XVelocity")->findBlock("NameTag") );
      
      Expr::Tag yVelTag = Expr::Tag();
      if (keSpec->findBlock("YVelocity"))
        yVelTag = parse_nametag( keSpec->findBlock("YVelocity")->findBlock("NameTag") );
      
      Expr::Tag zVelTag = Expr::Tag();
      if (keSpec->findBlock("ZVelocity"))
        zVelTag = parse_nametag( keSpec->findBlock("ZVelocity")->findBlock("NameTag") );
      
      bool totalKE=false;
      keSpec->getAttribute("total",totalKE);
      if (totalKE) {
        typedef typename TotalKineticEnergy<XVolField, YVolField, ZVolField>::Builder Builder;
        builder = scinew Builder(tag, xVelTag, yVelTag, zVelTag);
      } else {
        typedef typename KineticEnergy<SVolField, XVolField, YVolField, ZVolField>::Builder Builder;
        builder = scinew Builder(tag, xVelTag, yVelTag, zVelTag);
      }      
    }
    
    return builder;
  }
  
  //------------------------------------------------------------------
  
  template<typename FieldT>
  Expr::ExpressionBuilder*
  build_bc_expr( Uintah::ProblemSpecP params,
                 Uintah::ProblemSpecP wasatchSpec )
  {
    const Expr::Tag tag = parse_nametag( params->findBlock("NameTag") );
    
    Expr::ExpressionBuilder* builder = NULL;

    const TagNames& tagNames = TagNames::self();
    
    std::string exprType;
    Uintah::ProblemSpecP valParams = params->get("value",exprType);
    if( params->findBlock("Constant") ){
      double val;  params->get("Constant",val);
      typedef typename ConstantBC<FieldT>::Builder Builder;
      builder = scinew Builder( tag, val );
    }
    
    else if( params->findBlock("LinearFunction") ){
      double slope, intercept;
      Uintah::ProblemSpecP valParams = params->findBlock("LinearFunction");
      valParams->getAttribute("slope",slope);
      valParams->getAttribute("intercept",intercept);
      const Expr::Tag indepVarTag = parse_nametag( valParams->findBlock("NameTag") );
      typedef typename LinearBC<FieldT>::Builder Builder;
      builder = scinew Builder( tag, indepVarTag, slope, intercept );
    }
    
    else if ( params->findBlock("ParabolicFunction") ) {
      double a=0.0, b=0.0, c=0.0, x0=0.0, f0=0.0, h=0.0;
      Uintah::ProblemSpecP valParams = params->findBlock("ParabolicFunction");
      
      const Expr::Tag indepVarTag = parse_nametag( valParams->findBlock("NameTag") );
      
      std::string parabolaType;
      valParams->getAttribute("type", parabolaType);
      
      if (parabolaType.compare("CENTERED") == 0) {
        valParams = valParams->findBlock("Centered");
        valParams->getAttribute("x0",x0);
        valParams->getAttribute("f0",f0);
        valParams->getAttribute("h",h);
        a = - f0/(h*h);
        b = 0.0;
        c = f0;
      } else if (parabolaType.compare("GENERAL") == 0) {
        valParams = valParams->findBlock("General");
        valParams->getAttribute("a",a);
        valParams->getAttribute("b",b);
        valParams->getAttribute("c",c);
      }
      
      typedef typename ParabolicBC<FieldT>::Builder Builder;
      builder = scinew Builder( tag, indepVarTag, a, b, c, x0);
    }
    
    else if ( params->findBlock("PowerLawFunction") ) {
      double x0, phic, R, n;
      Uintah::ProblemSpecP valParams = params->findBlock("PowerLawFunction");
      valParams->getAttribute("x0",x0);
      valParams->getAttribute("PhiCenter",phic);
      valParams->getAttribute("HalfHeight",R);
      valParams->getAttribute("n",n);
      const Expr::Tag indepVarTag = parse_nametag( valParams->findBlock("NameTag") );
      typedef typename PowerLawBC<FieldT>::Builder Builder;
      builder = scinew Builder( tag, indepVarTag,x0, phic, R, n);
    }
    
    else if ( params->findBlock("VarDenMMSVelocity") ){
      std::string side;
      Uintah::ProblemSpecP valParams = params->findBlock("VarDenMMSVelocity");
      valParams->getAttribute("side",side);
      
      typedef VarDen1DMMSVelocity<FieldT> VarDenMMSVExpr;
      SpatialOps::BCSide bcSide;
      if      (side == "PLUS" ) bcSide = SpatialOps::PLUS_SIDE;
      else if (side == "MINUS"  ) bcSide = SpatialOps::MINUS_SIDE;
      else {
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__ << std::endl
        << "ERROR: The boundary side " << side
        << " is not supported in VarDen1DMMSVelocity expression." << std::endl;
        throw std::invalid_argument( msg.str() );
      }
      builder = scinew typename VarDenMMSVExpr::Builder( tag, tagNames.time, bcSide );
    }

    else if ( params->findBlock("VarDenMMSMomentum") ){
      std::string side;
      double rho0=1.29985, rho1=0.081889;
      Uintah::ProblemSpecP valParams = params->findBlock("VarDenMMSMomentum");
      valParams->getAttribute("side",side);
      valParams->get("rho0",rho0);
      valParams->get("rho1",rho1);
      typedef VarDen1DMMSMomentum<FieldT> VarDenMMSMomExpr;
      SpatialOps::BCSide bcSide;
      if      (side == "PLUS" ) bcSide = SpatialOps::PLUS_SIDE;
      else if (side == "MINUS"  ) bcSide = SpatialOps::MINUS_SIDE;
      else {
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__ << std::endl
        << "ERROR: The boundary side " << side
        << " is not supported in VarDen1DMMSMomentum expression." << std::endl;
        throw std::invalid_argument( msg.str() );
      }
      builder = scinew typename VarDenMMSMomExpr::Builder( tag, tagNames.time, rho0, rho1, bcSide );
    }

    else if ( params->findBlock("VarDenMMSMixtureFraction") ){
      Uintah::ProblemSpecP valParams = params->findBlock("VarDenMMSMixtureFraction");      
      typedef VarDen1DMMSMixtureFraction<FieldT> VarDen1DMMSMixtureFractionExpr;
      builder = scinew typename VarDen1DMMSMixtureFractionExpr::Builder( tag, tagNames.time );
    }

    else if ( params->findBlock("VarDenMMSDensity") ){
      double rho0=1.29985, rho1=0.081889;
      Uintah::ProblemSpecP valParams = params->findBlock("VarDenMMSDensity");
      valParams->get("rho0",rho0);
      valParams->get("rho1",rho1);

      typedef VarDen1DMMSDensity<FieldT> VarDen1DMMSDensityExpr;
      builder = scinew typename VarDen1DMMSDensityExpr::Builder( tag, tagNames.time, rho0, rho1 );
    }

    else if ( params->findBlock("VarDenMMSSolnVar") ){
      double rho0=1.29985, rho1=0.081889;
      Uintah::ProblemSpecP valParams = params->findBlock("VarDenMMSSolnVar");
      valParams->get("rho0",rho0);
      valParams->get("rho1",rho1);

      typedef VarDen1DMMSSolnVar<FieldT> VarDen1DMMSSolnVarExpr;
      builder = scinew typename VarDen1DMMSSolnVarExpr::Builder( tag, tagNames.time, rho0, rho1 );
    }
    
    else if ( params->findBlock("TurbulentInlet") ) {
      std::string inputFileName;
      std::string velDir;
      int period=1;
      double timePeriod;
      Uintah::ProblemSpecP valParams = params->findBlock("TurbulentInlet");
      valParams->get("InputFile",inputFileName);
      valParams->getAttribute("component",velDir);
      
      bool hasPeriod = valParams->getAttribute("period",period);
      bool hasTimePeriod = valParams->getAttribute("timeperiod",timePeriod);
      if (hasTimePeriod) period = 0;
      
      if (hasPeriod && hasTimePeriod) {
        std::ostringstream msg;
        msg << "ERROR: When specifying a TurbulentInletBC, you cannot specify both timeperiod AND period. Please revise your input file." << std::endl;
        throw Uintah::ProblemSetupException( msg.str(), __FILE__, __LINE__ );
      }
      
      typedef typename TurbulentInletBC<FieldT>::Builder Builder;
      builder = scinew Builder(tag,inputFileName, velDir,period, timePeriod);
    }
    
    if (wasatchSpec->findBlock("VarDenCorrugatedMMS")) {
      Uintah::ProblemSpecP varDen2DMMSParams = wasatchSpec->findBlock("VarDenCorrugatedMMS");
      const TagNames& tagNames = TagNames::self();
      double rho0, rho1, d, w, a, b, k, uf, vf;
      varDen2DMMSParams->getAttribute("rho0",rho0);
      varDen2DMMSParams->getAttribute("rho1",rho1);
      varDen2DMMSParams->getAttribute("uf",uf);
      varDen2DMMSParams->getAttribute("vf",vf);
      varDen2DMMSParams->getAttribute("k",k);
      varDen2DMMSParams->getAttribute("w",w);
      varDen2DMMSParams->getAttribute("d",d);
      varDen2DMMSParams->getAttribute("a",a);
      varDen2DMMSParams->getAttribute("b",b);

      if ( params->findBlock("VDCorMMSVelocity") ) {
        typedef VarDenCorrugatedMMSVelocityBC<FieldT> VDCorMMSVelocityT;
        builder = scinew typename VDCorMMSVelocityT::Builder(tag, tagNames.xxvolcoord, tagNames.yxvolcoord, tagNames.time, rho0, rho1, d, w, k, a, b, uf, vf);
      } else if ( params->findBlock("VDCorMMSMomentum") ) {
        builder = scinew typename VarDenCorrugatedMMSMomBC<FieldT>::Builder(tag, tagNames.xxvolcoord, tagNames.yxvolcoord, tagNames.time, rho0, rho1, d, w, k, a, b, uf, vf);
      } else if ( params->findBlock("VDCorMMSMixtureFraction") ) {
        builder = scinew typename VarDenCorrugatedMMSMixFracBC<FieldT>::Builder(tag, tagNames.xsvolcoord, tagNames.ysvolcoord, tagNames.time, rho0, rho1, d, w, k, a, b, uf, vf);
      } else if ( params->findBlock("VDCorMMSRhof") ) {
        builder = scinew typename VarDenCorrugatedMMSRhofBC<FieldT>::Builder(tag, tagNames.xsvolcoord, tagNames.ysvolcoord, tagNames.time, rho0, rho1, d, w, k, a, b, uf, vf);
      } else if ( params->findBlock("VDCorMMSYMomentum") ) {
        builder = scinew typename VarDenCorrugatedMMSyMomBC<FieldT>::Builder(tag, tagNames.xyvolcoord, tagNames.yyvolcoord, tagNames.time, rho0, rho1, d, w, k, a, b, uf, vf);
      } else if ( params->findBlock("VDCorMMSRho") ) {
        builder = scinew typename VarDenCorrugatedMMSRho<FieldT>::Builder(tag, tagNames.xsvolcoord, tagNames.ysvolcoord, tagNames.time, rho0, rho1, d, w, k, a, b, uf, vf);
      }
    }
    
    return builder;
  }
  
  //------------------------------------------------------------------
  
  void
  create_expressions_from_input( Uintah::ProblemSpecP parser,
                                GraphCategories& gc )
  {
    Expr::ExpressionBuilder* builder = NULL;
    
    //___________________________________
    // parse and build basic expressions
    for( Uintah::ProblemSpecP exprParams = parser->findBlock("BasicExpression");
        exprParams != 0;
        exprParams = exprParams->findNextBlock("BasicExpression") ){
      
      std::string fieldType;
      exprParams->getAttribute("type",fieldType);
      
      switch( get_field_type(fieldType) ){
        case SVOL : builder = build_basic_expr< SVolField >( exprParams );  break;
        case XVOL : builder = build_basic_expr< XVolField >( exprParams );  break;
        case YVOL : builder = build_basic_expr< YVolField >( exprParams );  break;
        case ZVOL : builder = build_basic_expr< ZVolField >( exprParams );  break;
        case PARTICLE : builder = build_basic_particle_expr< ParticleField >( exprParams );  break;
        default:
          std::ostringstream msg;
          msg << "ERROR: unsupported field type '" << fieldType << "'" << std::endl;
          throw Uintah::ProblemSetupException( msg.str(), __FILE__, __LINE__ );
      }
      
      const Category cat = parse_tasklist(exprParams,false);
      gc[cat]->exprFactory->register_expression( builder );
    }
    
    //________________________________________
    // parse and build Taylor-Green Vortex MMS
    for( Uintah::ProblemSpecP exprParams = parser->findBlock("TaylorVortexMMS");
        exprParams != 0;
        exprParams = exprParams->findNextBlock("TaylorVortexMMS") ){
      
      std::string fieldType;
      exprParams->getAttribute("type",fieldType);
      
      switch( get_field_type(fieldType) ){
        case SVOL : builder = build_taylor_vortex_mms_expr< SVolField >( exprParams );  break;
        case XVOL : builder = build_taylor_vortex_mms_expr< XVolField >( exprParams );  break;
        case YVOL : builder = build_taylor_vortex_mms_expr< YVolField >( exprParams );  break;
        case ZVOL : builder = build_taylor_vortex_mms_expr< ZVolField >( exprParams );  break;
        default:
          std::ostringstream msg;
          msg << "ERROR: unsupported field type '" << fieldType << "'" << std::endl;
          throw Uintah::ProblemSetupException( msg.str(), __FILE__, __LINE__ );
      }
      
      const Category cat = parse_tasklist(exprParams,false);
      gc[cat]->exprFactory->register_expression( builder );
    }
    
    //___________________________________________________
    // parse and build physical coefficients expressions
    for( Uintah::ProblemSpecP exprParams = parser->findBlock("PrecipitationBasicExpression");
        exprParams != 0;
        exprParams = exprParams->findNextBlock("PrecipitationBasicExpression") ){
      
      std::string fieldType;
      exprParams->getAttribute("type",fieldType);
      
      switch( get_field_type(fieldType) ){
        case SVOL : builder = build_precipitation_expr< SVolField >( exprParams , parser);  break;
        case XVOL : builder = build_precipitation_expr< XVolField >( exprParams , parser);  break;
        case YVOL : builder = build_precipitation_expr< YVolField >( exprParams , parser);  break;
        case ZVOL : builder = build_precipitation_expr< ZVolField >( exprParams , parser);  break;
        default:
          std::ostringstream msg;
          msg << "ERROR: unsupported field type '" << fieldType << "'" << std::endl;
          throw Uintah::ProblemSetupException( msg.str(), __FILE__, __LINE__ );
      }
      
      const Category cat = parse_tasklist(exprParams,false);
      gc[cat]->exprFactory->register_expression( builder );
    }
    
    //___________________________________________________
    // parse and build post-processing expressions
    for( Uintah::ProblemSpecP exprParams = parser->findBlock("PostProcessingExpression");
        exprParams != 0;
        exprParams = exprParams->findNextBlock("PostProcessingExpression") ){
      
      std::string fieldType;
      exprParams->getAttribute("type",fieldType);
      
      switch( get_field_type(fieldType) ){
        case SVOL : builder = build_post_processing_expr< SVolField >( exprParams );  break;
        case XVOL : builder = build_post_processing_expr< XVolField >( exprParams );  break;
        case YVOL : builder = build_post_processing_expr< YVolField >( exprParams );  break;
        case ZVOL : builder = build_post_processing_expr< ZVolField >( exprParams );  break;
        default:
          std::ostringstream msg;
          msg << "ERROR: unsupported field type '" << fieldType << "'. Postprocessing expressions are setup with SVOLFields as destination fields only." << std::endl
          << "You were trying to register a postprocessing expression with a non cell centered destination field. Please revise you input file." << std::endl;
          throw Uintah::ProblemSetupException( msg.str(), __FILE__, __LINE__ );
      }

      const Category cat = parse_tasklist(exprParams,false);
      gc[cat]->exprFactory->register_expression( builder );
    }
    
    //___________________________________________________
    // parse and build boundary condition expressions
    for( Uintah::ProblemSpecP exprParams = parser->findBlock("BCExpression");
        exprParams != 0;
        exprParams = exprParams->findNextBlock("BCExpression") ) {
      
      std::string fieldType;
      exprParams->getAttribute("type",fieldType);
      
      // get the list of tasks
      std::string taskNames;
      exprParams->require("TaskList", taskNames);
      std::stringstream ss(taskNames);
      std::istream_iterator<std::string> begin(ss);
      std::istream_iterator<std::string> end;
      std::vector<std::string> taskNamesList(begin,end);
      std::vector<std::string>::iterator taskNameIter = taskNamesList.begin();
      
      // iterate through the list of tasks to which this expression is to be added
      while (taskNameIter != taskNamesList.end()) {
        std::string taskName = *taskNameIter;
        
        switch( get_field_type(fieldType) ){
          case SVOL : builder = build_bc_expr< SVolField >( exprParams, parser );  break;
          case XVOL : builder = build_bc_expr< XVolField >( exprParams, parser );  break;
          case YVOL : builder = build_bc_expr< YVolField >( exprParams, parser );  break;
          case ZVOL : builder = build_bc_expr< ZVolField >( exprParams, parser );  break;
          default:
            std::ostringstream msg;
            msg << "ERROR: unsupported field type '" << fieldType << "' while trying to register BC expression.." << std::endl;
            throw Uintah::ProblemSetupException( msg.str(), __FILE__, __LINE__ );
        }
        
        Category cat = INITIALIZATION;
        if     ( taskName == "initialization"   )   cat = INITIALIZATION;
        else if( taskName == "advance_solution" )   cat = ADVANCE_SOLUTION;
        else{
          std::ostringstream msg;
          msg << "ERROR: unsupported task list '" << taskName << "' while parsing BCExpression." << std::endl;
          throw Uintah::ProblemSetupException( msg.str(), __FILE__, __LINE__ );
        }
        
        GraphHelper* const graphHelper = gc[cat];
        graphHelper->exprFactory->register_expression( builder );
        
        ++taskNameIter;
      }
    }
    
    // This is a special parser for turbulent inlets
    for( Uintah::ProblemSpecP exprParams = parser->findBlock("TurbulentInlet");
        exprParams != 0;
        exprParams = exprParams->findNextBlock("TurbulentInlet") ) {
      
      std::string inputFileName, velDir, baseName;
      int period=1;
      double timePeriod;
      exprParams->get("InputFile",inputFileName);
      exprParams->get("BaseName",baseName);
      
      bool hasPeriod = exprParams->getAttribute("period",period);
      bool hasTimePeriod = exprParams->getAttribute("timeperiod",timePeriod);
      if (hasTimePeriod) period = 0;
      
      if (hasPeriod && hasTimePeriod) {
        std::ostringstream msg;
        msg << "ERROR: When specifying a TurbulentInletBC, you cannot specify both timeperiod AND period. Please revise your input file." << std::endl;
        throw Uintah::ProblemSetupException( msg.str(), __FILE__, __LINE__ );
      }
      
      Expr::Tag xVelTag("x-" + baseName, Expr::STATE_NONE);
      typedef TurbulentInletBC<XVolField>::Builder xBuilder;
      
      Expr::Tag yVelTag("y-" + baseName, Expr::STATE_NONE);
      typedef TurbulentInletBC<YVolField>::Builder yBuilder;
      
      Expr::Tag zVelTag("z-" + baseName, Expr::STATE_NONE);
      typedef TurbulentInletBC<ZVolField>::Builder zBuilder;
      
      GraphHelper* const initGraphHelper = gc[INITIALIZATION];
      initGraphHelper->exprFactory->register_expression( scinew xBuilder(xVelTag, inputFileName, "X", period, timePeriod) );
      initGraphHelper->exprFactory->register_expression( scinew yBuilder(yVelTag, inputFileName, "Y", period, timePeriod) );
      initGraphHelper->exprFactory->register_expression( scinew zBuilder(zVelTag, inputFileName, "Z", period, timePeriod) );
      
      GraphHelper* const slnGraphHelper = gc[ADVANCE_SOLUTION];
      slnGraphHelper->exprFactory->register_expression( scinew xBuilder(xVelTag, inputFileName, "X", period, timePeriod) );
      slnGraphHelper->exprFactory->register_expression( scinew yBuilder(yVelTag, inputFileName, "Y", period, timePeriod) );
      slnGraphHelper->exprFactory->register_expression( scinew zBuilder(zVelTag, inputFileName, "Z", period, timePeriod) );
    }

    //_________________________________________________
    // This is a special parser for variable density MMS
    for( Uintah::ProblemSpecP exprParams = parser->findBlock("VarDenOscillatingMMS");
        exprParams != 0;
        exprParams = exprParams->findNextBlock("VarDenOscillatingMMS") ) {
      
      const TagNames& tagNames = TagNames::self();

      double rho0, rho1, uf, vf, k, w, d;
      exprParams->getAttribute("rho0",rho0);
      exprParams->getAttribute("rho1",rho1);
      exprParams->getAttribute("uf",uf);
      exprParams->getAttribute("vf",vf);
      exprParams->getAttribute("k",k);
      exprParams->getAttribute("w",w);
      exprParams->getAttribute("d",d);
      
      std::string x1, x2;
      exprParams->getWithDefault("x1",x1,"X");
      exprParams->getWithDefault("x2",x2,"Y");

      Expr::Tag x1Tag, x2Tag;

      if      (x1 == "X")  x1Tag = tagNames.xsvolcoord;
      else if (x1 == "Y")  x1Tag = tagNames.ysvolcoord;
      else if (x1 == "Z")  x1Tag = tagNames.zsvolcoord;

      if      (x2 == "X")  x2Tag = tagNames.xsvolcoord;
      else if (x2 == "Y")  x2Tag = tagNames.ysvolcoord;
      else if (x2 == "Z")  x2Tag = tagNames.zsvolcoord;
      
      GraphHelper* const initGraphHelper = gc[INITIALIZATION];

      std::string mixFracName;
      exprParams->get("Scalar", mixFracName);
      const Expr::Tag mixFracTag( mixFracName, Expr::STATE_NONE );
      typedef VarDenOscillatingMMSMixFrac<SVolField>::Builder MixFracBuilder;
      initGraphHelper->exprFactory->register_expression( scinew MixFracBuilder( mixFracTag, x1Tag, x2Tag, tagNames.time, rho0, rho1, w, k, uf, vf ) );

      const Expr::Tag diffCoefTag = parse_nametag(exprParams->findBlock("DiffusionCoefficient")->findBlock("NameTag"));
      const Expr::Tag densityTag = parse_nametag( parser->findBlock("Density")->findBlock("NameTag") );
      typedef DiffusiveConstant<SVolField>::Builder diffCoefBuilder;
      gc[ADVANCE_SOLUTION]->exprFactory->register_expression( scinew diffCoefBuilder( diffCoefTag, densityTag, d ) );
    }
    
    //_________________________________________________
    // This is a special parser for variable density MMS
    for( Uintah::ProblemSpecP exprParams = parser->findBlock("VarDenCorrugatedMMS");
        exprParams != 0;
        exprParams = exprParams->findNextBlock("VarDenCorrugatedMMS") ) {
      
      const TagNames& tagNames = TagNames::self();
      
      double rho0, rho1, uf, vf, k, w, d, a, b;
      exprParams->getAttribute("rho0",rho0);
      exprParams->getAttribute("rho1",rho1);
      exprParams->getAttribute("uf",uf);
      exprParams->getAttribute("vf",vf);
      exprParams->getAttribute("k",k);
      exprParams->getAttribute("w",w);
      exprParams->getAttribute("d",d);
      exprParams->getAttribute("a",a);
      exprParams->getAttribute("b",b);
      
      std::string x1, x2;
      exprParams->getWithDefault("x1",x1,"X");
      exprParams->getWithDefault("x2",x2,"Y");
      
      Expr::Tag x1Tag, x2Tag;
      
      if      (x1 == "X")  x1Tag = tagNames.xsvolcoord;
      else if (x1 == "Y")  x1Tag = tagNames.ysvolcoord;
      else if (x1 == "Z")  x1Tag = tagNames.zsvolcoord;
      
      if      (x2 == "X")  x2Tag = tagNames.xsvolcoord;
      else if (x2 == "Y")  x2Tag = tagNames.ysvolcoord;
      else if (x2 == "Z")  x2Tag = tagNames.zsvolcoord;
      
      GraphHelper* const initGraphHelper = gc[INITIALIZATION];
      
      std::string mixFracName;
      exprParams->get("Scalar", mixFracName);
      const Expr::Tag mixFracTag( mixFracName, Expr::STATE_NONE );
      typedef VarDenCorrugatedMMSMixFrac<SVolField>::Builder MixFracBuilder;
      initGraphHelper->exprFactory->register_expression( scinew MixFracBuilder( mixFracTag, x1Tag, x2Tag, tagNames.time, rho0, rho1, d, w, k, a, b, uf, vf ) );
      
      const Expr::Tag diffCoefTag = parse_nametag(exprParams->findBlock("DiffusionCoefficient")->findBlock("NameTag"));
      const Expr::Tag densityTag = parse_nametag( parser->findBlock("Density")->findBlock("NameTag") );
      typedef DiffusiveConstant<SVolField>::Builder diffCoefBuilder;
      gc[ADVANCE_SOLUTION]->exprFactory->register_expression( scinew diffCoefBuilder( diffCoefTag, densityTag, d ) );
    }

    
    //___________________________________________________
    // parse and build initial conditions for moment transport
    int nEnv = 0;
    if (parser->findBlock("MomentTransportEquation")) {
      parser->findBlock("MomentTransportEquation")->get( "NumberOfEnvironments", nEnv );
    }
    int nEqs = 2*nEnv; // we need the number of equations so that we only build the necessary number of moments for initialization
    for( Uintah::ProblemSpecP exprParams = parser->findBlock("MomentInitialization");
        exprParams != 0;
        exprParams = exprParams->findNextBlock("MomentInitialization") ){
      
      std::string populationName;
      
      exprParams->get("PopulationName", populationName);
      std::vector<double> initialMoments;
      
      exprParams->get("Values", initialMoments, nEqs); // get only the first nEqs moments
      
      assert( (int) initialMoments.size() == nEqs );
      
      const int nMoments = initialMoments.size();
      typedef Expr::ConstantExpr<SVolField>::Builder Builder;
      GraphHelper* const graphHelper = gc[INITIALIZATION];
      for (int i=0; i<nMoments; i++) {
        double val = initialMoments[i];
        std::stringstream ss;
        ss << i;
        Expr::Tag thisMomentTag("m_" + populationName + "_" + ss.str(), Expr::STATE_DYNAMIC);
        graphHelper->exprFactory->register_expression( scinew Builder( thisMomentTag, val ) );
      }
    }
  }
  
  //------------------------------------------------------------------
  
}
