
/*
 *  GuageWidget.h
 *
 *  Written by:
 *   James Purciful
 *   Department of Computer Science
 *   University of Utah
 *   Aug. 1994
 *
 *  Copyright (C) 1994 SCI Group
 */


#include <Widgets/GuageWidget.h>
#include <Constraints/DistanceConstraint.h>
#include <Constraints/SegmentConstraint.h>
#include <Constraints/RatioConstraint.h>
#include <Geom/Cylinder.h>
#include <Geom/Sphere.h>

const Index NumCons = 4;
const Index NumVars = 6;
const Index NumGeoms = 6;
const Index NumMatls = 4;
const Index NumSchemes = 3;
const Index NumPcks = 6;

enum { ConstLine, ConstDist, ConstSDist, ConstRatio };
enum { GeomPointL, GeomPointR, GeomShaft, GeomSlider,
       GeomResizeL, GeomResizeR };
enum { PickSphL, PickSphR, PickCyl, PickSlider,
       PickResizeL, PickResizeR };

GuageWidget::GuageWidget( Module* module, CrowdMonitor* lock, double widget_scale )
: BaseWidget(module, lock, NumVars, NumCons, NumGeoms, NumMatls, NumPcks, widget_scale),
  oldaxis(1, 0, 0)
{
   Real INIT = 1.0*widget_scale;
   // Scheme3 is for resizing.
   variables[PointLVar] = new PointVariable("PntL", solve, Scheme1, Point(0, 0, 0));
   variables[PointRVar] = new PointVariable("PntR", solve, Scheme1, Point(INIT, 0, 0));
   variables[SliderVar] = new PointVariable("Slider", solve, Scheme2, Point(INIT/2.0, 0, 0));
   variables[DistVar] = new RealVariable("Dist", solve, Scheme1, INIT);
   variables[SDistVar] = new RealVariable("SDistVar", solve, Scheme2, INIT/2.0);
   variables[RatioVar] = new RealVariable("Ratio", solve, Scheme1, 0.5);
   
   constraints[ConstLine] = new SegmentConstraint("ConstLine",
						  NumSchemes,
						  variables[PointLVar],
						  variables[PointRVar],
						  variables[SliderVar]);
   constraints[ConstLine]->VarChoices(Scheme1, 2, 2, 2);
   constraints[ConstLine]->VarChoices(Scheme2, 2, 2, 2);
   constraints[ConstLine]->VarChoices(Scheme3, 2, 2, 2);
   constraints[ConstLine]->Priorities(P_Default, P_Default, P_Highest);
   constraints[ConstDist] = new DistanceConstraint("ConstDist",
						   NumSchemes,
						   variables[PointLVar],
						   variables[PointRVar],
						   variables[DistVar]);
   constraints[ConstDist]->VarChoices(Scheme1, 1, 0, 1);
   constraints[ConstDist]->VarChoices(Scheme2, 1, 0, 1);
   constraints[ConstDist]->VarChoices(Scheme3, 2, 2, 2);
   constraints[ConstDist]->Priorities(P_Highest, P_Highest, P_Default);
   constraints[ConstSDist] = new DistanceConstraint("ConstSDist",
						       NumSchemes,
						       variables[PointLVar],
						       variables[SliderVar],
						       variables[SDistVar]);
   constraints[ConstSDist]->VarChoices(Scheme1, 1, 1, 1);
   constraints[ConstSDist]->VarChoices(Scheme2, 2, 2, 2);
   constraints[ConstSDist]->VarChoices(Scheme3, 1, 1, 1);
   constraints[ConstSDist]->Priorities(P_Lowest, P_Default, P_Default);
   constraints[ConstRatio] = new RatioConstraint("ConstRatio",
						 NumSchemes,
						 variables[SDistVar],
						 variables[DistVar],
						 variables[RatioVar]);
   constraints[ConstRatio]->VarChoices(Scheme1, 0, 0, 0);
   constraints[ConstRatio]->VarChoices(Scheme2, 2, 2, 2);
   constraints[ConstRatio]->VarChoices(Scheme3, 0, 0, 0);
   constraints[ConstRatio]->Priorities(P_Highest, P_Highest, P_Highest);

   materials[PointMatl] = PointWidgetMaterial;
   materials[EdgeMatl] = EdgeWidgetMaterial;
   materials[SliderMatl] = SliderWidgetMaterial;
   materials[HighMatl] = HighlightWidgetMaterial;

   geometries[GeomPointL] = new GeomSphere;
   picks[PickSphL] = new GeomPick(geometries[GeomPointL], module);
   picks[PickSphL]->set_highlight(materials[HighMatl]);
   picks[PickSphL]->set_cbdata((void*)PickSphL);
   GeomMaterial* sphlm = new GeomMaterial(picks[PickSphL], materials[PointMatl]);
   geometries[GeomPointR] = new GeomSphere;
   picks[PickSphR] = new GeomPick(geometries[GeomPointR], module);
   picks[PickSphR]->set_highlight(materials[HighMatl]);
   picks[PickSphR]->set_cbdata((void*)PickSphR);
   GeomMaterial* sphrm = new GeomMaterial(picks[PickSphR], materials[PointMatl]);
   GeomGroup* resizes = new GeomGroup;
   geometries[GeomResizeL] = new GeomCappedCylinder;
   picks[PickResizeL] = new GeomPick(geometries[GeomResizeL], module);
   picks[PickResizeL]->set_highlight(materials[HighMatl]);
   picks[PickResizeL]->set_cbdata((void*)PickResizeL);
   resizes->add(picks[PickResizeL]);
   geometries[GeomResizeR] = new GeomCappedCylinder;
   picks[PickResizeR] = new GeomPick(geometries[GeomResizeR], module);
   picks[PickResizeR]->set_highlight(materials[HighMatl]);
   picks[PickResizeR]->set_cbdata((void*)PickResizeR);
   resizes->add(picks[PickResizeR]);
   GeomMaterial* resizesm = new GeomMaterial(resizes, SpecialWidgetMaterial);
   geometries[GeomShaft] = new GeomCylinder;
   picks[PickCyl] = new GeomPick(geometries[GeomShaft], module);
   picks[PickCyl]->set_highlight(materials[HighMatl]);
   picks[PickCyl]->set_cbdata((void*)PickCyl);
   GeomMaterial* cylm = new GeomMaterial(picks[PickCyl], materials[EdgeMatl]);
   geometries[GeomSlider] = new GeomCappedCylinder;
   picks[PickSlider] = new GeomPick(geometries[GeomSlider], module);
   picks[PickSlider]->set_highlight(materials[HighMatl]);
   picks[PickSlider]->set_cbdata((void*)PickSlider);
   GeomMaterial* sliderm = new GeomMaterial(picks[PickSlider], materials[SliderMatl]);

   GeomGroup* w = new GeomGroup;
   w->add(sphlm);
   w->add(sphrm);
   w->add(resizesm);
   w->add(cylm);
   w->add(sliderm);

   SetEpsilon(widget_scale*1e-6);

   FinishWidget(w);
}


GuageWidget::~GuageWidget()
{
}


void
GuageWidget::widget_execute()
{
   ((GeomSphere*)geometries[GeomPointL])->move(variables[PointLVar]->point(),
						  1*widget_scale);
   ((GeomSphere*)geometries[GeomPointR])->move(variables[PointRVar]->point(),
						  1*widget_scale);
   ((GeomCylinder*)geometries[GeomShaft])->move(variables[PointLVar]->point(),
						variables[PointRVar]->point(),
						0.5*widget_scale);
   ((GeomCappedCylinder*)geometries[GeomResizeL])->move(variables[PointLVar]->point(),
							variables[PointLVar]->point()
							- (GetAxis() * 1.5 * widget_scale),
							0.5*widget_scale);
   ((GeomCappedCylinder*)geometries[GeomResizeR])->move(variables[PointRVar]->point(),
							variables[PointRVar]->point()
							+ (GetAxis() * 1.5 * widget_scale),
							0.5*widget_scale);
   ((GeomCappedCylinder*)geometries[GeomSlider])->move(variables[SliderVar]->point()
						       - (GetAxis() * 0.3 * widget_scale),
						       variables[SliderVar]->point()
						       + (GetAxis() * 0.3 * widget_scale),
						       1.1*widget_scale);

   SetEpsilon(widget_scale*1e-6);

   Vector v(GetAxis()), v1, v2;
   v.find_orthogonal(v1,v2);
   for (Index geom = 0; geom < NumPcks; geom++) {
      if ((geom == PickSlider) || (geom == PickResizeL) || (geom == PickResizeR))
	 picks[geom]->set_principal(v);
      else
	 picks[geom]->set_principal(v, v1, v2);
   }
}


void
GuageWidget::geom_moved( int /* axis */, double /* dist */, const Vector& delta,
			 void* cbdata )
{
   ((DistanceConstraint*)constraints[ConstSDist])->SetDefault(GetAxis());
   ((DistanceConstraint*)constraints[ConstDist])->SetDefault(GetAxis());

   switch((int)cbdata){
   case PickSphL:
      variables[PointLVar]->SetDelta(delta);
      break;
   case PickSphR:
      variables[PointRVar]->SetDelta(delta);
      break;
   case PickResizeL:
      variables[PointLVar]->SetDelta(delta, Scheme3);
      break;
   case PickResizeR:
      variables[PointRVar]->SetDelta(delta, Scheme3);
      break;
   case PickSlider:
      variables[SliderVar]->SetDelta(delta);
      break;
   case PickCyl:
      MoveDelta(delta);
      break;
   }
}


void
GuageWidget::MoveDelta( const Vector& delta )
{
   variables[PointLVar]->MoveDelta(delta);
   variables[PointRVar]->MoveDelta(delta);
   variables[SliderVar]->MoveDelta(delta);
}


Point
GuageWidget::ReferencePoint() const
{
   return (variables[PointLVar]->point()
	   + (variables[PointRVar]->point()
	      -variables[PointLVar]->point())/2.0);
}


void
GuageWidget::SetRatio( const Real ratio )
{
   ASSERT((ratio>=0.0) && (ratio<=1.0));
   variables[RatioVar]->Set(ratio);
   
   execute();
}


Real
GuageWidget::GetRatio() const
{
   return (variables[RatioVar]->real());
}


void
GuageWidget::SetEndpoints( const Point& end1, const Point& end2 )
{
   variables[PointLVar]->Move(end1);
   variables[PointRVar]->Set(end2);
   
   execute();
}


void
GuageWidget::GetEndpoints( Point& end1, Point& end2 ) const
{
   end1 = variables[PointLVar]->point();
   end2 = variables[PointRVar]->point();
}


const Vector&
GuageWidget::GetAxis()
{
   Vector axis(variables[PointRVar]->point() - variables[PointLVar]->point());
   if (axis.length2() <= 1e-6)
      return oldaxis;
   else 
      return (oldaxis = axis.normal());
}


