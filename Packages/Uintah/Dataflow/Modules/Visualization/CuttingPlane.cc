//static char *id="@(#) $Id$";

/*
 *  CuttingPlane.cc:  
 *
 *  Written by:
 *   Colette Mullenhoff
 *   Department of Computer Science
 *   University of Utah
 *   May 1995
 *
 *  Copyright (C) 1995 SCI Group
 */

#include <Core/Containers/Array1.h>
#include <Dataflow/Network/Module.h>
#include <Dataflow/Ports/GeometryPort.h>
#include <Dataflow/Ports/FieldPort.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldAlgo.h>
#include <Dataflow/Ports/ColorMapPort.h>
#include <Core/Geom/GeomGrid.h>
#include <Core/Geom/GeomGroup.h>
#include <Core/Geom/GeomLine.h>
#include <Core/Geom/Material.h>
#include <Core/Geometry/Point.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/MiscMath.h>
#include <Core/Malloc/Allocator.h>
#include <Core/GuiInterface/GuiVar.h>
#include <Core/Thread/CrowdMonitor.h>
#include <Core/Datatypes/LatticeVol.h>
#include <Dataflow/Widgets/ScaledFrameWidget.h>
#include <Packages/Uintah/Core/Datatypes/DispatchScalar1.h>
#include <Packages/Uintah/Core/Datatypes/LevelField.h>
#include <iostream>
using std::cerr;
using std::cin;
using std::endl;

#define CP_PLANE 0
#define CP_SURFACE 1
#define CP_CONTOUR 2

using namespace SCIRun;

namespace Uintah {

#define EPS 1.e-6

/**************************************
CLASS
   CuttingPlane
       CuttingPlane interpolates a planar slice through the
       unstructured data and maps data values to colors on a
       semitransparent surface.

GENERAL INFORMATION

   CuttingPlane
  
   Author:  Colette Mullenhoff<br>
            Department of Computer Science<br>
            University of Utah
   Date:    June 1995
   
   C-SAFE
   
   Copyright <C> 1995 SCI Group

KEYWORDS
   Visualization, Widget, GenColormap

DESCRIPTION
   CuttingPlane interpolates a planar slice through the
   unstructured data and maps data values to colors on a
   semitransparent surface.  The plane can be manipulated with a
   3D widget to allow the user to look at different
   cross sections of the electric potential.

WARNING
   None

****************************************/

class CuttingPlane : public Module {
   FieldIPort *infield;
   ColorMapIPort *inColorMap;
   GeometryOPort* ogeom;
   CrowdMonitor widget_lock;
   int init;
   int widget_id;
   ScaledFrameWidget *widget;
   virtual void widget_moved(int last);
   GuiInt cutting_plane_type;
   GuiInt num_contours;   
   GuiDouble offset;
   GuiDouble scale;
   GuiDouble where;
   GuiInt need_find;
   GuiInt localMinMaxGUI;
   GuiInt fullRezGUI;
   GuiInt exhaustiveGUI;
   MaterialHandle outcolor;
   int grid_id;
   string msg;
public:
 
        // GROUP:  Constructors:
        ///////////////////////////
        //
        // Constructs an instance of class CuttingPlane
        //
        // Constructor taking
        //    [in] id as an identifier
        //
   CuttingPlane(const string& id);

        // GROUP:  Destructor:
        ///////////////////////////
        // Destructor
   virtual ~CuttingPlane();

        // GROUP:  Access functions:
        ///////////////////////////
        //
        // execute() - execution scheduled by scheduler
   virtual void execute();
  
  bool interpolate(FieldHandle f, const Point& p, double& val);

  bool get_gradient(FieldHandle f, const Point& p, Vector& g);

  void get_minmax(FieldHandle f);
  template<class F> void dispatch_minmax(F *f) 
  { success_ = field_minmax(*f, minmax_); }

  bool get_dimensions(FieldHandle f, int& nx, int& ny, int& nz);
  template <class Mesh>  bool get_dimensions(Mesh, int&, int&, int&);
        //////////////////////////
        //
        // tcl_commands - overides tcl_command in base class Module, takes:
        //                                  findxy, findyz, findxz, 
        //                                  plusx, minusx, plusy, minusy, plusz, minusz,
        //                                  connectivity
   virtual void tcl_command(TCLArgs&, void*);
private:
  Point iPoint_;
  double result_;
  Vector gradient_;
  bool success_;
  IntVector dim_;
  pair<double, double> minmax_;
};

extern "C" Module* make_CuttingPlane(const string& id) {
  return new CuttingPlane(id);
}

//static string module_name("CuttingPlane");
static string widget_name("CuttingPlane Widget");

CuttingPlane::CuttingPlane(const string& id) :
  Module("CuttingPlane", id, Filter),
  widget_lock("Cutting plane widget lock"),
  cutting_plane_type("cutting_plane_type",id, this),
  num_contours("num_contours", id, this), 
  offset("offset", id, this), scale("scale", id, this), 
  where("where", id, this),
  need_find("need_find",id,this),
  localMinMaxGUI("localMinMaxGUI", id, this), 
  fullRezGUI("fullRezGUI", id, this),
  exhaustiveGUI("exhaustiveGUI", id, this)
{
    // Create the input ports
    // Need a scalar field and a ColorMap
    infield = scinew FieldIPort( this, "ScalarField",
					FieldIPort::Atomic);
    add_iport( infield);
    inColorMap = scinew ColorMapIPort( this, "ColorMap",
				     ColorMapIPort::Atomic);
    add_iport( inColorMap);
					
    // Create the output port
    ogeom = scinew GeometryOPort(this, "Geometry", 
			      GeometryIPort::Atomic);
    add_oport(ogeom);
    init = 1;
    float INIT(.1);

    widget = scinew ScaledFrameWidget(this, &widget_lock, INIT);
    grid_id=0;

    need_find.set(1);
    
    outcolor=scinew Material(Color(0,0,0), Color(0,0,0), Color(0,0,0), 0);
    msg="";
}

CuttingPlane::~CuttingPlane()
{
}

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1184
#endif

void CuttingPlane::get_minmax(FieldHandle field) {
  uintah_dispatch_scalar1(field, dispatch_minmax);
}


void CuttingPlane::execute()
{
    int old_grid_id = grid_id;
    static int find = -1;
    int cmapmin, cmapmax;

    // get the scalar field and ColorMap...if you can
    FieldHandle field;
    if (!infield->get( field ))
	return;

    if(!field->is_scalar()){
      cerr<<"Not a scalar field\n ";
      return;
    }
      
    ColorMapHandle cmap;
    if (!inColorMap->get( cmap ))
	return;
    cmapmin=(int)cmap->getMin();
    cmapmax=(int)cmap->getMax();

    if (init == 1) 
    {
	init = 0;
	GeomObj *w = widget->GetWidget() ;
	widget_id = ogeom->addObj( w, widget_name, &widget_lock );
	widget->Connect( ogeom );
	widget->SetRatioR( 0.4 );
	widget->SetRatioD( 0.4 );
    }
    if (need_find.get() != find)
    {
      BBox box;
      Point min, max;
      box = field->mesh()->get_bounding_box();
      min = box.min(); max = box.max();
      Point center = min + (max-min)/2.0;
	double max_scale;
	double wh=where.get();
	if (need_find.get() == 1)
	{   // Find the field and put in optimal place
	    // in xy plane with reasonable frame thickness
	    center.z(min.z()*wh+max.z()*(1-wh));
	    Point right( max.x(), center.y(), center.z());
	    Point down( center.x(), min.y(), center.z());
	    max_scale = Max( (max.x() - min.x()), (max.y() - min.y()) );
	    widget->SetScale( max_scale/30. );
	    widget->SetPosition( center, right, down);
	}
	else if (need_find.get() == 2)
	{   // Find the field and put in optimal place
	    // in yz plane with reasonable frame thickness
	    center.x(min.x()*wh+max.x()*(1-wh));
	    Point right( center.x(), center.y(), max.z());
	    Point down( center.x(), min.y(), center.z());	    
	    max_scale = Max( (max.z() - min.z()), (max.y() - min.y()) );
	    widget->SetScale( max_scale/30. );
	    widget->SetPosition( center, right, down);
	}
	else
	{   // Find the field and put in optimal place
	    // in xz plane with reasonable frame thickness
	    center.y(min.y()*wh+max.y()*(1-wh));
	    Point right( max.x(), center.y(), center.z());
	    Point down( center.x(), center.y(), min.z());	    
	    max_scale = Max( (max.x() - min.x()), (max.z() - min.z()) );
	    widget->SetScale( max_scale/30. );
	    widget->SetPosition( center, right, down);
	}
	find = need_find.get();
    }

    // advance or decrement along x, y, or z
    if (msg != "" &&
	(field->get_type_name(0) == "LatticeVol" ||
	 field->get_type_name(0) == "LevelField")) {
      int nx, ny, nz;
      get_dimensions(field, nx,ny,nz);
      if(field->data_at() == Field::CELL){ nx--; ny--; nz--; }
      Point center, right, down;
      widget->GetPosition(center, right, down);
      BBox box;
      Point min, max;
      Vector diag;
      box = field->mesh()->get_bounding_box();
      min = box.min(); max = box.max();
      diag=max-min;
      diag.x(diag.x()/(nx-1));
      diag.y(diag.y()/(ny-1));
      diag.z(diag.z()/(nz-1));
      if (msg=="plusx") {
	center.x(center.x()+diag.x());
	right.x(right.x()+diag.x());
	down.x(down.x()+diag.x());
      } else if (msg=="minusx") {
	center.x(center.x()-diag.x());
	right.x(right.x()-diag.x());
	down.x(down.x()-diag.x());
      } else if (msg=="plusy") {
	center.y(center.y()+diag.y());
	right.y(right.y()+diag.y());
	down.y(down.y()+diag.y());
      } else if (msg=="minusy") {
	center.y(center.y()-diag.y());
	right.y(right.y()-diag.y());
	down.y(down.y()-diag.y());
      } else if (msg=="plusz") {
	center.z(center.z()+diag.z());
	right.z(right.z()+diag.z());
	down.z(down.z()+diag.z());
      } else { // if (msg=="minusz")
	center.z(center.z()-diag.z());
	right.z(right.z()-diag.z());
	down.z(down.z()-diag.z());
      }
      widget->SetPosition( center, right, down);
      msg="";
    }

    // get the position of the frame widget
    Point 	corner, center, R, D;
    widget->GetPosition( center, R, D);
    Vector v1 = R - center,
      v2 = D - center;
         
    // calculate the corner and the
    // u and v vectors of the cutting plane
    
    corner = (center - v1) - v2;
    Vector u = v1 * 2.0,
           v = v2 * 2.0;

    int cptype = cutting_plane_type.get();
    
    // create the grid for the cutting plane
    double u_fac = widget->GetRatioR(),
           v_fac = widget->GetRatioD(),
           scale_fac = scale.get(),
           offset_fac = offset.get();

    int u_num = (int) (u_fac * 500),
        v_num = (int) (v_fac * 500);


    if (fullRezGUI.get()) {
      if(field->get_type_name(0) != "LevelField" &&
	 field->get_type_name(0) != "LatticeVol"){
	cerr << "Error - not a regular grid... can't use Full Resolution!\n";
      } else {
	int nx, ny, nz;
	get_dimensions(field, nx,ny,nz);
	if(field->data_at() == Field::CELL){ nx--; ny--; nz--; }
	int most=Max(Max(nx, ny), nz);
	u_num=v_num=most;
      }
    }
    //    cout << "u fac = " << u_fac << "\nv fac = " << v_fac << endl;
    
    int localMinMax=localMinMaxGUI.get();

    //int exhaustive=exhaustiveGUI.get();

    // Get the scalar values and corresponding
    // colors to put in the cutting plane
    if (cptype != CP_CONTOUR)
      {
#if 0
	double alpha=1.0;
#endif
	//	GeomGrid::Format format=GeomGrid::WithMaterials;
	//	if (cptype == CP_SURFACE)
	//	    format=GeomGrid::WithNormAndMatl;
	GeomGrid* grid = new GeomGrid( u_num, v_num, corner, u, v,
				       0);
        Vector unorm=u.normal();
        Vector vnorm=v.normal();
        Vector N(Cross(unorm, vnorm));
#if 0
	if (cptype != CP_SURFACE)
	  grid->texture();
#endif

	//int ix = 0;
	int i, j;

	int haveval=0;
	double minval, maxval, invrange;

	if (localMinMax) {
	    // get the min and max values from this slice
	    for (i = 0; i < u_num; i++)
		for (j = 0; j < v_num; j++) {
		    Point p = corner + u * ((double) i/(u_num-1)) + 
			v * ((double) j/(v_num-1));
		    double sval;
// 		    if (sfield->interpolate( p, sval, ix) || 
// 			(ix=0) || 
// 			sfield->
// 			  interpolate( p, sval, ix, EPS, EPS, exhaustive)) {
		    if( interpolate( field, p, sval) ){
		      if (!haveval) { minval=maxval=sval; haveval=1; }
		      else if (sval<minval) minval=sval;
		      else if (sval>maxval) maxval=sval;
		    }

		}
	    invrange=(cmapmax-cmapmin)/(maxval-minval);
	}

	for (i = 0; i < u_num; i++)
	  for (j = 0; j < v_num; j++) {
	    Point p = corner + u * ((double) i/(u_num-1)) + 
	      v * ((double) j/(v_num-1));
	    double sval;
	    
	    // get the color from cmap for p 	    
	    MaterialHandle matl;
// 	    if (sfield->interpolate( p, sval, ix) 
// 		|| (ix=0) || 
// 		sfield->interpolate( p, sval, ix, EPS, EPS, exhaustive)) {
	    if( interpolate( field, p, sval) ){
		if (localMinMax)	// use local min/max to scale
		    sval=(sval-minval)*invrange+cmapmin;
		matl = cmap->lookup( sval);
#if 0
		alpha = 0.8;
#endif
	    } else {
	      matl = outcolor;
	      sval = 0;
#if 0
	      alpha=0.0;
#endif
	    }
	    
#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 255
#endif

	    // put the color into the cutting plane (grid) at i, j
	    Vector G;
	    if (cptype == CP_SURFACE && get_gradient(field, p, G))
	      {
		double h = sval;
		double umag=Dot(unorm, G)*scale_fac;
		double vmag=Dot(vnorm, G)*scale_fac;
		Vector normal(N-unorm*umag-vnorm*vmag);
		grid->set(i, j, ((h*scale_fac) + offset_fac),
			  normal, matl/*,alpha*/);
	      }
	    else  			// if (cptype == CP_PLANE)
	      grid->set(i, j, 0, matl/*,alpha*/);
	  }
	// delete the old grid/cutting plane
	if (old_grid_id != 0) {
	  ogeom->delObj( old_grid_id );
	}
	grid_id = ogeom->addObj(grid, "Cutting Plane");
	old_grid_id = grid_id;
	
      }

    else
    {
	double minval, maxval, invrange;
	get_minmax(field);
	minval = minmax_.first;
	maxval = minmax_.second;

	if (localMinMax) {
	    // get the min and max values from this slice
	    // int ix = 0;
	    int i, j;
	    int haveval=0;
 
	    for (i = 0; i < u_num; i++)
		for (j = 0; j < v_num; j++) {
		    Point p = corner + u * ((double) i/(u_num-1)) + 
			v * ((double) j/(v_num-1));
		    double sval;
// 		    if (sfield->interpolate( p, sval, ix) || 
// 			(ix=0) || 
// 			sfield->
// 			interpolate( p, sval, ix, EPS, EPS, exhaustive)) {
		    if( interpolate( field, p, sval) ){
			if (!haveval) { minval=maxval=sval; haveval=1; }
			else if (sval<minval) minval=sval;
			else if (sval>maxval) maxval=sval;
		    }
		}
	    invrange=(cmapmax-cmapmin)/(maxval-minval);
	} else {
	    // get the min and max values from the field
	  get_minmax(field);
	  minval = minmax_.first;
	  maxval = minmax_.second;
	}
	// get the number of contours
	int contours = num_contours.get();
	if (contours >= 2)
	{
	    GeomGroup *cs = new GeomGroup;
	    Array1<GeomLines *> col_group( contours);
	    Array1<GeomMaterial *> colrs( contours);
	    Array1<double> values( contours);

	    // get the "contour" number of values from 
	    // the field, find corresponding colors,	
	    // and add a Material and it's group to the tree
	    int i;
	    for (i = 0; i < contours; i++)
	    {
		if (localMinMax) values[i]=i/(contours-1.)*
				     (cmapmax-cmapmin)+cmapmin;
		else values[i]=((double)i/(contours - 1)) * (maxval - minval) + minval;
		MaterialHandle matl;
		matl = cmap->lookup( values[i] );
		col_group[i] = new GeomLines;
		colrs[i] = new GeomMaterial( col_group[i], matl );
		cs->add( colrs[i] );
	    }
	    //int ix=0;
	    // look at areas in the plane to find the contours
	    for (i = 0; i < u_num-1; i++)
		for (int j = 0; j < v_num-1; j++)
		{
		    Point p1 = corner + u * ((double) i/(u_num-1)) + 
			                v * ((double) j/(v_num-1));
		    Point p2 = corner + u * ((double) (i+1)/(u_num-1)) + 
			                v * ((double) j/(v_num-1));
		    Point p3 = corner + u * ((double) (i+1)/(u_num-1)) + 
			                v * ((double) (j+1)/(v_num-1));
		    Point p4 = corner + u * ((double) i/(u_num-1)) + 
			                v * ((double) (j+1)/(v_num-1));
		    double sval1;
		    double sval2;
		    double sval3;
		    double sval4;

		    // get the value from the field for each point	    
// 		    if ( (sfield->interpolate( p1,sval1,ix) || (ix=0) || 
// 			  sfield->interpolate( p1, sval1, ix, EPS, EPS, exhaustive)) && 
// 			(sfield->interpolate( p2,sval2,ix) || (ix=0) || 
// 			 sfield->interpolate( p2, sval2, ix, EPS, EPS, exhaustive)) && 
// 			(sfield->interpolate( p3,sval3, ix) || (ix=0) || 
// 			 sfield->interpolate( p3, sval3, ix, EPS, EPS, exhaustive)) && 
// 			(sfield->interpolate( p4, sval4,ix) || (ix=0) || 
// 			 sfield->interpolate( p4, sval4, ix, EPS, EPS, exhaustive)))
		    if (interpolate(field, p1, sval1) &&
			interpolate(field, p2, sval2) &&
			interpolate(field, p3, sval3) &&
			interpolate(field, p4, sval4))
		    {
			if (localMinMax) {	// use local min/max to scale
			    sval1=(sval1-minval)*invrange+cmapmin;
			    sval2=(sval2-minval)*invrange+cmapmin;
			    sval3=(sval3-minval)*invrange+cmapmin;
			    sval4=(sval4-minval)*invrange+cmapmin;
			}			    
			// find the indices of the values array between smin & smax
			double smin = Min( Min( sval1, sval2), Min( sval3, sval4)),
			       smax = Max( Max( sval1, sval2), Max( sval3, sval4));
			int i1, i2;
			if (localMinMax) {
			    i1 = RoundUp((smin-cmapmin)*(contours - 1)/(cmapmax-cmapmin));
			    i2 = RoundDown((smax-cmapmin)*(contours - 1)/(cmapmax-cmapmin));
			} else {
			    i1=RoundUp((smin-minval)*(contours - 1)/(maxval - minval));
			    i2=RoundDown((smax-minval)*(contours-1)/(maxval - minval));
			}
			if(!localMinMax && (smin < minval || smax > maxval))
			{
			    cerr << "OOPS: " << endl;
			    cerr << "smin, smax=" << smin << " " << smax << endl;
			    cerr << "minval, maxval=" << minval << " " << maxval << endl;
			    continue;
			}
			
			// find and add the contour lines if they exist in this area
			for (int k = i1; k <= i2; k++)
			{
			    int found = 0;
			    Point x1, x2;

			    if ((sval1 <= values[k] && values[k] < sval2) ||
				(sval2 < values[k] && values[k] <= sval1))
			    {
				x1 = Interpolate( p1, p2, (values[k]-sval1)/(sval2-sval1));
				++found;
			    }
			    if ((sval2 <= values[k] && values[k] < sval3) ||
				(sval3 < values[k] && values[k] <= sval2))
			    {
				x2 = Interpolate( p2, p3, (values[k]-sval2)/(sval3-sval2));
				if (!found)
				    x1 = x2;
				++found;
			    }
			    if (((sval3 <= values[k] && values[k] < sval4) ||
				 (sval4 < values[k] && values[k] <= sval3)) && found < 2)
			    {
				x2 = Interpolate( p3, p4, (values[k]-sval3)/(sval4-sval3));
				if (!found)
				    x1 = x2;
				++found;
			    }
			    if (((sval1 < values[k] && values[k] <= sval4) ||
				 (sval4 <= values[k] && values[k] < sval1)) && found < 2)
			    {
				x2 = Interpolate( p1, p4, (values[k]-sval1)/(sval4-sval1));
				++found;
			    }
			    // did we find two points to draw a line with?
			    if (found == 2)
				col_group[k]->add(x1, x2);
			}
		    }
		}
    // delete the old grid/cutting plane
	    if (old_grid_id != 0) {
	      ogeom->delObj( old_grid_id );

	    }
	    grid_id =  ogeom->addObj( cs, "Contour Plane");
	    old_grid_id = grid_id;
		      
	}
	
    }

}

void CuttingPlane::widget_moved(int last)
{
    if(last && !abort_flag)
    {
	abort_flag=1;
	want_to_execute();
    }
}


void CuttingPlane::tcl_command(TCLArgs& args, void* userdata)
{
    if(args.count() < 2)
    {
	args.error("Streamline needs a minor command");
	return;
    }
    if(args[1] == "findxy")
    {
	need_find.set(1);
	want_to_execute();
    }
    else if(args[1] == "findyz")
    {
	need_find.set(2);
	want_to_execute();
    }
    else if(args[1] == "findxz")
    {
	need_find.set(3);
	want_to_execute();
    }
    else if(args[1] == "plusx")
    {
	msg="plusx";
	want_to_execute();
    }
    else if(args[1] == "minusx")
    {
	msg="minusx";
	want_to_execute();
    }
    else if(args[1] == "plusy")
    {
	msg="plusy";
	want_to_execute();
    }
    else if(args[1] == "minusy")
    {
	msg="minusy";
	want_to_execute();
    }
    else if(args[1] == "plusz")
    {
	msg="plusz";
	want_to_execute();
    }
    else if(args[1] == "minusz")
    {
	msg="minusz";
	want_to_execute();
    }
    else if(args[1] == "connectivity")
    {
	
    }
    else
    {
	Module::tcl_command(args, userdata);
    }
}

template <class Mesh>
bool 
CuttingPlane::get_dimensions(Mesh, int&, int&, int&)
  {
    return false;
  }

template<> 
bool CuttingPlane::get_dimensions(LevelMeshHandle m,
				 int& nx, int& ny, int& nz)
  {
    nx = m->get_nx();
    ny = m->get_ny();
    nz = m->get_nz();
    return true;
  }
template<> 
bool CuttingPlane::get_dimensions(LatVolMeshHandle m,
				 int& nx, int& ny, int& nz)
  {
    nx = m->get_nx();
    ny = m->get_ny();
    nz = m->get_nz();
    return true;
  }

bool
CuttingPlane::get_dimensions(FieldHandle texfld_,  int& nx, int& ny, int& nz)
{
  const string type = texfld_->get_type_name(1);
  if( texfld_->get_type_name(0) == "LevelField" ){
    LevelMeshHandle mesh_;
    if (type == "double") {
      LevelField<double> *fld =
	dynamic_cast<LevelField<double>*>(texfld_.get_rep());
      mesh_ = fld->get_typed_mesh();
    } else if (type == "float") {
      LevelField<float> *fld =
	dynamic_cast<LevelField<float>*>(texfld_.get_rep());
      mesh_ = fld->get_typed_mesh();
    } else if (type == "long") {
      LevelField<long> *fld =
	dynamic_cast<LevelField<long>*>(texfld_.get_rep());
      mesh_ = fld->get_typed_mesh();
    } else {
      cerr << "CuttingPlane error - unknown LevelField type: " << type << endl;
      return false;
    }
    return get_dimensions( mesh_, nx, ny, nz );
  } else if(texfld_->get_type_name(0) == "LatticeVol"){
    LatVolMeshHandle mesh_;
    if (type == "double") {
      LatticeVol<double> *fld =
	dynamic_cast<LatticeVol<double>*>(texfld_.get_rep());
      mesh_ = fld->get_typed_mesh();
    } else if (type == "float") {
      LatticeVol<float> *fld =
	dynamic_cast<LatticeVol<float>*>(texfld_.get_rep());
      mesh_ = fld->get_typed_mesh();
    } else if (type == "long") {
      LatticeVol<long> *fld =
	dynamic_cast<LatticeVol<long>*>(texfld_.get_rep());
      mesh_ = fld->get_typed_mesh();
    } else if (type == "int") {
      LatticeVol<int> *fld =
	dynamic_cast<LatticeVol<int>*>(texfld_.get_rep());
      mesh_ = fld->get_typed_mesh();
    } else {
      cerr << "CuttingPlane error - unknown LatticeVol type: " << type << endl;
      return false;
    }
    return get_dimensions( mesh_, nx, ny, nz );
  } else {
    return false;
  }
}


bool  
CuttingPlane::get_gradient(FieldHandle texfld_, const Point& p, Vector& g)
{
  const string field_type = texfld_->get_type_name(0);
  const string type = texfld_->get_type_name(1);
  if( texfld_->get_type_name(0) == "LevelField" ){
    if (type == "double") {
      LevelField<double> *fld =
	dynamic_cast<LevelField<double>*>(texfld_.get_rep());
      return fld->get_gradient(g,p);
    } else if (type == "float") {
      LevelField<float> *fld =
	dynamic_cast<LevelField<float>*>(texfld_.get_rep());
      return fld->get_gradient(g,p);   
    } else if (type == "long") {
      LevelField<long> *fld =
	dynamic_cast<LevelField<long>*>(texfld_.get_rep());
      return fld->get_gradient(g,p);    
    } else {
      cerr << "CuttingPlane::get_gradient:: error - unimplemented Field type: " << type << endl;
      return false;
    }
  } else if(texfld_->get_type_name(0) == "LatticeVol"){
    if (type == "double") {
      LatticeVol<double> *fld =
	dynamic_cast<LatticeVol<double>*>(texfld_.get_rep());
      return fld->get_gradient(g,p);
    } else if (type == "float") {
      LatticeVol<float> *fld =
	dynamic_cast<LatticeVol<float>*>(texfld_.get_rep());
      return fld->get_gradient(g,p);   
    } else if (type == "long") {
      LatticeVol<long> *fld =
	dynamic_cast<LatticeVol<long>*>(texfld_.get_rep());
      return fld->get_gradient(g,p);    
    } else if (type == "int") {
      LatticeVol<int> *fld =
	dynamic_cast<LatticeVol<int>*>(texfld_.get_rep());
      return fld->get_gradient(g,p);    
    } else {
      cerr << "CuttingPlane::get_gradient:: error - unimplemented Field type: " << type << endl;
      return false;
    }
  } else {
    return false;
  }
}

bool  
CuttingPlane::interpolate(FieldHandle texfld_, const Point& p, double& val)
{
  const string field_type = texfld_->get_type_name(0);
  const string type = texfld_->get_type_name(1);
  if( texfld_->get_type_name(0) == "LevelField" ){
    if (type == "double") {
      LevelField<double> *fld =
	dynamic_cast<LevelField<double>*>(texfld_.get_rep());
      return fld->interpolate(val ,p);
    } else if (type == "float") {
      float result;
      bool success;
      LevelField<float> *fld =
	dynamic_cast<LevelField<float>*>(texfld_.get_rep());
      success = fld->interpolate(result ,p);   
      val = (double)result;
      return success;
    } else if (type == "long") {
      long result;
      bool success;
      LevelField<long> *fld =
	dynamic_cast<LevelField<long>*>(texfld_.get_rep());
      success =  fld->interpolate(result,p);
      val = (double)result;
      return success;
    } else {
      cerr << "Uintah::CuttingPlane::interpolate:: error - unimplemented Field type: " << type << endl;
      return false;
    }
  } else if( texfld_->get_type_name(0) == "LatticeVol" ){
    if (type == "double") {
      LatticeVol<double> *fld =
	dynamic_cast<LatticeVol<double>*>(texfld_.get_rep());
      return fld->interpolate(val ,p);
    } else if (type == "float") {
      float result;
      bool success;
      LatticeVol<float> *fld =
	dynamic_cast<LatticeVol<float>*>(texfld_.get_rep());
      success = fld->interpolate(result ,p);   
      val = (double)result;
      return success;
    } else if (type == "long") {
      long result;
      bool success;
      LatticeVol<long> *fld =
	dynamic_cast<LatticeVol<long>*>(texfld_.get_rep());
      success =  fld->interpolate(result,p);
      val = (double)result;
      return success;
    } else if (type == "int") {
      int result;
      bool success;
      LatticeVol<int> *fld =
	dynamic_cast<LatticeVol<int>*>(texfld_.get_rep());
      success =  fld->interpolate(result,p);
      val = (double)result;
      return success;
    } else {
      cerr << "Uintah::CuttingPlane::interpolate:: error - unimplemented Field type: " << type << endl;
      return false;
    }
  } else {
    return false;
  }
}


} // End namespace Uintah
