
/*
 *  BldTransform.cc:  Build a 4x4 transformation matrix
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   March 1999
 *
 *  Copyright (C) 1999 SCI Group
 */

#include <Dataflow/Ports/MatrixPort.h>
#include <Dataflow/Ports/GeometryPort.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Geom/Switch.h>
#include <Core/Geometry/BBox.h>
#include <Core/Geometry/Transform.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Expon.h>
#include <Core/Math/MusilRNG.h>
#include <Core/Math/Trig.h>
#include <Core/TclInterface/TCLvar.h>
#include <Core/Thread/CrowdMonitor.h>
#include <Dataflow/Widgets/ScaledBoxWidget.h>
#include <iostream>
using std::cerr;
#include <stdio.h>

namespace SCIRun {


class BldTransform : public Module {
    MatrixIPort* imatrix;
    MatrixOPort* omatrix;
    MatrixHandle omh;
    TCLdouble rx, ry, rz, th;
    TCLdouble tx, ty, tz;
    TCLdouble scale, scalex, scaley, scalez;
    TCLdouble sha, shb, shc, shd;
    TCLint xmapTCL;
    TCLint ymapTCL;
    TCLint zmapTCL;
    TCLint pre;
    TCLint whichxform;
    TCLint widgetShowResizeHandles;
    TCLdouble widgetScale;
    ScaledBoxWidget *boxWidget;
    GeomSwitch *widget_switch;
    CrowdMonitor widget_lock;
    GeometryOPort* ogeom;
    Transform composite, latestT, latestWidgetT, widget_trans, widget_pose_inv;
    Point widget_pose_center;
    int ignorechanges;
    int init;
public:
    BldTransform(const clString& id);
    virtual ~BldTransform();
    virtual void widget_moved(int last);
    virtual void execute();
    void tcl_command( TCLArgs&, void * );
};

extern "C" Module* make_BldTransform(const clString& id)
{
    return new BldTransform(id);
}

static clString module_name("BldTransform");
static clString widget_name("TransformWidget");

BldTransform::BldTransform(const clString& id)
: Module("BldTransform", id, Filter),
  rx("rx", id, this), ry("ry", id, this), rz("rz", id, this), 
  th("th", id, this),
  tx("tx", id, this), ty("ty", id, this), tz("tz", id, this),
  scalex("scalex", id, this), scaley("scaley", id, this), 
  scalez("scalez", id, this), 
  scale("scale", id, this), pre("pre", id, this),
  xmapTCL("xmapTCL", id, this), ymapTCL("ymapTCL", id, this), 
  zmapTCL("zmapTCL", id, this),
  whichxform("whichxform", id, this),
  sha("sha", id, this), shb("shb", id, this), shc("shc", id, this),
  shd("shd", id, this), widget_lock("BldTransform widget lock"),
  widgetScale("widgetScale", id, this), ignorechanges(1),
  init(0), widgetShowResizeHandles("widgetShowResizeHandles", id, this)
{
    // Create the input port
    imatrix=scinew MatrixIPort(this, "Matrix", MatrixIPort::Atomic);
    add_iport(imatrix);

    // Create the output ports
    omatrix=scinew MatrixOPort(this, "Matrix", MatrixIPort::Atomic);
    add_oport(omatrix);
    ogeom=scinew GeometryOPort(this, "Geometry", GeometryIPort::Atomic);
    add_oport(ogeom);
    boxWidget=scinew ScaledBoxWidget(this, &widget_lock, 0.2);
    widget_switch=boxWidget->GetWidget();
}

BldTransform::~BldTransform()
{
}

void BldTransform::execute()
{
    int wh=whichxform.get();
    if (!init) {
	Point C, R, D, I;
	boxWidget->GetPosition(C,R,D,I);
	C=Point(0,0,0); R=Point(1,0,0); D=Point(0,1,0), I=Point(0,0,1);
	widget_pose_center=C;
	boxWidget->SetPosition(C,R,D,I);
	boxWidget->SetCurrentMode(2);
	if (wh != 5) widget_switch->set_state(0);
	ogeom->addObj(widget_switch, widget_name, &widget_lock);
	ogeom->flushViews();
	init=1;
    }
    int i, j;

    // get the input matrix if there is one
    MatrixHandle imh;
    Matrix* im;
    Transform inT;
    if (imatrix->get(imh) && (im=imh.get_rep())) {
	double inV[16];
	double *p=&(inV[0]);
	for (i=0; i<4; i++)
	    for (j=0; j<4; j++)
		*p++=(*im)[i][j];
	inT.set(inV);
    }

    Transform locT;

    // get the "fixed point"
    double txx=tx.get();
    double tyy=ty.get();
    double tzz=tz.get();
    Vector t(txx, tyy, tzz);

    // switch on the message and build the local matrix accordingly
    if (wh==0) {			       // TRANSLATE
	locT.post_translate(t);
    } else if (wh==1) {                        // SCALE
	double new_scale=scale.get();
	double s=pow(10.,new_scale);
	double new_scalex=scalex.get();
	double sx=pow(10.,new_scalex)*s;
	double new_scaley=scaley.get();
	double sy=pow(10.,new_scaley)*s;
	double new_scalez=scalez.get();
	double sz=pow(10.,new_scalez)*s;
	Vector sc(sx, sy, sz);
	locT.post_translate(t);	
	locT.post_scale(sc);
	locT.post_translate(-t);
    } else if (wh==2) {			       // ROTATE
	Vector axis(rx.get(),ry.get(),rz.get());
	if (!axis.length2()) axis.x(1);
	axis.normalize();
	locT.post_translate(t);
	locT.post_rotate(th.get()*M_PI/180., axis);
	locT.post_translate(-t);
    } else if (wh==3) {      		       // SHEAR
	locT.post_shear(t, Plane(sha.get(), shb.get(), shc.get(), shd.get()));
    } else if (wh==4) {			       // PERMUTE
	locT.post_permute(xmapTCL.get(), ymapTCL.get(), zmapTCL.get());
    } else { // (wh==5)			       // WIDGET
	Point R, D, I, C;
	boxWidget->GetPosition(C, R, D, I);

	double ratio=widgetScale.get();
	widgetScale.set(1);
	R=C+(R-C)*ratio;
	D=C+(D-C)*ratio;
	I=C+(I-C)*ratio;
	boxWidget->SetPosition(C, R, D, I);

	// find the difference between widget_pose(_inv) and the current pose
	if (!ignorechanges) {
	    locT.load_frame(C,R-C,D-C,I-C);
	    locT.post_trans(widget_pose_inv);
	    locT.post_translate(-widget_pose_center.vector());
	    locT.pre_translate(C.vector());
	}
	locT.post_trans(latestWidgetT);
	latestWidgetT=locT;
	widget_pose_center=C;
	widget_pose_inv.load_frame(C,R-C,D-C,I-C);
	widget_pose_inv.invert();
    }
    DenseMatrix *dm=scinew DenseMatrix(4,4);
    omh=dm;
    
    // now either pre- or post-multiply the transforms and store in matrix
    if (pre.get()) {
	locT.post_trans(composite);
	latestT=locT;
	locT.post_trans(inT);
    } else {
	locT.pre_trans(composite);
	latestT=locT;
	locT.pre_trans(inT);
    }
    double finalP[16];
    locT.get(finalP);
    double *p=&(finalP[0]);
    int cnt=0;
    for (i=0; i<4; i++) 
	for (j=0; j<4; j++, cnt++)
	    (*dm)[i][j]=*p++;

    omatrix->send(omh);
}

void BldTransform::widget_moved(int last)
{
    // only re-execute if this was a widget-release event
    if (last) {
	want_to_execute();
    }
}

void BldTransform::tcl_command(TCLArgs& args, void* userdata) {
    if (args[1] == "hide_widget") {
	widget_switch->set_state(0);
	ogeom->flushViews();
    } else if (args[1] == "show_widget") {
	widget_switch->set_state(1);
	ogeom->flushViews();
    } else if (args[1] == "reset_widget" || args[1] == "reset" || 
	       args[1] == "composite") {
	if (args[1] == "reset")
	    composite.load_identity();
	else if (args[1] == "composite")
	    composite=latestT;
	latestT.load_identity();
	latestWidgetT.load_identity();
        boxWidget->SetPosition(Point(0,0,0), 
			       Point(1,0,0), Point(0,1,0), Point(0,0,1));
	widget_pose_center=Point(0,0,0);
	widget_pose_inv.load_identity();
	want_to_execute();
    } else if (args[1] == "change_ignore") {
	if (args[2] == "1") {	// start ignoring widget changes
	    ignorechanges=1;
	} else {		// stop ignoring widget changes
	    ignorechanges=0;
	}
    } else if (args[1] == "change_handles") {
	if (args[2] == "1") {	// start showing resize handles
	    boxWidget->SetCurrentMode(1);
	    ogeom->flushViews();
	} else {		// stop showing resize handles
	    boxWidget->SetCurrentMode(2);
	    ogeom->flushViews();
	}
    } else {
        Module::tcl_command(args, userdata);
    }
}

} // End namespace SCIRun
