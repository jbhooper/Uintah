
/*
 *  BitVisualize.cc:  IsoSurfaces a SFRG bitwise
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   February 1995
 *
 *  Copyright (C) 1994 SCI Group
 */

#include <config.h>
#undef SCI_ASSERTION_LEVEL_3
#define SCI_ASSERTION_LEVEL_2
#include <Tester/RigorousTest.h>
#include <Classlib/BitArray1.h>
#include <Classlib/HashTable.h>
#include <Classlib/NotFinished.h>
#include <Dataflow/Module.h>
#include <Datatypes/GeometryPort.h>
#include <Datatypes/Mesh.h>
#include <Datatypes/ScalarField.h>
#include <Datatypes/ScalarFieldRG.h>
#include <Datatypes/ScalarFieldUG.h>
#include <Datatypes/ScalarFieldPort.h>
#include <Datatypes/SurfacePort.h>
#include <Datatypes/TriSurface.h>
#include <Geometry/Point.h>
#include <Geom/Geom.h>
#include <Geom/Pt.h>
#include <Geom/Group.h>
#include <Geom/Material.h>
#include <Geom/Tri.h>
#include <Malloc/Allocator.h>
#include <Multitask/Task.h>
#include <TCL/TCLvar.h>
#include <iostream.h>
#include <strstream.h>

// just so I can see the proccess id...

#include <sys/types.h>
#include <unistd.h>


#define NUM_MATERIALS 5

class BitVisualize : public Module {
    ScalarFieldIPort* infield;

    GeometryOPort* ogeom;
    Array1<SurfaceOPort* > osurfs;
    Array1<SurfaceHandle> surf_hands;
    Array1<TriSurface*> surfs;
    Array1<TCLint*> isovals;
    Array1<int> BitVisualize_id;
    Array1<int> calc_mat;
    int isoChanged;
    ScalarFieldRG* sfrg;
    Array1<MaterialHandle> matls;
    Array1<GeomGroup*> groups;
    Array1<int> geom_allocated;
    void vol_render1_grid();

    Point ov[9];
    Point v[9];
    Point maxP;
    Point minP;
    Point pmid;
    int bit;

    Mutex grouplock;
    int np;
    void vol_render1();
public:
    void parallel_vol_render1(int proc);
    BitVisualize(const clString& id);
    BitVisualize(const BitVisualize&, int deep);
    virtual ~BitVisualize();
    virtual Module* clone(int deep);
    virtual void execute();
};

#define FACE1 8
#define FACE2 4
#define FACE3 2
#define FACE4 1
#define ALLFACES (FACE1|FACE2|FACE3|FACE4)

#if 0
struct MCubeTable {
    int which_case;
    int permute[8];
    int nbrs;
};

#include "mcube.h"
#endif

extern "C" {
Module* make_BitVisualize(const clString& id)
{
    return scinew BitVisualize(id);
}
}

BitVisualize::BitVisualize(const clString& id)
: Module("BitVisualize", id, Filter)
{
    // Create the input ports
    infield=scinew ScalarFieldIPort(this, "Field", ScalarFieldIPort::Atomic);
    add_iport(infield);
    // Create the output port
    ogeom=scinew GeometryOPort(this, "Geometry", GeometryIPort::Atomic);
    add_oport(ogeom);
    groups.resize(NUM_MATERIALS);
    surfs.resize(NUM_MATERIALS);
    surf_hands.resize(NUM_MATERIALS);
    for (int i=0; i<NUM_MATERIALS; i++) {
	osurfs.add(scinew SurfaceOPort(this, "Surface", SurfaceIPort::Atomic));
	surf_hands[i]=surfs[i]=0;
	add_oport(osurfs[i]);
	geom_allocated.add(0);
	calc_mat.add(0);
	BitVisualize_id.add(0);
	surfs.add(scinew TriSurface);
	clString str;
	str = "iso" + to_string(i+1);
	isovals.add(scinew TCLint(str, id, this));
	isovals[i]->set(0);
	int r=(i+1)%2;		// 1 0 1 0 1 0
	int g=(((i+2)/3)%2);	// 0 1 1 1 0 0
	int b=((i/3)%2);	// 0 0 0 1 1 1
	matls.add(scinew Material(Color(.2,.2,.2), Color(r*.6, g*.6, b*.6), 
			       Color(r*.5, g*.5, b*.5), 20));
    }	
}

BitVisualize::BitVisualize(const BitVisualize& copy, int deep)
: Module(copy, deep)
{
    NOT_FINISHED("BitVisualize::BitVisualize");
}

BitVisualize::~BitVisualize()
{
}

Module* BitVisualize::clone(int deep)
{
    return scinew BitVisualize(*this, deep);
}

void BitVisualize::execute()
{
    ScalarFieldHandle field;
    if(!infield->get(field))
	return;
    sfrg=field->getRG();
    if(sfrg){
//	isoChanged=(int)sfrg->grid(0,0,0);
//	isoChanged=(1<<NUM_MATERIALS)-1;
	sfrg->grid(0,0,0)=0;
	int i;
	for (i=0; i<NUM_MATERIALS; i++) {
//	    bit=1<<i;
//	    if (isoChanged & bit) {		    // fld changed...
		if (isovals[i]->get()) {	    //   want this material...
		    if (geom_allocated[i]) {	    //     it's allocated...
			ogeom->delObj(BitVisualize_id[i]);//	      DELETE
		    } else {			    //     not allocated
			geom_allocated[i]=1; 	    // 	      MARK as ALLOCED
		    }
		    groups[i]=scinew GeomGroup;	    //	   ALLOCATE THEM
		    calc_mat[i]=1;		    //	   *Calculate mat.*
		} else {			    //   don't want material...
		    if (geom_allocated[i]) {	    //     it's allocated...
			ogeom->delObj(BitVisualize_id[i]);//	 DELETE
			geom_allocated[i]=0;	    //	      MARK as UNALLOCED
		    }
		    calc_mat[i]=0;		    //	   *Don't calc. mat.*
		}
//	    } else {				    // fld didn't change...
//		if (isovals[i]->get()) {	    //	 want this material...
//		    if (geom_allocated[i]) {	    //	   it's allocated...
//			calc_mat[i]=0;		    //	      *Don't calc. mat*
//		    } else {			    //	   not allocated
//			geom_allocated[i]=1;	    //	      MARK
//			groups[i]=scinew GeomGroup; //	      ALLOCATE THEM
//			geomPts[i]=scinew GeomGroup;
//			calc_mat[i]=1;		    //	      *Calcluate mat*
//		    }		
//		} else {			    //   don't want material
//		    calc_mat[i]=0;		    //     *Don't calc mat*
//		    if (geom_allocated[i]) {	    //	   it's allocated...
//			ogeom->delObj(BitVisualize_id[i]);//	      DELETE
//			geom_allocated[i]=0;	    //	      MARK as UNALLOCED
//		    }				    //	   not allocated
//		}				    //        NOTHING TO DO!
//	    }
	}
	vol_render1_grid();
	for (i=0; i<NUM_MATERIALS; i++) {
	    if (calc_mat[i]) {
		if (groups[i]->size()) {
		    GeomObj* topobj=scinew GeomMaterial(groups[i], matls[i]);
		    clString nm = "Material " + to_string(i+1);
		    BitVisualize_id[i] = ogeom->addObj(topobj, nm);
		} else {
		    geom_allocated[i]=0;
		}
	    }
	}
	ogeom->flushViews();
    } else {
	error("I can't BitVisualize this type of field...");
    }
}

static void do_parallel_vol_render1(void* obj, int proc)
{
  BitVisualize* module=(BitVisualize*)obj;
  module->parallel_vol_render1(proc);
}

void BitVisualize::vol_render1()
{
    np=Task::nprocessors();
    Task::multiprocess(np, do_parallel_vol_render1, this);
}


void BitVisualize::parallel_vol_render1(int proc) {
    int nx=sfrg->nx;
    int ny=sfrg->ny;
    int nz=sfrg->nz;
    double dx=(maxP.x()-minP.x())/(nx-1);
    double dy=(maxP.y()-minP.y())/(ny-1);
    double dz=(maxP.z()-minP.z())/(nz-1);
    int sx=proc*(nx-1)/np;
    int ex=(proc+1)*(nx-1)/np;
    Array1<GeomPts*> geomPts(NUM_MATERIALS);
    for (int m=0; m<NUM_MATERIALS; m++) {
	int bit=1<<m;
	if (calc_mat[m]) {
	    geomPts[m]=scinew GeomPts(10000, Vector(0,0,1));
	    double ii=minP.x()+sx*dx;
	    for (int i=sx; i<ex; i++, ii+=dx) {
		if (proc==0) update_progress(i,nx);
		double jj=minP.y();
		for (int j=0; j<ny; j++, jj+=dy) {
		    double kk=minP.z();
		    for (int k=0; k<nz; k++, kk+=dz) {
			int val=(int)(sfrg->grid(i,j,k));
			if ((val & bit) != 0) {
			    geomPts[m]->add(Point(ii,jj,kk));
			}
		    }
		}
	    }
	}	
    }
    grouplock.lock();
    for (m=0; m<NUM_MATERIALS; m++) {
	if (calc_mat[m]) {
	    if (geomPts[m]->pts.size())
		groups[m]->add(geomPts[m]);
//	    cerr << "adding material "<<m+1<<" with "<<geomPts[m]->pts.size()<<" points.\n";
	}
    }
    grouplock.unlock();
}

void BitVisualize::vol_render1_grid() {
    sfrg->get_bounds(minP, maxP);
    vol_render1();
}

#ifdef __GNUG__

#include <Classlib/Array1.cc>
template class Array1<SurfaceOPort*>;
template class Array1<SurfaceHandle>;
template class Array1<TriSurface*>;
template class Array1<TCLint*>;
template class Array1<GeomPts*>;
template class Array1<int>;
template class Array1<MaterialHandle>;
template class Array1<GeomGroup*>;

#endif

