
#ifndef TRI_H
#define TRI_H 1

#include <Packages/rtrt/Core/Object.h>
#include <Core/Geometry/Point.h>

namespace rtrt {

class Tri : public Object {
    Point p1, p2, p3;
    Vector n;
    double d;
    Vector e1p, e2p, e3p;
    Vector e1, e2, e3;
    double e1l, e2l, e3l;
    bool bad;
public:
    inline bool isbad() {
	return bad;
    }
    Tri(Material* matl, const Point& p1, const Point& p2, const Point& p3);
    virtual ~Tri();
    virtual void intersect(const Ray& ray, HitInfo& hit, DepthStats* st,
			   PerProcessorContext*);
    virtual void light_intersect(Light* light, const Ray& ray,
				 HitInfo& hit, double dist, Color& atten,
				 DepthStats* st, PerProcessorContext*);
    virtual Vector normal(const Point&, const HitInfo& hit);
    inline Vector normal() { return n; }
    virtual void compute_bounds(BBox&, double offset);

    Point centroid()
    {
	double one_third = 1./3.;

	return AffineCombination(p1,one_third,
				 p2,one_third,
				 p3,one_third);
    }
    Point pt(int i)
    {
	if (i==0)
	    return p1;
	else if (i==1)
	    return p2;
	else 
	    return p3;
    }
	       
};

} // end namespace rtrt

#endif
