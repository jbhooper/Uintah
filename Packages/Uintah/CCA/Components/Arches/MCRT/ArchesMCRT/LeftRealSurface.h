#ifndef LeftRealSurface_H
#define LeftRealSurface_H

#include "RealSurface.h"

class LeftRealSurface:public RealSurface{
public:
  //LeftRealSurface(int _surfaceIndex, int _xno);
  LeftRealSurface();
  ~LeftRealSurface();
  void setData(int _surfaceIndex, int _xno);
  virtual void set_n(double *nn);
  virtual void get_n();
  virtual void get_t1();
  virtual void get_t2();
  virtual void get_limits(double *VolTable, int &vIndex);
private:
  //double n[3], t1[3], t2[3];   
  double ylow, yup, zlow, zup;
  double xleft;
  int VolIndex;
  int xno;
};

#endif
  
