// 
// File:          framework_Component_Impl.hh
// Symbol:        framework.Component-v1.0
// Symbol Type:   class
// Babel Version: 0.7.0
// SIDL Created:  20020730 13:51:30 MST
// Generated:     20020730 13:51:36 MST
// Description:   Server-side implementation for framework.Component
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 

#ifndef included_framework_Component_Impl_hh
#define included_framework_Component_Impl_hh

#ifndef included_SIDL_cxx_hh
#include "SIDL_cxx.hh"
#endif
#ifndef included_framework_Component_IOR_h
#include "framework_Component_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_SIDL_BaseInterface_hh
#include "SIDL_BaseInterface.hh"
#endif
#ifndef included_framework_Component_hh
#include "framework_Component.hh"
#endif
#ifndef included_govcca_Services_hh
#include "govcca_Services.hh"
#endif


// DO-NOT-DELETE splicer.begin(framework.Component._includes)
// Put additional includes or other arbitrary code here...
// DO-NOT-DELETE splicer.end(framework.Component._includes)

namespace framework { 

  /**
   * Symbol "framework.Component" (version 1.0)
   */
  class Component_impl
  // DO-NOT-DELETE splicer.begin(framework.Component._inherits)
  // Put additional inheritance here...
  // DO-NOT-DELETE splicer.end(framework.Component._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    Component self;

    // DO-NOT-DELETE splicer.begin(framework.Component._implementation)
    // Put additional implementation details here...
    // DO-NOT-DELETE splicer.end(framework.Component._implementation)

  private:
    // private default constructor (required)
    Component_impl() {} 

  public:
    // SIDL constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    Component_impl( struct framework_Component__object * s ) : self(s,
      true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~Component_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

  public:


    /**
     * Obtain Services handle, through which the component communicates with the
     * framework. This is the one method that every CCA Component must implement. 
     */
    void
    setServices (
      /*in*/ govcca::Services svc
    )
    throw () 
    ;

  };  // end class Component_impl

} // end namespace framework

// DO-NOT-DELETE splicer.begin(framework.Component._misc)
// Put miscellaneous things here...
// DO-NOT-DELETE splicer.end(framework.Component._misc)

#endif
