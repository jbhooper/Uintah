// 
// File:          who_Com.hh
// Symbol:        who.Com-v1.0
// Symbol Type:   class
// Babel Version: 0.7.0
// SIDL Created:  20020730 17:57:40 MST
// Generated:     20020730 17:57:41 MST
// Description:   Client-side glue code for who.Com
// 
// WARNING: Automatically generated; changes will be lost
// 
// 

#ifndef included_who_Com_hh
#define included_who_Com_hh

// declare class before #includes
// (this alleviates circular #include guard problems)[BUG#393]
namespace who { 

  class Com;
} // end namespace who

// Some compilers need to define array template before the specializations
#ifndef included_SIDL_cxx_hh
#include "SIDL_cxx.hh"
#endif
namespace SIDL {
  template<>
  class array<who::Com>;
}

#ifndef included_SIDL_cxx_hh
#include "SIDL_cxx.hh"
#endif
#ifndef included_who_Com_IOR_h
#include "who_Com_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_SIDL_BaseInterface_hh
#include "SIDL_BaseInterface.hh"
#endif
#ifndef included_govcca_Services_hh
#include "govcca_Services.hh"
#endif

namespace who { 

  /**
   * Symbol "who.Com" (version 1.0)
   */
  class Com : public SIDL::StubBase {

    //////////////////////////////////////////////////
    // 
    // User Defined Methods
    // 

  public:

    /**
     * <p>
     * Add one to the intrinsic reference count in the underlying object.
     * Object in <code>SIDL</code> have an intrinsic reference count.
     * Objects continue to exist as long as the reference count is
     * positive. Clients should call this method whenever they
     * create another ongoing reference to an object or interface.
     * </p>
     * <p>
     * This does not have a return value because there is no language
     * independent type that can refer to an interface or a
     * class.
     * </p>
     */
    void
    addReference() throw ( SIDL::NullIORException ) 
    ;


    /**
     * Decrease by one the intrinsic reference count in the underlying
     * object, and delete the object if the reference is non-positive.
     * Objects in <code>SIDL</code> have an intrinsic reference count.
     * Clients should call this method whenever they remove a
     * reference to an object or interface.
     */
    void
    deleteReference() throw ( SIDL::NullIORException ) 
    ;


    /**
     * Return true if and only if <code>obj</code> refers to the same
     * object as this object.
     */
    bool
    isSame (
      /*in*/ SIDL::BaseInterface iobj
    )
    throw ( SIDL::NullIORException ) 
    ;



    /**
     * Check whether the object can support the specified interface or
     * class.  If the <code>SIDL</code> type name in <code>name</code>
     * is supported, then a reference to that object is returned with the
     * reference count incremented.  The callee will be responsible for
     * calling <code>deleteReference</code> on the returned object.  If
     * the specified type is not supported, then a null reference is
     * returned.
     */
    SIDL::BaseInterface
    queryInterface (
      /*in*/ std::string name
    )
    throw ( SIDL::NullIORException ) 
    ;



    /**
     * Return whether this object is an instance of the specified type.
     * The string name must be the <code>SIDL</code> type name.  This
     * routine will return <code>true</code> if and only if a cast to
     * the string type name would succeed.
     */
    bool
    isInstanceOf (
      /*in*/ std::string name
    )
    throw ( SIDL::NullIORException ) 
    ;



    /**
     * Obtain Services handle, through which the component communicates with the
     * framework. This is the one method that every CCA Component must implement. 
     */
    void
    setServices (
      /*in*/ govcca::Services svc
    )
    throw ( SIDL::NullIORException ) 
    ;



    //////////////////////////////////////////////////
    // 
    // End User Defined Methods
    // (everything else in this file is specific to
    //  Babel's C++ bindings)
    // 

  public:
    typedef struct who_Com__object ior_t;
    typedef struct who_Com__external ext_t;
    typedef struct who_Com__sepv sepv_t;

    // default constructor
    Com() : d_self(0), d_weak_reference(false) { }

    // static constructor
    static who::Com _create();

    // default destructor
    virtual ~Com ();

    // copy constructor
    Com ( const Com& original );

    // assignment operator
    Com& operator= ( const Com& rhs );

    // conversion from ior to C++ class
    Com ( Com::ior_t* ior );

    // Alternate constructor: does not call addReference()
    // (sets d_weak_reference=isWeak)
    // For internal use by Impls (fixes bug#275)
    Com ( Com::ior_t* ior, bool isWeak );

    // conversion from a StubBase
    Com ( const SIDL::StubBase& base );

    ior_t* _get_ior() { return d_self; }

    const ior_t* _get_ior() const { return d_self; }

    void _set_ior( ior_t* ptr ) { d_self = ptr; }

    bool _is_nil() const { return (d_self==0); }

    bool _not_nil() const { return (d_self!=0); }

    bool operator !() const { return (d_self==0); }

  protected:
    virtual void* _cast(const char* type) const;

  private:
    // Pointer to SIDL's IOR type (one per instance)
    ior_t * d_self;

    // Weak references (used by Impl's only) don't add/deleteRef()
    bool d_weak_reference;

    // Pointer to external (DLL loadable) symbols (shared among instances)
    static const ext_t * s_ext;

  public:
    static const ext_t * _get_ext() throw ( SIDL::NullIORException );

  }; // end class Com
} // end namespace who

// array specialization
namespace SIDL {
  template<>
  class array<who::Com> : public array_mixin< struct who_Com__array,
    struct who_Com__object*, who::Com>{
  public:
    /**
     * default constructor
     */
    array() : array_mixin< struct who_Com__array, struct who_Com__object*,
      who::Com>() {}

    /**
     * default destructor
     */
    virtual ~array() { if ( d_array ) { deleteReference(); } }

    /**
     * copy constructor
     */
    array( const array< who::Com >& original ) {
      d_array = original.d_array;
      if ( d_array ) { addReference(); }
    }

    /**
     * assignment operator
     */
    array< who::Com >& operator=( const array< who::Com >& rhs ) {
      if ( d_array != rhs.d_array ) {
        if ( d_array ) { deleteReference(); }
        d_array=rhs.d_array;
        if ( d_array ) { addReference(); }
      }
      return *this;
    }

  /**
   * conversion from ior to C++ class
   * (constructor/casting operator)
   */
  array( struct who_Com__array* src ) : array_mixin< struct who_Com__array,
    struct who_Com__object*, who::Com>(src) {}

  /**
   * static constructor: createRow
   */
  static array < who::Com >
  createRow(int32_t dimen, const int32_t lower[], const int32_t upper[]) {
    array < who::Com > a;
    a._set_ior( (struct who_Com__array*)
                SIDL_interface__array_createRow( dimen, lower, upper ) );
    return a;
  }

  /**
   * static constructor: createCol
   */
  static array < who::Com >
  createCol(int32_t dimen, const int32_t lower[], const int32_t upper[]) {
    array < who::Com > a;
    a._set_ior( (struct who_Com__array*)
                SIDL_interface__array_createCol( dimen, lower, upper ) );
    return a;
  }

  /**
   * static constructor: create1d
   */
  static array < who::Com >
  create1d(int32_t len) {
    array < who::Com > a;
    a._set_ior( (struct who_Com__array*)
                SIDL_interface__array_create1d(len));
    return a;
  }

  /**
   * static constructor: create2dCol
   */
  static array < who::Com >
  create2dCol(int32_t m, int32_t n) {
    array < who::Com > a;
    a._set_ior( (struct who_Com__array*)
                SIDL_interface__array_create2dCol(m,n));
    return a;
  }

  /**
   * static constructor: create2dRow
   */
  static array < who::Com >
  create2dRow(int32_t m, int32_t n) {
    array < who::Com > a;
    a._set_ior( (struct who_Com__array*)
                SIDL_interface__array_create2dRow(m,n));
    return a;
  }

  /**
   * constructor: slice
   */
  array < who::Com >
  slice( int dimen,
         const int32_t newElem[],
         const int32_t *srcStart = 0,
         const int32_t *srcStride = 0,
         const int32_t *newStart = 0) {
    array < who::Com > a;
    a._set_ior( (struct who_Com__array*)
                SIDL_interface__array_slice( (struct SIDL_interface__array *) 
      d_array,
                                   dimen, newElem, srcStart, srcStride,
      newStart) );
    return a;
  }

  /**
   * copy
   */
  void copy( const array<who::Com>& src ) {
    SIDL_interface__array_copy( (const struct SIDL_interface__array*) 
      src._get_ior(),(struct SIDL_interface__array*) _get_ior() );
  }

  /**
   * constructor: smart copy
   */
  void smartCopy() {
    struct SIDL_interface__array * p  = SIDL_interface__array_smartCopy(
      (struct SIDL_interface__array *) _get_ior() );
    if ( _not_nil() ) { deleteReference(); }
    _set_ior( (struct who_Com__array*) p);
  }

  /**
   * constructor: ensure
   */
  void ensure(int32_t dimen, int ordering ) {
    struct SIDL_interface__array * p  = SIDL_interface__array_ensure(
      (struct SIDL_interface__array *) d_array, dimen, ordering );
    if ( _not_nil() ) { deleteReference(); }
    _set_ior( (struct who_Com__array*) p);
  }
};
}

#endif
