/*
 * File:          whoc_IDPort_Impl.c
 * Symbol:        whoc.IDPort-v1.0
 * Symbol Type:   class
 * Babel Version: 0.7.4
 * SIDL Created:  20021110 23:39:22 MST
 * Generated:     20021110 23:39:25 MST
 * Description:   Server-side implementation for whoc.IDPort
 * 
 * WARNING: Automatically generated; only changes within splicers preserved
 * 
 * babel-version = 0.7.4
 * source-line   = 7
 * source-url    = file:/home/sparker/SCIRun/src/Babel/Components/whoc/whoc.sidl
 */

/*
 * DEVELOPERS ARE EXPECTED TO PROVIDE IMPLEMENTATIONS
 * FOR THE FOLLOWING METHODS BETWEEN SPLICER PAIRS.
 */

/*
 * Symbol "whoc.IDPort" (version 1.0)
 */

#include "whoc_IDPort_Impl.h"

/* DO-NOT-DELETE splicer.begin(whoc.IDPort._includes) */
/* Put additional includes or other arbitrary code here... */
#include "string.h"
/* DO-NOT-DELETE splicer.end(whoc.IDPort._includes) */

/*
 * Class constructor called when the class is created.
 */

#undef __FUNC__
#define __FUNC__ "impl_whoc_IDPort__ctor"

void
impl_whoc_IDPort__ctor(
  whoc_IDPort self)
{
  /* DO-NOT-DELETE splicer.begin(whoc.IDPort._ctor) */
  /* Insert the implementation of the constructor method here... */
  /* DO-NOT-DELETE splicer.end(whoc.IDPort._ctor) */
}

/*
 * Class destructor called when the class is deleted.
 */

#undef __FUNC__
#define __FUNC__ "impl_whoc_IDPort__dtor"

void
impl_whoc_IDPort__dtor(
  whoc_IDPort self)
{
  /* DO-NOT-DELETE splicer.begin(whoc.IDPort._dtor) */
  /* Insert the implementation of the destructor method here... */
  /* DO-NOT-DELETE splicer.end(whoc.IDPort._dtor) */
}

/*
 * Test prot. Return a string as an ID for Hello component
 */

#undef __FUNC__
#define __FUNC__ "impl_whoc_IDPort_getID"

char*
impl_whoc_IDPort_getID(
  whoc_IDPort self)
{
  /* DO-NOT-DELETE splicer.begin(whoc.IDPort.getID) */
  char *id="World (in C )";
  char *s=(char*)malloc(strlen(id)+1);
  strcpy(s, id);
  return s;
  /* DO-NOT-DELETE splicer.end(whoc.IDPort.getID) */
}
