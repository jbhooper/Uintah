#ifndef UINTAH_HOMEBREW_NCVARIABLE_H
#define UINTAH_HOMEBREW_NCVARIABLE_H

#include <Uintah/Grid/Array3.h>
#include <Uintah/Grid/NCVariableBase.h>
#include <Uintah/Grid/TypeDescription.h>
#include <SCICore/Exceptions/InternalError.h>
#include <SCICore/Geometry/Vector.h>
#include <Uintah/Exceptions/TypeMismatchException.h>
#include <Uintah/Grid/Region.h>

using namespace Uintah;

namespace Uintah {
   using SCICore::Exceptions::InternalError;
   using SCICore::Geometry::Vector;

   class TypeDescription;

/**************************************

CLASS
   NCVariable
   
GENERAL INFORMATION

   NCVariable.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
  
   Copyright (C) 2000 SCI Group

KEYWORDS
   NCVariable

DESCRIPTION
   Long description...
  
WARNING
  
****************************************/

   template<class T> class NCVariable : public Array3<T>, public NCVariableBase{
   public:
     
     NCVariable();
     NCVariable(const NCVariable<T>&);
     virtual ~NCVariable();
     
     //////////
     // Insert Documentation Here:
     static const TypeDescription* getTypeDescription();
     
     virtual void copyPointer(const NCVariableBase&);
     
     //////////
     // Insert Documentation Here:
     virtual NCVariable<T>* clone() const;
     
     //////////
     // Insert Documentation Here:
     virtual void allocate(const Region*);
     
     NCVariable<T>& operator=(const NCVariable<T>&);
     
     // Replace the values on the indicated face with value
     void fillFace(Region::FaceType face, Vector value) {};
     
     // Use to apply symmetry boundary conditions.  On the
     // indicated face, replace the component of the vector
     // normal to the face with 0.0
     void fillFaceNormal(Region::FaceType) {};
     
   private:
   };
   
   template<class T>
      const TypeDescription*
      NCVariable<T>::getTypeDescription()
      {
	 static TypeDescription* td;
	 if(!td)
	    td = new TypeDescription(false, TypeDescription::Node);
	 return td;
      }
   
   template<class T>
      NCVariable<T>::~NCVariable()
      {
      }
   
   template<class T>
      NCVariable<T>*
      NCVariable<T>::clone() const
      {
	 return new NCVariable<T>(*this);
      }
   
   template<class T>
      void
      NCVariable<T>::copyPointer(const NCVariableBase& copy)
      {
	 const NCVariable<T>* c = dynamic_cast<const NCVariable<T>* >(&copy);
	 if(!c)
	    throw TypeMismatchException("Type mismatch in NC variable");
	 *this = *c;
      }

   template<class T>
      NCVariable<T>&
      NCVariable<T>::operator=(const NCVariable<T>& copy)
      {
	 if(this != &copy){
	    Array3<T>::operator=(copy);
	 }
	 return *this;
      }
   
   template<class T>
      NCVariable<T>::NCVariable()
      {
      }
   
   template<class T>
      NCVariable<T>::NCVariable(const NCVariable<T>& copy)
      : Array3<T>(copy)
      {
      }
   
   template<class T>
      void
      NCVariable<T>::allocate(const Region* region)
      {
	 if(getWindow())
	    throw InternalError("Allocating an NCvariable that is apparently already allocated!");
	 IntVector res(region->getNCells());
	 resize(res.x()+1, res.y()+1, res.z()+1);
      }
   
} // end namespace Uintah

//
// $Log$
// Revision 1.12  2000/05/09 03:24:39  jas
// Added some enums for grid boundary conditions.
//
// Revision 1.11  2000/05/07 06:02:12  sparker
// Added beginnings of multiple patch support and real dependencies
//  for the scheduler
//
// Revision 1.10  2000/05/04 19:06:47  guilkey
// Added the beginnings of grid boundary conditions.  Functions still
// need to be filled in.
//
// Revision 1.9  2000/05/02 06:07:21  sparker
// Implemented more of DataWarehouse and SerialMPM
//
// Revision 1.8  2000/04/26 06:48:49  sparker
// Streamlined namespaces
//
// Revision 1.7  2000/04/20 18:56:30  sparker
// Updates to MPM
//
// Revision 1.6  2000/04/12 23:00:48  sparker
// Starting problem setup code
// Other compilation fixes
//
// Revision 1.5  2000/04/11 07:10:50  sparker
// Completing initialization and problem setup
// Finishing Exception modifications
//
// Revision 1.4  2000/03/21 02:22:57  dav
// few more updates to make it compile including moving Array3 stuff out of namespace as I do not know where it should be
//
// Revision 1.3  2000/03/16 22:07:59  dav
// Added the beginnings of cocoon docs.  Added namespaces.  Did a few other coding standards updates too
//
//

#endif
