#ifndef UINTAH_HOMEBREW_REFCOUNTED_H
#define UINTAH_HOMEBREW_REFCOUNTED_H

namespace Uintah {

/**************************************

CLASS
   RefCounted
   
   Short description...

GENERAL INFORMATION

   Task.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
  
   Copyright (C) 2000 SCI Group

KEYWORDS
   Reference_Counted, RefCounted, Handle

DESCRIPTION
   Long description...
  
WARNING
  
****************************************/

   class RefCounted {
   public:
      RefCounted();
      RefCounted(const RefCounted&);
      RefCounted& operator=(const RefCounted&);
      virtual ~RefCounted();
      
      //////////
      // Insert Documentation Here:
      void addReference();
      
      //////////
      // Insert Documentation Here:
      bool removeReference();

      int getReferenceCount() const {
	 return d_refCount;
      }
   private:
      //////////
      // Insert Documentation Here:
      int d_refCount;
      int d_lockIndex;
   };
   
} // end namespace Uintah

//
// $Log$
// Revision 1.4  2000/08/21 23:27:07  sparker
// Added getReferenceCount() method to RefCounted
// Correctly maintain ref counts on neighboring particle subsets in ParticleSubset
//
// Revision 1.3  2000/04/26 06:48:53  sparker
// Streamlined namespaces
//
// Revision 1.2  2000/03/16 22:08:00  dav
// Added the beginnings of cocoon docs.  Added namespaces.  Did a few other coding standards updates too
//
//

#endif

