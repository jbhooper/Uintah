#ifndef UINTAH_HOMEBREW_SimpleLoadBalancer_H
#define UINTAH_HOMEBREW_SimpleLoadBalancer_H

#include <Uintah/Interface/LoadBalancer.h>
#include <Uintah/Parallel/UintahParallelComponent.h>

namespace Uintah {
   
   /**************************************
     
     CLASS
       SimpleLoadBalancer
      
       Short Description...
      
     GENERAL INFORMATION
      
       SimpleLoadBalancer.h
      
       Steven G. Parker
       Department of Computer Science
       University of Utah
      
       Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
      
       Copyright (C) 2000 SCI Group
      
     KEYWORDS
       SimpleLoadBalancer
      
     DESCRIPTION
       Long description...
      
     WARNING
      
     ****************************************/
    
   class SimpleLoadBalancer : public LoadBalancer, public UintahParallelComponent {
   public:
      SimpleLoadBalancer(const ProcessorGroup* myworld);
      ~SimpleLoadBalancer();
      virtual void assignResources(TaskGraph& tg, const ProcessorGroup*);
      virtual int getPatchwiseProcessorAssignment(const Patch* patch,
						  const ProcessorGroup* resources);
   private:
      SimpleLoadBalancer(const SimpleLoadBalancer&);
      SimpleLoadBalancer& operator=(const SimpleLoadBalancer&);
      
   };
} // end namespace Uintah

//
// $Log$
// Revision 1.1  2000/09/20 16:00:28  sparker
// Added external interface to LoadBalancer (for per-processor tasks)
// Added message logging functionality. Put the tag <MessageLog/> in
//    the ups file to enable
//
//

#endif

