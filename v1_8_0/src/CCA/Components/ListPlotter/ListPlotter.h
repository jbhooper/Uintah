/*
  The contents of this file are subject to the University of Utah Public
  License (the "License"); you may not use this file except in compliance
  with the License.
  
  Software distributed under the License is distributed on an "AS IS"
  basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
  License for the specific language governing rights and limitations under
  the License.
  
  The Original Source Code is SCIRun, released March 12, 2001.
  
  The Original Source Code was developed by the University of Utah.
  Portions created by UNIVERSITY are Copyright (C) 2001, 1994 
  University of Utah. All Rights Reserved.
*/

/*
 *  ListPlotter.h
 *
 *  Written by:
 *   Keming Zhang 
 *   Department of Computer Science
 *   University of Utah
 *   March 2002
 *
 */

#ifndef SCIRun_Framework_ListPlotter_h
#define SCIRun_Framework_ListPlotter_h

#include <Core/CCA/spec/cca_sidl.h>
#include "ListPlotterForm.h"
//namespace SCIRun {
  
  class ImUIPort : public virtual gov::cca::ports::UIPort {

  public:
    void setServices(const gov::cca::Services::pointer &svc);
    virtual ~ImUIPort(){}
    virtual int ui();
  private:
    gov::cca::Services::pointer services;
  };

class ListPlotter : public gov::cca::Component{
                
  public:
    ListPlotter();
    virtual ~ListPlotter();

    virtual void setServices(const gov::cca::Services::pointer& svc);
  private:

    ListPlotter(const ListPlotter&);
    ListPlotter& operator=(const ListPlotter&);
    ImUIPort ui;
    gov::cca::Services::pointer services;
  };
//}




#endif
