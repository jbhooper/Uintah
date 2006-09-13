/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2004 Scientific Computing and Imaging Institute,
   University of Utah.

   License for the specific language governing rights and limitations under
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/


/*
 *  PortInstance.h:
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   October 2001
 *
 */

#ifndef Framework_PortInstance_h
#define Framework_PortInstance_h

#include <Core/CCA/spec/cca_sidl.h>
#include <string>

namespace SCIRun
{

/**
 * \class PortInstance
 *
 *
 *
 */
class PortInstance
{
public:
  PortInstance();
  virtual ~PortInstance();

  enum PortType {
    Uses = 0,
    Provides
  };

  /** */
  virtual bool connect(PortInstance*) = 0;

  /** */
  virtual PortType portType() = 0;

  /** */
  virtual std::string getType();

  /** */
  virtual std::string getModel();

  /** */
  virtual std::string getUniqueName() = 0;

  /** */
  virtual bool disconnect(PortInstance*) =0;

  /** */
  virtual bool canConnectTo(PortInstance *)=0;

  /** */
  virtual bool available();

  /** */
  virtual sci::cca::TypeMap::pointer getProperties() = 0;

  /** */
  virtual void setProperties(const sci::cca::TypeMap::pointer& tm) = 0;

  /** */
  virtual PortInstance* getPeer();

  // default keys
  static const std::string MAX_C;
  static const std::string MIN_C;
  static const std::string PROXY;
  static const std::string NAME;
  static const std::string TYPE;
  static const std::string MODEL;

private:
  PortInstance(const PortInstance&);
  PortInstance& operator=(const PortInstance&);
};

} // end namespace SCIRun

#endif
