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
 *  ColorMapReader.cc: Read a persistent colormap from a file
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   July 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#include <Dataflow/Ports/ColorMapPort.h>
#include <Dataflow/Modules/DataIO/GenericReader.h>
#include <Core/ImportExport/ColorMap/ColorMapIEPlugin.h>

namespace SCIRun {

template class GenericReader<ColorMapHandle>;

class ColorMapReader : public GenericReader<ColorMapHandle> {
protected:
  virtual bool call_importer(const string &filename);

public:
  ColorMapReader(GuiContext* ctx);
};

DECLARE_MAKER(ColorMapReader)

ColorMapReader::ColorMapReader(GuiContext* ctx)
  : GenericReader<ColorMapHandle>("ColorMapReader", ctx, "DataIO", "SCIRun")
{
}


bool
ColorMapReader::call_importer(const string &filename)
{
  ColorMapIEPluginManager mgr;
  ColorMapIEPlugin *pl = mgr.get_plugin("SomePlugin");
  if (pl)
  {
    handle_ = pl->filereader(this, filename.c_str());
    return handle_.get_rep();
  }
  return false;
}


} // End namespace SCIRun
