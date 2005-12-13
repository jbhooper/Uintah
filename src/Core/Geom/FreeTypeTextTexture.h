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
 *  FreeTypeTextTexture.h
 *
 *  Written by:
 *   McKay Davis
 *   School of Computing
 *   University of Utah
 *   November, 2005
 *
 *  Copyright (C) 2005 SCI Group
 */

#ifndef SCIRun_Dataflow_Modules_Render_FreeTypeTextTexture_h
#define SCIRun_Dataflow_Modules_Render_FreeTypeTextTexture_h

#include <Core/Geom/FreeType.h>

namespace SCIRun {

class NrrdTextureObj;

class FreeTypeTextTexture {
public:
  FreeTypeTextTexture(const string &text, FreeTypeFace *face);
  ~FreeTypeTextTexture();
  enum anchor_e { n, e, s, w, ne, se, sw, nw, c };
  void			set(const string &text);
  void			draw(double x, double y, anchor_e anchor = sw);
  int                   width();
  int                   height();
private:
  void			render_text_to_texture();
  NrrdTextureObj *	texture_;
  string		text_;
  FreeTypeFace *	face_;
  bool			dirty_;
};
  
}


#endif
