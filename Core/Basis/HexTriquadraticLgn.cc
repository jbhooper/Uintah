//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2004 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//  
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//  
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//  
//    File   : HexTriquadraticLgn.cc
//    Author : Martin Cole, Frank B. Sachse
//    Date   : Dec 04 2004

#include <Core/Basis/HexTriquadraticLgn.h>
 
namespace SCIRun {

double HexTriquadraticLgnUnitElement::unit_vertices[20][3] = {
  {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, 
  {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1},
  {0.5,0,0}, {1,0.5,0}, {.5,1,0}, {0,.5,0}, 
  {0,0,0.5}, {1,0,.5}, {0,1,.5}, {1,1,.5}, 
  {0.5,0,1}, {1,0.5,1}, {.5,1,1}, {0,.5,1}, 

};
int HexTriquadraticLgnUnitElement::unit_edges[12][2] = {{0,1}, {1,2}, {2,3}, {3,0},
				   {0,4}, {1,5}, {2,6}, {3,7},
				   {4,5}, {5,6}, {6,7}, {7,4}};
int HexTriquadraticLgnUnitElement::unit_faces[6][4] = {{0,1,2,3}, {0,1,5,4}, {1,2,6,5},
				  {2,3,7,6}, {3,0,4,7}, {4,5,6,7}};


} //namespace SCIRun

