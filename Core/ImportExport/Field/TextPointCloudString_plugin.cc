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
 *  TextPointCloudField_plugin.cc
 *
 *  Written by:
 *   Michael Callahan
 *   Department of Computer Science
 *   University of Utah
 *   December 2004
 *
 *  Copyright (C) 2004 SCI Group
 */

// Read in a .pcs.txt file.  The file should contain one entry per
// line where each entry specifies the x/y/z position separated by
// spaces and then the string.  Example:
//
// 0.0 0.0 0.0 Origin
// 1.0 0.0 0.0 Positive X Axis
// 0.0 1.0 0.0 Positive Y Axis
// 0.0 0.0 1.0 Positive Z Axis
// 1.2 0.5 0.2 Point of Interest


#include <Core/Datatypes/PointCloudField.h>
#include <Core/ImportExport/Field/FieldIEPlugin.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

using namespace std;
using namespace SCIRun;


static FieldHandle
reader(const char *filename)
{
  ifstream ptsstream(filename);

  PointCloudMesh *pcm = scinew PointCloudMesh();
  vector<string> strings;
  char buffer[1024];
  double x, y, z;
  
  while (!ptsstream.eof())
  {
      ptsstream >> x >> y >> z;
    pcm->add_point(Point(x, y, z));
    
    char c;
    ptsstream.get(c); // Eat up the one whitespace.
    ptsstream.getline(buffer, 1024);
    strings.push_back(string(buffer));
  }

  PointCloudField<string> *pc = 
    scinew PointCloudField<string>(pcm, Field::NODE);
  for (unsigned int i=0; i < strings.size(); i++)
  {
    pc->set_value(strings[i], PointCloudMesh::Node::index_type(i));
  }

  return FieldHandle(pc);
}    


static void
writer(FieldHandle field, const char *filename)
{
  ofstream ptsstream(filename);

  PointCloudField<string> *pcfs =
    dynamic_cast<PointCloudField<string> *>(field.get_rep());
  if (pcfs == 0)
  {
    // Handle error checking somehow.
    return;
  }
  PointCloudMeshHandle mesh = pcfs->get_typed_mesh();
  
  PointCloudMesh::Node::iterator itr, eitr;
  mesh->begin(itr);
  mesh->end(eitr);

  while (itr != eitr)
  {
    Point c;
    mesh->get_center(c, *itr);
    
    ptsstream << c.x() << " " << c.y() << " " << c.z() << " " 
	      << pcfs->value(*itr) << "\n";
    
    ++itr;
  }
}



static FieldIEPlugin
TextPointCloudString_plugin("TextPointCloudString",
			    "pcs.txt", "", reader, writer);
