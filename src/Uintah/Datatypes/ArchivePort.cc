
#include <Uintah/Datatypes/ArchivePort.h>
#include <Uintah/share/share.h>

namespace PSECore {
namespace Datatypes {

using namespace Uintah::Datatypes;


extern "C" {
UINTAHSHARE IPort* make_ArchiveIPort(Module* module, const clString& name) {
  return new SimpleIPort<ArchiveHandle>(module,name);
}
UINTAHSHARE OPort* make_ArchiveOPort(Module* module, const clString& name) {
  return new SimpleOPort<ArchiveHandle>(module,name);
}
}

template<> clString SimpleIPort<ArchiveHandle>::port_type("Archive");
template<> clString SimpleIPort<ArchiveHandle>::port_color("lightsteelblue4");


} // End namespace Datatypes
} // End namespace PSECore

