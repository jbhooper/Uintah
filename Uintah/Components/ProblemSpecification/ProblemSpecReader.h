#ifndef __PROBLEM_SPEC_READER_H__ 
#define __PROBLEM_SPEC_READER_H__

#include <Uintah/Interface/ProblemSpecInterface.h>
#include <string>

using Uintah::Interface::ProblemSpecP;
using Uintah::Interface::ProblemSpec;
using Uintah::Interface::ProblemSpecInterface;

class ProblemSpecReader : public ProblemSpecInterface {

public:
  ProblemSpecReader(const std::string& filename);
  ~ProblemSpecReader();
   
  virtual ProblemSpecP readInputFile();

private:
  ProblemSpecReader(const ProblemSpecReader&);
  ProblemSpecReader& operator=(const ProblemSpecReader&);
  
  std::string filename;
};

#endif // __PROBLEM_SPEC_READER_H__
