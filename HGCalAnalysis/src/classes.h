#ifdef __CINT_
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;
#pragma link C++ class Valid + ;
#pragma link C++ class vector < Valid> + ;
#pragma link C++ class vector <set<unsigned int> > + ;
#pragma link C++ class vector <vector<bool> > + ;
#endif /* __CINT__ */

#include "HGCalValidator/HGCalAnalysis/interface/Valid.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include <vector>
#include <set>

Valid vh;
std::vector<Valid> vvh;
edm::Wrapper<std::vector<Valid> > wvvh;

std::vector<std::set<unsigned int> > vsui;
std::vector<std::vector<bool> > vvb;
