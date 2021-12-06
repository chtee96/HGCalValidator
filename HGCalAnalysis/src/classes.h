#ifdef __CINT_
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;
#pragma link C++ class vector <set<unsigned int> > + ;
#pragma link C++ class vector <vector<bool> > + ;
#endif /* __CINT__ */

#include <vector>
#include <set>

std::vector<std::set<unsigned int> > vsui;
std::vector<std::vector<bool> > vvb;
