#ifndef SPECFEM2VTU_HPP
#define SPECFEM2VTU_HPP

#include <string>



class Parameters;



void specfem2vtu(const Parameters &param);

void read_dofs(const std::string &filename);

#endif // SPECFEM2VTU_HPP
