#include "specfem2vtu.hpp"
#include "parameters.hpp"
#include "point.hpp"

#include <vector>
#include <stdexcept>
#include <fstream>


void specfem2vtu(const Parameters &param)
{
  std::vector<Point> dofs; // coordinates of the degrees of freedom
  read_dofs(param._file_dofs);
}




void read_dofs(const std::string &filename)
{
  std::ifstream in(filename.c_str());
  if (!in)
  {
    throw std::runtime_error("File " + filename + " can't be opened");
  }
}





