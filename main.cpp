#include "parameters.hpp"
#include "specfem2vtu.hpp"




int main(int argc, char **argv)
{
  Parameters param(argc, argv);

  if (param._verbose > 1)
    param.print_parameters();

  param.check_parameters();

  specfem2vtu(param);

  return 0;
}

