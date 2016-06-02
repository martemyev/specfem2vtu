#include "parameters.hpp"
#include "specfem2vtu.hpp"

int main(int argc, char **argv)
{
  try
  {
    Parameters param(argc, argv);

    if (param._verbose > 1)
      param.print_parameters();

    param.check_parameters();

    specfem2vtu(param);
  }
  catch(const std::runtime_error &e)
  {
    std::cerr << "\n\nRuntime error caught:\n" << e.what() << std::endl;
    return 1;
  }
  catch(const std::exception &e)
  {
    std::cerr << "\n\nException caught:\n" << e.what() << std::endl;
    return 1;
  }
  catch(...)
  {
    std::cerr << "\n\nUnknown exception was caught\n" << std::endl;
    return 1;
  }

  return 0;
}

