#include "parameters.hpp"
#include "utilities.hpp"

#include <algorithm>



//==============================================================================
//
// Parameters class
//
//==============================================================================

std::string Parameters::DEFAULT_FILE_NAME = "no-file";
int Parameters::DEFAULT_PRINT_LEN = 10;




Parameters::Parameters(int argc, char **argv)
  : _file_dofs(DEFAULT_FILE_NAME),
    _file_solution_base(DEFAULT_FILE_NAME),
    _n_time_steps(0),
    _step_snapshot(0),
    _binary(false),
    _verbose(2),
    _longest_string_key_len(DEFAULT_PRINT_LEN),
    _longest_string_value_len(DEFAULT_PRINT_LEN)
{
  _parameters["-fdofs"] = std::unique_ptr<ParamBase>(new OneParam<std::string>("name of file with degrees of freedom", &_file_dofs));
  _parameters["-fbase"] = std::unique_ptr<ParamBase>(new OneParam<std::string>("base name of solution files", &_file_solution_base));
  _parameters["-nt"]    = std::unique_ptr<ParamBase>(new OneParam<int>("number of time steps", &_n_time_steps));
  _parameters["-step"]  = std::unique_ptr<ParamBase>(new OneParam<int>("step in snapshot export", &_step_snapshot));
  _parameters["-bin"]   = std::unique_ptr<ParamBase>(new OneParam<bool>("if the files are in binary format", &_binary));
  _parameters["-v"]     = std::unique_ptr<ParamBase>(new OneParam<int>("verbosity level", &_verbose));

  update_longest_string_key_len();

  if (argc == 1 || argcheck(argc, argv, "-help") || argcheck(argc, argv, "-h"))
  {
    print_options();
    exit(0);
  }

  read_command_line(argc, argv);

  update_longest_string_value_len();
}




void Parameters::read_command_line(int argc, char **argv)
{
  // all the command line entries
  std::vector<std::string> arguments(argv, argv+argc);

  // starting from the 1-st (not 0-th, because arguments[0] is the path to the
  // executable file of this program). Then we consider every second parameter,
  // since the command line goes like this: param0 value0 param1 value1 ...
  // and we need to consider only parameters.
  for (size_t ar = 1; ar < arguments.size(); ar += 2)
  {
    std::map<std::string, std::unique_ptr<ParamBase> >::const_iterator iter = _parameters.find(arguments[ar]);
    if (iter == _parameters.end())
    {
      std::cerr << "\nCommand line argument '" << arguments[ar]
                << "' wasn't found\n\n";
      exit(1);
    }
    if (ar+1 >= arguments.size())
    {
      std::cerr << "\nCommand line argument '" << arguments[ar] << "' doesn't "
                   "have any value\n\n";
      exit(1);
    }
    iter->second->read(arguments[ar+1]);
  }
}




std::ostream& Parameters::print_options(std::ostream &out) const
{
  std::map<std::string, std::unique_ptr<ParamBase> >::const_iterator iter = _parameters.begin();

  out << "\nAvailable options (default values in brackets)\n\n";

  for (; iter != _parameters.end(); ++iter)
  {
    const ParamBase *par = iter->second.get();
    out << add_space(iter->first, _longest_string_key_len+1)
        << par->_description
        << " [" << par->str() << "]\n";
  }

  out << "\n";

  return out;
}




std::ostream& Parameters::print_parameters(std::ostream &out) const
{
  std::map<std::string, std::unique_ptr<ParamBase> >::const_iterator iter = _parameters.begin();

  for (; iter != _parameters.end(); ++iter)
  {
    const ParamBase *par = iter->second.get();
    out << add_space(iter->first, _longest_string_key_len+1)
        << add_space(par->str(), _longest_string_value_len+1)
        << par->_description << "\n";
  }
  out << "\n";

  return out;
}




void Parameters::check_parameters() const
{
  if (_file_dofs.empty() || _file_dofs == DEFAULT_FILE_NAME)
  {
    std::cerr << "\nFile with degrees of freedom is empty or not defined\n\n";
    exit(1);
  }
  if (_file_solution_base.empty() || _file_solution_base == DEFAULT_FILE_NAME)
  {
    std::cerr << "\nBase of name of solution files is empty or not defined\n\n";
    exit(1);
  }
  if (_n_time_steps <= 0)
  {
    std::cerr << "\nNumber of time steps (" << _n_time_steps << ") must be "
                 ">0\n\n";
    exit(1);
  }
  if (_step_snapshot <= 0)
  {
    std::cerr << "\nStep of snapshot export (" << _step_snapshot << ") must be "
                 ">0\n\n";
    exit(1);
  }
}




void Parameters::update_longest_string_key_len()
{
  _longest_string_key_len = 0;

  std::map<std::string, std::unique_ptr<ParamBase> >::const_iterator iter = _parameters.begin();

  for (; iter != _parameters.end(); ++iter)
  {
    const int len_key_string = iter->first.size();

    if (len_key_string > _longest_string_key_len)
      _longest_string_key_len = len_key_string;
  }
}




void Parameters::update_longest_string_value_len()
{
  _longest_string_value_len = 0;

  std::map<std::string, std::unique_ptr<ParamBase> >::const_iterator iter = _parameters.begin();

  for (; iter != _parameters.end(); ++iter)
  {
    const ParamBase *par = iter->second.get();
    const int len_value_string = par->str().size();

    if (len_value_string > _longest_string_value_len)
      _longest_string_value_len = len_value_string;
  }
}


