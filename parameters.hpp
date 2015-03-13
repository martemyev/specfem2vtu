#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <string>
#include <map>
#include <memory>
#include <vector>
#include <ostream>
#include <iostream>
#include <sstream>




//==============================================================================
//
// Base abstract class for a one parameter
//
//==============================================================================
class ParamBase
{
public:

  /// Description of the parameter (appears when help is invoked, for example)
  std::string _description;

  /// Read the value of the parameter from a string
  /// @param from a string from which the value is read
  virtual void read(const std::string &from) = 0;

  /// Convert the value of the parameter into a string
  virtual std::string str() const = 0;

  /// Comparison of the two parameters by lexicographic order of their
  /// descriptions
  static bool compare_by_desc(const std::unique_ptr<ParamBase> &a,
                              const std::unique_ptr<ParamBase> &b)
  {
    return a->_description < b->_description;
  }
};




//==============================================================================
//
// Template class for a one parameter
//
//==============================================================================
template <typename T>
class OneParam : public ParamBase
{
public:

  /// Construator takes a description of the parameter and a pointer to a value
  OneParam(const std::string &desc,
           T* val)
    : _value(val)
  {
    _description = desc;
  }

  /// Value of the parameter is saved somewhere else, and here we keep the
  /// pointer to it
  T* _value;

  /// Read the value of the parameter from a string
  /// @param from a string from which the value is read
  virtual void read(const std::string &from)
  {
    std::istringstream is(from);
    is >> *(_value);
  }

  /// Convert the value of the parameter into a string
  virtual std::string str() const
  {
    std::ostringstream os;
    os << *(_value);
    return os.str();
  }
};




//==============================================================================
//
// Class which handles the parameters of this program
//
//==============================================================================
class Parameters
{
public:

  /// Constructor takes the command line
  Parameters(int argc, char **argv);

  /// Name of a file containing description of the coordinates of the degrees of
  /// freedom
  std::string _file_dofs;

  /// Base of the name (since for each time step the filename will be different,
  /// because contains the time step) of the file with the solution
  std::string _file_solution_base;

  /// How many time steps are in the simulation
  int _n_time_steps;

  /// Typically, the solution (the snapshots) are not exported on each time step
  /// of the simulation (because it's very time consuming), therefore this
  /// variable describes how many steps are between output of the snapshots. So,
  /// basically every (_step_snapshot)-th step is outputted
  int _step_snapshot;

  /// Verbosity level
  int _verbose;

  /// The map between the key word representing a parameters, and its value
  /// (and maybe other attributes such as description)
  std::map<std::string, std::unique_ptr<ParamBase> > _parameters;

  /// Read the values from the command line
  void read_command_line(int argc, char **argv);

  /// Print possible options of the program
  std::ostream& print_options(std::ostream &out = std::cout) const;

  /// Print the parameters with which the program is going to work
  std::ostream& print_parameters(std::ostream &out = std::cout) const;

  /// Check that the parameters make sense
  void check_parameters() const;

  /// Length of the longest string representing the key words of the parameters
  int _longest_string_key_len;

  /// Length of the longest string representing the values of the parameters
  int _longest_string_value_len;

  /// Default file name
  static std::string DEFAULT_FILE_NAME;

  /// Default length of strings for printing aligned key words and values of the
  /// parameters
  static int DEFAULT_PRINT_LEN;

protected:

  /// The name of the function is self-explaining
  void update_longest_string_key_len();
  /// The name of the function is self-explaining
  void update_longest_string_value_len();
};




//==============================================================================
//
// Auxiliary functions
//
//==============================================================================
/**
 * Check if there is a string arg in the array of strings argv of length argc.
 * This is used to determine if there is an argument in a command line.
 * @return The position of the arg in the array argv, so that the value of the
 * argument can be read at the next position
 */
int argcheck(int argc, char **argv, const char *arg);
/**
 * Add some empty space to a given string str up to the given length. It's used
 * to represent all options aligned.
 * @return A string extended by spaces
 */
std::string add_space(const std::string &str, int length);
/**
 * Get (extract) a file name from the given path.
 * @param path - a name of a file under interest including the path
 * @return a string representing a name of the file.
 *         For example:
 * @verbatim
   file_name("/home/user/file.dat") = "file.dat"
 * @endverbatim
 */
std::string file_name(const std::string &path);
/**
 * Extract a stem from a filename with a path.
 * @param path - a name of a file under interest including the path
 * @return a string which represents the name of the file without an extension -
 *         only a stem of the file.
 *         For example:
 * @verbatim
   stem("/home/user/file.dat") = "file"
 * @endverbatim
 */
std::string file_stem(const std::string &path);
/**
  Get (extract) a path of the given file.
  @param path - a name of a file under interest including the path
  @return a string representing the path to the file.
          For example:
  @verbatim
  file_name("/home/user/file.dat") = "/home/user/"
  @endverbatim
 */
std::string file_path(const std::string &path);

#endif // PARAMETERS_HPP
