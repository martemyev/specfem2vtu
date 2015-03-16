#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <stdexcept>
#include <sstream>

//==============================================================================
//
// d2s<T> - convert data of type T to a string
//
//==============================================================================
/**
 * Convert the data of any type which has oveloaded operator '<<' to string
 * @param data - the data
 * @param scientific - use scientific (exponential) format or not
 * @param precision - if scientific format is used, we can change the precision
 * @return data in string format
 */
template <typename T>
inline std::string d2s(T data,
                       bool scientific = false,
                       int precision = 6,
                       bool noperiod = false)
{
  std::ostringstream o;
  if (scientific) {
    o.setf(std::ios::scientific);
    o.precision(precision);
  }
  if (!(o << data))
    throw std::runtime_error("Bad conversion of data to string!");
  if (noperiod) { // eliminate a period in case of 'double' numbers
    std::string res = o.str(); // resulting string
    std::string::size_type pos = res.find('.');
    if (pos != std::string::npos)
      res.erase(pos, 1);
    return res;
  }
  return o.str();
}

//==============================================================================
//
// A set of functions on strings (file names, etc)
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
 * Add some empty space to the end of the given string 'str' up to the given
 * length. It's used to represent all options aligned.
 * @return A string extended by spaces
 */
std::string add_space(const std::string &str, int length);
/**
 * Add some zeros to the beginning of the given string 'str' up to the given
 * length.
 * @return A string extended by zeros
 */
std::string add_zeros(const std::string &str, int length);
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
 * Get a file extension.
 * @param path - a name of a file under interest including the path
 * @return a string representing an extension of the file with a dot before it.
 *         For example:
 * @verbatim
   file_extension("/home/user/file.dat") = ".dat"
 * @endverbatim
 */
std::string file_extension(const std::string &path);
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



#endif // UTILITIES_HPP
