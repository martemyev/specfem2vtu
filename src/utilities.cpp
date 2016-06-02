#include "utilities.hpp"

#include <cstring>


//==============================================================================
//
// Check if there is a key word in the array of command line arguments
//
//==============================================================================
int argcheck(int argc, char **argv, const char *arg)
{
  for(int i = 1; i < argc; ++i)
  {
    // strcmp returns 0 if the strings are equal
    if(strcmp(argv[i], arg) == 0)
      return(i);
  }

  return 0;
}

//==============================================================================
//
// Get a name of the given file
//
//==============================================================================
std::string file_name(const std::string &path)
{
  if (path == "") return path;

#if defined(__linux__) || defined(__APPLE__)
  // extract a filename
  const std::string fname = path.substr(path.find_last_of('/') + 1);
#elif defined(_WIN32)
  // extract a filename
  const std::string fname = path.substr(path.find_last_of('\\') + 1);
#endif

  return fname;
}

//==============================================================================
//
// Get an extension of the given file
//
//==============================================================================
std::string file_extension(const std::string &path)
{
  if (path == "") return path;

  // extract a file name from the path
  const std::string fname = file_name(path);

  // extract an extension and return it
  return fname.substr(fname.find_last_of('.'));
}

//==============================================================================
//
// Get a stem of the given file
//
//==============================================================================
std::string file_stem(const std::string &path)
{
  if (path == "") return path;

  // get a file name from the path
  const std::string fname = file_name(path);

  // extract a stem and return it
  return fname.substr(0, fname.find_last_of('.'));
}

//==============================================================================
//
// Get a path of the given file
//
//==============================================================================
std::string file_path(const std::string &path)
{
  if (path == "") return path;

#if defined(__linux__) || defined(__APPLE__)
  // extract a path
  const std::string path_ = path.substr(0, path.find_last_of('/') + 1);
#elif defined(_WIN32)
  // extract a path
  const std::string path_ = path.substr(0, path.find_last_of('\\') + 1);
#endif

  if (path_ == path)
    return ""; // no path delimeters = no path

  return path_;
}

//==============================================================================
//
// Extend the given string with empty space
//
//==============================================================================
std::string add_space(const std::string &str, int length)
{
  // how spaces need to add (if the string longer than the required length,
  // nothing is added)
  const int n_spaces = std::max(length - (int)str.size(), 0);
  return str + std::string(n_spaces, ' ');
}

//==============================================================================
//
// Fulfill the empty space of a string by zeros
//
//==============================================================================
std::string add_zeros(const std::string &str, int length)
{
  // how spaces need to add (if the string longer than the required length,
  // nothing is added)
  const int n_spaces = std::max(length - (int)str.size(), 0);
  return std::string(n_spaces, '0') + str;
}
