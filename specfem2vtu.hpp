#ifndef SPECFEM2VTU_HPP
#define SPECFEM2VTU_HPP

#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>


class Parameters;
class Point;
struct triangulateio;



void specfem2vtu(const Parameters &param);

void read_content(const std::string &filename,
                  bool binary_file,
                  std::vector<Point> &content);

void write_msh(const triangulateio &io,
               const std::string &filename);

void write_vtu(const std::string &filename,
               const triangulateio &io,
               const std::vector<Point> &U);

void zero_initialization(triangulateio &io);

#endif // SPECFEM2VTU_HPP
