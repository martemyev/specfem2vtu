#include "specfem2vtu.hpp"
#include "parameters.hpp"
#include "point.hpp"
#include "utilities.hpp"

#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>

// We use the Triangle project developed by J.R. Shewchuk to build the exact
// Delaunay triangulation based on the set (cloud) of points
#define REAL double
#include "triangle.h"

// Uncomment if you want to check how the Delaunay triangulation looks like in
// Gmsh
#define WRITE_MSH

// real type for the specfem output
#define SPEC_REAL float

// PI
const double PI = 3.14159265358979323846264338327950288419716939937510582;

// Convert the radians to degrees
#define toDegrees(x) (x*180./PI)

const int TRI_DIM = 2; // triangles dimenstion (used in getting point's
                       // coordinates)
const int N_TRI_VERTICES = 3; // number of vertices of a triangle
const int N_TRI_ANGLES   = 3; // number of angles of a triangle

// These angles are used to determine if a triangle is appropriate or not (if it
// has all angles within this range [MIN_ANGLE, MAX_ANGLE], it's good)
const double MIN_ANGLE = 1.0;   // in degrees
const double MAX_ANGLE = 160.0; // in degrees


//==============================================================================
//
// Main function performing the processing of the specfem results to the vtu
// files
//
//==============================================================================
void specfem2vtu(const Parameters &param)
{
  std::vector<Point> dofs; // coordinates of the degrees of freedom

  //----------------------------------------------------------------------------
  // read the dofs
  if (param._verbose > 0)
    std::cout << "Reading dofs..." << std::endl;

  read_content(param._file_dofs, param._binary, dofs);

  if (param._verbose > 0)
    std::cout << "Reading dofs is done.\nN dofs = " << dofs.size() << std::endl;
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  // create the Delaunay triangulation based on these dofs
  if (param._verbose > 0)
    std::cout << "Delaunay triangulation of the dofs..." << std::endl;

  char switches[256];
  strcpy(switches, "");
//  strcat(switches, "D"); // force to make conforming Delaunay triangulation
  strcat(switches, "z"); // numerate everything from 0
  strcat(switches, "B"); // no boundary markers in the output files
  strcat(switches, "P"); // no output .poly file
//  strcat(switches, "N"); // no output .node file
//  strcat(switches, "E"); // no output .ele file
  strcat(switches, "Y"); // no new vertices
  strcat(switches, "Y"); // assure it again
  if (param._verbose > 1)
    strcat(switches, "V"); // verbose mode
  else
    strcat(switches, "Q"); // quiet mode
//  strcat(switches, "q15.0"); // should be no angles smaller than 15 degrees

  if (param._verbose > 1)
    std::cout << "String of parameters for Triangle: " << switches << std::endl;

  // input and output info to work with the Triangle code
  triangulateio trin, trout;

  zero_initialization(trin);
  zero_initialization(trout);

  trin.numberofpoints = dofs.size();
  trin.pointlist = (REAL*)malloc(dofs.size() * 2 * sizeof(REAL));
  for (size_t i = 0; i < dofs.size(); ++i)
  {
    trin.pointlist[2*i + 0] = dofs[i].x();
    trin.pointlist[2*i + 1] = dofs[i].y();
  }

  if (param._verbose > 1)
    std::cout << "Start triangulation process" << std::endl;
  triangulate(switches, &trin, &trout, NULL);
  if (param._verbose > 1)
    std::cout << "Triangulation process is over" << std::endl;

  if (param._verbose > 1)
    std::cout << "Choose good triangles" << std::endl;
  std::vector<std::vector<int> > good_triangles;
  choose_triangles(MIN_ANGLE, MAX_ANGLE, trout, good_triangles);
  if (param._verbose > 1)
    std::cout << "Good triangles are chosen:\n  number of discarded triangles: "
              << trout.numberoftriangles - good_triangles.size() <<  std::endl;

#if defined(WRITE_MSH)
  if (param._verbose > 0)
    std::cout << "Writing to a .msh file..." << std::endl;
  write_msh(trout, good_triangles, "test.msh");
  if (param._verbose > 0)
    std::cout << "Writing to a .msh file is done" << std::endl;
#endif

  if (param._verbose > 0)
    std::cout << "Delaunay triangulation of the dofs is done" << std::endl;
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  // read the solutions at the specific time steps, and rewrite them in the
  // vtu format
  if (param._verbose > 0)
    std::cout << "Output in vtu format..." << std::endl;
  for (int tstep = 1; tstep <= param._n_time_steps; ++tstep)
  {
    // consider only every *-th time step, where * is set up as a variable
    // param._step_snapshot
    if (tstep % param._step_snapshot != 0) continue;

    // read the solution
    std::string solfile = file_path(param._file_dofs) +
                          param._file_solution_base +
                          add_zeros(d2s(tstep), 7) + "_01_000" +
                          (param._binary ? ".bin" : ".txt");

    //std::cout << "reading..." << std::flush;
    std::vector<Point> U; // solution
    read_content(solfile, param._binary, U);
    //std::cout << "done" << std::endl;

    // write the solution in the vtu format
    //std::cout << "writing..." << std::flush;
    solfile = param._file_solution_base + d2s(tstep) + ".vtu";
    write_vtu(solfile, trout, good_triangles, U);
    //std::cout << "done" << std::endl;
  }
  if (param._verbose > 0)
    std::cout << "Output in vtu format is done" << std::endl;
  //----------------------------------------------------------------------------

  free(trin.pointlist);
  free(trout.pointlist);
  free(trout.trianglelist);
}




//==============================================================================
//
// Read the points from the given file
//
//==============================================================================
void read_content(const std::string &filename,
                  bool binary_file,
                  std::vector<Point> &content)
{
  SPEC_REAL x, y; // Cartesian coordinates of the dofs or the components of the
                  // vector field

  if (binary_file)
  {
    std::ifstream in(filename.c_str(), std::ios::binary);
    if (!in)
      throw std::runtime_error("File " + filename + " can't be opened");

    in.seekg(0, in.end);
    int length = in.tellg(); // total length of the file in bytes

    int n_elements = length / sizeof(SPEC_REAL) / 2;

    content.resize(n_elements);

    in.seekg(0, in.beg);
    for (int i = 0; i < n_elements; ++i)
    {
      in.read((char*)&x, sizeof(SPEC_REAL));
      in.read((char*)&y, sizeof(SPEC_REAL));
      content[i] = Point(x, y);
    }

    in.close();
  }
  else
  {
    std::ifstream in(filename.c_str());
    if (!in)
      throw std::runtime_error("File " + filename + " can't be opened");

    content.clear();
    while (in >> x >> y)
      content.push_back(Point(x, y));

    in.close();
  }
}




//==============================================================================
//
// Write a .msh file in the ASCII 2.2 format which can be read by Gmsh
//
//==============================================================================
void write_msh(const triangulateio &io,
               const std::vector<std::vector<int> > triangles,
               const std::string &filename)
{
  std::ofstream out(filename.c_str());
  if (!out)
    throw std::runtime_error("File " + filename + " can't be opened");

  out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
  out << "$Nodes\n";
  out << io.numberofpoints << "\n";
  for (int i = 0; i < io.numberofpoints; ++i)
  {
    out << i + 1 << " "
        << io.pointlist[i*TRI_DIM + 0] << " "
        << io.pointlist[i*TRI_DIM + 1] << " 0\n";
  }
  out << "$EndNodes\n";
  out << "$Elements\n";
  const int n_triangles = triangles.size();
  out << n_triangles << "\n";
  for (int i = 0; i < n_triangles; ++i)
  {
    out << i + 1 << " 2 2 1 1 ";
    for (int j = 0; j < N_TRI_VERTICES; ++j)
      out << triangles[i][j] + 1 << " ";
    out << "\n";
  }
  out << "$EndElements\n";

  out.close();
}




//==============================================================================
//
// Write a .vtu file in the ASCII format which can be read by ParaView
//
//==============================================================================
void write_vtu(const std::string &filename,
               const triangulateio &io,
               const std::vector<std::vector<int> > triangles,
               const std::vector<Point> &U)
{
  std::ofstream out(filename.c_str()); // open the file for writing
  if (!out)
    throw std::runtime_error("File " + filename + " cannot be opened");

  const int n_triangles = triangles.size();

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
         "byte_order=\"LittleEndian\">\n";
  out << "  <UnstructuredGrid>\n";
  out << "    <Piece NumberOfPoints=\"" << io.numberofpoints
      << "\" NumberOfCells=\"" << n_triangles << "\">\n";
  out << "      <PointData Vectors=\"U\" Scalars=\"U_magnitude\">\n";

  out << "        <DataArray type=\"Float64\" Name=\"U\" format=\"ascii\""
         " NumberOfComponents=\"" << Point::N_COORD << "\">\n";

  for (int i = 0; i < io.numberofpoints; ++i)
  {
    for (int j = 0; j < Point::N_COORD; ++j)
      out << U[i].coord(j) << " ";
  }

  out << "\n";
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"Float64\" Name=\"U_magnitude\" "
         "format=\"ascii\" NumberOfComponents=\"1\">\n";

  for (int i = 0; i < io.numberofpoints; ++i)
  {
    double magnitude = 0.0;
    for (int j = 0; j < Point::N_COORD; ++j)
      magnitude += U[i].coord(j) * U[i].coord(j);
    magnitude = sqrt(magnitude);
    out << magnitude << " ";
  }

  out << "\n";
  out << "        </DataArray>\n";

  out << "      </PointData>\n";
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\""
      << Point::N_COORD << "\" format=\"ascii\">\n";

  for (int i = 0; i < io.numberofpoints; ++i)
  {
    for (int j = 0; j < TRI_DIM; ++j)
      out << io.pointlist[i*TRI_DIM + j] << " ";
    out << "0 ";
  }

  out << "\n";
  out << "        </DataArray>\n";
  out << "      </Points>\n";
  out << "      <Cells>\n";
  out << "        <DataArray type=\"Int32\" Name=\"connectivity\" "
         "format=\"ascii\">\n";

  for (int el = 0; el < n_triangles; ++el)
  {
    // write the dofs indices
    for (int v = 0; v < N_TRI_VERTICES; ++v)
      out << triangles[el][v] << " ";
  }

  out << "\n";
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Int32\" Name=\"offsets\" "
         "format=\"ascii\">\n";

  for (int el = 0; el < n_triangles; ++el)
    out << (el + 1) * N_TRI_VERTICES << " ";
  out << "\n";

  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";

  for (int el = 0; el < n_triangles; ++el)
    out << "5 ";
  out << "\n";

  out << "        </DataArray>\n";
  out << "      </Cells>\n";
  out << "    </Piece>\n";
  out << "  </UnstructuredGrid>\n";
  out << "</VTKFile>\n";
  out.close();
}




//==============================================================================
//
// Initialize the whole triangulateio structure with zero
//
//==============================================================================
void zero_initialization(triangulateio &io)
{
  io.pointlist = NULL;
  io.pointattributelist = NULL;
  io.pointmarkerlist = 0;
  io.numberofpoints = 0;
  io.numberofpointattributes = 0;

  io.trianglelist = NULL;
  io.triangleattributelist = NULL;
  io.trianglearealist = NULL;
  io.neighborlist = NULL;
  io.numberoftriangles = 0;
  io.numberofcorners = 0;
  io.numberoftriangleattributes = 0;

  io.segmentlist = NULL;
  io.segmentmarkerlist = NULL;
  io.numberofsegments = 0;

  io.holelist = NULL;
  io.numberofholes = 0;

  io.regionlist = NULL;
  io.numberofregions = 0;

  io.edgelist = NULL;
  io.edgemarkerlist = NULL;
  io.normlist = NULL;
  io.numberofedges = 0;
}




//==============================================================================
//
// Remove the triangle with angles smaller than the 'min_angle' and larger than
// the 'max_angle'.
//
//==============================================================================
void choose_triangles(double min_angle,
                      double max_angle,
                      triangulateio &io,
                      std::vector<std::vector<int> > &good_triangles)
{
  if (io.numberofcorners != N_TRI_VERTICES)
    throw std::runtime_error("Number of vertices of a triangle is unexpected");

  // vertices of a triangle
  std::vector<int> vert_indices(N_TRI_VERTICES);
  std::vector<Point> vertices(N_TRI_VERTICES);

  // angles of a triangle
  std::vector<double> angles(N_TRI_ANGLES);

  // make sure this vector is empty
  good_triangles.clear();

  // inspect every triangle
  for (int el = 0; el < io.numberoftriangles; ++el)
  {
    // get the vertices
    for (int v = 0; v < N_TRI_VERTICES; ++v)
    {
      const int vert = io.trianglelist[el*N_TRI_VERTICES + v]; // vertex index
      vert_indices[v] = vert; // save the index
      vertices[v].coord(0) = io.pointlist[vert*TRI_DIM + 0]; // x-coordinate
      vertices[v].coord(1) = io.pointlist[vert*TRI_DIM + 1]; // y-coordinate
    }

    // get the angles
    compute_angles(vertices, angles);

    // check if the angles are okay
    bool bad_triangle = false;
    for (int a = 0; a < N_TRI_ANGLES && !bad_triangle; ++a)
    {
      if (angles[a] < min_angle || angles[a] > max_angle)
        bad_triangle = true;
    }

    // save the triangle satisfying the requirements (we save the indices of the
    // vertices only)
    if (!bad_triangle)
      good_triangles.push_back(vert_indices);
  }
}




//==============================================================================
//
// Compute the angles (in degrees) of a triangle based on its vertices.
//
//==============================================================================
void compute_angles(const std::vector<Point> &vertices,
                    std::vector<double> &angles)
{
  if (angles.empty()) angles.resize(N_TRI_ANGLES);

  if ((int)vertices.size() != N_TRI_VERTICES)
    throw std::runtime_error("Unexpected number of vertices of a triangle");

  // lengths of the triangle's sides
  const double a = length(vertices[0], vertices[1]);
  const double b = length(vertices[0], vertices[2]);
  const double c = length(vertices[1], vertices[2]);

  // first angle
  angles[0] = toDegrees(acos((b*b + c*c - a*a) / (2.*b*c)));
  // second angle
  angles[1] = toDegrees(acos((c*c + a*a - b*b) / (2.*c*a)));
  // third angle
  angles[2] = 180. - angles[0] - angles[1];
}




//==============================================================================
//
// Compute the length of a line connecting the two vertices.
//
//==============================================================================
double length(const Point &a, const Point &b)
{
  double sum = 0.;
  for (int i = 0; i < Point::N_COORD; ++i)
  {
    const double diff = a.coord(i) - b.coord(i);
    sum += diff * diff;
  }

  return sqrt(sum);
}

