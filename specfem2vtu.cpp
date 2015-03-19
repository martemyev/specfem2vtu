#include "specfem2vtu.hpp"
#include "parameters.hpp"
#include "point.hpp"
#include "utilities.hpp"

#include <fstream>
#include <cassert>

// We use the Triangle project developed by J.R. Shewchuk to build the exact
// Delaunay triangulation based on the set (cloud) of points
#define REAL double
#include "triangle.h"

// Uncomment if you want to check how the Delaunay triangulation looks like in
// Gmsh
#define WRITE_MSH

// real type for the specfem output
#define SPEC_REAL float



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
  int p = 0;
//  switches[p++] = 'D'; // force to make conforming Delaunay triangulation
  switches[p++] = 'z'; // numerate everything from 0
  switches[p++] = 'B'; // no boundary markers in the output files
  switches[p++] = 'P'; // no output .poly file
//  switches[p++] = 'N'; // no output .node file
//  switches[p++] = 'E'; // no output .ele file
  switches[p++] = 'Y'; // no new vertices
  switches[p++] = 'Y'; // assure it again
  if (param._verbose > 1)
    switches[p++] = 'V'; // verbose mode
  else
    switches[p++] = 'Q'; // quiet mode
  switches[p++] = '\0';

  // input and output info to work with the Triangle code
  triangulateio trin, trout;

  zero_initialization(trin);
  zero_initialization(trout);

//  trin.numberofpointattributes = 0;
//  trin.numberofcorners = 0;
//  trin.numberofedges = 0;
//  trin.numberofholes = 0;
//  trin.numberofregions = 0;
//  trin.numberofsegments = 0;
//  trin.numberoftriangleattributes = 0;
//  trin.numberoftriangles = 0;

  trin.numberofpoints = dofs.size();
  trin.pointlist = (REAL*)malloc(dofs.size() * 2 * sizeof(REAL));
  for (size_t i = 0; i < dofs.size(); ++i)
  {
    trin.pointlist[2*i + 0] = dofs[i].x();
    trin.pointlist[2*i + 1] = dofs[i].y();
  }

//  trout.edgelist = NULL;
//  trout.edgemarkerlist = NULL;
//  trout.holelist = NULL;
//  trout.neighborlist = NULL;
//  trout.normlist = NULL;
//  trout.pointattributelist = NULL;
//  trout.pointlist = NULL;
//  trout.pointmarkerlist = NULL;
//  trout.regionlist = NULL;
//  trout.segmentlist = NULL;
//  trout.segmentmarkerlist = NULL;
//  trout.trianglearealist = NULL;
//  trout.triangleattributelist = NULL;
//  trout.trianglelist = NULL;

  if (param._verbose > 1)
    std::cout << "Start triangulation process" << std::endl;
  triangulate(switches, &trin, &trout, NULL);
  if (param._verbose > 1)
    std::cout << "Triangulation process is over" << std::endl;

#if defined(WRITE_MSH)
  if (param._verbose > 0)
    std::cout << "Writing to a .msh file..." << std::endl;
  write_msh(trout, "test.msh");
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
    write_vtu(solfile, trout, U);
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
        << io.pointlist[2*i + 0] << " "
        << io.pointlist[2*i + 1] << " 0\n";
  }
  out << "$EndNodes\n";
  out << "$Elements\n";
  out << io.numberoftriangles << "\n";
  for (int i = 0; i < io.numberoftriangles; ++i)
  {
    out << i + 1 << " 2 2 1 1 ";
    for (int j = 0; j < io.numberofcorners; ++j)
      out << io.trianglelist[i*io.numberofcorners + j] + 1 << " ";
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
               const std::vector<Point> &U)
{
  std::ofstream out(filename.c_str()); // open the file for writing
  if (!out)
    throw std::runtime_error("File " + filename + " cannot be opened");

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
         "byte_order=\"LittleEndian\">\n";
  out << "  <UnstructuredGrid>\n";
  out << "    <Piece NumberOfPoints=\"" << io.numberofpoints
      << "\" NumberOfCells=\"" << io.numberoftriangles << "\">\n";
  out << "      <PointData Vectors=\"U\">\n";

  out << "        <DataArray type=\"Float64\" Name=\"U\" format=\"ascii\""
         " NumberOfComponents=\"" << Point::N_COORD << "\">\n";

  for (int i = 0; i < io.numberofpoints; ++i)
  {
    for (int j = 0; j < Point::N_COORD; ++j)
      out << U[i].coord(j) << " ";
    out << "\n";
  }

  out << "        </DataArray>\n";
  out << "      </PointData>\n";
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\""
      << Point::N_COORD << "\" format=\"ascii\">\n";

  const int dim = 2;
  for (int i = 0; i < io.numberofpoints; ++i)
  {
    for (int j = 0; j < dim; ++j)
      out << io.pointlist[i*dim + j] << " ";
    out << "0 \n";
  }

  out << "        </DataArray>\n";
  out << "      </Points>\n";
  out << "      <Cells>\n";
  out << "        <DataArray type=\"Int32\" Name=\"connectivity\" "
         "format=\"ascii\">\n";

  for (int el = 0; el < io.numberoftriangles; ++el)
  {
    // write the dofs indices
    for (int d = 0; d < io.numberofcorners; ++d)
      out << io.trianglelist[el*io.numberofcorners + d] << " ";
    out << "\n";
  }

  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Int32\" Name=\"offsets\" "
         "format=\"ascii\">\n";

  for (int el = 0; el < io.numberoftriangles; ++el)
    out << (el + 1) * io.numberofcorners << " ";
  out << "\n";

  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";

  for (int el = 0; el < io.numberoftriangles; ++el)
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
