#include <iostream>
#include <MeshHoleFilling.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/Polygon_mesh_processing/border.h>

using Kernel = CGAL::Simple_cartesian<double>;
using Polyhedron = CGAL::Polyhedron_3<Kernel>;

int main(int argc, char* argv[])
{
    Polyhedron mesh;
    CGAL::IO::read_polygon_mesh("C:\\Dev\\MeshHoleFilling\\test\\L.ply", mesh);
    std::vector<Polyhedron::Halfedge_handle> border_edges;
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(mesh, std::back_inserter(border_edges));
    for(auto hh : border_edges)
    {
        triangulate_hole_w(mesh, hh);
    }
    CGAL::IO::write_polygon_mesh("C:\\Dev\\MeshHoleFilling\\test\\Lout.obj", mesh);
    return 0;
}