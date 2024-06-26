#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/graph_traits_HalfedgeDS_default.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

template <typename PolygonMesh, typename CGAL_NP_TEMPLATE_PARAMETERS>
void triangulate_hole_w(PolygonMesh& mesh, typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge, const CGAL_NP_CLASS& np = CGAL::parameters::default_values())
{
    CGAL_precondition(CGAL::is_valid_polygon_mesh(mesh));
    CGAL_precondition(CGAL::is_triangle_mesh(mesh));
    CGAL_precondition(CGAL::is_border(border_halfedge, mesh));

    using Face_output_iterator = typename CGAL::internal_np::Lookup_named_param_def<CGAL::internal_np::face_output_iterator_t,
                                                         CGAL_NP_CLASS,
                                                         CGAL::Emptyset_iterator>::type;
    Face_output_iterator face_out = CGAL::parameters::choose_parameter<CGAL::Emptyset_iterator>(CGAL::parameters::get_parameter(np, CGAL::internal_np::face_output_iterator));
    using Vertex_output_iterator = typename CGAL::internal_np::Lookup_named_param_def<CGAL::internal_np::vertex_output_iterator_t,
                                                         CGAL_NP_CLASS,
                                                         CGAL::Emptyset_iterator>::type;
    Vertex_output_iterator vertex_out = CGAL::parameters::choose_parameter<CGAL::Emptyset_iterator>(CGAL::parameters::get_parameter(np, CGAL::internal_np::vertex_output_iterator));

    using hHalfedge = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
    using hVertex = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
    using hFacet = typename boost::graph_traits<PolygonMesh>::face_descriptor;
    using VPmap = typename CGAL::GetVertexPointMap<PolygonMesh, CGAL_NP_CLASS>::type;
    using Point_3 = typename boost::property_traits<VPmap>::value_type;
    using Kernel = CGAL::GetGeomTraits<PolygonMesh, CGAL_NP_CLASS>::type;

    VPmap vpm = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np, CGAL::internal_np::vertex_point),
                             CGAL::get_property_map(boost::vertex_point, mesh));

    auto Source     = [&mesh](hHalfedge hh)->hVertex { return CGAL::source(hh, mesh); };
    auto Target     = [&mesh](hHalfedge hh)->hVertex { return CGAL::target(hh, mesh); };
    auto Next       = [&mesh](hHalfedge hh)->hHalfedge { return CGAL::next(hh, mesh); };
    auto Prev       = [&mesh](hHalfedge hh)->hHalfedge { return CGAL::prev(hh, mesh); };
    auto Opposite   = [&mesh](hHalfedge hh)->hHalfedge { return CGAL::opposite(hh, mesh); };
    auto Facet      = [&mesh](hHalfedge hh)->hFacet { return CGAL::face(hh, mesh); };
    auto Point      = [&vpm](hVertex hv)->Point_3 { return get(vpm, hv); };

    std::list<hHalfedge> bound_edges;
    std::vector<hVertex> bound_vertices;
    for(auto edge : CGAL::halfedges_around_face(border_halfedge, mesh))
    {
        bound_edges.push_back(edge);
        bound_vertices.push_back(Target(edge));
    }

    if(bound_edges.size() < 3)  // error
    {
        return;
    }

    double avg_length = 0.0;
    for(size_t i = 0; i < bound_vertices.size(); i++)
    {
        avg_length += std::sqrt(CGAL::squared_distance(Point(bound_vertices[i]), Point(bound_vertices[(i + 1) % bound_vertices.size()])));
    }
    avg_length /= bound_vertices.size();

    while(!bound_edges.empty())
    {
        if(bound_edges.size() == 3)
        {
            CGAL::Euler::fill_hole(bound_edges.front(), mesh);
            break;
        }
        auto min_angle_edge = bound_edges.end();
        double min_angle = 999;
        for(auto it = bound_edges.begin(); it != bound_edges.end(); it++)
        {
            auto curr_edge = *it;
            auto prev_edge = Prev(curr_edge);
            auto next_edge = Next(curr_edge);
            auto p = Point(Target(curr_edge));
            auto q = Point(Source(curr_edge));
            auto r = Point(Target(next_edge));
            auto angle = CGAL::approximate_angle(q, p, r);
            auto t0 = Point(Target(Opposite(prev_edge)));
            auto t1 = Point(Target(Next(Opposite(prev_edge))));
            auto t2 = Point(Target(Prev(Opposite(prev_edge))));
            auto t3 = Point(Target(Opposite(next_edge)));
            auto t4 = Point(Target(Next(Opposite(next_edge))));
            auto t5 = Point(Target(Prev(Opposite(next_edge))));
            auto n0 = CGAL::normal(t0, t1, t2);
            auto n1 = CGAL::normal(t3, t4, t5);

            if(CGAL::collinear(p, q, r))
            {
                if(CGAL::angle(q, p, r) == CGAL::ACUTE)
                {
                    angle = 0.0;
                }
                else
                {
                    angle = 180.0;
                }
            }
            else if(CGAL::scalar_product(CGAL::normal(p, q, r), n0) > 0.0 || CGAL::scalar_product(CGAL::normal(p, q, r), n1) > 0.0)
            {
                angle += 180.0;
            }
            if(angle < min_angle)
            {
                min_angle = angle;
                min_angle_edge = it;
            }
        }
        auto itnext = (std::next(min_angle_edge) == bound_edges.end() ? bound_edges.begin() : std::next(min_angle_edge));
        auto edge = *min_angle_edge;
        auto eprev = Prev(edge);
        auto enext = Next(edge);
        auto vmid = Target(edge);
        auto vprev = Source(edge);
        auto vnext = Target(enext);
        auto pmid =  Point(vmid);
        auto pprev = Point(vprev);
        auto pnext = Point(vnext);
        double new_edge_length = CGAL::to_double(std::sqrt(CGAL::squared_distance(pprev, pnext)));
        if(new_edge_length < 2 * avg_length)
        {
            auto new_edge = CGAL::Euler::add_face_to_border(eprev, enext, mesh);
            bound_edges.insert(min_angle_edge, Opposite(new_edge));
            bound_edges.erase(min_angle_edge);
            bound_edges.erase(itnext);
            *face_out = Facet(new_edge);
            face_out++;
        }
        else
        {
            auto edge_new = CGAL::Euler::add_vertex_and_face_to_border(eprev, enext, mesh);
            put(vpm, Target(edge_new), CGAL::midpoint(pprev, pnext));
            bound_edges.insert(min_angle_edge, Next(eprev));
            bound_edges.insert(min_angle_edge, Opposite(edge_new));
            bound_edges.erase(min_angle_edge);
            bound_edges.erase(itnext);
            CGAL::Polygon_mesh_processing::triangulate_face(Facet(edge_new), mesh, CGAL::parameters::face_output_iterator(face_out));
            *vertex_out = Target(edge_new);
            vertex_out++;
        }
    }
}