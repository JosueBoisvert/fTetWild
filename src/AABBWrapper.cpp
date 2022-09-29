#include <floattetwild/LocalOperations.h>
#include <floattetwild/TriangleInsertion.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/mesh/mesh_reorder.h>
#include "ftetwild.hpp"

namespace floatTetWild {
void AABBWrapper::init_b_mesh_and_tree(
  const std::vector<Eigen::Matrix<double, 3, 1>>& input_vertices,
  const std::vector<Eigen::Matrix<int, 3, 1>>&    input_faces,
  Mesh&                                           mesh)
{
    b_mesh.clear(false, false);
    std::vector<std::vector<int>>   conn_tris(input_vertices.size());
    std::vector<std::array<int, 2>> all_edges;
    all_edges.reserve(input_faces.size() * 3);
    for (int i = 0; i < input_faces.size(); i++) {
        for (int j = 0; j < 3; j++) {
            conn_tris[input_faces[i][j]].push_back(i);
            if (input_faces[i][j] < input_faces[i][(j + 1) % 3])
                all_edges.push_back({{input_faces[i][j], input_faces[i][(j + 1) % 3]}});
            else
                all_edges.push_back({{input_faces[i][(j + 1) % 3], input_faces[i][j]}});
        }
    }
    vector_unique(all_edges);

    std::vector<std::pair<std::array<int, 2>, std::vector<int>>> _;
    std::vector<std::array<int, 2>>                              b_edges;
    std::vector<bool>                                            _1;
    find_boundary_edges(input_vertices,
                        input_faces,
                        std::vector<bool>(input_faces.size(), true),
                        std::vector<bool>(input_faces.size(), true),
                        _,
                        _1,
                        b_edges);

    if (b_edges.empty()) {
        b_mesh.vertices.clear();
        b_mesh.vertices.create_vertices(1);
        b_mesh.vertices.point(0) = GEO::vec3(0, 0, 0);
        b_mesh.facets.clear();
        b_mesh.facets.create_triangles(1);
        b_mesh.facets.set_vertex(0, 0, 0);
        b_mesh.facets.set_vertex(0, 1, 0);
        b_mesh.facets.set_vertex(0, 2, 0);
    }
    else {
        b_mesh.vertices.clear();
        b_mesh.vertices.create_vertices((int)b_edges.size() * 2);
        int cnt = 0;
        for (auto& e : b_edges) {
            for (int j = 0; j < 2; j++) {
                GEO::vec3& p = b_mesh.vertices.point(cnt++);
                p[0]         = input_vertices[e[j]][0];
                p[1]         = input_vertices[e[j]][1];
                p[2]         = input_vertices[e[j]][2];
            }
        }
        b_mesh.facets.clear();
        b_mesh.facets.create_triangles((int)b_edges.size());
        for (int i = 0; i < b_edges.size(); i++) {
            b_mesh.facets.set_vertex(i, 0, i * 2);
            b_mesh.facets.set_vertex(i, 1, i * 2);
            b_mesh.facets.set_vertex(i, 2, i * 2 + 1);
        }
    }

    mesh_reorder(b_mesh, GEO::MESH_ORDER_MORTON);
    b_tree = std::make_shared<MeshFacetsAABBWithEps>(b_mesh);

    if (b_edges.empty())
        mesh.is_closed = true;

#ifdef NEW_ENVELOPE
    std::vector<Eigen::Matrix<double, 3, 1>> vs;
    std::vector<Eigen::Matrix<int, 3, 1>>    fs;
    if (b_edges.empty()) {
        vs.push_back(Eigen::Matrix<double, 3, 1>(0, 0, 0));
        fs.push_back(Eigen::Matrix<int, 3, 1>(0, 0, 0));
    }
    else {
        vs.resize(b_edges.size() * 2);
        fs.resize(b_edges.size());
        for (int i = 0; i < b_edges.size(); i++) {
            vs[i * 2]     = input_vertices[b_edges[i][0]];
            vs[i * 2 + 1] = input_vertices[b_edges[i][1]];
            fs[i]         = Eigen::Matrix<int, 3, 1>(i * 2, i * 2 + 1, i * 2 + 1);
        }
    }
    //    b_tree_exact = std::make_shared<fastEnvelope::FastEnvelope>(vs, fs, eps);
    b_tree_exact.init(vs, fs, mesh.params.eps);
#endif
}

void AABBWrapper::init_tmp_b_mesh_and_tree(
  const std::vector<Eigen::Matrix<double, 3, 1>>& input_vertices,
  const std::vector<Eigen::Matrix<int, 3, 1>>&    input_faces,
  const std::vector<std::array<int, 2>>&          b_edges1,
  const Mesh&                                     mesh,
  const std::vector<std::array<int, 2>>&          b_edges2)
{
    if (b_edges1.empty() && b_edges2.empty()) {
        tmp_b_mesh.vertices.clear();
        tmp_b_mesh.vertices.create_vertices(1);
        tmp_b_mesh.vertices.point(0) = GEO::vec3(0, 0, 0);
        tmp_b_mesh.facets.clear();
        tmp_b_mesh.facets.create_triangles(1);
        tmp_b_mesh.facets.set_vertex(0, 0, 0);
        tmp_b_mesh.facets.set_vertex(0, 1, 0);
        tmp_b_mesh.facets.set_vertex(0, 2, 0);
    }
    else {
        tmp_b_mesh.vertices.clear();
        tmp_b_mesh.vertices.create_vertices((int)(b_edges1.size() + b_edges2.size()) * 2);
        int cnt = 0;
        for (auto& e : b_edges1) {
            for (int j = 0; j < 2; j++) {
                GEO::vec3& p = tmp_b_mesh.vertices.point(cnt++);
                p[0]         = input_vertices[e[j]][0];
                p[1]         = input_vertices[e[j]][1];
                p[2]         = input_vertices[e[j]][2];
            }
        }
        for (auto& e : b_edges2) {
            for (int j = 0; j < 2; j++) {
                GEO::vec3& p = tmp_b_mesh.vertices.point(cnt++);
                p[0]         = mesh.tet_vertices[e[j]].pos[0];
                p[1]         = mesh.tet_vertices[e[j]].pos[1];
                p[2]         = mesh.tet_vertices[e[j]].pos[2];
            }
        }

        tmp_b_mesh.facets.clear();
        tmp_b_mesh.facets.create_triangles((int)b_edges1.size() + b_edges2.size());
        for (int i = 0; i < b_edges1.size(); i++) {
            tmp_b_mesh.facets.set_vertex(i, 0, i * 2);
            tmp_b_mesh.facets.set_vertex(i, 1, i * 2);
            tmp_b_mesh.facets.set_vertex(i, 2, i * 2 + 1);
        }
        for (int i = b_edges1.size(); i < b_edges1.size() + b_edges2.size(); i++) {
            tmp_b_mesh.facets.set_vertex(i, 0, i * 2);
            tmp_b_mesh.facets.set_vertex(i, 1, i * 2);
            tmp_b_mesh.facets.set_vertex(i, 2, i * 2 + 1);
        }
    }
    mesh_reorder(tmp_b_mesh, GEO::MESH_ORDER_MORTON);
    tmp_b_tree = std::make_shared<MeshFacetsAABBWithEps>(tmp_b_mesh);

#ifdef NEW_ENVELOPE
    std::vector<Eigen::Matrix<double, 3, 1>> vs;
    std::vector<Eigen::Matrix<int, 3, 1>>    fs;
    if (b_edges1.empty() && b_edges2.empty()) {
        vs.push_back(Eigen::Matrix<double, 3, 1>(0, 0, 0));
        fs.push_back(Eigen::Matrix<int, 3, 1>(0, 0, 0));
    }
    else {
        vs.resize((b_edges1.size() + b_edges2.size()) * 2);
        fs.resize(b_edges1.size() + b_edges2.size());
        for (int i = 0; i < b_edges1.size(); i++) {
            vs[i * 2]     = input_vertices[b_edges1[i][0]];
            vs[i * 2 + 1] = input_vertices[b_edges1[i][1]];
            fs[i]         = Eigen::Matrix<int, 3, 1>(i * 2, i * 2 + 1, i * 2 + 1);
        }
        for (int i = b_edges1.size(); i < b_edges1.size() + b_edges2.size(); i++) {
            vs[i * 2]     = mesh.tet_vertices[b_edges2[i - b_edges1.size()][0]].pos;
            vs[i * 2 + 1] = mesh.tet_vertices[b_edges2[i - b_edges1.size()][1]].pos;
            fs[i]         = Eigen::Matrix<int, 3, 1>(i * 2, i * 2 + 1, i * 2 + 1);
        }
    }
    //    tmp_b_tree_exact = std::make_shared<fastEnvelope::FastEnvelope>(vs, fs,
    //    mesh.params.eps_input);
    tmp_b_tree_exact.init(vs, fs, mesh.params.eps);
#endif
}
int AABBWrapper::get_nearest_face_sf(const Eigen::Matrix<double, 3, 1>& p) const
{
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    double    sq_dist = std::numeric_limits<double>::max();
    return sf_tree.nearest_facet(geo_p, nearest_p, sq_dist);
}
double AABBWrapper::get_sq_dist_to_sf(const Eigen::Matrix<double, 3, 1>& p) const
{
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    double    sq_dist = std::numeric_limits<double>::max();
    sf_tree.nearest_facet(geo_p, nearest_p, sq_dist);
    return sq_dist;
}
bool AABBWrapper::is_out_sf_envelope(const std::vector<GEO::vec3>& ps,
                                     const double                  eps_2,
                                     GEO::index_t                  prev_facet) const
{
    GEO::vec3 nearest_point;
    double    sq_dist = std::numeric_limits<double>::max();

    for (const GEO::vec3& current_point : ps) {
        if (prev_facet != GEO::NO_FACET) {
            get_point_facet_nearest_point(
              sf_mesh, current_point, prev_facet, nearest_point, sq_dist);
        }
        if (double(sq_dist) > eps_2) {
            sf_tree.facet_in_envelope_with_hint(
              current_point, eps_2, prev_facet, nearest_point, sq_dist);
        }
        if (double(sq_dist) > eps_2) {
            return true;
        }
    }

    return false;
}
bool AABBWrapper::is_out_b_envelope(const std::vector<GEO::vec3>& ps,
                                    const double                  eps_2,
                                    GEO::index_t                  prev_facet) const
{
    GEO::vec3 nearest_point;
    double    sq_dist = std::numeric_limits<double>::max();

    for (const GEO::vec3& current_point : ps) {
        if (prev_facet != GEO::NO_FACET) {
            get_point_facet_nearest_point(
              b_mesh, current_point, prev_facet, nearest_point, sq_dist);
        }
        if (double(sq_dist) > eps_2) {
            b_tree->facet_in_envelope_with_hint(
              current_point, eps_2, prev_facet, nearest_point, sq_dist);
        }
        if (double(sq_dist) > eps_2) {
            return true;
        }
    }

    return false;
}
bool AABBWrapper::is_out_tmp_b_envelope(const std::vector<GEO::vec3>& ps,
                                        const double                  eps_2,
                                        GEO::index_t                  prev_facet) const
{
    GEO::vec3 nearest_point;
    double    sq_dist = std::numeric_limits<double>::max();

    for (const GEO::vec3& current_point : ps) {
        if (prev_facet != GEO::NO_FACET) {
            get_point_facet_nearest_point(
              tmp_b_mesh, current_point, prev_facet, nearest_point, sq_dist);
        }
        if (double(sq_dist) > eps_2) {
            tmp_b_tree->facet_in_envelope_with_hint(
              current_point, eps_2, prev_facet, nearest_point, sq_dist);
        }
        if (double(sq_dist) > eps_2) {
            return true;
        }
    }

    return false;
}
bool AABBWrapper::is_out_sf_envelope(const Eigen::Matrix<double, 3, 1>& p,
                                     const double                       eps_2,
                                     GEO::index_t&                      prev_facet) const
{
    GEO::vec3 nearest_p;
    double    sq_dist;
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    prev_facet = sf_tree.facet_in_envelope(geo_p, eps_2, nearest_p, sq_dist);

    if (double(sq_dist) > eps_2)
        return true;
    return false;
}
bool AABBWrapper::is_out_sf_envelope(const GEO::vec3& geo_p,
                                     const double     eps_2,
                                     GEO::index_t&    prev_facet,
                                     double&          sq_dist,
                                     GEO::vec3&       nearest_p) const
{
    if (prev_facet != GEO::NO_FACET) {
        get_point_facet_nearest_point(sf_mesh, geo_p, prev_facet, nearest_p, sq_dist);
    }
    if (double(sq_dist) > eps_2) {
        sf_tree.facet_in_envelope_with_hint(geo_p, eps_2, prev_facet, nearest_p, sq_dist);
    }

    if (double(sq_dist) > eps_2)
        return true;
    return false;
}
bool AABBWrapper::is_out_sf_envelope(const Eigen::Matrix<double, 3, 1>& p,
                                     const double                       eps_2,
                                     GEO::index_t&                      prev_facet,
                                     double&                            sq_dist,
                                     GEO::vec3&                         nearest_p) const
{
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    return is_out_sf_envelope(geo_p, eps_2, prev_facet, sq_dist, nearest_p);
}
bool AABBWrapper::is_out_b_envelope(const Eigen::Matrix<double, 3, 1>& p,
                                    const double                       eps_2,
                                    GEO::index_t&                      prev_facet) const
{
    GEO::vec3 nearest_p;
    double    sq_dist;
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    prev_facet = b_tree->facet_in_envelope(geo_p, eps_2, nearest_p, sq_dist);

    if (double(sq_dist) > eps_2)
        return true;
    return false;
}
bool AABBWrapper::is_out_tmp_b_envelope(const Eigen::Matrix<double, 3, 1>& p,
                                        const double                       eps_2,
                                        GEO::index_t&                      prev_facet) const
{
    GEO::vec3 nearest_p;
    double    sq_dist;
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    prev_facet = tmp_b_tree->facet_in_envelope(geo_p, eps_2, nearest_p, sq_dist);

    if (double(sq_dist) > eps_2)
        return true;
    return false;
}
double AABBWrapper::dist_sf_envelope(const std::vector<GEO::vec3>& ps,
                                     const double                  eps_2,
                                     GEO::index_t                  prev_facet) const
{  /// only used for checking correctness
    GEO::vec3 nearest_point;
    double    sq_dist = std::numeric_limits<double>::max();

    for (const GEO::vec3& current_point : ps) {
        if (prev_facet != GEO::NO_FACET) {
            get_point_facet_nearest_point(
              sf_mesh, current_point, prev_facet, nearest_point, sq_dist);
        }
        if (double(sq_dist) > eps_2) {
            sf_tree.facet_in_envelope_with_hint(
              current_point, eps_2, prev_facet, nearest_point, sq_dist);
        }
        if (double(sq_dist) > eps_2) {
            return sq_dist;
        }
    }

    return 0;
}
double AABBWrapper::project_to_sf(Eigen::Matrix<double, 3, 1> &p) const {
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max(); //??
    sf_tree.nearest_facet(geo_p, nearest_p, sq_dist);
    p[0] = nearest_p[0];
    p[1] = nearest_p[1];
    p[2] = nearest_p[2];

    return sq_dist;
}
double AABBWrapper::project_to_b(Eigen::Matrix<double, 3, 1> &p) const {
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max(); //?
    b_tree->nearest_facet(geo_p, nearest_p, sq_dist);
    p[0] = nearest_p[0];
    p[1] = nearest_p[1];
    p[2] = nearest_p[2];

    return sq_dist;
}
double AABBWrapper::project_to_tmp_b(Eigen::Matrix<double, 3, 1> &p) const {
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max(); //?
    tmp_b_tree->nearest_facet(geo_p, nearest_p, sq_dist);
    p[0] = nearest_p[0];
    p[1] = nearest_p[1];
    p[2] = nearest_p[2];

    return sq_dist;
}
};  // namespace floatTetWild
