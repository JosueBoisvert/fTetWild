//
// Created by jboisvert on 2022-09-28.
//

#ifndef FLOATTETWILD_FTETWILD_HPP
#define FLOATTETWILD_FTETWILD_HPP

#include <Eigen/Dense>
#include <array>
#include <cassert>
#include <queue>
#include <random>
#include <sstream>
#include <unordered_set>
#include <vector>

#include <geogram/mesh/mesh_geometry.h>
#include "external/mesh_AABB.h"

namespace floatTetWild {

#define SCALAR_ZERO 1e-8
#define SCALAR_ZERO_2 1e-16
#define SCALAR_ZERO_3 1e-24
#define NOT_SURFACE SCHAR_MAX  // SHRT_MAX//INT_MAX
#define KNOWN_NOT_SURFACE (-SCHAR_MAX / 2)
#define KNOWN_SURFACE (SCHAR_MAX / 2)
#define NO_SURFACE_TAG 0
#define NOT_BBOX (-1)
#define OPP_T_ID_UNKNOWN (-2)
#define OPP_T_ID_BOUNDARY (-1)

class Parameters
{
  public:
    std::string log_path;
    std::string input_path;
    std::string output_path;
    std::string tag_path;
    std::string postfix;
    std::string envelope_log;
    std::string envelope_log_csv;

    bool not_sort_input              = false;
    bool correct_surface_orientation = false;
    bool is_quiet                    = false;
    int  log_level                   = 3;  // 2;

    bool smooth_open_boundary = false;
    bool manifold_surface     = false;
    bool disable_filtering    = false;
    bool use_floodfill        = false;
    bool use_general_wn       = false;
    bool use_input_for_wn     = false;
    bool coarsen              = false;

    bool            apply_sizing_field = false;
    Eigen::VectorXd V_sizing_field;
    Eigen::VectorXi T_sizing_field;
    Eigen::VectorXd values_sizing_field;
    // get sizing field value for a point
    std::function<double(const Eigen::Matrix<double, 3, 1>&)> get_sizing_field_value;

#ifdef NEW_ENVELOPE
    std::vector<double> input_epsr_tags;  // same length as the list of input faces
#endif

    // it decides the scale of the box, presents the deviation of the box from the model
    //( in % of  max((xmax-xmin), (ymax-ymin), (zmax-zmin)) of the input points)
    double box_scale = 1 / 15.0;

    // epsilon presents the tolerence permited (in % of the box diagonal)
    double eps_rel = 1e-3;

    // initial target edge length at every vertex(in % of the box diagonal)
    double ideal_edge_length_rel = 1 / 20.0;
    double min_edge_len_rel      = -1;

    int    max_its     = 80;
    double stop_energy = 10;

#ifdef NEW_ENVELOPE
    int stage = 1;
#else
    int stage = 2;
#endif

    unsigned int num_threads = std::numeric_limits<unsigned int>::max();

    int stop_p = -1;

    Eigen::Matrix<double, 3, 1> bbox_min;
    Eigen::Matrix<double, 3, 1> bbox_max;
    double                      bbox_diag_length;
    double                      ideal_edge_length;
    double                      ideal_edge_length_2;
    double                      eps_input;
    double                      eps;
    double                      eps_delta;
    double                      eps_2;
    double                      dd;
    double                      min_edge_length;

    double split_threshold;
    double collapse_threshold;
    double split_threshold_2;
    double collapse_threshold_2;

    double eps_coplanar;
    double eps_2_coplanar;
    double eps_simplification;
    double eps_2_simplification;
    double dd_simplification;

    bool init(double bbox_diag_l);
};

class MeshVertex
{
  public:
    std::vector<int>            conn_tets;
    Eigen::Matrix<double, 3, 1> pos;

    bool is_on_surface    = false;
    bool is_on_boundary   = false;
    bool is_on_cut        = false;
    int  on_boundary_e_id = -1;
    bool is_on_bbox       = false;
    bool is_outside       = false;

    bool is_removed = false;
    bool is_freezed = false;  // todo

    double sizing_scalar = 1;

    double scalar = 0;

  public:
    MeshVertex() = default;
    explicit MeshVertex(const Eigen::Matrix<double, 3, 1>& p)
        : pos(p) {};

    inline double& operator[](const int index)
    {
        assert(index >= 0 && index < 3);
        return pos[index];
    }

    inline double operator[](const int index) const
    {
        assert(index >= 0 && index < 3);
        return pos[index];
    }
};

class MeshTet
{
  public:
    Eigen::Matrix<int, 4, 1> indices;
    std::array<char, 4>      is_surface_fs = {{NOT_SURFACE, NOT_SURFACE, NOT_SURFACE, NOT_SURFACE}};
    std::array<char, 4>      is_bbox_fs    = {{NOT_BBOX, NOT_BBOX, NOT_BBOX, NOT_BBOX}};
    std::array<int, 4>       opp_t_ids     = {
      {OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN}};
    std::array<char, 4> surface_tags = {{0, 0, 0, 0}};

    double quality    = 0;
    double scalar     = 0;
    bool   is_removed = false;
    bool   is_outside = false;

  public:
    MeshTet() = default;
    explicit MeshTet(const Eigen::Matrix<int, 4, 1>& idx)
        : indices(idx) {};
    MeshTet(int v0, int v1, int v2, int v3)
        : indices(v0, v1, v2, v3) {};

    inline void reset()
    {
        is_surface_fs = {{NOT_SURFACE, NOT_SURFACE, NOT_SURFACE, NOT_SURFACE}};
        is_bbox_fs    = {{NOT_BBOX, NOT_BBOX, NOT_BBOX, NOT_BBOX}};
        opp_t_ids     = {{OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN}};
        surface_tags  = {{0, 0, 0, 0}};

        quality    = 0;
        scalar     = 0;
        is_removed = false;
        is_outside = false;
    }

    inline const Eigen::Matrix<int, 4, 1>& pts_indices() const { return indices; }

    inline int& operator[](const int index)
    {
        assert(index >= 0 && index < 4);
        return indices[index];
    }

    inline int operator[](const int index) const
    {
        assert(index >= 0 && index < 4);
        return indices[index];
    }

    inline int find(int ele) const
    {
        for (int j = 0; j < 4; j++)
            if (indices[j] == ele)
                return j;
        return -1;
    }

    inline int find_opp(int v0_id, int v1_id, int v2_id) const
    {
        for (int j = 0; j < 4; j++)
            if (indices[j] != v0_id && indices[j] != v1_id && indices[j] != v2_id)
                return j;
        return -1;
    }

    inline std::string to_string() const
    {
        std::ostringstream ss;
        ss << indices[0] << " " << indices[1] << " " << indices[2] << " " << indices[3];
        return ss.str();
    }
};

class Mesh
{
  public:
    std::vector<MeshVertex> tet_vertices;
    std::vector<MeshTet>    tets;
    Parameters              params;

    int  t_empty_start   = 0;
    int  v_empty_start   = 0;
    bool is_limit_length = true;
    bool is_closed       = true;

    bool is_input_all_inserted = false;
    bool is_coarsening         = false;

  public:
    void one_ring_vertex_coloring(std::vector<double>& colors) const;

    void one_ring_vertex_sets(const int                      threshold,
                              std::vector<std::vector<int>>& concurrent_sets,
                              std::vector<int>&              serial_set) const;

    void partition(const int n_parts, std::vector<std::vector<int>>& tets_id) const;

    static void one_ring_edge_set(const std::vector<std::array<int, 2>>&          edges,
                                  const std::vector<bool>&                        v_is_removed,
                                  const std::vector<bool>&                        f_is_removed,
                                  const std::vector<std::unordered_set<int>>&     conn_fs,
                                  const std::vector<Eigen::Matrix<double, 3, 1>>& input_vertices,
                                  std::vector<int>&                               safe_set);

    inline int t_empty_size() const
    {
        int cnt = 0;
        for (const auto& t : tets) {
            if (t.is_removed)
                cnt++;
        }
        return cnt;
    }

    inline int v_empty_size() const
    {
        int cnt = 0;
        for (const auto& v : tet_vertices) {
            if (v.is_removed)
                cnt++;
        }
        return cnt;
    }

    inline void reset_t_empty_start()
    {
        t_empty_start = tets.size();
        for (int i = 0; i < tets.size(); i++) {
            if (tets[i].is_removed) {
                t_empty_start = i;
                break;
            }
        }
    }

    inline void reset_v_empty_start()
    {
        v_empty_start = tet_vertices.size();
        for (int i = 0; i < tet_vertices.size(); i++) {
            if (tet_vertices[i].is_removed) {
                v_empty_start = i;
                break;
            }
        }
    }

    inline int get_v_num() const
    {
        int cnt = 0;
        for (int i = 0; i < tet_vertices.size(); i++) {
            if (!tet_vertices[i].is_removed)
                cnt++;
        }
        return cnt;
    }

    inline int get_t_num() const
    {
        int cnt = 0;
        for (int i = 0; i < tets.size(); i++) {
            if (!tets[i].is_removed)
                cnt++;
        }
        return cnt;
    }

    inline double get_max_energy() const
    {
        double max_energy = 0;
        int    cnt        = 0;
        for (auto& t : tets) {
            if (t.is_removed)
                continue;
            if (t.quality > max_energy)
                max_energy = t.quality;
        }
        return max_energy;
    }

    inline double get_avg_energy() const
    {
        double avg_energy = 0;
        int    cnt        = 0;
        for (auto& t : tets) {
            if (t.is_removed)
                continue;
            avg_energy += t.quality;
            cnt++;
        }
        avg_energy /= cnt;
        return avg_energy;
    }
};

class AABBWrapper
{
  public:
    GEO::Mesh        b_mesh;
    GEO::Mesh        tmp_b_mesh;
    const GEO::Mesh& sf_mesh;

    std::shared_ptr<MeshFacetsAABBWithEps> b_tree;
    std::shared_ptr<MeshFacetsAABBWithEps> tmp_b_tree;
    MeshFacetsAABBWithEps                  sf_tree;

  public:
    inline double get_sf_diag() const { return GEO::bbox_diagonal(sf_mesh); }
    explicit AABBWrapper(const GEO::Mesh& sf_mesh)
        : sf_mesh(sf_mesh)
        , sf_tree(sf_mesh) {};

#ifdef NEW_ENVELOPE
    fastEnvelope::FastEnvelope b_tree_exact;
    fastEnvelope::FastEnvelope tmp_b_tree_exact;
    fastEnvelope::FastEnvelope sf_tree_exact;
    fastEnvelope::FastEnvelope sf_tree_exact_simplify;

    inline void init_sf_tree(const std::vector<Eigen::Matrix<double, 3, 1>>& vs,
                             const std::vector<Eigen::Matrix<int, 3, 1>>&    fs,
                             double                                          eps)
    {
        sf_tree_exact.init(vs, fs, eps);
        sf_tree_exact_simplify.init(vs, fs, 0.8 * eps);
    }
    inline void init_sf_tree(const std::vector<Eigen::Matrix<double, 3, 1>>& vs,
                             const std::vector<Eigen::Matrix<int, 3, 1>>&    fs,
                             std::vector<double>&                            eps,
                             double                                          bbox_diag_length)
    {
        for (auto& e : eps)
            e *= bbox_diag_length;
        sf_tree_exact.init(vs, fs, eps);
        std::vector<double> eps_simplify = eps;
        for (auto& e : eps_simplify)
            e *= 0.8;
        sf_tree_exact_simplify.init(vs, fs, eps_simplify);
    }
#endif

    void init_b_mesh_and_tree(const std::vector<Eigen::Matrix<double, 3, 1>>& input_vertices,
                              const std::vector<Eigen::Matrix<int, 3, 1>>&    input_faces,
                              Mesh&                                           mesh);

    void init_tmp_b_mesh_and_tree(const std::vector<Eigen::Matrix<double, 3, 1>>& input_vertices,
                                  const std::vector<Eigen::Matrix<int, 3, 1>>&    input_faces,
                                  const std::vector<std::array<int, 2>>&          b_edges1,
                                  const Mesh&                                     mesh,
                                  const std::vector<std::array<int, 2>>&          b_edges2);

    double project_to_sf(Eigen::Matrix<double, 3, 1>& p) const;

    double project_to_b(Eigen::Matrix<double, 3, 1>& p) const;

    double project_to_tmp_b(Eigen::Matrix<double, 3, 1>& p) const;

    int get_nearest_face_sf(const Eigen::Matrix<double, 3, 1>& p) const;

    double get_sq_dist_to_sf(const Eigen::Matrix<double, 3, 1>& p) const;

    // envelope check - triangle
    bool is_out_sf_envelope(const std::vector<GEO::vec3>& ps,
                            const double                  eps_2,
                            GEO::index_t                  prev_facet = GEO::NO_FACET) const;

    bool is_out_b_envelope(const std::vector<GEO::vec3>& ps,
                           const double                  eps_2,
                           GEO::index_t                  prev_facet = GEO::NO_FACET) const;

    bool is_out_tmp_b_envelope(const std::vector<GEO::vec3>& ps,
                               const double                  eps_2,
                               GEO::index_t                  prev_facet = GEO::NO_FACET) const;

#ifdef NEW_ENVELOPE
    inline bool is_out_sf_envelope_exact(
      const std::array<Eigen::Matrix<double, 3, 1>, 3>& triangle) const
    {
        return sf_tree_exact.is_outside(triangle);
    }

    inline bool is_out_sf_envelope_exact_simplify(
      const std::array<Eigen::Matrix<double, 3, 1>, 3>& triangle) const
    {
        return sf_tree_exact_simplify.is_outside(triangle);
    }

    inline bool is_out_b_envelope_exact(
      const std::array<Eigen::Matrix<double, 3, 1>, 3>& triangle) const
    {
        return b_tree_exact.is_outside(triangle);
    }

    inline bool is_out_tmp_b_envelope_exact(
      const std::array<Eigen::Matrix<double, 3, 1>, 3>& triangle) const
    {
        return tmp_b_tree_exact.is_outside(triangle);
    }
#endif

    // envelope check - point
    bool is_out_sf_envelope(const Eigen::Matrix<double, 3, 1>& p,
                            const double                       eps_2,
                            GEO::index_t&                      prev_facet) const;

    bool is_out_sf_envelope(const Eigen::Matrix<double, 3, 1>& p,
                            const double                       eps_2,
                            GEO::index_t&                      prev_facet,
                            double&                            sq_dist,
                            GEO::vec3&                         nearest_p) const;
    bool is_out_sf_envelope(const GEO::vec3& geo_p,
                            const double     eps_2,
                            GEO::index_t&    prev_facet,
                            double&          sq_dist,
                            GEO::vec3&       nearest_p) const;

    bool is_out_b_envelope(const Eigen::Matrix<double, 3, 1>& p,
                           const double                       eps_2,
                           GEO::index_t&                      prev_facet) const;

    bool is_out_tmp_b_envelope(const Eigen::Matrix<double, 3, 1>& p,
                               const double                       eps_2,
                               GEO::index_t&                      prev_facet) const;

#ifdef NEW_ENVELOPE
    inline bool is_out_sf_envelope_exact(const Eigen::Matrix<double, 3, 1>& p) const
    {
        return sf_tree_exact.is_outside(p);
    }

    inline bool is_out_b_envelope_exact(const Eigen::Matrix<double, 3, 1>& p) const
    {
        return b_tree_exact.is_outside(p);
    }

    inline bool is_out_tmp_b_envelope_exact(const Eigen::Matrix<double, 3, 1>& p) const
    {
        return tmp_b_tree_exact.is_outside(p);
    }
#endif

    double dist_sf_envelope(const std::vector<GEO::vec3>& ps,
                            const double                  eps_2,
                            GEO::index_t                  prev_facet = GEO::NO_FACET) const;
};
}  // namespace floatTetWild

#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/Statistics.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/CSGTreeParser.hpp>
#include <floattetwild/Predicates.hpp>
#include <floattetwild/Types.hpp>

#endif  // FLOATTETWILD_FTETWILD_HPP
