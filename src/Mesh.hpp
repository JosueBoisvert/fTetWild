// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once
#ifndef FLOATTETWILD_MESH_HPP
#define FLOATTETWILD_MESH_HPP

//#include <floattetwild/Parameters.h>
//#include <floattetwild/Types.hpp>
//
//#include <vector>
//#include <array>
//#include <unordered_set>
//#include <cassert>
//#include <queue>
//#include <random>
//#include <sstream>


//namespace floatTetWild {
//
//#define NOT_SURFACE SCHAR_MAX//SHRT_MAX//INT_MAX
//#define KNOWN_NOT_SURFACE (-SCHAR_MAX/2)
//#define KNOWN_SURFACE (SCHAR_MAX/2)
//
//#define NO_SURFACE_TAG 0
//#define NOT_BBOX (-1)
//#define OPP_T_ID_UNKNOWN (-2)
//#define OPP_T_ID_BOUNDARY (-1)
//
//    class MeshVertex {
//    public:
//        explicit MeshVertex(const Eigen::Matrix<double, 3, 1>& p): pos(p){}
//        MeshVertex() = default;
//
//        Eigen::Matrix<double, 3, 1> pos;
//        //Eigen::Matrix<double, 3, 1> input_pos;
//
//        inline double &operator[](const int index) {
//            assert(index >= 0 && index < 3);
//            return pos[index];
//        }
//
//        inline double operator[](const int index) const {
//            assert(index >= 0 && index < 3);
//            return pos[index];
//        }
//
//        std::vector<int> conn_tets;
//
//        bool is_on_surface = false;
//        bool is_on_boundary = false;
//        bool is_on_cut = false;
//        int on_boundary_e_id = -1;
//        bool is_on_bbox = false;
//        bool is_outside = false;
//
//        bool is_removed = false;
//        bool is_freezed = false;//todo
//
//        double sizing_scalar = 1;
//
//        double scalar = 0;
//    };
//
//    class MeshTet {
//    public:
//        Eigen::Matrix<int, 4, 1> indices;
//
//        MeshTet() = default;
//        explicit MeshTet(const Eigen::Matrix<int, 4, 1> &idx) : indices(idx) {}
//        MeshTet(int v0, int v1, int v2, int v3) : indices(v0, v1, v2, v3) {}
//
//        inline void reset() {
//            is_surface_fs = {{NOT_SURFACE, NOT_SURFACE, NOT_SURFACE, NOT_SURFACE}};
//            is_bbox_fs = {{NOT_BBOX, NOT_BBOX, NOT_BBOX, NOT_BBOX}};
//            opp_t_ids = {{OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN}};
//            surface_tags = {{0, 0, 0, 0}};
//
//            quality = 0;
//            scalar = 0;
//            is_removed = false;
//            is_outside = false;
//        }
//
//        inline const Eigen::Matrix<int, 4, 1> &pts_indices() const { return indices; }
//
//        inline int &operator[](const int index) {
//            assert(index >= 0 && index < 4);
//            return indices[index];
//        }
//
//        inline int operator[](const int index) const {
//            assert(index >= 0 && index < 4);
//            return indices[index];
//        }
//
//        inline int find(int ele) const {
//            for (int j = 0; j < 4; j++)
//                if (indices[j] == ele)
//                    return j;
//            return -1;
//        }
//
//        inline int find_opp(int v0_id, int v1_id, int v2_id) const {
//            for (int j = 0; j < 4; j++)
//                if (indices[j] != v0_id && indices[j] != v1_id && indices[j] != v2_id)
//                    return j;
//            return -1;
//        }
//
//        inline std::string to_string() const
//        {
//            std::ostringstream ss;
//            ss << indices[0] << " " << indices[1] << " " << indices[2] << " " << indices[3];
//            return ss.str();
//        }
//
//        std::array<char, 4> is_surface_fs = {{NOT_SURFACE, NOT_SURFACE, NOT_SURFACE, NOT_SURFACE}};
//        std::array<char, 4> is_bbox_fs = {{NOT_BBOX, NOT_BBOX, NOT_BBOX, NOT_BBOX}};
//        std::array<int, 4> opp_t_ids = {{OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN}};
//        std::array<char, 4> surface_tags = {{0, 0, 0, 0}};
//
//        double quality = 0;
//        double scalar = 0;
//        bool is_removed = false;
//        bool is_outside = false;
//    };
//
//    class Mesh {
//    public:
//        std::vector<MeshVertex> tet_vertices;
//        std::vector<MeshTet> tets;
//        Parameters params;
//
//        int t_empty_start = 0;
//        int v_empty_start = 0;
//        bool is_limit_length = true;
//        bool is_closed = true;
//
//        bool is_input_all_inserted = false;
//        bool is_coarsening = false;
//
//        void one_ring_vertex_coloring(std::vector<double> &colors) const;
//
//        void one_ring_vertex_sets(const int threshold, std::vector<std::vector<int>> &concurrent_sets,
//                                  std::vector<int> &serial_set) const;
//
//        void partition(const int n_parts, std::vector<std::vector<int>> &tets_id) const;
//
//        static void one_ring_edge_set(const std::vector<std::array<int, 2>> &edges, const std::vector<bool> &v_is_removed,
//                          const std::vector<bool> &f_is_removed, const std::vector<std::unordered_set<int>> &conn_fs,
//                          const std::vector<Eigen::Matrix<double, 3, 1>> &input_vertices, std::vector<int> &safe_set);
//
//        inline int t_empty_size() const {
//            int cnt = 0;
//            for (const auto &t:tets) {
//                if (t.is_removed)
//                    cnt++;
//            }
//            return cnt;
//        }
//
//        inline int v_empty_size() const {
//            int cnt = 0;
//            for (const auto &v:tet_vertices) {
//                if (v.is_removed)
//                    cnt++;
//            }
//            return cnt;
//        }
//
//        inline void reset_t_empty_start() {
//            t_empty_start = tets.size();
//            for (int i = 0; i < tets.size(); i++) {
//                if (tets[i].is_removed) {
//                    t_empty_start = i;
//                    break;
//                }
//            }
//        }
//
//        inline void reset_v_empty_start() {
//            v_empty_start = tet_vertices.size();
//            for (int i = 0; i < tet_vertices.size(); i++) {
//                if (tet_vertices[i].is_removed) {
//                    v_empty_start = i;
//                    break;
//                }
//            }
//        }
//
//        inline int get_v_num() const {
//            int cnt = 0;
//            for (int i = 0; i < tet_vertices.size(); i++) {
//                if (!tet_vertices[i].is_removed)
//                    cnt++;
//            }
//            return cnt;
//        }
//
//        inline int get_t_num() const {
//            int cnt = 0;
//            for (int i = 0; i < tets.size(); i++) {
//                if (!tets[i].is_removed)
//                    cnt++;
//            }
//            return cnt;
//        }
//
//        inline double get_max_energy() const {
//            double max_energy = 0;
//            int cnt = 0;
//            for (auto &t: tets) {
//                if (t.is_removed)
//                    continue;
//                if (t.quality > max_energy)
//                    max_energy = t.quality;
//            }
//            return max_energy;
//        }
//
//        inline double get_avg_energy() const {
//            double avg_energy = 0;
//            int cnt = 0;
//            for (auto &t: tets) {
//                if (t.is_removed)
//                    continue;
//                avg_energy += t.quality;
//                cnt++;
//            }
//            avg_energy /= cnt;
//            return avg_energy;
//        }
//    };
//}

#endif  // FLOATTETWILD_MESH_HPP