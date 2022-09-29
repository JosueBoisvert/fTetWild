// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_INTERSECTIONS_H
#define FLOATTETWILD_INTERSECTIONS_H

//#include <floattetwild/Mesh.hpp>
#include "ftetwild.hpp"

namespace floatTetWild {
#define CUT_EDGE_0 0
#define CUT_EDGE_1 1
#define CUT_EDGE_2 2
#define CUT_FACE 3
#define CUT_COPLANAR 4
#define CUT_EMPTY -1

    int is_tri_tri_cutted(const std::array<Eigen::Matrix<double, 3, 1>, 3> &f_tri, const std::array<Eigen::Matrix<double, 3, 1>, 3> &f_tet,
                          const std::array<int, 3>& oris_tri);

    double seg_seg_squared_dist_3d(const std::array<Eigen::Matrix<double, 3, 1>, 2> &s1, const std::array<Eigen::Matrix<double, 3, 1>, 2> &s2);

    double p_seg_squared_dist_3d(const Eigen::Matrix<double, 3, 1> &p, const Eigen::Matrix<double, 3, 1> &a, const Eigen::Matrix<double, 3, 1> &b);
    double p_line_squared_dist_3d(const Eigen::Matrix<double, 3, 1> &p, const Eigen::Matrix<double, 3, 1> &a, const Eigen::Matrix<double, 3, 1> &b);

    bool is_p_inside_tri_2d(const Eigen::Matrix<double, 2, 1>& p, const std::array<Eigen::Matrix<double, 2, 1>, 3> &tri);
    bool is_seg_tri_cutted_2d(const std::array<Eigen::Matrix<double, 2, 1>, 2> &seg, const std::array<Eigen::Matrix<double, 2, 1>, 3> &tri);
    bool is_tri_tri_cutted_2d(const std::array<Eigen::Matrix<double, 2, 1>, 3> &p_tet, const std::array<Eigen::Matrix<double, 2, 1>, 3> &p_tri);

    bool seg_seg_intersection_2d(const std::array<Eigen::Matrix<double, 2, 1>, 2> &seg1, const std::array<Eigen::Matrix<double, 2, 1>, 2> &seg2, Scalar& t2);
    bool seg_line_intersection_2d(const std::array<Eigen::Matrix<double, 2, 1>, 2> &seg, const std::array<Eigen::Matrix<double, 2, 1>, 2> &line, Scalar& t_seg);
    bool seg_plane_intersection(const Eigen::Matrix<double, 3, 1> &p1, const Eigen::Matrix<double, 3, 1> &p2, const Eigen::Matrix<double, 3, 1> &a, const Eigen::Matrix<double, 3, 1> &n,
                                Eigen::Matrix<double, 3, 1> &p, double &d1);

    int get_t(const Eigen::Matrix<double, 3, 1> &p0, const Eigen::Matrix<double, 3, 1> &p1, const Eigen::Matrix<double, 3, 1> &p2);
    Eigen::Matrix<double, 2, 1> to_2d(const Eigen::Matrix<double, 3, 1> &p, int t);
    Eigen::Matrix<double, 2, 1> to_2d(const Eigen::Matrix<double, 3, 1> &p, const Eigen::Matrix<double, 3, 1>& n, const Eigen::Matrix<double, 3, 1>& pp, int t);

    bool is_crossing(int s1, int s2);

    int is_tri_tri_cutted(const Eigen::Matrix<double, 3, 1> &p1, const Eigen::Matrix<double, 3, 1> &p2, const Eigen::Matrix<double, 3, 1> &p3,//cutting tri
                          const Eigen::Matrix<double, 3, 1> &q1, const Eigen::Matrix<double, 3, 1> &q2, const Eigen::Matrix<double, 3, 1> &q3);//face of tet
    int is_tri_tri_cutted_hint(const Eigen::Matrix<double, 3, 1> &p1, const Eigen::Matrix<double, 3, 1> &p2, const Eigen::Matrix<double, 3, 1> &p3,//cutting tri
                               const Eigen::Matrix<double, 3, 1> &q1, const Eigen::Matrix<double, 3, 1> &q2, const Eigen::Matrix<double, 3, 1> &q3, int hint,
                               bool is_debug = false);//face of tet

    void get_bbox_face(const Eigen::Matrix<double, 3, 1>& p0, const Eigen::Matrix<double, 3, 1>& p1, const Eigen::Matrix<double, 3, 1>& p2, Eigen::Matrix<double, 3, 1>& min, Eigen::Matrix<double, 3, 1>& max, double eps = 0);
    void get_bbox_tet(const Eigen::Matrix<double, 3, 1>& p0, const Eigen::Matrix<double, 3, 1>& p1, const Eigen::Matrix<double, 3, 1>& p2, const Eigen::Matrix<double, 3, 1>& p3, Eigen::Matrix<double, 3, 1>& min, Eigen::Matrix<double, 3, 1>& max, double eps = 0);

    bool is_bbox_intersected(const Eigen::Matrix<double, 3, 1>& min1, const Eigen::Matrix<double, 3, 1>& max1, const Eigen::Matrix<double, 3, 1>& min2, const Eigen::Matrix<double, 3, 1>& max2);

    bool is_tri_inside_tet(const std::array<Eigen::Matrix<double, 3, 1>, 3>& ps,
                           const Eigen::Matrix<double, 3, 1>& p0t, const Eigen::Matrix<double, 3, 1>& p1t, const Eigen::Matrix<double, 3, 1>& p2t, const Eigen::Matrix<double, 3, 1>& p3t);
    bool is_point_inside_tet(const Eigen::Matrix<double, 3, 1>& p, const Eigen::Matrix<double, 3, 1>& p0t, const Eigen::Matrix<double, 3, 1>& p1t, const Eigen::Matrix<double, 3, 1>& p2t, const Eigen::Matrix<double, 3, 1>& p3t);
}

#endif //FLOATTETWILD_INTERSECTIONS_H
