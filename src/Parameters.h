// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

// To set the parameters related
#include <floattetwild/Types.hpp>

#include <array>
#include <vector>

#include <geogram/mesh/mesh.h>

namespace floatTetWild {
    class Parameters {
    public:
        std::string log_path = "";
        std::string input_path = "";
        std::string output_path = "";
        std::string tag_path = "";
        std::string postfix = "";

        std::string envelope_log = "";
        std::string envelope_log_csv = "";

        bool not_sort_input = false;
        bool correct_surface_orientation = false;
        bool is_quiet = false;
        int log_level = 3;//2;

        bool smooth_open_boundary = false;
        bool manifold_surface = false;
        bool disable_filtering = false;
        bool use_floodfill = false;
        bool use_general_wn = false;
        bool use_input_for_wn = false;
        bool coarsen = false;

        bool apply_sizing_field = false;
        Eigen::VectorXd V_sizing_field;
        Eigen::VectorXi T_sizing_field;
        Eigen::VectorXd values_sizing_field;
        std::function<double(const Eigen::Matrix<double, 3, 1>&)> get_sizing_field_value;//get sizing field value for an point

#ifdef NEW_ENVELOPE
        std::vector<double> input_epsr_tags;//same length as the list of input faces
#endif

        // it decides the scale of the box, presents the deviation of the box from the model
        //( in % of  max((xmax-xmin), (ymax-ymin), (zmax-zmin)) of the input points)
        double box_scale = 1 / 15.0;

        // epsilon presents the tolerence permited (in % of the box diagonal)
        double eps_rel = 1e-3;

        // initial target edge length at every vertex(in % of the box diagonal)
        double ideal_edge_length_rel = 1 / 20.0;
        double min_edge_len_rel = -1;

        int max_its = 80;
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
        double bbox_diag_length;
        double ideal_edge_length;
        double ideal_edge_length_2;
        double eps_input;
        double eps;
        double eps_delta;
        double eps_2;
        double dd;
        double min_edge_length;

        double split_threshold;
        double collapse_threshold;
        double split_threshold_2;
        double collapse_threshold_2;

        double eps_coplanar;
        double eps_2_coplanar;
        double eps_simplification;
        double eps_2_simplification;
        double dd_simplification;

        bool init(Scalar bbox_diag_l) {
            if (stage > 5)
                stage = 5;

            bbox_diag_length = bbox_diag_l;

            ideal_edge_length = bbox_diag_length * ideal_edge_length_rel;
            ideal_edge_length_2 = ideal_edge_length * ideal_edge_length;

            eps_input = bbox_diag_length * eps_rel;
            dd = eps_input;// / stage;
            dd /= 1.5;

#ifdef NEW_ENVELOPE
            double eps_usable = eps_input;
            eps_delta = eps_usable * 0.1;
            eps = eps_usable - eps_delta * (stage - 1);
#else
            double eps_usable = eps_input - dd / std::sqrt(3);
            eps_delta = eps_usable * 0.1;
            eps = eps_usable - eps_delta * (stage - 1);
#endif

//        dd /= 1.6;
//        eps_delta = dd / std::sqrt(3);
//        eps       = eps_input - eps_delta * stage;
//        dd /= 1.2;
            eps_2 = eps * eps;

            eps_coplanar = eps * 0.2;  // better to set it as eps-related
            if (eps_coplanar > bbox_diag_length * 1e-6)
                eps_coplanar = bbox_diag_length * 1e-6;
            eps_2_coplanar = eps_coplanar * eps_coplanar;

            eps_simplification = eps * 0.8;
            eps_2_simplification = eps_simplification * eps_simplification;
            dd_simplification = dd / eps * eps_simplification;
            //            dd_simplification = dd;

            if (min_edge_len_rel < 0)
                min_edge_len_rel = eps_rel;
            min_edge_length = bbox_diag_length * min_edge_len_rel;

            split_threshold = ideal_edge_length * (4 / 3.0);
            collapse_threshold = ideal_edge_length * (4 / 5.0);
            split_threshold_2 = split_threshold * split_threshold;
            collapse_threshold_2 = collapse_threshold * collapse_threshold;

            std::cout << "bbox_diag_length = " << bbox_diag_length << std::endl;
            std::cout << "ideal_edge_length = " << ideal_edge_length << std::endl;

            std::cout << "stage = " << stage << std::endl;
            std::cout << "eps_input = " << eps_input << std::endl;
            std::cout << "eps = " << eps << std::endl;
            std::cout << "eps_simplification = " << eps_simplification << std::endl;
            std::cout << "eps_coplanar = " << eps_coplanar << std::endl;

            std::cout << "dd = " << dd << std::endl;
            std::cout << "dd_simplification = " << dd_simplification << std::endl;

            return true;
        }
    };
}  // namespace floatTetWild
