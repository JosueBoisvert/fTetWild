// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include <floattetwild/Types.hpp>

#include <geogram/mesh/mesh.h>



#include <vector>
#include <map>
#include <string>

namespace floatTetWild {
class CSGTreeParser
{
    public:
    static void get_meshes(const nlohmann::json &csg_tree, std::vector<std::string> &meshes, nlohmann::json &csg_tree_with_ids)
    {
        int index = 0;
        meshes.clear();

        std::map<std::string, int> existings;

        get_meshes_aux(csg_tree, meshes, existings, index, csg_tree_with_ids);
    }

    static int get_max_id(const nlohmann::json &csg_tree_with_ids) {
        int max = -1;
        get_max_id_aux(csg_tree_with_ids, max);

        return max;
    }

    static bool keep_tet(const nlohmann::json &csg_tree_with_ids, const int t_id, const std::vector<Eigen::VectorXd> &w);

    static bool load_and_merge(const std::vector<std::string> &meshes, std::vector<Eigen::Matrix<double, 3, 1>> &V, std::vector<Eigen::Matrix<int, 3, 1>> &F, GEO::Mesh &sf_mesh, std::vector<int> &tags);
    static void merge(const std::vector<std::vector<Eigen::Matrix<double, 3, 1>>> &Vs, const std::vector<std::vector<Eigen::Matrix<int, 3, 1>>> &Fs, std::vector<Eigen::Matrix<double, 3, 1>> &V, std::vector<Eigen::Matrix<int, 3, 1>> &F, GEO::Mesh &sf_mesh, std::vector<int> &tags);

    private:
    static void get_meshes_aux(const nlohmann::json &csg_tree_node, std::vector<std::string> &meshes, std::map<std::string, int> &existings, int &index, nlohmann::json &current_node);
    static void get_max_id_aux(const nlohmann::json &csg_tree_node, int &max);
};

}
