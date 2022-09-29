// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "auto_table.hpp"

#include <cassert>

namespace floatTetWild {
	const std::vector<std::vector<Eigen::Matrix<int, 4, 1>>>& CutTable::get_tet_confs(const int idx) {
		static const std::array<std::vector<std::vector<Eigen::Matrix<int, 4, 1>>>, 64> table= {{

			{

				{
					Eigen::Matrix<int, 4, 1>(0, 1, 2, 3)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(4, 3, 2, 0),Eigen::Matrix<int, 4, 1>(4, 1, 2, 3)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(4, 3, 0, 1),Eigen::Matrix<int, 4, 1>(4, 2, 0, 3)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(1, 3, 5, 4),Eigen::Matrix<int, 4, 1>(4, 0, 3, 5),Eigen::Matrix<int, 4, 1>(5, 0, 3, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 3, 5, 4),Eigen::Matrix<int, 4, 1>(4, 2, 3, 5),Eigen::Matrix<int, 4, 1>(4, 0, 3, 2)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(4, 3, 1, 2),Eigen::Matrix<int, 4, 1>(4, 0, 1, 3)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(0, 4, 5, 3),Eigen::Matrix<int, 4, 1>(4, 5, 3, 1),Eigen::Matrix<int, 4, 1>(5, 2, 3, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 4, 5, 3),Eigen::Matrix<int, 4, 1>(4, 5, 3, 2),Eigen::Matrix<int, 4, 1>(4, 2, 3, 1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(2, 3, 5, 4),Eigen::Matrix<int, 4, 1>(4, 1, 3, 5),Eigen::Matrix<int, 4, 1>(5, 1, 3, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 3, 5, 4),Eigen::Matrix<int, 4, 1>(4, 0, 3, 5),Eigen::Matrix<int, 4, 1>(4, 1, 3, 0)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(4, 2, 1, 0),Eigen::Matrix<int, 4, 1>(4, 3, 1, 2)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(0, 2, 5, 4),Eigen::Matrix<int, 4, 1>(4, 1, 2, 5),Eigen::Matrix<int, 4, 1>(5, 1, 2, 3)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 2, 5, 4),Eigen::Matrix<int, 4, 1>(4, 3, 2, 5),Eigen::Matrix<int, 4, 1>(4, 1, 2, 3)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(4, 5, 0, 1),Eigen::Matrix<int, 4, 1>(4, 2, 0, 5),Eigen::Matrix<int, 4, 1>(4, 1, 3, 5),Eigen::Matrix<int, 4, 1>(4, 5, 3, 2)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(1, 7, 5, 4),Eigen::Matrix<int, 4, 1>(0, 7, 6, 4),Eigen::Matrix<int, 4, 1>(2, 3, 6, 7),Eigen::Matrix<int, 4, 1>(2, 7, 6, 0),Eigen::Matrix<int, 4, 1>(2, 7, 5, 3),Eigen::Matrix<int, 4, 1>(3, 7, 5, 1),Eigen::Matrix<int, 4, 1>(0, 4, 5, 7),Eigen::Matrix<int, 4, 1>(0, 5, 2, 7),Eigen::Matrix<int, 4, 1>(1, 4, 6, 7),Eigen::Matrix<int, 4, 1>(1, 6, 3, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 7, 5, 4),Eigen::Matrix<int, 4, 1>(0, 7, 6, 4),Eigen::Matrix<int, 4, 1>(2, 3, 6, 7),Eigen::Matrix<int, 4, 1>(2, 7, 6, 0),Eigen::Matrix<int, 4, 1>(2, 7, 5, 3),Eigen::Matrix<int, 4, 1>(3, 7, 5, 1),Eigen::Matrix<int, 4, 1>(0, 4, 2, 7),Eigen::Matrix<int, 4, 1>(4, 5, 2, 7),Eigen::Matrix<int, 4, 1>(1, 4, 6, 7),Eigen::Matrix<int, 4, 1>(1, 6, 3, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 7, 5, 4),Eigen::Matrix<int, 4, 1>(0, 7, 6, 4),Eigen::Matrix<int, 4, 1>(2, 3, 6, 7),Eigen::Matrix<int, 4, 1>(2, 7, 6, 0),Eigen::Matrix<int, 4, 1>(2, 7, 5, 3),Eigen::Matrix<int, 4, 1>(3, 7, 5, 1),Eigen::Matrix<int, 4, 1>(0, 4, 5, 7),Eigen::Matrix<int, 4, 1>(0, 5, 2, 7),Eigen::Matrix<int, 4, 1>(1, 4, 3, 7),Eigen::Matrix<int, 4, 1>(4, 6, 3, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 7, 5, 4),Eigen::Matrix<int, 4, 1>(0, 7, 6, 4),Eigen::Matrix<int, 4, 1>(2, 3, 6, 7),Eigen::Matrix<int, 4, 1>(2, 7, 6, 0),Eigen::Matrix<int, 4, 1>(2, 7, 5, 3),Eigen::Matrix<int, 4, 1>(3, 7, 5, 1),Eigen::Matrix<int, 4, 1>(0, 4, 2, 7),Eigen::Matrix<int, 4, 1>(4, 5, 2, 7),Eigen::Matrix<int, 4, 1>(1, 4, 3, 7),Eigen::Matrix<int, 4, 1>(4, 6, 3, 7)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(0, 4, 5, 1),Eigen::Matrix<int, 4, 1>(4, 5, 1, 2),Eigen::Matrix<int, 4, 1>(5, 3, 1, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 4, 5, 1),Eigen::Matrix<int, 4, 1>(4, 5, 1, 3),Eigen::Matrix<int, 4, 1>(4, 3, 1, 2)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(4, 0, 6, 5),Eigen::Matrix<int, 4, 1>(4, 5, 6, 2),Eigen::Matrix<int, 4, 1>(6, 2, 3, 1),Eigen::Matrix<int, 4, 1>(4, 2, 6, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 0, 6, 5),Eigen::Matrix<int, 4, 1>(4, 5, 6, 2),Eigen::Matrix<int, 4, 1>(6, 2, 3, 4),Eigen::Matrix<int, 4, 1>(4, 2, 3, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 0, 6, 5),Eigen::Matrix<int, 4, 1>(1, 5, 6, 2),Eigen::Matrix<int, 4, 1>(6, 2, 3, 1),Eigen::Matrix<int, 4, 1>(4, 5, 6, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 0, 6, 5),Eigen::Matrix<int, 4, 1>(4, 5, 3, 2),Eigen::Matrix<int, 4, 1>(6, 5, 3, 4),Eigen::Matrix<int, 4, 1>(4, 2, 3, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 0, 6, 5),Eigen::Matrix<int, 4, 1>(1, 5, 3, 2),Eigen::Matrix<int, 4, 1>(1, 5, 6, 3),Eigen::Matrix<int, 4, 1>(4, 5, 6, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 0, 6, 5),Eigen::Matrix<int, 4, 1>(1, 5, 3, 2),Eigen::Matrix<int, 4, 1>(4, 5, 6, 3),Eigen::Matrix<int, 4, 1>(4, 5, 3, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 0, 6, 5),Eigen::Matrix<int, 4, 1>(1, 4, 3, 7),Eigen::Matrix<int, 4, 1>(3, 7, 4, 6),Eigen::Matrix<int, 4, 1>(5, 7, 2, 6),Eigen::Matrix<int, 4, 1>(3, 6, 2, 7),Eigen::Matrix<int, 4, 1>(1, 7, 2, 5),Eigen::Matrix<int, 4, 1>(4, 7, 1, 5),Eigen::Matrix<int, 4, 1>(1, 7, 3, 2),Eigen::Matrix<int, 4, 1>(4, 5, 6, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 0, 6, 5),Eigen::Matrix<int, 4, 1>(1, 4, 6, 7),Eigen::Matrix<int, 4, 1>(1, 6, 3, 7),Eigen::Matrix<int, 4, 1>(5, 7, 3, 6),Eigen::Matrix<int, 4, 1>(3, 5, 2, 7),Eigen::Matrix<int, 4, 1>(4, 7, 1, 2),Eigen::Matrix<int, 4, 1>(4, 7, 2, 5),Eigen::Matrix<int, 4, 1>(1, 7, 3, 2),Eigen::Matrix<int, 4, 1>(4, 5, 6, 7)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(2, 7, 5, 4),Eigen::Matrix<int, 4, 1>(0, 5, 6, 7),Eigen::Matrix<int, 4, 1>(1, 7, 6, 3),Eigen::Matrix<int, 4, 1>(1, 0, 6, 7),Eigen::Matrix<int, 4, 1>(1, 3, 4, 7),Eigen::Matrix<int, 4, 1>(3, 2, 4, 7),Eigen::Matrix<int, 4, 1>(1, 4, 5, 7),Eigen::Matrix<int, 4, 1>(1, 5, 0, 7),Eigen::Matrix<int, 4, 1>(2, 7, 6, 5),Eigen::Matrix<int, 4, 1>(2, 7, 3, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 7, 5, 4),Eigen::Matrix<int, 4, 1>(0, 5, 6, 7),Eigen::Matrix<int, 4, 1>(1, 7, 6, 3),Eigen::Matrix<int, 4, 1>(1, 0, 6, 7),Eigen::Matrix<int, 4, 1>(1, 3, 4, 7),Eigen::Matrix<int, 4, 1>(3, 2, 4, 7),Eigen::Matrix<int, 4, 1>(1, 4, 0, 7),Eigen::Matrix<int, 4, 1>(4, 5, 0, 7),Eigen::Matrix<int, 4, 1>(2, 7, 6, 5),Eigen::Matrix<int, 4, 1>(2, 7, 3, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 7, 5, 4),Eigen::Matrix<int, 4, 1>(0, 5, 6, 7),Eigen::Matrix<int, 4, 1>(1, 7, 6, 3),Eigen::Matrix<int, 4, 1>(1, 0, 6, 7),Eigen::Matrix<int, 4, 1>(1, 3, 4, 7),Eigen::Matrix<int, 4, 1>(3, 2, 4, 7),Eigen::Matrix<int, 4, 1>(1, 4, 5, 7),Eigen::Matrix<int, 4, 1>(1, 5, 0, 7),Eigen::Matrix<int, 4, 1>(2, 7, 3, 5),Eigen::Matrix<int, 4, 1>(5, 7, 3, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 7, 5, 4),Eigen::Matrix<int, 4, 1>(0, 5, 6, 7),Eigen::Matrix<int, 4, 1>(1, 7, 6, 3),Eigen::Matrix<int, 4, 1>(1, 0, 6, 7),Eigen::Matrix<int, 4, 1>(1, 3, 4, 7),Eigen::Matrix<int, 4, 1>(3, 2, 4, 7),Eigen::Matrix<int, 4, 1>(1, 4, 0, 7),Eigen::Matrix<int, 4, 1>(4, 5, 0, 7),Eigen::Matrix<int, 4, 1>(2, 7, 3, 5),Eigen::Matrix<int, 4, 1>(5, 7, 3, 6)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(4, 0, 2, 1),Eigen::Matrix<int, 4, 1>(4, 3, 2, 0)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(1, 4, 5, 2),Eigen::Matrix<int, 4, 1>(4, 5, 2, 0),Eigen::Matrix<int, 4, 1>(5, 3, 2, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 4, 5, 2),Eigen::Matrix<int, 4, 1>(4, 5, 2, 3),Eigen::Matrix<int, 4, 1>(4, 3, 2, 0)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(1, 0, 5, 4),Eigen::Matrix<int, 4, 1>(4, 2, 0, 5),Eigen::Matrix<int, 4, 1>(5, 2, 0, 3)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 0, 5, 4),Eigen::Matrix<int, 4, 1>(4, 3, 0, 5),Eigen::Matrix<int, 4, 1>(4, 2, 0, 3)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 1),Eigen::Matrix<int, 4, 1>(4, 2, 6, 5),Eigen::Matrix<int, 4, 1>(6, 0, 3, 2),Eigen::Matrix<int, 4, 1>(4, 0, 6, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 1),Eigen::Matrix<int, 4, 1>(4, 2, 6, 5),Eigen::Matrix<int, 4, 1>(6, 4, 3, 2),Eigen::Matrix<int, 4, 1>(4, 0, 3, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 1),Eigen::Matrix<int, 4, 1>(0, 2, 6, 5),Eigen::Matrix<int, 4, 1>(6, 0, 3, 2),Eigen::Matrix<int, 4, 1>(4, 0, 6, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 1),Eigen::Matrix<int, 4, 1>(4, 2, 3, 5),Eigen::Matrix<int, 4, 1>(6, 4, 3, 5),Eigen::Matrix<int, 4, 1>(4, 0, 3, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 1),Eigen::Matrix<int, 4, 1>(0, 2, 3, 5),Eigen::Matrix<int, 4, 1>(0, 3, 6, 5),Eigen::Matrix<int, 4, 1>(4, 0, 6, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 1),Eigen::Matrix<int, 4, 1>(0, 2, 3, 5),Eigen::Matrix<int, 4, 1>(4, 3, 6, 5),Eigen::Matrix<int, 4, 1>(4, 0, 3, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 1),Eigen::Matrix<int, 4, 1>(0, 7, 3, 4),Eigen::Matrix<int, 4, 1>(3, 6, 4, 7),Eigen::Matrix<int, 4, 1>(5, 6, 2, 7),Eigen::Matrix<int, 4, 1>(3, 7, 2, 6),Eigen::Matrix<int, 4, 1>(0, 5, 2, 7),Eigen::Matrix<int, 4, 1>(4, 5, 0, 7),Eigen::Matrix<int, 4, 1>(0, 2, 3, 7),Eigen::Matrix<int, 4, 1>(4, 7, 6, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 1),Eigen::Matrix<int, 4, 1>(0, 7, 6, 4),Eigen::Matrix<int, 4, 1>(0, 7, 3, 6),Eigen::Matrix<int, 4, 1>(5, 6, 3, 7),Eigen::Matrix<int, 4, 1>(3, 7, 2, 5),Eigen::Matrix<int, 4, 1>(4, 2, 0, 7),Eigen::Matrix<int, 4, 1>(4, 5, 2, 7),Eigen::Matrix<int, 4, 1>(0, 2, 3, 7),Eigen::Matrix<int, 4, 1>(4, 7, 6, 5)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(4, 5, 1, 2),Eigen::Matrix<int, 4, 1>(4, 0, 1, 5),Eigen::Matrix<int, 4, 1>(4, 2, 3, 5),Eigen::Matrix<int, 4, 1>(4, 5, 3, 0)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(0, 4, 5, 7),Eigen::Matrix<int, 4, 1>(1, 4, 6, 7),Eigen::Matrix<int, 4, 1>(2, 3, 5, 7),Eigen::Matrix<int, 4, 1>(3, 0, 5, 7),Eigen::Matrix<int, 4, 1>(2, 7, 6, 3),Eigen::Matrix<int, 4, 1>(2, 1, 6, 7),Eigen::Matrix<int, 4, 1>(1, 7, 5, 4),Eigen::Matrix<int, 4, 1>(1, 7, 2, 5),Eigen::Matrix<int, 4, 1>(0, 7, 6, 4),Eigen::Matrix<int, 4, 1>(0, 7, 3, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 4, 5, 7),Eigen::Matrix<int, 4, 1>(1, 4, 6, 7),Eigen::Matrix<int, 4, 1>(2, 3, 5, 7),Eigen::Matrix<int, 4, 1>(3, 0, 5, 7),Eigen::Matrix<int, 4, 1>(2, 7, 6, 3),Eigen::Matrix<int, 4, 1>(2, 1, 6, 7),Eigen::Matrix<int, 4, 1>(1, 7, 2, 4),Eigen::Matrix<int, 4, 1>(4, 7, 2, 5),Eigen::Matrix<int, 4, 1>(0, 7, 6, 4),Eigen::Matrix<int, 4, 1>(0, 7, 3, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 4, 5, 7),Eigen::Matrix<int, 4, 1>(1, 4, 6, 7),Eigen::Matrix<int, 4, 1>(2, 3, 5, 7),Eigen::Matrix<int, 4, 1>(3, 0, 5, 7),Eigen::Matrix<int, 4, 1>(2, 7, 6, 3),Eigen::Matrix<int, 4, 1>(2, 1, 6, 7),Eigen::Matrix<int, 4, 1>(1, 7, 5, 4),Eigen::Matrix<int, 4, 1>(1, 7, 2, 5),Eigen::Matrix<int, 4, 1>(0, 7, 3, 4),Eigen::Matrix<int, 4, 1>(4, 7, 3, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 4, 5, 7),Eigen::Matrix<int, 4, 1>(1, 4, 6, 7),Eigen::Matrix<int, 4, 1>(2, 3, 5, 7),Eigen::Matrix<int, 4, 1>(3, 0, 5, 7),Eigen::Matrix<int, 4, 1>(2, 7, 6, 3),Eigen::Matrix<int, 4, 1>(2, 1, 6, 7),Eigen::Matrix<int, 4, 1>(1, 7, 2, 4),Eigen::Matrix<int, 4, 1>(4, 7, 2, 5),Eigen::Matrix<int, 4, 1>(0, 7, 3, 4),Eigen::Matrix<int, 4, 1>(4, 7, 3, 6)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(2, 7, 5, 4),Eigen::Matrix<int, 4, 1>(1, 7, 6, 4),Eigen::Matrix<int, 4, 1>(0, 3, 6, 7),Eigen::Matrix<int, 4, 1>(0, 7, 6, 1),Eigen::Matrix<int, 4, 1>(0, 7, 5, 3),Eigen::Matrix<int, 4, 1>(3, 7, 5, 2),Eigen::Matrix<int, 4, 1>(1, 4, 5, 7),Eigen::Matrix<int, 4, 1>(1, 5, 0, 7),Eigen::Matrix<int, 4, 1>(2, 4, 6, 7),Eigen::Matrix<int, 4, 1>(2, 6, 3, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 7, 5, 4),Eigen::Matrix<int, 4, 1>(1, 7, 6, 4),Eigen::Matrix<int, 4, 1>(0, 3, 6, 7),Eigen::Matrix<int, 4, 1>(0, 7, 6, 1),Eigen::Matrix<int, 4, 1>(0, 7, 5, 3),Eigen::Matrix<int, 4, 1>(3, 7, 5, 2),Eigen::Matrix<int, 4, 1>(1, 4, 0, 7),Eigen::Matrix<int, 4, 1>(4, 5, 0, 7),Eigen::Matrix<int, 4, 1>(2, 4, 6, 7),Eigen::Matrix<int, 4, 1>(2, 6, 3, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 7, 5, 4),Eigen::Matrix<int, 4, 1>(1, 7, 6, 4),Eigen::Matrix<int, 4, 1>(0, 3, 6, 7),Eigen::Matrix<int, 4, 1>(0, 7, 6, 1),Eigen::Matrix<int, 4, 1>(0, 7, 5, 3),Eigen::Matrix<int, 4, 1>(3, 7, 5, 2),Eigen::Matrix<int, 4, 1>(1, 4, 5, 7),Eigen::Matrix<int, 4, 1>(1, 5, 0, 7),Eigen::Matrix<int, 4, 1>(2, 4, 3, 7),Eigen::Matrix<int, 4, 1>(4, 6, 3, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 7, 5, 4),Eigen::Matrix<int, 4, 1>(1, 7, 6, 4),Eigen::Matrix<int, 4, 1>(0, 3, 6, 7),Eigen::Matrix<int, 4, 1>(0, 7, 6, 1),Eigen::Matrix<int, 4, 1>(0, 7, 5, 3),Eigen::Matrix<int, 4, 1>(3, 7, 5, 2),Eigen::Matrix<int, 4, 1>(1, 4, 0, 7),Eigen::Matrix<int, 4, 1>(4, 5, 0, 7),Eigen::Matrix<int, 4, 1>(2, 4, 3, 7),Eigen::Matrix<int, 4, 1>(4, 6, 3, 7)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(3, 2, 5, 4),Eigen::Matrix<int, 4, 1>(4, 0, 2, 5),Eigen::Matrix<int, 4, 1>(5, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(3, 2, 5, 4),Eigen::Matrix<int, 4, 1>(4, 1, 2, 5),Eigen::Matrix<int, 4, 1>(4, 0, 2, 1)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(1, 7, 6, 4),Eigen::Matrix<int, 4, 1>(3, 7, 6, 5),Eigen::Matrix<int, 4, 1>(0, 7, 4, 2),Eigen::Matrix<int, 4, 1>(0, 1, 4, 7),Eigen::Matrix<int, 4, 1>(0, 2, 5, 7),Eigen::Matrix<int, 4, 1>(2, 3, 5, 7),Eigen::Matrix<int, 4, 1>(2, 4, 6, 7),Eigen::Matrix<int, 4, 1>(2, 6, 3, 7),Eigen::Matrix<int, 4, 1>(0, 5, 6, 7),Eigen::Matrix<int, 4, 1>(0, 6, 1, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 7, 6, 4),Eigen::Matrix<int, 4, 1>(3, 7, 6, 5),Eigen::Matrix<int, 4, 1>(0, 7, 4, 2),Eigen::Matrix<int, 4, 1>(0, 1, 4, 7),Eigen::Matrix<int, 4, 1>(0, 2, 5, 7),Eigen::Matrix<int, 4, 1>(2, 3, 5, 7),Eigen::Matrix<int, 4, 1>(2, 4, 3, 7),Eigen::Matrix<int, 4, 1>(4, 6, 3, 7),Eigen::Matrix<int, 4, 1>(0, 5, 6, 7),Eigen::Matrix<int, 4, 1>(0, 6, 1, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 7, 6, 4),Eigen::Matrix<int, 4, 1>(3, 7, 6, 5),Eigen::Matrix<int, 4, 1>(0, 7, 4, 2),Eigen::Matrix<int, 4, 1>(0, 1, 4, 7),Eigen::Matrix<int, 4, 1>(0, 2, 5, 7),Eigen::Matrix<int, 4, 1>(2, 3, 5, 7),Eigen::Matrix<int, 4, 1>(2, 4, 6, 7),Eigen::Matrix<int, 4, 1>(2, 6, 3, 7),Eigen::Matrix<int, 4, 1>(0, 5, 1, 7),Eigen::Matrix<int, 4, 1>(5, 6, 1, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 7, 6, 4),Eigen::Matrix<int, 4, 1>(3, 7, 6, 5),Eigen::Matrix<int, 4, 1>(0, 7, 4, 2),Eigen::Matrix<int, 4, 1>(0, 1, 4, 7),Eigen::Matrix<int, 4, 1>(0, 2, 5, 7),Eigen::Matrix<int, 4, 1>(2, 3, 5, 7),Eigen::Matrix<int, 4, 1>(2, 4, 3, 7),Eigen::Matrix<int, 4, 1>(4, 6, 3, 7),Eigen::Matrix<int, 4, 1>(0, 5, 1, 7),Eigen::Matrix<int, 4, 1>(5, 6, 1, 7)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(0, 4, 5, 7),Eigen::Matrix<int, 4, 1>(3, 7, 6, 5),Eigen::Matrix<int, 4, 1>(1, 2, 4, 7),Eigen::Matrix<int, 4, 1>(1, 7, 4, 0),Eigen::Matrix<int, 4, 1>(1, 7, 6, 2),Eigen::Matrix<int, 4, 1>(2, 7, 6, 3),Eigen::Matrix<int, 4, 1>(2, 7, 5, 4),Eigen::Matrix<int, 4, 1>(2, 7, 3, 5),Eigen::Matrix<int, 4, 1>(0, 5, 6, 7),Eigen::Matrix<int, 4, 1>(0, 6, 1, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 4, 5, 7),Eigen::Matrix<int, 4, 1>(3, 7, 6, 5),Eigen::Matrix<int, 4, 1>(1, 2, 4, 7),Eigen::Matrix<int, 4, 1>(1, 7, 4, 0),Eigen::Matrix<int, 4, 1>(1, 7, 6, 2),Eigen::Matrix<int, 4, 1>(2, 7, 6, 3),Eigen::Matrix<int, 4, 1>(2, 7, 3, 4),Eigen::Matrix<int, 4, 1>(4, 7, 3, 5),Eigen::Matrix<int, 4, 1>(0, 5, 6, 7),Eigen::Matrix<int, 4, 1>(0, 6, 1, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 4, 5, 7),Eigen::Matrix<int, 4, 1>(3, 7, 6, 5),Eigen::Matrix<int, 4, 1>(1, 2, 4, 7),Eigen::Matrix<int, 4, 1>(1, 7, 4, 0),Eigen::Matrix<int, 4, 1>(1, 7, 6, 2),Eigen::Matrix<int, 4, 1>(2, 7, 6, 3),Eigen::Matrix<int, 4, 1>(2, 7, 5, 4),Eigen::Matrix<int, 4, 1>(2, 7, 3, 5),Eigen::Matrix<int, 4, 1>(0, 5, 1, 7),Eigen::Matrix<int, 4, 1>(5, 6, 1, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 4, 5, 7),Eigen::Matrix<int, 4, 1>(3, 7, 6, 5),Eigen::Matrix<int, 4, 1>(1, 2, 4, 7),Eigen::Matrix<int, 4, 1>(1, 7, 4, 0),Eigen::Matrix<int, 4, 1>(1, 7, 6, 2),Eigen::Matrix<int, 4, 1>(2, 7, 6, 3),Eigen::Matrix<int, 4, 1>(2, 7, 3, 4),Eigen::Matrix<int, 4, 1>(4, 7, 3, 5),Eigen::Matrix<int, 4, 1>(0, 5, 1, 7),Eigen::Matrix<int, 4, 1>(5, 6, 1, 7)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(1, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 0, 6, 5),Eigen::Matrix<int, 4, 1>(1, 0, 7, 5),Eigen::Matrix<int, 4, 1>(3, 5, 7, 6),Eigen::Matrix<int, 4, 1>(7, 2, 4, 5),Eigen::Matrix<int, 4, 1>(3, 2, 7, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 0, 6, 5),Eigen::Matrix<int, 4, 1>(1, 0, 7, 5),Eigen::Matrix<int, 4, 1>(3, 5, 7, 6),Eigen::Matrix<int, 4, 1>(7, 3, 4, 5),Eigen::Matrix<int, 4, 1>(3, 2, 4, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 0, 6, 5),Eigen::Matrix<int, 4, 1>(1, 0, 7, 5),Eigen::Matrix<int, 4, 1>(2, 5, 7, 6),Eigen::Matrix<int, 4, 1>(7, 2, 4, 5),Eigen::Matrix<int, 4, 1>(3, 2, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 0, 6, 5),Eigen::Matrix<int, 4, 1>(1, 0, 7, 5),Eigen::Matrix<int, 4, 1>(2, 8, 4, 3),Eigen::Matrix<int, 4, 1>(4, 7, 3, 8),Eigen::Matrix<int, 4, 1>(6, 7, 5, 8),Eigen::Matrix<int, 4, 1>(4, 8, 5, 7),Eigen::Matrix<int, 4, 1>(2, 6, 5, 8),Eigen::Matrix<int, 4, 1>(3, 6, 2, 8),Eigen::Matrix<int, 4, 1>(2, 5, 4, 8),Eigen::Matrix<int, 4, 1>(3, 8, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 1, 6, 5),Eigen::Matrix<int, 4, 1>(1, 0, 6, 5),Eigen::Matrix<int, 4, 1>(3, 5, 7, 6),Eigen::Matrix<int, 4, 1>(7, 2, 4, 5),Eigen::Matrix<int, 4, 1>(3, 2, 7, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 1, 6, 5),Eigen::Matrix<int, 4, 1>(1, 0, 6, 5),Eigen::Matrix<int, 4, 1>(3, 5, 7, 6),Eigen::Matrix<int, 4, 1>(7, 3, 4, 5),Eigen::Matrix<int, 4, 1>(3, 2, 4, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 1, 6, 5),Eigen::Matrix<int, 4, 1>(1, 0, 6, 5),Eigen::Matrix<int, 4, 1>(2, 5, 7, 6),Eigen::Matrix<int, 4, 1>(7, 2, 4, 5),Eigen::Matrix<int, 4, 1>(3, 2, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 1, 6, 5),Eigen::Matrix<int, 4, 1>(1, 0, 6, 5),Eigen::Matrix<int, 4, 1>(2, 8, 4, 3),Eigen::Matrix<int, 4, 1>(4, 7, 3, 8),Eigen::Matrix<int, 4, 1>(6, 7, 5, 8),Eigen::Matrix<int, 4, 1>(4, 8, 5, 7),Eigen::Matrix<int, 4, 1>(2, 6, 5, 8),Eigen::Matrix<int, 4, 1>(3, 6, 2, 8),Eigen::Matrix<int, 4, 1>(2, 5, 4, 8),Eigen::Matrix<int, 4, 1>(3, 8, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 0, 6, 5),Eigen::Matrix<int, 4, 1>(1, 0, 7, 4),Eigen::Matrix<int, 4, 1>(3, 5, 7, 6),Eigen::Matrix<int, 4, 1>(7, 2, 4, 5),Eigen::Matrix<int, 4, 1>(3, 2, 7, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 0, 6, 5),Eigen::Matrix<int, 4, 1>(1, 0, 7, 4),Eigen::Matrix<int, 4, 1>(3, 5, 7, 6),Eigen::Matrix<int, 4, 1>(7, 3, 4, 5),Eigen::Matrix<int, 4, 1>(3, 2, 4, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 0, 6, 5),Eigen::Matrix<int, 4, 1>(1, 0, 7, 4),Eigen::Matrix<int, 4, 1>(2, 5, 7, 6),Eigen::Matrix<int, 4, 1>(7, 2, 4, 5),Eigen::Matrix<int, 4, 1>(3, 2, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 0, 6, 5),Eigen::Matrix<int, 4, 1>(1, 0, 7, 4),Eigen::Matrix<int, 4, 1>(2, 8, 4, 3),Eigen::Matrix<int, 4, 1>(4, 7, 3, 8),Eigen::Matrix<int, 4, 1>(6, 7, 5, 8),Eigen::Matrix<int, 4, 1>(4, 8, 5, 7),Eigen::Matrix<int, 4, 1>(2, 6, 5, 8),Eigen::Matrix<int, 4, 1>(3, 6, 2, 8),Eigen::Matrix<int, 4, 1>(2, 5, 4, 8),Eigen::Matrix<int, 4, 1>(3, 8, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 5, 6, 4),Eigen::Matrix<int, 4, 1>(7, 1, 6, 4),Eigen::Matrix<int, 4, 1>(1, 0, 6, 5),Eigen::Matrix<int, 4, 1>(3, 5, 4, 6),Eigen::Matrix<int, 4, 1>(7, 3, 4, 6),Eigen::Matrix<int, 4, 1>(3, 2, 4, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 5, 6, 4),Eigen::Matrix<int, 4, 1>(7, 1, 6, 4),Eigen::Matrix<int, 4, 1>(1, 0, 6, 5),Eigen::Matrix<int, 4, 1>(2, 5, 4, 6),Eigen::Matrix<int, 4, 1>(2, 4, 7, 6),Eigen::Matrix<int, 4, 1>(3, 2, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 5, 6, 4),Eigen::Matrix<int, 4, 1>(7, 1, 6, 4),Eigen::Matrix<int, 4, 1>(1, 0, 6, 5),Eigen::Matrix<int, 4, 1>(2, 5, 4, 6),Eigen::Matrix<int, 4, 1>(3, 4, 7, 6),Eigen::Matrix<int, 4, 1>(3, 2, 4, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 5, 6, 4),Eigen::Matrix<int, 4, 1>(7, 1, 6, 4),Eigen::Matrix<int, 4, 1>(1, 0, 6, 5),Eigen::Matrix<int, 4, 1>(2, 8, 7, 3),Eigen::Matrix<int, 4, 1>(2, 8, 4, 7),Eigen::Matrix<int, 4, 1>(6, 7, 4, 8),Eigen::Matrix<int, 4, 1>(4, 8, 5, 6),Eigen::Matrix<int, 4, 1>(3, 5, 2, 8),Eigen::Matrix<int, 4, 1>(3, 6, 5, 8),Eigen::Matrix<int, 4, 1>(2, 5, 4, 8),Eigen::Matrix<int, 4, 1>(3, 8, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 6, 4),Eigen::Matrix<int, 4, 1>(0, 6, 7, 4),Eigen::Matrix<int, 4, 1>(1, 0, 7, 4),Eigen::Matrix<int, 4, 1>(3, 5, 4, 6),Eigen::Matrix<int, 4, 1>(7, 3, 4, 6),Eigen::Matrix<int, 4, 1>(3, 2, 4, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 6, 4),Eigen::Matrix<int, 4, 1>(0, 6, 7, 4),Eigen::Matrix<int, 4, 1>(1, 0, 7, 4),Eigen::Matrix<int, 4, 1>(2, 5, 4, 6),Eigen::Matrix<int, 4, 1>(2, 4, 7, 6),Eigen::Matrix<int, 4, 1>(3, 2, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 6, 4),Eigen::Matrix<int, 4, 1>(0, 6, 7, 4),Eigen::Matrix<int, 4, 1>(1, 0, 7, 4),Eigen::Matrix<int, 4, 1>(2, 5, 4, 6),Eigen::Matrix<int, 4, 1>(3, 4, 7, 6),Eigen::Matrix<int, 4, 1>(3, 2, 4, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 6, 4),Eigen::Matrix<int, 4, 1>(0, 6, 7, 4),Eigen::Matrix<int, 4, 1>(1, 0, 7, 4),Eigen::Matrix<int, 4, 1>(2, 8, 7, 3),Eigen::Matrix<int, 4, 1>(2, 8, 4, 7),Eigen::Matrix<int, 4, 1>(6, 7, 4, 8),Eigen::Matrix<int, 4, 1>(4, 8, 5, 6),Eigen::Matrix<int, 4, 1>(3, 5, 2, 8),Eigen::Matrix<int, 4, 1>(3, 6, 5, 8),Eigen::Matrix<int, 4, 1>(2, 5, 4, 8),Eigen::Matrix<int, 4, 1>(3, 8, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 6, 4),Eigen::Matrix<int, 4, 1>(1, 6, 7, 4),Eigen::Matrix<int, 4, 1>(1, 0, 6, 4),Eigen::Matrix<int, 4, 1>(3, 5, 4, 6),Eigen::Matrix<int, 4, 1>(7, 3, 4, 6),Eigen::Matrix<int, 4, 1>(3, 2, 4, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 6, 4),Eigen::Matrix<int, 4, 1>(1, 6, 7, 4),Eigen::Matrix<int, 4, 1>(1, 0, 6, 4),Eigen::Matrix<int, 4, 1>(2, 5, 4, 6),Eigen::Matrix<int, 4, 1>(2, 4, 7, 6),Eigen::Matrix<int, 4, 1>(3, 2, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 6, 4),Eigen::Matrix<int, 4, 1>(1, 6, 7, 4),Eigen::Matrix<int, 4, 1>(1, 0, 6, 4),Eigen::Matrix<int, 4, 1>(2, 5, 4, 6),Eigen::Matrix<int, 4, 1>(3, 4, 7, 6),Eigen::Matrix<int, 4, 1>(3, 2, 4, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 6, 4),Eigen::Matrix<int, 4, 1>(1, 6, 7, 4),Eigen::Matrix<int, 4, 1>(1, 0, 6, 4),Eigen::Matrix<int, 4, 1>(2, 8, 7, 3),Eigen::Matrix<int, 4, 1>(2, 8, 4, 7),Eigen::Matrix<int, 4, 1>(6, 7, 4, 8),Eigen::Matrix<int, 4, 1>(4, 8, 5, 6),Eigen::Matrix<int, 4, 1>(3, 5, 2, 8),Eigen::Matrix<int, 4, 1>(3, 6, 5, 8),Eigen::Matrix<int, 4, 1>(2, 5, 4, 8),Eigen::Matrix<int, 4, 1>(3, 8, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 8, 6, 1),Eigen::Matrix<int, 4, 1>(6, 7, 1, 8),Eigen::Matrix<int, 4, 1>(4, 7, 5, 8),Eigen::Matrix<int, 4, 1>(6, 8, 5, 7),Eigen::Matrix<int, 4, 1>(0, 4, 5, 8),Eigen::Matrix<int, 4, 1>(1, 4, 0, 8),Eigen::Matrix<int, 4, 1>(0, 5, 6, 8),Eigen::Matrix<int, 4, 1>(1, 8, 7, 4),Eigen::Matrix<int, 4, 1>(3, 5, 7, 6),Eigen::Matrix<int, 4, 1>(7, 2, 4, 5),Eigen::Matrix<int, 4, 1>(3, 2, 7, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 8, 6, 1),Eigen::Matrix<int, 4, 1>(6, 7, 1, 8),Eigen::Matrix<int, 4, 1>(4, 7, 5, 8),Eigen::Matrix<int, 4, 1>(6, 8, 5, 7),Eigen::Matrix<int, 4, 1>(0, 4, 5, 8),Eigen::Matrix<int, 4, 1>(1, 4, 0, 8),Eigen::Matrix<int, 4, 1>(0, 5, 6, 8),Eigen::Matrix<int, 4, 1>(1, 8, 7, 4),Eigen::Matrix<int, 4, 1>(3, 5, 7, 6),Eigen::Matrix<int, 4, 1>(7, 3, 4, 5),Eigen::Matrix<int, 4, 1>(3, 2, 4, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 8, 6, 1),Eigen::Matrix<int, 4, 1>(6, 7, 1, 8),Eigen::Matrix<int, 4, 1>(4, 7, 5, 8),Eigen::Matrix<int, 4, 1>(6, 8, 5, 7),Eigen::Matrix<int, 4, 1>(0, 4, 5, 8),Eigen::Matrix<int, 4, 1>(1, 4, 0, 8),Eigen::Matrix<int, 4, 1>(0, 5, 6, 8),Eigen::Matrix<int, 4, 1>(1, 8, 7, 4),Eigen::Matrix<int, 4, 1>(2, 5, 7, 6),Eigen::Matrix<int, 4, 1>(7, 2, 4, 5),Eigen::Matrix<int, 4, 1>(3, 2, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 8, 6, 1),Eigen::Matrix<int, 4, 1>(6, 7, 1, 8),Eigen::Matrix<int, 4, 1>(4, 7, 5, 8),Eigen::Matrix<int, 4, 1>(6, 8, 5, 7),Eigen::Matrix<int, 4, 1>(0, 4, 5, 8),Eigen::Matrix<int, 4, 1>(1, 4, 0, 8),Eigen::Matrix<int, 4, 1>(0, 5, 6, 8),Eigen::Matrix<int, 4, 1>(1, 8, 7, 4),Eigen::Matrix<int, 4, 1>(2, 9, 4, 3),Eigen::Matrix<int, 4, 1>(4, 7, 3, 9),Eigen::Matrix<int, 4, 1>(6, 7, 5, 9),Eigen::Matrix<int, 4, 1>(4, 9, 5, 7),Eigen::Matrix<int, 4, 1>(2, 6, 5, 9),Eigen::Matrix<int, 4, 1>(3, 6, 2, 9),Eigen::Matrix<int, 4, 1>(2, 5, 4, 9),Eigen::Matrix<int, 4, 1>(3, 9, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 8, 7, 1),Eigen::Matrix<int, 4, 1>(0, 8, 6, 7),Eigen::Matrix<int, 4, 1>(4, 7, 6, 8),Eigen::Matrix<int, 4, 1>(6, 8, 5, 4),Eigen::Matrix<int, 4, 1>(1, 5, 0, 8),Eigen::Matrix<int, 4, 1>(1, 4, 5, 8),Eigen::Matrix<int, 4, 1>(0, 5, 6, 8),Eigen::Matrix<int, 4, 1>(1, 8, 7, 4),Eigen::Matrix<int, 4, 1>(3, 5, 4, 6),Eigen::Matrix<int, 4, 1>(7, 3, 4, 6),Eigen::Matrix<int, 4, 1>(3, 2, 4, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 8, 7, 1),Eigen::Matrix<int, 4, 1>(0, 8, 6, 7),Eigen::Matrix<int, 4, 1>(4, 7, 6, 8),Eigen::Matrix<int, 4, 1>(6, 8, 5, 4),Eigen::Matrix<int, 4, 1>(1, 5, 0, 8),Eigen::Matrix<int, 4, 1>(1, 4, 5, 8),Eigen::Matrix<int, 4, 1>(0, 5, 6, 8),Eigen::Matrix<int, 4, 1>(1, 8, 7, 4),Eigen::Matrix<int, 4, 1>(2, 5, 4, 6),Eigen::Matrix<int, 4, 1>(2, 4, 7, 6),Eigen::Matrix<int, 4, 1>(3, 2, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 8, 7, 1),Eigen::Matrix<int, 4, 1>(0, 8, 6, 7),Eigen::Matrix<int, 4, 1>(4, 7, 6, 8),Eigen::Matrix<int, 4, 1>(6, 8, 5, 4),Eigen::Matrix<int, 4, 1>(1, 5, 0, 8),Eigen::Matrix<int, 4, 1>(1, 4, 5, 8),Eigen::Matrix<int, 4, 1>(0, 5, 6, 8),Eigen::Matrix<int, 4, 1>(1, 8, 7, 4),Eigen::Matrix<int, 4, 1>(2, 5, 4, 6),Eigen::Matrix<int, 4, 1>(3, 4, 7, 6),Eigen::Matrix<int, 4, 1>(3, 2, 4, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 8, 7, 1),Eigen::Matrix<int, 4, 1>(0, 8, 6, 7),Eigen::Matrix<int, 4, 1>(4, 7, 6, 8),Eigen::Matrix<int, 4, 1>(6, 8, 5, 4),Eigen::Matrix<int, 4, 1>(1, 5, 0, 8),Eigen::Matrix<int, 4, 1>(1, 4, 5, 8),Eigen::Matrix<int, 4, 1>(0, 5, 6, 8),Eigen::Matrix<int, 4, 1>(1, 8, 7, 4),Eigen::Matrix<int, 4, 1>(2, 9, 7, 3),Eigen::Matrix<int, 4, 1>(2, 9, 4, 7),Eigen::Matrix<int, 4, 1>(6, 7, 4, 9),Eigen::Matrix<int, 4, 1>(4, 9, 5, 6),Eigen::Matrix<int, 4, 1>(3, 5, 2, 9),Eigen::Matrix<int, 4, 1>(3, 6, 5, 9),Eigen::Matrix<int, 4, 1>(2, 5, 4, 9),Eigen::Matrix<int, 4, 1>(3, 9, 7, 6)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(4, 1, 0, 2),Eigen::Matrix<int, 4, 1>(4, 3, 0, 1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(4, 5, 2, 0),Eigen::Matrix<int, 4, 1>(4, 1, 2, 5),Eigen::Matrix<int, 4, 1>(4, 0, 3, 5),Eigen::Matrix<int, 4, 1>(4, 5, 3, 1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(2, 4, 5, 0),Eigen::Matrix<int, 4, 1>(4, 5, 0, 1),Eigen::Matrix<int, 4, 1>(5, 3, 0, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 4, 5, 0),Eigen::Matrix<int, 4, 1>(4, 5, 0, 3),Eigen::Matrix<int, 4, 1>(4, 3, 0, 1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(1, 7, 5, 4),Eigen::Matrix<int, 4, 1>(2, 5, 6, 7),Eigen::Matrix<int, 4, 1>(0, 3, 4, 7),Eigen::Matrix<int, 4, 1>(3, 1, 4, 7),Eigen::Matrix<int, 4, 1>(0, 7, 6, 3),Eigen::Matrix<int, 4, 1>(0, 2, 6, 7),Eigen::Matrix<int, 4, 1>(0, 4, 5, 7),Eigen::Matrix<int, 4, 1>(0, 5, 2, 7),Eigen::Matrix<int, 4, 1>(1, 7, 6, 5),Eigen::Matrix<int, 4, 1>(1, 7, 3, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 7, 5, 4),Eigen::Matrix<int, 4, 1>(2, 5, 6, 7),Eigen::Matrix<int, 4, 1>(0, 3, 4, 7),Eigen::Matrix<int, 4, 1>(3, 1, 4, 7),Eigen::Matrix<int, 4, 1>(0, 7, 6, 3),Eigen::Matrix<int, 4, 1>(0, 2, 6, 7),Eigen::Matrix<int, 4, 1>(0, 4, 2, 7),Eigen::Matrix<int, 4, 1>(4, 5, 2, 7),Eigen::Matrix<int, 4, 1>(1, 7, 6, 5),Eigen::Matrix<int, 4, 1>(1, 7, 3, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 7, 5, 4),Eigen::Matrix<int, 4, 1>(2, 5, 6, 7),Eigen::Matrix<int, 4, 1>(0, 3, 4, 7),Eigen::Matrix<int, 4, 1>(3, 1, 4, 7),Eigen::Matrix<int, 4, 1>(0, 7, 6, 3),Eigen::Matrix<int, 4, 1>(0, 2, 6, 7),Eigen::Matrix<int, 4, 1>(0, 4, 5, 7),Eigen::Matrix<int, 4, 1>(0, 5, 2, 7),Eigen::Matrix<int, 4, 1>(1, 7, 3, 5),Eigen::Matrix<int, 4, 1>(5, 7, 3, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 7, 5, 4),Eigen::Matrix<int, 4, 1>(2, 5, 6, 7),Eigen::Matrix<int, 4, 1>(0, 3, 4, 7),Eigen::Matrix<int, 4, 1>(3, 1, 4, 7),Eigen::Matrix<int, 4, 1>(0, 7, 6, 3),Eigen::Matrix<int, 4, 1>(0, 2, 6, 7),Eigen::Matrix<int, 4, 1>(0, 4, 2, 7),Eigen::Matrix<int, 4, 1>(4, 5, 2, 7),Eigen::Matrix<int, 4, 1>(1, 7, 3, 5),Eigen::Matrix<int, 4, 1>(5, 7, 3, 6)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(2, 1, 5, 4),Eigen::Matrix<int, 4, 1>(4, 0, 1, 5),Eigen::Matrix<int, 4, 1>(5, 0, 1, 3)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 1, 5, 4),Eigen::Matrix<int, 4, 1>(4, 3, 1, 5),Eigen::Matrix<int, 4, 1>(4, 0, 1, 3)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(0, 4, 5, 7),Eigen::Matrix<int, 4, 1>(2, 7, 6, 5),Eigen::Matrix<int, 4, 1>(1, 7, 4, 3),Eigen::Matrix<int, 4, 1>(3, 7, 4, 0),Eigen::Matrix<int, 4, 1>(1, 3, 6, 7),Eigen::Matrix<int, 4, 1>(1, 7, 6, 2),Eigen::Matrix<int, 4, 1>(1, 7, 5, 4),Eigen::Matrix<int, 4, 1>(1, 7, 2, 5),Eigen::Matrix<int, 4, 1>(0, 5, 6, 7),Eigen::Matrix<int, 4, 1>(0, 6, 3, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 4, 5, 7),Eigen::Matrix<int, 4, 1>(2, 7, 6, 5),Eigen::Matrix<int, 4, 1>(1, 7, 4, 3),Eigen::Matrix<int, 4, 1>(3, 7, 4, 0),Eigen::Matrix<int, 4, 1>(1, 3, 6, 7),Eigen::Matrix<int, 4, 1>(1, 7, 6, 2),Eigen::Matrix<int, 4, 1>(1, 7, 2, 4),Eigen::Matrix<int, 4, 1>(4, 7, 2, 5),Eigen::Matrix<int, 4, 1>(0, 5, 6, 7),Eigen::Matrix<int, 4, 1>(0, 6, 3, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 4, 5, 7),Eigen::Matrix<int, 4, 1>(2, 7, 6, 5),Eigen::Matrix<int, 4, 1>(1, 7, 4, 3),Eigen::Matrix<int, 4, 1>(3, 7, 4, 0),Eigen::Matrix<int, 4, 1>(1, 3, 6, 7),Eigen::Matrix<int, 4, 1>(1, 7, 6, 2),Eigen::Matrix<int, 4, 1>(1, 7, 5, 4),Eigen::Matrix<int, 4, 1>(1, 7, 2, 5),Eigen::Matrix<int, 4, 1>(0, 5, 3, 7),Eigen::Matrix<int, 4, 1>(5, 6, 3, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 4, 5, 7),Eigen::Matrix<int, 4, 1>(2, 7, 6, 5),Eigen::Matrix<int, 4, 1>(1, 7, 4, 3),Eigen::Matrix<int, 4, 1>(3, 7, 4, 0),Eigen::Matrix<int, 4, 1>(1, 3, 6, 7),Eigen::Matrix<int, 4, 1>(1, 7, 6, 2),Eigen::Matrix<int, 4, 1>(1, 7, 2, 4),Eigen::Matrix<int, 4, 1>(4, 7, 2, 5),Eigen::Matrix<int, 4, 1>(0, 5, 3, 7),Eigen::Matrix<int, 4, 1>(5, 6, 3, 7)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 2),Eigen::Matrix<int, 4, 1>(4, 0, 6, 5),Eigen::Matrix<int, 4, 1>(6, 1, 3, 0),Eigen::Matrix<int, 4, 1>(4, 1, 6, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 2),Eigen::Matrix<int, 4, 1>(4, 0, 6, 5),Eigen::Matrix<int, 4, 1>(6, 4, 3, 0),Eigen::Matrix<int, 4, 1>(4, 1, 3, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 2),Eigen::Matrix<int, 4, 1>(1, 0, 6, 5),Eigen::Matrix<int, 4, 1>(6, 1, 3, 0),Eigen::Matrix<int, 4, 1>(4, 1, 6, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 2),Eigen::Matrix<int, 4, 1>(4, 0, 3, 5),Eigen::Matrix<int, 4, 1>(6, 4, 3, 5),Eigen::Matrix<int, 4, 1>(4, 1, 3, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 2),Eigen::Matrix<int, 4, 1>(1, 0, 3, 5),Eigen::Matrix<int, 4, 1>(1, 3, 6, 5),Eigen::Matrix<int, 4, 1>(4, 1, 6, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 2),Eigen::Matrix<int, 4, 1>(1, 0, 3, 5),Eigen::Matrix<int, 4, 1>(4, 3, 6, 5),Eigen::Matrix<int, 4, 1>(4, 1, 3, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 2),Eigen::Matrix<int, 4, 1>(1, 7, 3, 4),Eigen::Matrix<int, 4, 1>(3, 6, 4, 7),Eigen::Matrix<int, 4, 1>(5, 6, 0, 7),Eigen::Matrix<int, 4, 1>(3, 7, 0, 6),Eigen::Matrix<int, 4, 1>(1, 5, 0, 7),Eigen::Matrix<int, 4, 1>(4, 5, 1, 7),Eigen::Matrix<int, 4, 1>(1, 0, 3, 7),Eigen::Matrix<int, 4, 1>(4, 7, 6, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 2),Eigen::Matrix<int, 4, 1>(1, 7, 6, 4),Eigen::Matrix<int, 4, 1>(1, 7, 3, 6),Eigen::Matrix<int, 4, 1>(5, 6, 3, 7),Eigen::Matrix<int, 4, 1>(3, 7, 0, 5),Eigen::Matrix<int, 4, 1>(4, 0, 1, 7),Eigen::Matrix<int, 4, 1>(4, 5, 0, 7),Eigen::Matrix<int, 4, 1>(1, 0, 3, 7),Eigen::Matrix<int, 4, 1>(4, 7, 6, 5)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(3, 4, 5, 1),Eigen::Matrix<int, 4, 1>(4, 5, 1, 0),Eigen::Matrix<int, 4, 1>(5, 2, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(3, 4, 5, 1),Eigen::Matrix<int, 4, 1>(4, 5, 1, 2),Eigen::Matrix<int, 4, 1>(4, 2, 1, 0)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(0, 7, 5, 4),Eigen::Matrix<int, 4, 1>(3, 5, 6, 7),Eigen::Matrix<int, 4, 1>(1, 2, 4, 7),Eigen::Matrix<int, 4, 1>(2, 0, 4, 7),Eigen::Matrix<int, 4, 1>(1, 7, 6, 2),Eigen::Matrix<int, 4, 1>(1, 3, 6, 7),Eigen::Matrix<int, 4, 1>(1, 4, 5, 7),Eigen::Matrix<int, 4, 1>(1, 5, 3, 7),Eigen::Matrix<int, 4, 1>(0, 7, 6, 5),Eigen::Matrix<int, 4, 1>(0, 7, 2, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 7, 5, 4),Eigen::Matrix<int, 4, 1>(3, 5, 6, 7),Eigen::Matrix<int, 4, 1>(1, 2, 4, 7),Eigen::Matrix<int, 4, 1>(2, 0, 4, 7),Eigen::Matrix<int, 4, 1>(1, 7, 6, 2),Eigen::Matrix<int, 4, 1>(1, 3, 6, 7),Eigen::Matrix<int, 4, 1>(1, 4, 3, 7),Eigen::Matrix<int, 4, 1>(4, 5, 3, 7),Eigen::Matrix<int, 4, 1>(0, 7, 6, 5),Eigen::Matrix<int, 4, 1>(0, 7, 2, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 7, 5, 4),Eigen::Matrix<int, 4, 1>(3, 5, 6, 7),Eigen::Matrix<int, 4, 1>(1, 2, 4, 7),Eigen::Matrix<int, 4, 1>(2, 0, 4, 7),Eigen::Matrix<int, 4, 1>(1, 7, 6, 2),Eigen::Matrix<int, 4, 1>(1, 3, 6, 7),Eigen::Matrix<int, 4, 1>(1, 4, 5, 7),Eigen::Matrix<int, 4, 1>(1, 5, 3, 7),Eigen::Matrix<int, 4, 1>(0, 7, 2, 5),Eigen::Matrix<int, 4, 1>(5, 7, 2, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 7, 5, 4),Eigen::Matrix<int, 4, 1>(3, 5, 6, 7),Eigen::Matrix<int, 4, 1>(1, 2, 4, 7),Eigen::Matrix<int, 4, 1>(2, 0, 4, 7),Eigen::Matrix<int, 4, 1>(1, 7, 6, 2),Eigen::Matrix<int, 4, 1>(1, 3, 6, 7),Eigen::Matrix<int, 4, 1>(1, 4, 3, 7),Eigen::Matrix<int, 4, 1>(4, 5, 3, 7),Eigen::Matrix<int, 4, 1>(0, 7, 2, 5),Eigen::Matrix<int, 4, 1>(5, 7, 2, 6)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(2, 4, 6, 7),Eigen::Matrix<int, 4, 1>(3, 5, 6, 7),Eigen::Matrix<int, 4, 1>(0, 1, 4, 7),Eigen::Matrix<int, 4, 1>(0, 7, 4, 2),Eigen::Matrix<int, 4, 1>(0, 7, 5, 1),Eigen::Matrix<int, 4, 1>(1, 7, 5, 3),Eigen::Matrix<int, 4, 1>(1, 7, 6, 4),Eigen::Matrix<int, 4, 1>(1, 7, 3, 6),Eigen::Matrix<int, 4, 1>(0, 7, 6, 5),Eigen::Matrix<int, 4, 1>(0, 7, 2, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 4, 6, 7),Eigen::Matrix<int, 4, 1>(3, 5, 6, 7),Eigen::Matrix<int, 4, 1>(0, 1, 4, 7),Eigen::Matrix<int, 4, 1>(0, 7, 4, 2),Eigen::Matrix<int, 4, 1>(0, 7, 5, 1),Eigen::Matrix<int, 4, 1>(1, 7, 5, 3),Eigen::Matrix<int, 4, 1>(1, 7, 3, 4),Eigen::Matrix<int, 4, 1>(4, 7, 3, 6),Eigen::Matrix<int, 4, 1>(0, 7, 6, 5),Eigen::Matrix<int, 4, 1>(0, 7, 2, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 4, 6, 7),Eigen::Matrix<int, 4, 1>(3, 5, 6, 7),Eigen::Matrix<int, 4, 1>(0, 1, 4, 7),Eigen::Matrix<int, 4, 1>(0, 7, 4, 2),Eigen::Matrix<int, 4, 1>(0, 7, 5, 1),Eigen::Matrix<int, 4, 1>(1, 7, 5, 3),Eigen::Matrix<int, 4, 1>(1, 7, 6, 4),Eigen::Matrix<int, 4, 1>(1, 7, 3, 6),Eigen::Matrix<int, 4, 1>(0, 7, 2, 5),Eigen::Matrix<int, 4, 1>(5, 7, 2, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 4, 6, 7),Eigen::Matrix<int, 4, 1>(3, 5, 6, 7),Eigen::Matrix<int, 4, 1>(0, 1, 4, 7),Eigen::Matrix<int, 4, 1>(0, 7, 4, 2),Eigen::Matrix<int, 4, 1>(0, 7, 5, 1),Eigen::Matrix<int, 4, 1>(1, 7, 5, 3),Eigen::Matrix<int, 4, 1>(1, 7, 3, 4),Eigen::Matrix<int, 4, 1>(4, 7, 3, 6),Eigen::Matrix<int, 4, 1>(0, 7, 2, 5),Eigen::Matrix<int, 4, 1>(5, 7, 2, 6)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(2, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 6, 0),Eigen::Matrix<int, 4, 1>(2, 4, 7, 0),Eigen::Matrix<int, 4, 1>(3, 6, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 5, 1),Eigen::Matrix<int, 4, 1>(3, 4, 7, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 6, 0),Eigen::Matrix<int, 4, 1>(2, 4, 7, 0),Eigen::Matrix<int, 4, 1>(3, 6, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 5, 3),Eigen::Matrix<int, 4, 1>(3, 4, 5, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 6, 0),Eigen::Matrix<int, 4, 1>(2, 4, 7, 0),Eigen::Matrix<int, 4, 1>(1, 6, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 5, 1),Eigen::Matrix<int, 4, 1>(3, 6, 7, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 6, 0),Eigen::Matrix<int, 4, 1>(2, 4, 7, 0),Eigen::Matrix<int, 4, 1>(1, 3, 5, 8),Eigen::Matrix<int, 4, 1>(5, 8, 3, 7),Eigen::Matrix<int, 4, 1>(6, 8, 4, 7),Eigen::Matrix<int, 4, 1>(5, 7, 4, 8),Eigen::Matrix<int, 4, 1>(1, 8, 4, 6),Eigen::Matrix<int, 4, 1>(3, 8, 1, 6),Eigen::Matrix<int, 4, 1>(1, 8, 5, 4),Eigen::Matrix<int, 4, 1>(3, 6, 7, 8)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 6, 2),Eigen::Matrix<int, 4, 1>(2, 4, 6, 0),Eigen::Matrix<int, 4, 1>(3, 6, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 5, 1),Eigen::Matrix<int, 4, 1>(3, 4, 7, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 6, 2),Eigen::Matrix<int, 4, 1>(2, 4, 6, 0),Eigen::Matrix<int, 4, 1>(3, 6, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 5, 3),Eigen::Matrix<int, 4, 1>(3, 4, 5, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 6, 2),Eigen::Matrix<int, 4, 1>(2, 4, 6, 0),Eigen::Matrix<int, 4, 1>(1, 6, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 5, 1),Eigen::Matrix<int, 4, 1>(3, 6, 7, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 6, 2),Eigen::Matrix<int, 4, 1>(2, 4, 6, 0),Eigen::Matrix<int, 4, 1>(1, 3, 5, 8),Eigen::Matrix<int, 4, 1>(5, 8, 3, 7),Eigen::Matrix<int, 4, 1>(6, 8, 4, 7),Eigen::Matrix<int, 4, 1>(5, 7, 4, 8),Eigen::Matrix<int, 4, 1>(1, 8, 4, 6),Eigen::Matrix<int, 4, 1>(3, 8, 1, 6),Eigen::Matrix<int, 4, 1>(1, 8, 5, 4),Eigen::Matrix<int, 4, 1>(3, 6, 7, 8)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 6, 0),Eigen::Matrix<int, 4, 1>(2, 5, 7, 0),Eigen::Matrix<int, 4, 1>(3, 6, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 5, 1),Eigen::Matrix<int, 4, 1>(3, 4, 7, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 6, 0),Eigen::Matrix<int, 4, 1>(2, 5, 7, 0),Eigen::Matrix<int, 4, 1>(3, 6, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 5, 3),Eigen::Matrix<int, 4, 1>(3, 4, 5, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 6, 0),Eigen::Matrix<int, 4, 1>(2, 5, 7, 0),Eigen::Matrix<int, 4, 1>(1, 6, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 5, 1),Eigen::Matrix<int, 4, 1>(3, 6, 7, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 6, 0),Eigen::Matrix<int, 4, 1>(2, 5, 7, 0),Eigen::Matrix<int, 4, 1>(1, 3, 5, 8),Eigen::Matrix<int, 4, 1>(5, 8, 3, 7),Eigen::Matrix<int, 4, 1>(6, 8, 4, 7),Eigen::Matrix<int, 4, 1>(5, 7, 4, 8),Eigen::Matrix<int, 4, 1>(1, 8, 4, 6),Eigen::Matrix<int, 4, 1>(3, 8, 1, 6),Eigen::Matrix<int, 4, 1>(1, 8, 5, 4),Eigen::Matrix<int, 4, 1>(3, 6, 7, 8)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 5, 6, 4),Eigen::Matrix<int, 4, 1>(7, 5, 6, 2),Eigen::Matrix<int, 4, 1>(2, 4, 6, 0),Eigen::Matrix<int, 4, 1>(3, 6, 5, 4),Eigen::Matrix<int, 4, 1>(7, 6, 5, 3),Eigen::Matrix<int, 4, 1>(3, 4, 5, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 5, 6, 4),Eigen::Matrix<int, 4, 1>(7, 5, 6, 2),Eigen::Matrix<int, 4, 1>(2, 4, 6, 0),Eigen::Matrix<int, 4, 1>(1, 6, 5, 4),Eigen::Matrix<int, 4, 1>(1, 6, 7, 5),Eigen::Matrix<int, 4, 1>(3, 6, 7, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 5, 6, 4),Eigen::Matrix<int, 4, 1>(7, 5, 6, 2),Eigen::Matrix<int, 4, 1>(2, 4, 6, 0),Eigen::Matrix<int, 4, 1>(1, 6, 5, 4),Eigen::Matrix<int, 4, 1>(3, 6, 7, 5),Eigen::Matrix<int, 4, 1>(3, 6, 5, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 5, 6, 4),Eigen::Matrix<int, 4, 1>(7, 5, 6, 2),Eigen::Matrix<int, 4, 1>(2, 4, 6, 0),Eigen::Matrix<int, 4, 1>(1, 3, 7, 8),Eigen::Matrix<int, 4, 1>(1, 7, 5, 8),Eigen::Matrix<int, 4, 1>(6, 8, 5, 7),Eigen::Matrix<int, 4, 1>(5, 6, 4, 8),Eigen::Matrix<int, 4, 1>(3, 8, 1, 4),Eigen::Matrix<int, 4, 1>(3, 8, 4, 6),Eigen::Matrix<int, 4, 1>(1, 8, 5, 4),Eigen::Matrix<int, 4, 1>(3, 6, 7, 8)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 6, 4),Eigen::Matrix<int, 4, 1>(0, 5, 7, 6),Eigen::Matrix<int, 4, 1>(2, 5, 7, 0),Eigen::Matrix<int, 4, 1>(3, 6, 5, 4),Eigen::Matrix<int, 4, 1>(7, 6, 5, 3),Eigen::Matrix<int, 4, 1>(3, 4, 5, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 6, 4),Eigen::Matrix<int, 4, 1>(0, 5, 7, 6),Eigen::Matrix<int, 4, 1>(2, 5, 7, 0),Eigen::Matrix<int, 4, 1>(1, 6, 5, 4),Eigen::Matrix<int, 4, 1>(1, 6, 7, 5),Eigen::Matrix<int, 4, 1>(3, 6, 7, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 6, 4),Eigen::Matrix<int, 4, 1>(0, 5, 7, 6),Eigen::Matrix<int, 4, 1>(2, 5, 7, 0),Eigen::Matrix<int, 4, 1>(1, 6, 5, 4),Eigen::Matrix<int, 4, 1>(3, 6, 7, 5),Eigen::Matrix<int, 4, 1>(3, 6, 5, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 6, 4),Eigen::Matrix<int, 4, 1>(0, 5, 7, 6),Eigen::Matrix<int, 4, 1>(2, 5, 7, 0),Eigen::Matrix<int, 4, 1>(1, 3, 7, 8),Eigen::Matrix<int, 4, 1>(1, 7, 5, 8),Eigen::Matrix<int, 4, 1>(6, 8, 5, 7),Eigen::Matrix<int, 4, 1>(5, 6, 4, 8),Eigen::Matrix<int, 4, 1>(3, 8, 1, 4),Eigen::Matrix<int, 4, 1>(3, 8, 4, 6),Eigen::Matrix<int, 4, 1>(1, 8, 5, 4),Eigen::Matrix<int, 4, 1>(3, 6, 7, 8)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 6, 4),Eigen::Matrix<int, 4, 1>(2, 5, 7, 6),Eigen::Matrix<int, 4, 1>(2, 5, 6, 0),Eigen::Matrix<int, 4, 1>(3, 6, 5, 4),Eigen::Matrix<int, 4, 1>(7, 6, 5, 3),Eigen::Matrix<int, 4, 1>(3, 4, 5, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 6, 4),Eigen::Matrix<int, 4, 1>(2, 5, 7, 6),Eigen::Matrix<int, 4, 1>(2, 5, 6, 0),Eigen::Matrix<int, 4, 1>(1, 6, 5, 4),Eigen::Matrix<int, 4, 1>(1, 6, 7, 5),Eigen::Matrix<int, 4, 1>(3, 6, 7, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 6, 4),Eigen::Matrix<int, 4, 1>(2, 5, 7, 6),Eigen::Matrix<int, 4, 1>(2, 5, 6, 0),Eigen::Matrix<int, 4, 1>(1, 6, 5, 4),Eigen::Matrix<int, 4, 1>(3, 6, 7, 5),Eigen::Matrix<int, 4, 1>(3, 6, 5, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 5, 6, 4),Eigen::Matrix<int, 4, 1>(2, 5, 7, 6),Eigen::Matrix<int, 4, 1>(2, 5, 6, 0),Eigen::Matrix<int, 4, 1>(1, 3, 7, 8),Eigen::Matrix<int, 4, 1>(1, 7, 5, 8),Eigen::Matrix<int, 4, 1>(6, 8, 5, 7),Eigen::Matrix<int, 4, 1>(5, 6, 4, 8),Eigen::Matrix<int, 4, 1>(3, 8, 1, 4),Eigen::Matrix<int, 4, 1>(3, 8, 4, 6),Eigen::Matrix<int, 4, 1>(1, 8, 5, 4),Eigen::Matrix<int, 4, 1>(3, 6, 7, 8)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 2, 6, 8),Eigen::Matrix<int, 4, 1>(6, 8, 2, 7),Eigen::Matrix<int, 4, 1>(5, 8, 4, 7),Eigen::Matrix<int, 4, 1>(6, 7, 4, 8),Eigen::Matrix<int, 4, 1>(0, 8, 4, 5),Eigen::Matrix<int, 4, 1>(2, 8, 0, 5),Eigen::Matrix<int, 4, 1>(0, 8, 6, 4),Eigen::Matrix<int, 4, 1>(2, 5, 7, 8),Eigen::Matrix<int, 4, 1>(3, 6, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 5, 1),Eigen::Matrix<int, 4, 1>(3, 4, 7, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 2, 6, 8),Eigen::Matrix<int, 4, 1>(6, 8, 2, 7),Eigen::Matrix<int, 4, 1>(5, 8, 4, 7),Eigen::Matrix<int, 4, 1>(6, 7, 4, 8),Eigen::Matrix<int, 4, 1>(0, 8, 4, 5),Eigen::Matrix<int, 4, 1>(2, 8, 0, 5),Eigen::Matrix<int, 4, 1>(0, 8, 6, 4),Eigen::Matrix<int, 4, 1>(2, 5, 7, 8),Eigen::Matrix<int, 4, 1>(3, 6, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 5, 3),Eigen::Matrix<int, 4, 1>(3, 4, 5, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 2, 6, 8),Eigen::Matrix<int, 4, 1>(6, 8, 2, 7),Eigen::Matrix<int, 4, 1>(5, 8, 4, 7),Eigen::Matrix<int, 4, 1>(6, 7, 4, 8),Eigen::Matrix<int, 4, 1>(0, 8, 4, 5),Eigen::Matrix<int, 4, 1>(2, 8, 0, 5),Eigen::Matrix<int, 4, 1>(0, 8, 6, 4),Eigen::Matrix<int, 4, 1>(2, 5, 7, 8),Eigen::Matrix<int, 4, 1>(1, 6, 7, 4),Eigen::Matrix<int, 4, 1>(7, 4, 5, 1),Eigen::Matrix<int, 4, 1>(3, 6, 7, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 2, 6, 8),Eigen::Matrix<int, 4, 1>(6, 8, 2, 7),Eigen::Matrix<int, 4, 1>(5, 8, 4, 7),Eigen::Matrix<int, 4, 1>(6, 7, 4, 8),Eigen::Matrix<int, 4, 1>(0, 8, 4, 5),Eigen::Matrix<int, 4, 1>(2, 8, 0, 5),Eigen::Matrix<int, 4, 1>(0, 8, 6, 4),Eigen::Matrix<int, 4, 1>(2, 5, 7, 8),Eigen::Matrix<int, 4, 1>(1, 3, 5, 9),Eigen::Matrix<int, 4, 1>(5, 9, 3, 7),Eigen::Matrix<int, 4, 1>(6, 9, 4, 7),Eigen::Matrix<int, 4, 1>(5, 7, 4, 9),Eigen::Matrix<int, 4, 1>(1, 9, 4, 6),Eigen::Matrix<int, 4, 1>(3, 9, 1, 6),Eigen::Matrix<int, 4, 1>(1, 9, 5, 4),Eigen::Matrix<int, 4, 1>(3, 6, 7, 9)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 2, 7, 8),Eigen::Matrix<int, 4, 1>(0, 7, 6, 8),Eigen::Matrix<int, 4, 1>(5, 8, 6, 7),Eigen::Matrix<int, 4, 1>(6, 5, 4, 8),Eigen::Matrix<int, 4, 1>(2, 8, 0, 4),Eigen::Matrix<int, 4, 1>(2, 8, 4, 5),Eigen::Matrix<int, 4, 1>(0, 8, 6, 4),Eigen::Matrix<int, 4, 1>(2, 5, 7, 8),Eigen::Matrix<int, 4, 1>(3, 6, 5, 4),Eigen::Matrix<int, 4, 1>(7, 6, 5, 3),Eigen::Matrix<int, 4, 1>(3, 4, 5, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 2, 7, 8),Eigen::Matrix<int, 4, 1>(0, 7, 6, 8),Eigen::Matrix<int, 4, 1>(5, 8, 6, 7),Eigen::Matrix<int, 4, 1>(6, 5, 4, 8),Eigen::Matrix<int, 4, 1>(2, 8, 0, 4),Eigen::Matrix<int, 4, 1>(2, 8, 4, 5),Eigen::Matrix<int, 4, 1>(0, 8, 6, 4),Eigen::Matrix<int, 4, 1>(2, 5, 7, 8),Eigen::Matrix<int, 4, 1>(1, 6, 5, 4),Eigen::Matrix<int, 4, 1>(1, 6, 7, 5),Eigen::Matrix<int, 4, 1>(3, 6, 7, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 2, 7, 8),Eigen::Matrix<int, 4, 1>(0, 7, 6, 8),Eigen::Matrix<int, 4, 1>(5, 8, 6, 7),Eigen::Matrix<int, 4, 1>(6, 5, 4, 8),Eigen::Matrix<int, 4, 1>(2, 8, 0, 4),Eigen::Matrix<int, 4, 1>(2, 8, 4, 5),Eigen::Matrix<int, 4, 1>(0, 8, 6, 4),Eigen::Matrix<int, 4, 1>(2, 5, 7, 8),Eigen::Matrix<int, 4, 1>(1, 6, 5, 4),Eigen::Matrix<int, 4, 1>(3, 6, 7, 5),Eigen::Matrix<int, 4, 1>(3, 6, 5, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 2, 7, 8),Eigen::Matrix<int, 4, 1>(0, 7, 6, 8),Eigen::Matrix<int, 4, 1>(5, 8, 6, 7),Eigen::Matrix<int, 4, 1>(6, 5, 4, 8),Eigen::Matrix<int, 4, 1>(2, 8, 0, 4),Eigen::Matrix<int, 4, 1>(2, 8, 4, 5),Eigen::Matrix<int, 4, 1>(0, 8, 6, 4),Eigen::Matrix<int, 4, 1>(2, 5, 7, 8),Eigen::Matrix<int, 4, 1>(1, 3, 7, 9),Eigen::Matrix<int, 4, 1>(1, 7, 5, 9),Eigen::Matrix<int, 4, 1>(6, 9, 5, 7),Eigen::Matrix<int, 4, 1>(5, 6, 4, 9),Eigen::Matrix<int, 4, 1>(3, 9, 1, 4),Eigen::Matrix<int, 4, 1>(3, 9, 4, 6),Eigen::Matrix<int, 4, 1>(1, 9, 5, 4),Eigen::Matrix<int, 4, 1>(3, 6, 7, 9)
				}
			},
			{},
			{},
			{},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(3, 0, 5, 4),Eigen::Matrix<int, 4, 1>(4, 1, 0, 5),Eigen::Matrix<int, 4, 1>(5, 1, 0, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(3, 0, 5, 4),Eigen::Matrix<int, 4, 1>(4, 2, 0, 5),Eigen::Matrix<int, 4, 1>(4, 1, 0, 2)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(1, 4, 5, 7),Eigen::Matrix<int, 4, 1>(3, 7, 6, 5),Eigen::Matrix<int, 4, 1>(0, 7, 4, 2),Eigen::Matrix<int, 4, 1>(2, 7, 4, 1),Eigen::Matrix<int, 4, 1>(0, 2, 6, 7),Eigen::Matrix<int, 4, 1>(0, 7, 6, 3),Eigen::Matrix<int, 4, 1>(0, 7, 5, 4),Eigen::Matrix<int, 4, 1>(0, 7, 3, 5),Eigen::Matrix<int, 4, 1>(1, 5, 6, 7),Eigen::Matrix<int, 4, 1>(1, 6, 2, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 4, 5, 7),Eigen::Matrix<int, 4, 1>(3, 7, 6, 5),Eigen::Matrix<int, 4, 1>(0, 7, 4, 2),Eigen::Matrix<int, 4, 1>(2, 7, 4, 1),Eigen::Matrix<int, 4, 1>(0, 2, 6, 7),Eigen::Matrix<int, 4, 1>(0, 7, 6, 3),Eigen::Matrix<int, 4, 1>(0, 7, 3, 4),Eigen::Matrix<int, 4, 1>(4, 7, 3, 5),Eigen::Matrix<int, 4, 1>(1, 5, 6, 7),Eigen::Matrix<int, 4, 1>(1, 6, 2, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 4, 5, 7),Eigen::Matrix<int, 4, 1>(3, 7, 6, 5),Eigen::Matrix<int, 4, 1>(0, 7, 4, 2),Eigen::Matrix<int, 4, 1>(2, 7, 4, 1),Eigen::Matrix<int, 4, 1>(0, 2, 6, 7),Eigen::Matrix<int, 4, 1>(0, 7, 6, 3),Eigen::Matrix<int, 4, 1>(0, 7, 5, 4),Eigen::Matrix<int, 4, 1>(0, 7, 3, 5),Eigen::Matrix<int, 4, 1>(1, 5, 2, 7),Eigen::Matrix<int, 4, 1>(5, 6, 2, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 4, 5, 7),Eigen::Matrix<int, 4, 1>(3, 7, 6, 5),Eigen::Matrix<int, 4, 1>(0, 7, 4, 2),Eigen::Matrix<int, 4, 1>(2, 7, 4, 1),Eigen::Matrix<int, 4, 1>(0, 2, 6, 7),Eigen::Matrix<int, 4, 1>(0, 7, 6, 3),Eigen::Matrix<int, 4, 1>(0, 7, 3, 4),Eigen::Matrix<int, 4, 1>(4, 7, 3, 5),Eigen::Matrix<int, 4, 1>(1, 5, 2, 7),Eigen::Matrix<int, 4, 1>(5, 6, 2, 7)
				}
			},
			{},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(2, 7, 6, 4),Eigen::Matrix<int, 4, 1>(3, 7, 6, 5),Eigen::Matrix<int, 4, 1>(0, 1, 4, 7),Eigen::Matrix<int, 4, 1>(1, 2, 4, 7),Eigen::Matrix<int, 4, 1>(0, 7, 5, 1),Eigen::Matrix<int, 4, 1>(0, 3, 5, 7),Eigen::Matrix<int, 4, 1>(0, 4, 6, 7),Eigen::Matrix<int, 4, 1>(0, 6, 3, 7),Eigen::Matrix<int, 4, 1>(1, 5, 6, 7),Eigen::Matrix<int, 4, 1>(1, 6, 2, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 7, 6, 4),Eigen::Matrix<int, 4, 1>(3, 7, 6, 5),Eigen::Matrix<int, 4, 1>(0, 1, 4, 7),Eigen::Matrix<int, 4, 1>(1, 2, 4, 7),Eigen::Matrix<int, 4, 1>(0, 7, 5, 1),Eigen::Matrix<int, 4, 1>(0, 3, 5, 7),Eigen::Matrix<int, 4, 1>(0, 4, 3, 7),Eigen::Matrix<int, 4, 1>(4, 6, 3, 7),Eigen::Matrix<int, 4, 1>(1, 5, 6, 7),Eigen::Matrix<int, 4, 1>(1, 6, 2, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 7, 6, 4),Eigen::Matrix<int, 4, 1>(3, 7, 6, 5),Eigen::Matrix<int, 4, 1>(0, 1, 4, 7),Eigen::Matrix<int, 4, 1>(1, 2, 4, 7),Eigen::Matrix<int, 4, 1>(0, 7, 5, 1),Eigen::Matrix<int, 4, 1>(0, 3, 5, 7),Eigen::Matrix<int, 4, 1>(0, 4, 6, 7),Eigen::Matrix<int, 4, 1>(0, 6, 3, 7),Eigen::Matrix<int, 4, 1>(1, 5, 2, 7),Eigen::Matrix<int, 4, 1>(5, 6, 2, 7)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 7, 6, 4),Eigen::Matrix<int, 4, 1>(3, 7, 6, 5),Eigen::Matrix<int, 4, 1>(0, 1, 4, 7),Eigen::Matrix<int, 4, 1>(1, 2, 4, 7),Eigen::Matrix<int, 4, 1>(0, 7, 5, 1),Eigen::Matrix<int, 4, 1>(0, 3, 5, 7),Eigen::Matrix<int, 4, 1>(0, 4, 3, 7),Eigen::Matrix<int, 4, 1>(4, 6, 3, 7),Eigen::Matrix<int, 4, 1>(1, 5, 2, 7),Eigen::Matrix<int, 4, 1>(5, 6, 2, 7)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(2, 4, 7, 5),Eigen::Matrix<int, 4, 1>(7, 1, 6, 4),Eigen::Matrix<int, 4, 1>(2, 1, 7, 4),Eigen::Matrix<int, 4, 1>(3, 4, 7, 6),Eigen::Matrix<int, 4, 1>(7, 0, 5, 4),Eigen::Matrix<int, 4, 1>(3, 0, 7, 4)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 4, 7, 5),Eigen::Matrix<int, 4, 1>(7, 1, 6, 4),Eigen::Matrix<int, 4, 1>(2, 1, 7, 4),Eigen::Matrix<int, 4, 1>(3, 4, 7, 6),Eigen::Matrix<int, 4, 1>(7, 3, 5, 4),Eigen::Matrix<int, 4, 1>(3, 0, 5, 4)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 4, 7, 5),Eigen::Matrix<int, 4, 1>(7, 1, 6, 4),Eigen::Matrix<int, 4, 1>(2, 1, 7, 4),Eigen::Matrix<int, 4, 1>(0, 4, 7, 6),Eigen::Matrix<int, 4, 1>(7, 0, 5, 4),Eigen::Matrix<int, 4, 1>(3, 0, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 4, 7, 5),Eigen::Matrix<int, 4, 1>(7, 1, 6, 4),Eigen::Matrix<int, 4, 1>(2, 1, 7, 4),Eigen::Matrix<int, 4, 1>(0, 8, 5, 3),Eigen::Matrix<int, 4, 1>(5, 7, 3, 8),Eigen::Matrix<int, 4, 1>(6, 7, 4, 8),Eigen::Matrix<int, 4, 1>(5, 8, 4, 7),Eigen::Matrix<int, 4, 1>(0, 6, 4, 8),Eigen::Matrix<int, 4, 1>(3, 6, 0, 8),Eigen::Matrix<int, 4, 1>(0, 4, 5, 8),Eigen::Matrix<int, 4, 1>(3, 8, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 4, 7, 5),Eigen::Matrix<int, 4, 1>(7, 2, 6, 4),Eigen::Matrix<int, 4, 1>(2, 1, 6, 4),Eigen::Matrix<int, 4, 1>(3, 4, 7, 6),Eigen::Matrix<int, 4, 1>(7, 0, 5, 4),Eigen::Matrix<int, 4, 1>(3, 0, 7, 4)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 4, 7, 5),Eigen::Matrix<int, 4, 1>(7, 2, 6, 4),Eigen::Matrix<int, 4, 1>(2, 1, 6, 4),Eigen::Matrix<int, 4, 1>(3, 4, 7, 6),Eigen::Matrix<int, 4, 1>(7, 3, 5, 4),Eigen::Matrix<int, 4, 1>(3, 0, 5, 4)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 4, 7, 5),Eigen::Matrix<int, 4, 1>(7, 2, 6, 4),Eigen::Matrix<int, 4, 1>(2, 1, 6, 4),Eigen::Matrix<int, 4, 1>(0, 4, 7, 6),Eigen::Matrix<int, 4, 1>(7, 0, 5, 4),Eigen::Matrix<int, 4, 1>(3, 0, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 4, 7, 5),Eigen::Matrix<int, 4, 1>(7, 2, 6, 4),Eigen::Matrix<int, 4, 1>(2, 1, 6, 4),Eigen::Matrix<int, 4, 1>(0, 8, 5, 3),Eigen::Matrix<int, 4, 1>(5, 7, 3, 8),Eigen::Matrix<int, 4, 1>(6, 7, 4, 8),Eigen::Matrix<int, 4, 1>(5, 8, 4, 7),Eigen::Matrix<int, 4, 1>(0, 6, 4, 8),Eigen::Matrix<int, 4, 1>(3, 6, 0, 8),Eigen::Matrix<int, 4, 1>(0, 4, 5, 8),Eigen::Matrix<int, 4, 1>(3, 8, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 4, 7, 5),Eigen::Matrix<int, 4, 1>(7, 1, 6, 4),Eigen::Matrix<int, 4, 1>(2, 1, 7, 5),Eigen::Matrix<int, 4, 1>(3, 4, 7, 6),Eigen::Matrix<int, 4, 1>(7, 0, 5, 4),Eigen::Matrix<int, 4, 1>(3, 0, 7, 4)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 4, 7, 5),Eigen::Matrix<int, 4, 1>(7, 1, 6, 4),Eigen::Matrix<int, 4, 1>(2, 1, 7, 5),Eigen::Matrix<int, 4, 1>(3, 4, 7, 6),Eigen::Matrix<int, 4, 1>(7, 3, 5, 4),Eigen::Matrix<int, 4, 1>(3, 0, 5, 4)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 4, 7, 5),Eigen::Matrix<int, 4, 1>(7, 1, 6, 4),Eigen::Matrix<int, 4, 1>(2, 1, 7, 5),Eigen::Matrix<int, 4, 1>(0, 4, 7, 6),Eigen::Matrix<int, 4, 1>(7, 0, 5, 4),Eigen::Matrix<int, 4, 1>(3, 0, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 4, 7, 5),Eigen::Matrix<int, 4, 1>(7, 1, 6, 4),Eigen::Matrix<int, 4, 1>(2, 1, 7, 5),Eigen::Matrix<int, 4, 1>(0, 8, 5, 3),Eigen::Matrix<int, 4, 1>(5, 7, 3, 8),Eigen::Matrix<int, 4, 1>(6, 7, 4, 8),Eigen::Matrix<int, 4, 1>(5, 8, 4, 7),Eigen::Matrix<int, 4, 1>(0, 6, 4, 8),Eigen::Matrix<int, 4, 1>(3, 6, 0, 8),Eigen::Matrix<int, 4, 1>(0, 4, 5, 8),Eigen::Matrix<int, 4, 1>(3, 8, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 4, 6, 5),Eigen::Matrix<int, 4, 1>(7, 2, 6, 5),Eigen::Matrix<int, 4, 1>(2, 1, 6, 4),Eigen::Matrix<int, 4, 1>(3, 4, 5, 6),Eigen::Matrix<int, 4, 1>(7, 3, 5, 6),Eigen::Matrix<int, 4, 1>(3, 0, 5, 4)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 4, 6, 5),Eigen::Matrix<int, 4, 1>(7, 2, 6, 5),Eigen::Matrix<int, 4, 1>(2, 1, 6, 4),Eigen::Matrix<int, 4, 1>(0, 4, 5, 6),Eigen::Matrix<int, 4, 1>(0, 5, 7, 6),Eigen::Matrix<int, 4, 1>(3, 0, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 4, 6, 5),Eigen::Matrix<int, 4, 1>(7, 2, 6, 5),Eigen::Matrix<int, 4, 1>(2, 1, 6, 4),Eigen::Matrix<int, 4, 1>(0, 4, 5, 6),Eigen::Matrix<int, 4, 1>(3, 5, 7, 6),Eigen::Matrix<int, 4, 1>(3, 0, 5, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(2, 4, 6, 5),Eigen::Matrix<int, 4, 1>(7, 2, 6, 5),Eigen::Matrix<int, 4, 1>(2, 1, 6, 4),Eigen::Matrix<int, 4, 1>(0, 8, 7, 3),Eigen::Matrix<int, 4, 1>(0, 8, 5, 7),Eigen::Matrix<int, 4, 1>(6, 7, 5, 8),Eigen::Matrix<int, 4, 1>(5, 8, 4, 6),Eigen::Matrix<int, 4, 1>(3, 4, 0, 8),Eigen::Matrix<int, 4, 1>(3, 6, 4, 8),Eigen::Matrix<int, 4, 1>(0, 4, 5, 8),Eigen::Matrix<int, 4, 1>(3, 8, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 4, 6, 5),Eigen::Matrix<int, 4, 1>(1, 6, 7, 5),Eigen::Matrix<int, 4, 1>(2, 1, 7, 5),Eigen::Matrix<int, 4, 1>(3, 4, 5, 6),Eigen::Matrix<int, 4, 1>(7, 3, 5, 6),Eigen::Matrix<int, 4, 1>(3, 0, 5, 4)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 4, 6, 5),Eigen::Matrix<int, 4, 1>(1, 6, 7, 5),Eigen::Matrix<int, 4, 1>(2, 1, 7, 5),Eigen::Matrix<int, 4, 1>(0, 4, 5, 6),Eigen::Matrix<int, 4, 1>(0, 5, 7, 6),Eigen::Matrix<int, 4, 1>(3, 0, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 4, 6, 5),Eigen::Matrix<int, 4, 1>(1, 6, 7, 5),Eigen::Matrix<int, 4, 1>(2, 1, 7, 5),Eigen::Matrix<int, 4, 1>(0, 4, 5, 6),Eigen::Matrix<int, 4, 1>(3, 5, 7, 6),Eigen::Matrix<int, 4, 1>(3, 0, 5, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 4, 6, 5),Eigen::Matrix<int, 4, 1>(1, 6, 7, 5),Eigen::Matrix<int, 4, 1>(2, 1, 7, 5),Eigen::Matrix<int, 4, 1>(0, 8, 7, 3),Eigen::Matrix<int, 4, 1>(0, 8, 5, 7),Eigen::Matrix<int, 4, 1>(6, 7, 5, 8),Eigen::Matrix<int, 4, 1>(5, 8, 4, 6),Eigen::Matrix<int, 4, 1>(3, 4, 0, 8),Eigen::Matrix<int, 4, 1>(3, 6, 4, 8),Eigen::Matrix<int, 4, 1>(0, 4, 5, 8),Eigen::Matrix<int, 4, 1>(3, 8, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 4, 6, 5),Eigen::Matrix<int, 4, 1>(2, 6, 7, 5),Eigen::Matrix<int, 4, 1>(2, 1, 6, 5),Eigen::Matrix<int, 4, 1>(3, 4, 5, 6),Eigen::Matrix<int, 4, 1>(7, 3, 5, 6),Eigen::Matrix<int, 4, 1>(3, 0, 5, 4)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 4, 6, 5),Eigen::Matrix<int, 4, 1>(2, 6, 7, 5),Eigen::Matrix<int, 4, 1>(2, 1, 6, 5),Eigen::Matrix<int, 4, 1>(0, 4, 5, 6),Eigen::Matrix<int, 4, 1>(0, 5, 7, 6),Eigen::Matrix<int, 4, 1>(3, 0, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 4, 6, 5),Eigen::Matrix<int, 4, 1>(2, 6, 7, 5),Eigen::Matrix<int, 4, 1>(2, 1, 6, 5),Eigen::Matrix<int, 4, 1>(0, 4, 5, 6),Eigen::Matrix<int, 4, 1>(3, 5, 7, 6),Eigen::Matrix<int, 4, 1>(3, 0, 5, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 4, 6, 5),Eigen::Matrix<int, 4, 1>(2, 6, 7, 5),Eigen::Matrix<int, 4, 1>(2, 1, 6, 5),Eigen::Matrix<int, 4, 1>(0, 8, 7, 3),Eigen::Matrix<int, 4, 1>(0, 8, 5, 7),Eigen::Matrix<int, 4, 1>(6, 7, 5, 8),Eigen::Matrix<int, 4, 1>(5, 8, 4, 6),Eigen::Matrix<int, 4, 1>(3, 4, 0, 8),Eigen::Matrix<int, 4, 1>(3, 6, 4, 8),Eigen::Matrix<int, 4, 1>(0, 4, 5, 8),Eigen::Matrix<int, 4, 1>(3, 8, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 8, 6, 2),Eigen::Matrix<int, 4, 1>(6, 7, 2, 8),Eigen::Matrix<int, 4, 1>(5, 7, 4, 8),Eigen::Matrix<int, 4, 1>(6, 8, 4, 7),Eigen::Matrix<int, 4, 1>(1, 5, 4, 8),Eigen::Matrix<int, 4, 1>(2, 5, 1, 8),Eigen::Matrix<int, 4, 1>(1, 4, 6, 8),Eigen::Matrix<int, 4, 1>(2, 8, 7, 5),Eigen::Matrix<int, 4, 1>(3, 4, 7, 6),Eigen::Matrix<int, 4, 1>(7, 0, 5, 4),Eigen::Matrix<int, 4, 1>(3, 0, 7, 4)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 8, 6, 2),Eigen::Matrix<int, 4, 1>(6, 7, 2, 8),Eigen::Matrix<int, 4, 1>(5, 7, 4, 8),Eigen::Matrix<int, 4, 1>(6, 8, 4, 7),Eigen::Matrix<int, 4, 1>(1, 5, 4, 8),Eigen::Matrix<int, 4, 1>(2, 5, 1, 8),Eigen::Matrix<int, 4, 1>(1, 4, 6, 8),Eigen::Matrix<int, 4, 1>(2, 8, 7, 5),Eigen::Matrix<int, 4, 1>(3, 4, 7, 6),Eigen::Matrix<int, 4, 1>(7, 3, 5, 4),Eigen::Matrix<int, 4, 1>(3, 0, 5, 4)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 8, 6, 2),Eigen::Matrix<int, 4, 1>(6, 7, 2, 8),Eigen::Matrix<int, 4, 1>(5, 7, 4, 8),Eigen::Matrix<int, 4, 1>(6, 8, 4, 7),Eigen::Matrix<int, 4, 1>(1, 5, 4, 8),Eigen::Matrix<int, 4, 1>(2, 5, 1, 8),Eigen::Matrix<int, 4, 1>(1, 4, 6, 8),Eigen::Matrix<int, 4, 1>(2, 8, 7, 5),Eigen::Matrix<int, 4, 1>(0, 4, 7, 6),Eigen::Matrix<int, 4, 1>(7, 0, 5, 4),Eigen::Matrix<int, 4, 1>(3, 0, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 8, 6, 2),Eigen::Matrix<int, 4, 1>(6, 7, 2, 8),Eigen::Matrix<int, 4, 1>(5, 7, 4, 8),Eigen::Matrix<int, 4, 1>(6, 8, 4, 7),Eigen::Matrix<int, 4, 1>(1, 5, 4, 8),Eigen::Matrix<int, 4, 1>(2, 5, 1, 8),Eigen::Matrix<int, 4, 1>(1, 4, 6, 8),Eigen::Matrix<int, 4, 1>(2, 8, 7, 5),Eigen::Matrix<int, 4, 1>(0, 9, 5, 3),Eigen::Matrix<int, 4, 1>(5, 7, 3, 9),Eigen::Matrix<int, 4, 1>(6, 7, 4, 9),Eigen::Matrix<int, 4, 1>(5, 9, 4, 7),Eigen::Matrix<int, 4, 1>(0, 6, 4, 9),Eigen::Matrix<int, 4, 1>(3, 6, 0, 9),Eigen::Matrix<int, 4, 1>(0, 4, 5, 9),Eigen::Matrix<int, 4, 1>(3, 9, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 8, 7, 2),Eigen::Matrix<int, 4, 1>(1, 8, 6, 7),Eigen::Matrix<int, 4, 1>(5, 7, 6, 8),Eigen::Matrix<int, 4, 1>(6, 8, 4, 5),Eigen::Matrix<int, 4, 1>(2, 4, 1, 8),Eigen::Matrix<int, 4, 1>(2, 5, 4, 8),Eigen::Matrix<int, 4, 1>(1, 4, 6, 8),Eigen::Matrix<int, 4, 1>(2, 8, 7, 5),Eigen::Matrix<int, 4, 1>(3, 4, 5, 6),Eigen::Matrix<int, 4, 1>(7, 3, 5, 6),Eigen::Matrix<int, 4, 1>(3, 0, 5, 4)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 8, 7, 2),Eigen::Matrix<int, 4, 1>(1, 8, 6, 7),Eigen::Matrix<int, 4, 1>(5, 7, 6, 8),Eigen::Matrix<int, 4, 1>(6, 8, 4, 5),Eigen::Matrix<int, 4, 1>(2, 4, 1, 8),Eigen::Matrix<int, 4, 1>(2, 5, 4, 8),Eigen::Matrix<int, 4, 1>(1, 4, 6, 8),Eigen::Matrix<int, 4, 1>(2, 8, 7, 5),Eigen::Matrix<int, 4, 1>(0, 4, 5, 6),Eigen::Matrix<int, 4, 1>(0, 5, 7, 6),Eigen::Matrix<int, 4, 1>(3, 0, 7, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 8, 7, 2),Eigen::Matrix<int, 4, 1>(1, 8, 6, 7),Eigen::Matrix<int, 4, 1>(5, 7, 6, 8),Eigen::Matrix<int, 4, 1>(6, 8, 4, 5),Eigen::Matrix<int, 4, 1>(2, 4, 1, 8),Eigen::Matrix<int, 4, 1>(2, 5, 4, 8),Eigen::Matrix<int, 4, 1>(1, 4, 6, 8),Eigen::Matrix<int, 4, 1>(2, 8, 7, 5),Eigen::Matrix<int, 4, 1>(0, 4, 5, 6),Eigen::Matrix<int, 4, 1>(3, 5, 7, 6),Eigen::Matrix<int, 4, 1>(3, 0, 5, 6)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 8, 7, 2),Eigen::Matrix<int, 4, 1>(1, 8, 6, 7),Eigen::Matrix<int, 4, 1>(5, 7, 6, 8),Eigen::Matrix<int, 4, 1>(6, 8, 4, 5),Eigen::Matrix<int, 4, 1>(2, 4, 1, 8),Eigen::Matrix<int, 4, 1>(2, 5, 4, 8),Eigen::Matrix<int, 4, 1>(1, 4, 6, 8),Eigen::Matrix<int, 4, 1>(2, 8, 7, 5),Eigen::Matrix<int, 4, 1>(0, 9, 7, 3),Eigen::Matrix<int, 4, 1>(0, 9, 5, 7),Eigen::Matrix<int, 4, 1>(6, 7, 5, 9),Eigen::Matrix<int, 4, 1>(5, 9, 4, 6),Eigen::Matrix<int, 4, 1>(3, 4, 0, 9),Eigen::Matrix<int, 4, 1>(3, 6, 4, 9),Eigen::Matrix<int, 4, 1>(0, 4, 5, 9),Eigen::Matrix<int, 4, 1>(3, 9, 7, 6)
				}
			},
			{},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 3),Eigen::Matrix<int, 4, 1>(4, 1, 6, 5),Eigen::Matrix<int, 4, 1>(6, 0, 2, 1),Eigen::Matrix<int, 4, 1>(4, 0, 6, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 3),Eigen::Matrix<int, 4, 1>(4, 1, 6, 5),Eigen::Matrix<int, 4, 1>(6, 4, 2, 1),Eigen::Matrix<int, 4, 1>(4, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 3),Eigen::Matrix<int, 4, 1>(0, 1, 6, 5),Eigen::Matrix<int, 4, 1>(6, 0, 2, 1),Eigen::Matrix<int, 4, 1>(4, 0, 6, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 3),Eigen::Matrix<int, 4, 1>(4, 1, 2, 5),Eigen::Matrix<int, 4, 1>(6, 4, 2, 5),Eigen::Matrix<int, 4, 1>(4, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 3),Eigen::Matrix<int, 4, 1>(0, 1, 2, 5),Eigen::Matrix<int, 4, 1>(0, 2, 6, 5),Eigen::Matrix<int, 4, 1>(4, 0, 6, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 3),Eigen::Matrix<int, 4, 1>(0, 1, 2, 5),Eigen::Matrix<int, 4, 1>(4, 2, 6, 5),Eigen::Matrix<int, 4, 1>(4, 0, 2, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 3),Eigen::Matrix<int, 4, 1>(0, 7, 2, 4),Eigen::Matrix<int, 4, 1>(2, 6, 4, 7),Eigen::Matrix<int, 4, 1>(5, 6, 1, 7),Eigen::Matrix<int, 4, 1>(2, 7, 1, 6),Eigen::Matrix<int, 4, 1>(0, 5, 1, 7),Eigen::Matrix<int, 4, 1>(4, 5, 0, 7),Eigen::Matrix<int, 4, 1>(0, 1, 2, 7),Eigen::Matrix<int, 4, 1>(4, 7, 6, 5)
				},
				{
					Eigen::Matrix<int, 4, 1>(4, 5, 6, 3),Eigen::Matrix<int, 4, 1>(0, 7, 6, 4),Eigen::Matrix<int, 4, 1>(0, 7, 2, 6),Eigen::Matrix<int, 4, 1>(5, 6, 2, 7),Eigen::Matrix<int, 4, 1>(2, 7, 1, 5),Eigen::Matrix<int, 4, 1>(4, 1, 0, 7),Eigen::Matrix<int, 4, 1>(4, 5, 1, 7),Eigen::Matrix<int, 4, 1>(0, 1, 2, 7),Eigen::Matrix<int, 4, 1>(4, 7, 6, 5)
				}
			},
			{},
			{},
			{},
			{},
			{},
			{},
			{}
		}};

		assert(!table[idx].empty());
		return table[idx];
	}


	const std::vector<std::vector<Eigen::Matrix<int, 2, 1>>>& CutTable::get_diag_confs(const int idx) {
		static const std::array<std::vector<std::vector<Eigen::Matrix<int, 2, 1>>>, 64> table= {{

			{

				{
				
				}
			},
			{

				{
				
				}
			},
			{

				{
				
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 4)
				}
			},
			{

				{
				
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(1, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 4)
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(1, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4)
				}
			},
			{},
			{

				{
				
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(1, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{

				{
				
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(1, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(2, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				}
			},
			{},
			{

				{
				
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(2, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				}
			},
			{

				{
				
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(2, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 4)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(2, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 7),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(2, 7)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 7),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(2, 7)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(2, 7),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(2, 7)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(2, 7)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 7),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(2, 7)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(2, 7),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(2, 7)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 7),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 7),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(2, 7)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(2, 7)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 7),Eigen::Matrix<int, 2, 1>(3, 5)
				}
			},
			{},
			{

				{
				
				}
			},
			{

				{
				
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(1, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(1, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(2, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 4),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(3, 5)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 4)
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 5),Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 5),Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(1, 7)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(1, 7)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{},
			{},
			{},
			{},
			{

				{
					Eigen::Matrix<int, 2, 1>(1, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 4)
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(1, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(2, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 5),Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{},
			{},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(2, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 5),Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 7)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 7)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 5),Eigen::Matrix<int, 2, 1>(2, 6),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 7),Eigen::Matrix<int, 2, 1>(1, 7),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(3, 4)
				}
			},
			{},
			{},
			{

				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 4),Eigen::Matrix<int, 2, 1>(1, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 4),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 6)
				},
				{
					Eigen::Matrix<int, 2, 1>(1, 4),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(2, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(2, 4),Eigen::Matrix<int, 2, 1>(2, 5)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 5),Eigen::Matrix<int, 2, 1>(1, 6),Eigen::Matrix<int, 2, 1>(2, 4)
				},
				{
					Eigen::Matrix<int, 2, 1>(0, 6),Eigen::Matrix<int, 2, 1>(1, 4),Eigen::Matrix<int, 2, 1>(2, 5)
				}
			},
			{},
			{},
			{},
			{},
			{},
			{},
			{}
		}};

		//assert(!table[idx].empty());
		return table[idx];
	}


	const std::vector<std::vector<std::array<bool, 4>>>& CutTable::get_surface_conf(const int idx) {
		static const std::array<std::vector<std::vector<std::array<bool, 4>>>, 64> table= {{

			{

				{
					{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}}
				}
			},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}}
				},
				{
					{{false, true, false, false}},{{false, false, false, false}},{{false, false, true, false}},{{false, false, false, false}}
				},
				{
					{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}}
				},
				{
					{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}}
				},
				{
					{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, true, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, true, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				}
			},
			{},
			{

				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{},
			{},
			{},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, true, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				}
			},
			{},
			{},
			{},
			{},
			{},
			{},
			{}
		}};

		assert(!table[idx].empty());
		return table[idx];
	}


	const std::vector<std::vector<Eigen::Matrix<int, 4, 1>>>& CutTable::get_face_id_conf(const int idx) {
		static const std::array<std::vector<std::vector<Eigen::Matrix<int, 4, 1>>>, 64> table= {{

			{

				{
					Eigen::Matrix<int, 4, 1>(0, 1, 2, 3)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(1, 3, 2, -1),Eigen::Matrix<int, 4, 1>(0, -1, 2, 3)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(2, 3, 0, -1),Eigen::Matrix<int, 4, 1>(1, -1, 0, 3)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 3, 2, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(1, 0, 3, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 3, 2, 0),Eigen::Matrix<int, 4, 1>(0, -1, 3, -1),Eigen::Matrix<int, 4, 1>(1, -1, 3, 2)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(0, 3, 1, -1),Eigen::Matrix<int, 4, 1>(2, -1, 1, 3)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 1, 2, 3),Eigen::Matrix<int, 4, 1>(-1, 2, 3, -1),Eigen::Matrix<int, 4, 1>(0, -1, 3, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 1, 2, 3),Eigen::Matrix<int, 4, 1>(1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(0, 2, 3, -1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 3, 0, 1),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(2, 1, 3, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 3, 0, 1),Eigen::Matrix<int, 4, 1>(1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(2, -1, 3, 0)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(3, 2, 1, -1),Eigen::Matrix<int, 4, 1>(0, -1, 1, 2)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 2, 3, 1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 3),Eigen::Matrix<int, 4, 1>(0, 1, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, 3, 1),Eigen::Matrix<int, 4, 1>(1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(0, -1, 2, 3)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(2, 3, -1, -1),Eigen::Matrix<int, 4, 1>(1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(2, -1, -1, 0),Eigen::Matrix<int, 4, 1>(1, 0, -1, -1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 2, 3, 1),Eigen::Matrix<int, 4, 1>(-1, 3, 1, -1),Eigen::Matrix<int, 4, 1>(0, -1, 1, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, 3, 1),Eigen::Matrix<int, 4, 1>(2, -1, 1, -1),Eigen::Matrix<int, 4, 1>(0, 3, 1, -1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(0, 2, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 2, 3, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, 1),Eigen::Matrix<int, 4, 1>(0, 2, 3, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(0, 2, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 2, 3, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, 1),Eigen::Matrix<int, 4, 1>(0, 2, 3, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, 3, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, 3, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(3, 0, 2, -1),Eigen::Matrix<int, 4, 1>(1, -1, 2, 0)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 0, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 3, 2, -1),Eigen::Matrix<int, 4, 1>(1, -1, 2, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, 3, 2),Eigen::Matrix<int, 4, 1>(0, -1, 2, -1),Eigen::Matrix<int, 4, 1>(1, 3, 2, -1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 0, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 0, 3),Eigen::Matrix<int, 4, 1>(1, 2, 0, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, 3, 2),Eigen::Matrix<int, 4, 1>(2, -1, 0, -1),Eigen::Matrix<int, 4, 1>(1, -1, 0, 3)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(0, 2, 3, -1),Eigen::Matrix<int, 4, 1>(0, -1, 3, -1),Eigen::Matrix<int, 4, 1>(1, 0, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 2, 3, -1),Eigen::Matrix<int, 4, 1>(0, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, 2),Eigen::Matrix<int, 4, 1>(1, -1, 3, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 2, 3, -1),Eigen::Matrix<int, 4, 1>(0, -1, 3, -1),Eigen::Matrix<int, 4, 1>(1, 0, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 2, 3, -1),Eigen::Matrix<int, 4, 1>(0, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, 2),Eigen::Matrix<int, 4, 1>(1, -1, 3, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 2, 3, -1),Eigen::Matrix<int, 4, 1>(0, -1, 3, 1),Eigen::Matrix<int, 4, 1>(0, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 2, 3, -1),Eigen::Matrix<int, 4, 1>(0, -1, 3, 1),Eigen::Matrix<int, 4, 1>(0, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 2, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 2, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(0, 3, -1, -1),Eigen::Matrix<int, 4, 1>(2, -1, -1, 3),Eigen::Matrix<int, 4, 1>(0, -1, -1, 1),Eigen::Matrix<int, 4, 1>(2, 1, -1, -1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 2, 1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1),Eigen::Matrix<int, 4, 1>(3, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, 1, 0),Eigen::Matrix<int, 4, 1>(0, -1, 2, -1),Eigen::Matrix<int, 4, 1>(3, -1, 2, 1)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 2, 1, -1),Eigen::Matrix<int, 4, 1>(3, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 2, 1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(3, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 1, -1),Eigen::Matrix<int, 4, 1>(3, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 2, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 2, 1, -1),Eigen::Matrix<int, 4, 1>(3, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 2, 1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(3, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 1, -1),Eigen::Matrix<int, 4, 1>(3, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 2, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 0, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 2, 1, -1),Eigen::Matrix<int, 4, 1>(3, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 0, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 2, 1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(3, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 0, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 1, -1),Eigen::Matrix<int, 4, 1>(3, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 2, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 0, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 0, 2),Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 0),Eigen::Matrix<int, 4, 1>(3, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 0, 2),Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 2, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 0, 2),Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 3),Eigen::Matrix<int, 4, 1>(-1, 2, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 0, 2),Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 0, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 0),Eigen::Matrix<int, 4, 1>(3, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 0, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 2, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 0, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 3),Eigen::Matrix<int, 4, 1>(-1, 2, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 0, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 0),Eigen::Matrix<int, 4, 1>(3, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 2, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 3),Eigen::Matrix<int, 4, 1>(-1, 2, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, 1, -1),Eigen::Matrix<int, 4, 1>(3, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, 1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(3, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 1, -1),Eigen::Matrix<int, 4, 1>(3, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 2, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 0),Eigen::Matrix<int, 4, 1>(3, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 2, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 3),Eigen::Matrix<int, 4, 1>(-1, 2, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(3, 1, 0, -1),Eigen::Matrix<int, 4, 1>(2, -1, 0, 1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(0, -1, -1, 3),Eigen::Matrix<int, 4, 1>(1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(0, 2, -1, -1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 3, 0, -1),Eigen::Matrix<int, 4, 1>(2, -1, 0, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(1, -1, 0, -1),Eigen::Matrix<int, 4, 1>(2, 3, 0, -1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 3),Eigen::Matrix<int, 4, 1>(2, 0, 1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(0, -1, 1, -1),Eigen::Matrix<int, 4, 1>(2, -1, 1, 3)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(2, 1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, 0),Eigen::Matrix<int, 4, 1>(2, -1, 3, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(2, 1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, 0),Eigen::Matrix<int, 4, 1>(2, -1, 3, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(1, 0, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1)
				}
			},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1),Eigen::Matrix<int, 4, 1>(-1, 2, 1, -1),Eigen::Matrix<int, 4, 1>(3, -1, 1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1),Eigen::Matrix<int, 4, 1>(0, -1, 1, -1),Eigen::Matrix<int, 4, 1>(3, 2, 1, -1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(2, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1),Eigen::Matrix<int, 4, 1>(3, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(2, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(3, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(2, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(3, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(2, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(2, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1),Eigen::Matrix<int, 4, 1>(3, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(2, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(3, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(2, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(3, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(2, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(2, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1),Eigen::Matrix<int, 4, 1>(3, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(2, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(3, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(2, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(3, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(2, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 1, 0, -1),Eigen::Matrix<int, 4, 1>(2, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 1, -1),Eigen::Matrix<int, 4, 1>(3, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 1, 0, -1),Eigen::Matrix<int, 4, 1>(2, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 3, 2, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 1, 0, -1),Eigen::Matrix<int, 4, 1>(2, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 3, 2, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 1, 0, -1),Eigen::Matrix<int, 4, 1>(2, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 1, -1),Eigen::Matrix<int, 4, 1>(3, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 3, 2, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 3, 2, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 1, -1),Eigen::Matrix<int, 4, 1>(3, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 3, 2, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 3, 2, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 2, 3, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1),Eigen::Matrix<int, 4, 1>(3, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(3, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(3, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 1, -1),Eigen::Matrix<int, 4, 1>(3, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 3, 2, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 3, 2, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1)
				}
			},
			{},
			{},
			{},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1),Eigen::Matrix<int, 4, 1>(-1, -1, 0, 2),Eigen::Matrix<int, 4, 1>(3, 1, 0, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1),Eigen::Matrix<int, 4, 1>(1, -1, 0, -1),Eigen::Matrix<int, 4, 1>(3, -1, 0, 2)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 3, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0)
				}
			},
			{},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 2, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0)
				}
			},
			{

				{
					Eigen::Matrix<int, 4, 1>(-1, 1, 3, -1),Eigen::Matrix<int, 4, 1>(2, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 0, 2, -1),Eigen::Matrix<int, 4, 1>(3, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 1, 3, -1),Eigen::Matrix<int, 4, 1>(2, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 0, 2, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(3, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 1, 3, -1),Eigen::Matrix<int, 4, 1>(2, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(3, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 1, 3, -1),Eigen::Matrix<int, 4, 1>(2, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(2, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 0, 2, -1),Eigen::Matrix<int, 4, 1>(3, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(2, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 0, 2, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(3, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(2, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(3, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(2, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(2, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 0, 2, -1),Eigen::Matrix<int, 4, 1>(3, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(2, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 0, 2, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(3, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(2, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(3, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(2, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 0),Eigen::Matrix<int, 4, 1>(2, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 0, 1),Eigen::Matrix<int, 4, 1>(3, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 0),Eigen::Matrix<int, 4, 1>(2, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 0),Eigen::Matrix<int, 4, 1>(2, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 1, 0),Eigen::Matrix<int, 4, 1>(2, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 0, 1),Eigen::Matrix<int, 4, 1>(3, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 0, 1),Eigen::Matrix<int, 4, 1>(3, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, -1, 3, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, 3, 0),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, -1),Eigen::Matrix<int, 4, 1>(3, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(3, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(3, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 0, 1),Eigen::Matrix<int, 4, 1>(3, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, 0, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1)
				}
			},
			{},
			{},
			{

				{
					Eigen::Matrix<int, 4, 1>(0, 1, 2, -1),Eigen::Matrix<int, 4, 1>(0, -1, 2, -1),Eigen::Matrix<int, 4, 1>(3, 0, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 1, 2, -1),Eigen::Matrix<int, 4, 1>(0, -1, 2, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, 1),Eigen::Matrix<int, 4, 1>(3, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 1, 2, -1),Eigen::Matrix<int, 4, 1>(0, -1, 2, -1),Eigen::Matrix<int, 4, 1>(3, 0, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 1, 2, -1),Eigen::Matrix<int, 4, 1>(0, -1, 2, -1),Eigen::Matrix<int, 4, 1>(-1, 0, -1, 1),Eigen::Matrix<int, 4, 1>(3, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 1, 2, -1),Eigen::Matrix<int, 4, 1>(0, -1, 2, 3),Eigen::Matrix<int, 4, 1>(0, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 1, 2, -1),Eigen::Matrix<int, 4, 1>(0, -1, 2, 3),Eigen::Matrix<int, 4, 1>(0, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, 2, 1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 1, 2, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1)
				},
				{
					Eigen::Matrix<int, 4, 1>(0, 1, 2, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, 1, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 0),Eigen::Matrix<int, 4, 1>(-1, 0, -1, -1),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 2),Eigen::Matrix<int, 4, 1>(-1, -1, -1, 3),Eigen::Matrix<int, 4, 1>(-1, -1, -1, -1)
				}
			},
			{},
			{},
			{},
			{},
			{},
			{},
			{}
		}};

		assert(!table[idx].empty());
		return table[idx];
	}


}
