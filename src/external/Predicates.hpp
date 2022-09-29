#pragma once

#include <floattetwild/Types.hpp>

namespace floatTetWild {

	class Predicates
	{
	public:
		static const int ORI_POSITIVE = 1;
        static const int ORI_ZERO = 0;
        static const int ORI_NEGATIVE = -1;
		static const int ORI_UNKNOWN = INT_MAX;

		static int orient_3d(const Eigen::Matrix<double, 3, 1>& p1, const Eigen::Matrix<double, 3, 1>& p2, const Eigen::Matrix<double, 3, 1>& p3, const Eigen::Matrix<double, 3, 1>& p4);
		static int orient_3d_tolerance(const Eigen::Matrix<double, 3, 1>& p1, const Eigen::Matrix<double, 3, 1>& p2, const Eigen::Matrix<double, 3, 1>& p3, const Eigen::Matrix<double, 3, 1>& p4);
		static double orient_3d_volume(const Eigen::Matrix<double, 3, 1>& p1, const Eigen::Matrix<double, 3, 1>& p2, const Eigen::Matrix<double, 3, 1>& p3, const Eigen::Matrix<double, 3, 1>& p4);

        static int orient_2d(const Eigen::Matrix<double, 2, 1>& p1, const Eigen::Matrix<double, 2, 1>& p2, const Eigen::Matrix<double, 2, 1>& p3);
	};

} // namespace floattetwild
