#include <floattetwild/Predicates.hpp>

//extern "C" floatTetWild::Scalar orient3d(const floatTetWild::Scalar *pa, const floatTetWild::Scalar *pb, const floatTetWild::Scalar *pc, const floatTetWild::Scalar *pd);
//extern "C" floatTetWild::Scalar orient2d(const floatTetWild::Scalar *pa, const floatTetWild::Scalar *pb, const floatTetWild::Scalar *pc);
#include <igl/predicates/predicates.h>

#include <geogram/delaunay/delaunay_3d.h>
namespace floatTetWild {
#define GEO_PREDICATES false
    const int Predicates::ORI_POSITIVE;
    const int Predicates::ORI_ZERO;
    const int Predicates::ORI_NEGATIVE;
    const int Predicates::ORI_UNKNOWN;

    int Predicates::orient_3d(const Eigen::Matrix<double, 3, 1>& p1, const Eigen::Matrix<double, 3, 1>& p2, const Eigen::Matrix<double, 3, 1>& p3, const Eigen::Matrix<double, 3, 1>& p4) {
#if GEO_PREDICATES
        const int result = -GEO::PCK::orient_3d(p1.data(), p2.data(), p3.data(), p4.data());
#else
//		const double result = orient3d(p1.data(), p2.data(), p3.data(), p4.data());
        igl::predicates::exactinit();
        auto res = igl::predicates::orient3d(p1, p2, p3, p4);
        double result;
        if(res == igl::predicates::Orientation::POSITIVE)
            result = 1;
        else if(res == igl::predicates::Orientation::NEGATIVE)
            result = -1;
        else
            result = 0;
#endif

//		if (result > SCALAR_ZERO)
//			return ORI_POSITIVE;
//		else if (result < -SCALAR_ZERO)
//			return ORI_NEGATIVE;
//		else
//			return ORI_ZERO;

        if (result > 0)
            return ORI_POSITIVE;
        else if (result < 0)
            return ORI_NEGATIVE;
        else
            return ORI_ZERO;
    }

    int Predicates::orient_3d_tolerance(const Eigen::Matrix<double, 3, 1>& p1, const Eigen::Matrix<double, 3, 1>& p2, const Eigen::Matrix<double, 3, 1>& p3, const Eigen::Matrix<double, 3, 1>& p) {
#if GEO_PREDICATES
        const int result = -GEO::PCK::orient_3d(p1.data(), p2.data(), p3.data(), p.data());
#else
//		const double result = orient3d(p1.data(), p2.data(), p3.data(), p.data());
        igl::predicates::exactinit();
        auto res = igl::predicates::orient3d(p1, p2, p3, p);
        double result;
        if(res == igl::predicates::Orientation::POSITIVE)
            result = 1;
        else if(res == igl::predicates::Orientation::NEGATIVE)
            result = -1;
        else
            result = 0;
#endif

        if (result == 0)
            return ORI_ZERO;

        Eigen::Matrix<double, 3, 1> n = ((p2 - p3).cross(p1 - p3)).normalized();
        double d = std::abs(n.dot(p - p1));
        if (d <= SCALAR_ZERO)
            return Predicates::ORI_ZERO;

        if (result > 0)
            return ORI_POSITIVE;
        else
            return ORI_NEGATIVE;
    }

    double Predicates::orient_3d_volume(const Eigen::Matrix<double, 3, 1>& p1, const Eigen::Matrix<double, 3, 1>& p2, const Eigen::Matrix<double, 3, 1>& p3, const Eigen::Matrix<double, 3, 1>& p4) {
#if GEO_PREDICATES
        const int ori = -GEO::PCK::orient_3d(p1.data(), p2.data(), p3.data(), p4.data());
#else
//		const double result = orient3d(p1.data(), p2.data(), p3.data(), p4.data());
        igl::predicates::exactinit();
        auto res = igl::predicates::orient3d(p1, p2, p3, p4);
        double ori;
        if(res == igl::predicates::Orientation::POSITIVE)
            ori = 1;
        else if(res == igl::predicates::Orientation::NEGATIVE)
            ori = -1;
        else
            ori = 0;
#endif
        if (ori <= 0)
            return ori;
        else
            return (p1 - p4).dot((p2 - p4).cross(p3 - p4)) / 6;
    }

    int Predicates::orient_2d(const Eigen::Matrix<double, 2, 1>& p1, const Eigen::Matrix<double, 2, 1>& p2, const Eigen::Matrix<double, 2, 1>& p3) {
#if GEO_PREDICATES
        const int result = -GEO::PCK::orient_2d(p1.data(), p2.data(), p3.data());
#else
//		const double result = orient2d(p1.data(), p2.data(), p3.data());
        igl::predicates::exactinit();
        auto res = igl::predicates::orient2d(p1, p2, p3);
        double result;
        if(res == igl::predicates::Orientation::POSITIVE)
            result = 1;
        else if(res == igl::predicates::Orientation::NEGATIVE)
            result = -1;
        else
            result = 0;
#endif
        if (result > 0)
            return ORI_POSITIVE;
        else if (result < 0)
            return ORI_NEGATIVE;
        else
            return ORI_ZERO;
    }

} // namespace floatTetWild
