#ifndef GFRD_POLYGON_VECTOR_OPERATION
#define GFRD_POLYGON_VECTOR_OPERATION

#include "Defs.hpp"
#include <ecell4/core/Real3.hpp>
#include <boost/math/quaternion.hpp>

template<typename coordT>
struct value_type_helper
{
    typedef typename coordT::value_type type;
};

template<typename coordT>
coordT rotation(const typename value_type_helper<coordT>::type angle,
                const coordT& axis, const coordT& target)
{
    typedef typename value_type_helper<coordT>::type valueT;
    typedef boost::math::quaternion<valueT> Quaternion;
    const valueT cos_t(cos(angle / 2));
    const valueT sin_t(sin(angle / 2));
    const valueT sin_normalize(sin_t / length(axis));

    const Quaternion Q(cos_t, axis[0] * sin_normalize, 
                              axis[1] * sin_normalize,
                              axis[2] * sin_normalize);
    const Quaternion P(0e0, target[0], target[1], target[2]);
    const Quaternion S(Q * P * conj(Q));

    return coordT(S.R_component_2(), S.R_component_3(), S.R_component_4());
}

template<typename coordT>
typename value_type_helper<coordT>::type
angle(const coordT& lhs, const coordT& rhs)
{
    typedef typename value_type_helper<coordT>::type valueT;
    const valueT lensq_l = length_sq(lhs);
    const valueT lensq_r = length_sq(rhs);
    const valueT inner = dot_product(lhs, rhs);
    return acos(inner / std::sqrt(lensq_l * lensq_r));
}

template<typename coordT>
bool is_same_vec(const coordT& lhs, const coordT& rhs,
         const typename value_type_helper<coordT>::type tol = 1e-10)
{
    return (std::abs(lhs[0] - rhs[0]) < tol) &&
           (std::abs(lhs[1] - rhs[1]) < tol) &&
           (std::abs(lhs[2] - rhs[2]) < tol);
}

template<typename coordT>
coordT reflect_plane(const coordT& begin, const coordT& end,
                     const coordT& normal, const coordT& plane)
{
    typedef typename value_type_helper<coordT>::type valueT;
    const valueT norm_b = dot_product((begin - plane), normal);
    const valueT norm_e = dot_product((end - plane), normal);
    if(norm_b == 0.0)
        throw std::invalid_argument("reflection: begin is on the plane");
    else if(std::abs(norm_e) < 1e-10) //end is just on the plane... add margin(1e-10)
        return (begin * 1e-10) + (end * (1.0 - 1e-10));
    else if(norm_b * norm_e >= 0.0) // same side of plane
        return end;
    else // reflect
        return end - (normal * (norm_e * 2.0));
}
#endif 
