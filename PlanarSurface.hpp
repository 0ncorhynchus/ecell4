#ifndef PLANAR_SURFACE_HPP
#define PLANAR_SURFACE_HPP

#include "Surface.hpp"
#include "Plane.hpp"

template<typename Ttraits_>
class PlanarSurface
    : public BasicRegionImpl<Ttraits_, Plane<typename Ttraits_::length_type> >
{
public:
    typedef BasicRegionImpl<Ttraits_, Plane<typename Ttraits_::length_type> > base_type;
    typedef typename base_type::traits_type traits_type;
    typedef typename base_type::identifier_type identifier_type;
    typedef typename base_type::shape_type shape_type;
    typedef typename base_type::rng_type rng_type;
    typedef typename base_type::position_type position_type;
    typedef typename base_type::length_type length_type;

    virtual position_type random_position(rng_type& rng) const
    {
        return ::random_position(base_type::shape(), boost::bind(&rng_type::uniform, rng, -1., 1.));
    }

    virtual position_type random_vector(length_type const& r, rng_type& rng) const
    {
        return multiply(
            normalize(
                add(
                    multiply(
                        base_type::shape().units()[0], rng.uniform(-1., 1.)),
                    multiply(
                        base_type::shape().units()[1], rng.uniform(-1., 1.)))), r);
    }

    virtual position_type bd_displacement(length_type const& r, rng_type& rng) const
    {
        length_type const x(rng.normal(0., r)), y(rng.normal(0., r));
        return add(
            multiply(base_type::shape().unit_x(), x),
            multiply(base_type::shape().unit_y(), y));
    }

    virtual length_type minimal_distance(length_type const& radius) const
    {
        return radius + traits_type::MINIMAL_SEPARATION_FACTOR;
    }

    PlanarSurface(identifier_type const& id, shape_type const& shape)
        : base_type(id, shape) {}
};


#endif /* PLANAR_SURFACE_HPP */
