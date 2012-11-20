#ifndef __PARTICLE_HPP
#define __PARTICLE_HPP

#include <map>

#include "types.hpp"
#include "Position3.hpp"
#include "Species.hpp"
#include "Identifier.hpp"


namespace ecell4
{

class Particle
{
public:

    Particle()
    {
        ;
    }

    Particle(
        Species const& sp, Position3 const& pos, Real const& radius,
        Real const& D)
        : species_(sp), position_(pos), radius_(radius), D_(D)
    {
        ;
    }

    Position3& position()
    {
        return position_;
    }

    Position3 const& position() const
    {
        return position_;
    }

    Real& radius()
    {
        return radius_;
    }

    Real const& radius() const
    {
        return radius_;
    }

    Real& D()
    {
        return D_;
    }

    Real const& D() const
    {
        return D_;
    }

    Species& species()
    {
        return species_;
    }

    Species const& species() const
    {
        return species_;
    }

private:

    Position3 position_;
    Real radius_, D_;
    Species species_;
};

struct ParticleID:
        public Identifier<ParticleID, unsigned long long, int>
{
    typedef Identifier<ParticleID, unsigned long long, int> base_type;

    ParticleID(value_type const& value = value_type(0, 0))
        : base_type(value)
    {
        ;
    }
};

} // ecell4

#endif /* __PARTICLE_HPP */
