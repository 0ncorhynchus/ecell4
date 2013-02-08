#ifndef __BD_WORLD_HPP
#define __BD_WORLD_HPP

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/extras.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/ParticleSpace.hpp>


namespace ecell4
{

namespace bd
{

struct ParticleInfo
{
    Real const radius;
    Real const D;
};

class BDWorld
{
public:

    typedef ParticleInfo particle_info_type;
    typedef ParticleSpace::particle_container_type particle_container_type;

public:

    BDWorld(
        Position3 const& edge_lengths,
        boost::shared_ptr<RandomNumberGenerator> rng)
        : ps_(new ParticleSpaceVectorImpl(edge_lengths)), rng_(rng)
    {
        ;
    }

    /**
     * create and add a new particle
     * @param p a particle
     * @return pid a particle id
     */
    ParticleID new_particle(Particle const& p)
    {
        ParticleID pid(pidgen_());
        // if (has_particle(pid))
        // {
        //     throw AlreadyExists("particle already exists");
        // }
        (*ps_).update_particle(pid, p);
        return pid;
    }

    /**
     * draw attributes of species and return it as a particle info.
     * @param sp a species
     * @return info a particle info
     */
    ParticleInfo get_particle_info(Species const& sp) const
    {
        const Real radius(std::atof(sp.get_attribute("radius").c_str()));
        const Real D(std::atof(sp.get_attribute("D").c_str()));
        ParticleInfo info = {radius, D};
        return info;
    }

    // SpaceTraits

    Real const& t() const
    {
        return (*ps_).t();
    }

    void set_t(Real const& t)
    {
        (*ps_).set_t(t);
    }

    // ParticleSpaceTraits

    Position3 const& edge_lengths() const
    {
        return (*ps_).edge_lengths();
    }

    Integer num_particles() const
    {
        return (*ps_).num_particles();
    }

    Integer num_particles(Species const& species) const
    {
        return (*ps_).num_particles(species);
    }

    bool has_particle(ParticleID const& pid) const
    {
        return (*ps_).has_particle(pid);
    }

    std::vector<std::pair<ParticleID, Particle> > list_particles() const
    {
        return (*ps_).list_particles();
    }

    std::vector<std::pair<ParticleID, Particle> >
    list_particles(Species const& species) const
    {
        return (*ps_).list_particles(species);
    }

    // ParticleSpace member functions

    bool update_particle(ParticleID const& pid, Particle const& p)
    {
        return (*ps_).update_particle(pid, p);
    }

    std::pair<ParticleID, Particle>
    get_particle(ParticleID const& pid) const
    {
        return (*ps_).get_particle(pid);
    }

    void remove_particle(ParticleID const& pid)
    {
        (*ps_).remove_particle(pid);
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        Position3 const& pos, Real const& radius) const
    {
        return (*ps_).list_particles_within_radius(pos, radius);
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        Position3 const& pos, Real const& radius, ParticleID const& ignore) const
    {
        return (*ps_).list_particles_within_radius(pos, radius, ignore);
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        Position3 const& pos, Real const& radius,
        ParticleID const& ignore1, ParticleID const& ignore2) const
    {
        return (*ps_).list_particles_within_radius(pos, radius, ignore1, ignore2);
    }

    inline Position3 periodic_transpose(
        Position3 const& pos1, Position3 const& pos2) const
    {
        return (*ps_).periodic_transpose(pos1, pos2);
    }

    inline Position3 apply_boundary(Position3 const& pos) const
    {
        return (*ps_).apply_boundary(pos);
    }

    inline Real distance_sq(Position3 const& pos1, Position3 const& pos2) const
    {
        return (*ps_).distance_sq(pos1, pos2);
    }

    inline Real distance(Position3 const& pos1, Position3 const& pos2) const
    {
        return (*ps_).distance(pos1, pos2);
    }

    // CompartmentSpaceTraits

    Integer num_molecules(Species const& sp) const
    {
        return num_particles(sp);
    }

    void add_molecules(Species const& sp, Integer const& num)
    {
        extras::throw_in_particles(*this, sp, num, *rng());
    }

    // Optional members

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return rng_;
    }

    particle_container_type const& particles() const
    {
        return (*ps_).particles();
    }

protected:

    boost::scoped_ptr<ParticleSpace> ps_;
    boost::shared_ptr<RandomNumberGenerator> rng_;

    SerialIDGenerator<ParticleID> pidgen_;
};

} // bd

} // ecell4

#endif /* __BD_WORLD_HPP */
