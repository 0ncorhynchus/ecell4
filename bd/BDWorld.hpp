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
    const Real radius;
    const Real D;
};

class BDWorld
{
public:

    typedef ParticleInfo particle_info_type;
    typedef ParticleSpace::particle_container_type particle_container_type;

public:

    BDWorld(
        const Position3& edge_lengths,
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
    ParticleID new_particle(const Particle& p)
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
    ParticleInfo get_particle_info(const Species& sp) const
    {
        const Real radius(std::atof(sp.get_attribute("radius").c_str()));
        const Real D(std::atof(sp.get_attribute("D").c_str()));
        ParticleInfo info = {radius, D};
        return info;
    }

    // SpaceTraits

    const Real& t() const
    {
        return (*ps_).t();
    }

    void set_t(const Real& t)
    {
        (*ps_).set_t(t);
    }

    // ParticleSpaceTraits

    const Position3& edge_lengths() const
    {
        return (*ps_).edge_lengths();
    }

    Integer num_particles() const
    {
        return (*ps_).num_particles();
    }

    Integer num_particles(const Species& species) const
    {
        return (*ps_).num_particles(species);
    }

    bool has_particle(const ParticleID& pid) const
    {
        return (*ps_).has_particle(pid);
    }

    std::vector<std::pair<ParticleID, Particle> > list_particles() const
    {
        return (*ps_).list_particles();
    }

    std::vector<std::pair<ParticleID, Particle> >
    list_particles(const Species& species) const
    {
        return (*ps_).list_particles(species);
    }

    // ParticleSpace member functions

    bool update_particle(const ParticleID& pid, const Particle& p)
    {
        return (*ps_).update_particle(pid, p);
    }

    std::pair<ParticleID, Particle>
    get_particle(const ParticleID& pid) const
    {
        return (*ps_).get_particle(pid);
    }

    void remove_particle(const ParticleID& pid)
    {
        (*ps_).remove_particle(pid);
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Position3& pos, const Real& radius) const
    {
        return (*ps_).list_particles_within_radius(pos, radius);
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Position3& pos, const Real& radius, const ParticleID& ignore) const
    {
        return (*ps_).list_particles_within_radius(pos, radius, ignore);
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Position3& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
    {
        return (*ps_).list_particles_within_radius(pos, radius, ignore1, ignore2);
    }

    inline Position3 periodic_transpose(
        const Position3& pos1, const Position3& pos2) const
    {
        return (*ps_).periodic_transpose(pos1, pos2);
    }

    inline Position3 apply_boundary(const Position3& pos) const
    {
        return (*ps_).apply_boundary(pos);
    }

    inline Real distance_sq(const Position3& pos1, const Position3& pos2) const
    {
        return (*ps_).distance_sq(pos1, pos2);
    }

    inline Real distance(const Position3& pos1, const Position3& pos2) const
    {
        return (*ps_).distance(pos1, pos2);
    }

    // CompartmentSpaceTraits

    Integer num_molecules(const Species& sp) const
    {
        return num_particles(sp);
    }

    void add_molecules(const Species& sp, const Integer& num)
    {
        extras::throw_in_particles(*this, sp, num, *rng());
    }

    // CompartmentSpace member functions

    void add_species(const Species& sp)
    {
        ; // just ignore
    }

    // Optional members

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return rng_;
    }

    const particle_container_type& particles() const
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
