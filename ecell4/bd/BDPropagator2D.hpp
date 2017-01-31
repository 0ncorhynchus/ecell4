#ifndef __ECELL4_BD_BD_PROPAGATOR_2D_HPP
#define __ECELL4_BD_BD_PROPAGATOR_2D_HPP

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>

#include "BDPropagator.hpp"
#include "BDPolygon.hpp"
#include "rotate_vector.hpp"

namespace ecell4
{

namespace bd
{

class BDPropagator2D
{
public:

    typedef ReactionInfo reaction_info_type;
    typedef BDPolygon polygon_type;
    typedef BDPolygon::face_type face_type;
    typedef BDPolygon::face_id_type face_id_type;

public:

    BDPropagator2D(
        Model& model, BDWorld& world, RandomNumberGenerator& rng, const Real& dt,
        std::vector<std::pair<ReactionRule, reaction_info_type> >& last_reactions)
        : model_(model), world_(world), poly_(world.container_2D().polygon()),
          rng_(rng), dt_(dt), last_reactions_(last_reactions), max_retry_count_(1)
    {
        queue_ = world.container_2D().list_particles();
        shuffle(rng_, queue_);
    }

    bool operator()();

    inline Real dt() const
    {
        return dt_;
    }

    inline RandomNumberGenerator& rng()
    {
        return rng_;
    }

    bool attempt_reaction(const ParticleID& pid, const Particle& particle,
                          const face_id_type& fid);
    bool attempt_reaction(
        const ParticleID& pid1, const Particle& particle1, const face_id_type& f1,
        const ParticleID& pid2, const Particle& particle2, const face_id_type& f2);

    void remove_particle(const ParticleID& pid);

    inline Real3 draw_displacement(const Particle& particle, const Real3& normal)
    {
        //TODO do more sophisticated way
        assert(length_sq(normal) - 1.0 < 1e-12);
        const Real sigma = std::sqrt(2 * particle.D() * dt());
        const Real x = rng_.gaussian(sigma);
        const Real y = rng_.gaussian(sigma);
        const Real leng  = std::sqrt(x * x + y * y);

        const Real theta = rng_.random() * 2. * M_PI;
        const Real3 dr = Real3(normal[1], -normal[0], 0) /
            std::sqrt(normal[1] * normal[1] + normal[0] * normal[0]);
        const Real3 disp = rotate(theta, normal, dr);

        assert(std::abs(dot_product(normal, disp)) < 1e-10);
        return disp * leng;
    }

    inline Real3 draw_ipv(const Real& sigma, const Real& t, const Real& D)
    {// XXX!
        return random_ipv_3d(rng(), sigma, t, D);
    }

private:

    class particle_finder
        : public std::unary_function<std::pair<ParticleID, Particle>, bool>
    {
    public:

        particle_finder(const ParticleID& pid)
            : pid_(pid)
        {
            ;
        }

        bool operator()(std::pair<ParticleID, Particle> pid_particle_pair)
        {
            return (pid_particle_pair.first == pid_);
        }

    protected:

        ParticleID pid_;
    };

protected:

    Model& model_;
    BDWorld& world_;
    const BDPolygon& poly_; // XXX additional
    RandomNumberGenerator& rng_;
    Real dt_;
    std::vector<std::pair<ReactionRule, reaction_info_type> >& last_reactions_;
    Integer max_retry_count_;

    BDWorld::particle_container_type queue_;
};

} // bd

} // ecell4

#endif /* __ECELL4_BD_BD_PROPAGATOR_HPP */
