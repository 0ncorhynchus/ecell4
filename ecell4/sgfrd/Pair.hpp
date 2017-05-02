#ifndef ECELL4_SGFRD_PAIR_DOMAIN
#define ECELL4_SGFRD_PAIR_DOMAIN
#include <ecell4/sgfrd/ShellID.hpp>
#include <ecell4/core/Particle.hpp>

namespace ecell4
{
namespace sgfrd
{

class Pair
{
  public:

    typedef ShellID    shell_id_type;
    typedef Particle   particle_type;
    typedef ParticleID particle_id_type;
    typedef std::pair<ParticleID, Particle>   particle_id_pair;
    typedef boost::array<particle_id_pair, 2> particle_array_type;

  public:
    Pair(): dt_(0.), last_time_(0.){}
    ~Pair(){}

    Pair(const Real dt, const Real last_time, shell_id_type const& sh)
        : dt_(dt), last_time_(last_time), shell_id_(sh)
    {}

    Pair(shell_id_type const& sh,
         particle_id_pair const& p0, particle_id_pair const& p1)
        : shell_id_(sh)
    {
        if(p0.second.D() < p1.second.D())
        {
            new(&particles_[0]) particle_id_pair(p0);
            new(&particles_[1]) particle_id_pair(p1);
        }
        else
        {
            new(&particles_[0]) particle_id_pair(p1);
            new(&particles_[1]) particle_id_pair(p0);
        }
    }

    shell_id_type&       shell_id()       {return shell_id_;}
    shell_id_type const& shell_id() const {return shell_id_;}

    Real& dt()       {return dt_;}
    Real  dt() const {return dt_;}
    Real& last_time()       {return last_time_;}
    Real  last_time() const {return last_time_;}

    particle_array_type&       particles()       {return particles_;}
    particle_array_type const& particles() const {return particles_;}

    std::size_t num_shells()   const {return 1;}
    std::size_t multiplicity() const {return 2;}

  private:

    Real dt_;
    Real last_time_;
    shell_id_type   shell_id_;
    particle_array_type particles_;
};

} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_PAIR_DOMAIN */
