#ifndef BD_PROPAGATOR_HPP
#define BD_PROPAGATOR_HPP

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/range/size.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/scoped_ptr.hpp>
#include "Defs.hpp"
#include "generator.hpp"
#include "exceptions.hpp"
#include "freeFunctions.hpp"
#include "utils.hpp"
#include "Logger.hpp"
#include "unassignable_adapter.hpp"

template<typename Ttraits_>
class BDPropagator
{
private:
    typedef typename Ttraits_::world_type world_type;
    typedef typename world_type::transaction_type transaction_type;
    typedef typename world_type::species_id_type species_id_type;
    typedef typename world_type::position_type position_type;
    typedef typename world_type::sphere_type sphere_type;
    typedef typename world_type::species_type species_type;
    typedef typename world_type::length_type length_type;
    typedef typename world_type::particle_id_type particle_id_type;
    typedef typename world_type::particle_type particle_type;
    typedef typename world_type::particle_id_pair particle_id_pair;
    typedef std::vector<particle_id_type> particle_id_vector_type;
    typedef typename world_type::particle_id_pair_generator particle_id_pair_generator;
    typedef typename world_type::particle_id_pair_list particle_id_pair_list;
    typedef typename Ttraits_::random_number_engine random_number_engine;
    typedef typename Ttraits_::network_rules_type network_rules_type;
    typedef typename network_rules_type::reaction_rules reaction_rules;
    typedef typename network_rules_type::reaction_rule_type reaction_rule_type;
    typedef unassignable_adapter<reaction_rule_type, get_default_impl::std::vector> reaction_rule_list_type;

public:
    typedef boost::iterator_range<typename reaction_rule_list_type::const_iterator> reaction_rules_range;

private:
    static const int max_retry_count = 100;

public:
    template<typename Tstate_, typename Trange_>
    BDPropagator(world_type const& world, Tstate_& state, transaction_type& tx, Trange_ const& particles)
        : world_(world),
          tx_(tx),
          rules_(state.get_network_rules()),
          dt_(state.get_dt()),
          rng_(state.get_rng()),
          queue_(),
          rejected_move_count_(0)
    {
        queue_.reserve(boost::size(particles));
        for (typename boost::range_const_iterator<Trange_>::type
                i(boost::begin(particles)),
                e(boost::end(particles)); i != e; ++i)
        {
            queue_.push_back(*i);
        }
    }

    bool operator()()
    {
        if (queue_.empty())
            return false;

        particle_id_type pid(queue_.back());
        queue_.pop_back();
        particle_id_pair pp(tx_.get_particle(pid));

        LOG_DEBUG(("propagating particle %s", boost::lexical_cast<std::string>(pp.first).c_str()));

        try
        {
            if (attempt_reaction(pp))
                return true;
        }
        catch (propagation_error const& reason)
        {
            log_.info("first-order reaction rejected (reason: %s)", reason.what());
            ++rejected_move_count_;
            return true;
        }

        const species_type species(tx_.get_species(pp.second.sid()));
        if (species.D() == 0.)
            return true;

        position_type new_pos(
            world_.apply_boundary(
                add(pp.second.position(), drawR_free(species.D()))));
        particle_id_pair particle_to_update(
                pp.first, particle_type(species.id(),
                sphere_type(new_pos, species.radius())));
        boost::scoped_ptr<particle_id_pair_list> overlapped(
            tx_.check_overlap(particle_to_update));
        if (overlapped)
        {
            switch (overlapped->size())
            {
            case 1:
                {
                    particle_id_pair const& closest(overlapped->at(0));
                    try
                    {
                        if (attempt_reaction(pp, closest))
                            return true;
                    }
                    catch (propagation_error const& reason)
                    {
                        log_.info("second-order reaction rejected (reason: %s)", reason.what());
                        ++rejected_move_count_;
                        return true;
                    }
                }
                break;

            default:
                log_.info("collision involving two or more particles; move rejected");
                ++rejected_move_count_;
                return true;
            }
        }
        tx_.update_particle(particle_to_update);
        return true;
    }

    reaction_rules_range get_reactions()
    {
        return reaction_rules_range(reactions_occurred_);
    }

    std::size_t get_rejected_move_count()
    {
        return rejected_move_count_;
    }

private:
    position_type drawR_free(Real D)
    {
        boost::normal_distribution<Real> dist(0.0, std::sqrt(2.0 * D * dt_));
        return position_type(dist(ur_), dist(ur_), dist(ur_));
    }

    Real getP_acct(Real k, Real D, Real sigma)
    {
        const Real p(k * dt_ / (I_bd(sigma, dt_, D) * 4.0 * M_PI));
        BOOST_ASSERT(p >= 0.);
        if (p >= 1.0)
        {
            throw propagation_error(
                "invalid acceptance ratio ("
                + boost::lexical_cast<std::string>(p)
                + ") for reaction rate "
                + boost::lexical_cast<std::string>(k)
                + ".");
        }
        return p;
    }

    bool attempt_reaction(particle_id_pair const& pp)
    {
        reaction_rule_type r(determine_reaction(pp.second.sid()));
        if (::valid(r))
        {
            typename reaction_rule_type::species_type_id_range products(
                    r.get_products());
            switch (boost::size(products))
            {
            case 0:
                remove_particle(pp.first);
                break;

            case 1:
                {
                    const species_type s0(tx_.get_species(products[0]));
                    const particle_id_pair new_p(
                        pp.first, particle_type(products[0],
                            sphere_type(pp.second.position(), s0.radius())));
                    if (boost::scoped_ptr<particle_id_pair_list>(tx_.check_overlap(new_p)))
                    {
                        throw propagation_error("no space");
                    }

                    tx_.update_particle(new_p);
                }
                break;

            case 2:
                {
                    const species_type s0(tx_.get_species(products[0])),
                            s1(tx_.get_species(products[1]));
                    const Real D01(s0.D() + s1.D());
                    const length_type r01(s0.radius() + s1.radius());
                    int i = max_retry_count;

                    for (;;)
                    {
                        if (--i < 0)
                        {
                            throw propagation_error("no space");
                        }

                        const Real rnd(ur_());
                        length_type pair_distance(
                            drawR_gbd(rnd, r01, dt_, D01));
                        const position_type m(random_unit_vector() * pair_distance);
                        const position_type np0(
                            world_.apply_boundary(pp.second.position()
                                + m * (s0.D() / D01)));
                        const position_type np1(
                            world_.apply_boundary(pp.second.position()
                                + m * (s1.D() / D01)));
                        boost::scoped_ptr<particle_id_pair_list> overlapped_s0(
                            tx_.check_overlap(
                                sphere_type(np0, s0.radius()),
                                pp.first));
                        boost::scoped_ptr<particle_id_pair_list> overlapped_s1(
                            tx_.check_overlap(
                                sphere_type(np1, s1.radius()),
                                pp.first));
                        if (!overlapped_s0 && !overlapped_s1)
                            break;
                    }
                }
                break;
            default:
                throw not_implemented("reactions that produces more than three products are not supported");
            }
            reactions_occurred_.push_back(r);
            return true;
        }
        return false;
    }

    bool attempt_reaction(particle_id_pair const& pp0, particle_id_pair const& pp1)
    {
        reaction_rule_type r(determine_reaction(pp0.second.sid(), pp1.second.sid()));
        if (::valid(r))
        {
            const species_type s0(tx_.get_species(pp0.second.sid())),
                    s1(tx_.get_species(pp1.second.sid()));
            const Real D01(s0.D() + s1.D());
            const length_type r01(s0.radius() + s1.radius());
            const Real p(getP_acct(r.k(), D01, r01));
            const Real rnd(ur_());
            if (p < rnd)
            {
                return false;
            }

            LOG_DEBUG(("fire reaction"));
            const typename reaction_rule_type::species_type_id_range products(
                r.get_products());
            BOOST_ASSERT(boost::size(products) == 1);
            const species_id_type product(products[0]);

            const position_type new_pos(
                world_.apply_boundary(
                    divide(
                        add(multiply(pp0.second.position(), s1.D()),
                            multiply(world_.cyclic_transpose(
                                pp1.second.position(),
                                pp0.second.position()), s0.D())),
                        D01)));
            if (boost::scoped_ptr<particle_id_pair_list>(tx_.check_overlap(sphere_type(new_pos, product))))
            {
                throw propagation_error("no space");
            }

            remove_particle(pp0.first);
            remove_particle(pp1.first);
            tx_.new_particle(product, new_pos);

            reactions_occurred_.push_back(r);
            return true;
        }
        return false;
    }

    reaction_rule_type const& determine_reaction(species_id_type const& sid)
    {
        const Real rnd(ur_() / dt_);
        Real prob = 0;

        reaction_rules const& rules(rules_.query_reaction_rule(sid));
        for (typename boost::range_const_iterator<reaction_rules>::type
                i(rules.begin()), e(rules.end()); i != e; ++i)
        {
            reaction_rule_type const& r(*i);
            prob += r.k();
            if (prob >= rnd)
                return r;
        }

        return empty_reaction_rule_;
    }

    reaction_rule_type const& determine_reaction(species_id_type const& sid0, species_id_type const& sid1)
    {
        const Real rnd(ur_() / dt_);
        Real prob = 0;

        reaction_rules const& rules(rules_.query_reaction_rule(sid0, sid1));
        for (typename boost::range_const_iterator<reaction_rules>::type
                i(rules.begin()), e(rules.end()); i != e; ++i)
        {
            reaction_rule_type const& r(*i);
            prob += r.k();
            if (prob >= rnd)
                return r;
        }

        return empty_reaction_rule_;
    }

    void remove_particle(particle_id_type const& pid)
    {
        LOG_DEBUG(("remove particle %s", boost::lexical_cast<std::string>(pid).c_str()));
        tx_.remove_particle(pid);
        typename particle_id_vector_type::iterator i(
            std::find(queue_.begin(), queue_.end(), pid));
        if (queue_.end() != i)
            queue_.erase(i);
    }

private:
    position_type random_unit_vector()
    {
        position_type v(ur_() - .5, ur_() - .5, ur_() - .5);
        return v / length(v);
    }

private:
    world_type const& world_;
    transaction_type& tx_;
    network_rules_type const& rules_;
    Real const dt_;
    random_number_engine& rng_;
    particle_id_vector_type queue_;
    reaction_rule_list_type reactions_occurred_;
    reaction_rule_type empty_reaction_rule_;
    int rejected_move_count_;
    static Logger& log_;
};

template<typename Ttraits_>
Logger& BDPropagator<Ttraits_>::log_(Logger::get_logger("BDPropagator"));

#endif /* BD_PROPAGATOR_HPP */
