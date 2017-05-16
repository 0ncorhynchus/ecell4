#ifndef ECELL4_SGFRD_SIMULATOR
#define ECELL4_SGFRD_SIMULATOR

#include <greens_functions/GreensFunction2DAbsSym.hpp>
#include <greens_functions/GreensFunction2DRefWedgeAbs.hpp>
#include <ecell4/core/SimulatorBase.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/geometry.hpp>
#include <ecell4/bd/BDPropagator.hpp>
#include "ShellContainer.hpp"
#include "ShellVisitorApplier.hpp"
#include "ShellVisitors.hpp"
#include "SGFRDEvent.hpp"
#include "SGFRDWorld.hpp"
#include <boost/make_shared.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/container/static_vector.hpp>
#include <iostream>

#ifndef SGFRD_NDEBUG
#define DUMP_MESSAGE( str ) std::cerr << str << std::endl
#else
#define DUMP_MESSAGE( str )
#endif //NDEBUG

namespace ecell4
{
namespace sgfrd
{

class SGFRDSimulator :
    public ecell4::SimulatorBase<ecell4::Model, SGFRDWorld>
{
  public:

    // polygon
    typedef polygon_traits polygon_traits_type;
    typedef ecell4::Polygon<polygon_traits_type>     polygon_type;
    typedef polygon_type::triangle_type     triangle_type;
    typedef polygon_type::face_id_type      face_id_type;
    typedef polygon_type::edge_id_type      edge_id_type;
    typedef polygon_type::vertex_id_type    vertex_id_type;
    typedef polygon_type::face_descripter   face_descripter;
    typedef polygon_type::edge_descripter   edge_descripter;
    typedef polygon_type::vertex_descripter vertex_descripter;
    typedef polygon_type::local_index_type  local_index_type;
    typedef polygon_type::barycentric_type  barycentric_type;

    // Event
    typedef SGFRDEvent                      event_type;
    typedef SGFRDEventScheduler             scheduler_type;
    typedef SGFRDEvent::domain_type         domain_type;
    typedef SGFRDEventScheduler::value_type event_id_pair_type;

    // Simulator
    typedef ecell4::SimulatorBase<ecell4::Model, SGFRDWorld> base_type;
    typedef base_type::world_type world_type;
    typedef base_type::model_type model_type;

    // ShellContainer
    typedef ecell4::SerialIDGenerator<ShellID> shell_id_generator_type;
    typedef ShellContainer<polygon_traits_type> shell_container_type;
    typedef shell_container_type::shell_type           shell_type;
    typedef shell_container_type::shell_id_pair_type   shell_id_pair_type;
    typedef shell_container_type::circle_type          circle_type;
    typedef shell_container_type::conical_surface_type conical_surface_type;
    typedef shell_container_type::circular_shell_type  circular_shell_type;
    typedef shell_container_type::conical_surface_shell_type
            conical_surface_shell_type;
    typedef mutable_shell_visitor_applier<polygon_traits_type>
            mutable_shell_visitor_applier_type;

    // reaction stuff
    typedef ecell4::ReactionRule     reaction_rule_type;
    typedef ecell4::bd::MoleculeInfo molecule_info_type;
    typedef ecell4::bd::ReactionInfo reaction_info_type;

  public:

    SGFRDSimulator(const boost::shared_ptr<world_type>& world,
                   const boost::shared_ptr<model_type>& model,
                   Real bd_dt_factor = 1e-5)
        : base_type(model, world), dt_(0), bd_dt_factor_(bd_dt_factor),
          rng_(*(world->rng())), shell_container_(*(world->polygon())),
          mut_sh_vis_applier(shell_container_)
    {}

    SGFRDSimulator(boost::shared_ptr<world_type> world, Real bd_dt_factor = 1e-5)
        : base_type(world), dt_(0), bd_dt_factor_(bd_dt_factor),
          rng_(*(world->rng())), shell_container_(*(world->polygon())),
          mut_sh_vis_applier(shell_container_)
    {}

    ~SGFRDSimulator(){}

    void initialize()
    {
        std::vector<std::pair<ParticleID, Particle> > const& ps =
            this->world_->list_particles();
        for(std::vector<std::pair<ParticleID, Particle> >::const_iterator
            iter = ps.begin(); iter != ps.end(); ++iter)
        {
            add_event(create_minimum_domain(create_minimum_shell(
                      iter->first, iter->second, this->get_face_id(iter->first)),
                  iter->first, iter->second));
        }
        return ;
    }
    void finalize()
    {
        std::vector<domain_type> tmp;
        while(scheduler_.size() != 0)
        {
            burst_event(*(this->scheduler_.pop().second), tmp);
        }
        return ;
    }

    void step()
    {
        this->set_time(this->scheduler_.next_time());
        this->fire_event(this->scheduler_.pop());
        DUMP_MESSAGE("now " << shell_container_.num_shells() << " shells exist.");
        return;
    }
    bool step(const Real& upto)
    {
        this->step();
        return this->time() < upto;
    }

    Real dt() const {return dt_;}

    bool check_reaction() const {return last_reactions_.size() > 0;}
    std::vector<std::pair<ReactionRule, reaction_info_type> > const&
    last_reactions() const {return last_reactions_;}

  private: // wrappers

    Real uniform_real(){return this->rng_.random();}

    world_type   const& world()   const {return *(this->world_);}
    polygon_type const& polygon() const {return *(this->world_->polygon());}

    bool update_particle(const ParticleID& pid, const Particle& p,
                         const face_id_type& fid)
    {return this->world_->update_particle(pid, p, fid);}

    shell_type&       get_shell(ShellID const& id)
    {return shell_container_.get_shell(id);}
    shell_type const& get_shell(ShellID const& id) const
    {return shell_container_.get_shell(id);}
    void remove_shell(ShellID const& id)
    {return shell_container_.remove_shell(id);}

    face_id_type      get_face_id(const ParticleID& pid) const
    {return this->world_->get_face_id(pid);}
    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const
    {return this->world_->get_particle(pid);}

    // consider scheduler.time
    Real time() const {return this->world_->t();}
    void set_time(const Real t) {return this->world_->set_t(t);}

  private:

    // visitors decl: algorithm specific part {{{

    struct domain_firer : boost::static_visitor<void>
    {
        domain_firer(SGFRDSimulator& s): sim(s){}
        void operator()(const Single&);
        void operator()(const Pair&);
        void operator()(const Multi&);
      private:
        SGFRDSimulator& sim;
    };

    struct domain_burster : boost::static_visitor<void>
    {
        domain_burster(SGFRDSimulator& s, std::vector<domain_type>& r)
            : results(r), sim(s)
        {}
        void operator()(const Single& dom)
        {
            single_burster burster(sim, dom, results);
            boost::apply_visitor(burster, sim.get_shell(dom.shell_id()));
            sim.remove_shell(dom.shell_id());
            return;
        }
        void operator()(const Pair& dom)
        {
            pair_burster burster(sim, dom, results);
            boost::apply_visitor(burster, sim.get_shell(dom.shell_id()));
            sim.remove_shell(dom.shell_id());
            return;
        }
        void operator()(const Multi& dom)
        {// TODO
            std::cerr << "bursting multi shell" << std::endl;
            return;
        }

        std::vector<domain_type>& results;

      private:
        SGFRDSimulator& sim;
    };
    struct single_burster : boost::static_visitor<void>
    {
        single_burster(SGFRDSimulator& s, const domain_type& d,
                       std::vector<domain_type>& r)
            : results(r), sim(s), dom(d)
        {}

        void operator()(const circular_shell_type& sh);
        void operator()(const conical_surface_shell_type& sh);

        std::vector<domain_type>& results;

      private:
        SGFRDSimulator& sim;
        const domain_type& dom;
    };
    struct pair_burster : boost::static_visitor<void>
    {
        pair_burster(SGFRDSimulator& s, const domain_type& d,
                     std::vector<domain_type>& r)
            : results(r), sim(s), dom(d)
        {}

        void operator()(const circular_shell_type& sh);
        void operator()(const conical_surface_shell_type& sh);

        std::vector<domain_type>& results;

      private:
        SGFRDSimulator& sim;
        const domain_type& dom;
    };

    struct single_escapement : boost::static_visitor<void>
    {
        typedef boost::container::static_vector<
            boost::tuple<ParticleID, Particle, face_id_type>, 1> executed_type;

        single_escapement(SGFRDSimulator& s, const Single& d)
            : sim(s), dom(d)
        {}

        void operator()(const circular_shell_type& sh);
        void operator()(const conical_surface_shell_type& sh);

        executed_type remnants;

      private:
        SGFRDSimulator& sim;
        Single const& dom;
    };

    struct single_reactor : boost::static_visitor<void>
    {
        typedef boost::container::static_vector<
            boost::tuple<ParticleID, Particle, face_id_type>, 2> executed_type;

        single_reactor(SGFRDSimulator& s, const Single& d)
            : sim(s), dom(d)
        {}

        void operator()(const circular_shell_type& sh);
        void operator()(const conical_surface_shell_type& sh);

        executed_type remnants;

      private:
        SGFRDSimulator& sim;
        Single const& dom;
    };
    // }}}

    //! make event from domain and push it into scheduler
    template<typename domainT>
    void add_event(const domainT& dom)
    {
        const DomainID did = scheduler_.add(boost::make_shared<event_type>(
                                            dom.begin_time() + dom.dt(), dom));
        domain_id_setter didset(did);
        mut_sh_vis_applier(didset, dom);
        return;
    }

    void fire_event(const event_id_pair_type& ev)
    {
        DUMP_MESSAGE("fire_event");
        domain_firer firer(*this);
        boost::apply_visitor(firer, ev.second->domain());
        return ;
    }

    void burst_event(const event_type& ev, std::vector<domain_type>& results)
    {
        domain_burster burster(*this, results);
        boost::apply_visitor(burster, ev.domain());
        return ;
    }

    void burst_non_multis(std::vector<domain_type>& domains)
    {
        std::cerr << "[WARNING]: burst_non_multis is not implemented yet"
                  << std::endl;
        return;
    }

    template<typename EventExecutor>
    void execute_event(EventExecutor& executor, const ShellID& sid)
    {
        boost::apply_visitor(executor, get_shell(sid));
        remove_shell(sid);

        ParticleID pid; Particle p; face_id_type fid;
        for(typename EventExecutor::executed_type::const_iterator
            iter(executor.remnants.begin()), end(executor.remnants.end());
            iter != end; ++iter)
        {
            boost::tie(pid, p, fid) = *iter;
            this->create_event(pid, p, fid);
        }
        return;
    }

    ShellID create_minimum_shell(
            const ParticleID& pid, const Particle& p, const face_id_type fid)
    {
        const ShellID sid(shell_id_gen());
        circular_shell_type sh(circle_type(p.radius(), p.position(),
                               this->polygon().triangle_at(fid).normal()), fid);
        shell_container_.add_shell(sid, sh, fid);
        return sid;
    }
    Single create_minimum_domain(
            const ShellID& sid, const ParticleID& pid, const Particle& p)
    {
        return Single(Single::ESCAPE, 0., this->time(), sid, std::make_pair(pid, p));
    }

    std::pair<ShellID, circle_type>
    create_single_circular_shell(const std::pair<Real3, face_id_type>& pos,
                                 const Real size)
    {
        DUMP_MESSAGE("create single circular shell");
        const ShellID id(shell_id_gen());
        const circle_type shape(size, pos.first,
                                polygon().triangle_at(pos.second).normal());
        shell_container_.add_shell(id, circular_shell_type(shape, pos.second),
                                   pos.second);
        return std::make_pair(id, shape);
    }

    std::pair<ShellID, conical_surface_type>
    create_single_conical_surface_shell(const vertex_id_type& vid,
                                        const Real size)
    {
        DUMP_MESSAGE("create single conical surface shell");
        const ShellID id(shell_id_gen());
        const conical_surface_type shape(polygon().vertex_at(vid).position,
                                         polygon().apex_angle(vid), size);
        shell_container_.add_shell(
                id, conical_surface_shell_type(shape, vid), vid);

        return std::make_pair(id, shape);
    }

    //! create domain and determine EventKind using GF
    Single create_single(const std::pair<ShellID, circle_type>& sh,
            const ParticleID& pid, const Particle& p)
    {
        //TODO consider single-reaction
        DUMP_MESSAGE("create single domain having circular shell");

        const greens_functions::GreensFunction2DAbsSym
            gf(/* D = */ p.D(),
               /* a = */ sh.second.size() - p.radius());
        const Real dt = gf.drawTime(uniform_real());
        DUMP_MESSAGE("delta t calculated: " << dt);

        return Single(Single::ESCAPE, dt, this->time(), sh.first,
                      std::make_pair(pid, p));
    }
    Single create_single(const std::pair<ShellID, conical_surface_type>& sh,
            const ParticleID& pid, const Particle& p)
    {
        //TODO consider single-reaction
        DUMP_MESSAGE("create single domain having conical shell");
        DUMP_MESSAGE("shell size = " << sh.second.size());
        DUMP_MESSAGE("D   = " << p.D());
        DUMP_MESSAGE("r0  = " << length(p.position() - sh.second.apex()));
        DUMP_MESSAGE("a   = " << sh.second.size() - p.radius());
        DUMP_MESSAGE("phi = " << sh.second.apex_angle());

        const greens_functions::GreensFunction2DRefWedgeAbs
            gf(/* D   = */ p.D(),
               /* r0  = */ length(p.position() - sh.second.apex()),
               /* a   = */ sh.second.size() - p.radius(),
               /* phi = */ sh.second.apex_angle());
        const Real dt = gf.drawTime(uniform_real());

        return Single(Single::ESCAPE, dt, this->time(), sh.first,
                      std::make_pair(pid, p));
    }

    //! make domain and call add_event
    void create_event(
            const ParticleID& pid, const Particle& p, const face_id_type fid);

    std::vector<std::pair<vertex_id_type, Real> >
    get_intrusive_vertices(const std::pair<Real3, face_id_type>& pos,
                           const Real radius) const
    {
        return polygon().list_vertices_within_radius(pos, radius);
    }

    std::vector<std::pair<DomainID, Real> >
    get_intrusive_domains(const std::pair<Real3, face_id_type>& pos,
                          const Real radius) const
    {
        const std::vector<std::pair<std::pair<ShellID, shell_type>, Real>
            > shells(shell_container_.list_shells_within_radius(pos, radius));

        std::vector<std::pair<DomainID, Real> > domains(shells.size());

        std::vector<std::pair<DomainID, Real> >::iterator dest = domains.begin();
        for(std::vector<std::pair<std::pair<ShellID, shell_type>, Real>
            >::const_iterator iter(shells.begin()), end(shells.end());
            iter != end; ++iter)
        {
            *dest = std::make_pair(
                    boost::apply_visitor(domain_id_getter(), iter->first.second),
                    iter->second);
            ++dest;
        }
        return domains;
    }

    std::vector<std::pair<DomainID, Real> >
    get_intrusive_domains(const vertex_id_type& vid, const Real radius) const
    {
        const std::pair<Real3, vertex_id_type> vpos = std::make_pair(
                polygon().vertex_at(vid).position, vid);
        const std::vector<std::pair<std::pair<ShellID, shell_type>, Real>
            > shells(shell_container_.list_shells_within_radius(vpos, radius));

        std::vector<std::pair<DomainID, Real> > domains(shells.size());

        std::vector<std::pair<DomainID, Real> >::iterator dest = domains.begin();
        for(std::vector<std::pair<std::pair<ShellID, shell_type>, Real>
            >::const_iterator iter(shells.begin()), end(shells.end());
            iter != end; ++iter)
        {
            *dest = std::make_pair(
                    boost::apply_visitor(domain_id_getter(), iter->first.second),
                    iter->second);
            ++dest;
        }
        return domains;
    }

    //! just a geometric restriction
    Real get_max_circle_size(const std::pair<Real3, face_id_type>& pos) const
    {
        Real lensq = std::numeric_limits<Real>::max();
        const boost::array<std::pair<Real3, Real3>, 6>& barrier =
            polygon().face_at(pos.second).barrier;

        for(std::size_t i=0; i<6; ++i)
        {
            const Real dist2 = distance_sq_to_segment(pos.first, barrier[i]);
            if(dist2 < lensq)
            {
                lensq = dist2;
            }
        }
        return std::sqrt(lensq);
    }
    Real get_max_cone_size(const vertex_id_type& vid) const
    {
        return polygon().vertex_at(vid).max_conical_shell_size * 0.5;
    }

    Real distance_sq_to_segment(const Real3& pos, const std::pair<Real3, Real3>& seg) const;

  private:

    static const Real single_circular_shell_factor;
    static const Real single_circular_shell_mergin;
    static const Real single_conical_surface_shell_factor;
    static const Real single_conical_surface_shell_mergin;

  private:

    // from SimulatorBase
    // boost::shared_ptr<model_type> model_;
    // boost::shared_ptr<world_type> world_;
    // Integer num_steps_;
    Real dt_;
    Real bd_dt_factor_;
    ecell4::RandomNumberGenerator&     rng_;
    scheduler_type                     scheduler_;
    shell_id_generator_type            shell_id_gen;
    shell_container_type               shell_container_;
    mutable_shell_visitor_applier_type mut_sh_vis_applier;
    std::vector<std::pair<reaction_rule_type, reaction_info_type> > last_reactions_;
};

} // sgfrd
} // ecell4

#endif // ECELL4_SGFRD_SIMULATOR
