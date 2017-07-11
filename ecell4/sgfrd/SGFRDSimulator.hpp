#ifndef ECELL4_SGFRD_SIMULATOR
#define ECELL4_SGFRD_SIMULATOR

#include <greens_functions/GreensFunction2DAbsSym.hpp>
#include <greens_functions/GreensFunction2DRefWedgeAbs.hpp>

#include <ecell4/core/SimulatorBase.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/geometry.hpp>

#include <ecell4/sgfrd/make_visitor.hpp>
#include <ecell4/sgfrd/ShellContainer.hpp>
#include <ecell4/sgfrd/ShellVisitorApplier.hpp>
#include <ecell4/sgfrd/ShellVisitors.hpp>
#include <ecell4/sgfrd/Informations.hpp>
#include <ecell4/sgfrd/SGFRDEvent.hpp>
#include <ecell4/sgfrd/SGFRDWorld.hpp>

#include <boost/make_shared.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/container/static_vector.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>

#include <iostream>

#ifndef SGFRD_NO_DEBUG
// #define SGFRD_LOG( sev, x ) BOOST_LOG_SEV(this->logger_, sev) << (x)
#define SGFRD_LOG( sev, x ) std::cerr << #sev << ':' << (x) << std::endl;
#endif // SGFRD_NO_DEBUG

namespace ecell4
{
namespace sgfrd
{

using boost::log::trivial::trace;
using boost::log::trivial::debug;
using boost::log::trivial::info;
using boost::log::trivial::warning;
using boost::log::trivial::error;
using boost::log::trivial::fatal;

class SGFRDSimulator :
    public ecell4::SimulatorBase<ecell4::Model, SGFRDWorld>
{
  public:
    typedef SGFRDSimulator self_type;

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
    typedef face_id_type   FaceID;// just for looks same as ParticleID
    typedef edge_id_type   EdgeID;
    typedef vertex_id_type VertexID;

    // Event & Domain
    typedef SGFRDEvent                      event_type;
    typedef SGFRDEvent::domain_type         domain_type;
    typedef EventID                         event_id_type;
    typedef DomainID                        domain_id_type;
    typedef SGFRDEventScheduler             scheduler_type;
    typedef SGFRDEventScheduler::value_type event_id_pair_type;

    // Simulator
    typedef ecell4::SimulatorBase<ecell4::Model, SGFRDWorld> base_type;
    typedef base_type::world_type world_type;
    typedef base_type::model_type model_type;
    typedef std::pair<ParticleID, Particle> particle_id_pair_type;
    typedef boost::tuple<ParticleID, Particle, FaceID> pid_p_fid_tuple_type;
    // length of this type may be tuned
    typedef boost::container::small_vector<pid_p_fid_tuple_type, 5> bursted_type;

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
    typedef shell_visitor_applier<shell_container_type>
            mutable_shell_visitor_applier_type;
    typedef shell_visitor_applier<const shell_container_type>
            immutable_shell_visitor_applier_type;

    // reaction
    typedef ecell4::ReactionRule           reaction_rule_type;
    typedef MoleculeInfo                   molecule_info_type;
    typedef ReactionInfo                   reaction_info_type;
    typedef std::pair<reaction_rule_type, reaction_info_type> reaction_log_type;
    typedef std::vector<reaction_log_type> reaction_archive_type;

  public:

    SGFRDSimulator(const boost::shared_ptr<world_type>& world,
                   const boost::shared_ptr<model_type>& model,
                   Real bd_dt_factor = 1e-5)
        : base_type(model, world), dt_(0), bd_dt_factor_(bd_dt_factor),
          rng_(*(world->rng())), shell_container_(*(world->polygon())),
          mut_sh_vis_applier(shell_container_), imm_sh_vis_applier(shell_container_)
    {}

    SGFRDSimulator(boost::shared_ptr<world_type> world, Real bd_dt_factor = 1e-5)
        : base_type(world), dt_(0), bd_dt_factor_(bd_dt_factor),
          rng_(*(world->rng())), shell_container_(*(world->polygon())),
          mut_sh_vis_applier(shell_container_), imm_sh_vis_applier(shell_container_)
    {}

    ~SGFRDSimulator(){}

    void initialize()
    {
        ParticleID pid; Particle p;
        BOOST_FOREACH(boost::tie(pid, p), this->world_->list_particles())
        {
            add_event(create_closely_fitted_domain(create_closely_fitted_shell(
                      pid, p, this->get_face_id(pid)), pid, p));
        }
        return ;
    }
    void finalize()
    {
        const Real tm(this->time());
        while(this->scheduler_.size() != 0)
        {
            this->burst_event(this->scheduler_.pop(), tm);
        }
        return ;
    }
    void step()
    {
        this->set_time(this->scheduler_.next_time());
        // fire event executes `create_event` inside.
        this->fire_event(this->scheduler_.pop());
        SGFRD_LOG(debug, boost::format("now %1% shells exist") %
                  shell_container_.num_shells());
        SGFRD_LOG(debug, boost::format("now %1% events exist") %
                  scheduler_.size());
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

  private:

    // simple wrappers to call member's member-method (e.g. world_->t()) {{{
    Real uniform_real(){return this->rng_.random();}

    world_type   const& world()   const {return *(this->world_);}
    polygon_type const& polygon() const {return *(this->world_->polygon());}

    bool update_particle(const ParticleID& pid, const Particle& p,
                         const FaceID& fid)
    {
        SGFRD_LOG(debug, boost::format("particle %1% is updated") % pid);
        return this->world_->update_particle(pid, p, fid);
    }

    shell_type&       get_shell(ShellID const& id)
    {
        SGFRD_LOG(debug, boost::format("searching shell %1%") % id);
        return shell_container_.get_shell(id);
    }
    shell_type const& get_shell(ShellID const& id) const
    {
        SGFRD_LOG(debug, boost::format("searching shell %1%") % id);
        return shell_container_.get_shell(id);
    }
    void remove_shell(ShellID const& id)
    {
        SGFRD_LOG(debug, boost::format("removing shell %1%") % id);
        return shell_container_.remove_shell(id);
    }
    template<typename shT, typename stridT>
    void update_shell(ShellID const& shid, shT const& sh, stridT strid)
    {
        this->shell_container_.update_shell(shid, sh, strid);
        return;
    }


    FaceID get_face_id(const ParticleID& pid) const
    {return this->world_->get_face_id(pid);}
    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const
    {return this->world_->get_particle(pid);}

    Real time() const {return this->world_->t();}
    void set_time(const Real t) {return this->world_->set_t(t);}

    boost::shared_ptr<event_type> get_event(const event_id_type& id)
    {
        SGFRD_LOG(debug, boost::format("getting event %1%") % id);
        return scheduler_.get(id);
    }
    boost::shared_ptr<event_type> pickout_event(const event_id_type& id)
    {
        SGFRD_LOG(debug, boost::format("picking out event %1%") % id);
        BOOST_AUTO(tmp, scheduler_.get(id));
        scheduler_.remove(id);
        return tmp;
    }
    void remove_event(const event_id_type id)
    {
        SGFRD_LOG(debug, boost::format("removing event %1%") % id);
        this->scheduler_.remove(id);
        return;
    }

    DomainID get_domain_id(Single const& dom) const
    {
        return boost::apply_visitor(domain_id_getter(), get_shell(dom.shell_id()));
    }
    DomainID get_domain_id(Pair const& dom) const
    {
        return boost::apply_visitor(domain_id_getter(), get_shell(dom.shell_id()));
    }
    DomainID get_domain_id(Multi const& dom) const
    {
        return boost::apply_visitor(domain_id_getter(),
                                    get_shell(dom.shell_ids().front()));
    }
    // }}}

  private:

//----------------------------------- single -----------------------------------

    /*! execute Event associated with Domain, remove shell, create next event. */
    void          fire_single(const Single& dom, DomainID did);
    bursted_type burst_single(const Single& dom, const Real tm);

    template<typename shellT>
    boost::tuple<ParticleID, Particle, FaceID>
    propagate_single(const shellT& sh, const Single& dom, const Real tm);

    template<typename shellT>
    boost::tuple<ParticleID, Particle, FaceID>
    escape_single(const shellT& sh, const Single& dom);

    template<typename shellT>
    boost::container::static_vector<pid_p_fid_tuple_type, 2>
    reaction_single(const shellT& sh, const Single& dom, const DomainID did);

    template<typename shellT>
    boost::container::static_vector<pid_p_fid_tuple_type, 2>
    attempt_reaction_single(const shellT& sh, const DomainID did,
            const ParticleID& pid, const Particle& p, const FaceID& fid);

    std::pair<ShellID, circle_type>
    create_single_circular_shell(
            const std::pair<Real3, FaceID>& pos, const Real size)
    {
        SGFRD_LOG(trace, "create single circular shell");
        const ShellID id(shell_id_gen());
        const circle_type shape(size, pos.first,
                                polygon().triangle_at(pos.second).normal());
        shell_container_.add_shell(id, circular_shell_type(shape, pos.second),
                                   pos.second);
        return std::make_pair(id, shape);
    }
    std::pair<ShellID, conical_surface_type>
    create_single_conical_surface_shell(
            const vertex_id_type& vid, const Real size)
    {
        SGFRD_LOG(trace, "create single conical surface shell");
        const ShellID id(shell_id_gen());
        const conical_surface_type shape(polygon().vertex_at(vid).position,
                                         polygon().apex_angle(vid), size);
        shell_container_.add_shell(
                id, conical_surface_shell_type(shape, vid), vid);

        return std::make_pair(id, shape);
    }

    Single create_single(const std::pair<ShellID, circle_type>& sh,
                         const ParticleID& pid, const Particle& p)
    {//TODO consider single-reaction
        SGFRD_LOG(trace, "create single domain having circular shell");

        const greens_functions::GreensFunction2DAbsSym
            gf(/* D = */ p.D(),
               /* a = */ sh.second.size() - p.radius());
        const Real dt = gf.drawTime(uniform_real());
        SGFRD_LOG(debug, boost::format("calculated delta t = %1%") % dt);

        return Single(Single::ESCAPE, dt, this->time(), sh.first,
                      std::make_pair(pid, p));
    }
    Single create_single(const std::pair<ShellID, conical_surface_type>& sh,
                         const ParticleID& pid, const Particle& p)
    {//TODO consider single-reaction
        SGFRD_LOG(trace, "create single domain having conical shell");
        SGFRD_LOG(debug, boost::format("shell size = %1%") % sh.second.size());
        SGFRD_LOG(debug, boost::format("D   = %1%") % p.D());
        SGFRD_LOG(debug, boost::format("r0  = %1%") %
                  length(p.position() - sh.second.apex()));
        SGFRD_LOG(debug, boost::format("a   = %1%") % (sh.second.size() - p.radius()));
        SGFRD_LOG(debug, boost::format("phi = %1%") % sh.second.apex_angle());

        const greens_functions::GreensFunction2DRefWedgeAbs
            gf(/* D   = */ p.D(),
               /* r0  = */ length(p.position() - sh.second.apex()),
               /* a   = */ sh.second.size() - p.radius(),
               /* phi = */ sh.second.apex_angle());
        const Real dt = gf.drawTime(uniform_real());

        return Single(Single::ESCAPE, dt, this->time(), sh.first,
                      std::make_pair(pid, p));
    }

//------------------------------------ pair ------------------------------------

    void fire_pair(const Pair& dom, DomainID did)
    {// TODO
        SGFRD_LOG(warning, "fire_pair has not been implemented yet.");
    }

    bursted_type burst_pair(const Pair& dom, const Real tm)
    {// TODO
        SGFRD_LOG(warning, "burst_pair has not been implemented yet.");
    }

//----------------------------------- multi ------------------------------------

    void fire_multi(Multi& dom, DomainID did)
    {
        SGFRD_LOG(trace, boost::format("fire multi for default dt(%1%)") % dom.dt());

        volume_clearer vc(did, dom, *this, this->imm_sh_vis_applier);
        dom.step(vc);
        switch(dom.eventkind())
        {
            case Multi::NONE:
            {
                SGFRD_LOG(trace, "nothing occurs");
                /* continuing multi domain: add this domain to scheduler */
                dom.begin_time() = this->time();
                this->add_event(dom);
                return;
            }
            case Multi::REACTION:
            {
                SGFRD_LOG(trace, "reaction occurs");
                std::copy(dom.last_reactions().begin(), dom.last_reactions().end(),
                          std::back_inserter(this->last_reactions_));
                ParticleID pid; Particle p; FaceID fid;
                BOOST_FOREACH(boost::tie(pid, p, fid),
                              this->remove_multi(dom))
                {
                    this->add_event(this->create_closely_fitted_domain(
                        this->create_closely_fitted_shell(pid, p, fid), pid, p));
                }
                SGFRD_LOG(trace, boost::format("multi domain (id = %1%) removed.") %
                          did);
                return;
            }
            case Multi::ESCAPE:
            {
                SGFRD_LOG(trace, "particle escapes");
                /* burst this domain! */
                ParticleID pid; Particle p; FaceID fid;
                BOOST_FOREACH(boost::tie(pid, p, fid),
                              this->remove_multi(dom))
                {
                    this->add_event(this->create_closely_fitted_domain(
                        this->create_closely_fitted_shell(pid, p, fid), pid, p));
                }
                SGFRD_LOG(trace, boost::format("multi domain (id = %1%) removed.") %
                          did);
                return;
            }
            default:
            {
                SGFRD_LOG(fatal, "Multi eventkind become invalid value!");
                throw std::logic_error("never reach here");
            }
        }
    }

    // simply remove all the shells. not add a domains for each particles.
    // particles are updated at each step, so here nothing is needed to
    // update world.
    bursted_type remove_multi(const Multi& dom)
    {
        SGFRD_LOG(trace, "remove multi called");
        bursted_type results;
        Particle p; ParticleID pid;
        BOOST_FOREACH(boost::tie(pid, p), dom.particles())
        {
            results.push_back(boost::make_tuple(pid, p, this->get_face_id(pid)));
        }
        SGFRD_LOG(trace, "particles are collected");
        BOOST_FOREACH(ShellID sid, dom.shell_ids())
        {
            this->remove_shell(sid);
        }
        SGFRD_LOG(trace, "shells are removed");
        return results;
    }

    // burst. step until(tm - dom.begin_time()).
    bursted_type burst_multi(Multi& dom, const Real tm)
    {
        SGFRD_LOG(trace, "burst_multi called");
        BOOST_AUTO(did, get_domain_id(dom));
        volume_clearer vc(did, dom, *this, this->imm_sh_vis_applier);
        dom.step(vc, tm - dom.begin_time());

        SGFRD_LOG(debug, boost::format("multi domain steps with delta_t = %1%") %
                  (tm - dom.begin_time()));

        if(dom.eventkind() == Multi::REACTION)
        {
            SGFRD_LOG(trace, "reaction occured");
            std::copy(dom.last_reactions().begin(), dom.last_reactions().end(),
                      std::back_inserter(this->last_reactions_));
        }
        return remove_multi(dom);
    }

// -----------------------------------------------------------------------------

    // XXX: second value of element of result_type is not a mere distance.
    //      - in Multi case, it is just a distance.
    //      - in bursted Single case, it become a distance between centers minus
    //        min_circular_shell of the particle
    std::vector<std::pair<DomainID, Real> >
    burst_and_shrink_non_multis(
            const ParticleID& pid, const Particle& p, const FaceID& fid,
            const std::vector<std::pair<DomainID, Real> >& intruders)
    {
        SGFRD_LOG(trace, "burst_and_shrink_non_multis called");
        const Real tm(this->time());
        std::vector<std::pair<DomainID, Real> > results;

        DomainID did; Real dist;
        BOOST_FOREACH(boost::tie(did, dist), intruders)
        {
            BOOST_AUTO(const& ev, get_event(did));
            if(ev->which_domain() == event_type::multi_domain)
            {
                SGFRD_LOG(trace, boost::format("domain %1% is multi") % did);
                results.push_back(std::make_pair(did, dist));
                continue;
            }
            SGFRD_LOG(trace, boost::format("domain %1% is single") % did);

            DomainID did_; ParticleID pid_; Particle p_; FaceID fid_;
            BOOST_FOREACH(boost::tie(pid_, p_, fid_),
                          burst_event(std::make_pair(did, ev), tm))
            {
                SGFRD_LOG(trace, boost::format(
                    "add closely-fitted domain to bursted particle %1%") % pid_);
                remove_event(did);
                did_ = add_event(create_closely_fitted_domain(
                    create_closely_fitted_shell(pid_, p_, fid_), pid_, p_));
                results.push_back(std::make_pair(did_, this->polygon().distance(
                    std::make_pair(p.position(),  fid),
                    std::make_pair(p_.position(), fid_)) -
                    calc_min_single_circular_shell_radius(p_)));
            }
            SGFRD_LOG(trace, boost::format("domain %1% is bursted and shrinked") % did);
        }

        std::sort(results.begin(), results.end(),
            ecell4::utils::pair_second_element_comparator<DomainID, Real>());
        SGFRD_LOG(trace, "results are sorted");
        return results;
    }

    // to clear volume. burst all the overlapping shells then add closely-fitted
    // shells to them. returns true if there are no overlapping particles.
    bool burst_and_shrink_overlaps(const Particle& p, const FaceID& fid);

    // form multi shell recursively
    DomainID form_multi(const ParticleID& pid, const Particle& p, const FaceID& fid,
                        const std::vector<std::pair<DomainID, Real> >& doms);

    // search intruder for multi, burst them and add them to multi if needed.
    void add_to_multi_recursive(Multi&);

    void merge_multi(Multi& from, Multi& to)
    {
        SGFRD_LOG(trace, "merge_multi called");
        // reset domain_id
        const domain_id_setter didset(this->get_domain_id(to));
        mut_sh_vis_applier(didset, from);

        // move particle
        ParticleID pid; Particle p; FaceID fid;
        BOOST_FOREACH(boost::tie(pid, p), from.particles())
        {
            const bool addp_result = to.add_particle(pid);
            assert(addp_result);
        }
        // move shell
        BOOST_FOREACH(ShellID sid, from.shell_ids())
        {
            const bool adds_result = to.add_shell(sid);
            assert(adds_result);
        }
        to.determine_reaction_length();
        to.determine_delta_t();

        remove_event(get_domain_id(from));
        SGFRD_LOG(trace, "multi from is removed");
        return;
    }

    struct volume_clearer
    {
        volume_clearer(domain_id_type d, const Multi& dom, SGFRDSimulator& s,
                       immutable_shell_visitor_applier_type& imm)
            : sim(s), did(d), domain(dom), applier(imm), logger_(sim.logger_)
        {}

        bool operator()(const Particle& p, const FaceID& fid)
        {
            SGFRD_LOG(trace, "volume clearer(p, fid) called");
            escaped_ = false;
            inside_checker is_inside(p.position(), fid, sim.polygon());
            if(applier(is_inside, domain)) return true;

            SGFRD_LOG(trace, "particle escaped");

            const bool no_overlap = sim.burst_and_shrink_overlaps(p, fid);
            escaped_ = no_overlap;
            SGFRD_LOG(trace, boost::format("overlap exist ? %1%") % no_overlap);
            return no_overlap;
        }
        bool operator()(const Particle& p, const FaceID& fid,
                        const ParticleID& ignore)
        {
            SGFRD_LOG(trace,
                    boost::format("volume clearer(p, fid, ignore = %1%) called")
                    % ignore);
            escaped_ = false;
            inside_checker is_inside(p.position(), fid, sim.polygon());
            if(applier(is_inside, domain)) return true;

            SGFRD_LOG(trace, "particle escaped");
            const bool no_overlap = sim.burst_and_shrink_overlaps(p, fid);
            escaped_ = no_overlap;
            SGFRD_LOG(trace, boost::format("overlap exist ? %1%") % no_overlap);
            return no_overlap;
        }
        bool operator()(const Particle& p, const FaceID& fid,
                        const ParticleID& ignore1, const ParticleID& ignore2)
        {
            SGFRD_LOG(trace,
                boost::format("volume clearer(p, fid, ignore1 = %1%, ignore2 = %2%) called")
                % ignore1 % ignore2);
            escaped_ = false;
            inside_checker is_inside(p.position(), fid, sim.polygon());
            if(applier(is_inside, domain)) return true;

            SGFRD_LOG(trace, "particle escaped");
            const bool no_overlap = sim.burst_and_shrink_overlaps(p, fid);
            escaped_ = no_overlap;
            SGFRD_LOG(trace, boost::format("overlap exist ? %1%") % no_overlap);
            return no_overlap;
        }

        bool escaped() const {return escaped_;}

      private:
        bool            escaped_;
        SGFRDSimulator& sim;
        domain_id_type  did;
        Multi const&    domain;
        immutable_shell_visitor_applier_type applier;
        boost::log::sources::severity_logger<boost::log::trivial::severity_level>
            logger_;
    };

//----------------------------------- event ------------------------------------

    //! make event from domain and push it into scheduler
    template<typename domainT>
    DomainID add_event(const domainT& dom)
    {
        SGFRD_LOG(trace, "add_event called");
        const DomainID did = scheduler_.add(
            boost::make_shared<event_type>(dom.begin_time() + dom.dt(), dom));
        SGFRD_LOG(trace, boost::format("domain ID = %1%") % did);
        domain_id_setter didset(did);
        mut_sh_vis_applier(didset, dom);
        SGFRD_LOG(trace, "return from add_event");
        return did;
    }

    // assuming the event is already poped
    void fire_event(event_id_pair_type ev)
    {
        SGFRD_LOG(trace, "fire_event");
        return boost::apply_visitor(make_visitor(
            resolve<Single const&, void>(boost::bind(
                    &self_type::fire_single, this, _1, ev.first)),
            resolve<Pair   const&, void>(boost::bind(
                    &self_type::fire_pair,   this, _1, ev.first)),
            resolve<Multi&,        void>(boost::bind(
                    &self_type::fire_multi,  this, _1, ev.first))
            ), ev.second->domain());
    }

    // assuming the event is already poped
    bursted_type burst_event(const event_id_pair_type& ev, Real tm)
    {
        SGFRD_LOG(trace, "burst_event");
        return boost::apply_visitor(make_visitor(
            resolve<Single const&, bursted_type>(boost::bind(
                    &self_type::burst_single, this, _1, tm)),
            resolve<Pair   const&, bursted_type>(boost::bind(
                    &self_type::burst_pair,  this, _1,  tm)),
            resolve<Multi&,        bursted_type>(boost::bind(
                    &self_type::burst_multi, this, _1,  tm))
            ), ev.second->domain());
    }

    ShellID create_closely_fitted_shell(
            const ParticleID& pid, const Particle& p, const FaceID fid)
    {
        SGFRD_LOG(trace,
            boost::format("creating closely fitted shell for particle %1%") % pid);
        const ShellID sid(shell_id_gen());
        circular_shell_type sh(circle_type(p.radius(), p.position(),
                               this->polygon().triangle_at(fid).normal()), fid);
        shell_container_.add_shell(sid, sh, fid);
        return sid;
    }
    Single create_closely_fitted_domain(
            const ShellID& sid, const ParticleID& pid, const Particle& p)
    {
        SGFRD_LOG(trace,
            boost::format("creating closely fitted domain for particle %1%") % pid);
        return Single(Single::ESCAPE, 0., this->time(), sid, std::make_pair(pid, p));
    }

    std::pair<ShellID, circle_type>
    create_minimum_single_shell(
            const ParticleID& pid, const Particle& p, const FaceID fid)
    {
        return create_single_circular_shell(std::make_pair(p.position(), fid),
            calc_min_single_circular_shell_radius(p));
    }

    Multi create_empty_multi()
    {
        return Multi(*this, *(this->world_), this->time());
    }

    //! make domain and call add_event
    DomainID create_event(const ParticleID&, const Particle&, const FaceID);

    std::vector<std::pair<vertex_id_type, Real> >
    get_intrusive_vertices(const std::pair<Real3, FaceID>& pos,
                           const Real radius) const
    {
        return polygon().list_vertices_within_radius(pos, radius);
    }

    std::vector<std::pair<DomainID, Real> >
    get_intrusive_domains(const std::pair<Real3, FaceID>& pos,
                          const Real radius) const
    {
        const std::vector<std::pair<std::pair<ShellID, shell_type>, Real>
            > shells(shell_container_.list_shells_within_radius(pos, radius));

        std::vector<std::pair<DomainID, Real> > domains;
        domains.reserve(shells.size());

        std::pair<ShellID, shell_type> shell_id_pair; Real dist;
        BOOST_FOREACH(boost::tie(shell_id_pair, dist), shells)
        {
            const DomainID did = boost::apply_visitor(
                    domain_id_getter(), shell_id_pair.second);

            if(std::find_if(domains.begin(), domains.end(),
                    ecell4::utils::pair_first_element_unary_predicator<
                    DomainID, Real>(did)) == domains.end())
            {
                domains.push_back(std::make_pair(did, dist));
            }
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

        std::vector<std::pair<DomainID, Real> > domains;
        domains.reserve(shells.size());

        std::pair<ShellID, shell_type> shell_id_pair; Real dist;
        BOOST_FOREACH(boost::tie(shell_id_pair, dist), shells)
        {
            const DomainID did = boost::apply_visitor(
                    domain_id_getter(), shell_id_pair.second);

            if(std::find_if(domains.begin(), domains.end(),
                    ecell4::utils::pair_first_element_unary_predicator<
                    DomainID, Real>(did)) == domains.end())
            {
                domains.push_back(std::make_pair(did, dist));
            }
        }
        return domains;
    }

    //! just a geometric restriction
    Real get_max_circle_size(const std::pair<Real3, FaceID>& pos) const
    {
        Real lensq = std::numeric_limits<Real>::max();
        const boost::array<ecell4::Segment, 6>& barrier =
            polygon().face_at(pos.second).barrier;

        for(std::size_t i=0; i<6; ++i)
        {
            const Real dist2 = ecell4::sgfrd::distance_sq(pos.first, barrier[i]);
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

    Real calc_min_single_circular_shell_radius(const Particle& p)
    {
        return p.radius() * single_circular_shell_factor;
    }
    Real calc_min_single_conical_shell_radius(const Particle& p)
    {
        return p.radius() * single_conical_surface_shell_factor;
    }

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
    ecell4::RandomNumberGenerator&       rng_;
    scheduler_type                       scheduler_;
    shell_id_generator_type              shell_id_gen;
    shell_container_type                 shell_container_;
    mutable_shell_visitor_applier_type   mut_sh_vis_applier;
    immutable_shell_visitor_applier_type imm_sh_vis_applier;
    std::vector<std::pair<reaction_rule_type, reaction_info_type> > last_reactions_;

    mutable
    boost::log::sources::severity_logger<boost::log::trivial::severity_level>
        logger_;
};

// XXX NOTE XXX:
// To avoid an error "specialization of template function after instantiation"
// these specialized functions are must be implemented before other member
// (maybe template) function directory calling it.
template<>
inline boost::tuple<ParticleID, Particle, SGFRDSimulator::FaceID>
SGFRDSimulator::propagate_single<SGFRDSimulator::circular_shell_type>(
        const circular_shell_type& sh, const Single& dom, const Real tm)
{
    SGFRD_LOG(trace, "propagating single circular shell");

    Particle   p   = dom.particle();
    ParticleID pid = dom.particle_id();

    greens_functions::GreensFunction2DAbsSym gf(p.D(), sh.size() - p.radius());

    const Real del_t = tm - dom.begin_time();
    const Real r     = gf.drawR(this->uniform_real(), del_t);
    const Real theta = this->uniform_real() * 2 * M_PI;
    SGFRD_LOG(debug, boost::format("r = %1%, theta = %2%") % r % theta);

    const FaceID         fid  = this->get_face_id(pid);
    const triangle_type& face = this->polygon().triangle_at(fid);
    const Real3 direction = rotate(theta, face.normal(), face.represent());
    const Real  len_direction = length(direction);
    SGFRD_LOG(debug, boost::format("direction = %1%, len = %2%")
                                   % direction % len_direction);

    std::pair<std::pair<Real3, FaceID>, Real3> state =
        std::make_pair(std::make_pair(p.position(), fid),
                       direction * r / len_direction);
    SGFRD_LOG(debug, boost::format("pos = %1%, fid = %2%")
              % state.first.first % state.first.second);

    unsigned int continue_count = 2;
    while(continue_count > 0)
    {
        state = this->polygon().move_next_face(state.first, state.second);
        const Real3& disp = state.second;
        if(disp[0] == 0. && disp[1] == 0. && disp[2] == 0.) break;
        --continue_count;
        SGFRD_LOG(debug, boost::format("pos = %1%, fid = %2%, dsp=%3%, count=%4%")
                  % state.first.first % state.first.second % disp % continue_count);
    }
    if(continue_count == 0)
    {
        SGFRD_LOG(warning, "moving on face: precision lost");
    }

    SGFRD_LOG(debug, boost::format("pos = %1%, fid = %2%, dsp = %3%, count = %4%")
          % state.first.first % state.first.second % state.second % continue_count);

    p.position() = state.first.first;
    this->update_particle(pid, p, state.first.second);

    SGFRD_LOG(trace, "propagated");
    return boost::make_tuple(pid, p, state.first.second);
}

template<>
inline boost::tuple<ParticleID, Particle, SGFRDSimulator::FaceID>
SGFRDSimulator::propagate_single<SGFRDSimulator::conical_surface_shell_type>(
        const conical_surface_shell_type& sh, const Single& dom, const Real tm)
{
    SGFRD_LOG(trace, "propagating single conical");
    Particle         p   = dom.particle();
    const ParticleID pid = dom.particle_id();
    const FaceID     fid = this->get_face_id(pid);
    SGFRD_LOG(debug, boost::format("pos = %1%, fid = %2%") % p.position() % fid);

    const Real r_max = sh.size() - p.radius();
    greens_functions::GreensFunction2DRefWedgeAbs
        gf(/* D   = */ p.D(),
           /* r0  = */ length(p.position() - sh.position()),
           /* a   = */ r_max,
           /* phi = */ sh.shape().apex_angle());

    const Real del_t = tm - dom.begin_time();
    const Real r     = gf.drawR(this->uniform_real(), del_t);
    const Real theta = gf.drawTheta(this->uniform_real(), r, del_t);

    SGFRD_LOG(debug, boost::format("r = %1%, theta = %2%") % r % theta);

    const std::pair<Real3, FaceID> state =
        this->polygon().rotate_around_vertex(std::make_pair(p.position(), fid),
                                           sh.structure_id(), r, theta);
    SGFRD_LOG(debug, boost::format("propagateed : pos = %1%, fid = %2%")
              % state.first % state.second);

    p.position() = state.first;
    this->update_particle(pid, p, state.second);
    return boost::make_tuple(pid, p, state.second);
}

template<>
inline boost::tuple<ParticleID, Particle, SGFRDSimulator::FaceID>
SGFRDSimulator::escape_single<SGFRDSimulator::circular_shell_type>(
        const circular_shell_type& sh, const Single& dom)
{
    SGFRD_LOG(trace, "single shell escapement circular shell");
    if(sh.size() == dom.particle().radius())
    {
        SGFRD_LOG(trace, "closely fitted shell. didnot move.");
        return boost::make_tuple(dom.particle_id(), dom.particle(),
                                 this->get_face_id(dom.particle_id()));
    }

    Particle   p   = dom.particle();
    ParticleID pid = dom.particle_id();

    const Real r   = sh.size() - p.radius();
    const Real theta = this->uniform_real() * 2.0 * M_PI;
    SGFRD_LOG(debug, boost::format("r = %1%, theta = %2%") % r % theta);
    const FaceID         fid  = this->get_face_id(pid);
    const triangle_type& face = this->polygon().triangle_at(fid);
    const Real3 direction = rotate(theta, face.normal(), face.represent());
    SGFRD_LOG(debug, boost::format("dir = %1%, len = %2%")
              % direction % length(direction));

    std::pair<std::pair<Real3, FaceID>, Real3> state =
        std::make_pair(/*position = */std::make_pair(p.position(), fid),
                   /*displacement = */direction * r / length(direction));

    SGFRD_LOG(debug, boost::format("pos = %1%, fid = %2%")
              % state.first.first % state.first.second);
    unsigned int continue_count = 2;
    while(continue_count > 0)
    {
        state = this->polygon().move_next_face(state.first, state.second);
        const Real3& disp = state.second;
        if(disp[0] == 0. && disp[1] == 0. && disp[2] == 0.) break;
        --continue_count;
        SGFRD_LOG(debug, boost::format("pos = %1%, fid = %2%, disp = %3%")
                  % state.first.first % state.first.second % disp);
    }
    if(continue_count == 0)
    {
        SGFRD_LOG(warning, "moving on face: precision lost");
    }
    SGFRD_LOG(trace, "escaped.");
    SGFRD_LOG(debug, boost::format("pos = %1%, fid = %2%")
              % state.first.first % state.first.second);

    p.position() = state.first.first;
    this->update_particle(pid, p, state.first.second);
    return boost::make_tuple(pid, p, state.first.second);
}

template<>
inline boost::tuple<ParticleID, Particle, SGFRDSimulator::FaceID>
SGFRDSimulator::escape_single<SGFRDSimulator::conical_surface_shell_type>(
        const conical_surface_shell_type& sh, const Single& dom)
{
    SGFRD_LOG(trace, "single shell escapement: conical shell");
    Particle           p   = dom.particle();
    const ParticleID   pid = dom.particle_id();
    const FaceID       fid = this->get_face_id(pid);
    SGFRD_LOG(debug, boost::format("pos = %1%, fid = %2%") % p.position() % fid);

    const Real r     = sh.size() - p.radius();
    greens_functions::GreensFunction2DRefWedgeAbs
        gf(p.D(), length(p.position() - sh.position()),
           r,     sh.shape().apex_angle());
    const Real theta = gf.drawTheta(this->uniform_real(), r, dom.dt());
    SGFRD_LOG(debug, boost::format("r = %1%, theta = %2%") % r % theta);

    const std::pair<Real3, FaceID> state =
        this->polygon().rotate_around_vertex(std::make_pair(p.position(), fid),
                                           sh.structure_id(), r, theta);

    SGFRD_LOG(trace, "escaped");
    SGFRD_LOG(debug, boost::format("pos = %1%, fid = %2%") % p.position() % fid);

    p.position() = state.first;
    this->update_particle(pid, p, state.second);
    SGFRD_LOG(trace, "particle updated");
    return boost::make_tuple(pid, p, state.second);
}

inline SGFRDSimulator::bursted_type
SGFRDSimulator::burst_single(const Single& dom, const Real tm)
{
    SGFRD_LOG(trace, "bursting single");
    const ShellID sid(dom.shell_id());
    bursted_type results;
    results.push_back(boost::apply_visitor(make_visitor(
        resolve<const circular_shell_type&,
                boost::tuple<ParticleID, Particle, FaceID> >(boost::bind(
            &self_type::propagate_single<circular_shell_type>,
            this, _1, dom, tm)),
        resolve<const conical_surface_shell_type&,
                boost::tuple<ParticleID, Particle, FaceID> >(boost::bind(
            &self_type::propagate_single<conical_surface_shell_type>,
            this, _1, dom, tm))
        ), get_shell(sid)));
    this->remove_shell(sid);
    SGFRD_LOG(trace, "shell removed");
    return results;
}

// calls propagate_single.
template<typename shellT>
boost::container::static_vector<
    boost::tuple<ParticleID, Particle, SGFRDSimulator::FaceID>, 2>
SGFRDSimulator::reaction_single(
        const shellT& sh, const Single& dom, const DomainID did)
{
    SGFRD_LOG(trace, "reaction_single called");
    ParticleID pid; Particle p; FaceID fid;
    boost::tie(pid, p, fid) = this->propagate_single(sh, dom, dom.dt());// XXX
    boost::container::static_vector<
        boost::tuple<ParticleID, Particle, FaceID>, 2> retval;
    return retval;
    // TODO make clear_volume() for single reaction
//         volume_clearer vc(did, dom, *this, this->imm_sh_vis_applier);
//         return this->attempt_reaction_single(sh, pid, p, fid, vc);
}


// TODO
// calls propagate_single.
template<typename shellT>
boost::container::static_vector<SGFRDSimulator::pid_p_fid_tuple_type, 2>
SGFRDSimulator::attempt_reaction_single(const shellT& sh, const DomainID did,
        const ParticleID& pid, const Particle& p, const FaceID& fid)
{
    SGFRD_LOG(trace,   "attempt_reaction_single called");
    SGFRD_LOG(warning, "attempt_reaction_single has not been implemented yet");
    boost::container::static_vector<pid_p_fid_tuple_type, 2> retval;
    return retval;
}

} // sgfrd
} // ecell4
#endif // ECELL4_SGFRD_SIMULATOR
