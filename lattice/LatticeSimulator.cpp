#include "LatticeSimulator.hpp"

//#include <set>

namespace ecell4
{

namespace lattice
{

void LatticeSimulator::initialize()
{
    scheduler_.clear();
    std::vector<Species> species(world_->list_species());
    for (std::vector<Species>::const_iterator itr(species.begin());
            itr != species.end(); ++itr)
    {
        register_step_event(*itr);
    }

    NetworkModel::reaction_rule_container_type rules(model_->reaction_rules());
    for (NetworkModel::reaction_rule_container_type::iterator itr(rules.begin());
            itr != rules.end(); ++itr)
    {
        if ((*itr).reactants().size() == 1)
        {
            const boost::shared_ptr<EventScheduler::Event> event(
                    create_first_order_reaction_event(*itr));
            scheduler_.add(event);
        }
    }

    is_initialized_ = true;
}

boost::shared_ptr<EventScheduler::Event> LatticeSimulator::create_step_event(
        const Species& species, const Real& t)
{
    boost::shared_ptr<EventScheduler::Event> event(new StepEvent(this, species, t));
    return event;
}

boost::shared_ptr<EventScheduler::Event>
LatticeSimulator::create_first_order_reaction_event(const ReactionRule& reaction_rule)
{
    boost::shared_ptr<EventScheduler::Event> event(new FirstOrderReactionEvent(
                this, reaction_rule));
    return event;
}

std::pair<bool, Reaction<Voxel> > LatticeSimulator::attempt_reaction_(
        LatticeWorld::particle_info& info, LatticeWorld::coordinate_type to_coord)
{
    const Species from_species(world_->get_molecular_type(info.first)->species());
    const MolecularTypeBase* to_mt(world_->get_molecular_type(to_coord));
    const Species to_species(to_mt->species());
    const std::vector<ReactionRule> rules(model_->query_reaction_rules(
                from_species, to_species));

    const Real factor((from_species == to_species) ? 2 : 1);
    const Real Da(boost::lexical_cast<Real>(from_species.get_attribute("D"))),
               Db(boost::lexical_cast<Real>(to_species.get_attribute("D")));
    const Real rnd(world_->rng()->uniform(0,1));
    Real accp(0.);
    for (std::vector<ReactionRule>::const_iterator itr(rules.begin());
            itr != rules.end(); ++itr)
    {
        const Real k((*itr).k());
        const Real P(k * factor/ (6 * sqrt(2.)
                    * (Da + Db) * world_->voxel_radius()));
        accp += P;
        if (accp > 1)
        {
            std::cerr << "accp : " << accp << std::endl;
        }
        if (accp >= rnd)
        {
            LatticeWorld::particle_info to_info(*(to_mt->find(to_coord)));
            return apply_reaction_(*itr, info, to_info);
        }
    }
    return std::pair<bool, Reaction<Voxel> >(false, Reaction<Voxel>());
}

/*
 * the Reaction between two molecules
 */
std::pair<bool, Reaction<Voxel> > LatticeSimulator::apply_reaction_(const ReactionRule& reaction_rule,
        LatticeWorld::particle_info& from_info, const LatticeWorld::particle_info& to_info)
{
    const ReactionRule::product_container_type& products(reaction_rule.products());
    Reaction<Voxel> reaction;
    reaction.rule = reaction_rule;
    reaction.reactants.push_back(Reaction<Voxel>::particle_type(
                from_info.second, Voxel(from_info.first,
                    world_->get_molecular_type(from_info.first)->species(), 0)));
    reaction.reactants.push_back(Reaction<Voxel>::particle_type(
                to_info.second, Voxel(to_info.first,
                    world_->get_molecular_type(to_info.first)->species(), 0)));

    if (products.size() == 0)
    {
        // Not test yet
        world_->remove_molecule(from_info.first);
        world_->remove_molecule(to_info.first);
        return std::pair<bool, Reaction<Voxel> >(true, reaction);
    }
    else if (products.size() == 1)
    {
        const LatticeWorld::private_coordinate_type coord(to_info.first);
        Species product_species(*(products.begin()));
        world_->remove_molecule(from_info.first);
        world_->remove_molecule(to_info.first);

        if (!world_->has_species(product_species))
            new_species_.push_back(product_species);
        std::pair<ParticleID, bool> new_mol(world_->place_voxel_private(product_species, coord));
        if (new_mol.second)
        {
            return std::pair<bool, Reaction<Voxel> >(false, reaction);
        }

        reaction.products.push_back(Reaction<Voxel>::particle_type(
                    new_mol.first, Voxel(coord, product_species, 0)));

        reactions_.push_back(reaction_rule);

        return std::pair<bool, Reaction<Voxel> >(true, reaction);
    }
    else if (products.size() == 2)
    {
        // Not test yet
        const LatticeWorld::private_coordinate_type from_coord(from_info.first);
        const LatticeWorld::private_coordinate_type to_coord(to_info.first);
        world_->remove_molecule(from_coord);
        world_->remove_molecule(to_coord);

        const Species product_species0(*(products.begin())),
                      product_species1(*(++(products.begin())));

        if (!world_->has_species(product_species0))
            new_species_.push_back(product_species0);
        std::pair<ParticleID, bool> new0(world_->place_voxel_private(product_species0, from_coord));
        reaction.products.push_back(Reaction<Voxel>::particle_type(
                    new0.first, Voxel(from_coord, product_species0, 0)));

        if (!world_->has_species(product_species1))
            new_species_.push_back(product_species1);
        std::pair<ParticleID, bool> new1(world_->place_voxel_private(product_species1, to_coord));
        reaction.products.push_back(Reaction<Voxel>::particle_type(
                    new1.first, Voxel(to_coord, product_species1, 0)));

        reactions_.push_back(reaction_rule);
        return std::pair<bool, Reaction<Voxel> >(true, reaction);
    }
    return std::pair<bool, Reaction<Voxel> >(false, reaction);
}

/*
 * the First Order Reaction
 */
std::pair<bool, Reaction<Voxel> > LatticeSimulator::apply_reaction_(
        const ReactionRule& reaction_rule, LatticeWorld::particle_info& info)
{
    Reaction<Voxel> reaction;
    reaction.rule = reaction_rule;
    reaction.reactants.push_back(Reaction<Voxel>::particle_type(info.second,
                Voxel(info.first,world_->get_molecular_type(info.first)->species(), 0)));

    if (reaction_rule.products().size() == 0)
    {
        world_->remove_molecule(info.first);
        return std::pair<bool, Reaction<Voxel> >(true, reaction);
    }
    else if (reaction_rule.products().size() == 1)
    {
        const Species product(*(reaction_rule.products().begin()));
        if (!world_->has_species(product))
                new_species_.push_back(product);
        const LatticeWorld::molecule_info_type
            minfo(world_->get_molecule_info(product));
        world_->update_voxel_private(
            ecell4::Voxel(product, info.first, minfo.radius, minfo.D));

        reactions_.push_back(reaction_rule);

        reaction.products.push_back(std::pair<ParticleID, Voxel>(
                    info.second, Voxel(world_->private2coord(info.first), product, 0)));
        return std::pair<bool, Reaction<Voxel> >(true, reaction);
    }
    else if (reaction_rule.products().size() == 2)
    {
        const Species product_species0(*(reaction_rule.products().begin())),
                      product_species1(*(++(reaction_rule.products().begin())));
        const Integer rnd(world_->rng()->uniform_int(0,11));
        const LatticeWorld::private_coordinate_type coord(info.first);
        std::pair<LatticeWorld::coordinate_type, bool> retval(world_->move_to_neighbor(info, rnd));
        if (retval.second)
        {
            if (!world_->has_species(product_species0))
                new_species_.push_back(product_species0);
            std::pair<ParticleID, bool> new_val(world_->place_voxel_private(product_species0, coord));

            if (!world_->has_species(product_species1))
                new_species_.push_back(product_species1);
            const LatticeWorld::molecule_info_type
                minfo1(world_->get_molecule_info(product_species1));
            world_->update_voxel_private(
                ecell4::Voxel(product_species1, info.first, minfo1.radius, minfo1.D));

            reactions_.push_back(reaction_rule);

            reaction.products.push_back(std::pair<ParticleID, Voxel>(
                        new_val.first, Voxel(info.first, product_species0, 0)));
            reaction.products.push_back(std::pair<ParticleID, Voxel>(
                        info.second, Voxel(retval.first, product_species1, 0)));
            return std::pair<bool, Reaction<Voxel> >(true, reaction);
        }
    }
    return std::pair<bool, Reaction<Voxel> >(false, reaction);
}

void LatticeSimulator::register_step_event(const Species& species)
{
    const boost::shared_ptr<EventScheduler::Event> event(
            create_step_event(species, world_->t()));
    scheduler_.add(event);
}

void LatticeSimulator::step()
{
    if (!is_initialized_)
    {
        initialize();
    }
    step_();
}

bool LatticeSimulator::step(const Real& upto)
{
    if (!is_initialized_)
    {
        initialize();
    }

    if (upto < t())
    {
        return false;
    }

    if (upto >= scheduler_.top().second->time())
    {
        step_();
        return true;
    }

    world_->set_t(upto);

    return false;
}

void LatticeSimulator::step_()
{
    reactions_.clear();
    new_species_.clear();

    EventScheduler::value_type top(scheduler_.pop());
    const Real time(top.second->time());
    top.second->fire(); // top.second->time_ is updated in fire()
    world_->set_t(time);
    EventScheduler::events_range events(scheduler_.events());
    for (EventScheduler::events_range::iterator itr(events.begin());
            itr != events.end(); ++itr)
    {
        (*itr).second->interrupt(time);
        scheduler_.update(*itr);
    }
    scheduler_.add(top.second);

    for (std::vector<Species>::const_iterator itr(new_species_.begin());
            itr != new_species_.end(); ++itr)
    {
        register_step_event(*itr);
    }
}

void LatticeSimulator::walk(const Species& species)
{
    boost::shared_ptr<GSLRandomNumberGenerator> rng(world_->rng());

    MolecularTypeBase* mt(world_->find_molecular_type(species));
    mt->shuffle(*rng);
    //std::set<ParticleID> pids;
    std::vector<ParticleID> pids;
    Integer i(0);
    Integer max(mt->size());
    pids.reserve(max);
    while(i < max)
    {
        LatticeWorld::particle_info& info(mt->at(i));
        //pids.insert(info.second);
        pids.push_back(info.second);
        const Integer nrnd(rng->uniform_int(0,11));
        std::pair<LatticeWorld::private_coordinate_type, bool> move(world_->move_to_neighbor(info, nrnd));
        if (!move.second)
        {
            const std::pair<bool, Reaction<Voxel> > retval(attempt_reaction_(info, move.first));
            if (retval.first)
            {
                const Reaction<Voxel> reaction(retval.second);
                for (std::vector<Reaction<Voxel>::particle_type>::const_iterator itr(reaction.reactants.begin());
                        itr != reaction.reactants.end(); ++itr)
                {
                    for (std::vector<ParticleID>::const_iterator pid_itr(pids.begin());
                            pid_itr != pids.end(); ++pid_itr)
                    {
                        if (*pid_itr == (*itr).first)
                        {
                            --i;
                            break;
                        }
                    }
                    /*
                    if (pids.find((*itr).first) != pids.end())
                    {
                        --i;
                    }
                    */
                }
            }
        }
        ++i;
        if (max > mt->size())
            max = mt->size();
    }
}

} // lattice

} // ecell4
