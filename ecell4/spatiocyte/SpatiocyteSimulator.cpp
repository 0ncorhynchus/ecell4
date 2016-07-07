#include "SpatiocyteSimulator.hpp"

#include <algorithm>
#include <iterator>
#include <ecell4/core/StructureType.hpp>

namespace ecell4
{

namespace spatiocyte
{

void SpatiocyteSimulator::initialize()
{
    scheduler_.clear();
    update_alpha_map();
    const std::vector<Species> species(world_->list_species());
    for (std::vector<Species>::const_iterator itr(species.begin());
        itr != species.end(); ++itr)
    {
        register_events(*itr);
    }


    const std::vector<ReactionRule>& rules(model_->reaction_rules());
    for (std::vector<ReactionRule>::const_iterator i(rules.begin());
        i != rules.end(); ++i)
    {
        const ReactionRule& rr(*i);
        if (rr.reactants().size() != 0)
        {
            continue;
        }
        const boost::shared_ptr<SpatiocyteEvent>
            zeroth_order_reaction_event(
                create_zeroth_order_reaction_event(rr, world_->t()));
        scheduler_.add(zeroth_order_reaction_event);
    }

    dt_ = scheduler_.next_time() - t();
}

void SpatiocyteSimulator::update_alpha_map()
{
    boost::shared_ptr<Model> model_(model());
    if (!model_ || !model_->is_static())
        return;

    const Model::reaction_rule_container_type reaction_rules(model_->reaction_rules());
    for (Model::reaction_rule_container_type::const_iterator itr(reaction_rules.begin());
            itr != reaction_rules.end(); ++itr)
    {
        const ReactionRule::reactant_container_type& reactants((*itr).reactants());
        if (reactants.size() != 2)
            continue;

        const Real alpha(calculate_alpha(*itr));
        for (int i(0); i < 2; ++i) {
            const Species& sp(reactants.at(i));
            alpha_map_type::iterator map_itr(alpha_map_.find(sp));
            if (map_itr == alpha_map_.end())
                alpha_map_.insert(alpha_map_type::value_type(sp, alpha));
            else if ((*map_itr).second > alpha)
                (*map_itr).second = alpha;
        }
    }
}

void SpatiocyteSimulator::register_events(const Species& sp)
{
    if (world_->has_molecule_pool(sp))
    {
        //TODO: Call steps only if sp is assigned not to StructureType.
        const boost::shared_ptr<SpatiocyteEvent> step_event(
                create_step_event(sp, world_->t()));
        scheduler_.add(step_event);
    }

    std::vector<ReactionRule> reaction_rules(model_->query_reaction_rules(sp));
    for (std::vector<ReactionRule>::const_iterator i(reaction_rules.begin());
        i != reaction_rules.end(); ++i)
    {
        const ReactionRule& rr(*i);
        const boost::shared_ptr<SpatiocyteEvent>
            first_order_reaction_event(
                create_first_order_reaction_event(rr, world_->t()));
        scheduler_.add(first_order_reaction_event);
    }
}

boost::shared_ptr<SpatiocyteEvent> SpatiocyteSimulator::create_step_event(
        const Species& species, const Real& t)
{
    double alpha(alpha_);
    alpha_map_type::const_iterator itr(alpha_map_.find(species));
    if (itr != alpha_map_.end() && (*itr).second < alpha)
        alpha = (*itr).second;

    boost::shared_ptr<SpatiocyteEvent> event(
        new StepEvent(this, species, t, alpha));
    return event;
}

boost::shared_ptr<SpatiocyteEvent>
SpatiocyteSimulator::create_zeroth_order_reaction_event(
    const ReactionRule& reaction_rule, const Real& t)
{
    boost::shared_ptr<SpatiocyteEvent> event(
            new ZerothOrderReactionEvent(world_, reaction_rule, t));
    return event;
}

boost::shared_ptr<SpatiocyteEvent>
SpatiocyteSimulator::create_first_order_reaction_event(
    const ReactionRule& reaction_rule, const Real& t)
{
    boost::shared_ptr<SpatiocyteEvent> event(new FirstOrderReactionEvent(
                world_, reaction_rule, t));
    return event;
}

void SpatiocyteSimulator::finalize()
{
    scheduler_type::events_range events(scheduler_.events());
    for (scheduler_type::events_range::iterator itr(events.begin());
            itr != events.end(); ++itr)
    {
        const Real queued_time((*itr).second->time() - (*itr).second->dt());
        StepEvent* step_event(dynamic_cast<StepEvent*>((*itr).second.get()));
        if (step_event != NULL && queued_time < t())
        {
            const Real alpha((t() - queued_time) / (*itr).second->dt());
            step_event->walk(alpha);
        }
    }
    initialize();
}

Real SpatiocyteSimulator::calculate_dimensional_factor(
    const VoxelPool* mt0, const VoxelPool* mt1) const
{
    const Species&
        speciesA(mt0->species()),
        speciesB(mt1->species());
    const Real
        D_A(mt0->D()),
        D_B(mt1->D());
    const Shape::dimension_kind
        dimensionA(mt0->get_dimension()),
        dimensionB(mt1->get_dimension());
    const Real Dtot(D_A + D_B);
    const Real gamma(pow(2 * sqrt(2.0) + 4 * sqrt(3.0) + 3 * sqrt(6.0) + sqrt(22.0), 2) /
        (72 * (6 * sqrt(2.0) + 4 * sqrt(3.0) + 3 * sqrt(6.0))));
    Real factor(0);
    if (dimensionA == Shape::THREE && dimensionB == Shape::THREE)
    {
        // if (speciesA != speciesB)
        //     factor = 1. / (6 * sqrt(2.0) * Dtot * world_->voxel_radius());
        // else
        //     factor = 1. / (6 * sqrt(2.0) * D_A * world_->voxel_radius());
        factor = 1. / (6 * sqrt(2.0) * Dtot * world_->voxel_radius());
    }
    else if (dimensionA == Shape::TWO && dimensionB == Shape::TWO)
    {
        // if (speciesA != speciesB)
        //     factor = gamma / Dtot;
        // else
        //     factor = gamma / D_A;
        factor = gamma / Dtot;
    }
    else if (dimensionA == Shape::THREE && dimensionB == Shape::TWO)
    {
        factor = sqrt(2.0) / (3 * D_A * world_->voxel_radius());
        if (mt1->is_structure()) // B is Surface
        {
            factor *= world_->unit_area();
        }
    }
    else if (dimensionA == Shape::TWO && dimensionB == Shape::THREE)
    {
        factor = sqrt(2.0) / (3 * D_B * world_->voxel_radius());
        if (mt0->is_structure()) // A is Surface
        {
            factor *= world_->unit_area();
        }
    }
    else
        throw NotSupported("The dimension of a structure must be two or three.");
    return factor;
}

Real SpatiocyteSimulator::calculate_alpha(const ReactionRule& rule) const
{
    const ReactionRule::reactant_container_type& reactants(rule.reactants());
    if (reactants.size() != 2)
        return 1.0;

    const Species species[2] = {reactants.at(0), reactants.at(1)};
    const MoleculeInfo info[2] = {
        world_->get_molecule_info(species[0]),
        world_->get_molecule_info(species[1])
    };
    VoxelPool* mt[2];
    bool is_created[2] = {false, false};
    for (int i(0); i < 2; ++i) {
        try
        {
            mt[i] = world_->find_voxel_pool(species[i]);
        }
        catch(NotFound e)
        {
            VoxelPool *location(&(VacantType::getInstance()));
            if (info[i].loc != "") {
                try
                {
                    location = world_->find_voxel_pool(Species(info[i].loc));
                }
                catch(NotFound e)
                {
                    ;
                }
            }
            mt[i] = new MolecularType(species[i], location, info[i].radius, info[i].D);
            is_created[i] = true;
        }
    }
    const Real factor(calculate_dimensional_factor(mt[0], mt[1]));
    for (int i(0); i < 2; ++i)
        if (is_created[i])
            delete mt[i];
    const Real alpha(1.0 / (factor * rule.k()));
    return alpha < 1.0 ? alpha : 1.0;
}

std::pair<SpatiocyteSimulator::attempt_reaction_result_type, SpatiocyteSimulator::reaction_type> SpatiocyteSimulator::attempt_reaction_(
    const SpatiocyteWorld::coordinate_id_pair_type& info, SpatiocyteWorld::coordinate_type to_coord,
    const Real& alpha)
{
    const VoxelPool* from_mt(
        world_->find_voxel_pool(info.coordinate));
    const VoxelPool* to_mt(
        world_->find_voxel_pool(to_coord));

    if (to_mt->is_vacant())
    {
        return std::make_pair(NO_REACTION, reaction_type());
    }

    const Species&
        speciesA(from_mt->species()),
        speciesB(to_mt->species());

    const std::vector<ReactionRule> rules(
        model_->query_reaction_rules(speciesA, speciesB));

    if (rules.empty())
    {
        return std::make_pair(NO_REACTION, reaction_type());
    }

    const Real factor(calculate_dimensional_factor(from_mt, to_mt));

    const Real rnd(world_->rng()->uniform(0,1));
    Real accp(0.0);
    for (std::vector<ReactionRule>::const_iterator itr(rules.begin());
        itr != rules.end(); ++itr)
    {
        const Real k((*itr).k());
        const Real P(k * factor * alpha);
        accp += P;
        if (accp > 1)
        {
            std::cerr << "The total acceptance probability [" << accp
                << "] exceeds 1 for '" << speciesA.serial()
                << "' and '" << speciesB.serial() << "'." << std::endl;
        }
        if (accp >= rnd)
        {
            std::pair<bool, SpatiocyteSimulator::reaction_type>
                retval = apply_second_order_reaction_(
                    *itr,
                    world_->make_pid_voxel_pair(from_mt, info),
                    world_->make_pid_voxel_pair(to_mt, to_coord));
            if (retval.first)
            {
                return std::make_pair(REACTION_SUCCEEDED, retval.second);
            }
            else
            {
                return std::make_pair(REACTION_FAILED, reaction_type());
            }
        }
    }
    return std::make_pair(REACTION_FAILED, reaction_type());
}

/*
 * the Reaction between two molecules
 */
std::pair<bool, SpatiocyteSimulator::reaction_type> SpatiocyteSimulator::apply_second_order_reaction_(
    const ReactionRule& reaction_rule,
    const SpatiocyteSimulator::reaction_info_type::particle_id_pair_type& p0,
    const SpatiocyteSimulator::reaction_info_type::particle_id_pair_type& p1)
{
    const ReactionRule::product_container_type&
        products(reaction_rule.products());

    std::pair<bool, reaction_type> retval;

    switch (products.size())
    {
        case 0:
            retval = apply_vanishment(reaction_rule, p0, p1);
            assert(retval.first);
            break;
        case 1:
            retval = apply_ab2c(reaction_rule, p0, p1, *(products.begin()));
            break;
        case 2:
            retval = apply_ab2cd(
                reaction_rule, p0, p1, *(products.begin()), *(++(products.begin())));
            break;
        default:
            return std::make_pair(false, reaction_type());
    }

    if (retval.first)
    {
        last_reactions_.push_back(retval.second);
    }
    return retval;
}

std::pair<bool, SpatiocyteSimulator::reaction_type> SpatiocyteSimulator::apply_vanishment(
    const ReactionRule& reaction_rule,
    const SpatiocyteSimulator::reaction_info_type::particle_id_pair_type& p0,
    const SpatiocyteSimulator::reaction_info_type::particle_id_pair_type& p1)
{
    reaction_info_type rinfo(world_->t());
    rinfo.add_reactant(p0);
    rinfo.add_reactant(p1);

    world_->remove_voxel(p0.second.coordinate());
    world_->remove_voxel(p1.second.coordinate());

    return std::make_pair(true, std::make_pair(reaction_rule, rinfo));
}

std::pair<bool, SpatiocyteSimulator::reaction_type> SpatiocyteSimulator::apply_ab2c(
    const ReactionRule& reaction_rule,
    const SpatiocyteSimulator::reaction_info_type::particle_id_pair_type& p0,
    const SpatiocyteSimulator::reaction_info_type::particle_id_pair_type& p1,
    const Species& product_species)
{
    // A and B (from_info and to_info) become C (product_species)
    const std::string location(world_->get_molecule_info(product_species).loc);
    const std::string fserial(get_serial(p0.second.coordinate()));
    const std::string floc(get_location(p0.second.coordinate()));
    const std::string tserial(get_serial(p1.second.coordinate()));
    const std::string tloc(get_location(p1.second.coordinate()));

    reaction_info_type rinfo(world_->t());

    if (tserial == location || tloc == location)
    {
        // B is on the location of C, or the location itself.
        // Place C at the coordinate of B, and remove A.
        rinfo.add_reactant(p0);
        rinfo.add_reactant(p1);

        if (tserial != location)
        {
            world_->remove_voxel(p1.second.coordinate());
        }

        world_->remove_voxel(p0.second.coordinate());
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
            world_->new_voxel(product_species, p1.second.coordinate()));

        rinfo.add_product(new_mol.first);
    }
    else if (fserial == location || floc == location)
    {
        // A is on the location of C, or the location itself.
        // Place C at the coordinate of A, and remove B.
        rinfo.add_reactant(p0);
        rinfo.add_reactant(p1);

        if (fserial != location)
        {
            world_->remove_voxel(p0.second.coordinate());
        }

        world_->remove_voxel(p1.second.coordinate());
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
            world_->new_voxel(product_species, p0.second.coordinate()));

        rinfo.add_product(new_mol.first);
    }
    else
    {
        return std::make_pair(false, reaction_type());
        // throw IllegalState(
        //     "no place for the product [" + product_species.serial() + "].");
    }
    return std::make_pair(true, std::make_pair(reaction_rule, rinfo));
}

// Not tested yet
std::pair<bool, SpatiocyteSimulator::reaction_type> SpatiocyteSimulator::apply_ab2cd(
    const ReactionRule& reaction_rule,
    const reaction_info_type::particle_id_pair_type& p0,
    const reaction_info_type::particle_id_pair_type& p1,
    const Species& product_species0,
    const Species& product_species1)
{
    const SpatiocyteWorld::coordinate_type from_coord(p0.second.coordinate());
    const SpatiocyteWorld::coordinate_type to_coord(p1.second.coordinate());
    const std::string aserial(get_serial(from_coord));
    const std::string aloc(get_location(from_coord));
    const std::string bserial(get_serial(to_coord));
    const std::string bloc(get_location(to_coord));
    const std::string cloc(world_->get_molecule_info(product_species0).loc);
    const std::string dloc(world_->get_molecule_info(product_species1).loc);

    if (aserial == cloc || aloc == cloc)
    {
        if (bserial == dloc || bloc == dloc)
        {
            if (aserial != cloc)
            {
                // Remove A once if A is not the location of C
                world_->remove_voxel(p0.second.coordinate());
            }
            if (bserial != dloc)
            {
                // Remove B once if B is not the location of D
                world_->remove_voxel(p1.second.coordinate());
            }
            return apply_ab2cd_in_order(
                reaction_rule, p0, p1, product_species0, product_species1,
                from_coord, to_coord);
        }
        else
        {
            std::pair<SpatiocyteWorld::coordinate_type, bool>
                neighbor(world_->check_neighbor(to_coord, dloc));

            if (neighbor.second)
            {
                world_->remove_voxel(p1.second.coordinate());
                if (aserial != cloc)
                {
                    // Remove A once if A is not the location of C
                    world_->remove_voxel(p0.second.coordinate());
                }
                return apply_ab2cd_in_order(
                    reaction_rule, p0, p1, product_species0, product_species1,
                    from_coord, neighbor.first);
            }
        }
    }
    else if (aserial == dloc || aloc == dloc)
    {
        if (bserial == cloc || bloc == dloc)
        {
            if (aserial != dloc)
            {
                // Remove A once if A is not the location of D
                world_->remove_voxel(p0.second.coordinate());
            }
            if (bserial != cloc)
            {
                // Remove B once if B is not the location of C
                world_->remove_voxel(p1.second.coordinate());
            }
            return apply_ab2cd_in_order(
                reaction_rule, p0, p1, product_species0, product_species1,
                to_coord, from_coord);
        }
        else
        {
            std::pair<SpatiocyteWorld::coordinate_type, bool>
                neighbor(world_->check_neighbor(to_coord, cloc));

            if (neighbor.second)
            {
                world_->remove_voxel(p1.second.coordinate());
                if (aserial != dloc)
                {
                    // Remove A once if A is not the location of D
                    world_->remove_voxel(p0.second.coordinate());
                }
                return apply_ab2cd_in_order(
                    reaction_rule, p0, p1, product_species0, product_species1,
                    neighbor.first, from_coord);
            }
        }
    }
    else if (bserial == cloc || bloc == cloc)
    {
        std::pair<SpatiocyteWorld::coordinate_type, bool>
            neighbor(world_->check_neighbor(to_coord, dloc));

        if (neighbor.second)
        {
            world_->remove_voxel(p0.second.coordinate());
            if (bserial != cloc)
            {
                // Remove B once if B is not the location of C
                world_->remove_voxel(p1.second.coordinate());
            }
            return apply_ab2cd_in_order(
                reaction_rule, p0, p1, product_species0, product_species1,
                to_coord, neighbor.first);
        }
    }
    else if (bserial == dloc || bloc == dloc)
    {
        std::pair<SpatiocyteWorld::coordinate_type, bool>
            neighbor(world_->check_neighbor(to_coord, dloc));

        if (neighbor.second)
        {
            world_->remove_voxel(p0.second.coordinate());
            if (bserial != dloc)
            {
                // Remove B once if B is not the location of D
                world_->remove_voxel(p1.second.coordinate());
            }
            return apply_ab2cd_in_order(
                reaction_rule, p0, p1, product_species0, product_species1,
                neighbor.first, to_coord);
        }
    }
    // else
    // {
    //     throw IllegalState("Not Supported.");
    // }
    return std::make_pair(false, reaction_type());
}

std::pair<bool, SpatiocyteSimulator::reaction_type> SpatiocyteSimulator::apply_ab2cd_in_order(
    const ReactionRule& reaction_rule,
    const SpatiocyteSimulator::reaction_info_type::particle_id_pair_type& p0,
    const SpatiocyteSimulator::reaction_info_type::particle_id_pair_type& p1,
    const Species& product_species0,
    const Species& product_species1,
    const SpatiocyteWorld::coordinate_type coord0,
    const SpatiocyteWorld::coordinate_type coord1)
{
    reaction_info_type rinfo(world_->t());
    rinfo.add_reactant(p0);
    rinfo.add_reactant(p1);

    std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
        world_->new_voxel(product_species0, coord0));
    if (!new_mol0.second)
    {
        throw IllegalState("no place for " + product_species0.serial());
    }
    std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
        world_->new_voxel(product_species1, coord1));
    if (!new_mol1.second)
    {
        throw IllegalState("no place for " + product_species1.serial());
    }

    rinfo.add_product(new_mol0.first);
    rinfo.add_product(new_mol1.first);
    return std::make_pair(true, std::make_pair(reaction_rule, rinfo));
}

void SpatiocyteSimulator::register_product_species(const Species& product_species)
{
    if (!world_->has_species(product_species))
    {
        new_species_.push_back(product_species);
    }
}

// void SpatiocyteSimulator::register_reactant_species(
//         const SpatiocyteWorld::coordinate_id_pair_type pinfo, reaction_type& reaction) const
// {
//     const VoxelPool* mtype(world_->find_voxel_pool(pinfo.first));
//     const std::string location(
//             mtype->location()->is_vacant() ? "" : mtype->location()->species().serial());
//     reaction.reactants.push_back(
//         reaction_type::particle_type(
//             pinfo.second,
//             Voxel(mtype->species(), pinfo.first,
//                   mtype->radius(), mtype->D(), location)));
// }

void SpatiocyteSimulator::step()
{
    step_();
    dt_ = scheduler_.next_time() - t();
}

bool SpatiocyteSimulator::step(const Real& upto)
{
    if (upto < t())
    {
        return false;
    }

    if (scheduler_.size() > 0 && upto >= scheduler_.top().second->time())
    {
        step_();
        dt_ = scheduler_.next_time() - t();
        return true;
    }

    world_->set_t(upto); //XXX: TODO
    last_reactions_.clear();
    new_species_.clear();
    dt_ = scheduler_.next_time() - t();
    return false;
}

void SpatiocyteSimulator::step_()
{
    last_reactions_.clear();
    new_species_.clear();

    scheduler_type::value_type top(scheduler_.pop());
    const Real time(top.second->time());
    world_->set_t(time);
    top.second->fire(); // top.second->time_ is updated in fire()
    if (!check_reaction())
        last_reactions_ = top.second->last_reactions();
    for (std::vector<reaction_type>::const_iterator itr(last_reactions_.begin());
            itr != last_reactions_.end(); ++itr)
        for (reaction_info_type::container_type::const_iterator
                product((*itr).second.products().begin());
                product != (*itr).second.products().end(); ++product)
            register_product_species((*product).second.species());

    scheduler_type::events_range events(scheduler_.events());
    for (scheduler_type::events_range::iterator itr(events.begin());
        itr != events.end(); ++itr)
    {
        (*itr).second->interrupt(time);
        scheduler_.update(*itr);
    }
    scheduler_.add(top.second);

    // update_alpha_map(); // may be performance cost
    for (std::vector<Species>::const_iterator itr(new_species_.begin());
        itr != new_species_.end(); ++itr)
    {
        register_events(*itr);
    }

    num_steps_++;
}


} // spatiocyte

} // ecell4
