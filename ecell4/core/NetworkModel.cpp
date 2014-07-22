#include <algorithm>

#include "exceptions.hpp"
#include "NetworkModel.hpp"


namespace ecell4
{

std::vector<ReactionRule> NetworkModel::query_reaction_rules(
    const Species& sp) const
{
    first_order_reaction_rules_map_type::const_iterator
        i(first_order_reaction_rules_map_.find(sp.serial()));

    std::vector<ReactionRule> retval;
    if (i != first_order_reaction_rules_map_.end())
    {
        retval.reserve((*i).second.size());
        for (first_order_reaction_rules_map_type::mapped_type::const_iterator
                 j((*i).second.begin()); j != (*i).second.end(); ++j)
        {
            retval.push_back(reaction_rules_[*j]);
        }
    }
    return retval;
}

std::vector<ReactionRule> NetworkModel::query_reaction_rules(
    const Species& sp1, const Species& sp2) const
{
    std::vector<ReactionRule> retval;
    const std::pair<Species::serial_type, Species::serial_type>
        key(sp1.serial() < sp2.serial()?
            std::make_pair(sp1.serial(), sp2.serial()):
            std::make_pair(sp2.serial(), sp1.serial()));

    second_order_reaction_rules_map_type::const_iterator
        i(second_order_reaction_rules_map_.find(key));
    if (i != second_order_reaction_rules_map_.end())
    {
        retval.reserve((*i).second.size());
        for (second_order_reaction_rules_map_type::mapped_type::const_iterator
                 j((*i).second.begin()); j != (*i).second.end(); ++j)
        {
            retval.push_back(reaction_rules_[*j]);
        }
    }
    return retval;
}

void NetworkModel::initialize()
{
    if (!dirty_)
    {
        return; // do nothing
    }

    species_cache_.clear();
    for (reaction_rule_container_type::const_iterator
        i(reaction_rules_.begin()); i != reaction_rules_.end(); ++i)
    {
        const ReactionRule::reactant_container_type&
            reactants((*i).reactants());
        const ReactionRule::product_container_type&
            products((*i).products());
        std::copy(reactants.begin(), reactants.end(),
                  std::back_inserter(species_cache_));
        std::copy(products.begin(), products.end(),
                  std::back_inserter(species_cache_));
    }
    std::sort(species_cache_.begin(), species_cache_.end());
    species_cache_.erase(
        std::unique(species_cache_.begin(), species_cache_.end()),
        species_cache_.end());

    dirty_ = false;
}

void NetworkModel::add_species_attribute(const Species& sp)
{
    if (has_species_attribute(sp))
    {
        throw AlreadyExists("species already exists");
    }
    species_attributes_.push_back(sp);

    dirty_ = true;
}

void NetworkModel::remove_species_attribute(const Species& sp)
{
    species_container_type::iterator i(
        std::find(species_attributes_.begin(), species_attributes_.end(), sp));
    if (i == species_attributes_.end())
    {
        std::ostringstream message;
        message << "Speices [" << sp.serial() << "] not found";
        throw NotFound(message.str()); // use boost::format if it's allowed
    }
    species_attributes_.erase(i);

    dirty_ = true;
}

bool NetworkModel::has_species_attribute(const Species& sp) const
{
    species_container_type::const_iterator i(
        std::find(species_attributes_.begin(), species_attributes_.end(), sp));
    return (i != species_attributes_.end());
}

void NetworkModel::add_reaction_rule(const ReactionRule& rr)
{
    reaction_rule_container_type::const_iterator
        i(std::find(reaction_rules_.begin(), reaction_rules_.end(), rr));
    if (i != reaction_rules_.end())
    {
        throw AlreadyExists("reaction rule already exists");
    }

    const reaction_rule_container_type::size_type idx(reaction_rules_.size());
    reaction_rules_.push_back(rr);

    if (rr.reactants().size() == 1)
    {
        first_order_reaction_rules_map_[rr.reactants()[0].serial()].push_back(idx);
    }
    else if (rr.reactants().size() == 2)
    {
        const Species::serial_type
            serial1(rr.reactants()[0].serial()),
            serial2(rr.reactants()[1].serial());
        const std::pair<Species::serial_type, Species::serial_type>
            key(serial1 < serial2?
                std::make_pair(serial1, serial2):
                std::make_pair(serial2, serial1));
        second_order_reaction_rules_map_[key].push_back(idx);
    }
    else
    {
        ;
    }

    dirty_ = true;
}

void NetworkModel::remove_reaction_rule(const ReactionRule& rr)
{
    // reaction_rule_container_type::iterator
    //     i(std::find(reaction_rules_.begin(), reaction_rules_.end(), rr));
    // if (i == reaction_rules_.end())
    // {
    //     throw NotFound("reaction rule not found");
    // }

    // reaction_rule_container_type::size_type const
    //     idx(i - reaction_rules_.begin()), last_idx(reaction_rules_.size() - 1);
    // reaction_rules_map_type::iterator
    //     j(reaction_rules_map_.find(rr.reactants()));
    // if (j == reaction_rules_map_.end())
    // {
    //     throw IllegalState("no corresponding map key found");
    // }
    // else if ((*j).second.erase(idx) == 0)
    // {
    //     throw IllegalState("no corresponding map value found");
    // }

    // if (idx < last_idx)
    // {
    //     reaction_rule_container_type::value_type const
    //         last_value(reaction_rules_[last_idx]);
    //     (*i) = last_value;
    //     j = reaction_rules_map_.find(last_value.reactants());
    //     if (j == reaction_rules_map_.end())
    //     {
    //         throw IllegalState("no corresponding map key for the last found");
    //     }
    //     else if ((*j).second.erase(last_idx) == 0)
    //     {
    //         throw IllegalState("no corresponding map value for the last found");
    //     }
    //     (*j).second.insert(idx);
    // }

    // reaction_rules_.pop_back();
    dirty_ = true;
}

bool NetworkModel::has_reaction_rule(const ReactionRule& rr) const
{
    reaction_rule_container_type::const_iterator
        i(std::find(reaction_rules_.begin(), reaction_rules_.end(), rr));
    return (i != reaction_rules_.end());
}

} // ecell4
