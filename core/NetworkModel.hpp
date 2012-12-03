#ifndef __NETWORK_MODEL_HPP
#define __NETWORK_MODEL_HPP

#include <set>

#include "get_mapper_mf.hpp"
#include "types.hpp"
#include "Species.hpp"
#include "ReactionRule.hpp"
#include "Model.hpp"


namespace ecell4
{

class NetworkModel
    : public Model
{
public:

    typedef std::vector<Species> species_container_type;
    typedef std::vector<ReactionRule> reaction_rule_container_type;

    NetworkModel()
        : species_(), reaction_rules_()
    {
        ;
    }

    species_container_type const& species() const
    {
        return species_;
    }

    reaction_rule_container_type const& reaction_rules() const
    {
        return reaction_rules_;
    }

    std::vector<ReactionRule> query_reaction_rules(Species const& sp) const;
    std::vector<ReactionRule> query_reaction_rules(
        Species const& sp1, Species const& sp2) const;

    void add_species(Species const& sp);
    void remove_species(Species const& sp);
    bool has_species(Species const& sp) const;

    void add_reaction_rule(ReactionRule const& rr);
    void remove_reaction_rule(ReactionRule const& rr);
    bool has_reaction_rule(ReactionRule const& rr) const;

protected:

    typedef utils::get_mapper_mf<
        ReactionRule::reactant_container_type,
        std::set<reaction_rule_container_type::size_type> >::type
    reaction_rules_map_type;

    species_container_type species_;
    reaction_rule_container_type reaction_rules_;
    reaction_rules_map_type reaction_rules_map_;
};

} // ecell4

#endif /* __NETWORK_MODEL_HPP */
