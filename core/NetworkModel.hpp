#ifndef __ECELL4_NETWORK_MODEL_HPP
#define __ECELL4_NETWORK_MODEL_HPP

// #include "get_mapper_mf.hpp"

#include <map>
#include <set>

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

protected:

    typedef std::map<ReactionRule::reactant_container_type,
                     std::set<reaction_rule_container_type::size_type> >
    reaction_rules_map_type;
    // typedef utils::get_mapper_mf<
    // ReactionRule::reactant_container_type,
    // std::set<reaction_rule_container_type::size_type> >::type
    // reaction_rules_map_type;

public:

    NetworkModel()
        : species_(), reaction_rules_()
    {
        ;
    }

    // ModelTraits

    std::vector<ReactionRule> query_reaction_rules(const Species& sp) const;
    std::vector<ReactionRule> query_reaction_rules(
        const Species& sp1, const Species& sp2) const;

    // NetworkModelTraits

    void add_species(const Species& sp);
    bool has_species(const Species& sp) const;
    void remove_species(const Species& sp);

    void add_reaction_rule(const ReactionRule& rr);
    void remove_reaction_rule(const ReactionRule& rr);
    bool has_reaction_rule(const ReactionRule& rr) const;

    // Optional functions

    const species_container_type& species() const
    {
        return species_;
    }

    const Species& species(const Species::serial_type& key) const
    {
        for (species_container_type::const_iterator i(species_.begin());
             i != species_.end(); ++i)
        {
            if ((*i).serial() == key)
            {
                return (*i);
            }
        }

        std::ostringstream message;
        message << "Speices [" << key << "] not found";
        throw NotFound(message.str()); // use boost::format if it's allowed
    }

    const reaction_rule_container_type& reaction_rules() const
    {
        return reaction_rules_;
    }

protected:

    species_container_type species_;
    reaction_rule_container_type reaction_rules_;
    reaction_rules_map_type reaction_rules_map_;
};

} // ecell4

#endif /* __ECELL4_NETWORK_MODEL_HPP */
