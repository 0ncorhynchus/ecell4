#ifndef __MODEL_HPP
#define __MODEL_HPP

#include "types.hpp"
#include "Species.hpp"
#include "ReactionRule.hpp"
#include "exceptions.hpp"


namespace ecell4
{

class Model
{
public:

    /**
     * a fundamental function to query unimolecular reaction rules from a reactant.
     * this must be overloaded by any sub classes of Model.
     * @param species Species of a reactant
     * @return the vector of ReactionRule(s)
     */
    virtual ReactionRuleVector query_reaction_rules(
        Species const& sp) const = 0;

    /**
     * a fundamental function to query bimolecular reaction rules from reactants.
     * this must be overloaded by any sub classes of Model.
     * @param species1 Species of the first reactant
     * @param species2 Species of the second reactant
     * @return the vector of ReactionRule(s)
     */
    virtual ReactionRuleVector query_reaction_rules(
        Species const& sp1, Species const& sp2) const = 0;

    /**
     * add a concrete species to the model.
     * this function is a part of the trait of NetworkModel.
     * @param species a new Species
     * @return if the species is not registered yet.
     */
    virtual bool add_species(Species const& sp)
    {
        throw NotSupported("add_species is not supported in this model class");
    }

    /**
     * add a concrete reaction rule to the model.
     * this function is a part of the trait of NetworkModel.
     * @param rr a new ReactionRule
     * @return if the reaction rule is not registered yet.
     */
    virtual bool add_reaction_rule(ReactionRule const& rr)
    {
        throw NotSupported(
            "add_reaction_rule is not supported in this model class");
    }
};

} // ecell4

#endif /* __MODEL_HPP */
