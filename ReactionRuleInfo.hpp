#ifndef REACTION_RULE_INFO_HPP
#define REACTION_RULE_INFO_HPP

#include <algorithm>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include "twofold_container.hpp"

template<typename Tid_, typename Tsid_, typename Trate_>
class ReactionRuleInfo
{
public:
    typedef Tsid_ species_id_type;

private:
    typedef std::vector<species_id_type> species_id_vector;

public:
    typedef species_id_vector species_id_range;
    typedef Tid_ identifier_type;
    typedef Trate_ rate_type;

    identifier_type const& id() const
    {
        return id_;
    }

    species_id_range const& get_products() const
    {
        return products_;
    }

    twofold_container<species_id_type> const& get_reactants() const
    {
        return reactants_;
    }

    rate_type k() const
    {
        return k_;
    }

    template<typename Tr1_, typename Tr2_>
    ReactionRuleInfo(identifier_type const& id, rate_type const& k,
            Tr1_ const& reactants, Tr2_ const& products)
        : id_(id), k_(k)
    {
        std::copy(boost::begin(reactants),
                boost::end(reactants),
                std::back_inserter(reactants_));
        std::copy(boost::begin(products),
                boost::end(products),
                std::back_inserter(products_));
    }

    ReactionRuleInfo(): id_(), k_(), reactants_(), products_() {}

private:
    identifier_type id_;
    rate_type k_;
    twofold_container<species_id_type> reactants_;
    species_id_vector products_;
};

#endif /* REACTION_RULE_INFO_HPP */
