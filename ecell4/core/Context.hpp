#ifndef ECELL4_CONTEXT_HPP
#define ECELL4_CONTEXT_HPP

#include "get_mapper_mf.hpp"
#include "Species.hpp"
#include "ReactionRule.hpp"
#include <boost/array.hpp>


namespace ecell4
{

namespace context
{

namespace rbex
{

inline bool is_wildcard(const std::string& name)
{
    return (name.size() > 0 && name[0] == '_');
}

inline bool is_unnamed_wildcard(const std::string& name)
{
    return name == "_";
}

inline bool is_pass_wildcard(const std::string& name)
{
    return name == "_0";
}

inline bool is_named_wildcard(const std::string& name)
{
    return (name.size() > 1 && name[0] == '_' && !is_pass_wildcard(name));
}

} // rbex

class MatchObject
{
public:

    typedef struct
    {
        typedef std::vector<Species::container_type::difference_type>
            iterator_container_type;
        typedef utils::get_mapper_mf<std::string, std::string>::type
            variable_container_type;

        iterator_container_type iterators;
        variable_container_type locals;
        variable_container_type globals;
    } context_type;

public:

    MatchObject(const UnitSpecies& pttrn)
        : pttrn_(pttrn)
    {
        ;
    }

    virtual ~MatchObject()
    {
        ;
    }

    std::pair<bool, context_type> match(
        const Species& sp, const context_type& ctx)
    {
        // target_ = sp;
        target_ = sp.units();
        itr_ = target_.begin();
        ctx_ = ctx;
        return next();
    }

    std::pair<bool, context_type> match(
        const std::vector<UnitSpecies>& target, const context_type& ctx)
    {
        target_ = target;
        itr_ = target_.begin();
        ctx_ = ctx;
        return next();
    }

    std::pair<bool, context_type> next();

protected:

    UnitSpecies pttrn_;
    // Species target_;
    std::vector<UnitSpecies> target_;
    Species::container_type::const_iterator itr_;
    context_type ctx_;
};

std::pair<bool, MatchObject::context_type>
uspmatch(const UnitSpecies& pttrn, const UnitSpecies& sp,
    const MatchObject::context_type& org);
bool __spmatch(
    Species::container_type::const_iterator itr,
    const Species::container_type::const_iterator& end,
    const Species& sp, const MatchObject::context_type& ctx);
bool spmatch(const Species& pttrn, const Species& sp);
Integer count_spmatches(const Species& pttrn, const Species& sp);
Integer count_spmatches(
    const Species& pttrn, const Species& sp,
    const MatchObject::context_type::variable_container_type& globals);

// ReactionRule create_reaction_rule_formatted(
//     const ReactionRule::reactant_container_type& reactants,
//     const ReactionRule::product_container_type& products, const Real k);

// inline ReactionRule create_reaction_rule_formatted(const ReactionRule& rr)
// {
//     return create_reaction_rule_formatted(rr.reactants(), rr.products(), rr.k());
// }

// class SpeciesExpressionMatcher
// {
// public:
// 
//     typedef MatchObject::context_type context_type;
// 
// public:
// 
//     SpeciesExpressionMatcher(const Species& pttrn)
//         : pttrn_(pttrn.units())
//     {
//         ;
//     }
// 
//     virtual ~SpeciesExpressionMatcher()
//     {
//         ;
//     }
// 
//     bool match(const Species& sp)
//     {
//         context_type::variable_container_type globals;
//         return match(sp, globals);
//     }
// 
//     bool match(
//         const Species& sp, const context_type::variable_container_type& globals)
//     {
//         matches_.clear();
//         for (Species::container_type::const_iterator i(pttrn_.begin());
//             i != pttrn_.end(); ++i)
//         {
//             matches_.push_back(MatchObject(*i));
//         }
// 
//         // target_ = sp;
//         target_ = sp.units();
//         itr_ = matches_.begin();
//         context_type ctx;
//         ctx.globals = globals;
//         return __match(ctx);
//     }
// 
//     bool __match(const context_type& ctx)
//     {
//         if (itr_ == matches_.end())
//         {
//             ctx_ = ctx;
//             return true;
//         }
// 
//         std::pair<bool, context_type> retval((*itr_).match(target_, ctx));
//         while (retval.first)
//         {
//             ++itr_;
//             const bool succeeded(__match(retval.second));
//             if (succeeded)
//             {
//                 return true;
//             }
//             --itr_;
//             retval = (*itr_).next();
//         }
//         return false;
//     }
// 
//     bool next()
//     {
//         if (itr_ != matches_.end())
//         {
//             return false;
//         }
//         else if (matches_.size() == 0)
//         {
//             return true;
//         }
// 
//         do
//         {
//             --itr_;
//             std::pair<bool, context_type> retval((*itr_).next());
//             while (retval.first)
//             {
//                 ++itr_;
//                 const bool succeeded(__match(retval.second));
//                 if (succeeded)
//                 {
//                     return true;
//                 }
//                 --itr_;
//                 retval = (*itr_).next();
//             }
//         }
//         while (itr_ != matches_.begin());
//         return false;
//     }
// 
//     Integer count(const Species& sp)
//     {
//         context_type::variable_container_type globals;
//         if (!match(sp, globals))
//         {
//             return 0;
//         }
//         Integer n(1);
//         while (next())
//         {
//             ++n;
//         }
//         return n;
//     }
// 
//     const context_type& context() const
//     {
//         return ctx_;
//     }
// 
// protected:
// 
//     // Species pttrn_;
//     // Species target_;
//     std::vector<UnitSpecies> pttrn_;
//     std::vector<UnitSpecies> target_;
//     std::vector<MatchObject> matches_;
//     std::vector<MatchObject>::iterator itr_;
//     context_type ctx_;
// };

class unit_species_comparerator
{
public:

    // typedef Species::container_type::size_type index_type;
    typedef unsigned int index_type;
    typedef std::pair<index_type, std::string> site_type;
    typedef utils::get_mapper_mf<std::string, std::vector<site_type> >::type
        connection_container_type;

public:

    unit_species_comparerator(const Species& sp)
        : root_(sp.units())
    {
        initialize();
    }

    const std::vector<UnitSpecies>& units() const
    {
        return root_;
    }

    void initialize()
    {
        connections_.clear();
        for (index_type idx(0); idx < root_.size(); ++idx)
        {
            const UnitSpecies usp(root_.at(idx));
            for (UnitSpecies::container_type::const_iterator i(usp.begin());
                 i != usp.end(); ++i)
            {
                if ((*i).second.second == "" || rbex::is_wildcard((*i).second.second))
                {
                    continue;
                }

                if (connections_.find((*i).second.second) == connections_.end())
                {
                    connections_.insert(std::make_pair(
                        (*i).second.second, std::vector<site_type>()));
                }
                connections_[(*i).second.second].push_back(
                    std::make_pair(idx, (*i).first));
            }
        }
    }

    int compare(const index_type& val1, const index_type& val2)
    {
        if (val1 == val2)
        {
            return 0;
        }

        const std::pair<index_type, index_type> pair_key((val1 < val2)?
            std::make_pair(val1, val2) : std::make_pair(val1, val2));
        if (std::binary_search(ignores_.begin(), ignores_.end(), pair_key))
        {
            return 0;
        }

        const UnitSpecies& lhs(root_.at(val1));
        const UnitSpecies& rhs(root_.at(val2));

        if (lhs.name() != rhs.name())
        {
            return (lhs.name() < rhs.name()? 1 : -1);
        }

        UnitSpecies::container_type::const_iterator
            i(lhs.begin()), j(rhs.begin());
        while (i != lhs.end() && j != rhs.end())
        {
            if ((*i).first != (*j).first)
            {
                // std::cout << "[1] " << lhs.serial() << "(" << val1 << ") vs "
                //     << rhs.serial() << "(" << val2 << ") -> " << (*i).first
                //     << " < " << (*j).first << std::endl;
                return ((*i).first < (*j).first? 1 : -1);
            }
            else if ((*i).second.first != (*j).second.first)
            {
                // std::cout << "[2] " << lhs.serial() << "(" << val1 << ") vs "
                //     << rhs.serial() << "(" << val2 << ")" << std::endl;
                return ((*i).second.first < (*j).second.first? 1 : -1);
            }
            else if (((*i).second.second == "") != ((*j).second.second == ""))
            {
                // std::cout << "[3] " << lhs.serial() << "(" << val1 << ") vs "
                //     << rhs.serial() << "(" << val2 << ") -> '"
                //     << (*i).second.second << "' < '" << (*j).second.second
                //     << "'" << std::endl;
                return ((*i).second.second == ""? 1 : -1);
            }

            ++i;
            ++j;
        }

        if (lhs.num_sites() != rhs.num_sites())
        {
            return (lhs.num_sites() < rhs.num_sites()? 1 : -1);
        }

        ignores_.insert(
            std::lower_bound(ignores_.begin(), ignores_.end(), pair_key),
            pair_key);
        i = lhs.begin();
        j = rhs.begin();
        while (i != lhs.end() && j != rhs.end())
        {
            if ((*i).second.second != "" && (*i).second.second != "")
            {
                const std::vector<site_type>&
                    pair1(connections_[(*i).second.second]);
                const std::vector<site_type>&
                    pair2(connections_[(*j).second.second]);
                const site_type& target1(
                    (pair1[0].first == val1 && pair1[0].second == (*i).first)?
                    pair1[1] : pair1[0]);
                const site_type& target2(
                    (pair2[0].first == val2 && pair2[0].second == (*j).first)?
                    pair2[1] : pair2[0]);
                if (target1.second != target2.second)
                {
                    ignores_.pop_back();
                    return (target1.second < target2.second? 1 : -1);
                }

                const int retval(compare(target1.first, target2.first));
                // std::cout << "[0] " << lhs.serial() << "(" << val1 << ") vs "
                //     << rhs.serial() << "(" << val2 << ") -> " << retval
                //     << std::endl;
                if (retval != 0)
                {
                    ignores_.pop_back();
                    return retval;
                }
            }

            ++i;
            ++j;
        }
        ignores_.pop_back();
        return 0;
    }

    bool operator()(const index_type& val1, const index_type& val2)
    {
        // return val1 < val2;
        ignores_.clear();
        return 0 < compare(val1, val2);
    }

    void reorder_units(
        std::vector<unsigned int>& unit_indices, const unsigned int& idx,
        unsigned int& stride)
    {
        if (unit_indices[idx] != root_.size())
        {
            return;
        }

        const UnitSpecies& usp(root_.at(idx));

        unit_indices[idx] = stride;
        ++stride;

        for (UnitSpecies::container_type::const_iterator i(usp.begin());
            i != usp.end(); ++i)
        {
            if ((*i).second.second == "" || rbex::is_wildcard((*i).second.second))
            {
                continue;
            }

            // const std::vector<unit_species_comparerator::site_type>&
            //     pair((*connections_.find((*i).second.second)).second);
            const std::vector<unit_species_comparerator::site_type>&
                pair(connections_[(*i).second.second]);
            const unit_species_comparerator::site_type&
                tgt((pair[0].first == idx && pair[0].second == (*i).first)?
                    pair[1] : pair[0]);

            reorder_units(unit_indices, tgt.first, stride);
        }
    }

protected:

    const std::vector<UnitSpecies> root_;
    connection_container_type connections_;
    std::vector<std::pair<index_type, index_type> > ignores_;
};

Species format_species(const Species& sp);

inline Species::serial_type unique_serial(const Species& sp)
{
    return format_species(sp).serial();
}

}  // context

namespace _context
{

template <typename T>
class rule_based_expression_matcher {};

template <>
class rule_based_expression_matcher<std::vector<UnitSpecies> >
{
public:

    typedef std::vector<UnitSpecies> value_type;
    typedef ecell4::context::MatchObject submatcher_type;

    typedef submatcher_type::context_type context_type;

public:

    rule_based_expression_matcher(const value_type& pttrn)
        : pttrn_(pttrn)
    {
        ;
    }

    bool match(const value_type& another)
    {
        context_type::variable_container_type globals;
        return match(another, globals);
    }

    bool match(
        const value_type& another, const context_type::variable_container_type& globals)
    {
        matches_.clear();
        for (value_type::const_iterator i(pttrn_.begin());
            i != pttrn_.end(); ++i)
        {
            matches_.push_back(submatcher_type(*i));
        }

        target_ = another;
        context_type ctx;
        ctx.globals = globals;
        return advance(0, ctx);
    }

    bool next()
    {
        if (matches_.size() == 0)
        {
            return true;
        }

        return __next();
    }

    Integer count(const value_type& another)
    {
        context_type::variable_container_type globals;
        if (!match(another, globals))
        {
            return 0;
        }

        Integer n(1);
        while (next())
        {
            ++n;
        }
        return n;
    }

    const context_type& context() const
    {
        return ctx_;
    }

protected:

    bool __next()
    {
        //XXX: Make sure match was already called..
        size_t pos = matches_.size();

        do
        {
            --pos;

            submatcher_type& matcher_at_pos(matches_[pos]);
            std::pair<bool, context_type> retval(matcher_at_pos.next());

            if (retval.first)
            {
                do
                {
                    if (this->advance(pos + 1, retval.second))
                    {
                        return true;
                    }
                    retval = matcher_at_pos.next();
                } while (retval.first);
            }
        }
        while (pos > 0);

        return false;
    }

    bool advance(const size_t pos, const context_type& ctx)
    {
        if (pos == matches_.size())
        {
            ctx_ = ctx;
            return true;
        }

        submatcher_type& matcher_at_pos(matches_[pos]);
        std::pair<bool, context_type> retval(matcher_at_pos.match(target_, ctx));

        if (!retval.first)
        {
            return false;
        }

        do
        {
            if (this->advance(pos + 1, retval.second))
            {
                return true;
            }
            retval = matcher_at_pos.next();
        } while (retval.first);

        return false;
    }

protected:

    value_type pttrn_;
    value_type target_;
    std::vector<submatcher_type> matches_;
    std::vector<submatcher_type>::iterator itr_;
    context_type ctx_;
};

template <>
class rule_based_expression_matcher<Species>
{
public:

    typedef Species value_type;
    typedef rule_based_expression_matcher<std::vector<UnitSpecies> > base_type;
    typedef base_type::context_type context_type;

    rule_based_expression_matcher(const value_type& pttrn)
        : base_(pttrn.units())
    {
        ;
    }

    bool match(const value_type& another)
    {
        return base_.match(another.units());
    }

    bool match(
        const value_type& another, const context_type::variable_container_type& globals)
    {
        return base_.match(another.units(), globals);
    }

    bool next()
    {
        return base_.next();
    }

    size_t count(const value_type& another)
    {
        return base_.count(another.units());
    }

    const context_type& context() const
    {
        return base_.context();
    }

protected:

    base_type base_;
};

template <>
class rule_based_expression_matcher<std::vector<Species> >
{
public:

    typedef std::vector<Species> value_type;
    typedef rule_based_expression_matcher<Species> submatcher_type;
    typedef submatcher_type::context_type context_type;

public:

    rule_based_expression_matcher(const value_type& pttrn)
        : pttrn_(pttrn)
    {
        ;
    }

    bool match(const value_type& another)
    {
        return __match(another);
    }

    bool match(
        const value_type& another, const context_type::variable_container_type& globals)
    {
        //XXX: Not implemented yet
        return __match(another);
    }

    bool next()
    {
        return __next();
    }

    //TODO: HERE
    context_type context() const
    {
        context_type ctx;
        if (matchers_.size() == 0)
        {
            return ctx;
        }

        ctx.globals = matchers_.back().context().globals;

        std::vector<unsigned int> strides;
        strides.reserve(another_.size());
        {
            unsigned int stride = 0;
            for (value_type::const_iterator i(another_.begin()); i != another_.end(); ++i)
            {
                strides.push_back(stride);
                stride += (*i).units().size();
            }
        }

        for (std::vector<submatcher_type>::const_iterator
            i(matchers_.begin()); i != matchers_.end(); ++i)
        {
            const unsigned int idx1 = std::distance(matchers_.begin(), i);  // a position in matcher_
            const unsigned int idx2 = permutation_[idx1];  // a position in reactants
            const unsigned int stride = strides[idx2];

            for (context_type::iterator_container_type::const_iterator
                j((*i).context().iterators.begin());
                j != (*i).context().iterators.end(); ++j)
            {
                // // a position in context.iterators
                // const unsigned int idx3 = std::distance((*i).context().iterators.begin(), j);
                const unsigned int idx4 = (*j);  // a position in units of a Species

                ctx.iterators.push_back(idx4 + stride);
            }
        }

        return ctx;
    }
    //TODO: THERE

protected:

    typedef std::vector<value_type::size_type> permutation_type;

    bool __match(const value_type& another)
    {
        permutation_type permutation;
        permutation.reserve(pttrn_.size());
        for (size_t i(0); i != pttrn_.size(); ++i)
        {
            permutation.push_back(i);
        }
        return __match(another, permutation);
    }

    bool __match(const value_type& another, const permutation_type& permutation)
    {
        assert(0 <= pttrn_.size() < 3);
        assert(pttrn_.size() == permutation.size());

        another_ = another; //XXX: copy?
        permutation_ = permutation;

        if (pttrn_.size() != another.size())
        {
            return false;
        }

        matchers_.clear();
        for (value_type::const_iterator i(pttrn_.begin()); i != pttrn_.end(); ++i)
        {
            matchers_.push_back(submatcher_type(*i));
        }

        context_type::variable_container_type globals;
        return this->advance(0, globals);
    }

    bool __next()
    {
        //XXX: Make sure match was already called..
        size_t pos = matchers_.size();

        do
        {
            --pos;

            submatcher_type& matcher_at_pos(matchers_[pos]);
            if (matcher_at_pos.next())
            {
                do
                {
                    if (this->advance(pos + 1, matcher_at_pos.context().globals))
                    {
                        return true;
                    }
                } while (matcher_at_pos.next());
            }
        }
        while (pos > 0);

        return false;
    }

    bool advance(const size_t pos, const context_type::variable_container_type& globals)
    {
        if (pos == matchers_.size())
        {
            return true;
        }

        submatcher_type& matcher_at_pos(matchers_[pos]);
        const value_type::value_type& another_at_pos(another_[pos]);

        if (!matcher_at_pos.match(another_at_pos, globals))
        {
            return false;
        }

        do
        {
            if (this->advance(pos + 1, matcher_at_pos.context().globals))
            {
                return true;
            }
        } while (matcher_at_pos.next());

        return false;
    }

protected:

    value_type pttrn_;
    value_type another_;
    permutation_type permutation_;
    std::vector<submatcher_type> matchers_;
};

} // _context

/**
 * New interfaces for the rule-based modeling
 */

inline Species format_species(const Species& sp)
{
    return context::format_species(sp);
}

struct SpeciesExpressionMatcher
{
    _context::rule_based_expression_matcher<Species> pttrn;

    SpeciesExpressionMatcher(const Species& pttrn)
        : pttrn(pttrn)
    {
        ;
    }

    bool match(const Species& sp)
    {
        return pttrn.match(sp);
    }

    bool next()
    {
        return pttrn.next();
    }

    size_t count(const Species& sp)
    {
        return pttrn.count(sp);
    }
};

std::vector<ReactionRule> generate_reaction_rules(
    const ReactionRule& pttrn,
    const ReactionRule::reactant_container_type& reactants);

} // ecell4

#endif /* ECELL4_CONTEXT_HPP */
