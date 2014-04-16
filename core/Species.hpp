#ifndef __ECELL4_SPECIES_HPP
#define __ECELL4_SPECIES_HPP

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>

#include "config.h"

#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include "get_mapper_mf.hpp"
#include "types.hpp"
#include "exceptions.hpp"
#include "UnitSpecies.hpp"


namespace ecell4
{

class Species
{
public:

    typedef std::string serial_type;
    typedef std::vector<UnitSpecies> container_type;

protected:

    typedef utils::get_mapper_mf<std::string, std::string>::type
    attributes_container_type;

public:

    Species()
        : units_()
    {
        ; // do nothing
    }

    Species(const std::string& name)
        : units_()
    {
        units_.push_back(UnitSpecies(name));
    }

    Species(
        const std::string& name, const std::string& D)
        : units_()
    {
        units_.push_back(UnitSpecies(name));
        set_attribute("D", D);
    }

    Species(
        const std::string& name, const std::string& radius, const std::string& D)
        : units_()
    {
        units_.push_back(UnitSpecies(name));
        set_attribute("radius", radius);
        set_attribute("D", D);
    }

    serial_type serial() const
    {
        if (units_.size() == 0)
        {
            return "";
        }

        container_type copied(units_.size());
        std::copy(units_.begin(), units_.end(), copied.begin());
        std::sort(copied.begin(), copied.end());

        container_type::const_iterator it(copied.begin());
        std::string retval((*it).name());
        ++it;
        for (; it != copied.end(); ++it)
        {
            retval += ".";
            retval += (*it).name();
        }
        return retval;
    }

    std::string name() const
    {
        if (units_.size() == 0)
        {
            return "";
        }

        container_type::const_iterator it(units_.begin());
        std::string retval((*it).name());
        ++it;
        for (; it != units_.end(); ++it)
        {
            retval += ".";
            retval += (*it).name();
        }
        return retval;
    }

    void add_unit(const UnitSpecies& usp)
    {
        units_.push_back(usp);
    }

    inline container_type::const_iterator begin() const
    {
        return units_.begin();
    }

    inline container_type::const_iterator end() const
    {
        return units_.end();
    }

    bool match(const Species& target) const
    {
        for (container_type::const_iterator i(units_.begin()); i != units_.end(); ++i)
        {
            if (std::count(target.begin(), target.end(), (*i))
                < std::count(units_.begin(), units_.end(), (*i))) //XXX:
            {
                return false;
            }
        }
        return true;
    }

    const attributes_container_type& attributes() const
    {
        return attributes_;
    }

    std::string get_attribute(const std::string& name_attr) const
    {
        attributes_container_type::const_iterator
            i(attributes_.find(name_attr));
        if (i == attributes_.end())
        {
            std::ostringstream message;
            message << "attribute [" << name_attr << "] not found";
            throw NotFound(message.str()); // use boost::format if it's allowed
        }

        return (*i).second;
    }

    void set_attribute(const std::string& name_attr, const std::string& value)
    {
        attributes_[name_attr] = value;
    }

    void set_attributes(const Species& sp)
    {
        attributes_ = sp.attributes();
    }

    void remove_attribute(const std::string& name_attr)
    {
        attributes_container_type::iterator
            i(attributes_.find(name_attr));
        if (i == attributes_.end())
        {
            std::ostringstream message;
            message << "attribute [" << name_attr << "] not found";
            throw NotFound(message.str()); // use boost::format if it's allowed
        }

        attributes_.erase(i);
    }

    bool has_attribute(const std::string& name_attr) const
    {
        return (attributes_.find(name_attr) != attributes_.end());
    }

    bool operator==(const Species& rhs) const;
    bool operator<(const Species& rhs) const;
    bool operator>(const Species& rhs) const;

protected:

    std::vector<UnitSpecies> units_;
    attributes_container_type attributes_;
};

} // ecell4

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std
{

namespace tr1
{
#elif defined(HAVE_STD_HASH)
namespace std
{
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost
{
#endif

template<>
struct hash<ecell4::Species>
{
    std::size_t operator()(const ecell4::Species& val) const
    {
        return hash<ecell4::Species::serial_type>()(val.serial());
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} // tr1

} // std
#elif defined(HAVE_STD_HASH)
} // std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // boost
#endif

#endif /* __ECELL4_SPECIES_HPP */
