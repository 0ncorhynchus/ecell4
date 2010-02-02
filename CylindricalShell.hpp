#ifndef CYLINDRICALSHELL_HPP
#define CYLINDRICALSHELL_HPP

#include <ostream>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include "Cylinder.hpp"

template<typename T_, typename Tdid_>
struct CylindricalShell
{
    typedef Cylinder<T_> cylinder_type;
    typedef Tdid_ domain_id_type;
    typedef typename cylinder_type::position_type position_type;
    typedef typename cylinder_type::length_type length_type;

    CylindricalShell(): cylinder_(), domain_id_() {}

    CylindricalShell(domain_id_type const& domain_id, cylinder_type const& cylinder)
        : cylinder_(cylinder), domain_id_(domain_id) {}

    length_type calculateDistanceToSelf(position_type pos)
    {
        return cylinder_.calculateDistanceToSelf(pos);
    }

    length_type calculateDistanceToSelfWithOffset(position_type pos, 
                                                  position_type offset)
    {
        return cylinder_.calculateDistanceToSelfWithOffset(pos, offset);
    }

    position_type& position()
    {
        return cylinder_.position();
    }

    position_type const& position() const
    {
        return cylinder_.position();
    }

    length_type& radius()
    {
        return cylinder_.radius();
    }

    length_type const& radius() const
    {
        return cylinder_.radius();
    }

    position_type& orientationZ()
    {
        return cylinder_.orientationZ();
    }

    position_type const& orientationZ() const
    {
        return cylinder_.orientationZ();
    }

    length_type& size()
    {
        return cylinder_.size();
    }

    length_type const& size() const
    {
        return cylinder_.size();
    }


    cylinder_type& as_cylinder()
    {
        return cylinder_;
    }

    cylinder_type const& as_cylinder() const
    {
        return cylinder_;
    }

    domain_id_type const& did() const
    {
        return domain_id_;
    }

    domain_id_type& did()
    {
        return domain_id_;
    }

    bool operator==(CylindricalShell const& rhs) const
    {
        return domain_id_ == rhs.did() && cylinder_ == rhs.as_cylinder();
    }

    bool operator!=(CylindricalShell const& rhs) const
    {
        return !operator==(rhs);
    }

private:
    cylinder_type cylinder_;
    domain_id_type domain_id_;
};

template<typename Tstrm_, typename Ttraits_, typename T_, typename Tdid_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm, const CylindricalShell<T_, Tdid_>& v)
{
    strm << "CylindricalShell(" << v.as_cylinder() << ", " << v.did() << ")";
    return strm;
}

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std { namespace tr1 {
#elif defined(HAVE_STD_HASH)
namespace std {
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost {
#endif

template<typename T_, typename Tdid_>
struct hash<CylindricalShell<T_, Tdid_> >
{
    typedef CylindricalShell<T_, Tdid_> argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<typename argument_type::position_type>()(val.position()) ^
            hash<typename argument_type::length_type>()(val.radius()) ^
            hash<typename argument_type::position_type>()(val.orientationZ()) ^
            hash<typename argument_type::length_type>()(val.size()) ^
            hash<typename argument_type::domain_id_type>()(val.did());
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} } // namespace std::tr1
#elif defined(HAVE_STD_HASH)
} // namespace std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // namespace boost
#endif

#endif /* SHELL_HPP */
