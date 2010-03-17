#ifndef SHELL_HPP
#define SHELL_HPP

#include <ostream>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

template<typename Tshape_, typename Tdid_>
struct Shell
{
    typedef Tshape_ shape_type;
    typedef Tdid_ domain_id_type;
    typedef typename shape_type::position_type position_type;
    typedef typename shape_type::length_type length_type;

    Shell(): sphere_(), domain_id_() {}

    Shell(domain_id_type const& domain_id, shape_type const& sphere)
        : sphere_(sphere), domain_id_(domain_id) {}

    position_type& position()
    {
        return sphere_.position();
    }

    position_type const& position() const
    {
        return sphere_.position();
    }

    length_type& radius()
    {
        return sphere_.radius();
    }

    length_type const& radius() const
    {
        return sphere_.radius();
    }

    shape_type& shape()
    {
        return sphere_;
    }

    shape_type const& shape() const
    {
        return sphere_;
    }

    domain_id_type const& did() const
    {
        return domain_id_;
    }

    domain_id_type& did()
    {
        return domain_id_;
    }

    bool operator==(Shell const& rhs) const
    {
        return domain_id_ == rhs.did() && sphere_ == rhs.shape();
    }

    bool operator!=(Shell const& rhs) const
    {
        return !operator==(rhs);
    }

private:
    shape_type sphere_;
    domain_id_type domain_id_;
};

template<typename Tstrm_, typename Ttraits_, typename Tshape_, typename Tdid_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm, const Shell<Tshape_, Tdid_>& v)
{
    strm << "Shell(" << v.shape() << ", " << v.did() << ")";
    return strm;
}

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std { namespace tr1 {
#elif defined(HAVE_STD_HASH)
namespace std {
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost {
#endif

template<typename Tshape_, typename Tdid_>
struct hash<Shell<Tshape_, Tdid_> >
{
    typedef Shell<Tshape_, Tdid_> argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<typename argument_type::position_type>()(val.position()) ^
            hash<typename argument_type::length_type>()(val.radius()) ^
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
