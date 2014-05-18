#ifndef POSITION3_TRAITS_HPP
#define POSITION3_TRAITS_HPP

// This file containing the reference to the ecell4::Position3 and 
//   some template traits of ecell4::Position3.
//

#include <ostream>
#include <iomanip>
#include <functional>
#include <algorithm>

#include <ecell4/core/Position3.hpp>

#include <boost/array.hpp>
#include "utils/array_traits.hpp"
#include "linear_algebra.hpp"

template <std::size_t N_>
struct is_vector<ecell4::Position3, N_>: public boost::mpl::true_ {};

template <>
struct element_type_of<ecell4::Position3>
{
    typedef ecell4::Position3::value_type type;
};

#endif /* POSITION3_TRAITS_HPP */
