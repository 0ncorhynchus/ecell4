#ifndef ECELL4_SPATIOCYTE_UTILS_HPP
#define ECELL4_SPATIOCYTE_UTILS_HPP

#include "SpatiocyteWorld.hpp"

namespace ecell4
{

namespace spatiocyte
{

const Real calculate_dimensional_factor(
    const Species& speciesA,
    const Species& speciesB,
    boost::shared_ptr<SpatiocyteWorld> world);

} // spatiocyte

} // ecell4

#endif /* ECELL4_SPATIOCYTE_UTILS_HPP */
