#include "utils.hpp"


namespace ecell4
{

namespace spatiocyte
{

const Real calculate_dimensional_factor(
    const Species& speciesA, const Species& speciesB,
    boost::shared_ptr<SpatiocyteWorld> world)
{
    const MoleculeInfo minfoA(world->get_molecule_info(speciesA));
    const MoleculeInfo minfoB(world->get_molecule_info(speciesB));

    const Real sqrt2(sqrt(2.0));

    if (minfoA.dimension == Shape::THREE && minfoB.dimension == Shape::THREE)
    {
        return 1. / (6 * sqrt2 * (minfoA.D + minfoB.D) * world->voxel_radius());
    }

    if (minfoA.dimension == Shape::TWO && minfoB.dimension == Shape::TWO)
    {
        const Real sqrt3(sqrt(3.0));
        const Real sqrt6(sqrt(6.0));
        const Real gamma(pow(2 * sqrt2 + 4 * sqrt3 + 3 * sqrt6 + sqrt(22.0), 2) /
            (72 * (6 * sqrt2 + 4 * sqrt3 + 3 * sqrt6)));
        return gamma / (minfoA.D + minfoB.D);
    }

    if (minfoA.dimension == Shape::THREE && minfoB.dimension == Shape::TWO)
    {
        const Real factor = sqrt2 / (3 * minfoA.D * world->voxel_radius());
        if (minfoB.is_structure) // B is Surface
        {
            return factor * world->unit_area();
        }
        return factor;
    }

    if (minfoA.dimension == Shape::TWO && minfoB.dimension == Shape::THREE)
    {
        const Real factor = sqrt2 / (3 * minfoB.D * world->voxel_radius());
        if (minfoA.is_structure) // A is Surface
        {
            return factor * world->unit_area();
        }
        return factor;
    }

    throw NotSupported("The dimension of a structure must be two or three.");
}

} // spatiocyte

} // ecell4
