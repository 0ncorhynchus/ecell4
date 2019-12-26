#include "SpatiocyteWorld.hpp"
#include "Voxel.hpp"

namespace ecell4
{

namespace spatiocyte
{

Voxel Voxel::get_neighbor_randomly(
    const boost::shared_ptr<RandomNumberGenerator> &rng,
    Shape::dimension_kind dimension) const
{
    std::vector<Voxel> neighbors;
    for (Integer idx = 0; idx < num_neighbors(); ++idx)
    {
        const Voxel neighbor = get_neighbor(idx);
        if (world->get_dimension(neighbor.get_voxel_pool()->species()) >
            dimension)
        {
            continue;
        }
        neighbors.push_back(neighbor);
    }

    const Integer idx(rng->uniform_int(0, neighbors.size() - 1));
    return neighbors.at(idx);
}

} // namespace spatiocyte

} // namespace ecell4
