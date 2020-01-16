#include <fstream>
#include <stdexcept>

#include "SpatiocyteWorld.hpp"

namespace ecell4
{

namespace spatiocyte
{

void SpatiocyteWorld::add_space(std::unique_ptr<VoxelSpaceBase> space)
{
    for (coordinate_type i(0); i < space->size(); ++i)
    {
        const Real3 position(space->coordinate2position(i));
        const coordinate_type nearest(
            get_root()->position2coordinate(position));

        for (Integer j(0); j < get_root()->num_neighbors(nearest); ++j)
        {
            const coordinate_type neighbor(
                get_root()->get_neighbor(nearest, j));
            if (length(get_root()->coordinate2position(neighbor) - position) <
                voxel_radius() * 2)
                interfaces_.add(neighbor, i + size_);
        }
    }

    for (const auto &interface : interfaces_)
    {
        std::vector<coordinate_type> neighbors;
        for (Integer i(0); i < get_root()->num_neighbors(interface.first); ++i)
        {
            const coordinate_type neighbor(
                get_root()->get_neighbor(interface.first, i));
            if (!interfaces_.find(neighbor))
                neighbors.push_back(neighbor);
        }

        for (const auto &adjoining : interface.second)
        {
            neighbors_.extend(adjoining, neighbors);
        }
    }

    size_ += space->size();
    spaces_.push_back(space_type(space.release()));
}

void SpatiocyteWorld::set_value(const Species &sp, const Real value)
{
    const Integer num1 = static_cast<Integer>(value);
    const Integer num2 = num_molecules_exact(sp);
    if (num1 > num2)
    {
        add_molecules(sp, num1 - num2);
    }
    else if (num1 < num2)
    {
        remove_molecules(sp, num2 - num1);
    }
}

std::vector<std::pair<ParticleID, Particle>>
SpatiocyteWorld::list_structure_particles() const
{
    const std::vector<Species> structure_species(list_structure_species());

    typedef std::vector<std::vector<std::pair<ParticleID, Particle>>> tmp_type;
    tmp_type tmp_vector(structure_species.size());
    Integer num_elements(0);

    for (const auto &species : structure_species)
    {
        std::vector<std::pair<ParticleID, Particle>> tmp(
            list_particles(species));
        tmp_vector.push_back(tmp);
        num_elements += tmp.size();
    }

    std::vector<std::pair<ParticleID, Particle>> retval;
    retval.reserve(num_elements);
    for (const auto &tmp : tmp_vector)
    {
        retval.insert(retval.end(), tmp.begin(), tmp.end());
    }

    return retval;
}

std::vector<std::pair<ParticleID, Particle>>
SpatiocyteWorld::list_non_structure_particles() const
{
    const std::vector<Species> non_structure_species(
        list_non_structure_species());

    typedef std::vector<std::vector<std::pair<ParticleID, Particle>>> tmp_type;
    tmp_type tmp_vector(non_structure_species.size());
    Integer num_elements(0);

    for (const auto &species : non_structure_species)
    {
        std::vector<std::pair<ParticleID, Particle>> tmp(
            list_particles(species));
        tmp_vector.push_back(tmp);
        num_elements += tmp.size();
    }

    std::vector<std::pair<ParticleID, Particle>> retval;
    retval.reserve(num_elements);
    for (const auto &tmp : tmp_vector)
    {
        retval.insert(retval.end(), tmp.begin(), tmp.end());
    }

    return retval;
}

std::vector<Species> SpatiocyteWorld::list_non_structure_species() const
{
    std::vector<Species> retval;
    for (const auto &species : list_species())
    {
        if (!find_voxel_pool(species)->is_structure())
            retval.push_back(species);
    }
    return retval;
}

std::vector<Species> SpatiocyteWorld::list_structure_species() const
{
    std::vector<Species> retval;
    for (const auto &species : list_species())
    {
        if (find_voxel_pool(species)->is_structure())
            retval.push_back(species);
    }
    return retval;
}

bool SpatiocyteWorld::add_molecules(const Species &sp, const Integer &num)
{
    if (num < 0)
    {
        throw std::invalid_argument(
            "The number of molecules must be positive.");
    }

    const MoleculeInfo info(get_molecule_info(sp));

    Integer count(0);
    while (count < num)
    {
        const Voxel voxel(coordinate2voxel(rng()->uniform_int(0, size() - 1)));

        if (voxel.get_voxel_pool()->species().serial() != info.loc)
        {
            continue;
        }
        else if (new_particle(sp, voxel))
        {
            ++count;
        }
    }
    return true;
}

bool SpatiocyteWorld::add_molecules(const Species &sp, const Integer &num,
                                    const boost::shared_ptr<const Shape> shape)
{
    if (num < 0)
    {
        throw std::invalid_argument(
            "The number of molecules must be positive.");
    }

    const MoleculeInfo info(get_molecule_info(sp));

    Integer count(0);
    while (count < num)
    {
        const Real3 pos(shape->draw_position(rng_));
        const Voxel voxel(get_voxel_nearby(pos));

        if (voxel.get_voxel_pool()->species().serial() != info.loc)
        {
            continue;
        }
        else if (new_particle(sp, voxel))
        {
            ++count;
        }
    }
    return true;
}

Integer
SpatiocyteWorld::add_structure(const Species &sp,
                               const boost::shared_ptr<const Shape> shape)
{
    const MoleculeInfo info(get_molecule_info(sp));
    get_root()->make_structure_type(sp, info.loc);

    if (shape->dimension() != info.dimension)
    {
        throw IllegalArgument("The dimension mismatch occurred between a given "
                              "species and shape");
    }

    switch (shape->dimension())
    {
    case Shape::THREE:
        return add_structure3(sp, info.loc, shape);
    case Shape::TWO:
        return add_structure2(sp, info.loc, shape);
    case Shape::ONE:
    case Shape::UNDEF:
        break;
    }

    throw NotSupported("The dimension of a shape must be two or three.");
}

Integer
SpatiocyteWorld::add_structure3(const Species &sp, const std::string &location,
                                const boost::shared_ptr<const Shape> shape)
{
    Integer count(0);
    for (coordinate_type coord(0); coord < size(); ++coord)
    {
        const Voxel voxel(coordinate2voxel(coord));

        if (!this->is_inside(coord) || shape->is_inside(voxel.position()) > 0)
        {
            continue;
        }

        if (voxel.get_voxel_pool()->species().serial() != location)
        {
            throw NotSupported("Mismatch in the location. Failed to place '" +
                               sp.serial() + "' to '" +
                               voxel.get_voxel_pool()->species().serial() +
                               "'. " + "'" + location + "' is expected.");
            continue;
        }

        if (new_voxel_structure(sp, voxel))
            ++count;
    }
    return count;
}

Integer
SpatiocyteWorld::add_structure2(const Species &sp, const std::string &location,
                                const boost::shared_ptr<const Shape> shape)
{
    Integer count(0);
    for (coordinate_type coord(0); coord < size(); ++coord)
    {
        const Voxel voxel(coordinate2voxel(coord));
        if (!this->is_inside(coord) || !is_surface_voxel(voxel, shape))
        {
            continue;
        }

        if (voxel.get_voxel_pool()->species().serial() != location)
        {
            throw NotSupported("Mismatch in the location. Failed to place '" +
                               sp.serial() + "' to '" +
                               voxel.get_voxel_pool()->species().serial() +
                               "'. " + "'" + location + "' is expected.");
            continue;
        }

        if (new_voxel_structure(sp, voxel))
            ++count;
    }
    return count;
}

bool SpatiocyteWorld::is_surface_voxel(
    const Voxel &voxel, const boost::shared_ptr<const Shape> shape) const
{
    const Real L(shape->is_inside(voxel.position()));
    if (L > 0 || L < -2 * voxel_radius())
        return false;

    for (Integer i(0); i < num_neighbors(voxel); ++i)
        if (shape->is_inside(get_neighbor(voxel, i).position()) > 0)
            return true;

    return false;
}

void SpatiocyteWorld::remove_molecules(const Species &sp, const Integer &num)
{
    if (num < 0)
    {
        throw std::invalid_argument(
            "The number of molecules must be positive.");
    }

    boost::shared_ptr<const MoleculePool> mtype(find_molecule_pool(sp));
    if (mtype->size() < num)
    {
        throw std::invalid_argument(
            "The number of molecules cannot be negative.");
    }

    Integer count(0);
    while (count < num)
    {
        const Integer idx(rng_->uniform_int(0, mtype->size() - 1));
        if (coordinate2voxel(mtype->at(idx).coordinate).clear())
        {
            ++count;
        }
    }
}

boost::optional<Voxel> SpatiocyteWorld::check_neighbor(const Voxel &voxel,
                                                       const std::string &loc)
{
    const std::size_t num(num_neighbors(voxel));

    std::vector<Voxel> tmp;
    tmp.reserve(num);

    for (unsigned int rnd(0); rnd < num; ++rnd)
    {
        const Voxel neighbor(get_neighbor(voxel, rnd));
        boost::shared_ptr<const VoxelPool> mt(neighbor.get_voxel_pool());
        const std::string serial(mt->is_vacant() ? "" : mt->species().serial());
        if (serial == loc)
        {
            tmp.push_back(neighbor);
        }
    }

    if (tmp.size() == 0)
    {
        return boost::none;
    }

    return tmp[rng()->uniform_int(0, tmp.size() - 1)];
}

const Voxel SpatiocyteWorld::get_neighbor_randomly(const Voxel &voxel) const
{
    const Integer idx(rng()->uniform_int(0, num_neighbors(voxel) - 1));
    return get_neighbor(voxel, idx);
}

const Voxel
SpatiocyteWorld::get_neighbor_randomly(const Voxel &voxel,
                                       Shape::dimension_kind dimension) const
{
    std::vector<Voxel> neighbors;
    for (Integer idx = 0; idx < num_neighbors(voxel); ++idx)
    {
        const Voxel neighbor = get_neighbor(voxel, idx);
        if (get_dimension(neighbor.get_voxel_pool()->species()) > dimension)
        {
            continue;
        }
        neighbors.push_back(neighbor);
    }

    const Integer idx(rng()->uniform_int(0, neighbors.size() - 1));
    return neighbors.at(idx);
}

} // namespace spatiocyte

} // namespace ecell4
