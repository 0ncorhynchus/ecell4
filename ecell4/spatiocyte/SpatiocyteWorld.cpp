#include <stdexcept>
#include <fstream>

#include "SpatiocyteWorld.hpp"

namespace ecell4
{

namespace spatiocyte
{

void SpatiocyteWorld::add_space(VoxelSpaceBase *space)
{
    for (std::size_t i(0); i < space->size(); ++i)
    {
        const Real3 position(space->coordinate2position(i));
        const coordinate_type nearest(get_root()->position2coordinate(position));

        for (Integer j(0); j < get_root()->num_neighbors(nearest); ++j)
        {
            const coordinate_type neighbor(get_root()->get_neighbor(nearest, j));
            if (length(get_root()->coordinate2position(neighbor) - position) < voxel_radius() * 2)
                interfaces_.add(neighbor, i + size_);
        }
    }

    for (OneToManyMap<coordinate_type>::const_iterator itr(interfaces_.begin());
         itr != interfaces_.end(); ++itr)
    {
        std::vector<coordinate_type> neighbors;
        for (Integer i(0); i < get_root()->num_neighbors((*itr).first); ++i)
        {
            const coordinate_type neighbor(get_root()->get_neighbor((*itr).first, i));
            if (! interfaces_.find(neighbor))
                neighbors.push_back(neighbor);
        }

        for (std::vector<coordinate_type>::const_iterator jtr((*itr).second.begin());
             jtr != (*itr).second.end(); ++jtr)
            neighbors_.extend(*jtr, neighbors);
    }

    spaces_.push_back(space_type(space, size_));

    size_ += space->size();
}

void SpatiocyteWorld::set_value(const Species& sp, const Real value)
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

std::vector<std::pair<ParticleID, Particle> >
SpatiocyteWorld::list_structure_particles() const
{
    const std::vector<Species> structure_species(list_structure_species());

    typedef std::vector<std::vector<std::pair<ParticleID, Particle> > > tmp_type;
    tmp_type tmp_vector(structure_species.size());
    Integer num_elements;

    for (std::vector<Species>::const_iterator itr(structure_species.begin());
            itr != structure_species.end(); ++itr)
    {
        std::vector<std::pair<ParticleID, Particle> > tmp(list_particles(*itr));
        tmp_vector.push_back(tmp);
        num_elements += tmp.size();
    }

    std::vector<std::pair<ParticleID, Particle> > retval;
    retval.reserve(num_elements);
    for (tmp_type::const_iterator itr(tmp_vector.begin());
            itr != tmp_vector.end(); ++itr)
    {
        retval.insert(retval.end(), (*itr).begin(), (*itr).end());
    }

    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
SpatiocyteWorld::list_non_structure_particles() const
{
    const std::vector<Species> non_structure_species(list_non_structure_species());

    typedef std::vector<std::vector<std::pair<ParticleID, Particle> > > tmp_type;
    tmp_type tmp_vector(non_structure_species.size());
    Integer num_elements;

    for (std::vector<Species>::const_iterator itr(non_structure_species.begin());
            itr != non_structure_species.end(); ++itr)
    {
        std::vector<std::pair<ParticleID, Particle> > tmp(list_particles(*itr));
        tmp_vector.push_back(tmp);
        num_elements += tmp.size();
    }

    std::vector<std::pair<ParticleID, Particle> > retval;
    retval.reserve(num_elements);
    for (tmp_type::const_iterator itr(tmp_vector.begin());
            itr != tmp_vector.end(); ++itr)
    {
        retval.insert(retval.end(), (*itr).begin(), (*itr).end());
    }

    return retval;
}

std::vector<Species> SpatiocyteWorld::list_non_structure_species() const
{
    const std::vector<Species> species(list_species());
    std::vector<Species> retval;
    for (std::vector<Species>::const_iterator itr(species.begin());
            itr != species.end(); ++itr)
    {
        if (!find_voxel_pool(*itr)->is_structure())
            retval.push_back(*itr);
    }
    return retval;
}

std::vector<Species> SpatiocyteWorld::list_structure_species() const
{
    const std::vector<Species> species(list_species());
    std::vector<Species> retval;
    for (std::vector<Species>::const_iterator itr(species.begin());
            itr != species.end(); ++itr)
    {
        if (find_voxel_pool(*itr)->is_structure())
            retval.push_back(*itr);
    }
    return retval;
}

bool SpatiocyteWorld::add_molecules(const Species& sp, const Integer& num)
{
    if (num < 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
    }

    const SpatiocyteWorld::molecule_info_type info(get_molecule_info(sp));

    Integer count(0);
    while (count < num)
    {
        const coordinate_type coord(rng()->uniform_int(0, size() - 1));

        if (get_voxel_pool_at(Voxel(coord))->species().serial() != info.loc)
        {
            continue;
        }
        else if (new_voxel(sp, Voxel(coord)))
        {
            ++count;
        }
    }
    return true;
}

bool SpatiocyteWorld::add_molecules(
    const Species& sp, const Integer& num, const boost::shared_ptr<const Shape> shape)
{
    if (num < 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
    }

    const SpatiocyteWorld::molecule_info_type info(get_molecule_info(sp));

    Integer count(0);
    while (count < num)
    {
        const Real3 pos(shape->draw_position(rng_));
        const Voxel voxel(position2voxel(pos));

        if (get_voxel_pool_at(voxel)->species().serial() != info.loc)
        {
            continue;
        }
        else if (new_voxel(sp, voxel))
        {
            ++count;
        }
    }
    return true;
}

Integer SpatiocyteWorld::add_structure(
    const Species& sp, const boost::shared_ptr<const Shape> shape)
{
    const SpatiocyteWorld::molecule_info_type info(get_molecule_info(sp));
    spaces_.at(0).make_structure_type(sp, shape->dimension(), info.loc);

    switch (shape->dimension())
    {
    case Shape::THREE:
        return add_structure3(sp, shape);
    case Shape::TWO:
        return add_structure2(sp, shape);
    case Shape::ONE:
    case Shape::UNDEF:
        break;
    }

    throw NotSupported("The dimension of a shape must be two or three.");
}

Integer SpatiocyteWorld::add_structure3(const Species& sp, const boost::shared_ptr<const Shape> shape)
{
    const SpatiocyteWorld::molecule_info_type info(get_molecule_info(sp));
    Integer count(0);
    for (coordinate_type coord(0); coord < size(); ++coord) {
        const Real L(shape->is_inside(voxel2position(Voxel(coord))));
        if (L > 0)
            continue;

        const ParticleVoxel v(sp, coord, info.radius, info.D, info.loc);
        if (get_voxel_pool_at(Voxel(coord))->species().serial() != info.loc)
        {
            continue;
        }

        if (new_voxel_structure(v).second)
            ++count;
    }
    return count;
}

Integer
SpatiocyteWorld::add_structure2(
        const Species& sp,
        const boost::shared_ptr<const Shape> shape)
{
    const SpatiocyteWorld::molecule_info_type info(get_molecule_info(sp));
    Integer count(0);
    for (coordinate_type coord(0); coord < size(); ++coord) {
        if (!is_surface_voxel(coord, shape))
            continue;

        const ParticleVoxel v(sp, coord, info.radius, info.D, info.loc);
        if (get_voxel_pool_at(Voxel(coord))->species().serial() != info.loc)
        {
            continue;
        }

        if (new_voxel_structure(v).second)
            ++count;
    }
    return count;
}

Integer
SpatiocyteWorld::add_interface(const Species& sp)
{
    const SpatiocyteWorld::molecule_info_type info(get_molecule_info(sp));
    //XXX: set the dimension properly
    spaces_.at(0).make_interface_type(sp, Shape::UNDEF, info.loc);
    return 0;  //XXX: dummpy
}

bool
SpatiocyteWorld::is_surface_voxel(
        const Voxel& voxel,
        const boost::shared_ptr<const Shape> shape) const
{
    const Real L(shape->is_inside(voxel2position(voxel)));
    if (L > 0 || L < -2 * voxel_radius())
        return false;

    for (Integer i(0); i < get_space(voxel)->num_neighbors(voxel.coordinate); ++i)
        if (shape->is_inside(voxel2position(get_neighbor(voxel, i))) > 0)
            return true;

    return false;
}

// TODO
Integer SpatiocyteWorld::add_neighbors(const Species& sp,
    const SpatiocyteWorld::coordinate_type center)
{
    Integer count(0);
    const SpatiocyteWorld::molecule_info_type info(get_molecule_info(sp));
    for (Integer i(0); i < 12; ++i)
    {
        if (new_voxel(sp, get_neighbor(Voxel(center), i)))
        {
            ++count;
        }
        else
        {
            throw "Error in add_neighbors()";
        }
    }
    return count;
}

void SpatiocyteWorld::remove_molecules(const Species& sp, const Integer& num)
{
    if (num < 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
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
        if (remove_voxel(Voxel(mtype->at(idx).coordinate)))
        {
            ++count;
        }
    }
}

boost::optional<Voxel>
SpatiocyteWorld::check_neighbor(const Voxel& voxel, const std::string& loc)
{
    std::vector<coordinate_type> tmp;
    const std::size_t num_neighbors(get_space(voxel)->num_neighbors(voxel.coordinate));
    tmp.reserve(num_neighbors);
    for (unsigned int rnd(0); rnd < num_neighbors; ++rnd)
    {
        const Voxel neighbor(get_neighbor(voxel, rnd));
        boost::shared_ptr<const VoxelPool> mt(get_voxel_pool_at(neighbor));
        const std::string serial(mt->is_vacant() ? "" : mt->species().serial());
        if (serial == loc)
        {
            tmp.push_back(neighbor.coordinate);
        }
    }

    if (tmp.size() == 0)
    {
        return boost::none;
    }

    return Voxel(tmp[rng()->uniform_int(0, tmp.size()-1)]);
}

} // spatiocyte

} // ecell4
