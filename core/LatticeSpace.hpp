#ifndef __ECELL4_LATTICE_SPACE_HPP
#define __ECELL4_LATTICE_SPACE_HPP

#include <vector>
#include <set>
#include <map>
#include <stdexcept>
#include "Space.hpp"
#include "MolecularType.hpp"
#include "VacantType.hpp"
#include "Global.hpp"
#include "LatticeSpaceHDF5Writer.hpp"

namespace ecell4
{

class LatticeSpace
    : public Space
{
protected:

    typedef std::map<Species, MolecularType> spmap;
    typedef std::vector<MolecularTypeBase*> voxel_container;


public:

    LatticeSpace(const Position3& edge_lengths);
    ~LatticeSpace();

    /*
     * Space APIs
     *
     * using ParticleID, Species and Posision3
     */
    const Position3& edge_lengths() const;

    Integer num_species() const;
    Integer num_molecules(const Species& sp)const;
    Integer num_particles() const;
    Integer num_particles(const Species& sp) const;

    bool has_species(const Species& sp) const;
    bool has_particle(const ParticleID& pid) const;

    std::vector<std::pair<ParticleID, Particle> >
        list_particles() const;
    std::vector<std::pair<ParticleID, Particle> >
        list_particles(const Species& sp) const;

    bool update_particle(const ParticleID& pid, const Particle& p);

    const Real t() const
    {
        return this->t_;
    }

    /*
     * for Simulator
     *
     * using Species and Coord
     */
    std::vector<Species> list_species() const;
    std::vector<Coord> list_coords(const Species& sp) const;
    MolecularTypeBase* get_molecular_type(const Species& sp);
    MolecularTypeBase* get_molecular_type(Coord coord) const;
    bool add_species(const Species& sp);
    bool add_molecule(const Species& sp, Coord coord, const ParticleID& pid)
        throw(std::out_of_range);
    bool move(Coord from, Coord to) throw(std::out_of_range);
    bool react(Coord coord, const Species& species) throw(std::out_of_range);

    Real normalized_voxel_radius() const
    {
        return theNormalizedVoxelRadius;
    }

    inline Integer num_col() const
    {
        return col_size_ - 2;
    }

    inline Integer num_row() const
    {
        return row_size_ - 2;
    }

    inline Integer num_layer() const
    {
        return layer_size_ - 2;
    }

    inline Integer num_colrow() const
    {
        return num_col() * num_row();
    }

    inline Integer size() const
    {
        return num_col() * num_row() * num_layer();
    }

    /*
     * Coordinate transformations
     */
    Coord global2coord(const Global& global) const;
    const Global coord2global(Coord coord) const;

    /*
     * HDF5 Save
     */
    void save(H5::H5File* fout, const std::string& hdf5path) const
    {
        LatticeSpaceHDF5Writer<LatticeSpace> writer(*this);
        writer.save(fout, hdf5path);
    }
protected:

    void set_lattice_properties();
    Coord get_coord(const ParticleID& pid) const;
    const Particle particle_at(Coord coord) const;
    bool is_in_range(Coord coord) const;

    /*
     * Coordinate transformations
     */
    Coord global2coord(const Global& global,
            Integer col_size, Integer row_size, Integer layer_size) const;
    const Global coord2global(Coord coord,
            Integer col_size, Integer row_size, Integer layer_size) const;

    const Position3 global2position(const Global& global) const;
    const Global position2global(const Position3& pos) const;

    const Position3 coord2position(Coord coord) const;
    Coord position2coord(const Position3& pos) const;

    Coord inner2general(Coord inner_cood) const;
    Coord general2inner(Coord genetal_coord) const;

protected:

    Real theNormalizedVoxelRadius;
    Real HCP_L, HCP_X, HCP_Y;

    Real t_;
    spmap spmap_;
    voxel_container voxels_;

    MolecularTypeBase* vacant_;
    MolecularTypeBase* border_;

    Position3 edge_lengths_;
    Integer row_size_, layer_size_, col_size_;

};

}

#endif
