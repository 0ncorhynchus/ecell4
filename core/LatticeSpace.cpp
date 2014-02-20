#include "LatticeSpace.hpp"

namespace ecell4
{

LatticeSpace::LatticeSpace(const Position3& edge_lengths) : t_(0)
{
    this->edge_lengths_ = edge_lengths;
    vacant_ = new VacantType();
    border_ = new MolecularType(Species("Border", "0"));
    set_lattice_properties();
}

LatticeSpace::~LatticeSpace()
{
    delete vacant_;
    delete border_;
}

/*
 * derived from SpatiocyteStepper::setLatticeProperties()
 */
void LatticeSpace::set_lattice_properties()
{
    theNormalizedVoxelRadius = 2.5e-9;

    HCP_L = theNormalizedVoxelRadius/sqrt(3);
    HCP_X = theNormalizedVoxelRadius*sqrt(8.0/3); //Lx
    HCP_Y = theNormalizedVoxelRadius*sqrt(3); //Ly

    const Real lengthX = edge_lengths_[0];
    const Real lengthY = edge_lengths_[1];
    const Real lengthZ = edge_lengths_[2];

    row_size_ = (Integer)rint((lengthZ/2)/theNormalizedVoxelRadius) + 2;
    layer_size_ = (Integer)rint(lengthY/HCP_Y) + 2;
    col_size_ = (Integer)rint(lengthX/HCP_X) + 2;

    for (Coord coord(0); coord < row_size_ * layer_size_ * col_size_; ++coord)
    {
        Global global(coord2global(coord, col_size_, row_size_, layer_size_));
        if (global.col == 0 || global.col == col_size_-1 ||
                global.row == 0 || global.row == row_size_-1 ||
                global.layer == 0 || global.layer == layer_size_-1)
        {
            voxels_.push_back(border_);
        }
        else
        {
            voxels_.push_back(vacant_);
        }
    }
}

Integer LatticeSpace::num_species() const
{
    return spmap_.size();
}

bool LatticeSpace::has_species(const Species& sp) const
{
    return spmap_.find(sp) != spmap_.end();
}

Integer LatticeSpace::num_molecules(const Species& sp) const
{
    return num_particles(sp) * 1;
}

const Position3& LatticeSpace::edge_lengths() const
{
    return edge_lengths_;
}

Integer LatticeSpace::num_particles() const
{
    Integer count(0);
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        count += (*itr).second.voxels().size();
    }

    return count;
}

Integer LatticeSpace::num_particles(const Species& sp) const
{
    spmap::const_iterator itr(spmap_.find(sp));
    Integer count(0);
    if (itr != spmap_.end())
    {
        count = (*itr).second.voxels().size();
    }
    return count;
}

bool LatticeSpace::has_particle(const ParticleID& pid) const
{
    bool flg(false);
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const MolecularTypeBase& mt((*itr).second);
        if (mt.is_vacant())
        {
            return false;
        }
        for (MolecularTypeBase::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
        {
            if ((*vitr).second == pid)
            {
                flg = true;
                break;
            }
        }
    }

    return flg;
}

std::vector<std::pair<ParticleID, Particle> >
    LatticeSpace::list_particles() const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const MolecularTypeBase& mt((*itr).second);
        if (mt.is_vacant())
        {
            continue;
        }
        for (MolecularTypeBase::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
        {
            retval.push_back(std::pair<ParticleID, Particle>(
                        (*vitr).second, particle_at((*vitr).first)));
        }
    }

    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
    LatticeSpace::list_particles(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    spmap::const_iterator itr(spmap_.find(sp));
    if (itr != spmap_.end())
    {
        const MolecularTypeBase& mt((*itr).second);
        if (mt.is_vacant())
        {
            return retval;
        }
        for (MolecularTypeBase::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
        {
            retval.push_back(std::pair<ParticleID, Particle>(
                        (*vitr).second, particle_at((*vitr).first)));
        }
    }

    return retval;
}

bool LatticeSpace::update_particle(const ParticleID& pid, const Particle& p)
{
    Coord coord(position2coord(p.position()));
    MolecularTypeBase* dest_mt = get_molecular_type(p.species());

    if (!is_in_range(coord))
    {
        return false;
    }
    if (has_particle(pid))
    {
        Coord coord2(get_coord(pid));
        MolecularTypeBase* src_ptr_mt(voxels_.at(coord));
        src_ptr_mt->removeVoxel(coord2);
        voxel_container::iterator itr(voxels_.begin() + coord2);
        voxels_.erase(itr);
        voxels_.insert(itr, NULL);
    }
    dest_mt->addVoxel(MolecularTypeBase::particle_info(coord, pid));

    voxel_container::iterator itr(voxels_.begin() + coord);
    (*itr) = dest_mt;

    return true;
}

/*
 * original methods
 */


std::vector<Species> LatticeSpace::list_species() const
{
    std::vector<Species> keys;
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        keys.push_back((*itr).first);
    }
    return keys;
}

std::vector<Coord> LatticeSpace::list_coords(const Species& sp) const
{
    std::vector<Coord> retval;
    spmap::const_iterator itr(spmap_.find(sp));
    if (itr == spmap_.end())
    {
        return retval;
    }

    const MolecularTypeBase* mt(&((*itr).second));

    for (MolecularTypeBase::container_type::const_iterator itr(mt->begin());
            itr != mt->end(); ++itr)
    {
        retval.push_back(general2inner((*itr).first));
    }

    return retval;
}

MolecularTypeBase* LatticeSpace::get_molecular_type(const Species& sp)
{
    spmap::iterator itr(spmap_.find(sp));
    if (itr == spmap_.end())
    {
        MolecularType mt(sp);
        std::pair<spmap::iterator, bool> result = spmap_.insert(spmap::value_type(sp, mt));
        if (result.second)
            itr = result.first;
        else
            throw "insert error";
    }
    return &((*itr).second);
}


/*
 * Coordinate transformations
 */

Coord LatticeSpace::global2coord(const Global& global) const
{
    return global2coord(global, num_col(), num_row(), num_layer());
}

const Global LatticeSpace::coord2global(Coord coord) const
{
    return coord2global(coord, num_col(), num_row(), num_layer());
}

/*
 * Protected functions
 */

Coord LatticeSpace::get_coord(const ParticleID& pid) const
{
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const MolecularTypeBase& mt((*itr).second);
        if (mt.is_vacant())
        {
            return false;
        }
        for (MolecularTypeBase::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
        {
            if ((*vitr).second == pid)
            {
                return (*vitr).first;
            }
        }
    }
    throw "Exception: Not in lattice";
}

MolecularTypeBase* LatticeSpace::get_molecular_type(Coord coord) const
{
    return voxels_.at(inner2general(coord));
}


bool LatticeSpace::add_species(const Species& sp)
{
    if (has_species(sp))
    {
        return false;
    }
    MolecularType mt(sp);
    std::pair<spmap::iterator, bool> retval = spmap_.insert(spmap::value_type(sp, mt));
    return retval.second;
}

bool LatticeSpace::add_molecule(const Species& sp, Coord coord, const ParticleID& pid)
        throw(std::out_of_range)
{
    if (!is_in_range(coord))
    {
        throw std::out_of_range("");
        return false;
    }

    MolecularTypeBase* mt(get_molecular_type(sp));
    if (mt->is_vacant())
    {
        std::cerr << "[" << sp.name() << " is vacant]";
        return false;
    }

    const MolecularTypeBase* mt_at(get_molecular_type(coord));
    if (!mt_at->is_vacant())
    {
        std::cerr << "[" << mt_at->species().name()  << " at " << coord << " is not vacant]"; //DEBUG
        return false;
    }
    const Coord general_coord(inner2general(coord));
    mt->addVoxel(std::pair<Coord, ParticleID>(general_coord, pid));

    voxel_container::iterator itr(voxels_.begin() + general_coord);
    (*itr) = mt;

    return true;
}

bool LatticeSpace::move(Coord from, Coord to) throw(std::out_of_range)
{
    if (!is_in_range(from) || !is_in_range(to))
    {
        throw std::out_of_range("");
        return false;
    }

    const MolecularTypeBase* to_mt(get_molecular_type(to));

    if (!to_mt->is_vacant()) {
        return false;
    }

    MolecularTypeBase* from_mt(get_molecular_type(from));

    const Coord general_from(inner2general(from));
    const Coord general_to(inner2general(to));

    if (!from_mt->is_vacant()) {
        MolecularTypeBase::container_type::iterator itr(from_mt->find(general_from));
        (*itr).first = general_to;
        voxel_container::iterator fitr(voxels_.begin() + general_from);
        (*fitr) = vacant_;
        voxel_container::iterator titr(voxels_.begin() + general_to);
        (*titr) = from_mt;
        std::cout << "[from: " << from << ", to: " << to << "] "; // DEBUG
    }

    return true;
}

bool LatticeSpace::react(Coord coord, const Species& species) throw(std::out_of_range)
{
    if (!is_in_range(coord))
    {
        throw std::out_of_range("");
        return false;
    }

    MolecularTypeBase* old_mt(get_molecular_type(coord));
    if (old_mt->is_vacant())
    {
        return false;
    }

    MolecularTypeBase* new_mt(get_molecular_type(species));

    const Coord general_coord(inner2general(coord));
    MolecularTypeBase::particle_info info(*old_mt->find(general_coord));
    old_mt->removeVoxel(general_coord);
    new_mt->addVoxel(info);
    voxel_container::iterator itr(voxels_.begin() + general_coord);
    (*itr) = new_mt;
    return true;
}

const Particle LatticeSpace::particle_at(Coord coord) const
{
    const MolecularTypeBase* ptr_mt(voxels_.at(coord));
    const Species& sp = ptr_mt->species();
    const Position3& pos = coord2position(coord);
    const Real& radius = 0;
    const Real& D = 0;
    Particle particle(sp, pos, radius, D);
    return particle;
}

/*
 * Coordinate transformations
 */

bool LatticeSpace::is_in_range(Coord coord) const
{
    coord = inner2general(coord);
    return coord >= 0 && coord < row_size_ * layer_size_ * col_size_;
}

Coord LatticeSpace::global2coord(const Global& global,
        Integer col_size, Integer row_size, Integer layer_size) const
{
    return (global.row) +
        row_size * (global.layer) +
        row_size * layer_size * (global.col);
}

const Global LatticeSpace::coord2global(Coord coord,
        Integer col_size, Integer row_size, Integer layer_size) const
{
    Global retval;
    retval.col = coord / (row_size * layer_size);
    retval.layer = (coord % (row_size * layer_size)) / row_size;
    retval.row = (coord % (row_size * layer_size)) % row_size;
    return retval;
}

const Position3 LatticeSpace::global2position(const Global& global) const
{
    //the center point of a voxel
    Position3 position;
    position[0] = global.col * HCP_X;
    position[1] = (global.col % 2) * HCP_L + HCP_Y * global.layer;
    position[2] = (global.row * 2 + (global.layer + global.col) % 2)
        * theNormalizedVoxelRadius;
    return position;
}

const Global LatticeSpace::position2global(const Position3& pos) const
{
    Global global;
    global.col = (Integer)(pos[0] / HCP_X);
    global.layer = (Integer)((pos[1] - (global.col % 2) * HCP_L) / HCP_Y);
    global.row = (Integer)((pos[2] - (((global.layer + global.col) % 2) *
                    theNormalizedVoxelRadius)) /
            theNormalizedVoxelRadius / 2);
    return global;
}

const Position3 LatticeSpace::coord2position(Coord coord) const
{
    return global2position(coord2global(coord, col_size_, row_size_, layer_size_));
}

Coord LatticeSpace::position2coord(const Position3& pos) const
{
    return global2coord(position2global(pos), col_size_, row_size_, layer_size_);
}


Coord LatticeSpace::inner2general(Coord inner_coord) const
{
    Global global(coord2global(inner_coord));
    global.col += 1;
    global.row += 1;
    global.layer += 1;
    return global2coord(global, col_size_, row_size_, layer_size_);
}

Coord LatticeSpace::general2inner(Coord general_coord) const
{
    Global global(coord2global(general_coord,
                col_size_, row_size_, layer_size_));
    global.col -= 1;
    global.row -= 1;
    global.layer -= 1;
    return global2coord(global);
}

} // ecell4
