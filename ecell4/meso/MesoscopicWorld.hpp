#ifndef __ECELL4_MESO_MESOSCOPIC_WORLD_HPP
#define __ECELL4_MESO_MESOSCOPIC_WORLD_HPP

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/SubvolumeSpace.hpp>
#include <ecell4/core/Model.hpp>

namespace ecell4
{

namespace meso
{

struct MoleculeInfo
{
    const Real D;
};

class MesoscopicWorld
    : public Space
{
public:

    typedef SubvolumeSpace::coordinate_type coordinate_type;

public:

    MesoscopicWorld(const Position3& edge_lengths,
        const Integer& cx, const Integer& cy, const Integer& cz,
        boost::shared_ptr<RandomNumberGenerator> rng)
        : cs_(new SubvolumeSpaceVectorImpl(edge_lengths, cx, cy, cz)), rng_(rng)
    {
        ;
    }

    MesoscopicWorld(const Position3& edge_lengths,
        const Integer& cx, const Integer& cy, const Integer& cz)
        : cs_(new SubvolumeSpaceVectorImpl(edge_lengths, cx, cy, cz))
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    virtual ~MesoscopicWorld()
    {
        ;
    }

    void bind_to(boost::shared_ptr<Model> model)
    {
        if (boost::shared_ptr<Model> bound_model = lock_model())
        {
            if (bound_model.get() != model.get())
            {
                std::cerr << "Warning: Model already bound to GillespieWorld."
                    << std::endl;
            }
        }
        this->model_ = model;
    }

    boost::shared_ptr<Model> lock_model() const
    {
        return model_.lock();
    }

    inline const boost::shared_ptr<RandomNumberGenerator>& rng()
    {
        return rng_;
    }

    MoleculeInfo get_molecule_info(const Species& sp) const;

    const Real& t() const;
    void set_t(const Real& t);
    const Integer num_subvolumes() const;
    const Real subvolume() const;
    const Real volume() const;
    const Position3 subvolume_edge_lengths() const;
    const Position3& edge_lengths() const;

    coordinate_type global2coord(const Global& g) const;
    Global coord2global(const coordinate_type& c) const;

    Real get_value(const Species& sp) const;
    Real get_value_exact(const Species& sp) const;
    Integer num_molecules(const Species& sp) const;
    Integer num_molecules_exact(const Species& sp) const;
    Integer num_molecules(const Species& sp, const coordinate_type& c) const;
    Integer num_molecules_exact(const Species& sp, const coordinate_type& c) const;
    void add_molecules(const Species& sp, const Integer& num, const coordinate_type& c);
    void remove_molecules(const Species& sp, const Integer& num, const coordinate_type& c);

    Integer num_molecules(const Species& sp, const Global& g) const
    {
        return cs_->num_molecules(sp, g);
    }

    Integer num_molecules_exact(const Species& sp, const Global& g) const
    {
        return cs_->num_molecules_exact(sp, g);
    }

    void add_molecules(const Species& sp, const Integer& num, const Global& g)
    {
        cs_->add_molecules(sp, num, g);
    }

    void remove_molecules(const Species& sp, const Integer& num, const Global& g)
    {
        cs_->remove_molecules(sp, num, g);
    }

    const std::vector<Species>& species() const;
    std::vector<Species> list_species() const;

    std::vector<std::pair<ParticleID, Particle> > list_particles() const;
    std::vector<std::pair<ParticleID, Particle> > list_particles_exact(const Species& sp) const;
    std::vector<std::pair<ParticleID, Particle> > list_particles(const Species& sp) const;

private:

    boost::scoped_ptr<SubvolumeSpace> cs_;
    boost::shared_ptr<RandomNumberGenerator> rng_;

    boost::weak_ptr<Model> model_;
};

} // meso

} // ecell4

#endif /* __ECELL4_MESO_MESOSCOPIC_WORLD_HPP */
