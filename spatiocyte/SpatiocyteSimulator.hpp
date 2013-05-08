#ifndef __ECELL4_SPATIOCYTE_SPATIOCYTE_SIMULATOR_HPP
#define __ECELL4_SPATIOCYTE_SPATIOCYTE_SIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Simulator.hpp>

#include "SpatiocyteWorld.hpp"

#include <H5Cpp.h>

namespace ecell4
{

namespace spatiocyte
{

class SpatiocyteSimulator
    : public Simulator
{
public:

    SpatiocyteSimulator(
        boost::shared_ptr<NetworkModel> model,
        boost::shared_ptr<SpatiocyteWorld> world)
        : model_(model), world_(world), num_steps_(0)
    {
        const NetworkModel::species_container_type&
            species((*model_).species());
        for (NetworkModel::species_container_type::const_iterator
                 i(species.begin()); i != species.end(); ++i)
        {
            (*world_).add_species(*i);
        }

        const NetworkModel::reaction_rule_container_type&
            reaction_rules((*model_).reaction_rules());
        for (NetworkModel::reaction_rule_container_type::const_iterator
                 i(reaction_rules.begin()); i != reaction_rules.end(); ++i)
        {
            (*world_).add_reaction_rule(*i);
        }

        initialize();
    }

    // SimulatorTraits

    Real t() const
    {
        return (*world_).t();
    }

    Real dt() const
    {
        return (*world_).dt();
    }

    Integer num_steps() const
    {
        return num_steps_;
    }

    void step();
    bool step(const Real& upto);

    // Optional members

    void initialize()
    {
        (*world_).initialize();
    }

    SpatiocyteStepper* spatiocyte_stepper() const
    {
        return (*world_).spatiocyte_stepper();
    }

    void save_hdf5_init(std::string);
    void save_hdf5(void);

protected:

    boost::shared_ptr<NetworkModel> model_;
    boost::shared_ptr<SpatiocyteWorld> world_;

    H5::H5File *file_;
    typedef struct h5_lattice {
        int h5_lattice_id;
        char h5_species_id[32];
        // double h5_particle_position[3];
    } h5_lattice;

    /**
     * the protected internal state of SpatiocyteSimulator.
     * they are needed to be saved/loaded with Visitor pattern.
     */

    Integer num_steps_;
};

} // spatiocyte

} // ecell4

#endif /* __ECELL4_SPATIOCYTE_SPATIOCYTE_SIMULATOR_HPP */
