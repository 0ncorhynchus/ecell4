#include <iostream>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/Species.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include "../../ODESimulator.hpp"

using namespace ecell4;
using namespace ecell4::ode;

/**
 * main function
 */
int main(int argc, char** argv)
{
    const Real volume(1e-18);

    Species sp1("A"), sp2("B"), sp3("C");
    ReactionRule rr1;
    rr1.set_k(1.0);
    rr1.add_reactant(sp1);
    rr1.add_product(sp2);
    rr1.add_product(sp3);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species(sp1);
    model->add_species(sp2);
    model->add_species(sp3);
    model->add_reaction_rule(rr1);

    boost::shared_ptr<ODEWorld> world(new ODEWorld(volume));
    world->add_species(sp1);
    world->add_species(sp2);
    world->add_species(sp3);
    world->add_molecules(sp1, 60);

    ODESimulator target(model, world);

	// ecell4_hdf5_manager<ODEWorld, double> hdf("dissociation.hdf5", model, world, "ODEWorld");

    Real next_time(0.0), dt(0.01);
	// target.save_hdf5_init(std::string("hogehoge.hdf5"));
	// ecell4_hdf5_manager<ODEWorld, double> hdf5_mng(std::string("ode_test.hdf5"), model, world, "MyOdeWorld");

    std::cout << target.t()
              << "\t" << world->num_molecules(sp1)
              << "\t" << world->num_molecules(sp2)
              << "\t" << world->num_molecules(sp3)
              << std::endl;
    for (unsigned int i(0); i < 200; ++i)
    {
        next_time += dt;
        target.step(next_time);
        std::cout << target.t()
                  << "\t" << world->num_molecules(sp1)
                  << "\t" << world->num_molecules(sp2)
                  << "\t" << world->num_molecules(sp3)
                  << std::endl;
    }
}
