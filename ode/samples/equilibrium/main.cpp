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
    Real const volume(1e-18);
    Real const N(60);
    Real const ka(0.1), U(0.5);

    Species sp1("A"), sp2("B"), sp3("C");
    ReactionRule rr1, rr2;
    rr1.set_k(ka);
    rr1.add_reactant(sp1);
    rr1.add_product(sp2);
    rr1.add_product(sp3);
    Real const kd(ka * volume * (1 - U) / (U * U * N));
    rr2.set_k(kd);
    rr2.add_reactant(sp2);
    rr2.add_reactant(sp3);
    rr2.add_product(sp1);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species(sp1);
    model->add_species(sp2);
    model->add_species(sp3);
    model->add_reaction_rule(rr1);
    model->add_reaction_rule(rr2);

    boost::shared_ptr<ODEWorld> world(new ODEWorld(volume));
    world->add_species(sp1);
    world->add_species(sp2);
    world->add_species(sp3);
    world->add_molecules(sp1, N);

    ODESimulator target(model, world);

    Real next_time(0.0), dt(0.01);
    std::cout << target.t()
              << "\t" << world->num_molecules(sp1)
              << "\t" << world->num_molecules(sp2)
              << "\t" << world->num_molecules(sp3)
              << std::endl;
    for (unsigned int i(0); i < 1000; ++i)
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
