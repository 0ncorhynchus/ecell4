#include <iostream>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/Species.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/ode/ODESimulator.hpp>

#include <ecell4/ode/ODENetworkModel.hpp>
#include <ecell4/ode/ODEReactionRule.hpp>
#include <ecell4/ode/ODESimulator_New.hpp>

#include <boost/format.hpp>


using namespace ecell4;
using namespace ecell4::ode;

/**
 * main function
 */
int main(int argc, char** argv)
{
    const Real L(1e-6);
    const Real3 edge_lengths(L, L, L);
    const Real volume(L * L * L);
    const Real N(60);
    const Real ka(0.1), U(0.5);

    Species sp1("A"), sp2("B"), sp3("C");
    ReactionRule rr1, rr2;
    rr1.set_k(ka);
    rr1.add_reactant(sp1);
    rr1.add_product(sp2);
    rr1.add_product(sp3);

    const Real kd(ka * volume * (1 - U) / (U * U * N));
    rr2.set_k(kd);
    rr2.add_reactant(sp2);
    rr2.add_reactant(sp3);
    rr2.add_product(sp1);

    ReactionRule rr3;
    rr3.add_reactant(sp1);
    rr3.add_product(sp3);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species_attribute(sp1);
    model->add_species_attribute(sp2);
    model->add_species_attribute(sp3);
    model->add_reaction_rule(rr1);
    model->add_reaction_rule(rr2);
    model->add_reaction_rule(rr3);
    //model->dump_reactions();
    boost::shared_ptr<ODENetworkModel>  ode_model(new ODENetworkModel(model) );//XXX OLD one!

    boost::shared_ptr<ODEWorld> world(new ODEWorld(edge_lengths));
    boost::shared_ptr<ODEWorld_New> new_world(new ODEWorld_New(edge_lengths));

    world->add_molecules(sp1, N);
    new_world->add_molecules(sp1, N);

    //ODESimulator target(model, world, RUNGE_KUTTA_CASH_KARP54);
    ODESimulator target(model, world, EULER);
    target.initialize();

    //ODESimulator_New new_target(model, new_world, RUNGE_KUTTA_CASH_KARP54);
    ODESimulator_New new_target(model, new_world, EULER);
    new_target.initialize();

    Real next_time(0.0), dt(0.01);
    std::cout << boost::format("%6s  %6s %6s %6s     %6s %6s %6s %s") % "time" % "A_old" % "B_old" % "C_old" % "A_new" % "B_new" % "C_new" % "BOOL" << std::endl;
    std::cout << boost::format("%6.f  %6d %6d %6d     %6d %6d %6d %s") % 
        target.t() % 
        world->num_molecules(sp1) % world->num_molecules(sp2) % world->num_molecules(sp3) %
        new_world->num_molecules(sp1) % new_world->num_molecules(sp2) % new_world->num_molecules(sp3) % 
        ( world->num_molecules(sp1) == new_world->num_molecules(sp1) && 
          world->num_molecules(sp2) == new_world->num_molecules(sp2) && 
          world->num_molecules(sp3) == new_world->num_molecules(sp3) )   ;

    //std::cout << target.t()
    //          << "\t" << world->num_molecules(sp1)
    //          << "\t" << world->num_molecules(sp2)
    //          << "\t" << world->num_molecules(sp3)
    //          << std::endl;
    
    bool ok_flag = true;
    for (unsigned int i(0); i < 1000; ++i)
    {
        next_time += dt;
        target.step(next_time);
        new_target.step(next_time);
        if (target.t() != new_target.t() ) {
            throw;
        }
        //std::cout << target.t()
        //          << "\t" << world->num_molecules(sp1)
        //          << "\t" << world->num_molecules(sp2)
        //          << "\t" << world->num_molecules(sp3)
        //          << std::endl;
        bool flag = (world->num_molecules(sp1) == new_world->num_molecules(sp1)) && ( world->num_molecules(sp2) == new_world->num_molecules(sp2) ) && ( world->num_molecules(sp3) == new_world->num_molecules(sp3) );
        std::cout << boost::format("%6f  %6d %6d %6d     %6d %6d %6d %s\n") % 
            target.t() % 
            world->num_molecules(sp1) % world->num_molecules(sp2) % world->num_molecules(sp3) %
            new_world->num_molecules(sp1) % new_world->num_molecules(sp2) % new_world->num_molecules(sp3) % (flag == true ? "TRUE" : "FALSE" );
        if (flag != true) {
            ok_flag = false;
        }
    }
    if (ok_flag) {
        std::cout << "Exactly the same result with new and old ODESimulator" << std::endl;
    }
}
