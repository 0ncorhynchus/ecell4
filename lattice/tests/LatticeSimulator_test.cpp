#define BOOST_TEST_MODULE "LatticeSimulator_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../LatticeSimulator.hpp"

using namespace ecell4;
using namespace ecell4::lattice;

BOOST_AUTO_TEST_CASE(LatticeSimulator_test_constructor)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(1e-8);

    const std::string D("1e-12"), radius("2.5e-13");

    ecell4::Species sp1("A", radius, D),
        sp2("B", radius, D),
        sp3("C", radius, D);
    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    (*model).add_species(sp1);
    (*model).add_species(sp2);
    (*model).add_species(sp3);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, rng));

    LatticeSimulator sim(model, world);
}

BOOST_AUTO_TEST_CASE(LatticeSimulator_test_hdf5_save)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(1e-8);
    const Integer N(60);

    const std::string D("1e-12"), radius("2.5e-9");

    ecell4::Species sp1("A", radius, D);
    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    (*model).add_species(sp1);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, rng));

    world->add_molecules(sp1, N / 2);

    BOOST_ASSERT(world->num_molecules(sp1) == N / 2);

    LatticeSimulator sim(model, world);

    world->add_molecules(sp1, N / 2);
    BOOST_ASSERT(world->num_molecules(sp1) == N);

    H5::H5File fout("data.h5", H5F_ACC_TRUNC);
    const std::string hdf5path("/");
    world->save(&fout, hdf5path);
}

BOOST_AUTO_TEST_CASE(LatticeSimulator_test_step_with_single_species)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(1e-8);
    const Integer N(60);

    const std::string D("1e-12"), radius("2.5e-13");

    ecell4::Species sp1("A", radius, D);
    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    (*model).add_species(sp1);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, rng));

    world->add_molecules(sp1, N / 2);

    BOOST_ASSERT(world->num_molecules(sp1) == N / 2);

    LatticeSimulator sim(model, world);

    world->add_molecules(sp1, N / 2);
    BOOST_ASSERT(world->num_molecules(sp1) == N);

    sim.step();
}

BOOST_AUTO_TEST_CASE(LatticeSimulator_test_reaction)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(1e-8);
    const std::string radius("2.5e-13");
    const ecell4::Species sp1("A", radius, "1.0e-12"),
          sp2("B", radius, "1.1e-12");

    //const ecell4::ReactionRule rr1();

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species(sp1);
    model->add_species(sp2);
    //model->add_reaction_rule(rr1);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, rng));

    LatticeSimulator sim(model, world);
}

BOOST_AUTO_TEST_CASE(LattiecSimulator_test_scheduler)
{
    std::cout << " <<LattiecSimulator_test_scheduler>> ";
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(1e-8);

    const std::string D1("1.0e-12"),
          D2("1.1e-12"),
          D3("1.2e-12"),
          radius("2.5e-13");

    const ecell4::Species sp1("A", radius, D1),
        sp2("B", radius, D2),
        sp3("C", radius, D3);
    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    (*model).add_species(sp1);
    (*model).add_species(sp2);
    (*model).add_species(sp3);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, rng));

    Coord c1(world->global2coord(Global(40,34,56))),
          c2(world->global2coord(Global(32,50,24))),
          c3(world->global2coord(Global(60,36,89)));
    BOOST_CHECK(world->add_molecule(sp1, c1));
    BOOST_CHECK(world->add_molecule(sp2, c2));
    BOOST_CHECK(world->add_molecule(sp3, c3));

    LatticeSimulator sim(model, world);

    sim.initialize();

    const MolecularTypeBase
        *mt1(world->get_molecular_type(sp1)),
        *mt2(world->get_molecular_type(sp2)),
        *mt3(world->get_molecular_type(sp3));
    std::vector<std::pair<Coord, ParticleID> >::const_iterator
        itr1(mt1->begin()),
        itr2(mt2->begin()),
        itr3(mt3->begin());

    BOOST_ASSERT(itr1 != mt1->end());
    BOOST_ASSERT(itr2 != mt2->end());
    BOOST_ASSERT(itr3 != mt3->end());

    c1 = (*itr1).first;
    c2 = (*itr2).first;
    c3 = (*itr3).first;

    sim.step();
    itr1 = mt1->begin();
    itr2 = mt2->begin();
    itr3 = mt3->begin();
    std::cout << "<itr1: " << (*itr1).first << "> ";
    std::cout << "<itr2: " << (*itr2).first << "> ";
    std::cout << "<itr3: " << (*itr3).first << "> ";
    BOOST_ASSERT((*itr1).first == c1);
    BOOST_ASSERT((*itr2).first == c2);
    BOOST_ASSERT((*itr3).first != c3);
    c3 = (*itr3).first;

    sim.step();
    itr1 = mt1->begin();
    itr2 = mt2->begin();
    itr3 = mt3->begin();
    BOOST_ASSERT((*itr1).first == c1);
    BOOST_ASSERT((*itr2).first != c2);
    BOOST_ASSERT((*itr3).first == c3);
    c2 = (*itr2).first;

    sim.step();
    itr1 = mt1->begin();
    itr2 = mt2->begin();
    itr3 = mt3->begin();
    BOOST_ASSERT((*itr1).first != c1);
    BOOST_ASSERT((*itr2).first == c2);
    BOOST_ASSERT((*itr3).first == c3);
    c1 = (*itr1).first;

}

