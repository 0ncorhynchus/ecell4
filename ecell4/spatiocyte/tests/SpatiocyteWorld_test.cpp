#define BOOST_TEST_MODULE "SpatiocyteWorld_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Sphere.hpp>

#include "../SpatiocyteWorld.hpp"

using namespace ecell4;
using namespace ecell4::spatiocyte;

const Real3 EDGE_LENGTHS(1e-6, 1e-6, 1e-6);
const Real VOXEL_RADIUS = 1e-8;

BOOST_AUTO_TEST_CASE(ConstructorTest)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SpatiocyteWorld world(EDGE_LENGTHS, VOXEL_RADIUS, rng);
}

BOOST_AUTO_TEST_CASE(TimeTest)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SpatiocyteWorld world(EDGE_LENGTHS, VOXEL_RADIUS, rng);
    BOOST_CHECK_EQUAL(world.t(), 0);
    world.set_t(23.4);
    BOOST_CHECK_EQUAL(world.t(), 23.4);
}

BOOST_AUTO_TEST_CASE(ModelTest)
{
    boost::shared_ptr<GSLRandomNumberGenerator> rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<NetworkModel> model(new NetworkModel());

    SpatiocyteWorld world(EDGE_LENGTHS, VOXEL_RADIUS, rng);
    world.bind_to(model);

    model->add_species_attribute(Species("A"));
    MoleculeInfo minfoA = world.get_molecule_info(Species("A"));
    BOOST_CHECK_EQUAL(minfoA.radius, VOXEL_RADIUS);
    BOOST_CHECK_EQUAL(minfoA.D, 0.0);
    BOOST_CHECK_EQUAL(minfoA.loc, "");
    BOOST_CHECK_EQUAL(minfoA.dimension, Shape::THREE);
    BOOST_CHECK(!minfoA.is_structure);

    model->add_species_attribute(Species("B", "2.5e-9", "1e-12"));
    MoleculeInfo minfoB = world.get_molecule_info(Species("B"));
    BOOST_CHECK_EQUAL(minfoB.radius, 2.5e-9);
    BOOST_CHECK_EQUAL(minfoB.D, 1e-12);
    BOOST_CHECK_EQUAL(minfoB.loc, "");
    BOOST_CHECK_EQUAL(minfoB.dimension, Shape::THREE);
    BOOST_CHECK(!minfoB.is_structure);

    Species speciesC("C");
    speciesC.set_attribute<Integer>("dimension", 2);
    model->add_species_attribute(speciesC);
    Species speciesD("D");
    speciesD.set_attribute("D", 1e-12);
    speciesD.set_attribute("location", "C");
    model->add_species_attribute(speciesD);
    MoleculeInfo minfoD = world.get_molecule_info(Species("D"));
    BOOST_CHECK_EQUAL(minfoD.radius, VOXEL_RADIUS);
    BOOST_CHECK_EQUAL(minfoD.D, 1e-12);
    BOOST_CHECK_EQUAL(minfoD.loc, "C");
    BOOST_CHECK_EQUAL(minfoD.dimension, Shape::TWO);
    BOOST_CHECK(!minfoD.is_structure);

    MoleculeInfo minfoC = world.get_molecule_info(Species("C"));
    BOOST_CHECK_EQUAL(minfoC.radius, VOXEL_RADIUS);
    BOOST_CHECK_EQUAL(minfoC.D, 0.0);
    BOOST_CHECK_EQUAL(minfoC.loc, "");
    BOOST_CHECK_EQUAL(minfoC.dimension, Shape::TWO);
    BOOST_CHECK(minfoC.is_structure);
}

BOOST_AUTO_TEST_CASE(SpeciesTest)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SpatiocyteWorld world(EDGE_LENGTHS, VOXEL_RADIUS, rng);

    Species sp(std::string("A"));

    BOOST_CHECK_EQUAL(world.list_species().size(), 0);
    BOOST_CHECK(!world.has_species(sp));
}

BOOST_AUTO_TEST_CASE(ListParticlesTest)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SpatiocyteWorld world(EDGE_LENGTHS, VOXEL_RADIUS, rng);
    std::vector<std::pair<ParticleID, Particle> > particles(world.list_particles());
    BOOST_CHECK_EQUAL(particles.size(), 0);
}

BOOST_AUTO_TEST_CASE(UpdateParticlesTest)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SerialIDGenerator<ParticleID> sidgen;
    SpatiocyteWorld world(EDGE_LENGTHS, VOXEL_RADIUS, rng);

    ParticleID pid(sidgen());
    Species sp(std::string("A"));
    const Real3 pos(2e-7, 1e-7, 0);
    Real r(0);
    Real d(0);
    Particle p(sp, pos, r, d);

    world.update_particle(pid, p);

    BOOST_CHECK(world.has_species(sp));
    BOOST_CHECK(world.has_particle(pid));
    BOOST_CHECK_EQUAL(world.list_particles().size(), 1);
    BOOST_CHECK_EQUAL(world.list_particles(sp).size(), 1);
}

BOOST_AUTO_TEST_CASE(AddMoleculeTest)
{
    const Real voxel_radius(2.5e-9);
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SpatiocyteWorld world(EDGE_LENGTHS, voxel_radius, rng);

    Species sp(std::string("TEST"));
    sp.set_attribute("radius", "2.5e-9");
    sp.set_attribute("D", "1e-12");

    const Voxel voxel(world.position2voxel(EDGE_LENGTHS / 2.0));
    // BOOST_CHECK(world.place_voxel(sp, coord).second);
    BOOST_CHECK(world.new_voxel(sp, voxel));
    BOOST_CHECK_EQUAL(world.num_particles(sp), 1);

    boost::shared_ptr<const VoxelPool> mt(voxel.get_voxel_pool());
    BOOST_CHECK(!mt->is_vacant());
}

BOOST_AUTO_TEST_CASE(AddMultipleMoleculesTest)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SpatiocyteWorld world(EDGE_LENGTHS, VOXEL_RADIUS, rng);

    Species sp(std::string("TEST"));
    sp.set_attribute("radius", "2.5e-9");
    sp.set_attribute("D", "1e-12");
    const Integer N(60);

    BOOST_CHECK(world.add_molecules(sp, N));
    BOOST_CHECK_EQUAL(world.num_particles(sp), N);
}

BOOST_AUTO_TEST_CASE(NeighborTest)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SpatiocyteWorld world(EDGE_LENGTHS, VOXEL_RADIUS, rng);

    const Voxel voxel(world.position2voxel(EDGE_LENGTHS / 2.0));
    const Real3 cp(voxel.position());

    Species sp(std::string("TEST"));
    sp.set_attribute("radius", "2.5e-9");
    sp.set_attribute("D", "1e-12");
    for (Integer i(0); i < voxel.num_neighbors(); ++i)
    {
        world.new_voxel(sp, voxel.get_neighbor(i));
    }
    std::vector<std::pair<ParticleID, Particle> > particles(world.list_particles());
    for (std::vector<std::pair<ParticleID, Particle> >::iterator itr(
                particles.begin()); itr != particles.end(); ++itr)
    {
        Real3 pos((*itr).second.position());
        BOOST_ASSERT(length(pos-cp) < VOXEL_RADIUS*2.1);
    }

#ifdef WITH_HDF5
    world.save("neighbor.h5");
#endif
}

BOOST_AUTO_TEST_CASE(AddShapeTest)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SpatiocyteWorld world(EDGE_LENGTHS, VOXEL_RADIUS, rng);

    Species sp(std::string("TEST"));
    sp.set_attribute("radius", "2.5e-9");
    sp.set_attribute("D", "1e-12");

    boost::shared_ptr<const Sphere> sphere(new Sphere(Real3(5e-7, 5e-7, 5e-7), 5e-7*1.5));

    const Integer n(world.add_structure(sp, sphere));
    BOOST_ASSERT(n > 0);
    BOOST_CHECK_EQUAL(world.num_particles(sp), n);

#ifdef WITH_HDF5
    world.save("sphere.h5");
#endif
}

BOOST_AUTO_TEST_CASE(MoveTest)
{
    const Real voxel_radius(2.5e-9);
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SpatiocyteWorld world(EDGE_LENGTHS, voxel_radius, rng);

    Species sp(std::string("TEST"));
    sp.set_attribute("radius", "2.5e-9");
    sp.set_attribute("D", "1e-12");

    const Voxel from(world.position2voxel(Real3(0.3e-6, 0.5e-6, 0.5e-6)));
    const Voxel to(world.position2voxel(Real3(0.5e-6, 0.5e-6, 0.5e-6)));

    BOOST_CHECK(world.new_voxel(sp, from));
    BOOST_CHECK(world.move(from, to));

    boost::shared_ptr<const VoxelPool> mt(to.get_voxel_pool());
    BOOST_CHECK(!mt->is_vacant());

    BOOST_CHECK(world.move(from, to));
}

BOOST_AUTO_TEST_CASE(StructureTest)
{
    const Real3 edge_lengths(5e-7, 5e-7, 5e-7);
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SpatiocyteWorld world(edge_lengths, VOXEL_RADIUS, rng);

    Species membrane("Membrane", "2.5e-9", "0");

    Species sp("SpeciesA", "2.5e-9", "1e-12");
    sp.set_attribute("location", "Membrane");

    boost::shared_ptr<const Sphere> sphere(new Sphere(Real3(2.5e-7, 2.5e-7, 2.5e-7), 2e-7));

    BOOST_CHECK(world.add_structure(membrane, sphere) == 5892);
    BOOST_CHECK(!world.new_particle(Particle(sp, Real3(2.5e-7, 2.5e-7, 4.5e-7), 2.5e-9, 1e-12)));
    BOOST_CHECK(world.new_particle(Particle(sp, Real3(2.5e-7, 2.5e-7, 4.5e-7 - VOXEL_RADIUS * 2), 2.5e-9, 1e-12)));

#ifdef WITH_HDF5
    world.save("structure.h5");
#endif
}
