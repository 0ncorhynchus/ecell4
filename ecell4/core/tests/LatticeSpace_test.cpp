#define BOOST_TEST_MODULE "LatticeSpace_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include <ecell4/core/MoleculePool.hpp>
#include <ecell4/core/VacantType.hpp>
#include <ecell4/core/LatticeSpaceVectorImpl.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>

using namespace ecell4;

struct Fixture
{
    const Real3 edge_lengths;
    const Real voxel_radius;
    LatticeSpaceVectorImpl space;
    SerialIDGenerator<ParticleID> sidgen;
    const Real D, radius;
    const Species sp;
    Fixture() :
        edge_lengths(2.5e-8, 2.5e-8, 2.5e-8),
        voxel_radius(2.5e-9),
        space(edge_lengths, voxel_radius, false),
        sidgen(), D(1e-12), radius(2.5e-9),
        sp("A", "2.5e-9", "1e-12")
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(suite, Fixture)

BOOST_AUTO_TEST_CASE(LatticeSpace_test_constructor)
{
    ;
}

BOOST_AUTO_TEST_CASE(CheckVacantSize)
{
    BOOST_CHECK_EQUAL(space.actual_size(), space.vacant()->size());
}

BOOST_AUTO_TEST_CASE(GetVoxel)
{
    const Real3 position(1.25e-8, 1.25e-8, 1.25e-8);
    const Integer coordinate(space.position2coordinate(position));

    {
        std::pair<ParticleID, Species> voxel(space.get_voxel_at(coordinate));
        BOOST_CHECK_EQUAL(voxel.first, ParticleID());
        BOOST_CHECK_EQUAL(voxel.second, space.vacant()->species());
    }

    ParticleID id(sidgen());
    BOOST_CHECK(space.update_voxel(id, ParticleVoxel(sp, coordinate, radius, D)));

    {
        std::pair<ParticleID, Species> voxel(space.get_voxel_at(coordinate));
        BOOST_CHECK_EQUAL(voxel.first, id);
        BOOST_CHECK_EQUAL(voxel.second, sp);
    }

    {
        boost::optional<ParticleVoxel> voxel(space.find_voxel(id));
        BOOST_ASSERT(voxel);
        BOOST_CHECK_EQUAL(voxel->species, sp);
    }
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_num_species)
{
    BOOST_CHECK_EQUAL(space.num_species(), 0);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_has_species)
{
    BOOST_CHECK(!space.has_species(sp));
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_update_particle)
{
    ParticleID id(sidgen());

    Real3 pos(2e-8, 1.7e-8, 1.5e-8);
    Real r(1.0);
    Real d(2.3);
    // Particle particle(sp, pos, r, d);
    ParticleVoxel v(sp, space.position2coordinate(pos), r, D);

    // BOOST_CHECK(space.update_particle(id, particle));
    BOOST_CHECK(space.update_voxel(id, v));
    BOOST_CHECK(space.has_species(sp));
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_num_particles)
{
    ParticleID id(sidgen());
    Real3 pos(2e-8, 1.7e-8, 1.5e-8);
    Real r(1.0), d(2.3);
    // Particle particle(sp, pos, r, d);
    ParticleVoxel v(sp, space.position2coordinate(pos), r, D);

    ParticleID a_id(sidgen());
    Species a(std::string("ANOTHER"));
    Real3 pos1(1e-8, 2e-8, 1e-9);
    Real r1(1.1);
    Real d1(4.3);
    // Particle another(a, pos1, r1, d1);
    ParticleVoxel another(a, space.position2coordinate(pos1), r1, d1);

    BOOST_CHECK(space.update_voxel(id, v));
    BOOST_CHECK(space.update_voxel(a_id, another));
    // BOOST_CHECK(space.update_particle(id, particle));
    // BOOST_CHECK(space.update_particle(a_id, another));
    BOOST_CHECK_EQUAL(space.num_voxels(sp), 1);
    BOOST_CHECK_EQUAL(space.num_voxels(), 2);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_list_particles)
{
    ParticleID id(sidgen());
    Real3 pos(2e-8, 1.7e-8, 1.5e-8);
    Real r(1.0), d(2.3);
    // Particle particle(sp, pos, r, d);
    ParticleVoxel v(sp, space.position2coordinate(pos), r, D);

    ParticleID a_id(sidgen());
    Species a(std::string("ANOTHER"));
    Real3 pos1(1e-8, 2e-8, 1e-9);
    Real r1(1.1);
    Real d1(4.3);
    // Particle another(a, pos1, r1, d1);
    ParticleVoxel another(a, space.position2coordinate(pos1), r1, d1);

    BOOST_CHECK(space.update_voxel(id, v));
    BOOST_CHECK(space.update_voxel(a_id, another));
    // BOOST_CHECK(space.update_particle(id, particle));
    // BOOST_CHECK(space.update_particle(a_id, another));

    BOOST_CHECK_EQUAL(space.list_voxels().size(), 2);
    BOOST_CHECK_EQUAL(space.list_voxels(sp).size(), 1);
}

// BOOST_AUTO_TEST_CASE(LatticeSpace_test_register_species)
// {
//     BOOST_CHECK(space.register_species(sp));
//     BOOST_CHECK(space.has_species(sp));
// 
//     std::vector<Species> list;
//     list.push_back(sp);
// 
//     BOOST_CHECK(list == space.list_species());
// }

BOOST_AUTO_TEST_CASE(LatticeSpace_test_coordinate_global_translation)
{
    for (Integer coord(0); coord < space.size(); ++coord)
    {
        const Integer3 global(space.coordinate2global(coord));
        VoxelSpaceBase::coordinate_type created_coord(space.global2coordinate(global));
        BOOST_CHECK_EQUAL(coord, created_coord);
    }
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_coordinate_position_translation)
{
    const Real3 origin_pos(space.coordinate2position(0));
    BOOST_ASSERT(origin_pos[0] < 0);
    BOOST_ASSERT(origin_pos[1] < 0);
    BOOST_ASSERT(origin_pos[2] < 0);

    const VoxelSpaceBase::coordinate_type origin(
                (space.col_size() + 3) * (space.row_size() + 2) + 1);
    const Real3 origin_p(space.coordinate2position(origin));
    BOOST_ASSERT(origin_p[0] == 0);
    BOOST_ASSERT(origin_p[1] == 0);
    BOOST_ASSERT(origin_p[2] == 0);

    BOOST_ASSERT(space.num_neighbors(origin) == 12);
    for (Integer i(0); i < 12; ++i)
    {
        const Real3 neighbor(
                space.coordinate2position(space.get_neighbor(origin, i)));
        BOOST_CHECK(origin_p != neighbor);
        const VoxelSpaceBase::coordinate_type coord(
                space.position2coordinate(origin_p * 0.7 + neighbor * 0.3));
        BOOST_CHECK_EQUAL(origin, coord);
    }

    Integer size(
            (space.col_size()+2) * (space.layer_size() + 2) * (space.row_size() + 2));
    for (VoxelSpaceBase::coordinate_type coord(0); coord < size; ++coord)
    {
        const Real3 pos(space.coordinate2position(coord));
        const Integer3 global(space.position2global(pos));
        const VoxelSpaceBase::coordinate_type created_coord(
                space.position2coordinate(pos));
        BOOST_CHECK_EQUAL(coord, created_coord);
    }
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_add_remove_molecule)
{
    const VoxelSpaceBase::coordinate_type coord(
            space.global2coordinate(Integer3(3,4,5)));
    ParticleID pid(sidgen());
    BOOST_CHECK(space.update_voxel(
        pid, ParticleVoxel(sp, coord, radius, D)));
    BOOST_CHECK_EQUAL(space.num_voxels(sp), 1);

    boost::shared_ptr<const VoxelPool> mt(space.get_voxel_pool_at(coord));
    BOOST_CHECK(!mt->is_vacant());

    BOOST_CHECK(space.remove_voxel(coord));
    boost::shared_ptr<const VoxelPool> vacant(space.get_voxel_pool_at(coord));
    BOOST_CHECK(vacant->is_vacant());
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_move)
{
    const Integer3 global0(2,3,4);
    const VoxelSpaceBase::coordinate_type coord(
            space.global2coordinate(global0));

    ParticleID pid(sidgen());
    BOOST_CHECK(space.update_voxel(
        pid, ParticleVoxel(sp, coord, radius, D)));

    boost::shared_ptr<VoxelPool> from_mt(space.get_voxel_pool_at(coord));
    BOOST_CHECK(!from_mt->is_vacant());

    const Integer3 global1(2,4,4);
    const VoxelSpaceBase::coordinate_type to_coord(
            space.global2coordinate(global1));

    BOOST_CHECK(space.move(coord, to_coord));

    boost::shared_ptr<VoxelPool> mt(space.get_voxel_pool_at(to_coord));
    BOOST_CHECK(!mt->is_vacant());

    BOOST_CHECK(space.update_voxel(
        sidgen(), ParticleVoxel(sp, coord, radius, D)));
    BOOST_CHECK(!space.move(coord, to_coord));
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_update_molecule)
{
    Species reactant(std::string("Reactant")),
            product(std::string("Product"));

    const Integer3 global(3,4,5);
    const VoxelSpaceBase::coordinate_type coord(
            space.global2coordinate(global));

    ParticleID pid(sidgen());
    BOOST_CHECK(space.update_voxel(
        pid, ParticleVoxel(reactant, coord, radius, D)));
    // space.update_voxel(
    //     ParticleVoxel(product, coord, radius, D));
    BOOST_CHECK(space.remove_voxel(coord));
    BOOST_CHECK(space.update_voxel(
        pid, ParticleVoxel(product, coord, radius, D)));

    boost::shared_ptr<const VoxelPool> mt(space.get_voxel_pool_at(coord));
    BOOST_ASSERT(mt->species() == product);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_update_voxel)
{
    const ParticleID pid(sidgen());
    for (VoxelSpaceBase::coordinate_type coord(0); coord < space.size(); ++coord)
    {
        if (!space.is_inside(coord))
        {
            continue;
        }

        space.update_voxel(pid, ParticleVoxel(sp, coord, radius, D));
        BOOST_CHECK_EQUAL(space.num_voxels(), 1);

        std::pair<ParticleID, ParticleVoxel> pair(space.list_voxels().at(0));
        BOOST_CHECK_EQUAL(pid, pair.first);
        BOOST_CHECK_EQUAL(coord, pair.second.coordinate);
        BOOST_CHECK_EQUAL(radius, pair.second.radius);
        BOOST_CHECK_EQUAL(D, pair.second.D);
    }
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_lattice_structure)
{
    for (VoxelSpaceBase::coordinate_type coord(0); coord < space.size(); ++coord)
    {
        if (space.is_inside(coord))
        {
            ParticleID pid(sidgen());
            BOOST_CHECK(space.update_voxel(pid, ParticleVoxel(sp, coord, radius, D)));
        }
    }
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_neighbor)
{
    for (VoxelSpaceBase::coordinate_type coord(0); coord < space.size(); ++coord)
    {
        if (!space.is_inside(coord))
        {
            continue;
        }

        BOOST_ASSERT(space.num_neighbors(coord) == 12);
        const Real3 center(space.coordinate2position(coord));
        for (int i(0); i < 12; ++i)
        {
            VoxelSpaceBase::coordinate_type neighbor(space.get_neighbor(coord, i));
            if (!space.is_inside(neighbor))
            {
                continue;
            }

            Real3 pos(space.coordinate2position(neighbor));
            Real3 vec((pos-center)/voxel_radius/2);
            Real r_ratio(length(pos-center)/voxel_radius/2);
            BOOST_ASSERT(r_ratio < 1.0001);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

struct PeriodicFixture
{
    const Real3 edge_lengths;
    const Real voxel_radius;
    LatticeSpaceVectorImpl space;
    SerialIDGenerator<ParticleID> sidgen;
    const Real D, radius;
    const Species sp;
    PeriodicFixture() :
        edge_lengths(2.5e-8, 2.5e-8, 2.5e-8),
        voxel_radius(2.5e-9),
        space(edge_lengths, voxel_radius, true),
        sidgen(), D(1e-12), radius(2.5e-9),
        sp(std::string("A"), "2.5e-9", "1e-12")
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(periodic_suite, PeriodicFixture)

BOOST_AUTO_TEST_CASE(LatticeSpace_test_periodic_col)
{
    std::cerr << " < periodic_col > ";
    const int col_size(space.col_size()),
              row_size(space.row_size()),
              layer_size(space.layer_size());
    for (int i(0); i < row_size; ++i)
        for (int j(0); j < layer_size; ++j)
        {
            const VoxelSpaceBase::coordinate_type
                coord(space.global2coordinate(Integer3(0, i, j)));

            BOOST_CHECK(space.update_voxel(sidgen(), ParticleVoxel(sp, coord, radius, D)));
        }

    // from 0 to col_size-1
    for (int i(0); i < row_size; ++i)
        for (int j(0); j < layer_size; ++j)
        {
            const VoxelSpaceBase::coordinate_type
                coord(space.global2coordinate(Integer3(0, i, j)));

            const Integer nrnd((j&1)==1?2:3);
            const VoxelSpaceBase::coordinate_type
                neighbor(space.get_neighbor(coord, nrnd));

            BOOST_CHECK_EQUAL(space.coordinate2global(neighbor).col, col_size-1);
            BOOST_CHECK(space.move(coord, neighbor));
        }

    // from col_size-1 to 0
    for (int i(0); i < row_size; ++i)
        for (int j(0); j < layer_size; ++j)
        {
            const VoxelSpaceBase::coordinate_type
                coord(space.global2coordinate(Integer3(col_size-1, i, j)));

            const Integer nrnd((j&1)==1?4:5);
            const VoxelSpaceBase::coordinate_type
                neighbor(space.get_neighbor(coord, nrnd));

            BOOST_CHECK_EQUAL(space.coordinate2global(neighbor).col, 0);
            BOOST_CHECK(space.move(coord, neighbor));
        }
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_periodic_row)
{
    const int col_size(space.col_size()),
              row_size(space.row_size()),
              layer_size(space.layer_size());
    for (int layer(0); layer < layer_size; ++layer)
        for (int col(0); col < col_size; ++col)
        {
            const VoxelSpaceBase::coordinate_type
                coord(space.global2coordinate(Integer3(col, 0, layer)));

            BOOST_CHECK(space.update_voxel(sidgen(), ParticleVoxel(sp, coord, radius, D)));
        }

    // from 0 to row_size-1
    for (int layer(0); layer < layer_size; ++layer)
        for (int col(0); col < col_size; ++col)
        {
            const VoxelSpaceBase::coordinate_type
                coord(space.global2coordinate(Integer3(col, 0, layer)));

            const Integer nrnd(0);
            const VoxelSpaceBase::coordinate_type
                neighbor(space.get_neighbor(coord, nrnd));

            BOOST_CHECK_EQUAL(space.coordinate2global(neighbor).row, row_size-1);
            BOOST_CHECK(space.move(coord, neighbor));
        }
    // from row_size-1 to 0
    for (int layer(0); layer < layer_size; ++layer)
        for (int col(0); col < col_size; ++col)
        {
            const VoxelSpaceBase::coordinate_type coord(
                    space.global2coordinate(Integer3(col, row_size-1, layer)));
            const Integer nrnd(1);
            const VoxelSpaceBase::coordinate_type
                neighbor(space.get_neighbor(coord, nrnd));

            BOOST_CHECK_EQUAL(space.coordinate2global(neighbor).row, 0);
            BOOST_CHECK(space.move(coord, neighbor));
        }
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_periodic_layer)
{
    const int col_size(space.col_size()),
              row_size(space.row_size()),
              layer_size(space.layer_size());

    for (int row(0); row < row_size; ++row)
        for (int col(0); col < col_size; ++col)
        {
            const VoxelSpaceBase::coordinate_type
                coord(space.global2coordinate(Integer3(col, row, 0)));

            BOOST_CHECK(space.update_voxel(sidgen(), ParticleVoxel(sp, coord, radius, D)));
        }

    // from 0 to layer_size-1
    for (int row(0); row < row_size; ++row)
        for (int col(0); col < col_size; ++col)
        {
            const VoxelSpaceBase::coordinate_type
                coord(space.global2coordinate(Integer3(col, row, 0)));

            const Integer nrnd((col&1)==1?8:9);
            const VoxelSpaceBase::coordinate_type
                neighbor(space.get_neighbor(coord, nrnd));

            BOOST_CHECK_EQUAL(space.coordinate2global(neighbor).layer, layer_size-1);
            BOOST_CHECK(space.move(coord, neighbor));
        }

    // from layer_size-1 to 0
    for (int row(0); row < row_size; ++row)
        for (int col(0); col < col_size; ++col)
        {
            const VoxelSpaceBase::coordinate_type coord(
                    space.global2coordinate(Integer3(col, row, layer_size-1)));
            const Integer nrnd((col&1)==1?10:11);
            const VoxelSpaceBase::coordinate_type
                neighbor(space.get_neighbor(coord, nrnd));

            BOOST_CHECK_EQUAL(space.coordinate2global(neighbor).layer, 0);
            BOOST_CHECK(space.move(coord, neighbor));
        }
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_coordinates2)
{
    const Integer3 g1(4, 4, 4);
    const VoxelSpaceBase::coordinate_type pc1(space.global2coordinate(g1));
    const Integer3 g3(space.coordinate2global(pc1));

    BOOST_CHECK(g1.col == g3.col && g1.row == g3.row && g1.layer == g3.layer);

    const Real3 p1(space.global2position(g1));
    const Integer3 g4(space.position2global(p1));

    BOOST_CHECK(g1.col == g4.col && g1.row == g4.row && g1.layer == g4.layer);
    BOOST_CHECK_EQUAL(pc1, space.position2coordinate(p1));

    const Real3 p2(space.coordinate2position(pc1));
    BOOST_CHECK_EQUAL(pc1, space.position2coordinate(p2));
}

BOOST_AUTO_TEST_SUITE_END()
