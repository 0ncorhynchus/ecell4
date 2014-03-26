#ifndef __ECELL4_PARTICLE_SPACE_HDF5_WRITER_HPP
#define __ECELL4_PARTICLE_SPACE_HDF5_WRITER_HPP

#include <cstring>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

#include <hdf5.h>
#include <H5Cpp.h>

#include "types.hpp"
#include "get_mapper_mf.hpp"
#include "Species.hpp"
#include "Particle.hpp"


namespace ecell4
{

template<typename Tspace_>
class ParticleSpaceHDF5Writer
{
public:

    typedef Tspace_ space_type;

protected:

    struct h5_species_struct {
        uint32_t id;
        char serial[32]; // species' serial may exceed the limit
    };

    struct h5_particle_struct {
        int lot;
        int serial;
        uint32_t sid;
        double posx;
        double posy;
        double posz;
        double radius;
        double D;
    };

public:

    ParticleSpaceHDF5Writer(const space_type& space)
        : space_(space)
    {
        ;
    }

    virtual ~ParticleSpaceHDF5Writer()
    {
        ;
    }

    void save(H5::Group* root) const
    {
        using namespace H5;

        typedef std::vector<std::pair<ParticleID, Particle> >
        particle_container_type;
        const particle_container_type& particles(space_.list_particles());
        const unsigned int num_particles(particles.size());

        std::vector<Species> species;
        typedef utils::get_mapper_mf<Species::serial_type, unsigned int>::type
        species_id_map_type;
        species_id_map_type species_id_map;

        boost::scoped_array<h5_particle_struct>
            h5_particle_table(new h5_particle_struct[num_particles]);
        for (unsigned int i(0); i < num_particles; ++i)
        {
            species_id_map_type::const_iterator
                it(species_id_map.find(particles[i].second.species().serial()));
            if (it == species_id_map.end())
            {
                species.push_back(particles[i].second.species());
                it = species_id_map.insert(
                    std::make_pair(particles[i].second.species().serial(),
                                   species.size())).first;
            }

            h5_particle_table[i].lot = particles[i].first.lot();
            h5_particle_table[i].serial = particles[i].first.serial();
            h5_particle_table[i].sid = (*it).second;
            h5_particle_table[i].posx = particles[i].second.position()[0];
            h5_particle_table[i].posy = particles[i].second.position()[1];
            h5_particle_table[i].posz = particles[i].second.position()[2];
            h5_particle_table[i].radius = particles[i].second.radius();
            h5_particle_table[i].D = particles[i].second.D();
        }

        boost::scoped_array<h5_species_struct>
            h5_species_table(new h5_species_struct[species.size()]);
        for (unsigned int i(0); i < species.size(); ++i)
        {
            h5_species_table[i].id = i + 1;
            std::strcpy(h5_species_table[i].serial,
                        species[i].serial().c_str());
        }

        CompType h5_particle_comp_type(sizeof(h5_particle_struct));
        h5_particle_comp_type.insertMember(
            std::string("lot"), HOFFSET(h5_particle_struct, lot),
            PredType::NATIVE_INT);
        h5_particle_comp_type.insertMember(
            std::string("serial"), HOFFSET(h5_particle_struct, serial),
            PredType::NATIVE_INT);
        h5_particle_comp_type.insertMember(
            std::string("sid"), HOFFSET(h5_particle_struct, sid),
            PredType::STD_I32LE);
        h5_particle_comp_type.insertMember(
            std::string("posx"), HOFFSET(h5_particle_struct, posx),
            PredType::NATIVE_DOUBLE);
        h5_particle_comp_type.insertMember(
            std::string("posy"), HOFFSET(h5_particle_struct, posy),
            PredType::NATIVE_DOUBLE);
        h5_particle_comp_type.insertMember(
            std::string("posz"), HOFFSET(h5_particle_struct, posz),
            PredType::NATIVE_DOUBLE);
        h5_particle_comp_type.insertMember(
            std::string("radius"), HOFFSET(h5_particle_struct, radius),
            PredType::NATIVE_DOUBLE);
        h5_particle_comp_type.insertMember(
            std::string("D"), HOFFSET(h5_particle_struct, D),
            PredType::NATIVE_DOUBLE);

        // const hsize_t dims[] = {3};
        // mtype.insertMember(
        //     MEMBER3, HOFFSET(h5_particles, h5_particle_position),
        //     ArrayType(PredType::NATIVE_DOUBLE, 1, dims));

        CompType h5_species_comp_type(sizeof(h5_species_struct));
        h5_species_comp_type.insertMember(
            std::string("id"), HOFFSET(h5_species_struct, id),
            PredType::STD_I32LE);
        h5_species_comp_type.insertMember(
            std::string("serial"), HOFFSET(h5_species_struct, serial),
            StrType(PredType::C_S1, 32));

        const int RANK = 1;
        hsize_t dim1[] = {num_particles};
        DataSpace dataspace1(RANK, dim1);
        boost::scoped_ptr<DataSet> dataset1(new DataSet(
            root->createDataSet("particles", h5_particle_comp_type, dataspace1)));

        hsize_t dim2[] = {species.size()};
        DataSpace dataspace2(RANK, dim2);
        boost::scoped_ptr<DataSet> dataset2(new DataSet(
            root->createDataSet("species", h5_species_comp_type, dataspace2)));

        dataset1->write(h5_particle_table.get(), h5_particle_comp_type);
        dataset2->write(h5_species_table.get(), h5_species_comp_type);

        const double t = space_.t();
        Attribute attr_t(
            root->createAttribute(
                "t", PredType::IEEE_F64LE, DataSpace(H5S_SCALAR)));
        attr_t.write(PredType::IEEE_F64LE, &t);

        const Position3 edge_lengths = space_.edge_lengths();
        const hsize_t dims[] = {3};
        const ArrayType lengths_type(PredType::NATIVE_DOUBLE, 1, dims);
        Attribute attr_lengths(
            root->createAttribute(
                "edge_lengths", lengths_type, DataSpace(H5S_SCALAR)));
        double lengths[] = {edge_lengths[0], edge_lengths[1], edge_lengths[2]};
        attr_lengths.write(lengths_type, lengths);
    }

protected:

    const space_type& space_;
};

} // ecell4

#endif /*  __ECELL4_PARTICLE_SPACE_HDF5_WRITER_HPP */
