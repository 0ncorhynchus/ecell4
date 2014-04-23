#ifndef __ECELL4_ODE_ODE_WORLD_HPP
#define __ECELL4_ODE_ODE_WORLD_HPP

#include <ecell4/core/Species.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/Space.hpp>
#include <ecell4/core/CompartmentSpaceHDF5Writer.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

namespace ecell4
{

namespace ode
{

class ODEWorld
{
protected:

    typedef std::vector<Real> num_molecules_container_type;
    typedef std::vector<Species> species_container_type;
    typedef utils::get_mapper_mf<
        Species, num_molecules_container_type::size_type>::type species_map_type;

public:

    ODEWorld(const Position3& edge_lengths)
        : volume_(1.0), t_(0.0)
    {
        const Real volume(edge_lengths[0] * edge_lengths[1] * edge_lengths[2]);
        set_volume(volume);
    }

    // SpaceTraits

    const Real& t() const
    {
        return t_;
    }

    void set_t(const Real& t)
    {
        if (t < 0.0)
        {
            throw std::invalid_argument("the time must be positive.");
        }
        t_ = t;
    }

    // CompartmentSpaceTraits

    const Real& volume() const
    {
        return volume_;
    }

    Real num_molecules(const Species& sp) const
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            // throw NotFound("Species not found");
            return 0.0;
        }

        return num_molecules_[(*i).second];
    }

    std::vector<Species> list_species() const
    {
        return species_;
    }

    // CompartmentSpace member functions

    void set_volume(const Real& volume)
    {
        if (volume <= 0.0)
        {
            throw std::invalid_argument("The volume must be positive.");
        }

        volume_ = volume;
    }

    void add_molecules(const Species& sp, const Real& num)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            reserve_species(sp);
            i = index_map_.find(sp);
        }

        num_molecules_[(*i).second] += num;
    }

    void remove_molecules(const Species& sp, const Real& num)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            throw NotFound("Species not found");
        }

        num_molecules_[(*i).second] -= num;
    }

    // Optional members

    void set_num_molecules(const Species& sp, const Real& num)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            throw NotFound("Species not found");
        }

        num_molecules_[(*i).second] = num;
    }

    void save(const std::string& filename) const
    {
        boost::scoped_ptr<H5::H5File>
            fout(new H5::H5File(filename, H5F_ACC_TRUNC));
        boost::scoped_ptr<H5::Group>
            group(new H5::Group(fout->createGroup("CompartmentSpace")));
        save_compartment_space<ODEWorld, H5DataTypeTraits_double>(*this, group.get());

        const uint32_t space_type = static_cast<uint32_t>(Space::ELSE);
        group->openAttribute("type").write(H5::PredType::STD_I32LE, &space_type);
    }

    void load(const std::string& filename)
    {
        clear();
        boost::scoped_ptr<H5::H5File>
            fin(new H5::H5File(filename, H5F_ACC_RDONLY));
        const H5::Group group(fin->openGroup("CompartmentSpace"));
        load_compartment_space<ODEWorld, H5DataTypeTraits_double>(group, this);
    }

    bool has_species(const Species& sp)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        return (i != index_map_.end());
    }

    void reserve_species(const Species& sp)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i != index_map_.end())
        {
            throw AlreadyExists("Species already exists");
        }

        index_map_.insert(std::make_pair(sp, num_molecules_.size()));
        species_.push_back(sp);
        num_molecules_.push_back(0);
    }

    void release_species(const Species& sp)
    {
        species_map_type::iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            throw NotFound("Species not found");
        }

        species_map_type::mapped_type
            idx((*i).second), last_idx(num_molecules_.size() - 1);
        if (idx != last_idx)
        {
            const species_container_type::size_type
                idx_(static_cast<species_container_type::size_type>(idx)),
                last_idx_(
                    static_cast<species_container_type::size_type>(last_idx));
            const Species& last_sp(species_[last_idx_]);
            species_[idx_] = last_sp;
            num_molecules_[idx] = num_molecules_[last_idx];
            index_map_[last_sp] = idx;
        }

        species_.pop_back();
        num_molecules_.pop_back();
        index_map_.erase(sp);
    }

    void bind_to(boost::shared_ptr<Model> model)
    {
        if (boost::shared_ptr<Model> bound_model = this->model_.lock())
        {
            if (bound_model.get() != model.get())
            {
                std::cerr << "Warning: Model already bound to ODEWorld."
                    << std::endl;
            }
        }
        this->model_ = model;
    }

protected:

    void clear()
    {
        index_map_.clear();
        num_molecules_.clear();
        species_.clear();
    }

protected:

    Real volume_;
    Real t_;

    num_molecules_container_type num_molecules_;
    species_container_type species_;
    species_map_type index_map_;

    boost::weak_ptr<Model> model_;
};

} // ode

} // ecell4

#endif /* __ECELL4_ODE_ODE_WORLD_HPP */
