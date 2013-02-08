#ifndef __BD_SIMULATOR_HPP
#define __BD_SIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/Model.hpp>
#include <ecell4/core/Simulator.hpp>

#include "BDWorld.hpp"
#include "BDPropagator.hpp"


namespace ecell4
{

namespace bd
{

class BDSimulator
    : public Simulator
{
public:

    BDSimulator(boost::shared_ptr<Model> model, boost::shared_ptr<BDWorld> world)
        : model_(model), world_(world), num_steps_(0), dt_(0)
    {
        ;
    }

    // SimulatorTraits

    Real t() const
    {
        return (*world_).t();
    }

    Real dt() const
    {
        return dt_;
    }

    Integer num_steps() const
    {
        return num_steps_;
    }

    void step();
    bool step(Real const& upto);

    // Optional members

    void set_t(Real const& t)
    {
        (*world_).set_t(t);
    }

    void set_dt(Real const& dt)
    {
        if (dt <= 0)
        {
            throw std::invalid_argument("The step size must be positive.");
        }
        dt_ = dt;
    }

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return (*world_).rng();
    }

protected:

    boost::shared_ptr<Model> model_;
    boost::shared_ptr<BDWorld> world_;

    /**
     * the protected internal state of BDSimulator.
     * they are needed to be saved/loaded with Visitor pattern.
     */
    Real dt_;
    Integer num_steps_;
};

} // bd

} // ecell4

#endif /* __BD_SIMULATOR_HPP */
