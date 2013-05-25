#include "ODESimulator.hpp"

#include <boost/numeric/odeint.hpp>
namespace odeint = boost::numeric::odeint;

namespace ecell4
{

namespace ode
{

bool ODESimulator::step(const Real& upto)
{
    if (upto <= t())
    {
        return false;
    }

    initialize();

    const std::vector<Species> species(model_->list_species());
    ODESystem::state_type x(species.size());

    {
        ODESystem::state_type::size_type i(0);
        for (NetworkModel::species_container_type::const_iterator
                 it(species.begin()); it != species.end(); ++it)
        {
            x[i] = static_cast<double>(world_->num_molecules(*it));
            ++i;
        }
    }

    typedef odeint::runge_kutta_cash_karp54<ODESystem::state_type>
        error_stepper_type;
    typedef odeint::controlled_runge_kutta<error_stepper_type>
        controlled_stepper_type;

    ODESystem func_obj(model_, world_->volume());
    StateAndTimeBackInserter::state_container_type x_vec;
    StateAndTimeBackInserter::time_container_type times;

    const Real dt(upto - t());

    // size_t steps(odeint::integrate(
    //                  func_obj, x, t(), upto, dt,
    //                  StateAndTimeBackInserter(x_vec, times)));

    const double abs_err(1e-10), rel_err(1e-6), a_x(1.0), a_dxdt(1.0);
    controlled_stepper_type controlled_stepper(
        odeint::default_error_checker<
            double, odeint::range_algebra, odeint::default_operations>(
                abs_err, rel_err, a_x, a_dxdt));
    const size_t steps(odeint::integrate_adaptive(
                           controlled_stepper, func_obj, x, t(), upto, dt,
                           StateAndTimeBackInserter(x_vec, times)));

    {
        ODESystem::state_type::size_type i(0);
        for (NetworkModel::species_container_type::const_iterator
                 it(species.begin()); it != species.end(); ++it)
        {
            world_->set_num_molecules(*it, static_cast<Real>(x_vec[steps][i]));
            ++i;
        }
    }

    set_t(upto);
    // set_t(times[steps]);
    ++num_steps_;
    return false;
}

} // ode

} // ecell4
