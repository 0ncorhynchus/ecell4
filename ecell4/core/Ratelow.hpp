#ifndef __ECELL4_RATELOW_HPP
#define __ECELL4_RATELOW_HPP

#include "types.hpp"
#include "exceptions.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/variant.hpp>
#include <vector>
#include <algorithm>

// For Jacobi
#include <boost/numeric/ublas/matrix.hpp>

namespace ecell4
{

class Ratelow
{
public:
    // The order of the species must be the same as
    //  reactants' container of ReactionRule object.
    //
    //  state_container_type must be resized when called 
    //   jacobi_func and deriv_func.
    typedef std::vector<Real> state_container_type;
    typedef boost::numeric::ublas::matrix<double> matrix_type;
public:
    virtual bool is_available() const = 0;
    virtual Real operator()(state_container_type const &state_array, Real volume) = 0;
    virtual Real deriv_func(state_container_type const &state_array, Real volume) = 0;
    virtual void jacobi_func(matrix_type &jacobian, state_container_type const &state_array, Real const volume) = 0;
};

class RatelowCppCallback : public Ratelow
{
    /** Function object for calculate ratelow called by C++
     *  This class must not be exposed to Cython interface.
     */
public:
    typedef double (*Ratelow_Callback)(state_container_type const &, double);
public:
    RatelowCppCallback(Ratelow_Callback func) : func_(func), h_(1.0e-8) {;}
    RatelowCppCallback() : func_(0), h_(1.0e-8) {}
    virtual ~RatelowCppCallback(){;}

    virtual bool is_available() const
    {
        return this->func_ != 0;
    }
    virtual Real operator()(state_container_type const &state_array, Real volume)
    {
        return this->deriv_func(state_array, volume);
    }
    virtual Real deriv_func(state_container_type const &state_array, Real volume)
    {
        if (!is_available())
        {
            throw IllegalState("Callback Function has not been registerd");
        }
        return this->func_(state_array, volume);
    }
    virtual void jacobi_func(
            matrix_type &jacobian, state_container_type const &state_array, Real volume)
    {
        Real h(this->h_); //XXX  1.0e-8. Should be fixed
        std::fill(jacobian.data().begin(), jacobian.data().end(), Real(0.0));
        Real flux( this->deriv_func(state_array, volume) );
        double num_reactants(state_array.size());
        for(int i(0); i < num_reactants; i++) 
        {
            //XXX For now, we are using FORWARD difference method.
            state_container_type h_shift(state_array);
            h_shift[i] += h;
            double deriv = ((this->deriv_func(h_shift, volume)) - flux) / h;
            for(int j(0); j < jacobian.size1() ; j++)
            {
                if (j < num_reactants) 
                {
                    jacobian(j, i) -= deriv;
                }
                else
                {
                    jacobian(j, i) += deriv;
                }
            }
        }
    }

    Ratelow_Callback get_callback() const
    {
        return this->func_;
    }
    Ratelow_Callback set_callback(Ratelow_Callback new_func)
    {
        if (new_func == 0)
        {
            throw std::invalid_argument("Ratelow Callback must not be 0");
        }
        Ratelow_Callback prev = get_callback();
        this->func_ = new_func;
        return prev;
    }
private:
    Real h_;
    Ratelow_Callback func_;
};

class RatelowMassAction : public Ratelow
{
public:
    RatelowMassAction(Real k = 0.0, std::size_t num_reactant = 0) 
        : k_(k), num_reactant_(num_reactant) {}
    virtual ~RatelowMassAction(){;}
    virtual bool is_available() const
    {
        return true;    // always true
    }
    virtual Real operator()(state_container_type const &state_array, Real volume)
    {
        // Forward to deriv_func()
        Real flux( this->deriv_func(state_array, volume) );
        return flux;
    }
    virtual Real deriv_func(state_container_type const &state_array, Real volume)
    {
        // The 1st argument 'state_array' must be resized when calling.
        Real flux(this->k_ * volume);
        for(state_container_type::const_iterator it(state_array.begin());
                it != state_array.end(); it++) 
        {
            flux *= (*it) / volume;
        }
        return flux;
    }
    virtual void jacobi_func(
            matrix_type &jacobian, 
            state_container_type const &state_array, Real const volume)
    {
        // The size of the argument 'state_array' must be resized to the number of reactants.
        // The size of the argument 'jacobian' must be resized to the number of (reactants + products)
        std::fill(jacobian.data().begin(), jacobian.data().end(), Real(0.0));
        Real flux( this->deriv_func(state_array, volume) );
        if (flux == Real(0.0))
        {
            return;
        }
        int num_reactants(state_array.size());
        for(int i(0); i < num_reactants; i++) 
        {
            Real partial(flux / state_array[i]);
            for(int j(0); j < jacobian.size1(); j++)
            {
                if (j < num_reactants)
                {
                    jacobian(j, i) -= partial;
                }
                else
                {
                    jacobian(j, i) += partial;
                }
            }
        }
    }
    void set_k(Real k)
    {
        this->k_ = k;
    }
    Real get_k(Real k) const
    {
        return this->k_;
    }
private:
    Real k_;
    std::size_t num_reactant_;
};

} // ecell4

#endif  //__ECELL4_RATELOW_HPP
