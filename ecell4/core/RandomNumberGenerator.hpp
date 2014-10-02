#ifndef __ECELL4_RANDOM_NUMBER_GENERATOR_HPP
#define __ECELL4_RANDOM_NUMBER_GENERATOR_HPP

#include <ctime>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <hdf5.h>
#include <H5Cpp.h>

#include "types.hpp"
#include "Position3.hpp"


namespace ecell4
{

class RandomNumberGenerator
{
public:

    virtual ~RandomNumberGenerator()
    {
        ;
    }

    virtual Real uniform(Real min, Real max) = 0;
    virtual Integer uniform_int(Integer min, Integer max) = 0;
    virtual Real gaussian(Real mean, Real sigma) = 0;
    virtual Integer binomial(Real p, Integer n) = 0;
    virtual Position3 direction3d(Real length) = 0;

    virtual Real normal(Real mean, Real sigma)  = 0;

    virtual void dir_2d(Real *x, Real *y) = 0;
    virtual void dir_3d(Real *x, Real *y, Real *z) = 0;
    virtual double operator() () = 0;
    virtual void seed(Integer val) = 0;
    virtual void seed() = 0;

    virtual void save(H5::CommonFG* root) const = 0;
    virtual void load(const H5::CommonFG& root) = 0;
};

template<typename Telem_>
inline void shuffle(RandomNumberGenerator& rng, std::vector<Telem_>& cont)
{
    typedef std::vector<Telem_> container_type;
    for (typename container_type::size_type i(cont.size()); i > 0;)
    {
        --i;
        typename container_type::size_type const j(rng.uniform_int(0, i));
        std::swap(cont[i], cont[j]);
    }
}

class GSLRandomNumberGenerator
    : public RandomNumberGenerator
{
public:

    typedef boost::shared_ptr<gsl_rng> rng_handle;

public:

    Real uniform(Real min, Real max)
    {
        return gsl_rng_uniform(rng_.get()) * (max - min) + min;
    }

    Integer uniform_int(Integer min, Integer max)
    {
        return gsl_rng_uniform_int(rng_.get(), max - min + 1) + min;
    }

    Real normal(Real loc, Real scale)
    {   // This function is implecated for comatible for epdp::GSLRandomNumberGenerator.
        // This function is the same as uniform().
        return gsl_ran_gaussian(rng_.get(), scale) + loc;

    }

    Real gaussian(Real mean, Real sigma)
    {
        return gsl_ran_gaussian(rng_.get(), sigma) + mean;
    }

    Integer binomial(Real p, Integer n)
    {
        return gsl_ran_binomial(rng_.get(), p, n);
    }

    Position3 direction3d(Real length)
    {
        double x, y, z;
        gsl_ran_dir_3d(rng_.get(), &x, &y, &z);
        return Position3(x * length, y * length, z * length);
    }

    void seed(Integer val)
    {
        gsl_rng_set(rng_.get(), val);
    }

    void seed()
    {
        gsl_rng_set(rng_.get(), unsigned(std::time(0)));
    }

    void save(H5::CommonFG* root) const;
    void load(const H5::CommonFG& root);

    GSLRandomNumberGenerator(rng_handle hdl)
        : rng_(hdl)
    {
        ;
    }

    GSLRandomNumberGenerator(gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937))
        : rng_(rng, &gsl_rng_free)
    {
        ;
    }

    /** for epdp
     */

    Integer get_raw()
    {
        return gsl_rng_get(rng_.get());
    }

    void dir_2d(Real *x, Real *y)
    {
        gsl_ran_dir_2d(rng_.get(), x, y);
    }

    void dir_3d(Real *x, Real *y, Real *z)
    {
        gsl_ran_dir_3d(rng_.get(), x, y, z);
    }

    Real operator()()
    {
        return gsl_rng_uniform(rng_.get());
    }

protected:

    rng_handle rng_;
};

} // ecell4

#endif /* __ECELL4_RANDOM_NUMBER_GENERATOR_HPP */
