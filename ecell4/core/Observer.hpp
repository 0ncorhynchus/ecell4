#ifndef __ECELL4_OBSERVER_HPP
#define __ECELL4_OBSERVER_HPP

#include "types.hpp"
#include "Space.hpp"


namespace ecell4
{

class Observer
{
public:

    Observer(const bool e)
        : every_(e)
    {
        ;
    }

    virtual ~Observer()
    {
        ; // do nothing
    }

    virtual const Real next_time() const
    {
        return inf;
    }

    virtual void initialize(const Space* space)
    {
        ;
    }

    virtual void fire(const Space* space) = 0;

    bool every()
    {
        return every_;
    }

private:

    const bool every_;
};

class FixedIntervalObserver
    : public Observer
{
public:

    typedef Observer base_type;

public:

    FixedIntervalObserver(const Real& dt)
        : base_type(false), tnext_(0.0), dt_(dt), num_steps_(0)
    {
        ;
    }

    virtual ~FixedIntervalObserver()
    {
        ;
    }

    const Real next_time() const
    {
        return tnext_;
    }

    virtual void initialize(const Space* space)
    {
        tnext_ = space->t();
        num_steps_ = 0;
    }

    virtual void fire(const Space* space)
    {
        tnext_ += dt_;
        ++num_steps_;
    }

protected:

    Real tnext_, dt_;
    Integer num_steps_;
};

struct NumberLogger
{
    typedef std::vector<std::vector<Real> > data_container_type;
    typedef std::vector<Species> species_container_type;

    NumberLogger(const std::vector<std::string>& species)
    {
        targets.reserve(species.size());
        for (std::vector<std::string>::const_iterator i(species.begin());
            i != species.end(); ++i)
        {
            targets.push_back(Species(*i));
        }
    }

    ~NumberLogger()
    {
        ;
    }

    void initialize()
    {
        data.clear();
    }

    void log(const Space* space)
    {
        data_container_type::value_type tmp;
        tmp.push_back(space->t());
        for (species_container_type::const_iterator i(targets.begin());
            i != targets.end(); ++i)
        {
            tmp.push_back(space->num_molecules(*i));
        }
        data.push_back(tmp);
    }

    data_container_type data;
    species_container_type targets;
};

class FixedIntervalNumberObserver
    : public FixedIntervalObserver
{
public:

    typedef FixedIntervalObserver base_type;

public:

    FixedIntervalNumberObserver(const Real& dt, const std::vector<std::string>& species)
        : base_type(dt), logger_(species)
    {
        ;
    }

    virtual ~FixedIntervalNumberObserver()
    {
        ;
    }

    virtual void initialize(const Space* space)
    {
        base_type::initialize(space);
        logger_.initialize();
    }

    virtual void fire(const Space* space)
    {
        logger_.log(space);
        base_type::fire(space);
    }

    NumberLogger::data_container_type data() const
    {
        return logger_.data;
    }

protected:

    NumberLogger logger_;
};

class NumberObserver
    : public Observer
{
public:

    typedef Observer base_type;

public:

    NumberObserver(const std::vector<std::string>& species)
        : base_type(true), logger_(species)
    {
        ;
    }

    virtual ~NumberObserver()
    {
        ;
    }

    virtual void initialize(const Space* space)
    {
        logger_.initialize();
    }

    virtual void fire(const Space* space)
    {
        logger_.log(space);
    }

    NumberLogger::data_container_type data() const
    {
        return logger_.data;
    }

protected:

    NumberLogger logger_;
};

} // ecell4

#endif /* __ECELL4_OBSEVER_HPP */
