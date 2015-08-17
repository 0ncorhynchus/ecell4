#ifndef __ECELL4_SHAPE_OPERATORS
#define __ECELL4_SHAPE_OPERATORS

#include "Shape.hpp"

namespace ecell4
{

struct Union
    : public Shape
{
public:

    Union(const boost::shared_ptr<const Shape>& a,
          const boost::shared_ptr<const Shape>& b)
        : a_(a), b_(b)
    {
        ;
    }

    Union(const Union& other)
        : a_(other.a_), b_(other.b_)
    {
        ;
    }

    ~Union()
    {
        ; // do nothing
    }

    virtual dimension_kind dimension() const
    {
        return a_->dimension();
    }

    virtual Real is_inside(const Real3& coord) const
    {
        const Real retval1 = a_->is_inside(coord);
        const Real retval2 = b_->is_inside(coord);
        return std::min(retval1, retval2);
    }

    virtual Real3 draw_position(
        boost::shared_ptr<RandomNumberGenerator>& rng) const
    {
        throw NotImplemented("not implemented yet");
    }

    virtual bool test_AABB(const Real3& l, const Real3& u) const
    {
        return (a_->test_AABB(l, u) || b_->test_AABB(l, u));
    }

    virtual void bounding_box(
        const Real3& edge_lengths, Real3& lower, Real3& upper) const
    {
        a_->bounding_box(edge_lengths, lower, upper);

        Real3 l, u;
        b_->bounding_box(edge_lengths, l, u);
        for (unsigned int dim(0); dim < 3; ++dim)
        {
            lower[dim] = std::min(lower[dim], l[dim]);
            upper[dim] = std::max(upper[dim], u[dim]);
        }
    }

protected:

    const boost::shared_ptr<const Shape> a_;
    const boost::shared_ptr<const Shape> b_;
};

struct Complement
    : public Shape
{
public:

    Complement(const boost::shared_ptr<const Shape>& a,
               const boost::shared_ptr<const Shape>& b)
        : a_(a), b_(b)
    {
        ;
    }

    Complement(const Complement& other)
        : a_(other.a_), b_(other.b_)
    {
        ;
    }

    ~Complement()
    {
        ; // do nothing
    }

    virtual dimension_kind dimension() const
    {
        return a_->dimension();
    }

    virtual Real is_inside(const Real3& coord) const
    {
        if (b_->is_inside(coord) > 0)
        {
            return a_->is_inside(coord);
        }
        else
        {
            return inf;
        }
    }

    virtual Real3 draw_position(
        boost::shared_ptr<RandomNumberGenerator>& rng) const
    {
        throw NotImplemented("not implemented yet");
    }

    virtual bool test_AABB(const Real3& l, const Real3& u) const
    {
        return (a_->test_AABB(l, u) && !b_->test_AABB(l, u));
    }

    virtual void bounding_box(
        const Real3& edge_lengths, Real3& lower, Real3& upper) const
    {
        return a_->bounding_box(edge_lengths, lower, upper);
    }

protected:

    const boost::shared_ptr<const Shape> a_;
    const boost::shared_ptr<const Shape> b_;
};

}

#endif /* __ECELL4_SHAPE_OPERATORS */
