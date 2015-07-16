#ifndef __ECELL4_MESH_HPP
#define __ECELL4_MESH_HPP

#include "Shape.hpp"

#include <vtkSmartPointer.h>
#include <vtkCellLocator.h>
#include <vtkSTLReader.h>
#include <vtkOBBTree.h>
#include <vtkPolyData.h>

namespace ecell4
{

struct MeshSurface
    : public Shape
{
public:

    MeshSurface(const std::string filename, const Real3& edge_lengths);
    MeshSurface(const MeshSurface& rhs);

    std::string filename() const
    {
        return filename_;
    }

    Real3 edge_lengths() const
    {
        return edge_lengths_;
    }

    virtual dimension_kind dimension() const
    {
        return TWO;
    }

    virtual Real is_inside(const Real3& pos) const;
    virtual Real3 draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const;
    virtual bool test_AABB(const Real3& l, const Real3& u) const;

    // virtual void bounding_box(
    //     const Real3& edge_lengths, Real3& lower, Real3& upper) const
    // {
    //     double bounds[6];
    //     reader_->GetOutput()->GetBounds(bounds);

    //     lower = Real3(std::max(0.0, bounds[0]), std::max(0.0, bounds[2]), std::max(0.0, bounds[4]));
    //     upper = Real3(std::min(edge_lengths[0], bounds[1]), std::min(edge_lengths[1], bounds[3]), std::min(edge_lengths[2], bounds[5]));
    // }

protected:

    std::string filename_;
    Real3 edge_lengths_;

    Real ratio_;
    Real3 shift_;

    vtkSmartPointer<vtkSTLReader> reader_;
    vtkSmartPointer<vtkOBBTree> tree_;
};

} // ecell4

#endif /* __ECELL4_MESH_HPP */
