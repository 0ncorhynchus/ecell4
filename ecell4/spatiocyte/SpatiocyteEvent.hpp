#ifndef ECELL4_SPATIOCYTE_EVENT_HPP
#define ECELL4_SPATIOCYTE_EVENT_HPP

#include "SpatiocyteReactions.hpp"
#include "SpatiocyteWorld.hpp"
#include <ecell4/core/EventScheduler.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/ReactionRule.hpp>

namespace ecell4
{

namespace spatiocyte
{

struct SpatiocyteEvent : public Event
{
public:
    typedef std::pair<ReactionRule, ReactionInfo> reaction_type;

    SpatiocyteEvent(Real const &time) : Event(time) {}
    virtual ~SpatiocyteEvent() {}

    const std::vector<reaction_type> &reactions() const { return reactions_; }

    virtual void fire()
    {
        reactions_.clear();
        fire_();
    }

protected:
    virtual void fire_() = 0;

    void push_reaction(const reaction_type &reaction)
    {
        reactions_.push_back(reaction);
    }

    std::vector<reaction_type> reactions_;
};

struct StepEvent : SpatiocyteEvent
{
    StepEvent(boost::shared_ptr<Model> model,
              boost::shared_ptr<SpatiocyteWorld> world, const Species &species,
              const Real &t, const Real alpha = 1.0);
    virtual ~StepEvent() {}

    Species const &species() const { return mpool_->species(); }

    Real const &alpha() const { return alpha_; }

    void fire_()
    {
        walk(alpha_);
        time_ += dt_;
    }

    void walk(const Real &alpha)
    {
        if (alpha < 0 || alpha > 1)
        {
            return; // INVALID ALPHA VALUE
        }

        MoleculePool::container_type voxels;
        copy(mpool_->begin(), mpool_->end(), back_inserter(voxels));

        std::size_t idx(0);
        for (const auto &info : voxels)
        {
            const Voxel voxel(space_, info.coordinate);

            if (voxel.get_voxel_pool() != mpool_)
            {
                // should skip if a voxel is not the target species.
                // when reaction has occured before, a voxel can be changed.
                continue;
            }

            const Voxel neighbor(
                world_->get_neighbor_randomly(voxel, dimension()));

            if (world_->can_move(voxel, neighbor))
            {
                if (world_->rng()->uniform(0, 1) <= alpha)
                    world_->move(voxel, neighbor, /*candidate=*/idx);
            }
            else
            {
                attempt_reaction_(info, neighbor, alpha);
            }

            ++idx;
        }
    }

    virtual const Shape::dimension_kind dimension() const = 0;

protected:
    void attempt_reaction_(const SpatiocyteWorld::coordinate_id_pair_type &info,
                           const Voxel &dst, const Real &alpha);

protected:
    boost::shared_ptr<Model> model_;
    boost::shared_ptr<SpatiocyteWorld> world_;
    boost::weak_ptr<VoxelSpaceBase> space_;
    boost::shared_ptr<MoleculePool> mpool_;

    const Real alpha_;
};

struct StepEvent3D : StepEvent
{
    StepEvent3D(boost::shared_ptr<Model> model,
                boost::shared_ptr<SpatiocyteWorld> world,
                const Species &species, const Real &t, const Real alpha = 1.0);

    const Shape::dimension_kind dimension() const { return Shape::THREE; }
};

struct StepEvent2D : StepEvent
{
    StepEvent2D(boost::shared_ptr<Model> model,
                boost::shared_ptr<SpatiocyteWorld> world,
                const Species &species, const Real &t, const Real alpha = 1.0);

    const Shape::dimension_kind dimension() const { return Shape::TWO; }
};

struct ZerothOrderReactionEvent : SpatiocyteEvent
{
    ZerothOrderReactionEvent(boost::shared_ptr<SpatiocyteWorld> world,
                             const ReactionRule &rule, const Real &t);

    virtual ~ZerothOrderReactionEvent() {}
    virtual void fire_();

    Real draw_dt();
    virtual void interrupt(Real const &t) { time_ = t + draw_dt(); }

protected:
    boost::shared_ptr<SpatiocyteWorld> world_;
    ReactionRule rule_;
};

struct FirstOrderReactionEvent : SpatiocyteEvent
{
    FirstOrderReactionEvent(boost::shared_ptr<SpatiocyteWorld> world,
                            const ReactionRule &rule, const Real &t);

    virtual ~FirstOrderReactionEvent() {}
    virtual void fire_();

    Real draw_dt();
    virtual void interrupt(Real const &t) { time_ = t + draw_dt(); }

protected:
    ReactionInfo::Item choice()
    {
        const Species &species(rule_.reactants().at(0));
        if (const auto space_and_molecule_pool =
                world_->find_space_and_molecule_pool(species))
        {
            const auto space = space_and_molecule_pool->first;
            const auto molecule_pool = space_and_molecule_pool->second;

            const auto i =
                rng_.lock()->uniform_int(0, molecule_pool->size() - 1);
            const auto &info = molecule_pool->at(i);

            return ReactionInfo::Item(info.pid, species,
                                      Voxel(space, info.coordinate));
        }
        throw "MoleculePool is not found";
    }

    boost::shared_ptr<SpatiocyteWorld> world_;
    boost::weak_ptr<RandomNumberGenerator> rng_;
    ReactionRule rule_;
};

} // namespace spatiocyte

} // namespace ecell4

#endif /* ECELL4_SPATIOCYTE_EVENT_HPP */
