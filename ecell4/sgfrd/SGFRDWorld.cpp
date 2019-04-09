#include <ecell4/sgfrd/SGFRDWorld.hpp>

namespace ecell4
{
namespace sgfrd
{

std::pair<std::pair<ParticleID, Particle>, bool>
SGFRDWorld::new_particle(const Particle& p)
{
    const std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        overlap3 = list_particles_within_radius(p.position(), p.radius());
    if(!overlap3.empty())
    {
        return std::make_pair(std::make_pair(pidgen_(), p), false);
    }
    const ParticleID pid = pidgen_();
    return std::make_pair(std::make_pair(pid, p), update_particle(pid, p));
}

std::pair<std::pair<ParticleID, Particle>, bool>
SGFRDWorld::new_particle(const Particle& p, const FaceID& fid)
{
    const ParticleID pid = pidgen_();
    // now this consider only 2D particles
    const std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        overlap2d(list_particles_within_radius(
            std::make_pair(p.position(), fid), p.radius()));

    if(!overlap2d.empty())
    {
//         std::cout << "overlapped particle = {";
//         std::cout << '{' << overlap2d.front().first.first
//                   << ", distance = " << overlap2d.front().second << '}';
//         for(std::size_t i=1; i<overlap2d.size(); ++i)
//         {
//             std::cout << ", {" << overlap2d.at(i).first.first
//                       << ", distance = " << overlap2d.at(i).second << '}';
//         }
//         std::cout << '}' << std::endl;
        return std::make_pair(std::make_pair(pid, p), false);
    }
    return std::make_pair(std::make_pair(pid, p), update_particle(pid, p, fid));
}

std::pair<std::pair<ParticleID, Particle>, bool>
SGFRDWorld::throw_in_particle(const Species& sp)
{
    const Real r = sp.get_attribute_as<Real>("radius");
    const Real D = sp.get_attribute_as<Real>("D");

    Real3 pos; FaceID fid;
    pos = this->polygon_->draw_position(this->rng_, fid);
    const Particle p(sp, pos, r, D);
    return this->new_particle(p, fid);
}

void SGFRDWorld::add_molecule(const Species& sp, const std::size_t N)
{
    for(std::size_t i=0; i<N; ++i)
    {
        while(this->throw_in_particle(sp).second == false)
        {
            /*do nothing*/
        }
    }
    return;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
SGFRDWorld::list_particles_within_radius(
        const std::pair<Real3, FaceID>& pos, const Real& radius) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>> retval;

    // look particles on the same face
    for(const auto& pid : this->list_particleIDs(pos.second))
    {
        const std::pair<ParticleID, Particle> pp = ps_->get_particle(pid);
        const Real dist = length(pos.first - pp.second.position()) -
                          pp.second.radius();
        if(dist < radius)
        {
            retval.push_back(std::make_pair(pp, dist));
        }
    }

    // look particles around
    for(const auto& fid : polygon_->neighbor_faces_of(pos.second))
    {
        for(const auto& pid : this->list_particleIDs(fid))
        {
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(pid);
            const Real dist = ecell4::polygon::distance(*polygon_,
                pos, std::make_pair(pp.second.position(), get_face_id(pp.first))
                ) - pp.second.radius();
            if(dist < radius)
            {
                retval.push_back(std::make_pair(pp, dist));
            }
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ParticleID, Particle>, Real>());
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
SGFRDWorld::list_particles_within_radius(
        const std::pair<Real3, FaceID>& pos, const Real& radius,
        const ParticleID& ignore) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>> retval;

    // look particles on the same face
    for(const auto& pid : this->list_particleIDs(pos.second))
    {
        if(pid == ignore)
        {
            continue;
        }
        const std::pair<ParticleID, Particle> pp = ps_->get_particle(pid);
        const Real dist = length(pos.first - pp.second.position()) -
                          pp.second.radius();

        if(dist < radius)
        {
            retval.push_back(std::make_pair(pp, dist));
        }
    }

    for(const auto& fid : polygon_->neighbor_faces_of(pos.second))
    {
        for(const auto& pid : this->list_particleIDs(fid))
        {
            if(pid == ignore)
            {
                continue;
            }
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(pid);
            const Real dist = ecell4::polygon::distance(*polygon_,
                pos, std::make_pair(pp.second.position(), get_face_id(pp.first))
                ) - pp.second.radius();
            if(dist < radius)
            {
                retval.push_back(std::make_pair(pp, dist));
            }
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ParticleID, Particle>, Real>());
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
SGFRDWorld::list_particles_within_radius(
        const std::pair<Real3, FaceID>& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    {// same face
        for(const auto& pid : this->list_particleIDs(pos.second))
        {
            if(pid == ignore1 || pid == ignore2)
            {
                continue;
            }
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(pid);
            const Real dist = length(pos.first - pp.second.position()) -
                              pp.second.radius();
            if(dist < radius)
            {
                retval.push_back(std::make_pair(pp, dist));
            }
        }
    }

    for(const auto& fid : polygon_->neighbor_faces_of(pos.second))
    {
        for(const auto& pid : this->list_particleIDs(fid))
        {
            if(pid == ignore1 || pid == ignore2)
            {
                continue;
            }
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(pid);
            const Real dist = ecell4::polygon::distance(*polygon_,
                pos, std::make_pair(pp.second.position(), get_face_id(pp.first))
                ) - pp.second.radius();
            if(dist < radius)
            {
                retval.push_back(std::make_pair(pp, dist));
            }
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ParticleID, Particle>, Real>());
    return retval;
}


bool SGFRDWorld::check_no_overlap(
        const std::pair<Real3, FaceID>& pos, const Real& radius) const
{
    {// same face
        const std::vector<ParticleID>& ids = this->list_particleIDs(pos.second);
        for(std::vector<ParticleID>::const_iterator
            i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = length(pos.first - pp.second.position()) -
                              pp.second.radius();
            if(dist < radius) {return false;}
        }
    }

    std::vector<FaceID> const& neighbors =
        polygon_->neighbor_faces_of(pos.second);
    for(std::vector<FaceID>::const_iterator
        iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
        const std::vector<ParticleID>& ids = registrator_.elements_over(*iter);
        for(std::vector<ParticleID>::const_iterator
                i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = ecell4::polygon::distance(*polygon_, pos,
                std::make_pair(pp.second.position(), get_face_id(pp.first))) -
                pp.second.radius();
            if(dist < radius) {return false;}
        }
    }
    return true;
}

bool SGFRDWorld::check_no_overlap(
        const std::pair<Real3, FaceID>& pos, const Real& radius,
        const ParticleID& ignore) const
{
    {// same face
        const std::vector<ParticleID>& ids = this->list_particleIDs(pos.second);
        for(std::vector<ParticleID>::const_iterator
            i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            if(*i == ignore) {continue;}
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = length(pos.first - pp.second.position()) -
                              pp.second.radius();
            if(dist < radius) {return false;}
        }
    }

    std::vector<FaceID> const& neighbors =
        polygon_->neighbor_faces_of(pos.second);
    for(std::vector<FaceID>::const_iterator
        iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
        const std::vector<ParticleID>& ids = this->list_particleIDs(*iter);
        for(std::vector<ParticleID>::const_iterator
                i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            if(*i == ignore) {continue;}
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = ecell4::polygon::distance(*polygon_, pos,
                std::make_pair(pp.second.position(), get_face_id(pp.first))) -
                pp.second.radius();
            if(dist < radius) {return false;}
        }
    }
    return true;
}

bool SGFRDWorld::check_no_overlap(
        const std::pair<Real3, FaceID>& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
{
    {// same face
        const std::vector<ParticleID>& ids = this->list_particleIDs(pos.second);
        for(std::vector<ParticleID>::const_iterator
            i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            if(*i == ignore1 || *i == ignore2) {continue;}
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = length(pos.first - pp.second.position()) -
                              pp.second.radius();
            if(dist < radius) {return false;}
        }
    }

    std::vector<FaceID> const& neighbors =
        polygon_->neighbor_faces_of(pos.second);
    for(std::vector<FaceID>::const_iterator
        iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
        const std::vector<ParticleID>& ids = this->list_particleIDs(*iter);
        for(std::vector<ParticleID>::const_iterator
                i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            if(*i == ignore1 || *i == ignore2) {continue;}
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = ecell4::polygon::distance(*polygon_, pos,
                std::make_pair(pp.second.position(), get_face_id(pp.first))) -
                pp.second.radius();
            if(dist < radius) {return false;}
        }
    }
    return true;
}

}// sgfrd
}// ecell4
