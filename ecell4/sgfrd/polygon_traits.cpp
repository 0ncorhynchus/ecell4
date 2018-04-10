#include "polygon_traits.hpp"
#include <cmath>
#include <cassert>

namespace ecell4
{

namespace sgfrd
{

polygon::vertex_descripter
make_vertex_information(const ecell4::Polygon<polygon_traits>& poly,
                        const polygon::vertex_id_type vid)
{
    polygon::vertex_descripter vtx;
    const std::vector<polygon::face_id_type>& faces = poly.connecting_faces(vid);
    vtx.neighbor_faces.reserve(faces.size()*2);
    vtx.neighbor_vertices.reserve(faces.size()+1);

    vtx.position[0] = std::numeric_limits<Real>::max();
    vtx.max_conical_shell_size = std::numeric_limits<Real>::max();
    for(std::vector<polygon::face_id_type>::const_iterator
        iter = faces.begin(); iter != faces.end(); ++iter)
    {
        const ecell4::Triangle& tri = poly.triangle_at(*iter);
        const boost::array<polygon::vertex_id_type, 3>& vtxs =
            poly.connecting_vertices(*iter);
        std::size_t vidx;
             if(vtxs[0] == vid) vidx = 0;
        else if(vtxs[1] == vid) vidx = 1;
        else if(vtxs[2] == vid) vidx = 2;
        else assert(false);

        if(vtx.position[0] == std::numeric_limits<Real>::max())
        {
            vtx.position = tri.vertex_at(vidx);
        }
        else
        {
            assert(length_sq(tri.vertex_at(vidx) - vtx.position) < 1e-4);
        }

        const Real length_of_perpendicular = length_sq(tri.edge_at(vidx)) -
                std::pow(dot_product(tri.edge_at(vidx),
                                     tri.edge_at(vidx==2?0:vidx+1)*(-1.0)), 2) /
                length_sq(tri.edge_at(vidx==2?0:vidx+1));
        vtx.max_conical_shell_size = std::min(vtx.max_conical_shell_size,
                                              length_of_perpendicular);
    }
    vtx.max_conical_shell_size = std::sqrt(vtx.max_conical_shell_size);

    for(std::vector<polygon::face_id_type>::const_iterator
        iter = faces.begin(); iter != faces.end(); ++iter)
    {
        const boost::array<polygon::vertex_id_type, 3>& vtxs =
            poly.connecting_vertices(*iter);
        uniquely_add(vtx.neighbor_vertices, vtxs[0]);
        uniquely_add(vtx.neighbor_vertices, vtxs[1]);
        uniquely_add(vtx.neighbor_vertices, vtxs[2]);
    }
    /* remove itself from neighbor_vertices */
    {
        const std::vector<polygon::vertex_id_type>::iterator self =
            std::find(vtx.neighbor_vertices.begin(),
                      vtx.neighbor_vertices.end(), vid);
        vtx.neighbor_vertices.erase(self);
    }


    for(std::vector<polygon::face_id_type>::const_iterator
        iter = faces.begin(); iter != faces.end(); ++iter)
    {
        const boost::array<polygon::face_id_type, 3>& adj
            = poly.adjacent_faces(*iter);
        uniquely_add(vtx.neighbor_faces,  *iter);
        uniquely_add(vtx.neighbor_faces, adj[0]);
        uniquely_add(vtx.neighbor_faces, adj[1]);
        uniquely_add(vtx.neighbor_faces, adj[2]);
    }

    return vtx;
}

polygon::face_descripter
make_face_information(const ecell4::Polygon<polygon_traits>& poly,
                      const polygon::face_id_type fid)
{
    polygon::face_descripter face;
    face.neighbor_faces = poly.neighbor_faces(fid);
    face.neighbor_vertices.reserve(6);

    for(std::vector<polygon::face_id_type>::const_iterator
        iter = face.neighbor_faces.begin(); iter != face.neighbor_faces.end(); ++iter)
    {
        const boost::array<polygon::vertex_id_type, 3>& vtxs =
            poly.connecting_vertices(*iter);
        uniquely_add(face.neighbor_vertices, vtxs[0]);
        uniquely_add(face.neighbor_vertices, vtxs[1]);
        uniquely_add(face.neighbor_vertices, vtxs[2]);
    }

    std::size_t idx=0;
    const Triangle& target = poly.triangle_at(fid);
    const boost::array<polygon::vertex_id_type, 3>& vtxs = // vtx of this face
        poly.connecting_vertices(fid);
    const boost::array<polygon::face_id_type, 3>& adjs = // id of adjacent face
        poly.adjacent_faces(fid);
    for(boost::array<polygon::face_id_type, 3>::const_iterator
        iter = adjs.begin(); iter != adjs.end(); ++iter)
    {
        const boost::array<polygon::vertex_id_type, 3>& vtxs_of_adj =
            poly.connecting_vertices(*iter); // vtx ids of adjacent face

        std::size_t vtx_idx = std::numeric_limits<std::size_t>::max();
        for(std::size_t i=0; i<3; ++i)
        {
            if(std::find(vtxs.begin(), vtxs.end(), vtxs_of_adj[i]) == vtxs.end())
            {
                assert(vtx_idx == std::numeric_limits<std::size_t>::max());
                vtx_idx = i;
            }
        }

        const Triangle& tri = poly.triangle_at(*iter);
        const Real3 vpos   = tri.vertex_at(vtx_idx);
        const Real  tilt   = ecell4::calc_angle(target.normal(), tri.normal());
        const Real3 axis   = tri.edge_at(vtx_idx==2?0:vtx_idx+1);
        const Real3 origin = tri.vertex_at(vtx_idx==2?0:vtx_idx+1);
        const Real3 developped_vtx = origin + ecell4::rotate(tilt, axis, vpos - origin);

        face.barrier.at(idx) = ecell4::Segment(
                developped_vtx, tri.vertex_at(vtx_idx==2?0:vtx_idx+1));
        idx++;
        face.barrier.at(idx) = ecell4::Segment(
                developped_vtx, tri.vertex_at(vtx_idx==0?2:vtx_idx-1));
        idx++;
    }

    return face;
}

void setup_descriptors(ecell4::Polygon<polygon_traits>& poly)
{
    const std::vector<polygon::face_id_type> fids = poly.list_face_id();
    for(std::vector<polygon::face_id_type>::const_iterator
        iter = fids.begin(); iter != fids.end(); ++iter)
    {
        poly.at(*iter) = make_face_information(poly, *iter);
    }

    const std::vector<polygon::vertex_id_type> vids = poly.list_vertex_id();
    for(std::vector<polygon::vertex_id_type>::const_iterator
        iter = vids.begin(); iter != vids.end(); ++iter)
    {
        poly.at(*iter) = make_vertex_information(poly, *iter);
    }

    return ;
}

} // sgfrd
} // ecell4
