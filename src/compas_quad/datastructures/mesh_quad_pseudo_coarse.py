from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from compas.geometry import Polyline
from compas.geometry import discrete_coons_patch

# from compas.datastructures import meshes_join_and_weld

from compas.datastructures import Mesh
from compas_quad.datastructures import PseudoQuadMesh, CoarseQuadMesh

from compas.itertools import pairwise, linspace
from compas.tolerance import Tolerance
geometric_key = Tolerance().geometric_key


__all__ = ['CoarsePseudoQuadMesh']


class CoarsePseudoQuadMesh(PseudoQuadMesh, CoarseQuadMesh):

    def __init__(self):
        super(CoarsePseudoQuadMesh, self).__init__()

    def densification_old(self, edges_to_curves=None):
        """Generate a denser quad mesh from the coarse quad mesh and its strip densities.

        Returns
        -------
        QuadMesh
            A denser quad mesh.
        edges_to_curves : dict, optional
            A dictionary with edges (u, v) pointing to curve for densification. The curves are lists of XYZ points.

        """

        edge_strip = {}
        for strip, edges in self.strips(data=True):
            for u, v in edges:
                edge_strip[u, v] = strip
                edge_strip[v, u] = strip

        pole_map = [geometric_key(self.vertex_coordinates(pole))
                    for pole in self.poles()]

        meshes = []
        for fkey in self.faces():
            polylines = []
            for u, v in self.face_halfedges(fkey):
                d = self.strip_density(edge_strip[u, v])

                if edges_to_curves:
                    polyline = []
                    if (u, v) in edges_to_curves:
                        curve = Polyline(edges_to_curves[u, v])
                        polyline = [curve.point_at(t) for t in linspace(0, 1, d)]
                    else:
                        curve = Polyline(edges_to_curves[v, u])
                        polyline = [curve.point_at(t) for t in linspace(0, 1, d)]
                        polyline[:] = polyline[::-1]
                else:
                    polyline = []
                    curve = Polyline(
                        [self.vertex_coordinates(u), self.vertex_coordinates(v)])
                    for i in range(0, d + 1):
                        point = curve.point_at(float(i) / float(d))
                        polyline.append(point)

                polylines.append(polyline)

            if self.is_face_pseudo_quad(fkey):
                pole = self.attributes['face_pole'][fkey]
                idx = self.face_vertices(fkey).index(pole)
                polylines.insert(idx, None)

            ab, bc, cd, da = polylines

            if cd:
                dc = cd[::-1]
            else:
                dc = None

            if da:
                ad = da[::-1]
            else:
                ad = None

            vertices, faces = discrete_coons_patch(ab, bc, dc, ad)
            faces = [[u for u, v in pairwise(
                face + face[:1]) if u != v] for face in faces]
            mesh = PseudoQuadMesh.from_vertices_and_faces_with_face_poles(
                vertices, faces)
            meshes.append(mesh)

        face_pole_map = {}
        for mesh in meshes:
            for fkey in mesh.faces():
                for u, v in pairwise(mesh.face_vertices(fkey) + mesh.face_vertices(fkey)[: 1]):
                    if geometric_key(mesh.vertex_coordinates(u)) in pole_map and geometric_key(mesh.vertex_coordinates(u)) == geometric_key(mesh.vertex_coordinates(v)):
                        face_pole_map[geometric_key(mesh.face_center(fkey))] = geometric_key(
                            mesh.vertex_coordinates(u))
                        break
        
        dense_mesh = PseudoQuadMesh()
        for mesh in meshes:
            dense_mesh.join(mesh)
        dense_mesh.weld()
        self.dense_mesh(dense_mesh)
        # self.dense_mesh(meshes_join_and_weld(meshes))

        face_pole = {}
        for fkey in self.dense_mesh().faces():
            if len(self.dense_mesh().face_vertices(fkey)) < 3:
                print(self.dense_mesh().face_vertices(fkey))
            if geometric_key(self.dense_mesh().face_center(fkey)) in face_pole_map:
                for vkey in self.dense_mesh().face_vertices(fkey):
                    if geometric_key(self.dense_mesh().vertex_coordinates(vkey)) == face_pole_map[geometric_key(self.dense_mesh().face_center(fkey))]:
                        face_pole[fkey] = vkey
                        break

        self.dense_mesh().attributes['face_pole'] = face_pole
        return self.dense_mesh()

    def densification(self, edges_to_curves=None):
        """Generate a denser quad mesh from the coarse quad mesh and its strip densities.

        Returns
        -------
        QuadMesh
            A denser quad mesh.
        edges_to_curves : dict, optional
            A dictionary with edges (u, v) pointing to curve for densification. The curves are lists of XYZ points.

        """

        edge_strip = {}
        for strip, edges in self.strips(data=True):
            for u, v in edges:
                edge_strip[u, v] = strip
                edge_strip[v, u] = strip

        parent2child = {}
        face_meshes = {}
        for fkey in self.faces():
            polylines = []
            for u, v in self.face_halfedges(fkey):
                d = self.strip_density(edge_strip[u, v])

                if edges_to_curves:
                    polyline = []
                    if (u, v) in edges_to_curves:
                        curve = Polyline(edges_to_curves[u, v])
                        polyline = [curve.point_at(t) for t in linspace(0, 1, d)]
                    else:
                        curve = Polyline(edges_to_curves[v, u])
                        polyline = [curve.point_at(t) for t in linspace(0, 1, d)]
                        polyline[:] = polyline[::-1]
                else:
                    polyline = []
                    curve = Polyline(
                        [self.vertex_coordinates(u), self.vertex_coordinates(v)])
                    for i in range(0, d + 1):
                        point = curve.point_at(float(i) / float(d))
                        polyline.append(point)

                polylines.append(polyline)

            if self.is_face_pseudo_quad(fkey):
                pole = self.attributes['face_pole'][fkey]
                idx = self.face_vertices(fkey).index(pole)
                polylines.insert(idx, None)

            ab, bc, cd, da = polylines

            dc = cd[::-1] if cd else None
            ad = da[::-1] if da else None
            
            vertices, faces = discrete_coons_patch(ad, dc, bc, ab)
            faces = [[u for u, v in pairwise(face + face[:1]) if u != v] for face in faces]

            mesh = PseudoQuadMesh.from_vertices_and_faces_with_face_poles(vertices, faces)
            face_meshes[fkey] = mesh
            fvkeys = [edge[0] for edge in self.face_halfedges(fkey)]
            if len(fvkeys) == 4:
                a, b, c, d = fvkeys
            elif len(fvkeys) == 3:
                u, v, w = fvkeys
                if ab is None:
                    a, b, c, d = u, u, v, w
                elif bc is None:
                    a, b, c, d = u, v, v, w
                elif cd is None:
                    a, b, c, d = u, v, w, w
                elif da is None:
                    a, b, c, d = u, v, w, u
            n = self.strip_density(edge_strip[(c, d)]) + 1 if ab is None else self.strip_density(edge_strip[(a, b)]) + 1
            m = self.strip_density(edge_strip[(d, a)]) + 1 if bc is None else self.strip_density(edge_strip[(b, c)]) + 1
            parent2child[(a, b)] = (fkey, list(range(n)))
            parent2child[(b, c)] = (fkey, [n - 1 + k * n for k in range(m)])
            parent2child[(c, d)] = (fkey, list(reversed(list(range(len(vertices) - n, len(vertices))))))
            parent2child[(d, a)] = (fkey, list(reversed([k * n for k in range(m)])))
        
        vkeys_flatt = {}
        vcount = 0
        for fkey, fmesh in face_meshes.items():
            for i, vkey in enumerate(fmesh.vertices()):
                vkeys_flatt[(fkey, vkey)] = vcount + i
            vcount += fmesh.number_of_vertices()
        
        weld_map = {}
        for (u, v), (fkey, vkeys) in parent2child.items():
            for i, vkey in enumerate(vkeys):
                fkey_s, vkey_s = fkey, vkey
                if u == v or (i == 0 or i == len(vkeys) - 1):
                    end = v if i == len(vkeys) - 1 else u
                    min_fkey = min(self.vertex_faces(end))
                    nbr = self.face_vertex_descendant(min_fkey, end)
                    fkey2, vkeys2 = parent2child[(end, nbr)]
                    (fkey_s, vkey_s) = (fkey2, vkeys2[0])
                elif not self.is_edge_on_boundary((u, v)):
                    alt_fkey, alt_vkeys = parent2child[(v, u)]
                    if alt_fkey < fkey:
                        (fkey_s, vkey_s) = (alt_fkey, alt_vkeys[len(alt_vkeys) - 1 - i])
                weld_map[vkeys_flatt[(fkey, vkey)]] = (vkeys_flatt[(fkey_s, vkey_s)], face_meshes[fkey_s].vertex_attributes(vkey_s)) 

        for (fkey, vkey), vkey2 in vkeys_flatt.items():
            if vkey2 not in weld_map:
                weld_map[vkey2] = (vkey2, face_meshes[fkey].vertex_attributes(vkey))

        dense_mesh = PseudoQuadMesh()
        for fkey, fmesh in face_meshes.items():
            for vkey in fmesh.vertices():
                vkey2, attr = weld_map[vkeys_flatt[(fkey, vkey)]]
                if not dense_mesh.has_vertex(vkey2):
                    dense_mesh.add_vertex(key=vkey2, attr_dict=attr)
        for fkey, fmesh in face_meshes.items():
            for subfkey in fmesh.faces():
                dense_mesh.add_face([weld_map[vkeys_flatt[(fkey, vkey)]][0] for vkey in fmesh.face_vertices(subfkey)])
        self.dense_mesh(dense_mesh)

        return self.dense_mesh()
