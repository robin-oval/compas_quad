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

        pole_map = [geometric_key(self.vertex_coordinates(pole))
                    for pole in self.poles()]

        parent2child = {}
        # parent_poles = {pole: [] for pole in self.poles()}
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

            if cd:
                dc = cd[::-1]
            else:
                dc = None

            if da:
                ad = da[::-1]
            else:
                ad = None

            
            # vertices, faces = discrete_coons_patch(ab, bc, dc, ad)
            vertices, faces = discrete_coons_patch(ad, dc, bc, ab)
            faces = [[u for u, v in pairwise(face + face[:1]) if u != v] for face in faces]

            mesh = PseudoQuadMesh.from_vertices_and_faces_with_face_poles(vertices, faces)
            face_meshes[fkey] = mesh
            fvkeys = [edge[0] for edge in self.face_halfedges(fkey)]
            if len(fvkeys) == 4:
                u, v, w, x = fvkeys
                n, m = self.strip_density(edge_strip[(u, v)]) + 1, self.strip_density(edge_strip[(v, w)]) + 1
                parent2child[(u, v)] = (fkey, list(range(n)))
                parent2child[(v, w)] = (fkey, [n - 1 + k * n for k in range(m)])
                parent2child[(w, x)] = (fkey, list(reversed(list(range(len(vertices) - n, len(vertices))))))
                parent2child[(x, u)] = (fkey, list(reversed([k * n for k in range(m)])))
            elif len(fvkeys) == 3:
                u, v, w = fvkeys
                if ab is None:
                    n, m = self.strip_density(edge_strip[(v, w)]) + 1, self.strip_density(edge_strip[(u, v)]) + 1
                    parent2child[(u, u)] = (fkey, list(range(n)))
                    # parent_poles[u].append((fkey, list(range(n))))
                    parent2child[(u, v)] = (fkey, [n - 1 + k * n for k in range(m)])
                    parent2child[(v, w)] = (fkey, list(reversed(list(range(len(vertices) - n, len(vertices))))))
                    parent2child[(w, u)] = (fkey, list(reversed([k * n for k in range(m)])))
                elif bc is None:
                    n, m = self.strip_density(edge_strip[(u, v)]) + 1, self.strip_density(edge_strip[(w, u)]) + 1
                    parent2child[(u, v)] = (fkey, list(range(n)))
                    parent2child[(v, v)] = (fkey, [n - 1 + k * n for k in range(m)])
                    # parent_poles[v].append((fkey, [n - 1 + k * n for k in range(m)]))
                    parent2child[(v, w)] = (fkey, list(reversed(list(range(len(vertices) - n, len(vertices))))))
                    parent2child[(w, u)] = (fkey, list(reversed([k * n for k in range(m)])))
                elif cd is None:
                    n, m = self.strip_density(edge_strip[(u, v)]) + 1, self.strip_density(edge_strip[(v, w)]) + 1
                    parent2child[(u, v)] = (fkey, list(range(n)))
                    parent2child[(v, w)] = (fkey, [n - 1 + k * n for k in range(m)])
                    parent2child[(w, w)] = (fkey, list(reversed(list(range(len(vertices) - n, len(vertices))))))
                    # parent_poles[w].append((fkey, list(reversed(list(range(len(vertices) - n, len(vertices)))))))
                    parent2child[(w, u)] = (fkey, list(reversed([k * n for k in range(m)])))
                elif da is None:
                    n, m = self.strip_density(edge_strip[(u, v)]) + 1, self.strip_density(edge_strip[(v, w)]) + 1
                    parent2child[(u, v)] = (fkey, list(range(n)))
                    parent2child[(v, w)] = (fkey, [n - 1 + k * n for k in range(m)])
                    parent2child[(w, u)] = (fkey, list(reversed(list(range(len(vertices) - n, len(vertices))))))
                    parent2child[(u, u)] = (fkey, list(reversed([k * n for k in range(m)])))
                    # parent_poles[u].append((fkey, list(reversed([k * n for k in range(m)]))))
                else:
                    print('missing collapsed edge', polylines)
            else:
                print('wrong number of vertices')
        
        # face_pole_map = {}
        # for fkey, mesh in face_meshes.items():
        #     for fkey in mesh.faces():
        #         for u, v in pairwise(mesh.face_vertices(fkey) + mesh.face_vertices(fkey)[: 1]):
        #             if geometric_key(mesh.vertex_coordinates(u)) in pole_map and geometric_key(mesh.vertex_coordinates(u)) == geometric_key(mesh.vertex_coordinates(v)):
        #                 face_pole_map[geometric_key(mesh.face_center(fkey))] = geometric_key(
        #                     mesh.vertex_coordinates(u))
        #                 break
        
        # dense_mesh = PseudoQuadMesh()
        # for mesh in meshes:
        #     dense_mesh.join(mesh)
        # dense_mesh.weld()
        # self.dense_mesh(dense_mesh)
        # # self.dense_mesh(meshes_join_and_weld(meshes))


        vkeys_flatt = {}
        vcount = 0
        for fkey, fmesh in face_meshes.items():
            for i, vkey in enumerate(fmesh.vertices()):
                vkeys_flatt[(fkey, vkey)] = vcount + i
            vcount += fmesh.number_of_vertices()
        # print('vkeys_flatt', vkeys_flatt)
        
        weld_map = {}
        for (u, v), (fkey, vkeys) in parent2child.items():
            if u != v:
                for i, vkey in enumerate(vkeys):
                    if i == 0:
                        min_fkey = min(self.vertex_faces(u))
                        nbr = self.face_vertex_descendant(min_fkey, u)
                        fkey2, vkeys2 = parent2child[(u, nbr)]
                        weld_map[vkeys_flatt[(fkey, vkey)]] = (vkeys_flatt[(fkey2, vkeys2[0])], face_meshes[fkey2].vertex_attributes(vkeys2[0]))
                    elif i == len(vkeys) - 1:
                        min_fkey = min(self.vertex_faces(v))
                        nbr = self.face_vertex_descendant(min_fkey, v)
                        fkey2, vkeys2 = parent2child[(v, nbr)]
                        weld_map[vkeys_flatt[(fkey, vkey)]] = (vkeys_flatt[(fkey2, vkeys2[0])], face_meshes[fkey2].vertex_attributes(vkeys2[0]))
                    else:
                        if self.is_edge_on_boundary((u, v)):
                            weld_map[vkeys_flatt[(fkey, vkey)]] = (vkeys_flatt[(fkey, vkey)], face_meshes[fkey].vertex_attributes(vkey))
                        else:
                            alt_fkey, alt_vkeys = parent2child[(v, u)]
                            if alt_fkey > fkey:
                                weld_map[vkeys_flatt[(fkey, vkey)]] = (vkeys_flatt[(fkey, vkey)], face_meshes[fkey].vertex_attributes(vkey))
                            else:
                                weld_map[vkeys_flatt[(fkey, vkey)]] = (vkeys_flatt[(alt_fkey, alt_vkeys[len(alt_vkeys) - 1 - i])], face_meshes[alt_fkey].vertex_attributes(alt_vkeys[len(alt_vkeys) - 1 - i]))
            if u == v:
                for i, vkey in enumerate(vkeys):
                    min_fkey = min(self.vertex_faces(u))
                    nbr = self.face_vertex_descendant(min_fkey, u)
                    fkey2, vkeys2 = parent2child[(u, nbr)]
                    weld_map[vkeys_flatt[(fkey, vkey)]] = (vkeys_flatt[(fkey2, vkeys2[0])], face_meshes[fkey2].vertex_attributes(vkeys2[0]))

        # print('weld_map', weld_map)
        # print(parent_poles)
        # for pole, data in parent_poles.items():
        #     new_vkeys = [vkeys_flatt[(fkey, vkey)] for fkey, vkeys in data for vkey in vkeys]
        #     new_vkey = min(new_vkeys)
        #     for fkey, vkeys in data:
        #         for vkey in vkeys:
        #             # print('attr', )
        #             weld_map[vkeys_flatt[(fkey, vkey)]] = (new_vkey, self.vertex_attributes(pole))
        
        # print('weld_map', weld_map)
        for (fkey, vkey), vkey2 in vkeys_flatt.items():
            if vkey2 not in weld_map:
                weld_map[vkey2] = (vkey2, face_meshes[fkey].vertex_attributes(vkey))
        # print('weld_map', weld_map)
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

        # face_pole = {}
        # for fkey

        # face_pole = {}
        # for fkey in self.dense_mesh().faces():
        #     if len(self.dense_mesh().face_vertices(fkey)) < 3:
        #         print(self.dense_mesh().face_vertices(fkey))
        #     if geometric_key(self.dense_mesh().face_center(fkey)) in face_pole_map:
        #         for vkey in self.dense_mesh().face_vertices(fkey):
        #             if geometric_key(self.dense_mesh().vertex_coordinates(vkey)) == face_pole_map[geometric_key(self.dense_mesh().face_center(fkey))]:
        #                 face_pole[fkey] = vkey
        #                 break

        # self.dense_mesh().attributes['face_pole'] = face_pole
        return self.dense_mesh()
