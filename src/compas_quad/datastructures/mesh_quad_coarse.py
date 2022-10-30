from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from copy import deepcopy
from math import floor, ceil

from compas.topology import adjacency_from_edges
from compas.topology import connected_components

from compas.geometry import discrete_coons_patch
from compas.geometry import vector_average

from compas.datastructures import meshes_join_and_weld

from compas_quad.datastructures import Mesh, QuadMesh

from compas.utilities import pairwise


__all__ = ['CoarseQuadMesh']


class CoarseQuadMesh(QuadMesh):

    def __init__(self):
        super(CoarseQuadMesh, self).__init__()
        self.attributes['strips_density'] = {}
        self.attributes['vertex_coarse_to_dense'] = {}
        self.attributes['edge_coarse_to_dense'] = {}
        self.attributes['dense_mesh'] = None
        self.attributes['pattern'] = None

    # --------------------------------------------------------------------------
    # constructors
    # --------------------------------------------------------------------------

    @classmethod
    def from_quad_mesh(cls, quad_mesh, collect_strips=True, collect_polyedges=True, attribute_density=True):
        """Build coarse quad mesh from quad mesh with density and child-parent element data.

        Parameters
        ----------
        quad_mesh : QuadMesh
            A quad mesh.
        attribute_density : bool, optional
            Keep density data of dense quad mesh and inherit it as aatribute.

        Returns
        ----------
        coarse_quad_mesh : CoarseQuadMesh
            A coarse quad mesh with density data.
        """
        polyedges = quad_mesh.singularity_polyedge_decomposition()

        # vertex data
        vertices = {vkey: quad_mesh.vertex_coordinates(
            vkey) for vkey in quad_mesh.vertices()}
        coarse_vertices_children = {
            vkey: vkey for polyedge in polyedges for vkey in [polyedge[0], polyedge[-1]]}
        coarse_vertices = {vkey: quad_mesh.vertex_coordinates(
            vkey) for vkey in coarse_vertices_children}

        # edge data
        coarse_edges_children = {
            (polyedge[0], polyedge[-1]): polyedge for polyedge in polyedges}
        singularity_edges = [(x, y) for polyedge in polyedges for u, v in pairwise(
            polyedge) for x, y in [(u, v), (v, u)]]

        # face data
        faces = {fkey: quad_mesh.face_vertices(
            fkey) for fkey in quad_mesh.faces()}
        adj_edges = {(f1, f2) for f1 in quad_mesh.faces() for f2 in quad_mesh.face_neighbors(
            f1) if f1 < f2 and quad_mesh.face_adjacency_halfedge(f1, f2) not in singularity_edges}
        coarse_faces_children = {}
        for i, connected_faces in enumerate(connected_components(adjacency_from_edges(adj_edges))):
            mesh = Mesh.from_vertices_and_faces(
                vertices, [faces[face] for face in connected_faces])
            coarse_faces_children[i] = [vkey for vkey in reversed(
                mesh.vertices_on_boundaries()[0]) if mesh.vertex_degree(vkey) == 2]

        coarse_quad_mesh = cls.from_vertices_and_faces(
            coarse_vertices, coarse_faces_children)

        # attribute relation child-parent element between coarse and dense quad meshes
        coarse_quad_mesh.attributes['vertex_coarse_to_dense'] = coarse_vertices_children
        coarse_quad_mesh.attributes['edge_coarse_to_dense'] = {
            u: {} for u in coarse_quad_mesh.vertices()}
        for (u, v), polyedge in coarse_edges_children.items():
            coarse_quad_mesh.attributes['edge_coarse_to_dense'][u][v] = polyedge
            coarse_quad_mesh.attributes['edge_coarse_to_dense'][v][u] = list(
                reversed(polyedge))

        # collect strip and polyedge attributes
        if collect_strips:
            coarse_quad_mesh.collect_strips()
        if collect_polyedges:
            coarse_quad_mesh.collect_polyedges()

        # store density attribute from input dense quad mesh
        if attribute_density:
            coarse_quad_mesh.strips_density(1)
            for skey in coarse_quad_mesh.strips():
                u, v = coarse_quad_mesh.strip_edges(skey)[0]
                d = len(coarse_edges_children.get(
                    (u, v), coarse_edges_children.get((v, u), [])))
                coarse_quad_mesh.strip_density(skey, d)

        # store quad mesh and use as polygonal mesh
        coarse_quad_mesh.dense_mesh(quad_mesh)
        coarse_quad_mesh.pattern(deepcopy(quad_mesh))

        return coarse_quad_mesh

    # --------------------------------------------------------------------------
    # meshes getters and setters
    # --------------------------------------------------------------------------

    def dense_mesh(self, quad_mesh=None):
        if quad_mesh:
            self.attributes['dense_mesh'] = quad_mesh
        return self.attributes['dense_mesh']

    def pattern(self, pattern=None):
        if pattern:
            self.attributes['pattern'] = pattern
        return self.attributes['pattern']

    # --------------------------------------------------------------------------
    # element child-parent relation getters
    # --------------------------------------------------------------------------

    def coarse_edge_dense_edges(self, u, v):
        """Return the child edges, or polyedge, in the dense quad mesh from a parent edge in the coarse quad mesh."""
        return self.attributes['edge_coarse_to_dense'][u][v]

    # --------------------------------------------------------------------------
    # density getters and setters
    # --------------------------------------------------------------------------

    def strip_density(self, skey, d=None):
        """Get or set the densty of one strip.

        Parameters
        ----------
        skey : hashable
            A strip key.
        d : int
            A density parameter.

        Returns
        -------
        int
            The density parameter.
        """
        if d:
            self.attributes['strips_density'][skey] = d
        return self.attributes['strips_density'][skey]

    def strips_density(self, d=None):
        """Get or set the density of all the strips.

        Parameters
        ----------
        d : int, optional
            A density parameter.

        Returns
        -------
        dict

        """
        if d:
            for skey in self.strips():
                self.strip_density(skey, d)
        return self.attributes['strips_density']

    def set_strip_density_target(self, skey, t):
        """Set the strip densities based on a target length and the average length of the strip edges.

        Parameters
        ----------
        skey : hashable
            A strip key.
        t : float
            A target length.
        """
        self.strip_density(skey, int(ceil(vector_average(
            [self.edge_length(u, v) for u, v in self.strip_edges(skey) if u != v]) / t)))

    def set_strip_density_func(self, skey, func, func_args):
        """Set the strip densities based on a function.

        Parameters
        ----------
        skey : hashable
            A strip key.
        """
        self.strip_density(skey, int(func(skey, func_args)))

    def set_strips_density_target(self, t, skeys=None):
        """Set the strip densities based on a target length and the average length of the strip edges.

        Parameters
        ----------
        t : float
            A target length.
        skeys : list, None
            The keys of strips to set density. If is None, all strips are considered.
        """
        if skeys is None:
            skeys = self.strips()
        for skey in skeys:
            self.set_strip_density_target(skey, t)

    def set_strips_density_func(self, func, func_args, skeys=None):
        """Set the strip densities based on a function.

        Parameters
        ----------
        skeys : list, None
            The keys of strips to set density. If is None, all strips are considered.
        """
        if skeys is None:
            skeys = self.strips()
        for skey in skeys:
            self.set_strip_density_func(skey, func, func_args)

    def set_mesh_density_face_target(self, nb_faces):
        """Set equal strip densities based on a target number of faces.

        Parameters
        ----------
        nb_faces : int
            The target number of faces.
        """
        n = (nb_faces / self.number_of_faces()) ** .5
        if ceil(n) - n > n - floor(n):
            n = int(floor(n))
        else:
            n = int(ceil(n))
        self.strips_density(n)

    # --------------------------------------------------------------------------
    # densification
    # --------------------------------------------------------------------------

    def densification(self):
        """Generate a denser quad mesh from the coarse quad mesh and its strip densities.
        """
        edge_strip = {}
        for skey, edges in self.strips(data=True):
            for edge in edges:
                edge_strip[edge] = skey
                edge_strip[tuple(reversed(edge))] = skey

        face_meshes = {}
        for fkey in self.faces():
            ab, bc, cd, da = [[self.edge_point(u, v, float(i) / float(self.strip_density(edge_strip[(u, v)])))
                               for i in range(0, self.strip_density(edge_strip[(u, v)]) + 1)] for u, v in self.face_halfedges(fkey)]
            vertices, faces = discrete_coons_patch(
                ab, bc, list(reversed(cd)), list(reversed(da)))
            face_meshes[fkey] = QuadMesh.from_vertices_and_faces(
                vertices, faces)

        self.dense_mesh(meshes_join_and_weld(list(face_meshes.values())))

    # def geometrical_densification(self):
    #   """Generate a denser quad mesh from the coarse quad mesh and its strip densities.

    #   WIP!

    #   Returns
    #   -------
    #   QuadMesh
    #       A denser quad mesh.

    #   """

    #   if self.quad_mesh is None:
    #       self.quad_mesh = self.copy()
    #       self.polygonal_mesh = self.copy()
    #       self.vertex_to_vertex = {vkey: vkey for vkey in self.vertices()}

    #       self.edge_to_polyedge = {vkey: {} for vkey in self.vertices()}
    #       for u, v in self.edges():
    #           self.edge_to_polyedge[u][v] = (u, v)
    #           self.edge_to_polyedge[v][u] = (v, u)

    #   quad_mesh = self.quad_mesh

    #   new_edge_polyline = {}

    #   for u, v in self.edges():
    #       d =  self.strip_density(self.edge_strip((u, v)))
    #       old_polyline = Polyline([quad_mesh.vertex_coordinates(vkey) for vkey in self.edge_to_polyedge[u][v]])
    #       new_polyline = [old_polyline.point(float(i) / float(d)) for i in range(0, d + 1)]
    #       new_edge_polyline[(u, v)] = new_polyline
    #       new_edge_polyline[(v, u)] = list(reversed(new_polyline))

    #   meshes = []

    #   new_edge_polyedge = {}
    #   new_vertex_vertex = {}

    #   for i, fkey in enumerate(self.faces()):
    #       ab, bc, cd, da = [new_edge_polyline[edge] for edge in self.face_halfedges(fkey)]

    #       a, b, c, d = self.face_vertices(fkey)
    #       n, m = len(ab), len(bc)
    #       vertices, faces = discrete_coons_patch(ab, bc, list(reversed(cd)), list(reversed(da)))
    #       meshes.append(QuadMesh.from_vertices_and_faces(vertices, faces))
    #       ab, bc, cd, da = list(self.face_halfedges(fkey))

    #       new_edge_polyedge[ab] = [i, range(0, m)]
    #       new_edge_polyedge[tuple(reversed(ab))] = [new_edge_polyedge[ab][0], list(reversed(new_edge_polyedge[ab][1]))]

    #       new_edge_polyedge[bc] = [i, range(m - 1, n * m, m)]
    #       new_edge_polyedge[tuple(reversed(bc))] = [new_edge_polyedge[bc][0], list(reversed(new_edge_polyedge[bc][1]))]

    #       new_edge_polyedge[cd] = [i, list(reversed(range((n - 1) * m, n * m)))]
    #       new_edge_polyedge[tuple(reversed(cd))] = [new_edge_polyedge[cd][0], list(reversed(new_edge_polyedge[cd][1]))]

    #       new_edge_polyedge[da] = [i, list(reversed(range(0, (n - 1) * m + 1, m)))]
    #       new_edge_polyedge[tuple(reversed(da))] = [new_edge_polyedge[da][0], list(reversed(new_edge_polyedge[da][1]))]

    #       new_vertex_vertex[a] = (new_edge_polyedge[ab][0], new_edge_polyedge[ab][1][0])
    #       new_vertex_vertex[b] = (new_edge_polyedge[ab][0], new_edge_polyedge[ab][1][-1])
    #       new_vertex_vertex[c] = (new_edge_polyedge[cd][0], new_edge_polyedge[cd][1][0])
    #       new_vertex_vertex[d] = (new_edge_polyedge[cd][0], new_edge_polyedge[cd][1][-1])

    #   self.quad_mesh, old_to_new_vertices = meshes_join_and_weld(meshes, data = True)

    #   self.vertex_to_vertex = {vkey: old_to_new_vertices[tuple(new_vertex_vertex[vkey])] for vkey in self.vertices()}

    #   for u, v in self.edges():
    #       i, vkeys = new_edge_polyedge[(u, v)]
    #       new_polyedge = [old_to_new_vertices[(i, vkey)] for vkey in vkeys]

    #       if self.vertex_to_vertex[u] == new_polyedge[0]:
    #           self.edge_to_polyedge[u][v] = new_polyedge
    #           self.edge_to_polyedge[v][u] = list(reversed(new_polyedge))

    #       elif self.vertex_to_vertex[u] == new_polyedge[-1]:
    #           self.edge_to_polyedge[u][v] = list(reversed(new_polyedge))
    #           self.edge_to_polyedge[v][u] = new_polyedge

    #       else:
    #           pass

    #   return self.quad_mesh
