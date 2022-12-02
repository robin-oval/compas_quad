from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from compas.datastructures import mesh_substitute_vertex_in_faces

__all__ = ['add_strip_lizard']


def add_strip_lizard(mesh, lizard, movements):
    """ Add a strip in mesh following a sequence of movements.

    Inputs
    ------
    mesh : Mesh
        The mesh to modify.
    lizard : tuple
        The initial position of the lizard as tuple of a triplet of vertex keys.
    movements : str
        A string of characters that dictate the movement of the lizard to add the strip: 't' for turns and 'p' for pivots.
    
    """

    def pivot(mesh, u, v):
        nbrs = mesh.vertex_neighbors(u, ordered=True)
        i = nbrs.index(v)
        w = nbrs[i + 1 - len(nbrs)]
        return u, w

    def turn(mesh, u, v):
        nbrs = mesh.vertex_neighbors(v, ordered=True)
        i = nbrs.index(u)
        w = nbrs[i + 1 - len(nbrs)]
        return v, w

    def vertex_neighbor_after(mesh, vkey, nbr):
        nbrs = mesh.vertex_neighbors(vkey, ordered=True)
        i = nbrs.index(nbr)
        return nbrs[i + 1 - len(nbrs)]

    def vertex_neighbor_before(mesh, vkey, nbr):
        nbrs = mesh.vertex_neighbors(vkey, ordered=True)
        i = nbrs.index(nbr)
        return nbrs[i - 1]

    def vertex_neighbors_between(mesh, vkey, from_nbr, to_nbr):
        # include from_nbr and exclude to_nbr
        nbrs = mesh.vertex_neighbors(vkey, ordered=True)
        i = nbrs.index(from_nbr)
        j = nbrs.index(to_nbr)
        return nbrs[j:] + nbrs[:i] if i < j else nbrs[j:i]

    def add_vertex_on_edge(mesh, u, v, t=0.5):
        xyz = mesh.edge_point(u, v, t=t)
        return mesh.add_vertex(attr_dict={key: value for key, value in zip('xyz', xyz)})

    def list_replace_item(items, old_item, new_items):
        i = items.index(old_item)
        del items[i]
        for item in reversed(new_items):
            items.insert(i, item)
        return items

    # to add: modify last tri, possibility to startinside with tri, closing strip

    tail, body, head = lizard

    for operation in movements:

        # pivot
        if operation == 'p':
            body, head = pivot(mesh, body, head)

        # turn
        if operation == 't':

            # store next position
            next_body, next_head = turn(mesh, body, head)

            # add new vertices on left and right sides
            body_left = add_vertex_on_edge(mesh, body, vertex_neighbor_before(mesh, body, head), t=0.1)
            body_right = add_vertex_on_edge(mesh, body, vertex_neighbor_after(mesh, body, head), t=0.1)

            # add tri face
            mesh.add_face([head, body_left, body_right])

            # sort left and right neighbors and faces except for previous quad face
            faces_left = [mesh.halfedge[nbr][body] for nbr in vertex_neighbors_between(mesh, body, head, tail)[1:]]            
            faces_right = [mesh.halfedge[nbr][body] for nbr in vertex_neighbors_between(mesh, body, tail, head)]
            
            # substitute vertices in left and right faces
            mesh_substitute_vertex_in_faces(mesh, body, body_left, fkeys=[fkey for fkey in faces_left if fkey is not None])
            mesh_substitute_vertex_in_faces(mesh, body, body_right, fkeys=[fkey for fkey in faces_right if fkey is not None])
            
            # convert previous tri face into quad face
            fkey = mesh.halfedge[tail][body]
            face_vertices = list_replace_item(mesh.face_vertices(fkey).copy(), body, [body_right, body_left])
            mesh.delete_face(fkey)
            mesh.add_face(face_vertices, fkey=fkey)

            # delete old vertex
            mesh.delete_vertex(body)
            
            # update lizard
            tail, body, head = body_right, next_body, next_head


if __name__ == '__main__':

    from compas_quad.datastructures import Mesh

    from compas_quad.datastructures import QuadMesh, CoarseQuadMesh

    from compas_view2.app import App

    # vertices = [
    #     [0.0, 0.0, 0.0],
    #     [0.5, 0.0, 0.0],
    #     [1.5, 0.0, 0.0],
    #     [2.0, 0.0, 0.0],
    #     [0.0, 1.0, 0.0],
    #     [1.0, 1.0, 0.0],
    #     [2.0, 1.0, 0.0],
    #     [0.0, 2.0, 0.0],
    #     [1.0, 2.0, 0.0],
    #     [2.0, 2.0, 0.0],
    # ]

    # faces = [
    #     [0, 1, 5, 4],
    #     [1, 2, 5],
    #     [2, 3, 6, 5],
    #     [4, 5, 8, 7],
    #     [5, 6, 9, 8],
    # ]

    # mesh = Mesh.from_vertices_and_faces(vertices, faces)


    vertices = [[0.5, 0.5, 0.0], [-0.5, 0.5, 0.0],
                [-0.5, -0.5, 0.0], [0.5, -0.5, 0.0]]
    faces = [[0, 1, 2, 3]]
    coarse = CoarseQuadMesh.from_vertices_and_faces(vertices, faces)

    # denser mesh
    coarse.collect_strips()
    coarse.strips_density(2)
    coarse.densification()
    mesh = coarse.dense_mesh()
    mesh.collect_strips()

    for vkey in mesh.vertices_on_boundary():
    if mesh.vertex_degree(vkey) == 2:
        body = vkey
        tail, head = [nbr for nbr in mesh0.vertex_neighbors(vkey) if mesh0.is_vertex_on_boundary(nbr)]
    break
    
    add_strip_lizard(mesh, (2, 5, 8), 'ttttptpt')

    # mesh.to_json('C:/Users/robin/OneDrive/Bureau/tmp.json')

    viewer = App(width=1600, height=900)
    viewer.add(mesh)
    viewer.show()
