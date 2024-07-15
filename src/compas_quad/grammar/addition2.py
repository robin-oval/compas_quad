from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from compas_quad.datastructures import mesh_substitute_vertex_in_faces

__all__ = ['lizard_atp', 'add_strip_lizard', 'add_strip_lizard_2']


def lizard_atp(mesh, lizard, movements):

    adding = False
    to_add_movements = ''

    while len(movements) > 0:

        operation, movements = movements[0], movements[1:]

        if operation == 'a':
            if adding:
                lizard = add_strip_lizard_2(mesh, lizard, to_add_movements)
                to_add_movements = ''
            adding = not adding
            
        else:
            if adding:
                to_add_movements += operation
            else:
                if operation == 'p':
                    lizard = lizard_pivot(mesh, *lizard)
                elif operation == 't':
                    lizard = lizard_turn(mesh, *lizard)
    
    return lizard


def lizard_pivot(mesh, tail, body, head):
    nbrs = mesh.vertex_neighbors(body, ordered=True)
    i = nbrs.index(head)
    new_head = nbrs[i + 1 - len(nbrs)]
    return tail, body, new_head


def lizard_turn(mesh, tail, body, head):
    nbrs = mesh.vertex_neighbors(head, ordered=True)
    i = nbrs.index(body)
    new_head = nbrs[i + 1 - len(nbrs)]
    return body, head, new_head


def add_strip_lizard_2(mesh, lizard, movements):
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
        xyz = mesh.edge_point((u, v), t=t)
        return mesh.add_vertex(attr_dict={key: value for key, value in zip('xyz', xyz)})

    def list_replace_item(items, old_item, new_items):
        i = items.index(old_item)
        del items[i]
        for item in reversed(new_items):
            items.insert(i, item)
        return items

    # to add: modify last tri, possibility to startinside with tri, closing strip

    tail, body, head = lizard

    new_fkeys = []

    movements_copy = movements

    while len(movements) > 0:

        operation, movements = movements[0], movements[1:]

        # pivot
        if operation == 'p':
            body, head = pivot(mesh, body, head)

        # turn
        if operation == 't':

            # store next position
            next_body, next_head = turn(mesh, body, head)
            # add new vertices on left and right sides
            if body != head:
                t1, t2 = 0.1, 0.1
            else:
                t1, t2 = 0.1, 0.2
            body_left = add_vertex_on_edge(mesh, body, vertex_neighbor_before(mesh, body, head), t=0.1)
            body_right = add_vertex_on_edge(mesh, body, vertex_neighbor_after(mesh, body, head), t=0.1)

            # add tri face
            new_fkeys.append(mesh.add_face([head, body_left, body_right]))

            # sort left and right neighbors and faces except for previous quad face
            if tail != head:
                faces_left = [mesh.halfedge[nbr][body] for nbr in vertex_neighbors_between(mesh, body, head, tail)[1:]]            
                faces_right = [mesh.halfedge[nbr][body] for nbr in vertex_neighbors_between(mesh, body, tail, head)]
            else:       
                fkeys = mesh.vertex_faces(body).copy()
                if body in mesh.halfedge[tail]:
                    fkey = mesh.halfedge[tail][body]
                    if fkey is not None:
                        fkeys.remove(fkey)
                faces_left, faces_right = fkeys, []

            # substitute vertices in left and right faces
            mesh_substitute_vertex_in_faces(mesh, body, body_left, fkeys=[fkey for fkey in faces_left if fkey is not None])
            mesh_substitute_vertex_in_faces(mesh, body, body_right, fkeys=[fkey for fkey in faces_right if fkey is not None])

            # convert previous tri face into quad face
            if body in mesh.halfedge[tail]:
                fkey = mesh.halfedge[tail][body]
                face_vertices = list_replace_item(mesh.face_vertices(fkey).copy(), body, [body_right, body_left])
                mesh.delete_face(fkey)
                mesh.add_face(face_vertices, fkey=fkey)

            # delete old vertex
            mesh.delete_vertex(body)
            
            # update lizard
            tail, body, head = body_right, next_body, next_head
        
    if len(new_fkeys) > 0 and mesh.is_vertex_on_boundary(body):
        
        # close strip if possible
        if mesh.halfedge[head][body] == new_fkeys[0]:
            # print('!')
            body_left, body_right = head, body
            fkey = new_fkeys[-1]
            face_vertices = list_replace_item(mesh.face_vertices(fkey).copy(), body, [body_right, body_left])
            mesh.delete_face(fkey)
            mesh.add_face(face_vertices, fkey=fkey)

        else:
            nbrs = mesh.vertex_neighbors(body, ordered=True)
            body_left = add_vertex_on_edge(mesh, body, nbrs[-1], t=0.1)
            body_right = add_vertex_on_edge(mesh, body, nbrs[0], t=0.1)
            i = nbrs.index(tail)
            faces_left = [mesh.halfedge[nbr][body] for nbr in nbrs[i:][1:]]            
            faces_right = [mesh.halfedge[nbr][body] for nbr in nbrs[:i]]
            mesh_substitute_vertex_in_faces(mesh, body, body_left, fkeys=[fkey for fkey in faces_left if fkey is not None])
            mesh_substitute_vertex_in_faces(mesh, body, body_right, fkeys=[fkey for fkey in faces_right if fkey is not None])        
            fkey = new_fkeys[-1]
            face_vertices = list_replace_item(mesh.face_vertices(fkey).copy(), body, [body_right, body_left])
            mesh.delete_face(fkey)
            mesh.add_face(face_vertices, fkey=fkey)
            mesh.delete_vertex(body)

    return tail, body_right, head


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
        xyz = mesh.edge_point((u, v), t=t)
        return mesh.add_vertex(attr_dict={key: value for key, value in zip('xyz', xyz)})

    def list_replace_item(items, old_item, new_items):
        i = items.index(old_item)
        del items[i]
        for item in reversed(new_items):
            items.insert(i, item)
        return items

    # to add: modify last tri, possibility to startinside with tri, closing strip

    tail, body, head = lizard

    new_fkeys = []

    movements_copy = movements

    while len(movements) > 0:

        operation, movements = movements[0], movements[1:]

        # pivot
        if operation == '0':
            body, head = pivot(mesh, body, head)

        # turn
        if operation == '1':

            # store next position
            next_body, next_head = turn(mesh, body, head)
            # add new vertices on left and right sides
            if body != head:
                t1, t2 = 0.1, 0.1
            else:
                t1, t2 = 0.1, 0.2
            body_left = add_vertex_on_edge(mesh, body, vertex_neighbor_before(mesh, body, head), t=0.1)
            body_right = add_vertex_on_edge(mesh, body, vertex_neighbor_after(mesh, body, head), t=0.1)

            # add tri face
            new_fkeys.append(mesh.add_face([head, body_left, body_right]))

            # sort left and right neighbors and faces except for previous quad face
            if tail != head:
                faces_left = [mesh.halfedge[nbr][body] for nbr in vertex_neighbors_between(mesh, body, head, tail)[1:]]            
                faces_right = [mesh.halfedge[nbr][body] for nbr in vertex_neighbors_between(mesh, body, tail, head)]
            else:       
                fkeys = mesh.vertex_faces(body).copy()
                if body in mesh.halfedge[tail]:
                    fkey = mesh.halfedge[tail][body]
                    if fkey is not None:
                        fkeys.remove(fkey)
                faces_left, faces_right = fkeys, []

            # substitute vertices in left and right faces
            mesh_substitute_vertex_in_faces(mesh, body, body_left, fkeys=[fkey for fkey in faces_left if fkey is not None])
            mesh_substitute_vertex_in_faces(mesh, body, body_right, fkeys=[fkey for fkey in faces_right if fkey is not None])

            # convert previous tri face into quad face
            if body in mesh.halfedge[tail]:
                fkey = mesh.halfedge[tail][body]
                face_vertices = list_replace_item(mesh.face_vertices(fkey).copy(), body, [body_right, body_left])
                mesh.delete_face(fkey)
                mesh.add_face(face_vertices, fkey=fkey)

            # delete old vertex
            mesh.delete_vertex(body)
            
            # update lizard
            tail, body, head = body_right, next_body, next_head
        
    if len(new_fkeys) > 0 and mesh.is_vertex_on_boundary(body):
        
        # close strip if possible
        if mesh.halfedge[head][body] == new_fkeys[0]:
            # print('!')
            body_left, body_right = head, body
            fkey = new_fkeys[-1]
            face_vertices = list_replace_item(mesh.face_vertices(fkey).copy(), body, [body_right, body_left])
            mesh.delete_face(fkey)
            mesh.add_face(face_vertices, fkey=fkey)

        else:
            # print('!!')
            nbrs = mesh.vertex_neighbors(body, ordered=True)
            body_left = add_vertex_on_edge(mesh, body, nbrs[-1], t=0.1)
            body_right = add_vertex_on_edge(mesh, body, nbrs[0], t=0.1)
            i = nbrs.index(tail)
            faces_left = [mesh.halfedge[nbr][body] for nbr in nbrs[i:][1:]]            
            faces_right = [mesh.halfedge[nbr][body] for nbr in nbrs[:i]]
            mesh_substitute_vertex_in_faces(mesh, body, body_left, fkeys=[fkey for fkey in faces_left if fkey is not None])
            mesh_substitute_vertex_in_faces(mesh, body, body_right, fkeys=[fkey for fkey in faces_right if fkey is not None])        
            fkey = new_fkeys[-1]
            face_vertices = list_replace_item(mesh.face_vertices(fkey).copy(), body, [body_right, body_left])
            mesh.delete_face(fkey)
            mesh.add_face(face_vertices, fkey=fkey)
            mesh.delete_vertex(body)

    return tail, body, head


if __name__ == '__main__':

    from math import pi, cos, sin

    from itertools import product

    from compas_quad.datastructures import Mesh

    from compas_quad.datastructures import QuadMesh, CoarseQuadMesh, CoarsePseudoQuadMesh

    from compas.numerical import fd_numpy

    from compas_view2.app import App


    def postprocessing(mesh):

        key2index = mesh.key_index()
        index2key = mesh.index_key()

        # map boundary to circle
        fixed = [key2index[key] for key in mesh.vertices_on_boundary()[:-1]]
        n = len(fixed)
        for i, vidx in enumerate(fixed):
            vkey = index2key[vidx]
            attr = mesh.vertex[vkey]
            attr['x'] = 0.5 * cos(i / n * 2 * pi)
            attr['y'] = 0.5 * sin(i / n * 2 * pi)
            attr['z'] = 0

        # force density method
        vertices = [mesh.vertex_coordinates(vkey) for vkey in mesh.vertices()]
        edges = [(key2index[u], key2index[v]) for u, v in mesh.edges()]
        q = [1.0] * len(edges)
        loads = [[0.0, 0.0, 0.0]] * len(vertices)

        xyz, q, f, l, r = fd_numpy(vertices, edges, fixed, q, loads)

        for i, (x, y, z) in enumerate(xyz):
            vkey = index2key[i]
            attr = mesh.vertex[vkey]
            attr['x'] = x
            attr['y'] = y
            attr['z'] = z

    vertices = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [2.0, 1.0, 0.0],
        [0.0, 2.0, 0.0],
        [1.0, 2.0, 0.0],
        [2.0, 2.0, 0.0],
    ]

    faces = [
        [0, 1, 4, 3],
        [1, 2, 5, 4],
        [3, 4, 7, 6],
        [4, 5, 8, 7],
    ]

    mesh0 = Mesh.from_vertices_and_faces(vertices, faces)

    for vkey in mesh0.vertices_on_boundary():
        if mesh0.vertex_degree(vkey) == 2:
            body = vkey
            tail, head = [nbr for nbr in mesh0.vertex_neighbors(vkey) if mesh0.is_vertex_on_boundary(nbr)]
        break
    
    lizard = tail, body, head
    print('lizard', lizard)

    for string in product('01', repeat=5):
        print('string', string)
        
        mesh = mesh0.copy()
        tail, body, head = add_strip_lizard(mesh, lizard, string)

        print('manifold', mesh_is_manifold(mesh))
        print('quad', mesh.is_quadmesh())

        # if mesh.is_quadmesh():
        # vertices, faces = mesh.to_vertices_and_faces()
        # mesh = CoarseQuadMesh.from_vertices_and_faces(vertices, faces)
        # mesh = CoarsePseudoQuadMesh.from_vertices_and_faces_with_poles(vertices, faces, poles=[pole])
        # mesh.collect_strips()
        # mesh.strips_density(3)
        # mesh.densification()
        # mesh = mesh.dense_mesh()
        # postprocessing(mesh)

    mesh.to_json('C:/Users/robin/OneDrive/Bureau/tmp.json')

    viewer = App(width=1600, height=900)
    viewer.add(mesh)
    viewer.show()
