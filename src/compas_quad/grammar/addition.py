from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from compas.topology import breadth_first_paths

from compas_quad.datastructures import mesh_substitute_vertex_in_faces

from compas.itertools import pairwise
from compas_quad.utilities import sublist_from_to_items_in_closed_list

__all__ = [
    'add_strip',
    'add_strips'
]


def add_strips(mesh, polyedges, callback=None, callback_args=None):
    to_add = polyedges[:]
    while len(to_add) > 0:
        polyedge = to_add.pop()
        add_strip(mesh, polyedge)
        # update polyedges
        if callback:
            if callable(callback):
                callback(mesh, callback_args)


def add_strip(mesh, polyedge):
    full_updated_polyedge = []
    # store data
    left_polyedge = []
    right_polyedge = []
    new_faces = []

    # exception if closed
    is_closed = polyedge[0] == polyedge[-1]
    if is_closed:
        polyedge.pop()

    k = -1
    count = len(polyedge) * 2
    while count and len(polyedge) > 0:
        k += 1
        count -= 1

        # select u, v, w if not closed
        if not is_closed:
            # u
            if len(new_faces) != 0:
                u1, u2 = left_polyedge[-1], right_polyedge[-1]
            else:
                u1, u2 = None, None
            # v
            v = polyedge.pop(0)
            full_updated_polyedge.append(v)
            # w
            if len(polyedge) != 0:
                w = polyedge[0]
            else:
                w = None

        # select u, v, w if closed
        else:
            # u
            if len(new_faces) != 0:
                u1, u2 = left_polyedge[-1], right_polyedge[-1]
            else:
                u1 = polyedge[-1]  # artificial u1
            # v
            v = polyedge.pop(0)
            full_updated_polyedge.append(v)
            # w
            if len(polyedge) != 0:
                w = polyedge[0]
            else:
                w = left_polyedge[0]

        # add new vertices
        faces = sort_faces(mesh, u1, v, w)
        v1, v2 = mesh.add_vertex(attr_dict=mesh.vertex[v]), mesh.add_vertex(
            attr_dict=mesh.vertex[v])

        if type(faces[0]) == list:
            faces_1, faces_2 = faces
        else:
            # exception necessary for U-turns
            if faces[0] in mesh.vertex_faces(left_polyedge[-2]):
                faces_1 = faces
                faces_2 = []
            else:
                faces_1 = []
                faces_2 = faces
        mesh_substitute_vertex_in_faces(mesh, v, v1, faces_1)
        mesh_substitute_vertex_in_faces(mesh, v, v2, faces_2)
        mesh.delete_vertex(v)
        left_polyedge.append(v1)
        right_polyedge.append(v2)

        # add new faces, different if at the start, end or main part of the polyedge
        if len(new_faces) == 0:
            if not is_closed:
                new_faces.append(mesh.add_face([v1, w, v2]))
            else:
                new_faces.append(mesh.add_face([v1, v2, u1]))
                new_faces.append(mesh.add_face([v1, w, v2]))
        elif len(polyedge) == 0:
            if not is_closed:
                u1, u2 = left_polyedge[-2], right_polyedge[-2]
                face = new_faces.pop()
                mesh.delete_face(face)
                new_faces.append(mesh.add_face([u1, v1, v2, u2]))
            else:
                u1, u2 = left_polyedge[-2], right_polyedge[-2]
                face = new_faces.pop()
                mesh.delete_face(face)
                new_faces.append(mesh.add_face([v1, u1, u2, v2]))
                face = new_faces.pop(0)
                mesh.delete_face(face)
                u1, u2 = left_polyedge[0], right_polyedge[0]
                new_faces.append(mesh.add_face([v1, u1, u2, v2]))

                mesh_substitute_vertex_in_faces(mesh, v, v1)
                mesh_substitute_vertex_in_faces(mesh, v, v2)
        else:
            face = new_faces.pop()
            mesh.delete_face(face)
            new_faces.append(mesh.add_face([u1, v1, v2, u2]))
            new_faces.append(mesh.add_face([v1, w, v2]))

        # update
        updated_polyedge = []
        via_vkeys = [v1, v2]
        for i, vkey in enumerate(polyedge):
            if vkey != v:
                updated_polyedge.append(vkey)
            else:
                from_vkey = polyedge[i - 1]
                to_vkey = polyedge[i + 1]
                updated_polyedge += polyedge_from_to_via_vertices(
                    mesh, from_vkey, to_vkey, via_vkeys)[1:-1]
        polyedge = updated_polyedge

    # include pseudo closed polyedges

    old_vkeys_to_new_vkeys = {u0: (u1, u2) for u0, u1, u2 in zip(
        full_updated_polyedge, left_polyedge, right_polyedge)}

    # for fkey in mesh.faces():
    #    print(mesh.face_vertices(fkey))
    n = update_strip_data(mesh, full_updated_polyedge, old_vkeys_to_new_vkeys)
    # print(left_polyedge, right_polyedge)
    return n, old_vkeys_to_new_vkeys


def add_strip_old(mesh, polyedge):
    """Add a strip along a mesh polyedge.

    Parameters
    ----------
    mesh : Mesh
        A mesh.
    polyedge : list
        List of vertex keys forming path.

    Returns
    -------
    new_skey, left_polyedge, right_polyedge : tuple
        The key of the new strip, the new strip vertices on the left, the new strip vertices on the right.

    """

    # close or open status
    closed = polyedge[0] == polyedge[-1]

    # store transversal strips to update later
    update = {mesh.edge_strip(edge): i for i,
              edge in enumerate(pairwise(polyedge))}
    transverse_strips = set(update.keys())

    # list faces on the left and right of the polyedge
    left_faces = [mesh.halfedge[u][v] for u, v in pairwise(polyedge)]
    right_faces = [mesh.halfedge[v][u] for u, v in pairwise(polyedge)]

    # add extremities for looping on data
    if closed:
        left_faces = [left_faces[-1]] + left_faces + [left_faces[0]]
        right_faces = [right_faces[-1]] + right_faces + [right_faces[0]]
    else:
        left_faces = [None] + left_faces + [None]
        right_faces = [None] + right_faces + [None]

    # remove duplicat extremity
    if closed:
        polyedge.pop()

    # duplicate polyedge
    left_polyedge = [mesh.add_vertex(
        attr_dict=mesh.vertex[vkey]) for vkey in polyedge]
    right_polyedge = [mesh.add_vertex(
        attr_dict=mesh.vertex[vkey]) for vkey in polyedge]

    # store changes to apply all at once later
    to_substitute = {vkey: [] for vkey in polyedge}

    all_left_faces = []
    all_right_faces = []
    # collect all faces to update along polyedge with corresponding new vertex
    for i, vkey in enumerate(polyedge):
        vertex_faces = mesh.vertex_faces(vkey, ordered=True, include_none=True)
        # on the left
        faces = sublist_from_to_items_in_closed_list(
            vertex_faces, left_faces[i], left_faces[i + 1])
        all_left_faces += faces
        to_substitute[vkey].append((left_polyedge[i], faces))
        # on the right
        faces = sublist_from_to_items_in_closed_list(
            vertex_faces, right_faces[i + 1], right_faces[i])
        all_right_faces += faces
        to_substitute[vkey].append((right_polyedge[i], faces))

    all_left_faces = list(set(all_left_faces))
    all_right_faces = list(set(all_right_faces))
    left_strips = list(set([skey for fkey in all_left_faces if fkey is not None for skey in mesh.face_strips(
        fkey) if skey not in transverse_strips]))
    right_strips = list(set([skey for fkey in all_right_faces if fkey is not None for skey in mesh.face_strips(
        fkey) if skey not in transverse_strips]))

    # apply changes
    for key, substitutions in to_substitute.items():
        for substitution in substitutions:
            new_key, faces = substitution
            mesh_substitute_vertex_in_faces(
                mesh, key, new_key, [face for face in faces if face is not None])

    # delete old vertices
    for vkey in polyedge:
        mesh.delete_vertex(vkey)

    # add strip faces
    if closed:
        polyedge.append(polyedge[0])
        left_polyedge.append(left_polyedge[0])
        right_polyedge.append(right_polyedge[0])
    for i in range(len(polyedge) - 1):
        mesh.add_face([right_polyedge[i], right_polyedge[i + 1],
                       left_polyedge[i + 1], left_polyedge[i]])

    # update transverse strip data
    for skey, i in update.items():
        mesh.attributes['strips'][skey] = mesh.collect_strip(
            *list(pairwise(left_polyedge))[i])

    # add new strip data
    new_skey = list(mesh.strips())[-1] + 1
    mesh.attributes['strips'][new_skey] = mesh.collect_strip(
        left_polyedge[0], right_polyedge[0])

    # update adjacent strips
    for i in range(len(polyedge)):
        old, left, right = polyedge[i], left_polyedge[i], right_polyedge[i]
        mesh.substitute_vertex_in_strips(old, left, left_strips)
        mesh.substitute_vertex_in_strips(old, right, right_strips)

    return new_skey, left_polyedge, right_polyedge


def add_element_start(mesh, u, v):
    pass


def add_element_main(mesh, u, v):
    pass


def add_element_end(mesh, u, v):
    pass


def update_strip_data(mesh, full_updated_polyedge, old_vkeys_to_new_vkeys):
    # orthogonal strips
    orth_to_update = {}
    # orth_skeys = []
    for old_u, old_v in pairwise(full_updated_polyedge):
        new_u = old_vkeys_to_new_vkeys[old_u][0]
        new_v = old_vkeys_to_new_vkeys[old_v][0]
        skey = mesh.edge_strip((old_u, old_v))
        edges = mesh.collect_strip(new_u, new_v)
        orth_to_update[skey] = edges
    for skey, edges in orth_to_update.items():
        mesh.attributes['strips'][skey] = edges

    # parallel strips
    paral_to_update = {}
    for skey, edges in mesh.attributes['strips'].items():
        if skey not in orth_to_update:
            for u, v in edges:
                if u in old_vkeys_to_new_vkeys:
                    u, v = v, u
                elif v in old_vkeys_to_new_vkeys:
                    u, v = u, v
                else:
                    continue
                new_v = [vkey for vkey in old_vkeys_to_new_vkeys[v]
                         if vkey in mesh.halfedge[u]][0]
                new_edges = mesh.collect_strip(u, new_v)
                paral_to_update[skey] = new_edges
                break
    for skey, edges in paral_to_update.items():
        mesh.attributes['strips'][skey] = edges

    # self strip
    n = max(mesh.attributes['strips']) + 1
    strip_edges = [tuple(old_vkeys_to_new_vkeys[vkey])
                   for vkey in full_updated_polyedge]
    mesh.attributes['strips'][n] = strip_edges

    return n


def update_polyedge(polyedge, old_vkey_to_new_vkey):

    return [old_vkey_to_new_vkey.get(vkey, vkey) for vkey in polyedge]


def sort_faces(mesh, u, v, w):

    sorted_faces = [[], []]
    k = 0
    vertex_faces = mesh.vertex_faces(v, ordered=True, include_none=True)
    f0 = mesh.halfedge[w][v] if w is not None else None
    i0 = vertex_faces.index(f0)
    vertex_faces = vertex_faces[i0:] + vertex_faces[:i0]

    for face in vertex_faces:

        if face is not None:
            sorted_faces[k].append(face)

        if u is None:
            if face is None:
                k = 1 - k
        elif face == mesh.halfedge[v][u]:
            k = 1 - k

    # indeterminate exception if u == w
    if u == w:
        return [face for faces in sorted_faces for face in faces]
    else:
        return sorted_faces


def adjacency_from_to_via_vertices(mesh, from_vkey, to_vkey, via_vkeys):
    # get mesh adjacency constraiend to from_vkey and via_keys, via_keys and via_keys, and via_vkeys and to_vkey

    all_vkeys = set([from_vkey, to_vkey] + via_vkeys)
    adjacency = {}
    for vkey, nbrs in mesh.adjacency.items():
        if vkey not in all_vkeys:
            continue
        else:
            sub_adj = {}
            for nbr, face in nbrs.items():
                if nbr not in all_vkeys or (vkey == from_vkey and nbr == to_vkey) or (vkey == to_vkey and nbr == from_vkey):
                    continue
                else:
                    sub_adj.update({nbr: face})
            adjacency.update({vkey: sub_adj})
    return adjacency


def polyedge_from_to_via_vertices(mesh, from_vkey, to_vkey, via_vkeys):
    # return shortest polyedge from_vkey to_vkey via_vkeys

    adjacency = adjacency_from_to_via_vertices(
        mesh, from_vkey, to_vkey, via_vkeys)
    return next(breadth_first_paths(adjacency, from_vkey, to_vkey))


def is_polyedge_valid_for_strip_addition(mesh, polyedge):
    if len(polyedge) > 2:
        if polyedge[0] == polyedge[-1] or (mesh.is_vertex_on_boundary(polyedge[0]) and mesh.is_vertex_on_boundary(polyedge[-1])):
            return True
    return False
