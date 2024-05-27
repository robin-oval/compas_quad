from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from compas.datastructures import Mesh


__all__ = ['Mesh',
           'mesh_substitute_vertex_in_faces'
           ]


class Mesh(Mesh):

    def __init__(self):
        super(Mesh, self).__init__()

    def move(self, vector):
        x, y, z = vector
        for vkey, attr in self.vertices(data=True):
            attr['x'] += x
            attr['y'] += y
            attr['z'] += z


def mesh_substitute_vertex_in_faces(mesh, old_vkey, new_vkey, fkeys=None):
    """Substitute in a mesh a vertex by another one.
    In all faces by default or in a given set of faces.

    Parameters
    ----------
    old_vkey : int
        The old vertex key.
    new_vkey : int
        The new vertex key.
    fkeys : list[int], optional
        List of face keys where to subsitute the old vertex by the new one.
        Default is to subsitute in all faces.

    Returns
    -------
    list[int]
        The list of modified faces.

    """

    # apply to all faces if there is none chosen
    if fkeys is None:
        fkeys = list(mesh.faces())

    # substitute vertices
    for fkey in fkeys:
        face_vertices = [new_vkey if key == old_vkey else key for key in mesh.face_vertices(fkey)]
        mesh.delete_face(fkey)
        mesh.add_face(face_vertices, fkey)

    return fkeys