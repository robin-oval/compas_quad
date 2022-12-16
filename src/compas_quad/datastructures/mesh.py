from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from compas.datastructures import Mesh


__all__ = ['Mesh']


class Mesh(Mesh):

    def __init__(self):
        super(Mesh, self).__init__()

    def move(self, vector):
        x, y, z = vector
        for vkey, attr in self.vertices(data=True):
            attr['x'] += x
            attr['y'] += y
            attr['z'] += z
