from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from .addition import add_strip
from .deletion import delete_strip

from itertools import product
from random import random


__all__ = ['Lizard']


class Lizard:

    def __init__(self, quad_mesh):
        self.lizard = None
        self.grow = False
        self.mesh = quad_mesh

    def initiate(self, tail=None, head=None):
        if tail in self.mesh.vertex and head in self.mesh.vertex and tail in self.mesh.vertex[head]:
            self.lizard = [head, tail]
        else:
            tail = list(self.mesh.vertices())[0]
            head = self.mesh.vertex_neighbors(tail)[0]
            self.lizard = [head, tail]

    def turn(self):
        nbrs = self.mesh.vertex_neighbors(self.lizard[0], ordered=True)
        i = nbrs.index(self.lizard[1])
        new_head = nbrs[i + 1 - len(nbrs)]
        self.lizard.insert(0, new_head)
        if not self.grow:
            del self.lizard[-1]

    def pivot(self):
        nbrs = self.mesh.vertex_neighbors(self.lizard[1], ordered=True)
        i = nbrs.index(self.lizard[0])
        new_head = nbrs[i + 1 - len(nbrs)]
        self.lizard[0] = new_head

    def add(self):
        if self.grow:
            n, old_vkeys_to_new_vkeys = add_strip(self.mesh, self.lizard[1:])
            head = old_vkeys_to_new_vkeys[self.lizard[0]
                                          ] if self.lizard[0] in old_vkeys_to_new_vkeys else self.lizard[0]
            tail = old_vkeys_to_new_vkeys[self.lizard[1]]
            self.lizard = [head, tail]

        self.grow = not self.grow

    def delete(self):
        if self.grow:
            # use update polyedge
            print('Not Implemented')
            pass

        skey = self.mesh.edge_strip((self.lizard[0], self.lizard[1]))
        self.pivot()  # check that pivot ends along strip
        old_vkeys_to_new_vkeys = delete_strip(self.mesh, skey)
        # specify one when old_vkeys_to_new_vkeys[i] is a tuple
        self.lizard = [old_vkeys_to_new_vkeys[i] for i in self.lizard]

    def from_vector_to_string(self, vector):
        string = []
        for i in range(len(vector) // 2):
            x, y = vector[2 * i: 2 * i + 2]
            if (x, y) == (0, 0):
                string.append('t')
            elif (x, y) == (0, 1):
                string.append('p')
            elif (x, y) == (1, 0):
                string.append('a')
            elif (x, y) == (1, 1):
                string.append('d')
        return string

    def from_string_to_vector(self, string):
        vector = []
        for k in string:
            if k == 't':
                vector += [0, 0]
            elif k == 'p':
                vector += [0, 1]
            elif k == 'a':
                vector += [1, 0]
            elif k == 'd':
                vector += [1, 1]
        return vector

    def from_string_to_rules(self, string):
        for k in string:
            if k == 't':
                self.turn()
            elif k == 'p':
                self.pivot()
            elif k == 'a':
                self.add()
            elif k == 'd':
                self.delete()


def string_generation_brute(characters, length):
    for letters in product(characters, repeat=length):
        yield ''.join(letters)


def string_generation_random(characters, number, length, ratios=None):
    if not ratios:
        ratios = [1.0 / len(characters)] * len(characters)

    for _ in range(number):
        string = ''
        for _ in range(length):
            x = random()
            for i in range(len(characters)):
                if x < sum(ratios[:i + 1]):
                    string += characters[i]
                    break
        yield string


def string_generation_structured(characters, number, length):

    for _ in range(number):
        string = ''
        polyedge_length = 0  # number of 't' between odd and even pari of 'a'

        for i in range(length):

            # add strip before end of string if polyedge is being collected
            if i == length - 1 and polyedge_length != 0:
                string += 'a'
                polyedge_length = 0
                continue

            x = random()

            # NO polyedge being collected - uniform chances
            if polyedge_length == 0:
                if x < 0.33:
                    string += 't'
                elif x < 0.67:
                    string += 'p'
                else:
                    string += 'a'
                    polyedge_length += 1

            # polyedge being collected - more chances to obtain 't' or 'p'
            
            # do not add if polyedge has only one vertex
            elif polyedge_length == 1:
                if x < 0.5:
                    string += 't'
                    polyedge_length += 1
                else:
                    string += 'p'

            else:
                if x < 0.4:
                    string += 't'
                    polyedge_length += 1
                elif x < 0.8:
                    string += 'p'
                else:
                    string += 'a'
                    polyedge_length = 0

        yield string
