from __future__ import print_function
from __future__ import absolute_import
from __future__ import division


__all__ = [
    'is_adjacency_two_colorable'
]


def is_adjacency_two_colorable(adjacency):
    """Try to color a data of adjacency with two colors only withtout any element adjacent to each other having the same colour.

    Parameters
    ----------
    adjacency : dict
        Dictionary of adjacency between elements, each elements points to the list of adjacent elements.

    Returns
    -------
    key_to_color : dict, None
        A dictionary with vertex keys pointing to colors, if two-colorable.
        None if not two-colorable.

    """

    # store color status of network vertices (-1 means no color)
    key_to_color = {vkey: -1 for vkey in adjacency.keys()}

    # start from any vertex, color it and propagate to neighbors
    key_0 = next(iter(adjacency.keys()))
    sources = [key_0]

    count = len(adjacency.keys()) * 2

    # propagate until all vertices are colored or two adjacent vertices have the same color
    while count > 0 and sources:
        count -= 1
        key = sources.pop()
        nbr_colors = set([key_to_color[nbr] for nbr in adjacency[key]])

        # if two colors already exist in the neighbors, the network is not two-colourable
        if 0 in nbr_colors and 1 in nbr_colors:
            return None
        # otherwise, color with an available color
        else:
            if 0 not in nbr_colors:
                key_to_color[key] = 0
            elif 1 not in nbr_colors:
                key_to_color[key] = 1
            # add uncolored neighbors to sources
            sources += [nbr for nbr in adjacency[key]
                        if key_to_color[nbr] == -1 and nbr not in sources]

    return key_to_color
