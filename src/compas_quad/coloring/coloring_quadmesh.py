from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from compas_quad.datastructures import CoarseQuadMesh

from compas.topology import adjacency_from_edges, vertex_coloring

from compas_quad.coloring import is_adjacency_two_colorable

from compas.geometry import centroid_points

from compas.utilities import pairwise


__all__ = [
    'quad_mesh_strip_2_coloring',
    'quad_mesh_strip_n_coloring',
    'quad_mesh_polyedge_2_coloring',
    'quad_mesh_polyedge_n_coloring',
    'quad_mesh_polyedge_subcolor',
    'dense_quad_mesh_polyedge_2_coloring'
]


def quad_mesh_strip_2_coloring(quad_mesh):
    """Try to color the strips of a quad mesh with two colors only without overlapping strips with the same color.

    Parameters
    ----------
    quad_mesh : QuadMesh
        A quad mesh.

    Returns
    -------
    dict, None
        A dictionary with strip keys pointing to colors, if two-colorable.
        None if not two-colorable.
    """

    vertices, edges = quad_mesh.strip_graph()
    return is_adjacency_two_colorable(adjacency_from_edges(edges))


def quad_mesh_strip_n_coloring(quad_mesh):
    """Color the strips of a quad mesh with a minimum number of colors without overlapping strips with the same color.

    Parameters
    ----------
    quad_mesh : QuadMesh
        A quad mesh.

    Returns
    -------
    dict
        A dictionary with strip keys pointing to colors.
    """

    vertices, edges = quad_mesh.strip_graph()
    return vertex_coloring(adjacency_from_edges(edges))


def quad_mesh_polyedge_2_coloring(quad_mesh, skip_singularities=True, edge_output=False):
    """Try to color the polyedges of a quad mesh with two colors only without overlapping polyedges with the same color.
    Polyedges connected by their extremities, which are singularities, do not count as overlapping.

    Parameters
    ----------
    quad_mesh : QuadMesh
        A quad mesh.
    skip_singularities : bool
            Boolean whether to discard (True) or consider (False) the connections between polyedges at singularities.
    edge_output : bool, optional
        Optional to output coloring per edge keys instead of polyedge key.

    Returns
    -------
    dict, None
        A dictionary with polyedge keys pointing to colors, if two-colorable. If edge_output, edge keys pointing to colors.
        None if not two-colorable.
    """

    vertices, edges = quad_mesh.polyedge_graph(skip_singularities)
    polyedge_coloring = is_adjacency_two_colorable(adjacency_from_edges(edges))

    if not polyedge_coloring or not edge_output:
        return polyedge_coloring

    else:
        edge_coloring = {}
        for pkey, group in polyedge_coloring.items():
            for u, v in pairwise(quad_mesh.attributes['polyedges'][pkey]):
                edge_coloring.update({(u, v): group, (v, u): group})
        return edge_coloring


def quad_mesh_polyedge_n_coloring(quad_mesh, edge_output=False):
    """Color the polyedges of a quad mesh with a minimum number of colors without overlapping polyedges with the same color.
    Polyedges connected by their extremities, which are singularities, do not count as overlapping.

    Parameters
    ----------
    quad_mesh : QuadMesh
        A quad mesh.
    edge_output : bool, optional
        Optional to output coloring per edge keys instead of polyedge key.

    Returns
    -------
    dict
        A dictionary with polyedge keys pointing to colors. If edge_output, edge keys pointing to colors.
    """

    vertices, edges = quad_mesh.polyedge_graph()
    polyedge_coloring = vertex_coloring(adjacency_from_edges(edges))
    if not polyedge_coloring or not edge_output:
        return polyedge_coloring
    else:
        edge_coloring = {}
        for pkey, group in polyedge_coloring.items():
            for u, v in pairwise(quad_mesh.attributes['polyedges'][pkey]):
                edge_coloring.update({(u, v): group, (v, u): group})
        return edge_coloring


def quad_mesh_polyedge_subcolor(quad_mesh, color=0):
    """Color the subset of quad-mesh polyedges that are from the same colour based on their adjacency through edges.

    Parameters
    ----------
    quad_mesh : QuadMesh
        A quad mesh.
    color : int, optional
        Select the group color to subcolor.

    Returns
    -------
    dict
        A dictionary with subset polyedge keys pointing to colors.
    """

    polyedge2color = quad_mesh_polyedge_n_coloring(quad_mesh)

    subnodes = []
    for pkey, polyedge in quad_mesh.polyedges(data=True):
        if polyedge2color[pkey] == 0:
            subnodes.append(centroid_points(polyedge))
            continue

    edge2polyedge = {}
    for pkey, polyedge in quad_mesh.polyedges(data=True):
        for u, v in pairwise(polyedge):
            edge2polyedge[(u, v)] = pkey
            edge2polyedge[(v, u)] = pkey

    subedges = set()
    for fkey in quad_mesh.faces():
        u, v, w, x = list(quad_mesh.face_vertices(fkey))
        if polyedge2color[edge2polyedge[(u, v)]] == color:
            pkey0, pkey1 = edge2polyedge[(u, v)], edge2polyedge[(w, x)]
            if (pkey0, pkey1) not in subedges and (pkey1, pkey0) not in subedges:
                subedges.add((pkey0, pkey1))
        elif polyedge2color[edge2polyedge[(v, w)]] == color:
            pkey0, pkey1 = edge2polyedge[(v, w)], edge2polyedge[(x, u)]
            if (pkey0, pkey1) not in subedges and (pkey1, pkey0) not in subedges:
                subedges.add((pkey0, pkey1))

    return vertex_coloring(adjacency_from_edges(subedges))


def dense_quad_mesh_polyedge_2_coloring(quad_mesh):
    # assume that strips and polyedges are collected in the quad mesh

    # get coarse quad mesh
    coarse_quad_mesh = CoarseQuadMesh.from_quad_mesh(quad_mesh)

    # get coarse strip color
    coarse_skey_to_color = quad_mesh_strip_2_coloring(coarse_quad_mesh)

    # get coarse edge color
    coarse_edge_to_color = {edge: coarse_skey_to_color[skey] for skey in coarse_quad_mesh.strips(
    ) for edge in coarse_quad_mesh.strip_edges(skey)}

    # get dense polyedge color
    dense_polyedge_to_color = {tuple(coarse_quad_mesh.attributes['edge_coarse_to_dense'][u][v]): color for (
        u, v), color in coarse_edge_to_color.items()}

    # get some dense edge color
    some_dense_edge_to_color = {edge: color for polyedge,
                                color in dense_polyedge_to_color.items() for edge in pairwise(polyedge)}

    # get strip color
    dense_strip_to_color = {}
    for skey in quad_mesh.strips():
        for u, v in quad_mesh.strip_edges(skey):
            color = some_dense_edge_to_color.get(
                (u, v), some_dense_edge_to_color.get((v, u), None))
            if color is not None:
                dense_strip_to_color[skey] = color
                break

    # get edge color
    all_dense_edge_to_color = {edge: color for skey, color in dense_strip_to_color.items(
    ) for edge in quad_mesh.strip_edges(skey)}

    # get polyedge color
    dense_polyedge_to_color = {}
    for pkey, polyedge in quad_mesh.polyedges(data=True):
        u, v = polyedge[:2]
        color = all_dense_edge_to_color.get(
            (u, v), all_dense_edge_to_color.get((v, u), None))
        dense_polyedge_to_color[pkey] = color

    return dense_polyedge_to_color
