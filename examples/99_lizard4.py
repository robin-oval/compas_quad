import os

from time import time

from math import pi, cos, sin

from compas.datastructures import Mesh

from compas_quad.datastructures import CoarsePseudoQuadMesh

from compas_quad.grammar.addition2 import add_strip_lizard, add_strip_lizard_2, lizard_atp

from compas_quad.grammar.lizard import string_generation_brute, string_generation_random, string_generation_structured, string_generation_evolution

from compas_fd.solvers import fd_numpy

from compas_viewer.viewer import Viewer

# custom postprocessing function
def postprocessing(mesh):

    key2index = {vkey: i for i, vkey in enumerate(mesh.vertices())}
    index2key = {i: vkey for i, vkey in enumerate(mesh.vertices())}

    # map mesh boundary vertices to circle
    fixed = mesh.vertices_on_boundary()[:-1]
    n = len(fixed)
    for i, vkey in enumerate(fixed):
        attr = mesh.vertex[vkey]
        attr['x'] = 0.5 * cos(i / n * 2 * pi)
        attr['y'] = 0.5 * sin(i / n * 2 * pi)
        attr['z'] = 0

    # form finding with force density method
    vertices = [mesh.vertex_coordinates(vkey) for vkey in mesh.vertices()]
    edges = [(key2index[u], key2index[v]) for u, v in mesh.edges()]
    fixed = [key2index[vkey] for vkey in fixed]
    q = [1.0] * len(edges)
    loads = [[0.0, 0.0, 1.0 / len(vertices)]] * len(vertices)
    result = fd_numpy(vertices=vertices, edges=edges, fixed=fixed, forcedensities=q, loads=loads)
    xyz = result.vertices
    for i, (x, y, z) in enumerate(xyz):
        vkey = index2key[i]
        attr = mesh.vertex[vkey]
        attr['x'] = x
        attr['y'] = y
        attr['z'] = z

### parameters ###

input_mesh_refinement = 1  # densify the input 1-face quad mesh
output_mesh_refinement = 3  # densify the ouput quad mesh

postprocess = True
densify = True
array = False
view = True
export_json = False

### intialise ###

# dummy mesh with a single quad face
vertices = [[0.5, 0.5, 0.0], [-0.5, 0.5, 0.0], [-0.5, -0.5, 0.0], [0.5, -0.5, 0.0]]
faces = [[0, 1, 2, 3]]
coarse = CoarsePseudoQuadMesh.from_vertices_and_faces(vertices, faces)

# denser mesh
coarse.collect_strips()
coarse.strips_density(input_mesh_refinement)
coarse.densification()
mesh0 = coarse.dense_mesh()
mesh0.collect_strips()

if view:
    viewer = Viewer()
    viewer.scene.add(mesh0)

### lizard - let's grooow! ###

# position lizard
for vkey in mesh0.vertices_on_boundary():
    if mesh0.vertex_degree(vkey) == 2:
        body = vkey
        tail, head = [nbr for nbr in mesh0.vertex_neighbors(vkey) if mesh0.is_vertex_on_boundary(nbr)]
    break
lizard = (tail, body, head)
print('lizard initial position', lizard)

# produce strings
# strings = ['t', 'tt', 'ttt', 'tttt']
# strings = ['ata', 'atta', 'attta', 'atttta']
strings = ['attatatattatta']
# strings = ['attatpatatptatta']

# apply
t0 = time()
mesh2string = {}
for k, string in enumerate(strings):
    print(string)

    # modifiy topology
    mesh = CoarsePseudoQuadMesh.from_vertices_and_faces(*mesh0.to_vertices_and_faces())
    # tail, body, head = add_strip_lizard_2(mesh, lizard, string)
    tail, body, head = lizard_atp(mesh, lizard, string)
    
    poles = []
    for fkey in mesh.faces():
        fv = mesh.face_vertices(fkey)
        if len(fv) == 3:
            if body in fv:
                poles.append(mesh.vertex_coordinates(body))
            else:
                # warn if pole missing
                'pbm identification pole'
                poles.append(mesh.vertex_coordinates(fv[0]))
    
    # warn if mesh not manifold
    if not mesh.is_manifold():
        print('mesh not manifold')
        continue

    # export JSON
    if export_json:
        HERE = os.path.dirname(__file__)
        FILE = os.path.join(HERE, 'data/{}_{}.json'.format(input_mesh_refinement, string))
        mesh_json = Mesh.from_vertices_and_faces(*mesh.to_vertices_and_faces())
        mesh_json.to_json('/Users/roval/Desktop/mesh4.json')

    # geometry and density processing
    if postprocess:
        postprocessing(mesh)
        if densify:
            mesh = CoarsePseudoQuadMesh.from_vertices_and_faces_with_poles(*mesh.to_vertices_and_faces(), poles=poles)
            mesh.collect_strips()
            mesh.strips_density(output_mesh_refinement)
            mesh.densification()
            mesh = mesh.dense_mesh()
            postprocessing(mesh)

    mesh2string[mesh] = string

# results
t1 = time()
print('computation time {}s for {} meshes'.format(round(t1 - t0, 3), len(mesh2string)))

if array:
    n = len(mesh2string)
    for k, mesh in enumerate(mesh2string):
        n2 = int(n ** 0.5)
        i = int(k / n2)
        j = int(k % n2)
        mesh.move([1.5 * (i + 1), 1.5 * (j + 1), 0.0])
else:
    for k, mesh in enumerate(mesh2string):
        mesh.move([1.5 * (k + 1), 0.0, 0.0])

if view:
    for mesh in mesh2string:
        viewer.scene.add(mesh)
    viewer.show()
