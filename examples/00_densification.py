import os

from compas_quad.datastructures import Mesh, CoarseQuadMesh

from compas.geometry import scale_vector

from compas_viewer.viewer import Viewer

HERE = os.path.dirname(__file__)
FILE = os.path.join(HERE, 'jsons/coarse_quad_mesh_british_museum.json')

# read input data
# mesh = Mesh.from_json(FILE)
# coarse = Mesh.from_vertices_and_faces(*mesh.to_vertices_and_faces())
# coarse = CoarseQuadMesh.from_json(FILE)
vertices = [[0,0,0],[1,0,0],[1,1,0],[0,1,0]]
faces = [[0,1,2,3]]
coarse = CoarseQuadMesh.from_vertices_and_faces(vertices, faces)
# box = coarse.bounding_box()
# vector = [1.1 * (box[1][0] - box[0][0]), 0.0, 0.0]
vector = [2,0,0]

# view coarse quad mesh
# viewer = App(width=1600, height=900)
viewer = Viewer()
viewer.scene.add(coarse)

# collect strip data
coarse.collect_strips()

# densification with uniform density
coarse.strips_density(3)
coarse.densification()

# plot dense quad mesh
dense = coarse.dense_mesh()
dense.move(vector)
viewer.scene.add(dense)

# densification with target length
coarse.set_strips_density_target(t=.5)
coarse.densification()

# plot dense quad mesh
dense = coarse.dense_mesh()
dense.move(scale_vector(vector, 2))
viewer.scene.add(dense)

# change density of one strip
skey = list(coarse.strips())[0]
coarse.strip_density(skey, 10)
coarse.densification()

# plot dense quad mesh
dense = coarse.dense_mesh()
dense.move(scale_vector(vector, 3))
viewer.scene.add(dense)
viewer.show()
