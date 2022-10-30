import os

from compas_quad.datastructures import QuadMesh, CoarseQuadMesh

from compas.geometry import scale_vector

from compas_view2.app import App

HERE = os.path.dirname(__file__)
FILE = os.path.join(HERE, 'jsons/dense_quad_mesh_british_museum.json')

# read input data
dense = QuadMesh.from_json(FILE)
coarse = CoarseQuadMesh.from_quad_mesh(dense)
box = coarse.bounding_box()
vector = [1.1 * (box[1][0] - box[0][0]), 0.0, 0.0]

# view coarse quad mesh
viewer = App(width=1600, height=900)
viewer.add(dense)

# densification with target length
coarse.set_strips_density_target(t=.33)
coarse.densification()
redense = coarse.dense_mesh()

# plot meshes
coarse.move(vector)
viewer.add(coarse)
redense.move(scale_vector(vector, 2))
viewer.add(redense)
viewer.show()
