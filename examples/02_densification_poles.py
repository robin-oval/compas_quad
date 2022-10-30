import os

from compas_quad.datastructures import CoarsePseudoQuadMesh

from compas_view2.app import App

HERE = os.path.dirname(__file__)
FILE = os.path.join(HERE, 'jsons/coarse_quad_mesh_british_museum_poles.json')

# read input data
coarse = CoarsePseudoQuadMesh.from_json(FILE)
box = coarse.bounding_box()
vector = [1.1 * (box[1][0] - box[0][0]), 0.0, 0.0]

# view coarse quad mesh
viewer = App(width=1600, height=900)
viewer.add(coarse)

# collect strip data
coarse.collect_strips()

# densification with target length
coarse.set_strips_density_target(t=.5)
coarse.densification()

# plot dense quad mesh
dense = coarse.dense_mesh()
dense.move(vector)
viewer.add(dense)
viewer.show()
