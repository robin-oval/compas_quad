import os

from compas_quad.datastructures import CoarseQuadMesh

from compas_plotters.plotter import Plotter

from compas_view2.app import App

HERE = os.path.dirname(__file__)
FILE = os.path.join(HERE, 'jsons/coarse_quad_mesh_british_museum.json')

# read input data
coarse_quad_mesh = CoarseQuadMesh.from_json(FILE)

# plot coarse quad mesh
# plotter = Plotter(figsize=(5, 5))
# plotter.vertex_color((1, 1, 1))
# plotter.add(coarse_quad_mesh)
# plotter.show()

viewer = App(width=1600, height=900)
viewer.add(coarse_quad_mesh)
viewer.show()

# # collect strip data
# coarse_quad_mesh.collect_strips()

# # densification with uniform density
# coarse_quad_mesh.set_strips_density(3)
# coarse_quad_mesh.densification()

# # plot dense quad mesh
# plotter = Plotter(figsize=(5, 5))
# plotter.add(coarse_quad_mesh.get_quad_mesh())
# plotter.show()

# # densification with target length
# coarse_quad_mesh.set_strips_density_target(t=.5)
# coarse_quad_mesh.densification()

# # plot dense quad mesh
# plotter = Plotter(figsize=(5, 5))
# plotter.add(coarse_quad_mesh.get_quad_mesh())
# plotter.show()

# # change density of one strip
# skey = list(coarse_quad_mesh.strips())[0]
# coarse_quad_mesh.set_strip_density(skey, 10)
# coarse_quad_mesh.densification()

# # plot dense quad mesh
# plotter = Plotter(figsize=(5, 5))
# plotter.add(coarse_quad_mesh.get_quad_mesh())
# plotter.show()
