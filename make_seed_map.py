# Make seed map

from distutils.ccompiler import gen_preprocess_options
import maup 
import geopandas 
from gerrychain import Graph, Partition 
import matplotlib.pyplot as plt

mggg_shapefile_path = '../data/tx_mggg/TX_vtds/TX_vtds.shp'
proposed_shapefile_path = 'zip://../data/TX_proposed/tx_cong_2021.zip'

# ex = geopandas.read_file('../Ex_data/tl_2012_25_sldu.shp')
# print("\nex:")
# print(ex)

districts = geopandas.read_file(proposed_shapefile_path)
units = geopandas.read_file(mggg_shapefile_path)
# print(units)
# print(districts)
assignment = maup.assign(units, districts)
# print(assignment)

units["District"] = assignment 
print(units) 

graph = Graph.from_file(mggg_shapefile_path, ignore_errors=True)
print("made graph")

# graph.join(units, columns=["District"], left_index="geometry", right_index="geometry")
graph.join(units, columns=["District"])
graph.to_json("tx_seed_graph.json")

real_life_plan = Partition(graph, "District")
real_life_plan.plot(units, figsize=(10, 10), cmap="tab20")
plt.axis('off')
plt.show()