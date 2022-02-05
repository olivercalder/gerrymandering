# Make seed map

from distutils.ccompiler import gen_preprocess_options
import maup
import geopandas
from gerrychain import Graph, Partition
import matplotlib.pyplot as plt

for orig, proposed, outfile in [
        ('data/mn_mggg/MN12_18/mn_precincts12_18.shp', 'zip://data/MN_proposed/c2103_0-shp.zip', "mn_proposed_map_C2103.json"),
        ('data/mn_mggg/MN12_18/mn_precincts12_18.shp', 'zip://data/MN_proposed/c2104_0-shp.zip', "mn_proposed_map_C2104.json"),
        ('data/mn_mggg/MN12_18/mn_precincts12_18.shp', 'zip://data/MN_proposed/c2105_0-shp.zip', "mn_proposed_map_C2105.json"),
        ('data/mn_mggg/MN12_18/mn_precincts12_18.shp', 'zip://data/MN_proposed/c2106_0-shp.zip', "mn_proposed_map_C2106.json"),
        ]:

    mggg_shapefile_path = orig
    proposed_shapefile_path = proposed

    districts = geopandas.read_file(proposed_shapefile_path)
    units = geopandas.read_file(mggg_shapefile_path)
    assignment = maup.assign(units, districts)

    units["District"] = assignment
    print(units)

    graph = Graph.from_file(mggg_shapefile_path, ignore_errors=True)
    print("Generated", outfile)

    # graph.join(units, columns=["District"], left_index="geometry", right_index="geometry")
    graph.join(units, columns=["District"])
    graph.to_json(outfile)

    real_life_plan = Partition(graph, "District")
    real_life_plan.plot(units, figsize=(30, 30))
    plt.axis('off')
    plt.savefig(outfile.replace('.json', '.png'))
    plt.close('all')
