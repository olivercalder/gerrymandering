import matplotlib.pyplot as plt
from gerrychain import (GeographicPartition, Partition, Graph, MarkovChain,
                        proposals, updaters, constraints, accept, Election, metrics)
from gerrychain.proposals import recom
from functools import partial
import pandas
import csv
import random
random.seed(1123)


def get_chain():
    """
    Get data, build graph, make updaters and constraints, run chain
    :return: chain
    """

    SHAPEFILE_PATH = "data/mn_mggg/MN16/mn_precincts16.shp"
    POPULATION_COL = "TOTPOP"
    ASSIGNMENT_COL = "CONGDIST"

    graph = Graph.from_file(SHAPEFILE_PATH)

    # Make updaters
    all_updaters = {"population": updaters.Tally(POPULATION_COL, alias="population")}
    elections = [
            Election('PRES16',  {'Democratic': 'PRES16D',   'Republican': 'PRES16R'}, alias='PRES16'),
            Election('USH16',   {'Democratic':  'USH16D',   'Republican':  'USH16R'}, alias='USH16'),
            Election('SSEN16',  {'Democratic': 'SSEN16D',   'Republican': 'SSEN16R'}, alias='SSEN16'),
            Election('SH16',    {'Democratic':   'SH16D',   'Republican':   'SH16R'}, alias='SH16'),
            ]
    all_updaters.update({election.name: election for election in elections})

    initial_partition = GeographicPartition(graph, assignment=ASSIGNMENT_COL, updaters=all_updaters)

    ideal_population = sum(initial_partition["population"].values()) / len(initial_partition)

    # We use functools.partial to bind the extra parameters (pop_col, pop_target, epsilon, node_repeats)
    # of the recom proposal.
    proposal = partial(recom,
                    pop_col=POPULATION_COL,
                    pop_target=ideal_population,
                    epsilon=0.02,
                    node_repeats=2
                    )

    pop_constraint = constraints.within_percent_of_ideal_population(initial_partition, 0.02)

    compactness_bound = constraints.UpperBound(
        lambda p: len(p["cut_edges"]),
        2*len(initial_partition["cut_edges"])
    )


    def acc_fn(partition):
        return random.uniform(0, 1)


    chain = MarkovChain(
        proposal=proposal,
        constraints=[
            pop_constraint,
            compactness_bound
        ],
        accept=accept.always_accept,
        initial_state=initial_partition,
        total_steps=1000
    )

    return chain


def get_chain_data(chain):
    return pandas.DataFrame(
            sorted(partition["USH16"].percents("Democratic"))
            for partition in chain.with_progress_bar())


def display_chain_data(data):
    fig, ax = plt.subplots(figsize=(8, 6))

    # Draw 50% line
    ax.axhline(0.5, color="#cccccc")

    # Draw boxplot
    data.boxplot(ax=ax, positions=range(len(data.columns)))

    # Draw initial plan's Democratic vote %s (.iloc[0] gives the first row)
    plt.plot(data.iloc[0], "ro")

    # Annotate
    ax.set_title("Comparing the 2011 plan to an ensemble")
    ax.set_ylabel("Democratic vote % (US House 2016)")
    ax.set_xlabel("Sorted districts")
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])

    plt.show()


def explore_chain(chain):
    print(partition.area)
    # print(partition.perimeter)
    # print(partition.population)
    # print(partition["SSEN16"].wins("Democratic"))
    # print(partition["SSEN16"].counts("Democratic"))



def out_csv(chain):
    with open("MN_data.csv", "w", newline="") as csvfile:
        f = csv.writer(csvfile)
        headings = ["State", "efficiency_gap", "partisan_bias", "d_seats", "r_seats"]
        num_dists = 8
        # add vote count headings
        for i in range(num_dists):
            headings.append("d_count_" + str(i))
            headings.append("r_count_" + str(i))
        f.writerow(headings)

        state = 1
        for partition in chain:

            # get election metrics, seats won and add to row
            e_g = partition["SSEN16"].efficiency_gap()
            p_b = partition["SSEN16"].partisan_bias()
            dem = partition["SSEN16"].wins("Democratic")
            rep = partition["SSEN16"].wins("Republican")
            row = [state, e_g, p_b, dem, rep]

            # get lists of vote counts and add to row
            dem_list = partition["SSEN16"].counts("Democratic")
            rep_list = partition["SSEN16"].counts("Republican")
            assert num_dists == len(dem_list)
            for i in range(num_dists):
                row.append(dem_list[i])
                row.append(rep_list[i])

            f.writerow(row)
            state += 1



if __name__ == "__main__":
    chain = get_chain()
    # data = get_chain_data(chain)
    # display_chain_data(data)

    # explore_chain(chain)

    out_csv(chain)


"""
Output we can get:
    Various election metrics
        Vote counts
        Mean-median
        Partisan bias
        etc.
    Partition characteristics
        Area
        Perimeter
        Cut edges
        Interior and exterior boundaries (?)
        All the node variables (see data readme):
            NH_WHITE
            NH_BLACK
            etc
"""
