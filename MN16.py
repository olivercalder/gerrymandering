import matplotlib.pyplot as plt
from gerrychain import (GeographicPartition, Partition, Graph, MarkovChain,
                        proposals, updaters, constraints, accept, Election, metrics)
from gerrychain.proposals import recom
from functools import partial
import pandas
import csv
import random
random.seed(1123)
from math import exp

partition_counter = 0


def compute_population_score(partition):
    '''
    Returns a real-valued score >= 0 related to population deviation.
    Lower is better.

    From Mattingly supplemental material.
    '''
    ideal_population = sum(partition.population.values()) / len(partition.population)
    population_deviation_sq = 0.0
    for district_pop in partition.population.values():
        population_deviation_sq += ((district_pop / ideal_population) - 1)**2
    population_deviation = population_deviation_sq**0.5
    return population_deviation


def compute_compactness_score(partition):
    '''
    Returns a real-valued score >= 0 related to isoparametric compactness.
    Lower is better.

    Similar to Polsby-Popper/Schwartzberg
    '''
    compactness_sum = 0.0
    for district in partition.perimeter:    # equivalently, in partition.area:
        compactness_sum += partition.perimeter[district]**2 / partition.area[district]
    return compactness_sum


def get_districts_in_county_sorted_by_dominance(county_info, partition):
    # Sort districts in given county by their dominance
    # Dominant = highest proportion of district's precincts in given county
    # county_info.contains{4, 8, 7} <- districts in county
    # Returns {8, 4, 7} <- same districts, sorted most-to-least dominant
    return sorted(county_info.contains, reverse=True, key=lambda district:
            len([precinct for precinct in partition.parts[district] if precinct in county_info.nodes])
            / len(partition.parts[district]))


def get_f_scores(county_info, partition):
    '''
    Returns a list of the f_n scores (n must not exceed the number of districts).
    ie. f_1 is the proportion of precincts in the county which fall within the
    in the dominant district.
    f_k is the proportion of precincts in the county which fall within the first
    k most dominant districts.

    f1 is at index 0 of output, fk is at index k-1 of output, etc.
    '''
    if county_info.split == 0:  #county is NOT_SPLIT
        return 1
    # Do we care more about newly-split counties than counties which were previously split?

    # Determine dominant district: with highest percentage of precincts in the given county
    # district with highest precinct percentage in the given county
    dominant_districts = get_districts_in_county_sorted_by_dominance(county_info, partition)

    # precinct_counts_per_district = [0, 0, 0] <- 3=#districts in county
    # Becomes precinct_counts_per_district = [54, 38, 13]
    #   <- number of precincts in those ditsricts (only those inside the county)
    precinct_counts_per_district = [0 for i in range(len(dominant_districts))]
    for precinct in county_info.nodes:
        for k in range(len(dominant_districts)):
            # k is the index of the district according to its dominance
            district = dominant_districts[k]
            # Count how many precincts are in the kth most dominant district
            if precinct in partition.parts[district]:
                precinct_counts_per_district[k] += 1

    f_scores = [0 for i in range(len(dominant_districts))]
    total_precincts_in_county = len(county_info.nodes)
    # first f score = # precincts in most dominant district / total precincts in county
    f_scores[0] = precinct_counts_per_district[0] / total_precincts_in_county
    for k in range(1, len(precinct_counts_per_district)):
        # Change to be count of precincts in the first k most dominant districts
        precinct_counts_per_district[k] += precinct_counts_per_district[k - 1]
        f_scores[k] = precinct_counts_per_district[k] / total_precincts_in_county
        # Number of precincts in k most dominant districts / total number of precincts in county
    return f_scores


def compute_county_split_score(partition):
    '''
    Returns a real-valued score >= 0 related to county splitting.
    Lower values mean fewer counties split, and splits are less dramatic.
    '''
    # TODO: make sure f and W scores aren't off by 1
    # hopefully the number of counties which are split is >= the maximum number
    # of districts which split any given county
    county_split_counts = [0 for i in range(len(partition.county_splits))] #[0, 0, 0, ..., 0] <- a 0 for every county
    w_scores = [0.0 for i in range(len(partition.county_splits))]
    # w_scores[k] will hold the w_{2+k} score, where w_{2+k} is computed as:
    # w_{2+k} = sum([(1 - f_{1+k}(county))**0.5 for county in counties split at least 2+k ways])
    max_splits = 0
    # loop over counties
    for county_info in partition.county_splits.values():
        splits = len(county_info.contains)
        max_splits = max(splits, max_splits)
        f_scores = get_f_scores(county_info, partition)
        # loop over number of county splits
        for k in range(splits):
            # if a county is split k=4 ways, add 1 to county_split_counts 1, 2, 3, and 4
            county_split_counts[k] += 1
            # add county component to W_k score
            w_scores[k] += (1 - f_scores[k])**0.5

    mc = 100    # TODO Change this!!!!!!!!
    # mc^k is the coefficient for the w_{2+k} score
    # or use mc*k
    j_c = 0.0   # County split score
    for k in range(max_splits):
        j_c += mc**k * county_split_counts[k] * w_scores[k]
    return j_c


def vra_score(partition):
    pass


def score_function(partition):
    pop_score_weight = 1
    comp_score_weight = 1   # TODO Change these weights
    county_score_weight = 1
    partition_score = 0.0
    partition_score += pop_score_weight * compute_population_score(partition)
    partition_score += comp_score_weight * compute_compactness_score(partition)
    partition_score += county_score_weight * compute_county_split_score(partition)
    return partition_score


def Q_func(partition1, partition2):
    # cut edges are all the edges from one part of the partition to another
    conflicted = len(partition1.cut_edges)
    part1 = partition1.cut_edges
    part2 = partition2.cut_edges
    # cross edges are edges that were cut in the 1st but are not in the 2nd
    cross = len(part2-part1)
    return (cross / conflicted) / 2


def acceptance_function(partition):
    # beta = 0 for first 10,000 accepted steps
    #   From 10,000 to 70,000, beta grows linearly to 1, only growing on accepted steps
    #   From 70,000 and up, beta = 1
    global partition_counter
    if partition_counter < 10000:
        beta = 0
    elif partition_counter < 70000:
        beta = (partition_counter - 10000) / 60000
    else:
        beta = 1

    Q2 = Q_func(partition.parent, partition)
    if Q2 == 0:
        return False
    Q1 = Q_func(partition, partition.parent)

    p = min(1, (Q1 / Q2) * exp(-beta * (score_function(partition) - score_function(partition.parent))))

    # accept the partition with probability p
    if random.random() < p:
        partition_counter += 1
        return True
    else:
        return False


def get_chain():
    """
    Get data, build graph, make updaters and constraints, run chain
    :return: chain
    """

    SHAPEFILE_PATH = "../data/mn_mggg/MN16/mn_precincts16.shp"
    POPULATION_COL = "TOTPOP"
    VAP_COL = 'VAP'
    BVAP_COL = 'BVAP'
    HVAP_COL = 'HVAP'
    ASSIGNMENT_COL = "CONGDIST"

    graph = Graph.from_file(SHAPEFILE_PATH)

    # Make updaters
    all_updaters = {
            'population': updaters.Tally(POPULATION_COL, alias='population'),
            'vap': updaters.Tally(VAP_COL, alias='vap'),
            'bvap': updaters.Tally(BVAP_COL, alias='bvap'),
            'hvap': updaters.Tally(HVAP_COL, alias='hvap'),
            'county_splits': updaters.county_splits('county_splits', 'COUNTYNAME'),
            #"county_splits": updaters.CountySplit('county_splits', 'COUNTYNAME'),  # TODO check on this one
            }

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

    partition_counter = 0

    chain = MarkovChain(
        proposal=proposal,
        constraints=[
            pop_constraint,
            compactness_bound
        ],
        accept=acceptance_function,
        # accept=accept.always_accept,
        initial_state=initial_partition,
        total_steps=100
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
    #print(partition.area)
    for partition in chain:
        print(type(partition.cut_edges))
        print(partition.cut_edges)
        print(partition.graph.edges)
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
    data = get_chain_data(chain)
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
