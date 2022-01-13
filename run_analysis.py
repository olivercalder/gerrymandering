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
from networkx import is_connected, connected_components

accepted_partition_counter = 0


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
        return [1]
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

    # mc = COUNTY_SPLIT_COEFFICIENT
    # mc^k is the coefficient for the w_{2+k} score
    # or use mc*k
    j_c = 0.0   # County split score
    for k in range(max_splits):
        j_c += COUNTY_SPLIT_COEFFICIENT**k * county_split_counts[k] * w_scores[k]
    return j_c


def compute_vra_score(partition):
    black_opportunity_districts = 0
    hisp_opportunity_districts = 0
    for key in partition.vap:
        if partition.bvap[key] / partition.vap[key] > OPPORTUNITY_THRESHOLD:
            black_opportunity_districts += 1
        if partition.hvap[key] / partition.vap[key] > OPPORTUNITY_THRESHOLD:
            hisp_opportunity_districts += 1
    return BLACK_OPP_WEIGHT * (black_opportunity_districts - BLACK_OPP_TARGET) ** 2 \
            + HISPANIC_OPP_WEIGHT * (hisp_opportunity_districts - HISPANIC_OPP_TARGET) ** 2


def score_function(partition):
    partition_score = 0.0
    partition_score += POPULATION_SCORE_WEIGHT * compute_population_score(partition)
    partition_score += COMPACTNESS_SCORE_WEIGHT * compute_compactness_score(partition)
    partition_score += COUNTY_SCORE_WEIGHT * compute_county_split_score(partition)
    partition_score += VRA_SCORE_WEIGHT * compute_vra_score(partition)
    return partition_score


def Q_func(partition1, partition2):
    # cut edges are all the edges from one part of the partition to another
    conflicted = len(partition1.cut_edges)
    # cross edges are edges that were cut in the 1st but are not in the 2nd
    cross = len(partition1.cut_edges - partition2.cut_edges)
    return (cross / conflicted) / 2


def acceptance_function(partition):
    # beta = 0 for first 10,000 accepted steps
    #   From 10,000 to 70,000, beta grows linearly to 1, only growing on accepted steps
    #   From 70,000 and up, beta = 1
    global accepted_partition_counter
    if accepted_partition_counter < 10000:
        beta = 0
    elif accepted_partition_counter < 70000:
        beta = (accepted_partition_counter - 10000) / 60000
    else:
        beta = 1

    Q2 = Q_func(partition.parent, partition)
    if Q2 == 0:
        return False
    Q1 = Q_func(partition, partition.parent)

    p = min(1, (Q1 / Q2) * exp(-beta * (score_function(partition) - score_function(partition.parent))))

    # accept the partition with probability p
    if random.random() < p:
        accepted_partition_counter += 1
        return True
    else:
        return False


def get_chain():
    '''
    Get data, build graph, make updaters and constraints, run chain
    :return: chain
    '''

    # NOTE: there are a few errors in some shapefiles; remove ignore_errors to see them
    graph = Graph.from_file(SHAPEFILE_PATH, ignore_errors=True)

    # delete islands (see gerrychain docs)
    components = list(connected_components(graph))
    biggest_component_size = max([len(c) for c in components])
    problem_components = [c for c in components if len(c) != biggest_component_size]
    for component in problem_components:
        for node in component:
            graph.remove_node(node)
            print('Removed a node from a smaller connected component\n')

    # Make updaters
    all_updaters = {
            'population': updaters.Tally(POPULATION_COL, alias='population'),
            'vap': updaters.Tally(VAP_COL, alias='vap'),
            'bvap': updaters.Tally(BVAP_COL, alias='bvap'),
            'hvap': updaters.Tally(HVAP_COL, alias='hvap'),
            'county_splits': updaters.county_splits('county_splits', COUNTY_COL),
            }

    elections = [
            Election(ELECTION_NAME,  {'Democratic': ELECTION_DEM_COL,   'Republican': ELECTION_REP_COL}, alias=ELECTION_NAME),
            ]

    all_updaters.update({election.name: election for election in elections})

    initial_partition = GeographicPartition(graph, assignment=DISTRICT_COL, updaters=all_updaters)

    ideal_population = sum(initial_partition['population'].values()) / len(initial_partition)

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
        lambda p: len(p['cut_edges']),
        2*len(initial_partition['cut_edges'])
    )

    global accepted_partition_counter
    accepted_partition_counter = 0

    chain = MarkovChain(
        proposal=proposal,
        constraints=[
            pop_constraint,
            compactness_bound
        ],
        accept=acceptance_function,
        # accept=accept.always_accept,
        initial_state=initial_partition,
        total_steps=TOTAL_STEPS
    )

    return chain


def get_chain_data(chain):
    return pandas.DataFrame(
            sorted(partition[ELECTION_NAME].percents('Democratic'))
            for partition in chain.with_progress_bar())


def display_chain_data(data):
    fig, ax = plt.subplots(figsize=(8, 6))

    # Draw 50% line
    ax.axhline(0.5, color='#cccccc')

    # Draw boxplot
    data.boxplot(ax=ax, positions=range(len(data.columns)))

    # Draw initial plan's Democratic vote %s (.iloc[0] gives the first row)
    plt.plot(data.iloc[0], 'ro')

    # Annotate
    ax.set_title('Comparing the 2011 plan to an ensemble')
    ax.set_ylabel('Democratic vote % (US House 2016)')
    ax.set_xlabel('Sorted districts')
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
    # print(partition['SSEN16'].wins('Democratic'))
    # print(partition['SSEN16'].counts('Democratic'))


def out_csv(chain):
    with open(OUTPUT_PATH, 'w', newline='') as csvfile:
        f = csv.writer(csvfile)
        headings = ['State', 'efficiency_gap', 'partisan_bias', 'd_seats', 'r_seats', 'black_opportunity', 'hisp_opportunity']
        # add vote count headings
        for i in range(TOTAL_DISTRICTS):
            headings.append(f'vap_{i}')
            headings.append(f'bvap_{i}')
            headings.append(f'hvap_{i}')
            headings.append(f'd_count_{i}')
            headings.append(f'r_count_{i}')
        f.writerow(headings)

        state = 1
        for partition in chain:

            # get election metrics, seats won and add to row
            e_g = partition[ELECTION_NAME].efficiency_gap()
            p_b = partition[ELECTION_NAME].partisan_bias()
            dem = partition[ELECTION_NAME].wins('Democratic')
            rep = partition[ELECTION_NAME].wins('Republican')
            b_opp_index = 5
            h_opp_index = 6
            row = [state, e_g, p_b, dem, rep, 0, 0]

            # get lists of vote counts and add to row
            dem_list = partition[ELECTION_NAME].counts('Democratic')
            rep_list = partition[ELECTION_NAME].counts('Republican')
            assert TOTAL_DISTRICTS == len(dem_list)
            black_opportunity_districts = 0
            hisp_opportunity_districts = 0
            for i in range(TOTAL_DISTRICTS):
                vap_key = str(i + 1)            # TODO make sure these correspond to correct district numbers
                if PAD_DISTRICT_NUMBERS:
                    vap_key = (len(str(TOTAL_DISTRICTS)) - len(vap_key)) * '0' + vap_key
                vap = partition.vap[vap_key]
                bvap = partition.bvap[vap_key]
                hvap = partition.hvap[vap_key]
                row.append(vap)
                row.append(bvap)
                row.append(hvap)
                row.append(dem_list[i])
                row.append(rep_list[i])
                if bvap / vap > OPPORTUNITY_THRESHOLD:
                    black_opportunity_districts += 1
                if hvap / vap > OPPORTUNITY_THRESHOLD:
                    hisp_opportunity_districts += 1
            row[b_opp_index] = black_opportunity_districts
            row[h_opp_index] = hisp_opportunity_districts
            f.writerow(row)
            state += 1


if __name__ == '__main__':
    chain = get_chain()
    data = get_chain_data(chain)
    # display_chain_data(data)

    # explore_chain(chain)

    out_csv(chain)


'''
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
'''
