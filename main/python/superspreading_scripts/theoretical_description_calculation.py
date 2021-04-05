import argparse
import csv
import multiprocessing
import numpy as np
import random
import time
import xml.etree.ElementTree as ET

from scipy import integrate
from scipy.stats import gamma

def get_contact_rates(pooltype, contact_rate_tree, maxAge):
    # Create matrix of zeroes
    contact_rates = []
    for age in range(maxAge + 1):
        contact_rates.append(0)

    # Get contact rates for the right pooltype from xml tree
    max_participant_age = 0
    pooltype_contacts = contact_rate_tree.find(pooltype)
    # This is for aggregated contacts by age (participants -> contacts -> contact -> age = all)
    for participant in pooltype_contacts.findall("participant"):
        participant_age = int(participant.find("age").text)
        if (participant_age > max_participant_age):
            max_participant_age = participant_age
        contact_rate = float(participant.find("contacts/contact/rate").text)
        contact_rates[participant_age] = contact_rate

    # Fill in contact rates for ages > max_participant_age
    # Set contact_rate[age] = contact_rate[max_participant_age] if age > max_participant_age
    for age in range(max_participant_age + 1, maxAge + 1):
        contact_rates[age] = contact_rates[max_participant_age]

    return contact_rates

def add_to_pool(pools, pool_id, person_id, age):
    if pool_id > 0:
        if pool_id in pools:
            pools[pool_id].append((person_id, age))
        else:
            pools[pool_id] = [(person_id, age)]

def sum_pools_and_contacts(function, person, all_pools, contact_rates, infectious_period_length,
    transmission_probability, mean_transmission_probability, overdispersion, estimated_mean=None):

    result = 0

    person_id = person["id"]
    age = person["age"]

    # Initialize gamma distribution if necessary
    if overdispersion is not None:
        shape = overdispersion
        scale = mean_transmission_probability / shape
        pdf_tp = gamma.pdf(transmission_probability, shape, scale=scale)
        cdf1 = gamma.cdf(1, shape, scale=scale)
        cdf0 = gamma.cdf(0, shape, scale=scale)

    # Iterate over contact pools this person belongs to
    for pool_type, pools in all_pools.items():
        pool_id = person[pool_type + "_id"]
        if pool_id > 0:
            pool_members = pools[pool_id]
            pool_size = len(pool_members)
            # Iterate over all members in this contact pool
            for member in pool_members:
                member_id = member[0]
                member_age = member[1]
                if member_id != person_id: # Check that this is not the same person
                    # This is for aggregated contacts by age (participants -> contacts -> contact -> age = all)
                    contact_rate1 = contact_rates[pool_type][age]
                    contact_rate2 = contact_rates[pool_type][member_age]
                    contact_probability1 = contact_rate1 / (pool_size - 1)
                    contact_probability2 = contact_rate2 / (pool_size - 1)
                    contact_probability = min(contact_probability1, contact_probability2)
                    # Households are assumed to be fully connected in Stride
                    if pool_type == "household":
                        contact_probability = 0.999
                    if contact_probability >= 1:
                        contact_probability = 0.999

                    # Function to sum over
                    if overdispersion is None:
                        if function == "mean":
                            result += (1 - (1 - (mean_transmission_probability * contact_probability))**infectious_period_length)
                        elif function == "variance":
                            result += ((1 - (mean_transmission_probability * contact_probability))**infectious_period_length) * (1 - (1 - mean_transmission_probability * contact_probability)**infectious_period_length)
                    else:
                        if function == "mean":
                            result += (1 - ((1 - (transmission_probability * contact_probability))**infectious_period_length) * (pdf_tp / (cdf1 - cdf0)))
                        elif function == "ev": # E[Var(Y | X)]
                            result += ((1 - (transmission_probability * contact_probability))**infectious_period_length) * (1 - (1 - transmission_probability * contact_probability)**infectious_period_length) * (pdf_tp / (cdf1 - cdf0))
                        elif function == "ve": # Var(E[Y | X])
                            result += (1 - (1 - transmission_probability * contact_probability)**infectious_period_length)

    # If we are calculating the variance of the expected value Var(E[Y | X])
    if function == "ve":
        result = ((result - estimated_mean)**2) * (pdf_tp / (cdf1 - cdf0))

    return (person_id, result)


def integrate_function(function, person, pools, contact_rates, infectious_period_length, mean_transmission_probability, overdispersion, estimated_mean=None):

    result = 0

    if function == "mean":
        func = lambda x: (sum_pools_and_contacts("mean", person, pools, contact_rates, infectious_period_length, x, mean_transmission_probability, overdispersion))[1]
        integral, upper_error = integrate.quad(func, 0, 1, epsabs=1.49e-2)
        result = integral
    elif function == "variance":
        func_ev = lambda x: (sum_pools_and_contacts("ev", person, pools, contact_rates, infectious_period_length, x, mean_transmission_probability, overdispersion))[1]
        func_ve = lambda x: (sum_pools_and_contacts("ve", person, pools, contact_rates, infectious_period_length, x, mean_transmission_probability, overdispersion, estimated_mean))[1]

        ev_integral, ev_upper_error = integrate.quad(func_ev, 0, 1, epsabs=1.49e-2)
        ve_integral, ve_upper_error = integrate.quad(func_ve, 0, 1, epsabs=1.49e-2)

        result = ev_integral + ve_integral

    return (person["id"], result)


def estimate_effective_contacts(population_file, contact_matrix_file, transmission_probabilities,
    infectious_period_length, overdispersion=None, person_ids=[]):
    maxAge = 111

    ############################################################################
    # Read contact matrices from file                                          #
    ############################################################################

    contact_rates = {
        "household": [],
        "work": [],
        "school": [],
        "primary_community": [],
        "secondary_community": []
    }

    # Parse xml file
    contact_matrix_tree = ET.parse(contact_matrix_file).getroot()
    contact_rates["household"] = get_contact_rates("household", contact_matrix_tree, maxAge)
    contact_rates["work"] = get_contact_rates("work", contact_matrix_tree, maxAge)
    contact_rates["school"] = get_contact_rates("school", contact_matrix_tree, maxAge)
    contact_rates["primary_community"] = get_contact_rates("primary_community", contact_matrix_tree, maxAge)
    contact_rates["secondary_community"] = get_contact_rates("secondary_community", contact_matrix_tree, maxAge)

    ############################################################################
    # From population file, get age constitutions                              #
    # and sizes of different contact pools                                     #
    ############################################################################

    households = {}
    workplaces = {}
    schools = {}
    primary_communities = {}
    secondary_communities = {}

    population = []
    person_id = 1

    with open(population_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            age = int(row["age"])
            # Keep track of pool constitutions
            add_to_pool(households, int(float(row["household_id"])), person_id, age)
            add_to_pool(workplaces, int(float(row["work_id"])), person_id, age)
            add_to_pool(schools, int(float(row["school_id"])), person_id, age)
            add_to_pool(primary_communities, int(float(row["primary_community"])), person_id, age)
            add_to_pool(secondary_communities, int(float(row["secondary_community"])), person_id, age)
            # Add person to population list
            population.append({"id": person_id, "age": age,
                                "household_id": int(float(row["household_id"])),
                                "work_id": int(float(row["work_id"])),
                                "school_id": int(float(row["school_id"])),
                                "primary_community_id": int(float(row["primary_community"])),
                                "secondary_community_id": int(float(row["secondary_community"]))})
            person_id += 1

    pools = {
        "household": households,
        "work": workplaces,
        "school": schools,
        "primary_community": primary_communities,
        "secondary_community": secondary_communities
    }

    ############################################################################
    # For each transmission probability to be tested:                          #
    # iterate over all persons in the population,                              #
    # and calculate the number of secondary cases they would cause             #
    # if infected in a completely susceptible population                       #
    # (cfr. theoretical description by A. Torneri)                             #
    ############################################################################
    mean_effective_contacts_by_tp = {}
    variance_effective_contacts_by_tp = {}

    for tp in transmission_probabilities:
        if tp != 0:
            all_effective_contacts = []
            selected_population = population[:2]

            # Check if calculation has to happen for a given selection of persons
            # instead of for entire population.
            if len(person_ids) > 0:
                selected_population = [population[i] for i in person_ids]

            with multiprocessing.Pool(processes=4) as pool:
                if overdispersion is None: # Calculate without distribution for individual transmission probability
                    mean_effective_contacts = pool.starmap(sum_pools_and_contacts,
                                                        [("mean", person, pools, contact_rates, infectious_period_length, tp,
                                                        tp, overdispersion) for person in selected_population])
                    var_effective_contacts = pool.starmap(sum_pools_and_contacts,
                                                        [("variance", person, pools, contact_rates, infectious_period_length, tp,
                                                        tp, overdispersion) for person in selected_population])
                else:
                    mean_effective_contacts = pool.starmap(integrate_function,
                                                        [("mean", person, pools, contact_rates,
                                                        infectious_period_length, tp, overdispersion) for person in selected_population])
                    var_effective_contacts = pool.starmap(integrate_function,
                                                        [("variance", person, pools, contact_rates,
                                                        infectious_period_length, tp, overdispersion,
                                                        next(person_result[1] for person_result in mean_effective_contacts if person_result[0] == person["id"])) for person in selected_population])

                mean_effective_contacts_by_tp[tp] = np.mean([person_result[1] for person_result in mean_effective_contacts])
                variance_effective_contacts_by_tp[tp] = np.mean([person_result[1] for person_result in var_effective_contacts])
        else:
            mean_effective_contacts_by_tp[tp] = 0
            variance_effective_contacts_by_tp[tp] = 0

    return (mean_effective_contacts_by_tp, variance_effective_contacts_by_tp)

def main(population_file, contact_matrix_file, transmission_probabilities, infectious_period_length, overdispersion):
    begin = time.perf_counter()

    mean_effective_contacts, var_effective_contacts = estimate_effective_contacts(population_file, contact_matrix_file,
                                                                    transmission_probabilities,
                                                                    infectious_period_length,
                                                                    overdispersion=overdispersion)

    end = time.perf_counter()
    print("Run time = {} seconds".format(end - begin))

    for tp in mean_effective_contacts:
        print("{}: mean R0 estimated to be {}".format(tp, mean_effective_contacts[tp]))
    for tp in var_effective_contacts:
        print("{}: variance of R0 estimated to be {}".format(tp, var_effective_contacts[tp]))

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("population_file", type=str)
    parser.add_argument("contact_matrix_file", type=str)
    parser.add_argument("--transmission_probabilities", type=float, nargs="+", default=[0.025, 0.05, 0.075, 0.10])
    parser.add_argument("--infectious_period_length", type=int, default=7, help="Mean length of infectious period (in days)")
    parser.add_argument("--overdispersion", type=float, default=None)

    args = parser.parse_args()
    main(args.population_file, args.contact_matrix_file, args.transmission_probabilities, args.infectious_period_length, args.overdispersion)
