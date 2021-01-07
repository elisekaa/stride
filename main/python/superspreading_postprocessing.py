import argparse
import collections
import csv
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os
import statistics

# TODO offspring distributions
# TODO distributions of individual transmission probabilities?

def get_experiment_ids(output_dir, scenario_name):
    experiment_ids = []
    summary_file = os.path.join(output_dir, scenario_name, scenario_name + "_summary.csv")
    with open(summary_file) as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            exp_id = int(row["exp_id"])
            experiment_ids.append(exp_id)
    return experiment_ids

def get_secondary_cases_by_individual(output_dir, scenario_name, experiment_id, num_days):
    log_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(experiment_id), "event_log.txt")
    potential_infectors = {}
    with open(log_file) as f:
        for line in f:
            line = line.split(" ")
            tag = line[0]
            sim_day = int(line[6])
            if sim_day < num_days:
                if tag == "[PRIM]":
                    infected_id = int(float(line[1]))
                    if infected_id not in potential_infectors:
                        potential_infectors[infected_id] = 0
                elif tag == "[TRAN]":
                    infector_id = int(float(line[2]))
                    infected_id = int(float(line[1]))
                    # Add infected to potential infectors
                    if infected_id not in potential_infectors:
                        potential_infectors[infected_id] = 0
                    # Add infector to infectors (if not yet done)
                    # And add this infection to total secondary cases count
                    if infector_id not in potential_infectors:
                        potential_infectors[infector_id] = 1
                    else:
                        potential_infectors[infector_id] += 1
    return potential_infectors

def get_effective_r_per_day(output_dir, scenario_name, experiment_id, num_days):
    infected_by_day = {}
    potential_infectors = {}

    for day in range(num_days):
        infected_by_day[day] = []

    log_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(experiment_id), "event_log.txt")
    with open(log_file) as f:
        for line in f:
            line = line.split(" ")
            tag = line[0]
            if tag == "[PRIM]":
                infected_id = int(float(line[1]))
                sim_day = int(line[6])
                infected_by_day[sim_day].append(infected_id)
                if infected_id not in potential_infectors:
                    potential_infectors[infected_id] = 0
            elif tag == "[TRAN]":
                infected_id = int(float(line[1]))
                sim_day = int(line[6])
                infected_by_day[sim_day].append(infected_id)

                infector_id = int(float(line[2]))
                # Add infected to potential infectors
                if infected_id not in potential_infectors:
                    potential_infectors[infected_id] = 0
                # Add infector to infectors (if not yet done)
                # And add this infection to total secondary cases count
                if infector_id not in potential_infectors:
                    potential_infectors[infector_id] = 1
                else:
                    potential_infectors[infector_id] += 1

    rt_by_day = {}
    for day in infected_by_day:
        infected = infected_by_day[day]
        secondary_cases = []
        for infector_id in infected:
            secondary_cases.append(potential_infectors[infector_id])
        if len(secondary_cases) > 0:
            rt = sum(secondary_cases) / len(secondary_cases)
        else:
            rt = 0
        rt_by_day[day] = rt
    return rt_by_day

def get_new_cases_per_day(output_dir, scenario_name, experiment_id, num_days):
    log_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(experiment_id), "event_log.txt")
    new_cases_per_day = {}
    for day in range(num_days):
        new_cases_per_day[day] = 0
    with open(log_file) as f:
        for line in f:
            line = line.split(" ")
            tag = line[0]
            if tag == "[TRAN]":
                sim_day = int(line[6])
                new_cases_per_day[sim_day] += 1
    return new_cases_per_day

def get_herd_immunity_threshold(output_dir, scenario_name, experiment_id, num_days, extinction_threshold=0):
    rt_by_day = get_effective_r_per_day(output_dir, scenario_name, experiment_id, num_days)
    new_cases_per_day = get_new_cases_per_day(output_dir, scenario_name, experiment_id, num_days)
    cumulative_cases = 0.0
    herd_immunity_threshold = 0
    for day in range(num_days):
        cumulative_cases += new_cases_per_day[day]
        rt = rt_by_day[day]
        if rt < 1 and day > 50 and cumulative_cases > extinction_threshold:
            herd_immunity_threshold = cumulative_cases / 3000000
            break
    return herd_immunity_threshold

def get_p80(secondary_cases_by_individual, extinction_threshold=0):
    p80s = []
    for run in secondary_cases_by_individual:
        total_cases = sum(run.values())
        if total_cases >= extinction_threshold:
            # Order number of secondary cases by individual in decreasing order
            secondary_cases_sorted = list(run.values())
            secondary_cases_sorted.sort(reverse=True)

            num_cases_responsible = 0
            num_cases_caused = 0
            for s in secondary_cases_sorted:
                num_cases_caused += s
                num_cases_responsible += 1
                if num_cases_caused >= (total_cases * 0.80):
                    break
            p80s.append(num_cases_responsible / total_cases)
    return p80s

def plot_p80(p80s, scenario_display_names):
    plt.boxplot(p80s, labels=scenario_display_names)
    plt.xticks(rotation=25)
    plt.ylabel("P80")
    plt.tight_layout()
    plt.show()

def plot_secondary_cases_frequency(scenario_name, all_secondary_cases, extinction_threshold=0):
    for sim_run in all_secondary_cases:
        secondary_cases = sim_run.values()
        total_infections = sum(secondary_cases)
        # Count frequencies of number of secondary cases / person
        counter = collections.Counter(secondary_cases)
        # Scale frequencies to percentage of total infected cases
        if total_infections > 0:
            for num_secondary_cases in counter:
                counter[num_secondary_cases] = (counter[num_secondary_cases] / total_infections) * 100
        num_secondary_cases_sorted = list(counter.keys())
        num_secondary_cases_sorted.sort()
        # Plot results for this simulation run
        plt.plot(num_secondary_cases_sorted, [counter[x] for x in num_secondary_cases_sorted], "o")
    plt.xlabel("Number of secondary cases")
    plt.xlim(-2, 85)
    plt.ylabel("Percentage of infected individuals")
    plt.ylim(-2, 105)
    plt.title(scenario_name)
    plt.show()

def plot_extinction_probabilities(all_total_cases, scenario_display_names, extinction_threshold):
    extinction_probabilities = []
    for i in range(len(scenario_display_names)):
        total_cases = all_total_cases[i]
        extinction_probabilities.append((len([x for x in total_cases if x < extinction_threshold]) / len(total_cases)))

    plt.bar(range(len(scenario_display_names)), extinction_probabilities)
    print(extinction_probabilities)
    plt.xticks(range(len(scenario_display_names)), scenario_display_names, rotation=25)
    plt.ylabel("Extinction probability (threshold = {} cases)".format(extinction_threshold))
    plt.ylim(0, 1.1)
    plt.tight_layout()
    plt.show()

def plot_ar(all_total_cases, scenario_display_names, num_days, extinction_threshold):
    all_total_cases_prop_pop = []
    for scenario in all_total_cases:
        all_total_cases_prop_pop.append([total_cases / 3000000 for total_cases in scenario if total_cases >= extinction_threshold])

    plt.boxplot(all_total_cases_prop_pop, labels=scenario_display_names)
    plt.xticks(rotation=25)
    plt.ylabel("AR (after {} days)".format(num_days))
    plt.ylim(0.90, 0.93)
    plt.tight_layout()
    plt.show()

def plot_evolution(new_cases_per_day, num_days):
    for run in new_cases_per_day:
        plt.plot(range(num_days), [run[day] for day in range(num_days)])
    plt.xlabel("Simulation day")
    plt.ylabel("New cases")
    plt.tight_layout()
    plt.show()

def plot_effective_r_by_day(effective_r_by_day, num_days, extinction_threshold):
    effective_r_by_day_no_extinction = [run for run in effective_r_by_day if sum(run.values()) > extinction_threshold]

    mean_effective_rs = []
    median_effetive_rs = []
    lower_effective_rs = []
    upper_effective_rs = []
    for day in range(num_days):
        all_effective_rs_on_day = [run[day] for run in effective_r_by_day_no_extinction]
        mean = sum(all_effective_rs_on_day) / len(all_effective_rs_on_day)
        mean_effective_rs.append(mean)
        median = statistics.median(all_effective_rs_on_day)
        median_effetive_rs.append(median)
        ll = np.percentile(all_effective_rs_on_day, 2.5)
        lower_effective_rs.append(ll)
        hh = np.percentile(all_effective_rs_on_day, 97.5)
        upper_effective_rs.append(hh)

    plt.fill_between(range(num_days), lower_effective_rs, upper_effective_rs, color="lightgrey")
    plt.plot(range(num_days), mean_effective_rs)
    plt.plot(range(num_days), [1] * num_days)
    plt.xlabel("Simulation day")
    plt.ylabel("Rt")
    plt.ylim(0, 30)
    plt.show()

def plot_herd_immunity_threshold(herd_immunity_thresholds, scenario_display_names):
    plt.boxplot(herd_immunity_thresholds, labels=scenario_display_names)
    plt.xticks(rotation=25)
    plt.ylabel("Herd immunity threshold")
    plt.tight_layout()
    plt.show()

def main(output_dir, scenario_names, scenario_display_names):
    extinction_threshold = 40
    num_days = 180

    if scenario_display_names is None:
        scenario_display_names = scenario_names

    all_p80s = []
    all_total_cases = []
    all_total_cases_prop_pop = []
    all_herd_immunity_thresholds = []
    for scenario in scenario_names:
        print(scenario)
        experiment_ids = get_experiment_ids(output_dir, scenario)
        with multiprocessing.Pool(processes=4) as pool:
            secondary_cases = pool.starmap(get_secondary_cases_by_individual,
                                           [(output_dir, scenario, exp_id, num_days) for exp_id in experiment_ids])
            #effective_r_by_day = pool.starmap(get_effective_r_per_day,
            #                                [(output_dir, scenario, exp_id, num_days) for exp_id in experiment_ids])
            #herd_immunity_thresholds = pool.starmap(get_herd_immunity_threshold,
            #                                [(output_dir, scenario, exp_id, num_days, extinction_threshold) for exp_id in experiment_ids])
            #all_herd_immunity_thresholds.append([x for x in herd_immunity_thresholds if x > 0])
            #new_cases_per_day = pool.starmap(get_new_cases_per_day,
            #                                [(output_dir, scenario, exp_id, num_days) for exp_id in experiment_ids])

            #plot_evolution(new_cases_per_day, num_days)
            #plot_effective_r_by_day(effective_r_by_day, num_days, extinction_threshold)

            total_cases = [sum(x.values()) for x in secondary_cases]
            all_total_cases.append(total_cases)

            #total_cases.sort(reverse=True)
            #print(total_cases)

            #all_p80s.append(get_p80(secondary_cases, extinction_threshold=extinction_threshold))
            #plot_secondary_cases_frequency(scenario, secondary_cases, extinction_threshold=1)

    #plot_extinction_probabilities(all_total_cases, scenario_display_names, extinction_threshold)
    #plot_p80(all_p80s, scenario_display_names)
    plot_ar(all_total_cases, scenario_display_names, num_days, extinction_threshold)
    #plot_herd_immunity_threshold(all_herd_immunity_thresholds, scenario_display_names)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str)
    parser.add_argument("scenario_names", type=str, nargs="+")
    parser.add_argument("--scenario_display_names", type=str, nargs="+", default=None)
    args = parser.parse_args()
    main(args.output_dir, args.scenario_names, args.scenario_display_names)
