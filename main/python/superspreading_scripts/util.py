############################################################################ #
#  This file is part of the Stride software.
#  It is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or any
#  later version.
#  The software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  You should have received a copy of the GNU General Public License,
#  along with the software. If not, see <http://www.gnu.org/licenses/>.
#  see http://www.gnu.org/licenses/.
#
#
#  Copyright 2020, Willem L, Kuylen E & Broeckhove J
############################################################################ #

"""
    Some utility functions for superspreading postprocessing scripts.
"""

import csv
import matplotlib.pyplot as plt
import os

def get_experiment_ids(output_dir, scenario_name):
    experiment_ids = []
    summary_file = os.path.join(output_dir, scenario_name, scenario_name + "_summary.csv")
    with open(summary_file) as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            exp_id = int(row["exp_id"])
            experiment_ids.append(exp_id)
    return experiment_ids

def get_trans_prob_by_exp(output_dir, scenario_name):
    experiments = {}
    summary_file = os.path.join(output_dir, scenario_name, scenario_name + "_summary.csv")
    with open(summary_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            exp_id = int(row["exp_id"])
            transmission_probability = float(row["transmission_probability"])
            experiments[exp_id] = transmission_probability

    return experiments

def get_cumulative_cases(output_dir, scenario_name, experiment_id, num_days, include_index_cases=False):
    log_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(experiment_id), "event_log.txt")
    cumulative_cases = 0
    with open(log_file) as f:
        for line in f:
            line = line.split(" ")
            tag = line[0]
            if tag == "[PRIM]" and include_index_cases:
                sim_day = int(line[6])
                if sim_day < num_days:
                    cumulative_cases += 1
            elif tag == "[TRAN]":
                sim_day = int(line[6])
                if sim_day < num_days:
                    cumulative_cases += 1
    return cumulative_cases

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
                if sim_day < num_days:
                    new_cases_per_day[sim_day] += 1
    return new_cases_per_day

def get_secondary_cases_by_individual(output_dir, scenario_name, experiment_id, num_days):
    log_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(experiment_id), "event_log.txt")
    potential_infectors = {}
    with open(log_file) as f:
        for line in f:
            line = line.split(" ")
            tag = line[0]
            sim_day = int(line[6])
            if sim_day < (num_days):
                if tag == "[PRIM]": # Index case
                    infected_id = int(float(line[1]))
                    if infected_id not in potential_infectors:
                        potential_infectors[infected_id] = 0
                elif tag == "[TRAN]": # Transmission
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

def plot_ar(output_dir, fig_name, display_scenario_names, all_total_cases, num_days, pop_size, extinction_threshold=0, violin_plot=False):
    all_total_cases_prop_pop = []
    for scenario in all_total_cases:
        all_total_cases_prop_pop.append([total_cases / pop_size for total_cases in scenario if total_cases >= extinction_threshold])

    if violin_plot:
        plt.violinplot(all_total_cases_prop_pop)
        plt.xticks(range(1, len(all_total_cases_prop_pop) + 1), display_scenario_names, rotation=25)
    else:
        plt.boxplot(all_total_cases_prop_pop, labels=display_scenario_names)
        plt.xticks(rotation=25)

    plt.ylabel("AR (after {} days)".format(num_days))

    save_figure(output_dir, fig_name, extension="png")

def plot_day_last_infection(output_dir, fig_name, display_scenario_names, all_last_days_with_infections, violin_plot=False):
    if violin_plot:
        plt.violinplot(all_last_days_with_infections)
        plt.xticks(range(1, len(all_last_days_with_infections) + 1), display_scenario_names, rotation=25)
    else:
        plt.boxplot(all_last_days_with_infections, labels=display_scenario_names)
        plt.xticks(rotation=25)

    plt.plot(range(1, len(all_last_days_with_infections) + 1), [np.mean(x) for x in all_last_days_with_infections], linestyle="None", marker="o")

    plt.ylabel("Last day with new infections")
    save_figure(output_dir, fig_name, extension="png")


def plot_evolution(output_dir, fig_name, scenario_name, new_cases_per_day, num_days, y_max):
    for run in new_cases_per_day:
        plt.plot(range(num_days), [run[day] for day in range(num_days)])

    plt.xlabel("Simulation day")
    plt.ylabel("New cases")
    plt.ylim(0, y_max)

    save_figure(output_dir, fig_name + "_" + scenario_name)

def plot_final_size_frequencies(output_dir, fig_name, display_scenario_names, all_total_cases, num_days):
    """
        Plot frequency with which final sizes occur.
        Visualisation for extinction threshold.
    """

    plt.hist(all_total_cases, histtype="bar", stacked=True)

    plt.xlabel("Outbreak size after {} days".format(num_days))
    plt.ylabel("Frequency")
    plt.legend(display_scenario_names)

    save_figure(output_dir, fig_name)

def save_figure(output_dir, figure_name, extension="eps", dpi=1000):
    if not os.path.exists(os.path.join(output_dir, "fig")):
        os.mkdir(os.path.join(output_dir, "fig"))

    plt.savefig(os.path.join(output_dir, "fig", figure_name + "." + extension), format=extension, bbox_inches='tight', dpi=dpi)
    plt.clf()
