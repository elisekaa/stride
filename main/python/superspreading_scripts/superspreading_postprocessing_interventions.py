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
"""

import argparse
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os

from statsmodels.nonparametric.smoothers_lowess import lowess

from util import get_experiment_ids, get_new_cases_per_day, get_secondary_cases_by_individual, save_figure
from util import plot_ar, plot_day_last_infection, plot_evolution, plot_final_size_frequencies

def plot_num_cases_over_period(output_dir, fig_name, all_num_cases, display_scenario_names, y_label):
    plt.violinplot(all_num_cases)
    plt.plot(range(1, len(all_num_cases) + 1), [np.mean(scen) for scen in all_num_cases], marker="o", linestyle="None")
    plt.xticks(range(1, len(all_num_cases) + 1), scenario_display_names, rotation=25)
    plt.ylabel(y_label)
    save_figure(output_dir, fig_name, extension="png")

def plot_resurgence_probabilities(output_dir, fig_name, all_num_case_after_release, display_scenario_names, extinction_threshold):
    resurgence_probs = []
    for scenario in all_num_case_after_release:
        resurgence_probs.append(len([run for run in scenario if run >= extinction_threshold]) / len(scenario))

    plt.bar(range(len(resurgence_probs)), resurgence_probs)
    plt.xticks(range(len(resurgence_probs)), display_scenario_names, rotation=25)
    plt.ylabel("Resurgence probability")
    plt.ylim((0, 1))
    save_figure(output_dir, fig_name)

def main(output_dir, scenario_names, display_scenario_names):
    extinction_threshold = 0 # FIXME
    num_days = 1100 # FIXME
    pop_size = 3000000 # FIXME

    if len(display_scenario_names) < len(scenario_names):
        display_scenario_names = scenario_names

    all_total_cases = []

    all_num_cases_lockdown = []
    all_num_cases_lockdown_exclude_extinct_before_d30 = []
    all_num_cases_release = []
    all_num_cases_release_exclude_extinct_before_d30 = []

    all_last_days_with_infections = []

    for scenario in scenario_names:
        print(scenario)

        experiment_ids = get_experiment_ids(output_dir, scenario)
        with multiprocessing.Pool(processes=4) as pool:
            # Get number of secondary cases for each infected individual
            secondary_cases_by_individual = pool.starmap(get_secondary_cases_by_individual,
                                                [(output_dir, scenario, exp_id, num_days) for exp_id in experiment_ids])
            # Get number of new cases per day
            new_cases_per_day = pool.starmap(get_new_cases_per_day,
                                                [(output_dir, scenario, exp_id, num_days) for exp_id in experiment_ids])

            # Calculate the final total number of cases for each simulation
            # Note: index cases are not counted.
            total_cases = [sum(x.values()) for x in secondary_cases_by_individual]
            all_total_cases.append(total_cases)

            # Calculate last day with new infections
            last_day_with_infections = []
            for run in new_cases_per_day:
                for day in range(num_days - 1, -1, -1):
                    if run[day] > 0:
                        last_day_with_infections.append(day)
                    elif day == 0:
                        last_day_with_infections.append(day)
            all_last_days_with_infections.append(last_day_with_infections_exclude_extinction)

            # Calculate number of cases during lockdown and after release
            num_cases_lockdown = [sum([run[d] for d in range(30, 100)]) for run in new_cases_per_day]
            num_cases_lockdown_exclude_extinction = [sum([run[d] for d in range(30, 100)]) for run in new_cases_per_day if run[30] > 0]

            num_cases_release = [sum([run[d] for d in range(100, 1100)]) for run in new_cases_per_day]
            num_cases_release_exclude_extinction = [sum([run[d] for d in range(100, 1100)]) for run in new_cases_per_day if run[30] > 0]

            all_num_cases_lockdown.append(num_cases_lockdown)
            all_num_cases_lockdown_exclude_extinct_before_d30.append(num_cases_lockdown_exclude_extinction)
            all_num_cases_release.append(num_cases_release)
            all_num_cases_release_exclude_extinct_before_d30.append(num_cases_release_exclude_extinction)

            # Plot evolution of new cases per day
            plot_evolution(output_dir, "evolution", scenario, new_cases_per_day, num_days)

    plot_day_last_infection(output_dir, "day_of_last_infection", display_scenario_names, all_last_days_with_infections, violin_plot=True)
    plot_ar(output_dir, "ar", display_scenario_names, all_total_cases, num_days, pop_size, 0, violin_plot=True)
    plot_final_size_frequencies(output_dir, "final_size_frequencies", display_scenario_names, all_total_cases, num_days)

    plot_num_cases_over_period(output_dir, "num_cases_lockdown", all_num_cases_lockdown, display_scenario_names, "Number of cases during lockdown")
    plot_num_cases_over_period(output_dir, "num_cases_release", all_num_cases_release, display_scenario_names, "Number of cases after relaxing of social distancing")
    plot_num_cases_over_period(output_dir, "num_cases_lockdown_exclude_extinction", all_num_cases_lockdown_exclude_extinct_before_d30, display_scenario_names, "Number of cases during lockdown")
    plot_num_cases_over_period(output_dir, "num_cases_release_exclude_extinction", all_num_cases_release_exclude_extinct_before_d30, display_scenario_names, "Number of cases after relaxing social distancing")

    plot_resurgence_probabilities(output_dir, "resurgence_probabilities", all_num_cases_release_exclude_extinct_before_d30, display_scenario_names, 500)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str, help="Directory containing simulation results")
    parser.add_argument("scenario_names", type=str, nargs="+", help="Names of scenarios to be postprocessed")
    parser.add_argument("--display_scenario_names", type=str, nargs="+", default=[], help="Names for scenarios to be displayed on plots")

    args = parser.parse_args()
    main(args.output_dir, args.scenario_names, args.display_scenario_names)
