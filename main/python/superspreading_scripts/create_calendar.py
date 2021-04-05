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
    Script to create holidays_file for scenario with period of social distancing
    for superspreading project.
    Start date: November 1, 2020
    Next: 30 days without interventions,
            followed by 70 days with schools_closed,
            workplace_distancing and community_distancing.
"""

import csv
import time

from datetime import timedelta, date

def main():

    start_date = date(2020, 11, 1) # November 1, 2020

    start_interventions = 30
    interventions_duration = 70

    start_interventions_date = start_date + timedelta(start_interventions)

    with open("calendar_social_distancing_no_holidays.csv", "w") as csvfile:
        fieldnames = ["category", "date", "value", "type", "age"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()

        for d in range(interventions_duration):
            day = start_interventions_date + timedelta(d)
            # Schools closed ages 0-25
            for age in range(26):
                writer.writerow({
                    "category": "schools_closed",
                    "date": day,
                    "value": 1,
                    "type": "boolean",
                    "age": age
                })
            # Workplace distancing
            writer.writerow({
                "category": "workplace_distancing",
                "date": day,
                "value": 1,
                "type": "boolean",
                "age": "NA"
            })
            # Community distancing
            writer.writerow({
                "category": "community_distancing",
                "date": day,
                "value": 1,
                "type": "boolean",
                "age": "NA"
            })


if __name__=="__main__":
    main()
