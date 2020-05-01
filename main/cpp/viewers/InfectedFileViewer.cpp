/*
 *  This is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *  The software is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License
 *  along with the software. If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright 2018, Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file
 * Definition of Observer for SimEvents for commandline interface usage.
 */

#include "InfectedFileViewer.h"

#include "pop/Population.h"
#include "sim/Sim.h"
#include "sim/SimRunner.h"

using namespace std;
using namespace stride::sim_event;

namespace stride {
namespace viewers {

InfectedFileViewer::InfectedFileViewer(std::shared_ptr<SimRunner> runner, const std::string& output_prefix)
            : m_infected(), m_infected_file(output_prefix,"infected"),
			  m_exposed(), m_exposed_file(output_prefix,"exposed"),
			  m_infectious(),m_infectious_file(output_prefix,"infectious"),
			  m_symptomatic(),m_symptomatic_file(output_prefix,"symptomatic"),
			  m_infected_total(),m_infected_total_file(output_prefix,"cases"),
			  m_runner(std::move(runner))
        {
        }



void InfectedFileViewer::Update(const sim_event::Id id)
{
        switch (id) {
        case Id::AtStart:
        case Id::Stepped: {
                const auto pop = m_runner->GetSim()->GetPopulation();
                m_infected.push_back(pop->CountInfectedCases());
                m_exposed.push_back(pop->CountExposedCases());
                m_infectious.push_back(pop->CountInfectiousCases());
                m_symptomatic.push_back(pop->CountSymptomaticCases());
                m_infected_total.push_back(pop->GetTotalInfected());
                break;
        }
        case Id::Finished: {
                m_infected_file.Print(m_infected);
                m_exposed_file.Print(m_exposed);
                m_infectious_file.Print(m_infectious);
                m_symptomatic_file.Print(m_symptomatic);
                m_infected_total_file.Print(m_infected_total);
                break;
        }
        default: break;
        }
}

} // namespace viewers
} // namespace stride
