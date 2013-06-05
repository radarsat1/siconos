/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
#include "NonSmoothEvent.hpp"
#include "EventFactory.hpp"
#include "EventDriven.hpp"
#include "OneStepIntegrator.hpp"
#include "LagrangianDS.hpp"


using namespace EventFactory;

// Default constructor
NonSmoothEvent::NonSmoothEvent(): Event(0.0, NS_EVENT)
{}

NonSmoothEvent::NonSmoothEvent(double time, int): Event(time, NS_EVENT)
{}

NonSmoothEvent::~NonSmoothEvent()
{}

void NonSmoothEvent::process(Simulation& simulation)
{
  if (Type::value(simulation) != Type::EventDriven)
    RuntimeException::selfThrow("NonSmoothEvent::process failed; Simulation is not of EventDriven type.");

  if (!(simulation.oneStepNSProblems()->empty()))
  {
    EventDriven& eventDriven = static_cast<EventDriven&>(simulation);

    // Compute y[0], y[1] and update index sets. => already done
    // during advance to event ...

    //       simulation->updateOutput(0, 1);

    //       simulation->updateIndexSets();

    // Get the required index sets ...
    SP::InteractionsGraph indexSet0 = simulation.indexSet(0);
    SP::DynamicalSystemsGraph dsG = simulation.model()->nonSmoothDynamicalSystem()->topology()->dSG(0);
    DynamicalSystemsGraph::VIterator vi, viend;

    // Update all the index sets ...
    eventDriven.updateIndexSets();
    SP::InteractionsGraph indexSet1 = simulation.indexSet(1);
    SP::InteractionsGraph indexSet2 = simulation.indexSet(2);
    bool found = true;
    InteractionsGraph::VIterator ui, uiend;
    for (std11::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui)
    {
      found = indexSet2->is_vertex(indexSet1->bundle(*ui));
      if (!found) break;
    }
    // ---> solve impact LCP if IndexSet[1]\IndexSet[2] is not empty.
    if (!found)
    {
      // For Event-Driven algo., memories vectors are of size 2
      // (ie 2 interactionBlocks).  First interactionBlock (pos 0, last
      // in) for post-event values and last interactionBlock (pos 1,
      // first in) for pre-event values.

      simulation.saveInMemory();  // To save pre-impact values

      // solve the LCP-impact => y[1],lambda[1]
      eventDriven.computeOneStepNSProblem(SICONOS_OSNSP_ED_IMPACT); // solveLCPImpact();
      // compute p[1], post-impact velocity, y[1] and indexSet[2]
      simulation.update(1);
      // Update the corresponding index set ...
      eventDriven.updateIndexSets();
    }

    /*
    // Display the variable after processing NSEvent
    std::cout<< "-------After processing NS events---------" <<std::endl;
    for (std11::tie(ui, uiend)=indexSet0->vertices(); ui != uiend; ++ui)
      {
        SP::Interaction inter = indexSet0->bundle(*ui);
        std::cout << "Velocity at this Interaction: " << (*inter->y(1))(0) <<std::endl;
      }

    for (std11::tie(vi, viend) = dsG->vertices(); vi != viend; ++vi)
      {
        SP::DynamicalSystem ds = dsG->bundle(*vi);
        SP::LagrangianDS Lag_ds = std11::static_pointer_cast<LagrangianDS>(ds);
        std::cout << "Velocity of DS: " << (*Lag_ds->velocity())(0) <<std::endl;
      }
    //
    */
    //---> solve acceleration LCP if IndexSet[2] is not empty
    if (indexSet2->size() > 0)
    {
      // // Update the state of the DS
      // OSIIterator itOSI; for(itOSI = simulation->oneStepIntegrators()->begin();
      //                        itOSI!= simulation->oneStepIntegrators()->end() ;
      //                        ++itOSI)
      //                           (*itOSI)->updateState(2);

      // solve LCP-acceleration
      eventDriven.computeOneStepNSProblem(SICONOS_OSNSP_ED_SMOOTH_ACC); // solveLCPAcceleration();
      // update input of level 2, acceleration and output of level 2
      simulation.update(2);
      // for all index in IndexSets[2], update the index set according to y[2] and/or lambda[2] sign.
      eventDriven.updateIndexSetsWithDoubleCondition();
    }

    // Save results in memory
    simulation.saveInMemory();
  }
}

AUTO_REGISTER_EVENT(NS_EVENT, NonSmoothEvent)
