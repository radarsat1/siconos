/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
/*! \file SiconosKernelExport.hpp
Python export with pyplusplus
*/

#ifndef SiconosKernelExport_hpp
#define SiconosKernelExport_hpp

#include "Model.h"
#include "Utils.hpp"
#include "ModelingTools.hpp"
#include "SimulationTools.hpp"
//#include "ControlTools.hpp"
#include "XmlTools.hpp"
#include "PluginTypes.hpp"


namespace pyplusplus
{
namespace aliases
{

typedef SiconosSet<DynamicalSystem, int> DynamicalSystemsSet;
/*    typedef boost::shared_ptr<SiconosSet<DynamicalSystem,int > > SPtrDynamicalSystemsSet;
typedef boost::shared_ptr<SiconosSet<UnitaryRelation,double* > > SPtrUnitaryRelationsSet;*/
}
}



#endif