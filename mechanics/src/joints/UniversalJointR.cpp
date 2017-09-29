/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2017 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*! \file UniversalJointR.cpp
*/

#include "UniversalJointR.hpp"
#include <NewtonEulerDS.hpp>
#include <Interaction.hpp>
#include <boost/math/quaternion.hpp>
#include <BlockVector.hpp>

#include <iostream>

// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

// Calculate the angle between two axes, and finding a vector exactly
// 90 degrees from the first axis in the same plane.
static
SP::SiconosVector calcOrthogonal(SP::SiconosVector A1, SP::SiconosVector A2)
{
}

UniversalJointR::UniversalJointR()
  : Pivot2JointR()
{
  setSuspensionFree(false);
}

UniversalJointR::UniversalJointR(SP::SiconosVector P, SP::SiconosVector A1,
                                 SP::SiconosVector A2, bool absoluteRef,
                                 SP::NewtonEulerDS d1, SP::NewtonEulerDS d2)
  : Pivot2JointR(P, A1, A2, false, absoluteRef, d1, d2)
{}

void UniversalJointR::setBasePositions(SP::SiconosVector q1,
                                       SP::SiconosVector q2)
{
  // Set the second axis internally by calculating the angle between
  // the two axes and finding the second axis in the same plane and
  // same direction but at exactly 90 degrees.
  SP::SiconosVector _axis1 = calcOrthogonal(_axes[0], _axes[1]);
  Pivot2JointR::setBasePositions(q1, q2);
}
