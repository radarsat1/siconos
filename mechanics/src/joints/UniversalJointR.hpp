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

/*! \file UniversalJointR.hpp
*/

#ifndef UniversalJointRELATION_H
#define UniversalJointRELATION_H

#include <MechanicsFwd.hpp>
#include <SiconosFwd.hpp>
#include <Pivot2JointR.hpp>

/** \class UniversalJointR \brief This class implements a universal
 *  joint between one or two Newton/Euler Dynamical system.  It is
 *  exactly a Pivot2JointR with no suspension and orthogonal axes,
 *  i.e. two bodies connected through a cross-shaped coupler, without
 *  the need for an extra supporting body.
 */
class UniversalJointR : public Pivot2JointR
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(UniversalJointR);

public:

  /** Empty constructor. The relation may be initialized later by
   * setPoint, setAxis, setAbsolute, and setBasePositions. */
  UniversalJointR();

  /** Constructor based on one or two dynamical systems, a point and an axis.
   *  \param d1 first DynamicalSystem linked by the joint.
   *  \param d2 second DynamicalSystem linked by the joint, or NULL
   *            for absolute frame.
   *  \param P SiconosVector of size 3 that defines the point around
   *           which rotation is allowed.
   *  \param A1 SiconosVector of size 3 that defines the primary universal
   *            axis.  The second axis of rotation is defined as the
   *            vector orthogonal to it starting at P, in the plane A1-A2.
   *  \param A2 SiconosVector of size 3 that defines the plane A1-A2,
   *            from which the secondary universal axis is derived as
   *            90 degrees from A1.
   *  \param absoluteRef if true, P and A are in the absolute frame,
   *                     otherwise P and A are in d1 frame.
   */
  UniversalJointR(SP::SiconosVector P, SP::SiconosVector A1,
                  SP::SiconosVector A2, bool absoluteRef,
                  SP::NewtonEulerDS d1 = SP::NewtonEulerDS(),
                  SP::NewtonEulerDS d2 = SP::NewtonEulerDS());

  /** destructor
   */
  virtual ~UniversalJointR() {};

  /** Initialize the joint constants based on the provided base positions.
   * \param q1 A SiconosVector of size 7 indicating translation and
   *           orientation in inertial coordinates.
   * \param q2 An optional SiconosVector of size 7 indicating
   *           translation and orientation; if null, the inertial
   *           frame will be considered as the second base. */
  virtual void setBasePositions(SP::SiconosVector q1,
                                SP::SiconosVector q2 = SP::SiconosVector());

  /* Everything else is inherited from Pivot2, as this joint is
   * exactly just a Pivot2 joint with orthogonal axes. */
};
#endif  //UniversalJointRELATION_H
