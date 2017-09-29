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

/*! \file Pivot2JointR.hpp
*/

#ifndef Pivot2JointRELATION_H
#define Pivot2JointRELATION_H

#include <MechanicsFwd.hpp>
#include <SiconosFwd.hpp>
#include <NewtonEulerJointR.hpp>

/** \class Pivot2JointR
 *  \brief This class implements a pivot2 joint between one or
 *  two Newton/Euler Dynamical system.  It is similar to a
 *  PrismaticJointR but allows for rotation around the axis.
 *
 * From a given axis, we construct two unit othorgonal vectors to the
 *  axis V1 and V2 such that (axis,V1,V2) is an orthogonal frame
 */
class Pivot2JointR : public NewtonEulerJointR
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(Pivot2JointR);

  /** _V[0] and _V[1] each define two planes as a set of two unit
   * vectors orthogonal to the axes of rotation.
   */
  SP::SiconosVector _V[2][2];

  SP::SiconosVector _V0;

  double _cq2q101;
  double _cq2q102;
  double _cq2q103;
  double _cq2q104;

  /** P is the point defining the location of the line created by
   * _axis0.  It is stored in the q1 frame, i.e. the vector from
   * initial G1 to P, called _G1P0. */
  SP::SiconosVector _G1P0;

  /** _G2P0 is the vector from initial G2 to P0 */
  SP::SiconosVector _G2P0;

  /** _G1P1q1 is the vector from initial G1 to P1 in q1 frame */
  SP::SiconosVector _G1P1q1;

  /** _G2P2q2 is the vector from initial G2 to P2 in q2 frame */
  SP::SiconosVector _G2P2q2;

  /** Cumulative number of twists around the joint relative to initial
   * angular difference. */
  int _twistCount;    // TODO: Should be in a graph work vector?
  double _previousAngle; // Needed to track _twistCount, TODO: work vector?
  double _initialAngle;

  /** If true, the first axis is a free prismatic axis, i.e. a
   * cylindrical joint rather than a pivot.  Note: Changes
   * numberOfDoF()/numberOfConstraints()! */
  bool _suspensionfree;

public:

  /** Empty constructor. The relation may be initialized later by
   * setPoint, setAxis, setAbsolute, and setBasePositions. */
  Pivot2JointR();

  /** Constructor based on one or two dynamical systems, a point and an axis.
   *  \param d1 first DynamicalSystem linked by the joint.
   *  \param d2 second DynamicalSystem linked by the joint, or NULL
   *            for absolute frame.
   *  \param P1 SiconosVector of size 3 that defines the point around
   *            which rotation is allowed.
   *  \param A1 SiconosVector of size 3 that defines the first pivot axis.
   *  \param A2 SiconosVector of size 3 that defines the second pivot axis.
   *  \param suspensionFree True if A1 has prismatic freedom.
   *  \param absoluteRef if true, P and A are in the absolute frame,
   *                     otherwise P and A are in d1 frame.
   */
  Pivot2JointR(SP::SiconosVector P, SP::SiconosVector A1, SP::SiconosVector A2,
               bool suspensionFree, bool absoluteRef,
               SP::NewtonEulerDS d1 = SP::NewtonEulerDS(),
               SP::NewtonEulerDS d2 = SP::NewtonEulerDS());

  /** Initialize the joint constants based on the provided base positions.
   * \param q1 A SiconosVector of size 7 indicating translation and
   *           orientation in inertial coordinates.
   * \param q2 An optional SiconosVector of size 7 indicating
   *           translation and orientation; if null, the inertial
   *           frame will be considered as the second base. */
  virtual void setBasePositions(SP::SiconosVector q1,
                                SP::SiconosVector q2 = SP::SiconosVector());

  /** Set whether or not the suspension axis is free.  If true, the
   * first axis is a free prismatic DoF, like a cylindrical joint
   * rather than a pivot joint. Important: Changes
   * numberOfDoF()/numberOfConstraints(), and therefore the associated
   * NonSmoothLaw size must be changed to correspond! */
  virtual void setSuspensionFree(bool suspensionfree)
    { _suspensionfree = suspensionfree; }

  /** Get whether or not the suspension axis (first axis) is
   * prismatically free, see setSuspensionFree(). */
  virtual bool suspensionFree() { return _suspensionfree; }

  void computeRotationPlanes();

  /** destructor
   */
  virtual ~Pivot2JointR() {};

  virtual void computeJachq(double time, Interaction& inter, SP::BlockVector q0 );

  virtual void computeh(double time, BlockVector& q0, SiconosVector& y);

  /** Compute the vector of linear and angular positions of the free axes */
  virtual void computehDoF(double time, BlockVector& q0, SiconosVector& y,
                           unsigned int axis);

  /** Compute the jacobian of linear and angular DoF with respect to some q */
  virtual void computeJachqDoF(double time, Interaction& inter,
                               SP::BlockVector q0, SimpleMatrix& jachq,
                               unsigned int axis);

  void Jd1d2(
    double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
    double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  void Jd1(
    double X1, double Y1, double Z1, double q10, double q11, double q12, double q13);

  /** Get the number of constraints defined in the joint
      \return the number of constraints
   */
  virtual unsigned int numberOfConstraints();// { return 3 /*- _suspensionfree*1*/; }

  /** Return the number of degrees of freedom of this joint.
      \return the number of degrees of freedom (DoF)
   */
  virtual unsigned int numberOfDoF() { return 4 /* + _suspensionfree*1*/; }

  /** Return the type of a degree of freedom of this joint.
      \return the type of the degree of freedom (DoF)
  */
  virtual DoF_Type typeOfDoF(unsigned int axis) {
    if (axis==0) return DOF_TYPE_LINEAR;
    else if (axis==1) return DOF_TYPE_ANGULAR;
    else return DOF_TYPE_INVALID;
  }
};
#endif  //Pivot2JointRELATION_H
