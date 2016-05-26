/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
/*! \file FirstOrderLinearDS.hpp

 */
#ifndef FirstOrderLinearDS_H
#define FirstOrderLinearDS_H

#include "FirstOrderNonLinearDS.hpp"

class FirstOrderLinearDS;
namespace boost { namespace serialization {
    template<class Archive>
    inline void load_construct_data(Archive & ar, FirstOrderLinearDS * t, const unsigned int file_version);
  }}

typedef   void (*LDSPtrFunction)(double, unsigned int, double*, unsigned int, double*);


/** First Order Linear Systems - \f$M \dot x = A(t)x(t)+ b(t) + r, x(t_0)=x_0\f$.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 29, 2004
 *
 *
 *  This class represents first order linear systems of the form:
 * \f[
 * M \dot x = A(t)x(t)+ b(t) + r,
 *  x(t_0)=x_0
 * \f]
 * where
 *    - \f$x \in R^{n} \f$ is the state,
 *    - \f$r \in R^{n} \f$  the input due to the Non Smooth Interaction.
 *    - \f$M \in R^{n\times n} \f$ is an optional constant invertible matrix
 *
 *  The right-hand side is described by
 *    - \f$A \in R^{n\times n}\f$
 *    - \f$b \in R^{n} \f$
 *
 * Specific members of this class are \f$A\f$ and \f$b\f$.
 *
 * f is not set for such system and thus calls to computeF or other
 * related functions are forbidden.
 *
 *  Thus, the main steps for FirstOrderLinearDS handling consist in:
 *
 * - Construction: A and b are optional, and can be given as a
 *   matrix/vector or a plug-in.
 * - Initialization: compute values at time=t0 (rhs, jacobianfx, A
 *    ...), usually done when calling simulation->initialize.
 * - Computation at time t, by calling "compute" functions
 *      => computeA
 *      => computeb
 *      => computeRhs, compute \f$ \dot x = M^{-1}(Ax + b + r) \f$
 *
 * Any call to a plug-in requires that it has been set correctly
 * before simulation using one of the following:
 *   => setComputeAFunction
 *   => setComputeBFunction
 *
 **/
class FirstOrderLinearDS : public FirstOrderNonLinearDS
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderLinearDS);


  /** matrix specific to the FirstOrderLinearDS \f$ A \in R^{n \times n}  \f$*/
  SP::SiconosMatrix _A;

  /** FirstOrderLinearDS plug-in to compute A(t,z), id = "A"
  * @param time : current time
  * @param sizeOfA : size of square-matrix A
  * @param[in,out] A : pointer to the first element of A
  * @param size of vector z
  * @param[in,out] z a vector of user-defined parameters
  */
  SP::PluggedObject _pluginA;

  /** FirstOrderLinearDS plug-in to compute b(t,z), id = "b"
  * @param time : current time
  * @param sizeOfB : size of vector b
  * @param[in,out] b : pointer to the first element of b
  * @param size of vector z
  * @param[in,out] param  : a vector of user-defined parameters
  */
  SP::PluggedObject _pluginb;

  /** default constructor
   */
  FirstOrderLinearDS(): FirstOrderNonLinearDS() {};

public:

  /** === CONSTRUCTORS/DESTRUCTOR === */

  /** constructor from a set of data
   *  \param newX0 the initial state of this DynamicalSystem
   *  \param APlugin plugin for A
   *  \param bPlugin plugin for b
   */
  FirstOrderLinearDS(SP::SiconosVector newX0, const std::string& APlugin, const std::string& bPlugin);

  /** constructor from a set of data
   *  \param newX0 the initial state of this DynamicalSystem
   *  \param newA matrix A
   */
  FirstOrderLinearDS(SP::SiconosVector newX0, SP::SiconosMatrix newA);

  /** constructor from the minimum set of data
   *  \param newX0 the initial state of this DynamicalSystem
   */
  FirstOrderLinearDS(SP::SiconosVector newX0);

  /** constructor from a set of data
   *  \param newX0 the initial state of this DynamicalSystem
   *  \param newA matrix A
   *  \param newB b
   */
  FirstOrderLinearDS(SP::SiconosVector newX0, SP::SiconosMatrix newA, SP::SiconosVector newB);

  /** Copy constructor
   * \param FOLDS the original FirstOrderLinearDS we want to copy
   */
  FirstOrderLinearDS(const FirstOrderLinearDS & FOLDS);

  /** destructor */
  virtual ~FirstOrderLinearDS() {};

  /** check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  bool checkDynamicalSystem();

  /** say that the system is linear
   * \return true if the DynamicalSystem is linear.
   */
  virtual bool isLinear()
  {
    return true;
  }

  /** Initialization function for the rhs and its jacobian.
   *  \param time time of initialization.
   */
  void initRhs(double time) ;

  /** Call all plugged-function to initialize plugged-object values
      \param time the time to give to the plugin functions
   */
  virtual void updatePlugins(double time);

  // --- getter and setter ---

  /** get the matrix \f$A\f$
   *  \return pointer (SP) on a matrix
   */
  inline SP::SiconosMatrix A() const
  {
    return _A;
  }

  /** get the matrix \f$A\f$
   *  \return pointer (SP) on a matrix
   */
  virtual SP::SiconosMatrix jacobianfx() const
  {
    return _A;
  };
  /**  function to compute \f$ f: (x,t)\f$
   * \param time time instant used in the computation of \f$f\f$
   */
  virtual void computef(double time);

  /** function to compute \f$ f: (x,t)\f$ with x different from
      current saved state.
   * \param time time instant used in the computation of \f$f\f$
   * \param x2 the state vector used for the computation of \f$f\f$
   */
  virtual void computef(double time, SiconosVector& x2);

  /** set A to pointer newPtr
   *  \param newA the new A matrix
   */
  inline void setAPtr(SP::SiconosMatrix newA)
  {
    _A = newA;
  }

  /** set A to a new matrix
   * \param newA the new A matrix
   **/
  void setA(const SiconosMatrix& newA);

  // --- plugins related functions

  /** set a specified function to compute the matrix A => same action as setComputeJacobianfxFunction
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeAFunction(const std::string& pluginPath, const std::string& functionName);

  /** set a specified function to compute the matrix A
   *  \param fct a pointer on a function
   */
  void setComputeAFunction(LDSPtrFunction fct);

  /** set a specified function to compute the vector b
   *  \param pluginPath the complete path to the plugin file
   *  \param functionName the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputebFunction(const std::string& pluginPath, const std::string& functionName);

  /** set a specified function to compute the vector b
   *  \param fct a pointer on a function
   */
  void setComputebFunction(LDSPtrFunction fct);

  /** default function to compute matrix A => same action as
      computeJacobianfx
      \param time time instant used to compute A
   */
  void computeA(double time);

  /** default function to compute vector b
   * \param time time instant used to compute b
   */
  virtual void computeb(double time);

  /** Default function to the right-hand side term. This is used only by
   * LsodarOSI with EventDriven
   *  \param time current time
   *  \param isDSup flag to avoid recomputation of operators
   *  \warning the \f$z\f$ input is not taken into account when computing the RHS
   */
  void computeRhs(double time, bool isDSup = false);

  /** Default function to jacobian of the right-hand side term
      with respect to x
   *  \param time current time
   *  \param isDSup boolean to avoid recomputation of operators (unused)
   */
  void computeJacobianRhsx(double time, bool isDSup = false);

  /** data display on screen
   */
  void display() const;

  /** overload LagrangianDS corresponding function
   * \return a double, always zero.
   */
  double dsConvergenceIndicator()
  {
    return 1.0;
  }

  /** Get _pluginA
   * \return the plugin for A
   */
  inline SP::PluggedObject getPluginA() const
  {
    return _pluginA;
  };

  /** Get _pluginb
   * \return the plugin for b
   */
  inline SP::PluggedObject getPluginb() const
  {
    return _pluginb;
  };

  /** Set _pluginA
   * \param newPluginA the new plugin
   */
  inline void setPluginA(SP::PluggedObject newPluginA)
  {
    _pluginA = newPluginA;
  };

  /** Set _pluginb
   * \param newPluginB the new plugin
   */
  inline void setPluginB(SP::PluggedObject newPluginB)
  {
    _pluginb = newPluginB;
  };

  /** Reset all the plugins */
  virtual void zeroPlugin();

  friend class boost::serialization::access;
  template<class Archive> inline friend void boost::serialization::load_construct_data(Archive &ar, FirstOrderLinearDS *t, const unsigned int file_version);

  ACCEPT_STD_VISITORS();

};

namespace boost { namespace serialization {
    template<class Archive>
    inline void load_construct_data(
      Archive & ar, FirstOrderLinearDS * t, const unsigned int file_version
      ){
      ::new(t)FirstOrderLinearDS();
      t->_x.resize(0);
      t->_workspace.resize(0);
      t->_workMatrix.resize(0);
    };
  }}

TYPEDEF_SPTR(FirstOrderLinearDS)

#endif // FirstOrderLinearDS_H
