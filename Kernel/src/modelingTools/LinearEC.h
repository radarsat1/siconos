/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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
#ifndef LINEAREC_H
#define LINEAREC_H

#include "EqualityConstraint.h"
#include "LinearECXML.h"

/** \class LinearEC
 *  \brief Linear Equality Constraint
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.3.
 *  \date 17/01/2005
 *
 *
 */

class LinearEC : public EqualityConstraint
{
public:

  /** \fn LinearEC(void);
   * \brief default constructor
   */
  LinearEC();
  virtual ~LinearEC();

  LinearEC(EqualityConstraintXML*);

  /** \fn void createEqualityConstraint(LagrangianECXML * ecXML)
   *  \brief allows to create the EqualityConstraint with an xml file, or the needed data
   *  \param LagrangianECXML * : the XML object for this EqualityConstraint
   *  \exception RuntimeException
   */
  void createEqualityConstraint(EqualityConstraintXML * ecXML , int number = -1,
                                SiconosMatrix *G = NULL, std::vector<DSInputOutput*> *dsioVector = NULL);
};

#endif // LINEAREC_H

