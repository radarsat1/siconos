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
/** \class XMLException
*   \brief This class represent an exeption causing by an XML class of the platform
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.1.3.
*   \date (Creation) 05/25/2004
*
*
*
* XMLException must be throws when an error is find in an XML class
* This exception can be catched by "catch(XMLException)" or "catch(SiconosException)"
*
*
*/

#ifndef __XMLException__
#define __XMLException__

#include "SiconosException.h"
#include <iostream>

// --------------------------------------------------------------------------
class XMLException: public SiconosException
{
public:

  /**
   * \fn XMLException()
   * \brief constructor
   */
  XMLException();

  /**
   * \fn XMLException(const std::string& report)
   * \brief constructor with a report
   * \param std::string report : exception description
   */
  XMLException(const std::string& report);

  /**
   * \fn ~XMLException()
   * \brief destructor
   */
  ~XMLException();

  /**
   * \fn static void selfThrow()
   * \brief static function which throw a XMLException
   * \exception XMLException
   */
  static void selfThrow() ;

  /**
   * \fn static void selfThrow(const std::string& report)
   * \brief static function which throw a XMLException with a report
   * \param std::string report : exception description
   * \exception XMLException
   */
  static void selfThrow(const std::string& report) ;

};

#endif //__XMLException__
