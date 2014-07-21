//---------------------------------------------------------------------------//
//                                                                           //
//                                 MGDO                                      //
//                                                                           //
//    This code is the intellectual property of the Majorana and             //
//    GERDA Collaborations.                                                  //
//                                                                           //
//                        *********************                              //
//                                                                           //
//    Neither the authors of this software system, nor their employing       //
//    institutes, nor the agencies providing financial support for this      //
//    work  make  any representation or  warranty, express or implied,       //
//    regarding this software system or assume any liability for its use.    //
//    By copying, distributing or modifying the Program (or any work based   //
//    on on the Program) you indicate your acceptance of this statement,     //
//    and all its terms.                                                     //
//                                                                           //
//---------------------------------------------------------------------------//

#ifndef _MDOUtils_HH
#define _MDOUtils_HH

// Header description:
//
// This file contains utilities and definitions for use throughout the MGDO
// software package.


// Error reporting
#include <iostream>
#define MGDOerr std::cerr << "Error at " << __FILE__ << ":" << __LINE__ << " - " 
#define MGDOwarn std::cout << "Warning at " << __FILE__ << ":" << __LINE__ << " - " 

// External units and constants
#include <float.h> // this is to include DBL_MAX, DBL_MIN, etc
#include <limits.h> // this is to include LONG_MAX, LONG_MIN, etc
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"


// Local additional units and constants - add to CLHEP namespace
namespace CLHEP
{
  // new units
  static const double gigahertz = 1.e+9*hertz;
  static const double minute = 60.0*s;
  static const double hour = 60.0*minute;
  static const double day = 24.0*hour;
  static const double year = 365.2425*day;
  static const double tonne = 1000.0*kg;
  static const double meV = 1.e-3*eV;

  // aliases
  static const double us = microsecond;

  static const double Hz = hertz;
  static const double kHz = kilohertz;
  static const double MHz = megahertz;
  static const double GHz = gigahertz;
}

class MGDOUtils
{
  public:
    static const char* GetMGDOTag();
    static const char* GetMGDORevision();
    static const char* GetMGDODate();
};

#endif
