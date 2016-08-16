/*                                                                            */
/*                           ----  SPARTA   ----                              */
/*     Shifts Prediction from Analogue of Residue type and Torsion Angle      */
/*                           Yang Shen and Ad Bax                             */
/*                    J. Biomol. NMR, 38, 289-302 (2007)                      */
/*                 NIH, NIDDK, Laboratory of Chemical Physics                 */
/*                     version, 1.2 (build 2009.0928.17)                     */
/*                                                                            */
/*                      for any problem, please contact                       */
/*                          shenyang@niddk.nih.gov                            */
/*                                                                            */
/******************************************************************************/
///  modified for use inside CS-Rosetta  by Oliver Lange
///
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief Set of functions from Sparta class that don't actually manipulate
/// any state in the class and have thus been removed.
/// @author Oliver Lange
/// @author James Thompson

#ifndef protocols_sparta_SpartaUtil_hh
#define protocols_sparta_SpartaUtil_hh

#include <boost/unordered_map.hpp>

#include <protocols/sparta/GDB.hh>


#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>
#include <protocols/sparta/Sparta.hh>

namespace protocols {
namespace sparta {

void
calc_per_residue_scores(
	Sparta::SpartaLib::AtomNameList& atom_names,
	GDB & Pred_Sum,
	GDB & REF_CS_Tab,
	GDB & COMP_Tab,
	utility::vector1< float > & per_residue_scores
);

core::Real compareRef_fxn(
	Sparta::SpartaLib::AtomNameList&,
	GDB & Pred_Sum,
	GDB & REF_CS_Tab,
	GDB & COMP_Tab // pass-by-reference, will be obliterated
);

float getDiff( float ang1, float ang2 ); // calculate the different between two angles
float getAVG( utility::vector0< float > &v1 );
float getSTD( utility::vector0< float > &v1 );
float getRMS( utility::vector0< float > &v1, utility::vector0< float > &v2 );

int MKDIR(const char *dirName);
bool isDirExists(const std::string &Dir) ;

char * ftoa( float n, char *buff, char f='g', int prec=6 );
char * itoa( int n, char *buff, int base=10 );

} // sparta
} // protocols

#endif
