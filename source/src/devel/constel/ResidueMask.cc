// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief Definition of class ResidueMask
/// @author Andrea Bazzoli (bazzoli@ku.edu)

#include "ResidueMask.hh"
#include "Primitives.hh"
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
#include <fstream>

static basic::Tracer TR( "devel.constel.ResidueMask");

namespace devel {
namespace constel {


/// @brief: reads a residue mask from file
///
/// @param[in] ps pose over whose residues the mask is defined
/// @param[in] fname path to the input file
///
/// @details The input file lists the residues whose bit must be set to true.
///  The format of the file is the following:
///  I1 C1\n
///  ...
///  IN CN\n ,
///  where Ii and Ci are the residue index and the chain ID, respectively, of
///  the ith residue whose bit must be set to true (i=1,...,N).
///
ResidueMask::ResidueMask(Pose& ps, std::string const &fname) :
	mask(ps.size(), false) {

	std::ifstream ifs(fname.c_str());
	if ( !ifs ) {
		TR << "can't open " << fname << std::endl;
		throw utility::excn::EXCN_BadInput(fname);
	}

	int ri;
	char rc;
	while ( ifs >> ri >> rc ) {
		if ( rc == '_' ) {
			rc = ' ';
		}
		mask[get_pose_resnum(ri, rc, ps)] = true;
	}
}


/// @brief: prints the residue mask as a binary string
///
/// @param[out] os output stream
///
/// @details: the ith output digit equals 1 if the ith element in the mask is true;
///           it equals 0 if the ith element is false (i=1,...,SIZ).
///
void ResidueMask::print(std::ostream& os) const {

	Size const SIZ = mask.size();
	for ( Size i=1; i<=SIZ; ++i ) {
		os << (mask[i] ? "1" : "0");
	}
}


/// @brief: the mask's subscription operator
///
/// @param[in] i index in mask
///
/// @return value of the mask at that index
///
bool ResidueMask::operator[] (Size const i) {return mask[i];}


} // constel
} // devel

