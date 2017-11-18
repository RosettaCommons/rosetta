// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/polar_hydrogens/util.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/polar_hydrogens/util.hh>
#include <protocols/stepwise/modeler/polar_hydrogens/PolarHydrogenPacker.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.modeler.polar_hydrogens.util" );

using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace polar_hydrogens {


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
check_if_proton_chi_atom( pose::Pose const & pose, Size const rsd, Size const atomno ){
	return check_if_proton_chi_atom( pose.residue_type( rsd ), atomno );
}

Size
check_if_proton_chi_atom( chemical::ResidueType const & rsd_type, Size const atomno ) {
	for ( Size n = 1; n <= rsd_type.n_proton_chi(); n++ ) {
		Size chino = rsd_type.proton_chi_2_chi( n );
		Size const & proton_chi_atom = rsd_type.chi_atoms( chino )[4];
		if ( proton_chi_atom == atomno ) return n;
	}
	return 0;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pack_polar_hydrogens( pose::Pose & pose,
	bool allow_virtual_o2prime_hydrogens /* = false */ ) {
	PolarHydrogenPacker polar_hydrogen_packer;
	polar_hydrogen_packer.set_allow_virtual_o2prime_hydrogens( allow_virtual_o2prime_hydrogens );
	polar_hydrogen_packer.apply( pose );
}

} //polar_hydrogens
} //modeler
} //stepwise
} //protocols
