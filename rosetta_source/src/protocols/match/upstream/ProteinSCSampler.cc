// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/upstream/ProteinSCSampler.cc
/// @brief  Class declarations for base class and simple derived class for SC sampling
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/upstream/ProteinSCSampler.hh>

// Package headers
#include <protocols/match/upstream/OriginalScaffoldBuildPoint.hh>

// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>

#include <utility/string_util.hh>

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>


namespace protocols {
namespace match {
namespace upstream {

/// @details If there is no library, returns a one-element vector.
DunbrackSCSampler::DunbrackRotamerSampleDataVector
DunbrackSCSampler::samples(
	ScaffoldBuildPoint const & bb_conf,
	core::chemical::ResidueType const & restype
) const
{
	using namespace core::scoring;
	using namespace core::pack::dunbrack;

	RotamerLibrary const & rotlib( RotamerLibrary::get_instance() );
	SingleResidueRotamerLibraryCAP res_rotlib( rotlib.get_rsd_library( restype ) );

	if ( res_rotlib != 0 ) {

		SingleResidueDunbrackLibraryCAP dun_rotlib(
			dynamic_cast< SingleResidueDunbrackLibrary const * >
			( res_rotlib.get() ));

		if ( dun_rotlib == 0 ) {
			utility_exit_with_message( "ERROR: Failed to retrieve a Dunbrack rotamer library for AA: " +
				utility::to_string( restype.aa() ) +  " named " +  restype.name() );
		}

		{
		ProteinBackboneBuildPoint const * bbptr(
			dynamic_cast< ProteinBackboneBuildPoint const * >
			( & bb_conf ));
		if ( bbptr == 0 ) {
			utility_exit_with_message( "ERROR: DunbrackSCSampler expects a ProteinBackboneBuildPoint but"
				"was handed an incompatible type.  ScaffoldBuildPoint #" +
				utility::to_string( bb_conf.index() ) + " is of an incompatible type" );
		}
		}

		ProteinBackboneBuildPoint const & bb(
			static_cast< ProteinBackboneBuildPoint const & >
			( bb_conf ));

		return dun_rotlib->get_all_rotamer_samples( bb.phi(), bb.psi() );
	} else {
		/// No rotamer library.  Return one-element vector
		DunbrackRotamerSampleDataVector one_element_vector( 1 );
		return one_element_vector;
	}

}

}
}
}

