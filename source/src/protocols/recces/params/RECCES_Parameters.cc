// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/params/RECCES_Parameters.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/recces/params/RECCES_Parameters.hh>
#include <protocols/rna/denovo/secstruct_legacy/RNA_SecStructLegacyInfo.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.recces.params.RECCES_Parameters" );

namespace protocols {
namespace recces {
namespace params {

//Constructor
RECCES_Parameters::RECCES_Parameters( core::pose::Pose const & pose )
{
	using namespace rna::denovo::secstruct_legacy;

	if ( has_rna_secstruct_legacy( pose ) ) {
		std::string const & rna_secstruct_legacy( get_rna_secstruct_legacy_from_const_pose( pose ) );
		for ( Size n = 1; n <= pose.total_residue(); n++ ) {
			if ( rna_secstruct_legacy[ n - 1 ] == 'H' ) {
				bp_res_.push_back( n );
			} else {
				runtime_assert( rna_secstruct_legacy[ n - 1 ] == 'X' );
				dangling_res_.push_back( n );
			}
		}
	}

}

//Destructor
RECCES_Parameters::~RECCES_Parameters() = default;

} //params
} //recces
} //protocols
