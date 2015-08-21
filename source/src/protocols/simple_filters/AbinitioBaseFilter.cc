// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/filters/Filters.hh
/// @brief header file for Filters.cc
/// @details
///   Contains currently: Filters
///
///
/// @author Robert Vernon

// Unit Headers

// Project Headers
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <basic/options/keys/OptionKeys.hh>


#include <protocols/simple_filters/AbinitioBaseFilter.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>


// Utility headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>


//// C++ headers

static thread_local basic::Tracer tr( "protocols.simple_filters.AbinitioBaseFilter" );

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

namespace protocols {
namespace simple_filters {

AbinitioBaseFilter::AbinitioBaseFilter() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	sstype_ = "";
}

bool
AbinitioBaseFilter::apply( core::pose::Pose const & pose ) const {
	if ( sstype_ == "" ) sstype_=get_protein_sstype( pose );
	tr.Info << "apply filter: " << name() << std::endl;
	return true;
}

std::string
AbinitioBaseFilter::get_protein_sstype( core::pose::Pose const & pose ) const {
	using core::Size;

	utility::vector1< char > secstructs = read_psipred_ss2_file( pose );

	if ( secstructs.size() == 0 ) {
		tr.Error << "Warning: Needs psipred_ss2 to run filters" << std::endl;
		//  disable_all_filters_ = true;
		return "fail";
	}

	// float beta_ratio;
	int alpha = 0;
	// int
	beta_ = 0;
	int helix_length = 0;
	max_helix_length_ = 0;

	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( ! pose.residue(i).is_protein() ) continue;
		if ( secstructs[i] == 'H' ) {
			++alpha;
			++helix_length;
			if ( helix_length > max_helix_length_ ) {
				max_helix_length_ = helix_length;
			}
		} else if ( secstructs[i] == 'E' ) {
			++beta_;
			helix_length = 0;
		} else {
			helix_length = 0;
		}
	}
	max_helix_fraction_ = static_cast< core::Real >( max_helix_length_ ) / pose.total_residue();

	//car protein_ss_type
	std::string protein_sstype;
	beta_ratio_ = static_cast< float >( beta_ ) / ( alpha + beta_ );
	if ( beta_ratio_ >= 0.8 ) {
		tr.Trace << "Protein type: all beta  Fraction beta: " <<
			F( 6, 3, beta_ratio_ ) << std::endl;
		protein_sstype = 'b';
	} else if ( beta_ratio_ > 0.2 && beta_ >= 10 ) {
		tr.Trace << "Protein type: alpha/beta  Fraction beta: " <<
			F( 6, 3, beta_ratio_ ) << std::endl;
		protein_sstype = "ab";
	} else {
		tr.Trace << "Protein type: all alpha  Fraction beta: " <<
			F( 6, 3, beta_ratio_ ) << "nbeta" << I( 6, beta_ ) << std::endl;
		protein_sstype = 'a';
	}

	return protein_sstype;

}

} // filters
} // protocols


