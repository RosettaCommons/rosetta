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
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>

#include <utility/string_util.hh>

#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace upstream {

ProteinSCSampler::~ProteinSCSampler() {}

DunbrackSCSampler::DunbrackSCSampler() :
	desymmeterize_( false )
{}

void DunbrackSCSampler::set_desymmeterize( bool setting ) {
	desymmeterize_ = setting;
}

bool DunbrackSCSampler::desymmeterize() const {
	return desymmeterize_;
}

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
	SingleResidueRotamerLibraryCOP res_rotlib = rotlib.get_rsd_library( restype ).lock();

	if ( res_rotlib != 0 ) {

		SingleResidueDunbrackLibraryCOP dun_rotlib(
			utility::pointer::dynamic_pointer_cast< SingleResidueDunbrackLibrary const >
			( res_rotlib ));

		if ( dun_rotlib == 0 ) {
			utility_exit_with_message( "ERROR: Failed to retrieve a Dunbrack rotamer library for AA: " +
				utility::to_string( restype.aa() ) +  " named " +  restype.name() );
		}

		{
		if ( ! dynamic_cast< ProteinBackboneBuildPoint const * > ( & bb_conf ) ) {
			utility_exit_with_message( "ERROR: DunbrackSCSampler expects a ProteinBackboneBuildPoint but"
				"was handed an incompatible type.  ScaffoldBuildPoint #" +
				utility::to_string( bb_conf.index() ) + " is of an incompatible type" );
		}
		}

		ProteinBackboneBuildPoint const & bb(
			static_cast< ProteinBackboneBuildPoint const & >
			( bb_conf ));

		DunbrackRotamerSampleDataVector rotamers = dun_rotlib->get_all_rotamer_samples( bb.phi(), bb.psi() );
		if ( desymmeterize_ && rotamers.size() != 0 ) {
			using namespace core::chemical;
			AA aa = restype.aa();
			if ( aa == aa_asp || aa == aa_glu || aa == aa_phe || aa == aa_tyr ) {
				/// the odd elements will hold the original rotamers, the even elements will hold the shifted by 180 rotamers
				DunbrackRotamerSampleDataVector desymm_rots( rotamers.size() * 2 );

				/// initialize the odd elements
				for ( core::Size ii = 1; ii <= rotamers.size(); ++ii ) {
					desymm_rots[ 2*ii - 1 ] = rotamers[ ii ];
					desymm_rots[ 2*ii - 1  ].set_prob( rotamers[ ii ].probability() * 0.5 );
				}

				/// initialize the even elements
				core::Size const last_chi = rotamers[1].nchi();
				for ( core::Size ii = 1; ii <= rotamers.size(); ++ii ) {
					DunbrackRotamerSampleData iisample = desymm_rots[ 2*ii - 1 ];
					if ( iisample.chi_is_nonrotameric( last_chi )) {
						/// the ProteinUpstreamBuilder assumes that nrchi_lower_boundary < chi_mean < nrchi_upper_boundary,
						/// so just add 180 to everything.
						iisample.set_nrchi_lower_boundary( iisample.nrchi_lower_boundary() + 180 );
						iisample.set_nrchi_upper_boundary( iisample.nrchi_upper_boundary() + 180 );
					}
					iisample.set_chi_mean( last_chi, iisample.chi_mean()[ last_chi ] + 180 );
					desymm_rots[ 2*ii ] = iisample;
				}

				rotamers.swap( desymm_rots );
			}
		}
		return rotamers;
	} else {
		/// No rotamer library.  Return one-element vector
		DunbrackRotamerSampleDataVector one_element_vector( 1 );
		return one_element_vector;
	}

}

}
}
}

