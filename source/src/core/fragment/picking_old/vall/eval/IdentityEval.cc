// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/picking_old/vall/eval/IdentityEval.cc
/// @brief  scores a fragment based on secondary structure identity and sequence identity
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <core/fragment/picking_old/vall/eval/IdentityEval.hh>

// project headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>

// numeric headers
#include <numeric/random/random.hh>

#include <core/fragment/picking_old/concepts/Extent.hh>
#include <utility/vector1.hh>



namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace eval {


// static initialization

static thread_local basic::Tracer TR( "core.fragment.picking_old.vall.eval.IdentityEval" );


/// @brief default constructor
IdentityEval::IdentityEval() :
	Super(),
	ss_penalty_( 1.0 ),
	aa_penalty_( 1.0 ),
	randomize_( true )
{}


/// @brief full values constructor
/// @param ss secondary structure string to match against
/// @param aa amino acid structure string to match against
IdentityEval::IdentityEval(
	String const & ss,
	String const & aa,
	Real const ss_penalty,
	Real const aa_penalty,
	bool const randomize
) :
	Super(),
	ss_( ss ),
	aa_( aa ),
	ss_penalty_( ss_penalty ),
	aa_penalty_( aa_penalty ),
	randomize_( randomize )
{
	runtime_assert( ss_.length() == aa_.length() );
}


/// @brief secondary structure constructor
IdentityEval::IdentityEval(
	String const & ss,
	Real const ss_penalty,
	bool const randomize
) :
	Super(),
	ss_( ss ),
	aa_( ss.length(), '.' ),
	ss_penalty_( ss_penalty ),
	aa_penalty_( 0.0 ),
	randomize_( randomize )
{
	runtime_assert( ss_.length() == aa_.length() );
}


/// @brief default copy constructor
IdentityEval::IdentityEval( IdentityEval const & rval ) :
	Super( rval ),
	ss_( rval.ss_ ),
	aa_( rval.aa_ ),
	ss_penalty_( rval.ss_penalty_ ),
	aa_penalty_( rval.aa_penalty_ ),
	randomize_( rval.randomize_ )
{}


/// @brief default destructor
IdentityEval::~IdentityEval() {}


/// @brief copy assignment
IdentityEval & IdentityEval::operator =( IdentityEval const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );
		ss_ = rval.ss_;
		aa_ = rval.aa_;
		ss_penalty_ = rval.ss_penalty_;
		aa_penalty_ = rval.aa_penalty_;
		randomize_ = rval.randomize_;
	}
	return *this;
}


/// @brief clone this object
VallFragmentEvalOP IdentityEval::clone() const {
	return new IdentityEval( *this );
}


/// @brief for a fragment extent, evaluate and store results in a VallFragmentScore
/// @return true, so score is always stored during VallLibrarian::catalog()
bool IdentityEval::eval_impl(
	Extent const & extent,
	VallFragmentScore & fs
)
{
	// no runtime_asserts here, will slow down Librarian operation
	assert( extent.distance() == ss_.length() );
	assert( ss_.length() == aa_.length() );

	Size str_idx = 0;
	for ( VallResidueConstIterator i = extent.begin; i != extent.end; ++i, ++str_idx ) {
		assert( str_idx != ss_.length() );
		assert( str_idx != aa_.length() );

		switch ( ss_.at( str_idx ) ) {
			case 'D': {
				break; // degenerate sec.struct do nothing
			}
			default: {
				if ( i->ss() != ss_.at( str_idx ) ) {
					fs.score += ss_penalty_;
				}
				break;
			}
		}

		switch ( aa_.at( str_idx ) ) {
			case '.': {
				break; // char for degenerate a.a., do nothing
			}
			default: {
				if ( i->aa() != aa_.at( str_idx ) ) {
					fs.score += aa_penalty_;
				}
				break;
			}
		}

	} // foreach residue in extent

	// finalize scores
	if ( randomize_ ) {
		fs.score += ( numeric::random::rg().uniform() * 0.001 );
	}

	return true;
}


/// @brief operation to be perform before catalog() starts
void IdentityEval::pre_catalog_op( VallLibrary const & ) {
	TR << "ss = " << ss_ << "  |  " << "aa = " << aa_ << std::endl;
}


} // namespace eval
} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core
