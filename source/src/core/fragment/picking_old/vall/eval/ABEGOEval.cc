// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/vall/eval/ABEGOEval.cc
/// @brief  scores a fragment based on abego identity
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// unit headers
#include <core/fragment/picking_old/vall/eval/ABEGOEval.hh>
#include <core/sequence/ABEGOManager.hh>

// project headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility>
#include <utility/exit.hh>

// numeric headers
#include <numeric/random/random.hh>

#include <core/fragment/picking_old/concepts/Extent.hh>
#include <utility/vector1.hh>

#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/option.hh>

namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace eval {


// static initialization

static basic::Tracer TR( "core.fragment.picking_old.vall.eval.ABEGOEval" );


/// @brief default constructor
ABEGOEval::ABEGOEval() :
	Super(),
	penalty_( 1.0 ),
	randomize_( true )
{
	abego_.clear();
	am_ = ABEGOManagerOP( new ABEGOManager );
}


/// @brief full values constructor
/// @param abego abego string to match against
ABEGOEval::ABEGOEval(
	utility::vector1< String > const & abego,
	Real const penalty,
	bool const randomize
) :
	Super(),
	abego_( abego ),
	penalty_( penalty ),
	randomize_( randomize )
{
	am_ = ABEGOManagerOP( new ABEGOManager );
}


/// @brief default copy constructor
ABEGOEval::ABEGOEval( ABEGOEval const & /*rval*/ ) = default;


/// @brief default destructor
ABEGOEval::~ABEGOEval() = default;


/// @brief copy assignment
ABEGOEval & ABEGOEval::operator =( ABEGOEval const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );
		abego_ = rval.abego_;
		penalty_ = rval.penalty_;
		randomize_ = rval.randomize_;
		am_ = rval.am_;
	}
	return *this;
}


/// @brief clone this object
VallFragmentEvalOP ABEGOEval::clone() const {
	return VallFragmentEvalOP( new ABEGOEval( *this ) );
}


/// @brief for a fragment extent, evaluate and store results in a VallFragmentScore
/// @return true, so score is always stored during VallLibrarian::catalog()
bool ABEGOEval::eval_impl(
	Extent const & extent,
	VallFragmentScore & fs
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// no runtime_asserts here, will slow down Librarian operation
	debug_assert( extent.distance() == abego_.size() );

	Size pos = 0;
	Real default_A_phi = -63.8;
	Real default_A_psi = -41.1;

	for ( VallResidueConstIterator i = extent.begin; i != extent.end; ++i, ++pos ) {
		String::size_type index = abego_[ pos+1 ].find( "X" );
		if ( index == String::npos ) {
			bool flag( false );
			for ( char ii : abego_[ pos+1 ] ) {
				if ( am_->check_rama( ii, i->phi(), i->psi(), i->omega() ) ) {
					if ( (ii != 'A') || (!(option[OptionKeys::frags::ABEGO::phi_psi_range_A].user())) ) {
						flag = true;
						break;
					} else {
						Real range = option[OptionKeys::frags::ABEGO::phi_psi_range_A]();
						if ( (i->phi() > default_A_phi-range) && (i->phi() < default_A_phi+range) && (i->psi() > default_A_psi - range) && (i->psi() < default_A_psi + range) ) {
							flag = true;
							break;
						}
					}
				}
				if ( ii == 'D' ) {
					flag = true;
					break;
				}
			}
			if ( !flag ) fs.score += penalty_;
		}
	} // foreach residue in extent

	// finalize scores
	if ( randomize_ ) {
		fs.score += ( numeric::random::rg().uniform() * 0.001 );
	}

	return true;
}


/// @brief operation to be perform before catalog() starts
void ABEGOEval::pre_catalog_op( VallLibrary const & ) {

	std::ostringstream abego;
	for ( Size ii=1; ii<=abego_.size(); ++ii ) {
		Size length = abego_[ ii ].length();
		if ( length > 1 ) {
			std::ostringstream multi;
			multi << "[";
			for ( char jj : abego_[ ii ] ) {
				multi << jj;
			}
			multi << "]";
			abego << multi.str();
		} else {
			abego << abego_[ ii ].at( 0 );
		}
	}
	TR << "abego = " << abego.str() << std::endl;
}


} // namespace eval
} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core
