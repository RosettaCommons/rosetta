// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/loops/loop_closure/ccd/RamaCheck.cc
/// @brief  Implementation for RamaCheck classes (RamaCheckBase, RamaCheck1B, RamaCheck2B)
/// @author Brian D. Weitzner


// Unit Headers
#include <protocols/loops/loop_closure/ccd/RamaCheck.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/Ramachandran2B.hh>
#include <core/scoring/ScoringManager.hh>

// Utility Headers
#include <utility/assert.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

// Basic Headers
#include <basic/Tracer.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/option.hh>

// Numeric Headers
#include <numeric/random/random.hh>

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

static THREAD_LOCAL basic::Tracer TR( "protocols.loops.loop_closure.ccd.RamaCheck" );

// constants
core::Real const RamaCheckBase::BAD_SCORE = 500.0;

// Empty constructor
/// @details  Initialize the ReferenceCount and the private data
RamaCheckBase::RamaCheckBase() : Parent()
{
	init();
}

// Copy constructor
RamaCheckBase::RamaCheckBase( RamaCheckBase const & object_to_copy ) : Parent( object_to_copy )
{
	copy_data( * this, object_to_copy );
}

// Assignment operator
RamaCheckBase &
RamaCheckBase::operator=( RamaCheckBase const & object_to_copy )
{
	// Check for self-assignment.
	if ( this != &object_to_copy ) {
		Parent::operator=( object_to_copy );
		copy_data( * this, object_to_copy );
	}

	return * this;
}

// Destructor
RamaCheckBase::~RamaCheckBase() {}

/// @details Print the RamaCheck's basic information and details about the instance's configuration.
void RamaCheckBase::show( std::ostream & output ) const
{
	using std::endl;

	output << name() << ":" << endl;
	output <<
		"\nTemperature:                                   " << temperature_ <<
		"\nMaximum Ramachandran Score Increase Tolerated: " << max_rama_score_increase_ << endl;
}

/// @details Enable the RamaCheck classes to parse tags with options for "max_rama_score_increase" and "temperature".
void RamaCheckBase::parse_my_tag( utility::tag::TagCOP tag )
{
	using core::Real;

	// Ok, let's parse this shit.
	if ( tag->hasOption( "max_rama_score_increase" ) ) {
		max_rama_score_increase( tag->getOption< Real >( "max_rama_score_increase" ) );
	}

	if ( tag->hasOption( "temperature" ) ) {
		temperature( tag->getOption< Real >( "temperature" ) );
	}
}

/// @details Associate the RamaCheck classes with the -loops:ccd:max_rama_score_increase and -loops:ccd:temperature
/// options
void RamaCheckBase::register_options()
{
	using namespace basic::options;

	option.add_relevant( OptionKeys::loops::ccd::max_rama_score_increase );
	option.add_relevant( OptionKeys::loops::ccd::temperature );
}

void RamaCheckBase::initialize_starting_rama_scores( core::pose::Pose const & pose ) const
{
	if ( TR.Warning.visible() ) {
		TR.Trace << "Initializing starting_rama_scores with current conformation." << std::endl;
	}
	starting_rama_scores_ = RamaScoreVector( pose.size(), BAD_SCORE ); // init with high value
	for ( core::uint seqpos = 1; seqpos <= pose.size(); ++seqpos ) {
		if ( ! pose.residue( seqpos ).is_protein() ) { continue; }
		starting_rama_scores_[ seqpos ] = compute_rama_score( pose, seqpos, pose.phi( seqpos ), pose.psi( seqpos ) );
	}
}

/// @details New conformations are accepted using the Metropolis criterion. More specifically, if the new conformation
/// has an improved score, it is always accepted. If the score gets worse (increases), it is accepted with a penalty
/// based on the difference of energy between the two states. Even if a score is accepted with a penalty, it may be
/// rejected if it exceeds the maximum score increase cutoff, which can be accessed by max_rama_score_increase().
bool
RamaCheckBase::accept_new_conformation(
	core::pose::Pose const & pose,
	core::id::TorsionID const & torsion_id,
	core::Angle const alpha
) const

{
	using utility::vector1;
	using core::Angle;
	using core::Real;
	using core::id::phi_torsion;
	using core::id::psi_torsion;

	// Rama only depends on phi & psi.
	if ( torsion_id.torsion() != phi_torsion && torsion_id.torsion() != psi_torsion ) { return true; }

	core::uint const seqpos( torsion_id.rsd() );
	// Rama only makes sense for protein residues.
	if ( ! pose.residue( seqpos ).is_protein() ) { return true; }

	// Make sure starting_rama_scores_ has been initialized. If not, initialize it with the current conformation.
	if ( starting_rama_scores_.size() != pose.size() ) {
		if ( TR.Warning.visible() ) {
			TR.Warning << "Starting Rama scores have not been initialized!" << std::endl;
		}
		initialize_starting_rama_scores( pose );
	}

	vector1< Angle > delta_phi_psi( 2, 0.0 );
	delta_phi_psi[ torsion_id.torsion() ] = alpha;

	// Evaluate the rama score difference.
	Real const current_rama_score( compute_rama_score(pose, seqpos, pose.phi( seqpos ), pose.psi( seqpos ) ) );
	Real const proposed_rama_score( compute_rama_score(pose, seqpos,
		pose.phi( seqpos ) + delta_phi_psi[ phi_torsion ],
		pose.psi( seqpos ) + delta_phi_psi[ psi_torsion ] ) );

	if ( TR.Debug.visible() ) {
		TR.Debug  << "Current rama score at position " << seqpos << ": " << current_rama_score <<
			"; Proposed score: " << proposed_rama_score << std::endl;
	}

	return accept_new_conformation( current_rama_score, proposed_rama_score, starting_rama_scores_[ seqpos ]);
}

/// @details The total net change in Ramachandran score is simply the sum of the differences of the current Ramachandran
/// score and the initial Ramachandran score of each residue in the range <first_res>, <last_res>.
/// @note <first_res> must be a smaller value than <last_res> and all values in [<first_res>, <last_res>] must
/// correspond to valid residue numbers in <pose>.
core::Real
RamaCheckBase::total_net_change_in_rama_score_over_range(
	core::pose::Pose const & pose,
	core::uint const first_res,
	core::uint const last_res
) const
{
	using core::Real;

	debug_assert( first_res < last_res );
	debug_assert( last_res <= pose.size() );

	Real total_change_in_rama_score( 0.0 );

	// Make sure starting_rama_scores_ has been initialized. If not, initialize it with the current conformation.
	// This will result in this metric evaluating to exactly 0.
	if ( starting_rama_scores_.size() != pose.size() ) {
		if ( TR.Debug.visible() ) {
			TR.Debug << "Starting Rama scores have not been initialized!" << std::endl;
		}
		initialize_starting_rama_scores( pose );
	}

	for ( core::uint i = first_res; i <= last_res; ++i ) {
		Real const final_rama_score( compute_rama_score( pose, i, pose.phi( i ), pose.psi( i ) ) );
		total_change_in_rama_score += ( final_rama_score - starting_rama_scores_[ i ] );
	}
	return total_change_in_rama_score;
}

/// @details The average change in Ramachandran score is the total net change in Ramachandran score divided by the
/// number of residues in the range [<first_res>, <last_res>].
/// @note <first_res> must be a smaller value than <last_res> and all values in [<first_res>, <last_res>] must
/// correspond to valid residue numbers in <pose>.
core::Real
RamaCheckBase::average_change_in_rama_score_over_range(
	core::pose::Pose const & pose,
	core::uint first_res,
	core::uint last_res
) const
{
	using core::Real;

	debug_assert( first_res < last_res );
	Real length( last_res - first_res + 1 );

	return total_net_change_in_rama_score_over_range( pose, first_res, last_res ) / length;
}

void RamaCheckBase::copy_data( RamaCheckBase & to, RamaCheckBase const & from ) const
{
	to.temperature_ = from.temperature_;
	to.max_rama_score_increase_ = from.max_rama_score_increase_;
	to.starting_rama_scores_ = from.starting_rama_scores_;
}

void RamaCheckBase::init()
{
	temperature_ = 0.25;   // Should this be higher?
	max_rama_score_increase_ = 2.0;
	starting_rama_scores_ = RamaScoreVector();
}

void RamaCheckBase::init_options()
{
	using namespace basic::options;

	if ( option[ OptionKeys::loops::ccd::max_rama_score_increase ].user() ) {
		max_rama_score_increase( option[ OptionKeys::loops::ccd::max_rama_score_increase ]() );
	}

	if ( option[ OptionKeys::loops::ccd::temperature ].user() ) {
		temperature( option[ OptionKeys::loops::ccd::temperature ]() );
	}
}

// Apply the Metroplis criterion, but reject if the score exceeds the cutoff
bool
RamaCheckBase::accept_new_conformation(
	core::Real const current_score,
	core::Real const proposed_score,
	core::Real const starting_score
) const
{
	using std::max;
	using numeric::random::uniform;
	using core::Real;

	if ( proposed_score > current_score ) {
		Real const boltz_factor( ( current_score - proposed_score ) / temperature_ );
		Real const minimum_boltz_factor( -40.0f );
		Real const probability ( exp( max( minimum_boltz_factor, boltz_factor ) ) );

		if ( ( proposed_score - starting_score ) > max_rama_score_increase_ || uniform() >= probability ) {
			return false;
		}
	}

	return true;
}

// End of RamaCheckBase //

// RamaCheck1B //
// Empty constructor
/// @details Initialize RamaCheck1B with a pointer to the Ramachandran potential.
RamaCheck1B::RamaCheck1B() : RamaCheckBase(),
	rama_( core::scoring::ScoringManager::get_instance()->get_Ramachandran_ptr() )
{}

// Copy constructor
RamaCheck1B::RamaCheck1B( RamaCheck1B const & object_to_copy ) : RamaCheckBase( object_to_copy )
{
	copy_data( * this, object_to_copy );
}

// Assignment operator
RamaCheck1B &
RamaCheck1B::operator=( RamaCheck1B const & object_to_copy )
{
	// Check for self-assignment.
	if ( this != & object_to_copy ) {
		RamaCheckBase::operator=( object_to_copy );
		copy_data( *this, object_to_copy );
	}

	return * this;
}

// Destructor
RamaCheck1B::~RamaCheck1B() {}

std::string RamaCheck1B::name() const
{
	return "RamaCheck1B";
}

/// @details Return a fully configured copy of this instance upcasted to a RamaCheckOP
RamaCheckBaseOP
RamaCheck1B::clone() const
{
	return RamaCheckBaseOP( new RamaCheck1B( * this ) );
}

/// @details The neighbor-independent Ramachandran maps are used to compute the score for a residue in the conformation
/// defined by <phi> and <psi>. These values represent the conformation preferences of the residue in its average
/// sequence context.
core::Real
RamaCheck1B::compute_rama_score(
	core::pose::Pose const & pose,
	core::uint const seqpos,
	core::Real const phi,
	core::Real const psi
) const
{
	// if the AA is aa_unk, then try using backbone_aa. This is useful for ncaa.
	core::chemical::AA currAA;

	if ( pose.aa(seqpos) == core::chemical::aa_unk ) {
		currAA = pose.residue(seqpos).backbone_aa();
	} else {
		currAA = pose.aa(seqpos);
	}

	return rama_->eval_rama_score_residue( currAA, phi, psi );
}

void RamaCheck1B::copy_data( RamaCheck1B & to, RamaCheck1B const & from ) const
{
	// This is a shallow copy because we only want to point to the singe shared instance of the potential
	to.rama_ = from.rama_;
}

// End of RamaCheck1B //

// RamaCheck2B //
// Empty constructor
/// @details Initialize RamaCheck2B with a pointer to the Ramachandran2B potential.
RamaCheck2B::RamaCheck2B() : RamaCheckBase(),
	rama_( core::scoring::ScoringManager::get_instance()->get_Ramachandran2B_ptr() )
{}

// Copy constructor
RamaCheck2B::RamaCheck2B( RamaCheck2B const & object_to_copy ) : RamaCheckBase( object_to_copy )
{
	copy_data( * this, object_to_copy );
}

// Assignment operator
RamaCheck2B &
RamaCheck2B::operator=( RamaCheck2B const & object_to_copy )
{
	// Check for self-assignment.
	if ( this != & object_to_copy ) {
		RamaCheckBase::operator=( object_to_copy );
		copy_data( *this, object_to_copy );
	}

	return * this;
}

// Destructor
RamaCheck2B::~RamaCheck2B() {}

std::string RamaCheck2B::name() const
{
	return "RamaCheck2B";
}

/// @details Return a fully configured copy of this instance upcasted to a RamaCheckOP
RamaCheckBaseOP
RamaCheck2B::clone() const
{
	return RamaCheckBaseOP( new RamaCheck2B( * this ) );
}

/// @details The neighbor-dependent Ramachandran maps are used to compute the score for a residue in the conformation
/// defined by <phi> and <psi>. The neighboring residue identities are taken from the pose, and these three residues
/// are used to compute the score.  In cases where there are no neighboring residues (i.e. the first and last residues
/// in the pose), a single neighbor-dependent score is returned. This is a little strange because these reisudes are
/// actually missing a torsion angle (residue 1 is missing phi, residue Nres is missing psi). Accordingly, the current
/// behavior of Ramachandran2B returns a score of zero in these circumstances.
core::Real
RamaCheck2B::compute_rama_score(
	core::pose::Pose const & pose,
	core::uint const seqpos,
	core::Real const phi,
	core::Real const psi
) const
{
	// If we're looking at the first residue, we cannot compute the energy of the "triple" so we'll
	// return the neighbor-dependent ramachandran score considering only the downstream (C-terminal,
	// right, higher res no -- however you want to put it) neighbor
	if ( seqpos == 1 ) {
		return rama_->RamaE_Upper( phi, psi, pose.aa( seqpos ), pose.aa( seqpos + 1 ) );
	}

	// If we're looking at the last residue, we cannot compute the energy of the "triple" so we'll
	// return the neighbor-dependent ramachandran score considering only the upstream (N-terminal,
	// left, lower res no -- however you want to put it) neighbor
	if ( seqpos == pose.size() ) {
		return rama_->RamaE_Lower( phi, psi, pose.aa( seqpos ), pose.aa( seqpos - 1 ) );
	}

	// Compute the score of the "triple" - that is, consider the conformation of the current residue
	// and the identity of the residues to which it is bound. Since we're working in score space, we
	// can simply add the scores of the upper and lower neighbors and subtract away the score of the
	// current conformation of the central residue averaged over all neighbor types.
	return rama_->eval_rama_score_residue( phi, psi, pose.aa( seqpos ), pose.aa( seqpos - 1 ), pose.aa( seqpos + 1 ) );
}

void RamaCheck2B::copy_data( RamaCheck2B & to, RamaCheck2B const & from ) const
{
	// This is a shallow copy because we only want to point to the singe shared instance of the potential
	to.rama_ = from.rama_;
}

// End of RamaCheck2B //

// Helper methods /////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that RamaCheckBase can be "printed" in PyRosetta).
std::ostream &
operator<< ( std::ostream & os, RamaCheckBase const & rama_check )
{
	rama_check.show( os );
	return os;
}

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols
