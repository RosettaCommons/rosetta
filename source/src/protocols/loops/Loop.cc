// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/Loop.cc
/// @brief
/// @author Chu Wang
/// @author Mike Tyka
/// @author James Thompson

// Unit Headers
#include <protocols/loops/Loop.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/string.functions.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace loops {

/// @details Auto-generated virtual destructor
Loop::~Loop() {}

static THREAD_LOCAL basic::Tracer tr( "protocols.loops.Loop" );

/// @brief switch DOF_Type for residues in loop. id::CHI, id::BB --- don't use
/// with id::JUMP
void
Loop::switch_movemap(
	core::kinematics::MoveMap& movemap,
	core::id::TorsionType id,
	bool allow_moves
) const {
	for ( Size pos = start(); pos <= stop(); pos++ ) {
		movemap.set( core::kinematics::MoveMap::MoveMapTorsionID( pos, id ), allow_moves );
	}
}


/////////////////////////////////////////////////////////////////////
/// @details  Choose a cutpoint for the loop if one is not specified.
///    Allow any residue to serve as the cutpoint, but prefer those to the center.
void
Loop::choose_cutpoint( core::pose::Pose const & pose ) {
	using core::Size;

	Size const loop_size  ( stop_ - start_ + 1 );
	Size const nres( pose.total_residue() );

	// Special case if we are extending the loop????

	// set up a weight-map for picking cut_s, we want the cut to be internal to
	// the loop, so that the splice-rmsds are calculated between idealized atoms
	// so also only use chainbreak_overlap = 1
	Size const      n_cut_s ( loop_size - 1 );
	ObjexxFCL::FArray1D_float cut_weight( n_cut_s );
	core::Real     total_cut_weight( 0 );
	for ( Size i = 1; i <= n_cut_s; ++i ) {
		// eg: for a 7 rsd loop: 1 2 3 4 3 2 1
		core::Real const weight ( std::max( i, n_cut_s-i + 1 ) );
		total_cut_weight += weight;
		cut_weight(i) = total_cut_weight;
	}
	// choose the cut_:
	cut_ = 0;

	if ( start_ > 1 && stop_ < nres ) {
		//char ss;
		Size nfail( 0 );
		do {
			nfail++;
			core::Real const weight ( numeric::random::uniform()*total_cut_weight );
			for ( Size i = 1; i <= n_cut_s; ++i ) {
				if ( weight <= cut_weight(i) ) {
					cut_ = start_ + i - 1;
					break;
				}
			}
			char ss = pose.secstruct( cut_ );

			// First try to put cutpoint outside of secondary structure
			if ( nfail < 20 ) if ( ( ss == 'H' || ss == 'E' ) ) continue;
			// later only insist that cutpoint is not infront of proline
			if ( nfail < 40 ) if ( pose.residue(cut_+1).aa() == core::chemical::aa_pro ) continue;
			if ( nfail >= 40 ) tr.Error << "Cutpoint choice problem, setting cut_ = " << cut_ << std::endl;
		}while(false);

		if ( cut_ == 0 ) {
			cut_ = ( start_ + n_cut_s/2 );
			tr.Warning << "Cutpoint choice problem, setting cut_ = " << cut_ << std::endl;
		}
		runtime_assert ( cut_ >= start_ && cut_ <= stop_ );
	} else if ( start_ == 1 ) {
		cut_ = 1;
	} else if ( stop_ == nres ) {
		cut_ = nres;
	} else {
		utility_exit_with_message(
			"Somthing is wrong with the input loop: Pose has " +
			ObjexxFCL::string_of( nres ) +
			" residues, but loop definition runs from " +
			ObjexxFCL::string_of( start_ ) + " to " +
			ObjexxFCL::string_of( stop_ )
		);
	}

	tr.Info << "Autoset cut_ for loop " << start_ << " " << stop_
		<< " as " << cut_ << "." << std::endl;
}

void Loop::get_residues( utility::vector1< Size>& selection ) const {
	for ( core::Size i = start_; i<= stop_; i++ ) {
		selection.push_back( i );
	}
}


///////////////////////////////////////////////////////////////////////
/// @details  Detect a terminal loop, logic is more complicated for multi-chain
/// poses. Returns TRUE for terminal loops.
bool
Loop::is_terminal( core::pose::Pose const & pose ) const
{
	if ( start() == 1 )                  return true;     // loop start at first residue
	if ( stop() == pose.total_residue() ) return true;     // loop end at last residue
	if ( !pose.residue( start() -1 ).is_protein() )   return true;  // residue before start is not protein
	if ( !pose.residue( stop()   +1 ).is_protein() )   return true;  // residue after end is not protein
	if ( pose.chain( start() -1 ) != pose.chain( start() ) )  return true; // residues before start is other chain
	if ( pose.chain( stop()   +1 ) != pose.chain( stop() ) )    return true; // residues before start is other chain
	if ( pose.residue( start() ).is_lower_terminus() ) return true; // explicit terminus variant @ start of loop
	if ( pose.residue( stop() ).is_upper_terminus() ) return true; // explicit terminus variant @ end of loop
	return false;
}

/// @details Print the start, stop and cut residues along with the skip rate and the whether or not the loop
/// will be extended.
void
Loop::show( std::ostream & output ) const
{
	output << "LOOP start: " << start_ << "  stop: " << stop_ << "  cut: " << cut_ <<
		"  size: " << size() << "  skip rate: " << skip_rate_ << "  extended?: " <<
		( extended_ ? "True" : "False" ) <<  std::endl;
}

//////////////////////////////////////////////////////////////////////
std::ostream &
operator<<( std::ostream & os, Loop const & loop )
{
	loop.show( os );
	return os;
}

} // namespace loops
} // namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::loops::Loop::save( Archive & arc ) const {
	arc( CEREAL_NVP( start_ ) ); // core::Size
	arc( CEREAL_NVP( stop_ ) ); // core::Size
	arc( CEREAL_NVP( cut_ ) ); // core::Size
	arc( CEREAL_NVP( skip_rate_ ) ); // core::Real
	arc( CEREAL_NVP( extended_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::loops::Loop::load( Archive & arc ) {
	arc( start_ ); // core::Size
	arc( stop_ ); // core::Size
	arc( cut_ ); // core::Size
	arc( skip_rate_ ); // core::Real
	arc( extended_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::loops::Loop );
CEREAL_REGISTER_TYPE( protocols::loops::Loop )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_loops_Loop )
#endif // SERIALIZATION
