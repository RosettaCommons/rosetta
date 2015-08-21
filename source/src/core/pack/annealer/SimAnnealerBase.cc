// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/SimAnnealerBase.cc
/// @brief  Packer's simulated annealing base class implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/annealer/SimAnnealerBase.hh>


// Package Headers

// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.hh>

// NUmeric headers
#include <numeric/random/random.hh>

// C++ headers
#include <algorithm>

// option key includes

#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <basic/options/option.hh>


using namespace ObjexxFCL;

namespace core {
namespace pack {
namespace annealer {

const core::PackerEnergy SimAnnealerBase::hightemp = 100.0;
const core::PackerEnergy SimAnnealerBase::lowtemp = 0.3;
const core::PackerEnergy SimAnnealerBase::calc_freq_temp = 1.0;
//bool annealing_starts_at_low_temperature = false;

/// @brief constructor
SimAnnealerBase::SimAnnealerBase
(
	int num_rots_to_pack,
	FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	FArray1_int & current_rot_index,
	bool calc_rot_freq,
	FArray1D< core::PackerEnergy > & rot_freq
) :
	num_rots_to_pack_( num_rots_to_pack ),
	bestrotamer_at_seqpos_(bestrotamer_at_seqpos),
	bestenergy_(bestenergy),
	start_with_current_(start_with_current),
	current_rot_index_(current_rot_index),
	calc_rot_freq_(calc_rot_freq),
	rot_freq_(rot_freq),
	outeriterations_( 0 ),
	inneriterations_( 0 ),
	quench_( false ),
	//hightemp_( annealing_starts_at_low_temperature ? 10 : hightemp ),
	hightemp_( hightemp ),
	lowtemp_( lowtemp ),
	temperature_( hightemp_ ),
	jump_(0),
	outeriterations_scaling_(1.),
	inneriterations_scaling_(1.),
	low_temp_annealing_( false ),
	disallow_quench_( false )
{
	bestrotamer_at_seqpos_ = 0;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	outeriterations_scaling_ = option[ OptionKeys::packing::outeriterations_scaling ]();
	inneriterations_scaling_ = option[ OptionKeys::packing::inneriterations_scaling ]();
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
/// virtual destructor
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
////////////////////////////////////////////////////////////////////////////////
SimAnnealerBase:: ~SimAnnealerBase() {}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
/// get the iterations number for simulation
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
////////////////////////////////////////////////////////////////////////////////

int SimAnnealerBase::get_outeriterations() const
{
	return int (outeriterations_scaling_ * outeriterations_);
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
////////////////////////////////////////////////////////////////////////////////

int SimAnnealerBase::get_inneriterations() const
{
	if ( low_temp_annealing_ ) {
		return int (inneriterations_scaling_ *
			(get_temperature() > 3.0f ? inneriterations_ * 0.5 : inneriterations_ ));
	} else {
		return int (inneriterations_scaling_ * inneriterations_);
	}
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
////////////////////////////////////////////////////////////////////////////////
void SimAnnealerBase::setup_iterations()
{
	setup_iterations( num_rots_to_pack_ );
}


void SimAnnealerBase::setup_iterations( const int & num_of_state_changes )
{

	if ( start_with_current_ ) {
		inneriterations_ = std::max( num_of_state_changes, (int) current_rot_index_.size1() );
		outeriterations_ = 10;
	} else {
		if ( low_temp_annealing_ ) {
			inneriterations_ = std::max( 10 * num_of_state_changes, 2 * ((int) current_rot_index_.size1()) );
			outeriterations_ = 10;
		} else {
			inneriterations_ = std::max( 5 * num_of_state_changes, (int) current_rot_index_.size1() );
			outeriterations_ = 20;
		}
	}

	//std::cerr << "SimAnnealerBase::setup iterations(); inner = " << inneriterations_ << " outer: " << outeriterations_ << std::endl;
}

void SimAnnealerBase::setup_temperature( const int & nn )
{
	if ( nn == get_outeriterations() && !disallow_quench_ ) set_to_quench();
	temperature_ = ( hightemp_ - lowtemp_ )*std::exp(-core::PackerEnergy(jump_)) + lowtemp_;
	jump_++;
}

void SimAnnealerBase::setup_temperature(const FArray1D< core::PackerEnergy > & loopenergy,int nn)
{
	bool calc_rot_freq = get_calc_rot_freq();
	core::PackerEnergy avgloopE = 0.0f;
	if ( ( nn == get_outeriterations() )&&(!calc_rot_freq) && !disallow_quench_ ) {
		set_to_quench();
		temperature_ = lowtemp_;
	} else if ( jump_ > 3 ) {
		avgloopE = (loopenergy( nn - 4 ) + loopenergy( nn -3 ) + loopenergy( nn - 2 ))/3.0;
		//std::cout << "avgloopE: " << avgloopE << " vs " << loopenergy( nn - 1) << std::endl;
		if ( ( loopenergy( nn - 1) - avgloopE ) > -1.0 ) {
			//std::cout << "High temp!" << std::endl;
			temperature_ = hightemp_;
			jump_ = 1;
		} else {
			temperature_ = (hightemp_ - lowtemp_)*std::exp(-core::PackerEnergy(jump_)) + lowtemp_;
			jump_++;
			//std::cout << "Cool " << temperature_ << std::endl;
		}
	} else {
		// This is terrible.  If instructions are to be passed to the packer, they are to be passed through the PackerTask.
		// Do not insert code into the middle of a low level function like this such that it can go totally unnoticed
		// by people who rely upon it.
		//else if( pKa_mode::get_pKa_flag() && pKa_mode::pKa_packer_flags::packer_set_temp )
		//{
		// temperature_ = pKa_mode::pKa_packer_flags::packer_temp;
		//}
		temperature_ = (hightemp_ - lowtemp_)*std::exp(-core::PackerEnergy(jump_)) + lowtemp_;
		jump_++;


		if ( calc_rot_freq && (temperature_ < calc_freq_temp) ) {
			temperature_ = calc_freq_temp;
		}
	}
}

// SimAnnealerBase::clear
void
SimAnnealerBase::clear() {
	jump_ = 0;
}

void SimAnnealerBase::set_temperature( core::PackerEnergy new_temp ){ temperature_ = new_temp;}

core::PackerEnergy SimAnnealerBase::get_temperature() const { return temperature_; }
void SimAnnealerBase::set_to_quench(){ quench_ = true;}
void SimAnnealerBase::set_not_to_quench(){ quench_ = false;}
bool SimAnnealerBase::quench() const { return quench_; }
bool SimAnnealerBase::get_start_with_current() const { return start_with_current_; }
bool SimAnnealerBase::get_calc_rot_freq() const { return calc_rot_freq_; }

void SimAnnealerBase::num_rots_to_pack( Size setting ) { num_rots_to_pack_ = setting; }
FArray1D_int& SimAnnealerBase::bestrotamer_at_seqpos() { return bestrotamer_at_seqpos_; }
FArray1D_int const & SimAnnealerBase::bestrotamer_at_seqpos() const { return bestrotamer_at_seqpos_; }
core::PackerEnergy & SimAnnealerBase::bestenergy() { return bestenergy_;}
bool SimAnnealerBase::start_with_current() const { return start_with_current_;}
FArray1_int & SimAnnealerBase::current_rot_index() { return current_rot_index_;}
FArray1_int const & SimAnnealerBase::current_rot_index() const  { return current_rot_index_;}
bool SimAnnealerBase::calc_rot_freq() const { return calc_rot_freq_; }
FArray1D< core::PackerEnergy >& SimAnnealerBase::rot_freq() { return rot_freq_; }
FArray1D< core::PackerEnergy > const & SimAnnealerBase::rot_freq() const { return rot_freq_; }

void SimAnnealerBase::set_hightemp( core::PackerEnergy high) { hightemp_ = high; }
void SimAnnealerBase::set_lowtemp( core::PackerEnergy low) { lowtemp_ = low; }

void SimAnnealerBase::set_disallow_quench( bool const & setting ){ disallow_quench_ = setting; }

////////////////////////////////////////////////////////////////////////////////
///
///
/// @brief
/// accept or reject movement based on Metropolis criterion
/// if this is the first movement, accept by default.
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
////////////////////////////////////////////////////////////////////////////////
bool SimAnnealerBase::pass_metropolis( core::PackerEnergy delta_energy) const
{
	const core::PackerEnergy GOOD_ENERGY = 0.0;
	return pass_metropolis( GOOD_ENERGY, delta_energy );
}

bool SimAnnealerBase::pass_metropolis( core::PackerEnergy previous_energy, core::PackerEnergy delta_energy ) const
{

	// skip evaluating random number generator if quenching
	if ( quench_ ) {
		if ( delta_energy < 0 ) return true;
		else return false;
	}

	core::PackerEnergy lnprob = 0.0f;
	core::PackerEnergy beta = 1.0/temperature_;
	core::PackerEnergy probability = 0.0f;

	/// call this every time, for better numerical stability. Otherwise delta_energies ~ 0 can lead to
	/// instability in the number of calls to the random number generator

	core::PackerEnergy const rg_uniform( numeric::random::rg().uniform() );

	if ( delta_energy < 0 ) {
		return true;
	} else { //evaluate prob of substitution
		lnprob = beta * delta_energy;
		if ( previous_energy > 1.0 ) { // this is specific to FixbbSimAnnealer
			//if both previous energy and new energy are poor
			//increase the probability of accept
			lnprob /= previous_energy;
		}
		if ( lnprob < 10.0 ) {
			probability = std::exp(-lnprob);
			if ( probability > rg_uniform ) return true;
		}
	}
	return false;
}

} // namespace annealer
} // namespace pack
} // namespace core
