// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/netcharge_energy/NetChargeEnergySetupSetup.hh
/// @brief Headers for a helper for the EnergyMethod that penalizes deviation from a desired net charge.
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.
/// @author Vikram K. Mulligan (vmullig@uw.edu).



#ifndef INCLUDED_core_scoring_netcharge_energy_NetChargeEnergySetupSetup_hh
#define INCLUDED_core_scoring_netcharge_energy_NetChargeEnergySetupSetup_hh

// Unit headers
#include <core/scoring/netcharge_energy/NetChargeEnergySetup.fwd.hh>

// Package headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueProperty.hh>

// Project headers
#include <core/types.hh>
#include <map>
#include <string>
#include <utility/vector1.hh>
#include <math.h>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace netcharge_energy {

/// @brief The different types of behaviour for the penalty values past the defined range.
/// @details When values are added to this enum, they should also be added to the NetChargeEnergySetup::get_tailfunction_name() function.
enum TailFunction {
	//When adding an effect to this enum, add its name to the get_tailfunction_name() function.
	tf_linear = 1,
	tf_quadratic,
	tf_constant,
	tf_unknown, //Keep this second-to-last
	tf_end_of_list=tf_unknown //Keep this last
};

/// @brief NetChargeEnergySetup, a helper class for the NetChargeEnergy energy method
/// that stores all of its setup data.
class NetChargeEnergySetup : public utility::pointer::ReferenceCount {
public:

	/// @brief Default constructor for NetChargeEnergySetup.
	///
	NetChargeEnergySetup();

	/// @brief Default destructor for NetChargeEnergySetup.
	///
	virtual ~NetChargeEnergySetup();

	/// @brief Clone: create a copy of this object, and return an owning pointer
	/// to the copy.
	virtual NetChargeEnergySetupOP clone() const;

public:

	/**********************
	Public setup functions:
	**********************/

	/// @brief Reset all data in this data storage object.
	///
	void reset();

	/// @brief Initialize from a .charge file.
	///
	void initialize_from_file( std::string const &filename );

	/// @brief Initialize from a string in the format of a .charge file.
	/// @details Allows external code to initialize object without having it read
	/// directly from disk.
	void initialize_from_file_contents( std::string const &filecontents );


public:

	/***********************
	Public lookup functions:
	***********************/

	/// @brief Get tail function name from enum.
	///
	std::string get_tailfunction_name( TailFunction const tf ) const;

	/// @brief Get tail function enum from name.
	/// @details This is slow; it calls get_tailfunction_name repeatedly.  Intended only for use during setup.
	/// Returns tf_unknown if tail function couldn't be parsed.
	TailFunction get_tailfunction_from_name( std::string const &name ) const;

public:

	/*************************
	Public accessor functions:
	*************************/

	/// @brief Get the desired charge.
	inline unsigned long int desired_charge() const { return desired_charge_; }

	/// @brief Look up the penalty given the net charge observed.
	inline core::Real penalty( signed long const observed) const {
		if ( observed < charge_penalties_range_.first ) return out_of_bounds_func( observed, true );
		if ( observed > charge_penalties_range_.second ) return out_of_bounds_func( observed, false );
		return penalties_[ static_cast<core::Size>( observed - charge_penalties_range_.first + 1 ) ];
	}

	/// @brief Get a summary of the data stored in this object
	///
	std::string report() const;

private:

	/******************
	Private functions:
	******************/

	/// @brief Parse out penalty definition from a single block of lines from file.
	///
	void parse_a_penalty_definition( utility::vector1 < std::string > const &lines );

	/// @brief Do some final checks to ensure that data were loaded properly.
	///
	void check_data() const;

	/// @brief Calculate the out-of-bounds behaviour of the penalty function and return the penalty value.
	/// @details Inputs are:
	/// @param[in] observed The observed count.
	/// @param[in] before If true, this is beyond the low range of penalty values.  If false, we're beyond the high range.
	inline core::Real out_of_bounds_func(
		signed long const observed,
		bool const before
	) const {
		TailFunction tf( tf_unknown );
		if ( before ) tf=tailfunctions_.first;
		else tf=tailfunctions_.second;

		switch(tf) {
		case tf_constant :
			return const_out_of_bounds_func(before);
			//break; //cppcheck complains about this below a return.
		case tf_linear :
			return linear_out_of_bounds_func( before, observed );
			//break; //cppcheck complains about this below a return.
		case tf_quadratic :
			return quadratic_out_of_bounds_func( before, observed );
			//break; //cppcheck complains about this below a return.
		default :
			utility_exit_with_message( "Error in core::scoring::netcharge_energy::NetChargeEnergySetup::out_of_bounds_func(): Unknown out of bounds function." );
		}

		return 0.0;
	}

	/// @brief Return a constant value for an out-of-bounds behaviour.
	/// @details Inputs are:
	/// @param[in] before If true, this is beyond the low range of penalty values.  If false, we're beyond the high range.
	inline core::Real const_out_of_bounds_func(
		bool const before
	) const {
		if ( before ) {
			return penalties_[1];
		} else { //if after
			return penalties_[penalties_.size()];
		}
		return 0.0;
	}

	/// @brief Return a linear value for an out-of-bounds behaviour.
	/// @details This fits the first two or last two points to a straight line, then extends it.  Inputs are:
	/// @param[in] before If true, this is beyond the low range of penalty values.  If false, we're beyond the high range.
	/// @param[in] observed The observed count.
	inline core::Real linear_out_of_bounds_func(
		bool const before,
		signed long const observed
	) const {
		core::Real m(0.0); //slope
		core::Real b(0.0); //intercept

		if ( before ) {
			m = (penalties_[2] - penalties_[1]); //denominator is 1
			b = (penalties_[1] - m * static_cast<core::Real>(charge_penalties_range_.first) );
		} else { //if after
			m = (penalties_[ penalties_.size() ] - penalties_[ penalties_.size()-1 ]); //denominator is 1
			b = (penalties_[ penalties_.size() ] - m * static_cast<core::Real>(charge_penalties_range_.second) );
		}
		return m*static_cast<core::Real>( observed )+b;
	}

	/// @brief Return a quadratic value for an out-of-bounds behaviour.
	/// @details This fits the first two or the last two values to a parabola centred on the origin, then extends it.  Inputs are:
	/// @param[in] before If true, this is beyond the low range of penalty values.  If false, we're beyond the high range.
	/// @param[in] observed The observed count.
	inline core::Real quadratic_out_of_bounds_func(
		bool const before,
		signed long const observed
	) const {
		core::Real a(0.0); //coefficient of y=ax^2+b
		core::Real b(0.0); //intercept of y=ax^2+b

		signed long const delta_expected( observed - desired_charge_ );

		if ( before ) {
			signed long const rangefirst_sq( pow( charge_penalties_range_.first - desired_charge_, 2) );
			signed long const rangefirstplusone_sq( pow( charge_penalties_range_.first + 1 - desired_charge_, 2) );
			a = static_cast< core::Real >(penalties_[1] - penalties_[2]) / static_cast<core::Real>( rangefirst_sq - rangefirstplusone_sq );
			b = penalties_[1] - a*rangefirst_sq;
		} else { //if after
			signed long const rangeend_sq( pow( charge_penalties_range_.second - desired_charge_, 2) );
			signed long const rangeendminusone_sq( pow( charge_penalties_range_.second - 1 - desired_charge_, 2) );
			a = static_cast< core::Real >(penalties_[penalties_.size()] - penalties_[penalties_.size()-1]) / static_cast<core::Real>( rangeend_sq - rangeendminusone_sq );
			b = penalties_[penalties_.size()] - a*rangeend_sq;
		}
		return a*pow( static_cast<core::Real>(delta_expected), 2) + b;
	}

private:

	/******************
	Private variables:
	******************/

	/// @brief The desired net charge.
	signed long int desired_charge_;

	/// @brief The lower and upper ends of the charge range specified in the penalties.
	std::pair < signed long int, signed long int > charge_penalties_range_;

	/// @brief Penalties for each residue count (relative to expected).
	utility::vector1 < core::Real > penalties_;

	/// @brief The behaviours at the end of the user-defined range.
	/// @details By default, past the user-defined range, the energy increases quadratically with
	/// increasing deviations from ideal sequence composition.  The user can set this to be linear
	/// or constant, though.
	std::pair <TailFunction, TailFunction > tailfunctions_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // netcharge_energy
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_netcharge_energy_NetChargeEnergySetup )
#endif // SERIALIZATION

#endif // INCLUDED_core_scoring_EtableEnergy_HH
