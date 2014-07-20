// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerSamplerOneValue.hh
/// @brief Generate rotamer for one value. Abstract base class.
/// @author Rhiju Das (rhiju@stanford.edu)


#ifndef INCLUDED_protocols_rotamer_sampler_RotamerSamplerOneValue_HH
#define INCLUDED_protocols_rotamer_sampler_RotamerSamplerOneValue_HH

// Unit headers
#include <protocols/rotamer_sampler/RotamerSamplerOneValue.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/RotamerSamplerSized.hh>

typedef utility::vector1<core::Real> ValueList;

namespace protocols {
namespace rotamer_sampler {

class RotamerSamplerOneValue : public RotamerSamplerSized {
public:

	RotamerSamplerOneValue();

	RotamerSamplerOneValue(
		ValueList const & values
	);

	virtual ~RotamerSamplerOneValue();

	/// @brief Initialization
	virtual void init() {
		set_init( true );
		reset();
	}

	/// @brief Apply the current rotamer to pose
	virtual void apply( core::pose::Pose & ){}

	/// @brief Apply the i-th rotamer to pose
	virtual void apply( core::pose::Pose &, Size const ){}

	/// @brief Get the total number of rotamers in sampler
	virtual core::Size size() const {
		runtime_assert( is_init() );
		return values_.size();
	}

	/// @brief Get the value of current value
	virtual core::Real value() const {
		runtime_assert( is_init() );
		return values_[id_];
	}

	/// @brief Get the value of i-th value
	virtual core::Real value( core::Size const i ) const {
		runtime_assert( is_init() );
		return values_[i];
	}

	/// @brief Set the allowed values in sampler
	virtual void set_values( ValueList const & setting	) {
		set_and_reinit( values_, setting );
	}

	/// @brief Name of the class
	virtual std::string get_name() const { return "RotamerSamplerOneValue"; }

	/// @brief Type of class (see enum in RotamerSamplerTypes.hh)
	virtual RotamerSamplerType type() const { return ONE_VALUE; }

protected: // will be available to derived classes.
	ValueList values_;
};

} //rotamer_sampler
} //protocols

#endif
