// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/SamplerPlusPlus.hh
/// @brief Abstract Base Class for Sampler generator.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_sampler_plus_plus_SamplerPlusPlus_HH
#define INCLUDED_sampler_plus_plus_SamplerPlusPlus_HH

// Unit headers
#include <protocols/toolbox/SamplerPlusPlus.fwd.hh>
#include <protocols/toolbox/SamplerPlusPlusType.hh>

// Project headers
#include <protocols/moves/Mover.hh>
#include <core/types.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

/////////////////////////////////////////////////////////////////////////////////////////////
//
// @details
//
//    Special kind of sampler that controls DOFs in a pose and is
//     'incrementable' (through operator++). This allows for convenient
//     use in loops for
//
//      (1) enumeration (stepwise assembly),
//      (2) random sampling within defined torsion combinations (stepwise monte carlo), and
//      (3) monte carlo with no pose copies (recces).
//
//    See protocols/stepwise/sampler and protocols/recces/sampler/ for specific use cases.
//
/////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace toolbox {

class SamplerPlusPlus : public moves::Mover {
public:
	typedef utility::vector1<core::Real> TorsionList;

	SamplerPlusPlus():
		moves::Mover(),
		is_init_( false )
	{}

	virtual ~SamplerPlusPlus(){}

	/// @brief Initialization
	virtual void init() = 0;

	/// @brief Reset to the first (or random if random()) rotamer
	virtual void reset() = 0;

	/// @brief Move to next rotamer
	virtual void operator++() = 0;

	/// @brief Apply the current rotamer to pose
	virtual void apply( core::pose::Pose & ) = 0;

	/// @brief Name of the class
	virtual std::string get_name() const = 0;

	/// @brief Type of class (see enum in SamplerPlusPlusTypes.hh)
	virtual SamplerPlusPlusType type() const = 0;

	using moves::Mover::show;
	/// @brief output summary of class
	virtual
	void show( std::ostream & out, Size const indent) const {
		for ( Size n = 1; n <= indent; n++ ) out << ' ';
		out << get_name() << std::endl;
	}

	/// @brief return OP to the subsampler that controls exactly this torsion_id (assume only one).
	/// This is worked out in recces/sampler/MC_Sampler & derived classes -- wouldn't be too hard to
	///  generalize to stepwise/sampler/ classes, too.
	// virtual
	// SamplerPlusPlusOP
	// find( core::id::TorsionID const & torsion_id ) { return 0; } /* = 0; */

protected:
	/// @brief Check if the sampler has been initialized
	bool is_init() const { return is_init_; }

	/// @brief Set the initialization state
	void set_init( bool const setting ) { is_init_ = setting; }

	/// @brief Set value of a parameter, return to pre-init state
	template <class T>
	void set_and_reinit( T & params, T const & setting ) {
		params = setting;
		is_init_ = false;
	}

private:
	bool is_init_;
};

} //toolbox
} //protocols

#endif
