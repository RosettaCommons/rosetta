// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerSamplerBase.hh
/// @brief Abstract Base Class for RotamerSampler generator.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_rotamer_sampler_RotamerSamplerBase_HH
#define INCLUDED_protocols_rotamer_sampler_RotamerSamplerBase_HH

// Unit headers
#include <protocols/rotamer_sampler/RotamerSamplerBase.fwd.hh>
#include <protocols/rotamer_sampler/RotamerSamplerType.hh>

// Project headers
#include <protocols/moves/Mover.hh>
#include <core/types.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace rotamer_sampler {

class RotamerSamplerBase : public moves::Mover {
public:
	typedef utility::vector1<core::Real> TorsionList;

	RotamerSamplerBase():
		moves::Mover(),
		is_init_( false ),
		random_( true )
	{}

	virtual ~RotamerSamplerBase(){}

	/// @brief Initialization
	virtual void init() = 0;

	/// @brief Reset to the first (or random if random()) rotamer
	virtual void reset() = 0;

	/// @brief Move to next rotamer
	virtual void operator++() = 0;

	/// @brief Check if there are more rotamers available
	virtual bool not_end() const { return true; }

	/// @brief Check if is random sampling
	virtual bool random() const { return random_; }

	/// @brief Set the random sampling state
	virtual void set_random( bool const setting ) { random_ = setting; }

	/// @brief Apply the current rotamer to pose
	virtual void apply( core::pose::Pose & ) = 0;

	/// @brief Name of the class
	virtual std::string get_name() const = 0;

	/// @brief Type of class (see enum in RotamerSamplerTypes.hh)
	virtual RotamerSamplerType type() const = 0;

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
	bool is_init_, random_;
};

} //rotamer_sampler
} //protocols

#endif
