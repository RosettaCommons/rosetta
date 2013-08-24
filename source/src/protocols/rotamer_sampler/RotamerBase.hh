// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerBase.hh
/// @brief Abstract Base Class for Rotamer generator.
/// @detailed
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_rotamer_sampler_RotamerBase_HH
#define INCLUDED_protocols_rotamer_sampler_RotamerBase_HH

#include <protocols/rotamer_sampler/RotamerBase.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/types.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace rotamer_sampler {

class RotamerBase : public moves::Mover {
public:
	RotamerBase():
		moves::Mover(),
		is_init_( false ),
		random_( true )
	{}

	virtual ~RotamerBase(){}

	/// @brief Initialization
	virtual void init() = 0;

	/// @brief Reset to the first (or random if is_random()) rotamer
	virtual void reset() = 0;

	/// @brief Move to next rotamer
	virtual void operator++() = 0;

	/// @brief Check if there are more rotamers available
	virtual bool not_end() const { return true; }

	/// @brief Check if is random sampling
	virtual bool is_random() const { return random_; }

	/// @brief Set the random sampling state
	virtual void set_random( bool const setting ) { random_ = setting; }

	/// @brief Apply the current rotamer to pose
	virtual void apply( core::pose::Pose & ) = 0;

	/// @brief Name of the class
	virtual std::string get_name() const { return "RotamerBase"; }

protected:
	/// @brief Check if the sampler has been initialized
	virtual bool is_init() const { return is_init_; }

	/// @brief Set the initialization state
	virtual void set_init( bool const setting ) { is_init_ = setting; }

	/// @brief Set value of a parameter, return to pre-init state
	template <class T>
	void set_and_reinit( T & params, T const setting ) {
		params = setting;
		set_init( false );
	}

private:
	bool is_init_, random_;
};

}
}

#endif
