// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerSized.hh
/// @brief Abstract Base Class for Rotamer sampler with finite size.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_rotamer_sampler_RotamerSized_HH
#define INCLUDED_protocols_rotamer_sampler_RotamerSized_HH

// Unit headers
#include <protocols/rotamer_sampler/RotamerSized.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/RotamerBase.hh>

namespace protocols {
namespace rotamer_sampler {

class RotamerSized: public RotamerBase {

public:
	RotamerSized();

	virtual ~RotamerSized();

	/// @brief Get the total number of rotamers in sampler
	virtual core::Size size() const = 0;

	/// @brief Initialization
	virtual void init();

	/// @brief Reset to the first (or random if random()) rotamer.
	virtual void reset();

	/// @brief Move to next rotamer
	virtual void operator++();

	/// @brief Check if reach the end of rotamer list
	virtual bool not_end() const;

	/// @brief Apply the i-th rotamer to pose
	virtual void apply( core::pose::Pose&, core::Size const ) = 0;

	/// @brief Apply the current rotamer to pose
	void apply( core::pose::Pose & pose ) { apply( pose, id_ ); }

	/// @brief Name of the class
	virtual std::string get_name() const { return "RotamerSized"; }

	/// @brief Type of class (see enum in RotamerTypes.hh)
	virtual RotamerType type() const { return SIZED; }

	Size const & id() const { return id_; }

	virtual void fast_forward(){ id_ = size(); }

	virtual void set_id( Size const setting ){ id_ = setting; }

protected:

	core::Size id_;

};

}
}
#endif
