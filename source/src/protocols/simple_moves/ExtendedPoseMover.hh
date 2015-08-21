// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/ExtendedPoseMover.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_protocols_simple_moves_ExtendedPoseMover_HH
#define INCLUDED_protocols_simple_moves_ExtendedPoseMover_HH

// Unit header
#include <protocols/simple_moves/ExtendedPoseMover.fwd.hh>

// C/C++ headers
#include <string>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Package headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

class ExtendedPoseMover : public protocols::moves::Mover {
	typedef std::string string;

public:
	ExtendedPoseMover(const string& sequence = "",
		const string& residue_type_set_ = "centroid");

	/// @brief Creates an extended, idealized pose from the sequence and residue
	/// type set specified in the constructor.
	void apply(core::pose::Pose& pose);

	/// @brief Returns the name of this mover
	string get_name() const;

	/// @brief Returns true if this instance is valid (i.e. contains a sequence).
	/// Since RosettaScripts mandates the existence of a no-argument constructor,
	/// we lose the ability to reason about the validity of a particular instance.
	bool valid() const;

	// -- Accessors -- //
	/// @brief Returns the sequence
	const string& sequence() const;

	/// @brief Returns the residue type set
	const string& residue_type_set() const;

	// -- Mutators -- //
	/// @brief Updates the sequence to be used in calls to apply()
	void sequence(const string& sequence);

	/// @brief Updates the residue type set to be used in calls to apply()
	void residue_type_set(const string& residue_type_set);

	// -- RosettaScripts -- //
	/// @brief Creates a replica of this protocols::moves::Mover
	protocols::moves::MoverOP clone() const;

	/// @brief Creates a new instance by calling the no-argument constructor
	protocols::moves::MoverOP fresh_instance() const;

	/// @brief protocols::moves::Mover-specific parsing required by RosettaScripts
	void parse_my_tag(utility::tag::TagCOP tag,
		basic::datacache::DataMap& data,
		const protocols::filters::Filters_map& filters,
		const protocols::moves::Movers_map& movers,
		const core::pose::Pose& pose);

	string chain() const { return chain_; }

	void chain( string const& setting ) { chain_ = setting; }

private:
	string sequence_;
	string residue_type_set_;
	string chain_;
};

}  // namespace simple_moves
}  // namespace protocols

#endif // INCLUDED_protocols_simple_moves_ExtendedPoseMover_HH
