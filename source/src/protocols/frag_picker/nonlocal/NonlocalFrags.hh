// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/frag_picker/nonlocal/NonlocalFrags.hh
/// @author David Kim (dekim@u.washington.edu)

#ifndef INCLUDED_PROTOCOLS_FRAG_PICKER_NONLOCAL_NONLOCALFRAGS_HH
#define INCLUDED_PROTOCOLS_FRAG_PICKER_NONLOCAL_NONLOCALFRAGS_HH

// Unit header
#include <protocols/frag_picker/nonlocal/NonlocalFrags.fwd.hh>

// C/C++ headers
#include <string>
#include <map>

// Package headers

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace frag_picker {
namespace nonlocal {


class NonlocalFrags : public protocols::moves::Mover {

public:
	/// @brief Constructs a new mover
	NonlocalFrags();

	/// @brief Finds interacting non-local fragment pairs
	void apply(core::pose::Pose& pose);

	/// @brief Returns the name of this mover.
	std::string get_name() const;

	/// @brief Creates a copy of this instance
	protocols::moves::MoverOP clone() const;

	/// @brief Creates a new instance by calling the no-argument constructor
	protocols::moves::MoverOP fresh_instance() const;

	/// @brief Registers applicable options
	static void register_options();

private:

	void initialize();


	bool recover_checkpoint( std::string const & tag, core::pose::Pose& pose );

	void write_checkpoint( std::string const & tag, std::string const & data );

	void read_checkpoint_file();

	/* Members */

	bool single_chain_;
	bool relax_input_;
	bool relax_input_with_coordinate_constraints_;
	Size relax_frags_repeats_;

	std::string checkpointfile_;
	std::map< std::string, std::string > checkpoints_map_;

	utility::vector1<Size> frag_sizes_;

	Size min_seq_sep_;
	core::Real ca_dist_squared_;
	Size min_contacts_per_res_;
	core::Real max_rmsd_after_relax_;
	core::Real max_ddg_score_;

	bool output_frags_pdbs_;
	bool output_idealized_;
};

}  // namespace nonlocal
}  // namespace frag_picker
}  // namespace protocols

#endif  // PROTOCOLS_FRAG_PICKER_NONLOCAL_FRAGS_NONLOCALFRAGS_HH_
