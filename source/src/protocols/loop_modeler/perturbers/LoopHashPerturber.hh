// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// @author Xingjie Pan (xingjiepan@gmail.com)

#ifndef INCLUDED_protocols_loop_modeler_perturbers_LoopHashPerturber_HH
#define INCLUDED_protocols_loop_modeler_perturbers_LoopHashPerturber_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/loop_modeler/perturbers/LoopHashPerturber.fwd.hh>

// Protocol headers
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/BackboneDB.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>


namespace protocols {
namespace loop_modeler {
namespace perturbers {

/// @brief Sample backbone torsions using loop hash
class LoopHashPerturber : public protocols::kinematic_closure::perturbers::Perturber {

	typedef utility::vector1< std::pair< protocols::loophash::BackboneSegment, std::string > > BackboneSegments;

public:

	/// @brief constructor
	LoopHashPerturber(protocols::loophash::LoopHashLibraryOP lh_library);

	/// @copydoc Perturber::get_name
	std::string get_name() const { return "LoopHashPerturber"; }

	/// @copydoc Perturber::perturb_subset()
	void perturb_subset(
		core::pose::Pose const & pose,
		kinematic_closure::IndexList const & residues,
		kinematic_closure::ClosureProblemOP problem);

	/// @copydoc Perturber::perturb_subset_with_balance()
	void perturb_subset_with_balance(
		core::pose::Pose const & pose,
		kinematic_closure::IndexList const & residues,
		kinematic_closure::ClosureProblemOP problem);

	///@brief Find a backbone segment from the loop hash library
	void get_backbone_segments(core::pose::Pose const& pose,
		core::Size loophash_fragment_start,
		core::Size loophash_fragment_end);

	///@brief Find a random backbone segment from the loop hash library
	///@detail The selected fragment need not to match the leap
	void get_random_backbone_segments(core::pose::Pose const& pose,
		core::Size loophash_fragment_start,
		core::Size loophash_fragment_end);

	///@brief Extract the fragment information into a pair
	std::pair< protocols::loophash::BackboneSegment, std::string > extract_fragment(core::Size frag_index, core::Size loop_size);

	///@brief Set if use radial loopup
	void use_radial_lookup(bool value){
		use_radial_lookup_ = value;
	}

	///@brief Set if use random mode
	void random_mode(bool value){
		random_mode_ = value;
	}

	///@brief Set if perturb sequence
	void perturb_sequence(bool value){
		perturb_sequence_ = value;
	}

	///@brief Set the sequence positions that should not be mutated
	void seqposes_no_mutate_str(std::string value){
		seqposes_no_mutate_str_ = value;
	}

private:

	protocols::loophash::LoopHashLibraryOP lh_library_;

	bool use_radial_lookup_ = true;

	// In random mode, random fragments from the backbone database will be selected
	bool random_mode_ = false;

	bool perturb_sequence_ = false;

	// Sequence positions that do not mutate
	std::string seqposes_no_mutate_str_;

	BackboneSegments bb_segs_;
	numeric::geometry::hashing::Real6 last_loop_transform_;
};

}
}
}

#endif
