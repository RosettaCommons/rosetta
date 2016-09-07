// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/MatchPositionModifiers.hh
/// @brief  header file for MatchPositionModifiers
/// @author Florian Richter (floric@u.washington.edu ), may 2010

#ifndef INCLUDED_protocols_match_MatchPositionModifiers_hh
#define INCLUDED_protocols_match_MatchPositionModifiers_hh

// unit headers
#include <protocols/match/MatchPositionModifiers.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

//package headers
#include <protocols/match/MatcherTask.fwd.hh>

//project headers
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// utility headers
#include <utility/vector1.hh>

// C++ headers
#include <set>

namespace protocols {
namespace match {

/// @brief "factory" function to create the match position modifiers
MatchPositionModifierCOP
create_match_position_modifier(
	std::string const & mpm_name,
	core::Size geom_cst,
	utility::vector1< std::string > const & input_tokens );

/// @brief base class for objects that modify the match positions based
/// on some criterion
class MatchPositionModifier : public utility::pointer::ReferenceCount {

public:
	MatchPositionModifier();

	~MatchPositionModifier() override;

	/// @brief the positions in the vector1< Size > that is returned
	/// by this function will replace the match positions in the MatcherTask.
	virtual
	utility::vector1< core::Size >
	modified_match_positions(
		utility::vector1< core::Size > const & original_positions,
		core::pose::Pose const & match_pose,
		protocols::match::MatcherTaskCOP mtask
	) const = 0;

};

/// @brief removes positions at which the pose does not have the
/// desired secondary structure
class SecondaryStructureMPM : public MatchPositionModifier {

public:
	SecondaryStructureMPM( utility::vector1< std::string > const & input_tokens );

	~SecondaryStructureMPM() override;


	utility::vector1< core::Size >
	modified_match_positions(
		utility::vector1< core::Size > const & original_positions,
		core::pose::Pose const & match_pose,
		protocols::match::MatcherTaskCOP mtask
	) const override;

private:

	std::set< char > desired_ss_chars_;
	utility::vector1< std::string > ss_motifs_;

};

/// @brief removes positions whose numer of neighbors
/// below a 10A cutoff is not within the desired range.
/// if either min_neighbors_ or max_neighbors_ are unspecified (0),
/// this means that they won't be taken into account, i.e.
/// if min is 5 and max is 0, every position that has
/// more than 4 neighbors will be allowed.
/// also offers the possibility of combining the num neighbors
/// cutoff with the angle between the CA->CB vector of the residue
/// and the CA->protein_center_of_mass vector, for example to
/// only allow positions that point inward
class NumNeighborsMPM : public MatchPositionModifier {

public:
	NumNeighborsMPM( utility::vector1< std::string > const & input_tokens );

	~NumNeighborsMPM() override;


	utility::vector1< core::Size >
	modified_match_positions(
		utility::vector1< core::Size > const & original_positions,
		core::pose::Pose const & match_pose,
		protocols::match::MatcherTaskCOP mtask
	) const override;

	bool
	passes_com_vector_criterion(
		core::Size seqpos,
		core::pose::Pose const & pose,
		core::Vector const & com
	) const;

private:

	core::Size min_neighbors_, max_neighbors_;

	bool com_vector_criterion_;  //whether to use center of mass vector
	bool both_criteria_needed_to_pass_;  //if fullfilling the com vector criterion and the num neighbors criterion is necessary for a position to pass
	core::Real min_com_vector_ang_cos_, max_com_vector_ang_cos_;

};


/// @brief removes positions at which the bfactors for
/// c-alpha atoms are above a desired cutoff. bfactors
/// stored in the pose pdbinfo are taken.
/// if relative bfactors are used, all bfactors are divided
/// by the largest observed bfactor
class BfactorMPM : public MatchPositionModifier {

public:
	BfactorMPM( utility::vector1< std::string > const & input_tokens );

	~BfactorMPM() override;


	utility::vector1< core::Size >
	modified_match_positions(
		utility::vector1< core::Size > const & original_positions,
		core::pose::Pose const & match_pose,
		protocols::match::MatcherTaskCOP mtask
	) const override;

	utility::vector1< core::Real >
	get_ca_bfactors( core::pose::Pose const & pose ) const;

private:

	bool use_relative_bfactors_;
	mutable bool all_bfactors_zero_; //safeguard against pdbs that had their bfactors wiped
	core::Real max_bfactor_;

};

/// @brief MPM that returns a vector of all protein positions in the pose
/// i.e. allowing matching everywhere
class AddAllPositionsMPM : public MatchPositionModifier {

public:
	AddAllPositionsMPM();

	~AddAllPositionsMPM() override;


	utility::vector1< core::Size >
	modified_match_positions(
		utility::vector1< core::Size > const & original_positions,
		core::pose::Pose const & match_pose,
		protocols::match::MatcherTaskCOP mtask
	) const override;


private:

};

/// @brief added by olga and flo 1/2011
/// class to exclude positions at the n and c termini of proteins from matching
class RemoveNorCTermMPM : public MatchPositionModifier {

public:

	RemoveNorCTermMPM( utility::vector1< std::string > const & input_tokens );

	~RemoveNorCTermMPM() override;


	utility::vector1< core::Size >
	modified_match_positions(
		utility::vector1< core::Size > const & original_positions,
		core::pose::Pose const & match_pose,
		protocols::match::MatcherTaskCOP mtask
	) const override;


private:
	core::Size cterm_length_, nterm_length_;

};

/// @brief mpm that will get a task operation as specified in the tag
/// from the TaskOperationFactory, apply the task operation to the pose
/// and every residue that is then set to designing in the task will be
/// a match position
class TaskOperationMPM : public MatchPositionModifier {

public:

	TaskOperationMPM(
		core::Size geom_cst,
		utility::vector1< std::string > const & input_tokens
	);

	~TaskOperationMPM() override;


	utility::vector1< core::Size >
	modified_match_positions(
		utility::vector1< core::Size > const & original_positions,
		core::pose::Pose const & match_pose,
		protocols::match::MatcherTaskCOP mtask
	) const override;


private:
	core::Size which_geom_cst_;
	core::pack::task::operation::TaskOperationOP task_op_;
};


}
}

#endif
