// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/movers/DisulfidizeMover.hh
/// @brief The DisulfidizeMover Protocol
/// @detailed
/// @author Tom Linsky (tlinsky@uw.edu) -- Adapting code from remodelmover into a mover
/// @author Gabe Rocklin (grocklin@uw.edu) -- Disulfide code


#ifndef INCLUDED_protocols_denovo_design_movers_DisulfidizeMover_hh
#define INCLUDED_protocols_denovo_design_movers_DisulfidizeMover_hh

// Unit headers
#include <protocols/denovo_design/movers/DisulfidizeMover.fwd.hh>

// Protocol headers

// Package headers
#include <protocols/moves/Mover.hh>

// Core headers
#include <core/pack/task/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.fwd.hh>

// C++ headers
#include <set>

namespace protocols {
namespace denovo_design {

class DisulfidizeMover : public protocols::moves::Mover {
public:
	typedef utility::vector1< std::pair< core::Size, core::Size > > DisulfideList;
	typedef std::list< core::pose::PoseCOP > PoseList;

	/// @brief default constructor
	DisulfidizeMover();

  /// @brief virtual constructor to allow derivation
	virtual ~DisulfidizeMover();

  /// @brief Parses the DisulfidizerMoverTags
	void parse_my_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap & data,
			protocols::filters::Filters_map const & filters,
			protocols::moves::Movers_map const & movers,
			core::pose::Pose const & pose );

  /// @brief Return the name of this mover.
  virtual std::string get_name() const;

  /// @brief return a fresh instance of this class in an owning pointer
	virtual protocols::moves::MoverOP clone() const;

  /// @brief Apply the DisulfidizerMover. Overloaded apply function from mover base class.
	virtual void apply( core::pose::Pose & pose );

	// public methods
public:
	/// @brief populate the internally cached list of results for a given pose
	void generate_results( core::pose::Pose const & pose );

	/// @brief pushes the given pose to the accumulator
	void push_result( core::pose::PoseCOP pose );

	/// @brief pops the given pose off the accumulator
	core::pose::PoseCOP pop_result();

	/// @brief Function for recursively creating multiple disulfides
	utility::vector1< DisulfideList >
		recursive_multiple_disulfide_former(
				DisulfideList const & disulfides_formed,
				DisulfideList const & disulfides_possible ) const;

	/// @brief forms a disulfide between res1 and res2, optionally allowing backbone movement
	void make_disulfide(
			core::pose::Pose & pose,
			core::Size const res1,
			core::Size const res2,
			bool const relax_bb ) const;

	/// @brief creates disulfides given the list of pairs given
	void make_disulfides(
			core::pose::Pose & pose,
			DisulfideList const & disulf,
			bool const relax_bb ) const;

	/// @brief temporarily tries building a disulfide between the given positions, scores, and restores the pose
	core::Real build_and_score_disulfide(
			core::pose::Pose & blank_pose,
			core::scoring::ScoreFunctionOP sfxn,
			const bool relax_bb,
			core::Size const res1,
			core::Size const res2 ) const;

	/// @brief find disulfides in the given neighborhood
	DisulfideList find_disulfides_in_the_neighborhood(
			core::pose::Pose const & pose,
			core::pack::task::residue_selector::ResidueSubset const & residueset1,
			core::pack::task::residue_selector::ResidueSubset const & residueset2 ) const;

	/// @brief find disulfides in the given neighborhood between residues in set 1 and residues in set 2 
	DisulfideList find_disulfides_in_the_neighborhood(
			core::pose::Pose const & pose,
			std::set< core::Size > const & set1,
			std::set< core::Size > const & set2 ) const;

	/// @brief checks seqpos to ensure that min_loop is satisfied
	bool check_residue_type( core::pose::Pose const & pose, core::Size const res ) const;

	/// @brief checks seqpos to ensure that min_loop is satisfied
	bool check_disulfide_seqpos( core::Size const res1, core::Size const res2 ) const;

	/// @brief checks disulfide CB-CB distance
	bool check_disulfide_cb_distance(
			core::pose::Pose const & pose,
			core::Size const res1,
			core::Size const res2 ) const;

	/// @brief checks disulfide rosetta score
	bool check_disulfide_score(
			core::pose::Pose & pose,
			core::Size const res1,
			core::Size const res2,
			core::scoring::ScoreFunctionOP sfxn ) const;

	/// @brief checks disulfide match rt
	bool check_disulfide_match_rt(
			core::pose::Pose const & pose,
			core::Size const res1,
			core::Size const res2,
			core::scoring::disulfides::DisulfideMatchingPotential const & disulfPot ) const;

private:   // options
	core::Real match_rt_limit_;
	core::Real max_disulf_score_;
	core::Size min_loop_;
	core::Size min_disulfides_ ;
	core::Size max_disulfides_ ;
  bool include_current_ds_;
	bool keep_current_ds_;

private:   // other data
	/// @brief list of results
	PoseList accumulator_;
	/// @brief copy of the pose last seen by apply
	core::pose::PoseCOP last_pose_;

};

} // denovo_design
} // protocols

#endif
