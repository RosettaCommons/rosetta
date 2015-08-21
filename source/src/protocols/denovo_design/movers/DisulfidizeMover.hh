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
/// @details
/// @author Tom Linsky (tlinsky@uw.edu) -- Adapting code from remodelmover into a mover
/// @author Gabe Rocklin (grocklin@uw.edu) -- Disulfide code


#ifndef INCLUDED_protocols_denovo_design_movers_DisulfidizeMover_hh
#define INCLUDED_protocols_denovo_design_movers_DisulfidizeMover_hh

// Unit headers
#include <protocols/denovo_design/movers/DisulfidizeMover.fwd.hh>

// Protocol headers

// Package headers
#include <protocols/rosetta_scripts/MultiplePoseMover.hh>

// Core headers
#include <core/pack/task/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.fwd.hh>

// C++ headers
#include <set>

namespace protocols {
namespace denovo_design {

class DisulfidizeMover : public protocols::rosetta_scripts::MultiplePoseMover {
public:
	typedef utility::vector1< std::pair< core::Size, core::Size > > DisulfideList;
	typedef std::list< core::pose::PoseOP > PoseList;

	/// @brief Default constructor
	///
	DisulfidizeMover();

	/// @brief Copy constructor
	///
	DisulfidizeMover( DisulfidizeMover const &src );

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

	// public methods
public:
	/// @brief finds existing disulfides within a pose
	DisulfideList find_current_disulfides( core::pose::Pose const & pose ) const;

	/// @brief mutates the given disulfides to ALA
	void mutate_disulfides_to_ala(
			core::pose::Pose & pose,
			DisulfideList const & current_ds ) const;

	/// @brief Function for recursively creating multiple disulfides
	utility::vector1< DisulfideList >
		recursive_multiple_disulfide_former(
				DisulfideList const & disulfides_formed,
				DisulfideList const & disulfides_possible ) const;

	/// @brief creates a residue tags on disulfides to inform users that this disulfide was created by disulfidize
	void tag_disulfide(
			core::pose::Pose & pose,
			core::Size const res1,
			core::Size const res2 ) const;

	/// @brief creates a residue tags on disulfides to inform users that this disulfide was created by disulfidize
	void tag_disulfides(
			core::pose::Pose & pose,
			DisulfidizeMover::DisulfideList const & disulf ) const;

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
	DisulfideList find_possible_disulfides(
			core::pose::Pose const & pose,
			core::pack::task::residue_selector::ResidueSubset const & residueset1,
			core::pack::task::residue_selector::ResidueSubset const & residueset2 ) const;

	/// @brief find disulfides in the given neighborhood between residues in set 1 and residues in set 2
	DisulfideList find_possible_disulfides(
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
			core::scoring::disulfides::DisulfideMatchingPotential const & disulfPot,
			bool const mirror
	) const;
	
	/// @brief Returns true if this is a mixed D/L disulfide, false otherwise.
	///
	bool mixed_disulfide (
		core::pose::Pose const &pose,
		core::Size const res1,
		core::Size const res2
	) const;

public: //mutators
	/// @brief sets the selector for set 1 -- disulfides will connect residues in set 1 to residues in set 2
	void set_set1_selector( core::pack::task::residue_selector::ResidueSelectorCOP selector );

	/// @brief sets the selector for set 2 -- disulfides will connect residues in set 1 to residues in set 2
	void set_set2_selector( core::pack::task::residue_selector::ResidueSelectorCOP selector );

	/// @brief sets the min_loop value (number of residues between disulfide-joined residues) (default=8)
	void set_min_loop( core::Size const minloopval );

	/// @brief sets the maximum allowed per-disulfide dslf_fa13 score (default=0.0)
	void set_max_disulf_score( core::Real const maxscoreval );

	/// @brief sets the maximum allowed "match-rt-limit" (default=1.0)
	/// @details Changed from default=2.0 by VKM on 19 Aug 2015 after fixing major bug in
	/// DisulfideMatchPotential.
	void set_match_rt_limit( core::Real const matchrtval );

	/// @brief Set the types of cysteines that we design with:
	/// @details By default, we use only L-cysteine (not D-cysteine).
	void set_cys_types( bool const lcys, bool const dcys );

protected:
	/// @brief Identifies disulfides for a given input pose
	virtual bool process_pose(
			core::pose::Pose & pose,
			utility::vector1 < core::pose::PoseOP > & additional_poses );

private:   // options
	core::Real match_rt_limit_;
	core::Real max_disulf_score_;
	core::Size min_loop_;
	core::Size min_disulfides_ ;
	core::Size max_disulfides_ ;
	bool include_current_ds_;
	bool keep_current_ds_;
	bool score_or_matchrt_;

private:   // other data
	/// @brief disulfides connect residues from set1 to residues from set2
	core::pack::task::residue_selector::ResidueSelectorCOP set1_selector_;
	core::pack::task::residue_selector::ResidueSelectorCOP set2_selector_;

	/// @brief Can disulfides involve L-cystine?
	/// @details Default true
	bool allow_l_cys_;

	/// @brief Can disulfides involve D-cystine?
	/// @details Default false
	bool allow_d_cys_;

};

} // denovo_design
} // protocols

#endif
