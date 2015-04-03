// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/MatcherTask.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_match_MatcherTask_hh
#define INCLUDED_protocols_match_MatcherTask_hh

// Unit headers
#include <protocols/match/MatcherTask.fwd.hh>

// Project headers
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers
#include <numeric/geometry/BoundingBox.hh>
#include <numeric/xyzVector.hh>

// C++ headers
#include <string>
#include <list>
#include <map>

#include <core/id/AtomID.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {

class MatcherTask : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< MatcherTask >
{
public:
	typedef core::Real                               Real;
	typedef core::Size                               Size;
	typedef core::Vector                             Vector;
	typedef numeric::geometry::BoundingBox< Vector > BoundingBox;

public:
	/// Construction and Destruction
	MatcherTask();
	MatcherTask( MatcherTask const & other );
	MatcherTask const & operator = ( MatcherTask const & rhs );

	virtual ~MatcherTask();

	/// self pointers
	inline MatcherTaskCOP get_self_ptr() const { return shared_from_this(); }
	inline MatcherTaskOP get_self_ptr() { return shared_from_this(); }
	//inline MatcherTaskCAP get_self_weak_ptr() const { return MatcherTaskCAP( shared_from_this() ); }
	//inline MatcherTaskAP get_self_weak_ptr() { return MatcherTaskAP( shared_from_this() ); }

public:

	/// Setup
	void set_upstream_pose(
		core::pose::Pose const & pose
	);

	void
	set_downstream_pose(
		core::pose::Pose const & input_pose
	);

	void set_downstream_pose(
		core::pose::Pose const & pose,
		utility::vector1< core::id::AtomID > const & orientation_atoms
	);

	void
	set_downstream_orientation_atoms(
		utility::vector1< core::id::AtomID > const & orientation_atoms
	);

	void
	set_enumerate_ligand_rotamers( bool setting );

	void
	set_only_enumerate_non_match_redundant_ligand_rotamers( bool setting );

	void clear_downstream_orientation_atoms();

	void
	set_ignore_cmdline_for_build_points( bool setting );
	/// @brief Uniformly consider the same set of build points for each of the geometric constrains
	void set_original_scaffold_build_points( utility::vector1< Size > const & resids );

	/// @brief modify the match positions according to what is specified
	/// in the cstfile
	void
	modify_pose_build_resids_from_endes_input();

	/// @brief Set up the task so that it keeps different backbone build points for each
	/// geometric constraint -- the task then needs to know how many geometric constraints
	/// there are.
	void use_different_build_points_for_each_geometric_constraint(
		Size n_geometric_constraints
	);

	/// @brief Set the build point id's for a particular geometric constraint
	void set_original_scaffold_build_points_for_geometric_constraint(
		Size geom_cst_id,
		utility::vector1< Size > const & resids
	);

	void
	define_active_site_from_gridlig_file( std::string const & file_name );

	void
	define_active_site_from_residue_radii_list();

	void
	append_upstream_resiue_as_defining_active_site( Size resid, Real radius );

	void
	append_downstream_atom_to_active_site_required_list( core::id::AtomID atid );


	//std::list< core::id::AtomID > const & downstream_atoms_required_inside_active_site_grid() const;

	/// @brief Set the bounding box for the region of space that the OccupiedSpaceHash should accept hits
	/// in side of.  If the 3rd orientation atom of the downstream partner is outside of this bounding
	/// box, the hit will be rejected out of hand.
	void set_occupied_space_bounding_box( BoundingBox const & bb );
	/// @brief For the occupied space hash, set the euclidean-bin width to a uniform value for xy&z
	void set_hash_euclidean_bin_width( Real width );
	/// @brief For the occupeid space hash, set the euler-bin width to a uniform value for phi,psi&theta
	void set_hash_euler_bin_width(     Real width );

	/// @brief For collision detection, select the amount of collision that should
	/// be tolerated between heavy atoms.  This should be a positive value in Angstroms.
	void set_permitted_overlap( Real permitted_overlap );

	/// @brief Initialize many parameters from the command line options
	void initialize_from_command_line();

	/// @brief Matches may either be output as soon as they are generated, or they may be consolidated.
	/// When consolidating matches, the MatchProcessor (the MatchEvaluator)
	/// waits until all matches are seen before outputting any; it groups matches and selects the top N
	/// matches from each group to output.  This can require a lot of memory if there are very many
	/// possible matches.  MatchConsolidation is on by default.
	void consolidate_matches( bool setting );
	/// @brief For use with the match consolidator; specify the number of output matches that the
	/// consolidator should select for each group.
	void n_to_output_per_group( Size setting );
	/// @brief Add a filter by name to the set of filters being included.  If that filter requires
	/// extra data (as, for example, the UpstreamCollisionFilter) then the task should be expanded
	/// to include all the data necessary to create and initialize that filter.  No valid options
	/// currently.
	void add_filter( std::string const & filter_name );
	/// @brief Specify the name of the match-consolidator related match-grouper class.
	/// This class will group matches together; the consolidator will then pick the top N from
	/// each group for output.  Valid options include: SameChiBinComboGrouper,
	/// SameSequenceGrouper, and SameRotamerComboGrouper.
	void grouper_name( std::string const & setting );
	/// @brief Specify the name of the match-consolidator related match-evaluator class.
	/// This class will rank each of the matches so that the consolidator may pick the top N.
	/// Valid options include: DownstreamRMSEvaluator.  More evaluator options will be implemented shortly.
	void evaluator_name( std::string const & setting );
	/// @brief Specify the name of the class that will write the output.
	/// Valid options include: KinWriter.  More output options will be implemented shortly.
	void output_writer_name( std::string const & setting );
	/// @brief Indicate the name of the single output file to which the matches will be written
	void output_file_name( std::string const & setting );

	/// @brief Set the name of the single output file to which the scores of matches will be written
	void score_output_file_name( std::string const & setting );

	/// @brief Set the matcher-file input data.  The Matcher will read this data when initializing itself.
	void set_enz_input_data( toolbox::match_enzdes_util::EnzConstraintIOCOP data );

	void filter_upstream_residue_collisions( bool setting );
	void filter_upstream_collisions_by_score( bool setting );
	void upstream_residue_collision_tolerance( Real setting );
	void upstream_residue_collision_score_cutoff( Real setting );
	void upstream_residue_collision_Wfa_atr( Real setting );
	void upstream_residue_collision_Wfa_rep( Real setting );
	void upstream_residue_collision_Wfa_sol( Real setting );

	void filter_upstream_downstream_collisions( bool setting );
	void filter_upstream_downstream_collisions_by_score( bool setting );
	void upstream_downstream_atom_collision_tolerance( Real setting );
	void upstream_downstream_residue_collision_score_cutoff( Real setting );
	void upstream_downstream_residue_collision_Wfa_atr( Real setting );
	void upstream_downstream_residue_collision_Wfa_rep( Real setting );
	void upstream_downstream_residue_collision_Wfa_sol( Real setting );

	void define_match_by_single_downstream_positioning( bool setting );


public:  // Accessors

	core::pose::PoseCOP
	upstream_pose() const;

	core::pose::PoseCOP
	downstream_pose() const;

	utility::vector1< core::id::AtomID > const &
	downstream_orientation_atoms() const;

	bool
	enumerate_ligand_rotamers() const;

	bool
	only_enumerate_non_match_redundant_ligand_rotamers() const;

	utility::vector1< Size > const &
	upstream_pose_build_resids_for_geometric_constraint( Size cst_id ) const;

	std::map< core::Size, core::Size > const &
	upstream_only_geom_cst() const;

	/// @brief Define the active site through a gridlig file (true), or by listing residue/radii paris (false)?
	bool
	gridlig_active_site_definition() const;

	/// @brief Accessor for the file name containing the active-site definition in gridlig format
	std::string const &
	gridlig_file_name() const;

	/// @brief Accessor for the data defining the active site by-residue.  This data is only
	/// active if gridlig_active_site_definition() returns false.
	std::list< std::pair< Size, Real > > const &
	upstream_resids_and_radii_defining_active_site() const;

	std::list< core::id::AtomID > const &
	downstream_atoms_required_inside_active_site() const;

	utility::vector1< core::id::AtomID > const &
	relevant_downstream_atoms() const;

	BoundingBox const &
	occ_space_bounding_box() const;

	Vector euclidean_bin_widths() const;
	Vector euler_bin_widths() const;

	Real permitted_overlap() const;

	bool use_input_sc() const;
	bool dynamic_grid_refinement() const;
	bool consolidate_matches() const;

	Size n_to_output_per_group() const;

	std::list< std::string > const &
	filter_names() const;

	std::string const & upstream_pose_name() const;
	std::string const & cstfile_name() const;
	std::string const & grouper_name() const;
	std::string const & evaluator_name() const;
	std::string const & output_writer_name() const;
	std::string const & output_file_name() const;
	Real grouper_ds_rmsd() const;

	std::string const & score_output_file_name() const;
	bool output_scores() const;

	bool output_matchres_only() const;

	utility::vector1< core::Size > const & geom_csts_downstream_output() const;

	toolbox::match_enzdes_util::EnzConstraintIOCOP
	enz_input_data() const;

	bool filter_upstream_residue_collisions() const;
	bool filter_upstream_collisions_by_score() const;
	Real upstream_residue_collision_tolerance() const;
	Real upstream_residue_collision_score_cutoff() const;
	Real upstream_residue_collision_Wfa_atr() const;
	Real upstream_residue_collision_Wfa_rep() const;
	Real upstream_residue_collision_Wfa_sol() const;

	bool filter_upstream_downstream_collisions() const;
	bool filter_upstream_downstream_collisions_by_score() const;
	Real upstream_downstream_atom_collision_tolerance() const;
	Real upstream_downstream_residue_collision_score_cutoff() const;
	Real upstream_downstream_residue_collision_Wfa_atr() const;
	Real upstream_downstream_residue_collision_Wfa_rep() const;
	Real upstream_downstream_residue_collision_Wfa_sol() const;

	bool define_match_by_single_downstream_positioning() const;

	bool build_round1_hits_twice() const;

private:

	void
	validate_downstream_orientation_atoms() const;

	/// @brief Read the file describing the occupied space grid euclidean dimensions.
	/// The "details" tag for this function describes the file format used.
	void
	initialize_occupied_space_bounding_box_from_command_line();

	/// @brief Read one of two files given on the command line that defines the set
	/// of residues on the scaffold to consider as potential launch points for the
	/// scaffold's active site.  File formats are described in the "details" tag.
	void
	initialize_scaffold_active_site_residue_list_from_command_line();

	/// @brief in cases where the upstream pose to be matched already
	/// contains some of the desired interactions (as specified in the
	/// REMARK header, the match position list for every geomcst will
	/// be set to these positions
	void
	set_active_site_residue_list_to_preexisting_partial_match();

	/// @brief in case the upstream pose containts a copy
	/// of the downstream object (i.e. if a previously matched
	/// partial match is being read in again )
	void
	remove_downstream_object_from_upstream_pose();

	void
	initialize_enzdes_input_data_from_command_line();

	/// @brief queries the enzdes input for which atoms are relevant to the matcher,
	/// i.e. which atoms in the downstream object interact with any of the match residues
	void
	determine_all_match_relevant_downstream_atoms();

	void
	initialize_orientation_atoms_from_command_line();

	/// @brief Read the command line arguments specifying the subset of downstream
	/// partner atoms that are required to be in the active site, as well as a definition
	/// of the region called the active site.  The "details" tag for this function
	/// describes three file formats used in this function.
	void
	initialize_active_site_definition_from_command_line();

	void
	initialize_upstream_residue_collision_filter_data_from_command_line();

	void
	initialize_upstream_downstream_collision_filter_data_from_command_line();

	void
	initialize_output_options_from_command_line();

private:

	core::pose::PoseCOP upstream_pose_;

	core::pose::PoseCOP downstream_pose_;
	utility::vector1< core::id::AtomID > downstream_orientation_atoms_;
	utility::vector1< core::id::AtomID > relevant_downstream_atoms_;
	bool enumerate_ligand_rotamers_;
	bool only_enumerate_non_match_redundant_ligand_rotamers_;

	bool ignore_cmdline_for_build_points_;
	bool share_build_points_for_geomcsts_;
	utility::vector1< Size > generic_pose_build_resids_;
	utility::vector1< utility::vector1< Size > > per_cst_pose_build_resids_;

	//for upstream only constraints, this maps between the geom cst
	//and the geom cst of the upstream target
	std::map< core::Size, core::Size > upstream_only_geom_cst_;

	bool gridlig_active_site_definition_;
	std::string gridlig_fname_;
	std::string upstream_pose_name_;
	std::string cstfile_name_;
	std::list< std::pair< Size, Real > > upstream_resids_and_radii_defining_active_site_;

	std::list< core::id::AtomID > downstream_atoms_required_inside_active_site_;
	utility::vector1< Size > downstream_atom_inds_used_in_matching_;

	BoundingBox occ_space_bounding_box_;
	Vector euclidean_bin_widths_;
	Vector euler_bin_widths_;

	Real permitted_overlap_;

	bool use_input_sc_;
	bool dynamic_grid_refinement_;
	bool consolidate_matches_; /// MatchConsolidator vs MatchOutputter
	Size n_to_output_per_group_;
	std::list< std::string > filter_names_;
	std::string grouper_name_;
	std::string evaluator_name_;
	std::string output_writer_name_;
	std::string output_file_name_;
	Real grouper_ds_rmsd_;

	// Score output
	std::string score_output_file_name_;
	bool output_scores_;

	//some options for outputting
	bool output_matchres_only_;
	utility::vector1< core::Size > geom_csts_downstream_output_; //for which of the geometric constraints will the ligand be output

	toolbox::match_enzdes_util::EnzConstraintIOCOP enz_input_data_;

	bool filter_upstream_residue_collisions_;
	bool filter_upstream_collisions_by_score_;
	Real upstream_residue_collision_tolerance_;
	Real upstream_residue_collision_score_cutoff_;
	Real upstream_residue_collision_Wfa_atr_;
	Real upstream_residue_collision_Wfa_rep_;
	Real upstream_residue_collision_Wfa_sol_;

	bool filter_upstream_and_downstream_residue_collisions_;
	bool filter_upstream_and_downstream_collisions_by_score_;
	Real upstream_downstream_atom_collision_tolerance_;
	Real upstream_downstream_residue_collision_score_cutoff_;
	Real upstream_downstream_residue_collision_Wfa_atr_;
	Real upstream_downstream_residue_collision_Wfa_rep_;
	Real upstream_downstream_residue_collision_Wfa_sol_;

	bool define_match_by_single_downstream_positioning_;

	bool build_round1_hits_twice_;

};


}
}

#endif
