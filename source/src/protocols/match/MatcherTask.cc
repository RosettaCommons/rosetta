// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/MatcherTask.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/match/MatcherTask.hh>

// Project headers
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/enzdes/enzdes_util.hh>
#include <protocols/toolbox/match_enzdes_util/AllowedSeqposForGeomCst.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>

#include <protocols/match/MatchPositionModifiers.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>

#include <core/conformation/Residue.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/match.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/func/Func.hh>

#include <basic/Tracer.hh>

/// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>

// C++ headers
#include <sstream>

#include <utility/vector1.hh>


namespace protocols {
namespace match {

static THREAD_LOCAL basic::Tracer TR( "protocols.match.MatcherTask" );

MatcherTask::MatcherTask() :
	enumerate_ligand_rotamers_( false ),
	only_enumerate_non_match_redundant_ligand_rotamers_( false ),
	ignore_cmdline_for_build_points_(false),
	share_build_points_for_geomcsts_( true ),
	gridlig_active_site_definition_( false ),
	occ_space_bounding_box_( Vector( 0.0 ), Vector( 0.0 ) ),
	euclidean_bin_widths_( 1.0 ),
	euler_bin_widths_( 10.0 ),
	permitted_overlap_( 0.0 ),
	use_input_sc_( false ),
	dynamic_grid_refinement_( false ),
	consolidate_matches_( true ),
	n_to_output_per_group_( 5 ),
	grouper_name_( "SameChiBinComboGrouper" ),
	evaluator_name_( "DownstreamRMSEvaluator" ),
	output_writer_name_( "PDB" ),
	output_file_name_( "output_matches.kin" ),
	grouper_ds_rmsd_(1.5),
	score_output_file_name_( "matcher_score.sc" ),
	output_scores_( false ),
	output_matchres_only_(true),
	filter_upstream_residue_collisions_( false ),
	filter_upstream_collisions_by_score_( false ),
	upstream_residue_collision_tolerance_( 1.0 ),
	upstream_residue_collision_score_cutoff_( 1.0 ),
	upstream_residue_collision_Wfa_atr_( 1.0 ),
	upstream_residue_collision_Wfa_rep_( 1.0 ),
	upstream_residue_collision_Wfa_sol_( 1.0 ),
	filter_upstream_and_downstream_residue_collisions_( false ),
	filter_upstream_and_downstream_collisions_by_score_( false ),
	upstream_downstream_atom_collision_tolerance_( 1.0 ),
	upstream_downstream_residue_collision_score_cutoff_( 1.0 ),
	upstream_downstream_residue_collision_Wfa_atr_( 1.0 ),
	upstream_downstream_residue_collision_Wfa_rep_( 1.0 ),
	upstream_downstream_residue_collision_Wfa_sol_( 1.0 ),
	define_match_by_single_downstream_positioning_( false ),
	build_round1_hits_twice_( false )
{}

MatcherTask::MatcherTask( MatcherTask const & other ) :
	ReferenceCount(),
	utility::pointer::enable_shared_from_this< MatcherTask >()
{
	(*this) = other;
}

MatcherTask &
MatcherTask::operator = ( MatcherTask const & rhs )
{
	if ( this != & rhs ) {
		upstream_pose_ = rhs.upstream_pose_;
		downstream_pose_ = rhs.downstream_pose_;
		downstream_orientation_atoms_ = rhs.downstream_orientation_atoms_;
		relevant_downstream_atoms_ = rhs.relevant_downstream_atoms_;
		enumerate_ligand_rotamers_ = rhs.enumerate_ligand_rotamers_;
		only_enumerate_non_match_redundant_ligand_rotamers_ = rhs.only_enumerate_non_match_redundant_ligand_rotamers_;
		ignore_cmdline_for_build_points_ = rhs.ignore_cmdline_for_build_points_;
		share_build_points_for_geomcsts_ = rhs.share_build_points_for_geomcsts_;
		generic_pose_build_resids_ = rhs.generic_pose_build_resids_;
		per_cst_pose_build_resids_ = rhs.per_cst_pose_build_resids_;
		upstream_only_geom_cst_ = rhs.upstream_only_geom_cst_;
		gridlig_active_site_definition_ = rhs.gridlig_active_site_definition_;
		gridlig_fname_ = rhs.gridlig_fname_;
		upstream_pose_name_ = rhs.upstream_pose_name_;
		cstfile_name_ = rhs.cstfile_name_;
		upstream_resids_and_radii_defining_active_site_ =
			rhs.upstream_resids_and_radii_defining_active_site_;
		downstream_atoms_required_inside_active_site_ =
			rhs.downstream_atoms_required_inside_active_site_;
		occ_space_bounding_box_ = rhs.occ_space_bounding_box_;
		euclidean_bin_widths_ = rhs.euclidean_bin_widths_;
		euler_bin_widths_ = rhs.euler_bin_widths_;
		permitted_overlap_ = rhs.permitted_overlap_;
		use_input_sc_ = rhs.use_input_sc_;
		dynamic_grid_refinement_ = rhs.dynamic_grid_refinement_;
		consolidate_matches_ = rhs.consolidate_matches_;
		n_to_output_per_group_ = rhs.n_to_output_per_group_;
		filter_names_ = rhs.filter_names_;
		grouper_name_ = rhs.grouper_name_;
		evaluator_name_ = rhs.evaluator_name_;
		output_writer_name_ = rhs.output_writer_name_;
		output_file_name_ = rhs.output_file_name_;
		score_output_file_name_ = rhs.score_output_file_name_;
		output_scores_ = rhs.output_scores_;
		grouper_ds_rmsd_ = rhs.grouper_ds_rmsd_;
		output_matchres_only_ = rhs.output_matchres_only_;
		geom_csts_downstream_output_ = rhs.geom_csts_downstream_output_;
		enz_input_data_ = rhs.enz_input_data_;
		filter_upstream_residue_collisions_ = rhs.filter_upstream_residue_collisions_;
		filter_upstream_collisions_by_score_ = rhs.filter_upstream_collisions_by_score_;
		upstream_residue_collision_tolerance_ = rhs.upstream_residue_collision_tolerance_;
		upstream_residue_collision_score_cutoff_ = rhs.upstream_residue_collision_score_cutoff_;
		upstream_residue_collision_Wfa_atr_ = rhs.upstream_residue_collision_Wfa_atr_;
		upstream_residue_collision_Wfa_rep_ = rhs.upstream_residue_collision_Wfa_rep_;
		upstream_residue_collision_Wfa_sol_ = rhs.upstream_residue_collision_Wfa_sol_;
		define_match_by_single_downstream_positioning_ = rhs.define_match_by_single_downstream_positioning_;
		filter_upstream_and_downstream_residue_collisions_ = rhs.filter_upstream_and_downstream_residue_collisions_;
		filter_upstream_and_downstream_collisions_by_score_ = rhs.filter_upstream_and_downstream_collisions_by_score_;
		upstream_downstream_atom_collision_tolerance_ = rhs.upstream_downstream_atom_collision_tolerance_;
		upstream_downstream_residue_collision_score_cutoff_ = rhs.upstream_downstream_residue_collision_score_cutoff_;
		upstream_downstream_residue_collision_Wfa_atr_ = rhs.upstream_downstream_residue_collision_Wfa_atr_;
		upstream_downstream_residue_collision_Wfa_rep_ = rhs.upstream_downstream_residue_collision_Wfa_rep_;
		upstream_downstream_residue_collision_Wfa_sol_ = rhs.upstream_downstream_residue_collision_Wfa_sol_;
		build_round1_hits_twice_ = rhs.build_round1_hits_twice_;
	}
	return *this;
}

MatcherTask::~MatcherTask() = default;

void
MatcherTask::set_upstream_pose(
	core::pose::Pose const & input_pose
)
{
	upstream_pose_ = core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( input_pose ) ) );
	utility::vector1< std::string > path = utility::string_split( input_pose.pdb_info()->name(), '/' );
	upstream_pose_name_ = utility::string_split( path[ path.size() ], '.' )[1];
}

void
MatcherTask::set_downstream_pose(
	core::pose::Pose const & input_pose
)
{
	downstream_pose_ = core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( input_pose ) ) );

	if ( downstream_orientation_atoms_.size() != 0 ) {
		validate_downstream_orientation_atoms();
	}
}

void
MatcherTask::set_downstream_pose(
	core::pose::Pose const & input_pose,
	utility::vector1< core::id::AtomID > const & orientation_atoms
)
{
	downstream_pose_ = core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( input_pose ) ) );
	downstream_orientation_atoms_ = orientation_atoms;

	validate_downstream_orientation_atoms();
}

void
MatcherTask::set_downstream_orientation_atoms(
	utility::vector1< core::id::AtomID > const & orientation_atoms
)
{
	downstream_orientation_atoms_ = orientation_atoms;
	if ( downstream_pose_ ) {
		validate_downstream_orientation_atoms();
	}
}

void
MatcherTask::clear_downstream_orientation_atoms() {
	downstream_orientation_atoms_.clear();
}

void
MatcherTask::set_ignore_cmdline_for_build_points( bool setting ){
	ignore_cmdline_for_build_points_ = setting;
}

void
MatcherTask::set_original_scaffold_build_points( utility::vector1< Size > const & resids )
{
	if ( share_build_points_for_geomcsts_ ) {
		generic_pose_build_resids_ = resids;
	} else {
		for ( Size ii = 1; ii <= per_cst_pose_build_resids_.size(); ++ii ) {
			per_cst_pose_build_resids_[ ii ] = resids;
		}
	}
}


void
MatcherTask::modify_pose_build_resids_from_endes_input()
{

	runtime_assert( enz_input_data_ != nullptr );
	runtime_assert( upstream_pose_ != nullptr );
	bool switch_to_per_cst_resids_necessary( share_build_points_for_geomcsts_ );

	for ( core::Size i =1; i<= enz_input_data_->num_mcfi_lists(); ++i ) {

		std::map< std::string, utility::vector1< std::string > > const & alg_map = (*(enz_input_data_->mcfi_list( i )->mcfi( 1 ))).algorithm_inputs();
		auto map_it( alg_map.find( "match_positions" ) );
		if ( map_it == alg_map.end() ) continue;

		// silly: this could mean that different positions are used for each cst, so we might have to switch
		if ( switch_to_per_cst_resids_necessary ) {
			share_build_points_for_geomcsts_ = false;
			switch_to_per_cst_resids_necessary = false;
			per_cst_pose_build_resids_.resize( enz_input_data_->num_mcfi_lists() );
			for ( core::Size j =1; j <= per_cst_pose_build_resids_.size(); ++j ) per_cst_pose_build_resids_[j] = generic_pose_build_resids_;
			generic_pose_build_resids_.clear();
		}

		for ( core::Size mpm_string =1; mpm_string <= map_it->second.size(); ++mpm_string ) {
			utility::vector1< std::string > tokens;
			tokens.push_back(""); //weird util vect1 copy behavior makes this necessary
			tokens = utility::split( map_it->second[mpm_string] );
			MatchPositionModifierCOP mpm( create_match_position_modifier( tokens[1], i, tokens ) );
			if ( !mpm ) utility_exit_with_message("Could not create a MatchPositionModifier based on name "+tokens[1]+".");
			per_cst_pose_build_resids_[i] = mpm->modified_match_positions( per_cst_pose_build_resids_[ i ], *upstream_pose_, get_self_ptr() );
		}

		TR << "Match position modifiers changed match positions for geomcst " << i << " to the following: " << std::endl;
		for ( core::Size newpos = 1; newpos <= per_cst_pose_build_resids_[i].size(); ++newpos ) TR << per_cst_pose_build_resids_[i][newpos] << "+";
		TR << std::endl;

	} //loop over all mcfi lists
}

void
MatcherTask::set_enumerate_ligand_rotamers( bool setting )
{
	enumerate_ligand_rotamers_ = setting;
}

void
MatcherTask::set_only_enumerate_non_match_redundant_ligand_rotamers( bool setting )
{
	only_enumerate_non_match_redundant_ligand_rotamers_ = setting;
}


void MatcherTask::use_different_build_points_for_each_geometric_constraint(
	Size n_geometric_constraints
)
{
	share_build_points_for_geomcsts_ = false;
	per_cst_pose_build_resids_.clear();
	per_cst_pose_build_resids_.resize( n_geometric_constraints );
}

/// @brief Set the build point id's for a particular geometric constraint
void
MatcherTask::set_original_scaffold_build_points_for_geometric_constraint(
	Size geom_cst_id,
	utility::vector1< Size > const & resids
)
{
	runtime_assert( ! share_build_points_for_geomcsts_ ); // use_different_build_points_for_each_geometric_constraint() must be called first
	runtime_assert( geom_cst_id <= per_cst_pose_build_resids_.size() );
	per_cst_pose_build_resids_[ geom_cst_id ] = resids;
}

void
MatcherTask::define_active_site_from_gridlig_file( std::string const & file_name )
{
	if ( ! gridlig_active_site_definition_ && ! upstream_resids_and_radii_defining_active_site_.empty() ) {
		TR << "overriding res/radius active site definition by setting gridlig file " << file_name << std::endl;
	}
	gridlig_active_site_definition_ = true;
	gridlig_fname_ = file_name;
}

void
MatcherTask::define_active_site_from_residue_radii_list()
{
	gridlig_active_site_definition_ = false;
}

void
MatcherTask::append_upstream_resiue_as_defining_active_site(
	Size resid,
	Real radius
)
{
	if ( gridlig_active_site_definition_ ) {
		TR << "overriding gridlig_active_site_definition_ by appending res/radius pair: " << resid << " " << radius << std::endl;
		gridlig_active_site_definition_ = false;
	}
	upstream_resids_and_radii_defining_active_site_.push_back( std::make_pair( resid, radius ) );
}

void
MatcherTask::append_downstream_atom_to_active_site_required_list(
	core::id::AtomID atid
)
{
	downstream_atoms_required_inside_active_site_.push_back( atid );
}


void
MatcherTask::set_occupied_space_bounding_box( BoundingBox const & bb )
{
	occ_space_bounding_box_ = bb;
}

void
MatcherTask::set_hash_euclidean_bin_width( Real width )
{
	euclidean_bin_widths_ = width;
}

void
MatcherTask::set_hash_euler_bin_width(     Real width )
{
	euler_bin_widths_ = width;
}

void
MatcherTask::set_permitted_overlap( Real permitted_overlap )
{
	permitted_overlap_ = permitted_overlap;
}

void
MatcherTask::initialize_from_command_line()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::match;

	/// 1. Read the ligand grid header file and initialize the boundary for the active site.
	/// 2. Initialize the euclid and euler bin witdhs
	/// 3. Note the location of the output file
	/// 4. Read in the EnzConstraintIO file.
	/// 5. Read in the build position list

	/// 1. liggrid
	initialize_occupied_space_bounding_box_from_command_line();

	euclidean_bin_widths_ = option[ euclid_bin_size ];
	euler_bin_widths_ = option[ euler_bin_size ];

	initialize_enzdes_input_data_from_command_line();

	if ( !ignore_cmdline_for_build_points_ ) initialize_scaffold_active_site_residue_list_from_command_line();

	modify_pose_build_resids_from_endes_input();

	//note: this function only does something if the input pose
	//already contains a partial match
	set_active_site_residue_list_to_preexisting_partial_match();

	/// bump tolerance
	permitted_overlap_ = option[ bump_tolerance ];

	use_input_sc_ = option[ OptionKeys::packing::use_input_sc ];
	dynamic_grid_refinement_ = option[ OptionKeys::match::dynamic_grid_refinement ];

	initialize_orientation_atoms_from_command_line();

	consolidate_matches_ = option[ OptionKeys::match::consolidate_matches ];
	n_to_output_per_group_ = option[ output_matches_per_group ];

	initialize_active_site_definition_from_command_line();
	initialize_upstream_residue_collision_filter_data_from_command_line();
	initialize_upstream_downstream_collision_filter_data_from_command_line();

	initialize_output_options_from_command_line();

	enumerate_ligand_rotamers_ = option[ OptionKeys::match::enumerate_ligand_rotamers ];
	only_enumerate_non_match_redundant_ligand_rotamers_ =  option[ OptionKeys::match::only_enumerate_non_match_redundant_ligand_rotamers ];

	build_round1_hits_twice_ = option[ OptionKeys::match::build_round1_hits_twice ];
}

void
MatcherTask::consolidate_matches( bool setting )
{
	consolidate_matches_ = setting;
}

void MatcherTask::filter_upstream_residue_collisions( bool setting )
{
	filter_upstream_residue_collisions_ = setting;
}

void MatcherTask::filter_upstream_collisions_by_score( bool setting )
{
	assert( filter_upstream_residue_collisions_ ); // must set filter_upstream_residue_collisions_ first!
	filter_upstream_collisions_by_score_ = setting;
}

void MatcherTask::upstream_residue_collision_tolerance( Real setting )
{
	upstream_residue_collision_tolerance_ = setting;
}

void MatcherTask::upstream_residue_collision_score_cutoff( Real setting )
{
	upstream_residue_collision_score_cutoff_ = setting;
}

void MatcherTask::upstream_residue_collision_Wfa_atr( Real setting )
{
	upstream_residue_collision_Wfa_atr_ = setting;
}

void MatcherTask::upstream_residue_collision_Wfa_rep( Real setting )
{
	upstream_residue_collision_Wfa_rep_ = setting;
}

void MatcherTask::upstream_residue_collision_Wfa_sol( Real setting )
{
	upstream_residue_collision_Wfa_sol_ = setting;
}

void MatcherTask::filter_upstream_downstream_collisions( bool setting )
{
	filter_upstream_and_downstream_residue_collisions_= setting;
}

void MatcherTask::filter_upstream_downstream_collisions_by_score( bool setting )
{
	filter_upstream_and_downstream_collisions_by_score_ = setting;
}

void MatcherTask::upstream_downstream_atom_collision_tolerance( Real setting )
{
	upstream_downstream_atom_collision_tolerance_ = setting;
}

void MatcherTask::upstream_downstream_residue_collision_score_cutoff( Real setting )
{
	upstream_downstream_residue_collision_score_cutoff_ = setting;
}

void MatcherTask::upstream_downstream_residue_collision_Wfa_atr( Real setting )
{
	upstream_downstream_residue_collision_Wfa_atr_ = setting;
}

void MatcherTask::upstream_downstream_residue_collision_Wfa_rep( Real setting )
{
	upstream_downstream_residue_collision_Wfa_rep_ = setting;
}

void MatcherTask::upstream_downstream_residue_collision_Wfa_sol( Real setting )
{
	upstream_downstream_residue_collision_Wfa_sol_ = setting;
}

void MatcherTask::define_match_by_single_downstream_positioning( bool setting )
{
	define_match_by_single_downstream_positioning_ = setting;
}

void MatcherTask::n_to_output_per_group( Size setting )
{
	n_to_output_per_group_ = setting;
}

void
MatcherTask::add_filter( std::string const & filter_name )
{
	filter_names_.push_back( filter_name );
}

void
MatcherTask::grouper_name( std::string const & setting )
{
	grouper_name_ = setting;
}

void
MatcherTask::evaluator_name( std::string const & setting )
{
	evaluator_name_ = setting;
}

void
MatcherTask::output_writer_name( std::string const & setting )
{
	output_writer_name_ = setting;
}

void
MatcherTask::output_file_name( std::string const & setting )
{
	output_file_name_ = setting;
}

/// @details Sets the output score filename and turns on score output boolean
void
MatcherTask::score_output_file_name( std::string const & setting )
{
	score_output_file_name_ = setting;
	output_scores_ = true;
}

void
MatcherTask::set_enz_input_data( toolbox::match_enzdes_util::EnzConstraintIOCOP data )
{
	enz_input_data_ = data;
	determine_all_match_relevant_downstream_atoms();
}

core::pose::PoseCOP
MatcherTask::upstream_pose() const
{
	return upstream_pose_;
}

core::pose::PoseCOP
MatcherTask::downstream_pose() const
{
	return downstream_pose_;
}

utility::vector1< core::id::AtomID > const &
MatcherTask::downstream_orientation_atoms() const
{
	return downstream_orientation_atoms_;
}

bool
MatcherTask::enumerate_ligand_rotamers() const
{
	return enumerate_ligand_rotamers_;
}

bool
MatcherTask::only_enumerate_non_match_redundant_ligand_rotamers() const
{
	return only_enumerate_non_match_redundant_ligand_rotamers_;
}


utility::vector1< core::Size > const &
MatcherTask::upstream_pose_build_resids_for_geometric_constraint( Size cst_id ) const
{
	if ( share_build_points_for_geomcsts_ ) {
		return generic_pose_build_resids_;
	} else {
		runtime_assert( cst_id <= per_cst_pose_build_resids_.size() );
		return per_cst_pose_build_resids_[ cst_id ];
	}
}

std::map< core::Size, core::Size > const &
MatcherTask::upstream_only_geom_cst() const{
	return upstream_only_geom_cst_;}

/// @brief Define the active site through a gridlig file (true), or by listing residue/radii paris (false)?
bool
MatcherTask::gridlig_active_site_definition() const
{
	return gridlig_active_site_definition_;
}

/// @brief Accessor for the file name containing the active-site definition in gridlig format
std::string const &
MatcherTask::gridlig_file_name() const
{
	runtime_assert( gridlig_active_site_definition_ );
	return gridlig_fname_;
}


std::list< std::pair< MatcherTask::Size, MatcherTask::Real > > const &
MatcherTask::upstream_resids_and_radii_defining_active_site() const
{
	runtime_assert( ! gridlig_active_site_definition_ );
	return upstream_resids_and_radii_defining_active_site_;
}

std::list< core::id::AtomID > const &
MatcherTask::downstream_atoms_required_inside_active_site() const
{
	return downstream_atoms_required_inside_active_site_;
}

utility::vector1< core::id::AtomID > const &
MatcherTask::relevant_downstream_atoms() const
{
	return relevant_downstream_atoms_;
}


MatcherTask::BoundingBox const &
MatcherTask::occ_space_bounding_box() const
{
	return occ_space_bounding_box_;
}

MatcherTask::Vector
MatcherTask::euclidean_bin_widths() const
{
	return euclidean_bin_widths_;
}

MatcherTask::Vector
MatcherTask::euler_bin_widths() const
{
	return euler_bin_widths_;
}

MatcherTask::Real
MatcherTask::permitted_overlap() const
{
	return permitted_overlap_;
}

bool
MatcherTask::use_input_sc() const {
	return use_input_sc_;
}

bool
MatcherTask::dynamic_grid_refinement() const
{
	return dynamic_grid_refinement_;
}

bool
MatcherTask::consolidate_matches() const
{
	return consolidate_matches_;
}

MatcherTask::Size
MatcherTask::n_to_output_per_group() const
{
	return n_to_output_per_group_;
}


std::list< std::string > const &
MatcherTask::filter_names() const
{
	return filter_names_;
}

std::string const &
MatcherTask::upstream_pose_name() const
{
	return upstream_pose_name_;
}

std::string const &
MatcherTask::cstfile_name() const
{
	return cstfile_name_;
}

std::string const &
MatcherTask::grouper_name() const
{
	return grouper_name_;
}

std::string const &
MatcherTask::evaluator_name() const
{
	return evaluator_name_;
}

std::string const &
MatcherTask::output_writer_name() const
{
	return output_writer_name_;
}

std::string const &
MatcherTask::output_file_name() const
{
	return output_file_name_;
}

std::string const &
MatcherTask::score_output_file_name() const
{
	return score_output_file_name_;
}

bool
MatcherTask::output_scores() const
{
	return output_scores_;
}

MatcherTask::Real
MatcherTask::grouper_ds_rmsd() const
{
	return grouper_ds_rmsd_;
}


bool
MatcherTask::output_matchres_only() const
{
	return output_matchres_only_;
}

utility::vector1< core::Size > const &
MatcherTask::geom_csts_downstream_output() const
{
	return geom_csts_downstream_output_;
}

toolbox::match_enzdes_util::EnzConstraintIOCOP
MatcherTask::enz_input_data() const
{
	return enz_input_data_;
}

bool MatcherTask::filter_upstream_residue_collisions() const
{
	return filter_upstream_residue_collisions_;
}

bool MatcherTask::filter_upstream_collisions_by_score() const
{
	return filter_upstream_collisions_by_score_;
}

MatcherTask::Real
MatcherTask::upstream_residue_collision_tolerance() const
{
	return upstream_residue_collision_tolerance_;
}

MatcherTask::Real MatcherTask::upstream_residue_collision_score_cutoff() const
{
	return upstream_residue_collision_score_cutoff_;
}

MatcherTask::Real MatcherTask::upstream_residue_collision_Wfa_atr() const
{
	return upstream_residue_collision_Wfa_atr_;
}

MatcherTask::Real MatcherTask::upstream_residue_collision_Wfa_rep() const
{
	return upstream_residue_collision_Wfa_rep_;
}

MatcherTask::Real MatcherTask::upstream_residue_collision_Wfa_sol() const
{
	return upstream_residue_collision_Wfa_sol_;
}

bool MatcherTask::define_match_by_single_downstream_positioning() const
{
	return define_match_by_single_downstream_positioning_;
}

bool MatcherTask::filter_upstream_downstream_collisions() const
{
	return filter_upstream_and_downstream_residue_collisions_;
}

bool MatcherTask::filter_upstream_downstream_collisions_by_score() const
{
	return filter_upstream_and_downstream_collisions_by_score_;
}

MatcherTask::Real
MatcherTask::upstream_downstream_atom_collision_tolerance() const
{
	return upstream_downstream_atom_collision_tolerance_;
}

MatcherTask::Real
MatcherTask::upstream_downstream_residue_collision_score_cutoff() const
{
	return upstream_downstream_residue_collision_score_cutoff_;
}

MatcherTask::Real
MatcherTask::upstream_downstream_residue_collision_Wfa_atr() const
{
	return upstream_downstream_residue_collision_Wfa_atr_;
}

MatcherTask::Real
MatcherTask::upstream_downstream_residue_collision_Wfa_rep() const
{
	return upstream_downstream_residue_collision_Wfa_rep_;
}

MatcherTask::Real
MatcherTask::upstream_downstream_residue_collision_Wfa_sol() const
{
	return upstream_downstream_residue_collision_Wfa_sol_;
}

bool
MatcherTask::build_round1_hits_twice() const {
	return build_round1_hits_twice_;
}

void
MatcherTask::validate_downstream_orientation_atoms() const
{
	/// STUBBED OUT
}


/// @details This function reads a file specified on the command line through the
/// -match::grid_boundary <filename> flag.  The "grid boundary" specifies a bounding
/// box inside which the third orientation atom of the downstream partner must lie
/// for a downstream conformation to become a match.  The orientation atoms are
/// described in the match::orientation_atoms flag description.
/// The format for the file describing this bounding box is the first four lines
/// of the gridlig file format from rosetta++.  The same file may be listed as
/// both the grid_boundary file and the active_site_definition_by_gridlig file.
/// Though the active_site_definition_by_gridlig flag is not required for the matcher to
/// work, the grid_boundary flag is necessary.
///
/// Format:
/// The gridlig format has been coopted.  This file specifies voxels for a region of 3D.
/// Though the occupied space bounding box does not require a discritization for its
/// description, this gridlig file format does.
///
/// Gridlig format
/// The gridlig format specifies a discretization of euclidean space, and then gives
/// the 0 and 1 values for each voxel.
/// The first line in the liggrid file is "NAME: gridlig"
/// The second line describes the coordinate of the lower corner:
/// "BASE: <xcoord> <ycoord> <zcoord>
/// The third line describes the number of voxels in each dimension:
/// "SIZE: <xcount> <ycount> <zcount>
/// The fourth line describes the length of each dimension describing a voxel:
/// "LENGTH: <xlength> <ylength> <zlength>
/// Only the first four lines of the gridlig file are read.  Additional lines
/// are ignored.
/// Example:
/// <begin file>
/// NAME: gridlig
/// BASE:   27.712   15.827   39.155
/// SIZE:   39   39   34
/// LENGTH:    0.500    0.500    0.500
/// <end file>
///
/// NOTE: The discritization of the occupied space bounding box has no impact
/// on the hash size when detecting that two hits match.  The
/// -match::euclid_bin_size <Real> flag alone controls that behavior.
void
MatcherTask::initialize_occupied_space_bounding_box_from_command_line()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::match;

	if ( option[ grid_boundary ].user() ) {

		std::string filename = option[ grid_boundary ]();
		utility::io::izstream istr( filename.c_str() );
		std::string name, liggrid;
		istr >> name >> liggrid;
		runtime_assert( name == "NAME:" || name == "@NAME:" ); // for some odd reason, the scaffold library has extra "@"'s.
		runtime_assert( liggrid == "gridlig" );
		std::string base;
		Real xbase( 0.0 ), ybase( 0.0 ), zbase( 0.0 );
		istr >> base;
		runtime_assert( base == "BASE:" );
		istr >> xbase; runtime_assert( ! istr.bad() );
		istr >> ybase; runtime_assert( ! istr.bad() );
		istr >> zbase; runtime_assert( ! istr.bad() );

		std::string size;
		istr >> size;
		runtime_assert( size == "SIZE:" );
		Size xsize( 0 ), ysize( 0 ), zsize( 0 );
		istr >> xsize; runtime_assert( ! istr.bad() );
		istr >> ysize; runtime_assert( ! istr.bad() );
		istr >> zsize; runtime_assert( ! istr.bad() );

		runtime_assert( xsize != 0 );
		runtime_assert( ysize != 0 );
		runtime_assert( zsize != 0 );

		std::string length;
		istr >> length;
		runtime_assert( length == "LENGTH:");
		Real xwidth( 0.0 ), ywidth( 0.0 ), zwidth( 0.0 );

		istr >> xwidth; runtime_assert( ! istr.bad() );
		istr >> ywidth; runtime_assert( ! istr.bad() );
		istr >> zwidth; runtime_assert( ! istr.bad() );

		runtime_assert( xwidth != 0 );
		runtime_assert( ywidth != 0 );
		runtime_assert( zwidth != 0 );

		Vector lower_corner( xbase, ybase, zbase );
		Vector upper_corner( lower_corner );
		upper_corner.x() += xwidth * xsize;
		upper_corner.y() += ywidth * ysize;
		upper_corner.z() += zwidth * zsize;
		occ_space_bounding_box_.set_lower( lower_corner );
		occ_space_bounding_box_.set_upper( upper_corner );
	} else {
		if ( !upstream_pose_ ) utility_exit_with_message( "Grid boundary defintion file must be specified with the -match::grid_boundary <filename> flag, or the upstram_pose_ must be set in the MatcherTask." );

		TR << "WARNING WARNING WARNING: no grid file specified for option -grid_boundary. The bounding grid will be generated from the dimensions of the pose. This is experimental at the moment." << std::endl;

		Real lowx(upstream_pose_->residue(1).xyz(1).x());
		Real lowy(upstream_pose_->residue(1).xyz(1).y());
		Real lowz(upstream_pose_->residue(1).xyz(1).z());
		Real highx(lowx), highy(lowy), highz(lowz);

		for ( Size i =1; i <= upstream_pose_->size(); ++i ) {
			core::conformation::Residue const & cur_res( upstream_pose_->residue( i ) );

			if ( !cur_res.is_protein() ) continue;
			utility::vector1< Size > const & bb_atoms( cur_res.type().all_bb_atoms() );

			for ( Size j =1; j <= bb_atoms.size(); ++j ) {
				Real curx( cur_res.xyz( bb_atoms[j] ).x() ), cury( cur_res.xyz( bb_atoms[j] ).y() ), curz( cur_res.xyz( bb_atoms[j] ).z() );

				if (  curx  < lowx ) lowx = curx;
				else if ( curx > highx ) highx = curx;
				if (  cury < lowy ) lowy = cury;
				else if ( cury > highy ) highy = cury;
				if (  curz < lowz ) lowz = curz;
				else if ( curz > highz ) highz = curz;
			}
		}
		//lowx -= padding, lowy -= padding, lowz -= padding;
		//highx += padding, highy += padding, highz += padding;


		Vector lower_corner( lowx , lowy, lowz );
		Vector upper_corner( highx, highy, highz );
		TR << "Bounding box lower corner set to (" << lower_corner.x() << "," << lower_corner.y() << "," << lower_corner.z() << "), upper corner set to to (" << upper_corner.x() << "," << upper_corner.y() << "," << upper_corner.z() << ")." << std::endl;
		occ_space_bounding_box_.set_lower( lower_corner );
		occ_space_bounding_box_.set_upper( upper_corner );
	}
}

/// @details flo may '12, bulk of this function has been moved
/// to AllowedSeqposForGeomCst file
void
MatcherTask::initialize_scaffold_active_site_residue_list_from_command_line()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::match;

	if ( ! option[ scaffold_active_site_residues ].user() && ! option[ scaffold_active_site_residues_for_geomcsts ].user() && ( generic_pose_build_resids_.size() == 0 ) && (per_cst_pose_build_resids_.size() == 0 ) ) {
		utility_exit_with_message( "Scaffold build point defintion file must be specified with either "
			"the -match::scaffold_active_site_residues <filename> flag\n"
			"or the -match::scaffold_active_site_residues_for_geomcsts <filename> flag" );
	}

	protocols::toolbox::match_enzdes_util::AllowedSeqposForGeomCst scaf_seqpos;
	scaf_seqpos.initialize_from_command_line( upstream_pose_ );

	//now go through and initialize the MatcherTask data from the file reader
	//first check whether there's actually any positions to go from
	if ( scaf_seqpos.num_seqpos_lists() == 0 ) utility_exit_with_message("Apparently match positions couldn't get read from file...");

	else if ( scaf_seqpos.num_seqpos_lists() == 1 ) {
		share_build_points_for_geomcsts_ = true;
		generic_pose_build_resids_.clear();
		for ( core::Size i = 1; i <= scaf_seqpos.seqpos_for_geomcst( 1 ).size(); ++i ) {
			if ( scaf_seqpos.seqpos_for_geomcst( 1 )[i] <= upstream_pose_->size() ) {
				if ( upstream_pose_->residue_type(  scaf_seqpos.seqpos_for_geomcst( 1 )[i] ).is_protein() ) generic_pose_build_resids_.push_back( scaf_seqpos.seqpos_for_geomcst( 1 )[i] );
			} else std::cerr << "Warning: seqpos number " << scaf_seqpos.seqpos_for_geomcst( 1 )[i] << " is larger than size of pose, ignoring. Check if you're using the right match posfile.";

		}
	} else {
		if ( scaf_seqpos.num_seqpos_lists() != enz_input_data_->mcfi_lists_size() ) {
			utility_exit_with_message("#geometric constraints disagreement between file given for option match::scaffold_active_site_residues_for_geomcsts having " + utility::to_string(scaf_seqpos.num_seqpos_lists() ) + "constraitns and Enz .cst file: " + option[ geometric_constraint_file ]()() +
				" having (" + utility::to_string(enz_input_data_->mcfi_lists_size()) + ") constraints." );
		}

		share_build_points_for_geomcsts_ = false;
		per_cst_pose_build_resids_.clear();
		per_cst_pose_build_resids_.resize( scaf_seqpos.num_seqpos_lists() );

		for ( core::Size i =1; i <= scaf_seqpos.num_seqpos_lists(); ++i ) {
			for ( core::Size j = 1; j <= scaf_seqpos.seqpos_for_geomcst( i ).size(); ++j ) {
				if ( scaf_seqpos.seqpos_for_geomcst( i )[j] <= upstream_pose_->size() ) {
					if ( upstream_pose_->residue_type(  scaf_seqpos.seqpos_for_geomcst( i )[j] ).is_protein() ) per_cst_pose_build_resids_[i].push_back( scaf_seqpos.seqpos_for_geomcst( i )[j] );
				} else std::cerr <<  "Warning: seqpos number " << scaf_seqpos.seqpos_for_geomcst( i )[j] << " is larger than size of pose, ignoring. Check if you're using the right match posfile.";
			}
		}
	}
}

/// @details go through all of the remarks in the pose and see if any of them specify
/// already existing catalytic interactions
void
MatcherTask::set_active_site_residue_list_to_preexisting_partial_match()
{
	if ( !upstream_pose_ ) return;

	core::pose::PDBInfoCOP pose_pdbinfo = upstream_pose_->pdb_info();

	core::io::Remarks const & pose_remarks = pose_pdbinfo->remarks();

	bool change_build_point_lists( false );
	bool switch_to_different_build_points_required( share_build_points_for_geomcsts_ );
	utility::vector1< Size > generic_pose_build_resids_copy( generic_pose_build_resids_ );
	std::set< Size > cst_ids_already_present;
	Size n_geometric_constraints( enz_input_data_->mcfi_lists_size() );

	// Loop over all the remarks in pdb file
	for ( auto const & pose_remark : pose_remarks ) {

		std::string resA_type(""), resB_type("");
		int resA_num(0), resB_num(0);
		Size cst_block(0);
		std::string resA_chain ,resB_chain;
		core::Size ex_geom_id;
		if ( protocols::toolbox::match_enzdes_util::split_up_remark_line(pose_remark.value, resA_chain,
				resA_type, resA_num, resB_chain, resB_type,
				resB_num, cst_block, ex_geom_id)  ) {

			if ( !change_build_point_lists ) {
				change_build_point_lists = true;
				if ( switch_to_different_build_points_required ) use_different_build_points_for_each_geometric_constraint( n_geometric_constraints );
			}

			Size upstream_seqpos(  pose_pdbinfo->pdb2pose(resB_chain.c_str()[0],resB_num) );
			utility::vector1< Size > res_id;
			runtime_assert( upstream_pose_->residue_type( upstream_seqpos ).is_protein() );
			TR << "An interaction for geometric constraint " << cst_block << " already seems to be present in the pose at seqpos " << upstream_seqpos << ". Matching for this geomcst will only be done at this position." << std::endl;
			res_id.push_back( upstream_seqpos );
			set_original_scaffold_build_points_for_geometric_constraint(cst_block,res_id);
			cst_ids_already_present.insert( cst_block );
		}
	}
	if ( change_build_point_lists ) {
		//it is likely that the downstream object is part of the upstream_pose_
		//in this case, bc the starting structure already contained a partial match.
		remove_downstream_object_from_upstream_pose();

		if ( switch_to_different_build_points_required ) {
			for ( Size i = 1; i <= n_geometric_constraints; ++i ) {
				if ( cst_ids_already_present.find( i ) == cst_ids_already_present.end() ) {
					set_original_scaffold_build_points_for_geometric_constraint( i, generic_pose_build_resids_copy );
				}
			}
		}
	}
}  //set_active_site_residue_list_to_preexisting_partial_match()

void
MatcherTask::remove_downstream_object_from_upstream_pose()
{
	//right now this function only works for ligand downstream objects
	if ( (downstream_pose_->size() != 1 ) ) utility_exit_with_message("Can't remove a downstream pose containing more than one residue from the upstream pose.");

	utility::vector1< Size > seqpos_to_remove;
	for ( core::Size i = 1; i<= upstream_pose_->conformation().num_chains(); ++i ) {
		core::Size chain_begin( upstream_pose_->conformation().chain_begin( i ) );
		if ( (upstream_pose_->conformation().chain_end( i ) - chain_begin == 0 ) &&
				downstream_pose_->residue_type(1).name3() == upstream_pose_->residue_type( chain_begin ).name3() ) {
			seqpos_to_remove.push_back( chain_begin );
		}
	}
	if ( seqpos_to_remove.size() == 0 ) return;

	core::pose::PoseOP mod_up_pose( new core::pose::Pose( *upstream_pose_ ) );
	for ( core::Size i = seqpos_to_remove.size(); i >= 1; --i ) {
		TR << "Removing seqpos " << seqpos_to_remove[i] << " from upstream_pose_ because it is part of the downstream object." << std::endl;
		mod_up_pose->conformation().delete_residue_slow( seqpos_to_remove[i] );
	}
	upstream_pose_ = mod_up_pose;
}

void
MatcherTask::initialize_enzdes_input_data_from_command_line()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::match;

	using namespace protocols::toolbox::match_enzdes_util;
	using namespace core::chemical;
	/// create a local non-const version of the input data.
	EnzConstraintIOOP enz_input_data( new EnzConstraintIO(
		ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) ) );
	if ( ! option[ geometric_constraint_file ].user() ) {
		utility_exit_with_message( "Option match::geometric_constraint_file must be specified on the command line" );
	}
	cstfile_name_ =  option[ geometric_constraint_file ]();

	enz_input_data->read_enzyme_cstfile( cstfile_name_ );
	enz_input_data_ = enz_input_data;
	determine_all_match_relevant_downstream_atoms();

}


/// @details note: for the time being, this function sets all the residue numbers of
/// the AtomIDs in relevant_downstream_atoms_ to 1, i.e. it will only function properly
/// for cases where the downstream pose consists of one residue only (i.e. ligands )
void
MatcherTask::determine_all_match_relevant_downstream_atoms()
{
	runtime_assert( enz_input_data_ != nullptr );
	relevant_downstream_atoms_.clear();
	std::set< core::id::AtomID > seen_atoms;

	for ( core::Size i = 1; i <= enz_input_data_->num_mcfi_lists(); ++i ) {

		std::map< std::string, utility::vector1< std::string > > const &
			alg_info(  enz_input_data_->mcfi_list( i )->mcfi( 1 )->algorithm_inputs() );
		if ( alg_info.find( "match" ) != alg_info.end() ) {
			//make sure this is not upstream matching
			bool upstream_matching( false);
			utility::vector1< std::string > const & info( alg_info.find( "match" )->second );
			for ( Size ll = 1; ll <= info.size(); ++ll ) {
				std::string llstr = info[ ll ];
				std::istringstream llstream( llstr );
				std::string first, second;
				llstream >> first;
				if ( first == "SECONDARY_MATCH:" ) {
					llstream >> second;
					if ( second == "UPSTREAM_CST" ) {
						core::Size target_id(0);
						llstream >> target_id;
						upstream_matching = true;
						upstream_only_geom_cst_.insert( std::pair< core::Size, core::Size >( i, target_id ));
						break;
					}
				}
			}
			if ( upstream_matching ) continue;
		}
		//upstream matching check over

		for ( core::Size j = 1; j <= enz_input_data_->mcfi_list( i )->num_mcfis(); ++j ) {
			protocols::toolbox::match_enzdes_util::MatchConstraintFileInfoCOP cur_mcfi( enz_input_data_->mcfi_list( i )->mcfi( j ) );
			core::Size ds_matchres( cur_mcfi->downstream_res() );
			if ( cur_mcfi->allowed_restypes( ds_matchres ).size() == 0 ) {
				utility_exit_with_message("No allowed residue types seen for downstream residue for constraint " + utility::to_string( i ) + " block " + utility::to_string( j ) );
			}
			core::chemical::ResidueTypeCOP ds_restype( cur_mcfi->allowed_restypes( ds_matchres )[1] );

			for ( core::Size k = 1; k <= 3; ++k ) {

				for ( core::Size l = 1; l <= cur_mcfi->template_atom_inds( ds_matchres, k, *ds_restype ).size(); ++l ) {
					core::id::AtomID this_at( cur_mcfi->template_atom_inds( ds_matchres, k, *ds_restype )[l], 1 );
					if ( seen_atoms.find( this_at ) == seen_atoms.end() ) {
						seen_atoms.insert( this_at );
						relevant_downstream_atoms_.push_back( this_at );
						//std::cerr << " atom " << cur_mcfi->template_atom_inds( ds_matchres, k, *ds_restype )[l] << " is relevant for matcher" << std::endl;
					}
				} //loop over template atom inds
			} //loop over template atoms
		} //loop over mcifs
	} //loop over mcfi lists
}

void
MatcherTask::initialize_orientation_atoms_from_command_line()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::match;

	if ( option[ orientation_atoms ].user() ) {
		if ( downstream_pose_ ) {
			utility::vector1< std::string > names = option[ orientation_atoms ];
			if ( names.size() != 3 ) {
				std::cerr << "ERROR: expected exactly three atom names for the downstream pose, but read ";
				std::cerr << names.size() << " from the command line." << std::endl;
				for ( Size ii = 1; ii <= names.size(); ++ii ) {
					std::cerr << names[ ii ] << " ";
				}
				std::cerr << std::endl;
				utility_exit_with_message( "Invalid use of the flag -match::orientation_atoms" );
			}

			if ( downstream_pose_->size() != 1 ) {
				std::cerr << "ERROR: Cannot use the flag -match::orientation_atoms if the downstream pose has multiple residues" << std::endl;
				utility_exit_with_message( "Invalide use of the flag -match::orientation_atoms" );
			}
			utility::vector1< Size > ats( 3 );
			for ( Size ii = 1; ii <= 3; ++ii ) {
				if ( ! downstream_pose_->residue( 1 ).has( names[ ii ] ) ) {
					std::cerr << "Could not find atom named '" << names[ ii ] << "' in residue " << downstream_pose_->residue( 1 ).name() << std::endl;
					utility_exit_with_message ("Unrecognized atom name listed in flag -match::orientation_atoms" );
				} else {
					ats[ ii ] = downstream_pose_->residue( 1 ).atom_index( names[ ii ] );
					downstream_orientation_atoms_[ ii ] = core::id::AtomID( ats[ ii ], 1 );
				}
			}

			Real angle = numeric::constants::d::radians_to_degrees * numeric::angle_radians(
				downstream_pose_->residue( 1 ).xyz( ats[ 1 ] ),
				downstream_pose_->residue( 1 ).xyz( ats[ 2 ] ),
				downstream_pose_->residue( 1 ).xyz( ats[ 3 ] ) );
			if ( angle > 175 ) {
				std::cerr << "ERROR: The three atoms used to define the orientation of the downstream partner are nearly co-linear: angle= " << angle << "\n";
				std::cerr << "Please choose atoms with an angle that will not be so numerically sensitive." << std::endl;
				utility_exit_with_message( "Nearly-colinear points chosen as orientation atoms of the downstream partner" );
			}

		} else {
			std::cout << "Warning, Downstream Pose has not been set in MatcherTask when initialize_from_command_line\n";
			std::cout << "was invoked.  The parameters in the command-line flag -match::orientation_atoms are being ignored" << std::endl;
		}
	}

}

/// @details This function reads three types of files relating to the active site.
/// The file formats are described below.
/// 1. A file listing the atom names of those atoms in the downstream partner
/// which are required to be in the active site -- conformations of the downstream
/// partner in which any of these atoms lie outside of the active site are excluded.
/// This file format is only valid if the downstream partner is a single residue
/// (e.g. a ligand).
/// 2. A file defining the active site geometry by listing residue indices of the
/// scaffold and, for each residue i, a reach, r_i (a radius).  A point in space is
/// within the active site it is within r_i of any atom of residue i.
/// 3. A "gridlig" file from rosetta++ where a box in space whose voxels are labeled
/// with 0's and 1's; voxels with 1's are part of the active site, those with 0's are
/// not.
///
/// File formats:
/// 1.  Active site atoms.
/// The file lists on one or more lines the atom names for all atoms that are required
/// to be present in the active site.  The file format does not support comments.
/// Example:
/// <begin file>
/// C7 O2 O3 C22 C23 N5
/// <end file>
///
/// 2.  Residue/reach pairs.
/// Each line of the file should list a residue index (where residue numbering starts
/// from 1) for a residue on the scaffold and a reach (a radius).
/// Example.
/// <begin file>
/// 10 5.5
/// 12 7.2
/// 13 7.0
/// 23 10.0
/// 25 7.0
/// 58 5.5
/// 60 7.0
/// <end file>
///
/// 3. Gridlig format
/// The gridlig format specifies a discretization of euclidean space, and then gives
/// the 0 and 1 values for each voxel.
/// The first line in the liggrid file is "NAME: gridlig"
/// The second line describes the coordinate of the lower corner:
/// "BASE: <xcoord> <ycoord> <zcoord>
/// The third line describes the number of voxels in each dimension:
/// "SIZE: <xcount> <ycount> <zcount>
/// The fourth line describes the length of each dimension describing a voxel:
/// "LENGTH: <xlength> <ylength> <zlength>
/// The subsequent lines are arranged in blocks sweeping through the grid in
/// "row-major" order.  Each line has zcount space-separated 0's and 1's.  Lines
/// are grouped into ycount blocks; a single empty line separates blocks.
/// Example:
/// <begin file>
/// NAME: gridlig
/// BASE:   27.712   15.827   39.155
/// SIZE:   39   39   34
/// LENGTH:    0.500    0.500    0.500
/// 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
/// 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
/// 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
/// 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1
/// <the file contains 1563 lines total>
void
MatcherTask::initialize_active_site_definition_from_command_line()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::match;

	if ( option[ required_active_site_atom_names ].user() ) {
		if ( downstream_pose_->size() != 1 ) {
			utility_exit_with_message( "Flag match::required_active_site_atom_names may only be used when the downstream pose is a single residue" );
		}
		if (
				! option[ active_site_definition_by_residue ].user() &&
				! option[ active_site_definition_by_gridlig ].user() ) {
			utility_exit_with_message( "Flag match::required_active_site_atom_names must be used in combination with\n" "either the match::active_site_definition_by_residue flag or the match::active_site_definition_by_gridlig flag" );
		}
		std::string filename = option[ required_active_site_atom_names ]();
		utility::io::izstream istr( filename.c_str() );
		while ( istr ) {
			std::string atname;
			istr >> atname;
			if ( ! istr.bad() ) {
				if ( atname == "" ) break;
				if ( downstream_pose_->residue( 1 ).has( atname ) ) {
					TR << "Requiring downstream atom '" << atname << "' to reside in the scaffold's active site" << std::endl;
					downstream_atoms_required_inside_active_site_.emplace_back( downstream_pose_->residue( 1 ).atom_index( atname ), 1 );
				} else {
					std::cerr << "ERROR reading required_active_site_atom_names: atom named '" << atname << "' is not an atom of " << downstream_pose_->residue( 1 ).name() << std::endl;
					utility_exit_with_message( "Problem reading required_active_site_atom_names file " + filename );
				}
			}
		}
		if ( downstream_atoms_required_inside_active_site_.empty() ) {
			utility_exit_with_message( "No active site atoms read; empty required_active_site_atom_names file " + filename );
		}
	}

	if ( option[ active_site_definition_by_residue ].user() ) {
		if ( option[ active_site_definition_by_gridlig ].user() ) {
			std::cerr << "ERROR: found incompatible flags active_site_definition_by_residue and active_site_definition_by_gridlig on the command line" << std::endl;
			utility_exit_with_message( "Ambiguous command line." );
		}
		std::string filename = option[ active_site_definition_by_residue ]();
		utility::io::izstream istr( filename.c_str() );
		while ( istr ) {
			Size resid; Real radius;
			istr >> resid;
			if ( resid > upstream_pose_->size() ) {
				std::cerr << "ERROR reading active_site_definition: residue " << resid << " exceeds the number of residues in the scaffold pose." << std::endl;
				utility_exit_with_message( "Problem reading active_site_definition " + filename );
			}
			if ( ! istr.good() ) break;
			istr >> radius;
			if ( istr.bad() ) {
				std::cerr << "Error reading active_site_definition: expected to read a radius after reading resid " << resid << std::endl;
				utility_exit_with_message( "Problem reading active_site_definition " + filename );
			}
			TR << "Including sphere of radius " << radius << " surrounding scaffold residue " << resid << " in active site definition" << std::endl;
			upstream_resids_and_radii_defining_active_site_.push_back( std::make_pair( resid, radius ));
		}
		if ( upstream_resids_and_radii_defining_active_site_.empty() ) {
			if ( ! downstream_atoms_required_inside_active_site_.empty() ) {
				/// No matches possible as the active site is undefined, yet the
				/// there are atoms required to reside in the active site.  Quit.
				std::cerr << "ERROR reading active_site_definition: no active site residues were read!" << std::endl;
				utility_exit_with_message( "Problem reading active_site_definition " + filename );
			}
		}
	}

	if ( option[ active_site_definition_by_gridlig ].user() ) {
		gridlig_active_site_definition_ = true;
		gridlig_fname_ = option[ active_site_definition_by_gridlig ]();
	}
}

void
MatcherTask::initialize_upstream_residue_collision_filter_data_from_command_line()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::match;

	if ( option[ OptionKeys::match::filter_colliding_upstream_residues ] ) {
		filter_upstream_residue_collisions_ = true;
		if ( option[ OptionKeys::match::upstream_residue_collision_score_cutoff ].user() ) {
			filter_upstream_collisions_by_score_ = true;
			upstream_residue_collision_score_cutoff_ = option[ OptionKeys::match::upstream_residue_collision_score_cutoff ];
			upstream_residue_collision_Wfa_atr_ = option[ OptionKeys::match::upstream_residue_collision_Wfa_atr ];
			upstream_residue_collision_Wfa_rep_ = option[ OptionKeys::match::upstream_residue_collision_Wfa_rep ];
			upstream_residue_collision_Wfa_sol_ = option[ OptionKeys::match::upstream_residue_collision_Wfa_sol ];
		} else {
			filter_upstream_collisions_by_score_ = false;
			if ( option[ OptionKeys::match::upstream_residue_collision_tolerance ].user() ) {
				upstream_residue_collision_tolerance_ = option[ OptionKeys::match::upstream_residue_collision_tolerance ];
			} else {
				upstream_residue_collision_tolerance_ = option[ OptionKeys::match::bump_tolerance ];
			}
		}
	}
}

void
MatcherTask::initialize_upstream_downstream_collision_filter_data_from_command_line()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::match;

	if ( option[ OptionKeys::match::filter_upstream_downstream_collisions ] ) {
		filter_upstream_and_downstream_residue_collisions_ = true;
		if ( option[ OptionKeys::match::updown_residue_collision_score_cutoff ].user() ) {
			filter_upstream_and_downstream_collisions_by_score_ = true;
			upstream_downstream_residue_collision_score_cutoff_ = option[ OptionKeys::match::updown_residue_collision_score_cutoff ];
			upstream_downstream_residue_collision_Wfa_atr_ = option[ OptionKeys::match::updown_residue_collision_Wfa_atr ];
			upstream_downstream_residue_collision_Wfa_rep_ = option[ OptionKeys::match::updown_residue_collision_Wfa_rep ];
			upstream_downstream_residue_collision_Wfa_sol_ = option[ OptionKeys::match::updown_residue_collision_Wfa_sol ];
		} else {
			filter_upstream_and_downstream_collisions_by_score_ = false;
			if ( option[ OptionKeys::match::updown_collision_tolerance ].user() ) {
				upstream_downstream_atom_collision_tolerance_ = option[ OptionKeys::match::updown_collision_tolerance ];
			} else {
				upstream_downstream_atom_collision_tolerance_ = option[ OptionKeys::match::bump_tolerance ];
			}
		}
	}
}

void
MatcherTask::initialize_output_options_from_command_line()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys::match;

	output_writer_name_ = option[ OptionKeys::match::output_format ];
	grouper_name_ = option[ OptionKeys::match::match_grouper ];

	grouper_ds_rmsd_ = option[ OptionKeys::match::grouper_downstream_rmsd ];
	output_matchres_only_ =  option[ OptionKeys::match::output_matchres_only ];

	// Check if score output option is enabled
	if ( option[ OptionKeys::out::file::scorefile ].user() ) {
		std::string scorefile_path = option[ OptionKeys::out::file::scorefile ];
		if ( scorefile_path == "" ) {
			output_scores_ = true;
		} else {
			score_output_file_name( option[ OptionKeys::out::file::scorefile ] );
		}
	}

	utility::vector1< int > tempvec = option[ OptionKeys::match::geom_csts_downstream_output ];
	if ( tempvec.size() == 0 ) utility_exit_with_message("Bad user input: empty vector specified for option -match::geom_csts_downstream_output.");
	geom_csts_downstream_output_.clear();

	core::Size num_geom_csts( enz_input_data_->num_mcfi_lists() );

	for ( core::Size i = 1; i <= tempvec.size(); ++i ) {

		if ( tempvec[ i ]  < 1 ) {
			utility_exit_with_message("Retarded user input. Output for geom cst with id smaller than 1 requested.");
		}

		if ( (Size) tempvec[ i ] > num_geom_csts ) {
			utility_exit_with_message("Bad user input. Output for geom cst with id higher than the number of total geomcsts requested.");
		}

		geom_csts_downstream_output_.push_back( (core::Size ) tempvec[ i ] );
	}

	if ( option[ OptionKeys::match::define_match_by_single_downstream_positioning ] ) {
		if ( geom_csts_downstream_output_.size() != 1 ) utility_exit_with_message("Bad user input: values specified for options -match::define_match_by_single_downstream_positioning and -match::geom_csts_downstream_output are incompatible.");
		if ( ! option[ OptionKeys::match::consolidate_matches ] ) {
			define_match_by_single_downstream_positioning_ = true;
		} else {
			define_match_by_single_downstream_positioning_ = false;
		}
	} else {
		if ( option[ OptionKeys::match::output_format ]() == "PDB" &&
				! option[ OptionKeys::match::consolidate_matches ] && (geom_csts_downstream_output_.size() == 1 ) ) {
			TR << "MatcherTask initialised to treat matches as defined by single downstream positioning." << std::endl;
			define_match_by_single_downstream_positioning_ = true;
		} else {
			define_match_by_single_downstream_positioning_ = false;
		}
	}
}

}
}
