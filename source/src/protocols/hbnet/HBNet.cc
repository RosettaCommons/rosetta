// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/hbnet/HBNet.cc
/// @brief base class HBNet; to explicitly detect and design h-bond networks
/// @author Scott Boyken (sboyken@gmail.com)

// Header files for project
#include <protocols/hbnet/HBNet.hh>
#include <protocols/hbnet/HBNetCreator.hh>
#include <protocols/hbnet/HBNet_util.hh>

// Utility headers
#include <utility/io/ozstream.hh>
#include <utility/numbers.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>
#include <utility>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>

// Core headers
#include <core/chemical/AtomType.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>

#include <core/id/AtomID.hh>
#include <core/io/Remarks.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pack/hbonds/HBondGraph_util.hh>
#include <core/pack/hbonds/MCHBNetInteractionGraph.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/InteractionGraphFactory.hh>
#include <core/pack/interaction_graph/PDInteractionGraph.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/pose/PDBInfo.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/selection.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/graph/AtomLevelHBondGraph.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/util.hh>

#include <core/select/residue_selector/LayerSelector.hh>
#include <core/select/residue_selector/ResidueNameSelector.hh>
#include <core/select/residue_selector/ResiduePDBInfoHasLabelSelector.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/CachedResidueSubset.hh>
#include <core/select/util.hh>

// Protocols headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <protocols/pose_creation/MakePolyXMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/jd2/util.hh>
#include <protocols/residue_selectors/StoreResidueSubsetMover.hh>
#include <protocols/moves/mover_schemas.hh>

#include <numeric/random/random.hh>

//Currently Unused
//#include <core/chemical/ResidueTypeSet.hh>
//#include <core/conformation/ResidueFactory.hh>
//#include <core/pack/pack_rotamers.hh>
//#include <core/pack/rotamer_set/RotamerSetFactory.hh>
//#include <core/pack/rotamer_set/RotamerSet_.hh>
//#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSet_.hh>
//#include <core/scoring/facts/FACTSEnergy.hh>
//#include <core/scoring/hbonds/HBEvalTuple.hh>
//#include <protocols/enzdes/AddorRemoveCsts.hh>
//#include <protocols/scoring/Interface.hh>

using namespace core;
using namespace pose;
using namespace pack;
using namespace rotamer_set;
using namespace interaction_graph;
using namespace conformation;
using namespace scoring::hbonds;
using namespace core::pack::hbonds;
using namespace core::scoring::hbonds::graph;
using namespace utility;

namespace protocols {
namespace hbnet {

static basic::Tracer TR( "protocols.hbnet.HBNet" );

HBNet::HBNet( ) :
	moves::Mover( "HBNet" ),
	write_network_pdbs_(false),
	write_cst_files_(true),
	output_poly_ala_background_(false),
	store_network_scores_in_pose_( false ),
	find_native_(false),
	only_native_(false),
	keep_existing_networks_(false),
	extend_existing_networks_(false),
	only_extend_existing_(false),
	verbose_(false),
	symmetric_(false),
	multi_component_(false),
	show_task_(false),
	minimize_(true),
	//bridging_waters_(0),
	allow_no_hbnets_(false),
	tyr_hydroxyls_must_donate_(false),
	hydroxyls_must_donate_(false),
	use_pdb_numbering_(true),
	no_heavy_unsats_allowed_(true),
	keep_start_selector_rotamers_fixed_(false),
	min_network_size_(3),
	max_network_size_(15),
	min_unique_networks_(1),
	min_core_res_(0),
	min_boundary_res_(0),
	max_unsat_Hpol_(5),
	max_lig_unsat_(15),
	max_rep_(1),
	max_replicates_before_branch_(0),
	max_replicates_before_unsat_check_(1),
	//bridging_water_search_depth_(1),
	des_residues_("STRKHYWNQDE"),
	start_res_vec_( /* NULL */ ),
	network_vector_( /* NULL */ ),
	native_networks_( /* NULL */ ),
	merged_vecs_( /* NULL */ ),
	output_vector_(/* NULL */ ),
	hydrogen_bond_threshold_(-0.5),
	onebody_hb_threshold_(-0.4),
	charge_charge_rep_cutoff_(1.0),
	clash_threshold_(1.0),
	upper_score_limit_(15.0),
	min_percent_hbond_capacity_(0.5),
	task_factory_( /* NULL */ ),
	init_scorefxn_(core::scoring::ScoreFunctionFactory::create_score_function( "HBNet" )),
	scorefxn_( core::scoring::get_score_function() ),
	rotamer_sets_( /* NULL */ ),
	ig_( /* NULL */ ),
	rotamer_links_(nullptr),
	store_subnetworks_(false),
	secondary_search_(false),
	secondary_threshold_(-0.25),
	upweight_twobody_(1.0),
	start_selector_( /* NULL */ ),
	core_selector_( /* NULL */ ),
	boundary_selector_( /* NULL */ ),
	core_residues_(0),
	boundary_residues_(0),
	input_hbnet_info_residues_(0)
{
	set_monte_carlo_data_to_default();
}

HBNet::HBNet( std::string const name ) :
	moves::Mover( name ),
	write_network_pdbs_(false),
	write_cst_files_(true),
	output_poly_ala_background_(false),
	store_network_scores_in_pose_( false ),
	find_native_(false),
	only_native_(false),
	keep_existing_networks_(false),
	extend_existing_networks_(false),
	only_extend_existing_(false),
	verbose_(false),
	symmetric_(false),
	multi_component_(false),
	show_task_(false),
	minimize_(true),
	//bridging_waters_(0),
	allow_no_hbnets_(false),
	tyr_hydroxyls_must_donate_(false),
	hydroxyls_must_donate_(false),
	use_pdb_numbering_(true),
	no_heavy_unsats_allowed_(true),
	keep_start_selector_rotamers_fixed_(false),
	min_network_size_(3),
	max_network_size_(15),
	min_unique_networks_(1),
	min_core_res_(0),
	min_boundary_res_(0),
	max_unsat_Hpol_(5),
	max_lig_unsat_(15),
	max_rep_(1),
	max_replicates_before_branch_(0),
	max_replicates_before_unsat_check_(1),
	//bridging_water_search_depth_(1),
	des_residues_("STRKHYWNQDE"),
	start_res_vec_( /* NULL */ ),
	network_vector_( /* NULL */ ),
	native_networks_( /* NULL */ ),
	merged_vecs_( /* NULL */ ),
	output_vector_(/* NULL */ ),
	hydrogen_bond_threshold_(-0.5),
	onebody_hb_threshold_(-0.4),
	charge_charge_rep_cutoff_(1.0),
	clash_threshold_(1.0),
	upper_score_limit_(15.0),
	min_percent_hbond_capacity_(0.5),
	task_factory_( /* NULL */ ),
	init_scorefxn_(core::scoring::ScoreFunctionFactory::create_score_function( "HBNet" )),
	scorefxn_( core::scoring::get_score_function() ),
	rotamer_sets_( /* NULL */ ),
	ig_( /* NULL */ ),
	rotamer_links_(nullptr),
	store_subnetworks_(false),
	secondary_search_(false),
	secondary_threshold_(-0.25),
	upweight_twobody_(1.0),
	start_selector_( /* NULL */ ),
	core_selector_( /* NULL */ ),
	boundary_selector_( /* NULL */ ),
	core_residues_(0),
	boundary_residues_(0),
	input_hbnet_info_residues_(0)
{
	set_monte_carlo_data_to_default();
}

//constructor to be called from code, NEED TO CLEAN THIS UP!
HBNet::HBNet(
	core::scoring::ScoreFunctionCOP scorefxn,
	Size max_unsat_Hpol,
	Size min_network_size, /* 3 */
	Real hb_threshold, /* -0.75 */
	Size max_network_size, /* 15 */
	std::string des_residues, /* "STRKHYWNQDE" */
	bool find_native, /*false*/
	bool only_native, /*false*/
	bool keep_existing, /*false*/
	bool extend_existing, /*false*/
	bool only_extend /*false*/
	//bool bridging_waters, /*false*/
	//bool minimize /*false*/
) :
	moves::Mover( "HBNet" ),
	write_network_pdbs_(false),
	write_cst_files_(true),
	output_poly_ala_background_(false),
	store_network_scores_in_pose_( false ),
	find_native_(find_native),
	only_native_(only_native),
	keep_existing_networks_(keep_existing),
	extend_existing_networks_(extend_existing),
	only_extend_existing_(only_extend),
	verbose_(false),
	symmetric_(false),
	multi_component_(false),
	show_task_(false),
	minimize_(true),
	//bridging_waters_(0),
	allow_no_hbnets_(false),
	tyr_hydroxyls_must_donate_(false),
	hydroxyls_must_donate_(false),
	use_pdb_numbering_(true),
	no_heavy_unsats_allowed_(true),
	keep_start_selector_rotamers_fixed_(false),
	min_network_size_(min_network_size),
	max_network_size_(max_network_size),
	min_unique_networks_(1),
	min_core_res_(0),
	min_boundary_res_(0),
	max_unsat_Hpol_(max_unsat_Hpol),
	max_lig_unsat_(15),
	max_rep_(1),
	max_replicates_before_branch_(0),
	max_replicates_before_unsat_check_(1),
	//bridging_water_search_depth_(1),
	des_residues_(std::move(des_residues)),
	start_res_vec_( /* NULL */ ),
	network_vector_( /* NULL */ ),
	native_networks_( /* NULL */ ),
	merged_vecs_( /* NULL */ ),
	output_vector_(/* NULL */ ),
	hydrogen_bond_threshold_(hb_threshold),
	onebody_hb_threshold_(-0.4),
	charge_charge_rep_cutoff_(1.0),
	clash_threshold_(1.0),
	upper_score_limit_(15.0),
	min_percent_hbond_capacity_(0.5),
	task_factory_( /* NULL */ ),
	init_scorefxn_(core::scoring::ScoreFunctionFactory::create_score_function( "HBNet" )),
	scorefxn_( scorefxn->clone() ),
	rotamer_sets_( /* NULL */ ),
	ig_( /* NULL */ ),
	rotamer_links_(nullptr),
	store_subnetworks_(false),
	secondary_search_(false),
	secondary_threshold_(-0.25),
	upweight_twobody_(1.0),
	start_selector_( /* NULL */ ),
	core_selector_( /* NULL */ ),
	boundary_selector_( /* NULL */ ),
	core_residues_(0),
	boundary_residues_(0),
	input_hbnet_info_residues_(0)
{
	//lkb_ = new core::scoring::lkball::LK_BallEnergy(init_scorefxn_->energy_method_options());
	//lkb_(init_scorefxn_->energy_method_options());
	if ( only_native_ ) find_native_ = true;
	if ( only_extend_existing_ ) extend_existing_networks_ = true;

	set_monte_carlo_data_to_default();

}

//// TODO NEED: Copy constructor to only copy over key settings and parameters;

moves::MoverOP
HBNet::clone() const {
	return pointer::make_shared< HBNet >( *this );
}

moves::MoverOP
HBNet::fresh_instance() const {
	return pointer::make_shared< HBNet >();
}

//destructor
HBNet::~HBNet()= default;

void
HBNet::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &data, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	hydrogen_bond_threshold_ = tag->getOption<core::Real>( "hb_threshold" ,-0.5);
	//bw_cutoff_ = tag->getOption<core::PackerEnergy>("bw_cutoff",-0.5);
	onebody_hb_threshold_ = tag->getOption<core::PackerEnergy>("onebody_hb_threshold",-0.4);
	charge_charge_rep_cutoff_ = tag->getOption<core::PackerEnergy>("charge_charge_rep_cutoff",1.0);
	clash_threshold_ = tag->getOption<core::PackerEnergy>("clash_threshold",1.0);
	upper_score_limit_ = tag->getOption<core::Real>("upper_score_limit",15.0);
	find_native_ = tag->getOption<bool>( "find_native_networks", false );
	only_native_ = tag->getOption<bool>( "find_only_native_networks", false );
	keep_existing_networks_ = tag->getOption<bool>( "keep_existing_networks", false );
	extend_existing_networks_ = tag->getOption<bool>( "extend_existing_networks", false );
	only_extend_existing_ = tag->getOption<bool>( "only_extend_existing", false );
	minimize_ = tag->getOption< bool > ("minimize",true);
	store_subnetworks_ = tag->getOption< bool >("store_subnetworks",false);
	secondary_search_ = tag->getOption< bool >("secondary_search",false);
	secondary_threshold_ = tag->getOption< Real >("secondary_threshold",-0.25);
	write_network_pdbs_ = tag->getOption< bool >( "write_network_pdbs", false );
	write_cst_files_ = tag->getOption< bool >( "write_cst_files", true );
	output_poly_ala_background_ = tag->getOption< bool >( "output_poly_ala_background", false );
	store_network_scores_in_pose_ = tag->getOption< bool >( "store_network_scores_in_pose", false );
	max_rep_ = tag->getOption< Size >( "max_replicates", 1 );
	max_replicates_before_branch_ = tag->getOption< Size >( "max_replicates_before_branch", 0 );
	max_replicates_before_unsat_check_ = tag->getOption< Size >( "max_replicates_before_unsat_check", 1 );
	//    bridging_water_search_depth_ = tag->getOption< Size >( "bridging_water_search_depth", 1 );
	tyr_hydroxyls_must_donate_ = tag->getOption<bool>( "tyr_hydroxyls_must_donate", false ); // only for unsat check -- should move to filter
	hydroxyls_must_donate_ = tag->getOption<bool>( "hydroxyls_must_donate", false ); // only for unsat check -- should move to filter
	use_pdb_numbering_ = tag->getOption<bool>( "use_pdb_numbering", true );
	no_heavy_unsats_allowed_ = tag->getOption<bool>( "no_heavy_unsats_allowed", true );
	keep_start_selector_rotamers_fixed_ = tag->getOption<bool>( "use_only_input_rot_for_start_res", false );
	show_task_ = tag->getOption<bool>( "show_task", false );
	min_network_size_ = tag->getOption<Size>( "min_network_size", 3 );
	max_network_size_ = tag->getOption<Size>( "max_network_size", 15 );
	min_unique_networks_ = tag->getOption<Size>( "min_unique_networks", 1 );
	min_core_res_ = tag->getOption<Size>( "min_core_res", 0 );
	min_boundary_res_ = tag->getOption<Size>( "min_boundary_res", 0 );
	max_unsat_Hpol_ = tag->getOption<Size>( "max_unsat_Hpol", 5 );
	min_percent_hbond_capacity_ = tag->getOption<core::Real>( "min_percent_hbond_capacity", 0.5 );
	allow_no_hbnets_ = tag->getOption< bool >( "allow_no_hbnets", false);

	//MONTE CARLO/////////
	monte_carlo_branch_ = tag->getOption< bool >( "monte_carlo_branch", monte_carlo_branch_ );//deprecated
	monte_carlo_branch_ = tag->getOption< bool >( "monte_carlo", monte_carlo_branch_ );

	max_mc_nets_ = tag->getOption<Size>( "max_mc_nets", max_mc_nets_ );
	total_num_mc_runs_ = tag->getOption<Size>( "total_num_mc_runs", total_num_mc_runs_ );
	monte_carlo_seed_threshold_ = tag->getOption< Real >( "seed_hbond_threshold", monte_carlo_seed_threshold_ );
	monte_carlo_seed_must_be_buried_ = tag->getOption< bool >( "monte_carlo_seed_must_be_buried", monte_carlo_seed_must_be_buried_ );
	monte_carlo_seed_must_be_fully_buried_ = tag->getOption< bool >( "monte_carlo_seed_must_be_fully_buried", monte_carlo_seed_must_be_fully_buried_ );
	///////////////////////

	if ( only_native_ ) find_native_ = true;
	if ( only_extend_existing_ ) extend_existing_networks_ = true;

	verbose_ = tag->getOption<bool>( "verbose", false);

	if ( tag->hasOption("design_residues") ) {
		des_residues_ = tag->getOption<std::string>("design_residues","STRKHYWNQDE");
	}

	// get task operations
	if ( tag->hasOption("task_operations") ) {
		task_factory( rosetta_scripts::parse_task_operations( tag, data ) );
	}

	if ( tag->hasOption("start_selector") ) {
		if ( tag->hasOption("start_resnums" ) ) {
			TR.Warning << "cannot use both start_resnums and start_selector options; start_selector will be used" << std::endl;
		}
		std::string const selector_name ( tag->getOption< std::string >( "start_selector" ) );
		if ( TR.visible() ) TR << "Set selector name to " << selector_name << "." << std::endl;
		core::select::residue_selector::ResidueSelectorCOP selector;
		try {
			selector = data.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name );
		} catch ( utility::excn::Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from AddCompositionConstraintMover::parse_tag()\n" + e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_message );
		}
		runtime_assert( selector );
		start_selector_ = selector->clone();
	} else if ( tag->hasOption("start_resnums") ) {
		start_selector_ = core::pose::get_resnum_selector( tag, "start_resnums" );
	}
	if ( tag->hasOption("core_selector") ) {
		std::string const selector_name ( tag->getOption< std::string >( "core_selector" ) );
		if ( TR.visible() ) TR << "Set selector name to " << selector_name << "." << std::endl;
		core::select::residue_selector::ResidueSelectorCOP selector;
		try {
			selector = data.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name );
		} catch ( utility::excn::Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from AddCompositionConstraintMover::parse_tag()\n" + e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_message );
		}
		runtime_assert( selector );
		core_selector_ = selector->clone();
		//THIS NEEDS TO HAPPEN AT APPLY TIME, NOT PARSE TIME!!!!!
		//core_residues_ = selector->apply( pose );
	}
	if ( tag->hasOption("boundary_selector") ) {
		std::string const selector_name ( tag->getOption< std::string >( "boundary_selector" ) );
		if ( TR.visible() ) TR << "Set selector name to " << selector_name << "." << std::endl;
		core::select::residue_selector::ResidueSelectorCOP selector;
		try {
			selector = data.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name );
		} catch ( utility::excn::Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from AddCompositionConstraintMover::parse_tag()\n" + e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_message );
		}
		runtime_assert( selector );
		boundary_selector_ = selector->clone();
	}
	core::scoring::ScoreFunctionOP new_score_function( rosetta_scripts::parse_score_function( tag, data ) );
	if ( new_score_function == nullptr ) return;
	set_score_function( new_score_function );
}

char
HBNet::get_aa_for_state( Size const res, Size const rot ) const {
	return rotamer_sets_->rotamer_set_for_residue( (platform::uint)(res) )->rotamer(rot)->name1();
}
bool
HBNet::res_is_boundary( Size const res ) const {
	if ( res > boundary_residues_.size() ) return false;
	else if ( boundary_residues_[ res ] ) return true;
	return false;
}

bool
HBNet::res_is_core( Size const res ) const {
	if ( res > core_residues_.size() ) return false;
	else if ( core_residues_[ res ] ) return true;
	return false;
}

///@details Traverse IG to enumerate all possible h-bond networks;
/// networks are stored using the store_networks() function and pushed to the back of network_vector_
void
HBNet::traverse_IG( Real const hb_threshold ){

	// Start the IG traversal at pre-deterimined starting residues (e.g. interface residues if HBNetStapleInterface being called)
	for ( core::Size res1 : start_res_vec_ ) {

		platform::uint const first_ni = rotamer_sets_->resid_2_moltenres(res1); //returns 0 if res1 not molten (not design/repack-able)
		if ( first_ni == 0 ) continue;
		auto const first_node_ind = (int)(first_ni);

		// Loop over all edges in IG for each starting residue
		for ( ig_->reset_edge_list_iterator_for_node(first_node_ind); !ig_->edge_list_iterator_at_end(); ig_->increment_edge_list_iterator() ) {
			EdgeBase const & edge(ig_->get_edge());
			int const second_node_ind = edge.get_other_ind(first_node_ind);
			platform::uint const second_ni = edge.get_other_ind(first_node_ind);
			Size const res2 = rotamer_sets_->moltenres_2_resid(second_ni);
			Size res1_ind(res1),res2_ind(res2);
			// if pose is symmetric, get independent resnums and set score multiply factors
			// Symmetric IG only contains nodes/edges for independent residues; symm clone energies added to ind edges
			if ( symmetric_ ) {
				res1_ind = get_ind_res( *orig_pose_, res1);
				res2_ind = get_ind_res( *orig_pose_, res2);
			}
			Size scmult_2b((symmetric_) ? symm_info_->score_multiply(res1,res2) : 1);
			Size scmult_1b((symmetric_) ? symm_info_->score_multiply_factor() : 1);

			// for all rotamers of the starting position X all rotamers at the other position of edge, look at energies:
			for ( int ii = 1; ii <= ig_->get_num_states_for_node(first_node_ind); ++ii ) {

				// allows derived classes to apply their own logic for starting criteria
				if ( !(this->state_is_starting_aa_type( res1_ind, ii )) ) continue;

				// if symmetric_ need to check one-body energy to make sure rotamer doesn't clash with itself accross the interface
				//    SymmetricRotamerSets compute_energies() adds 2-body e to the 1-body e for cases where symm residues that interact with own clones
				if ( symmetric_ ) {
					Real one_body_i = (ig_->get_one_body_energy_for_node_state( first_node_ind, ii ))/scmult_1b;
					if ( one_body_i > clash_threshold_ ) {
						continue;
					}
					// check for residues that hydrogen bond with their symmetric clones; these are stored in 1-body in symmeric case!!
					if ( one_body_i < onebody_hb_threshold_ ) {
						utility::vector1< HBondResStructCOP > residues;
						residues.push_back(
							pointer::make_shared< HBondResStruct >(
							res1_ind,
							platform::uint( ii ),
							rotamer_sets_->rotamer_set_for_moltenresidue( first_ni )->rotamer( ii )->name1(),
							orig_pose_->pdb_info()->chain( res1_ind ),
							true,
							false,
							false
							)
						);
						store_network( residues, one_body_i, true, true, false, false ); //network complete, store it
					}
				}
				for ( int jj = 1; jj <= ig_->get_num_states_for_node(second_node_ind); ++jj ) {
					if ( !(this->pair_meets_starting_criteria( res1_ind, ii, res2_ind, jj )) ) continue;

					if ( symmetric_ ) {
						Real one_body_j = (ig_->get_one_body_energy_for_node_state( second_node_ind, jj ))/scmult_1b;
						if ( one_body_j > clash_threshold_ ) {
							continue;
						} else if ( one_body_j < onebody_hb_threshold_ ) {
							utility::vector1< HBondResStructCOP > residues;
							residues.push_back(
								pointer::make_shared< HBondResStruct >(
								res2_ind,
								platform::uint( jj ),
								rotamer_sets_->rotamer_set_for_moltenresidue( second_ni )->rotamer( jj )->name1(),
								orig_pose_->pdb_info()->chain( res2_ind ),
								true,
								false,
								false
								)
							);
							store_network( residues, one_body_j, true, true, false, false ); //network complete, store it
						}
					}
					auto const & pdedge =  static_cast< PDEdge const & >(edge);
					Real twobody( (first_node_ind < second_node_ind) ? (pdedge.get_two_body_energy(ii, jj))/scmult_2b : (pdedge.get_two_body_energy(jj, ii))/scmult_2b );
					twobody = twobody * ( this->upweight_starting_twobody_energy() );
					if ( twobody < 0.0 ) {
						twobody = this->scale_twobody_energy( twobody, rotamer_sets_->rotamer_set_for_moltenresidue(first_ni)->rotamer(ii)->name1(), rotamer_sets_->rotamer_set_for_moltenresidue(second_ni)->rotamer(jj)->name1() );
					}
					// if energy < threshold, we found an h-bond, and if doesn't clash with anything we add the rotamer to the current network
					if ( twobody < hb_threshold ) {
						core::Real init_sc = twobody;
						utility::vector1< HBondResStructCOP > residues;
						residues.push_back(
							pointer::make_shared< HBondResStruct >(
							res1_ind,
							platform::uint( ii ),
							rotamer_sets_->rotamer_set_for_moltenresidue( first_ni )->rotamer( ii )->name1(),
							orig_pose_->pdb_info()->chain( res1_ind ),
							true,
							false,
							false
							)
						);
						recursive_traverse( second_node_ind, jj, res2_ind, res1_ind, residues, 1, init_sc, hb_threshold ); //recursive traverse of IG
					}
				}
			}
		}
	}
}

///@details recursive call for depth-first traversal of IG
void
HBNet::recursive_traverse( int const new_node_ind, int const newstate, Size const newres, Size const prevres,
	utility::vector1< HBondResStructCOP > residues, Size network_rec_count, Real init_sc, Real const hb_threshold, bool const second_search /* false */ )
{
	// secondary search option allows dead-end cases to search again with a lowered threshold
	if ( !second_search ) {
		residues.push_back(
			pointer::make_shared< HBondResStruct >(
			newres,
			platform::uint( newstate ),
			rotamer_sets_->rotamer_set_for_moltenresidue( new_node_ind )->rotamer( newstate )->name1(),
			orig_pose_->pdb_info()->chain( newres ),
			true,
			false,
			false
			)
		);
	}
	Size prev_resnum = prevres;
	network_rec_count++; // numer of recursive calls; network_rec_count = # rotamers in the current h-bond network

	if ( network_rec_count >= max_network_size_ ) { //arbitrary limit for size of h-bond network to prevent explosion
		store_network( residues, init_sc, false, false, false, false );
		return;
	}
	Size instance_rec_call_cnt(0); // Keep track of how many times this particular instance makes new recursive calls

	// Continue depth-first traversal of IG starting from new_node_ind
	for ( auto edge_iter = ig_->get_node(new_node_ind)->edge_list_begin(); edge_iter != ig_->get_node(new_node_ind)->edge_list_end(); ++edge_iter ) {
		runtime_assert(newres == rotamer_sets_->moltenres_2_resid(new_node_ind));
		int const second_node_ind = (*edge_iter)->get_other_ind(new_node_ind);
		platform::uint const second_ni = (*edge_iter)->get_other_ind(new_node_ind);
		Size const res2 = rotamer_sets_->moltenres_2_resid(second_ni);

		if ( res2 == prev_resnum ) { //check for case where immediately return to previous residue (e.g. res_x -> res_y -> res_x)
			continue;
		}
		Size scmult_2b((symmetric_) ? symm_info_->score_multiply(newres,res2) : 1);
		Size scmult_1b((symmetric_) ? symm_info_->score_multiply_factor() : 1);

		// now search the IG for new rotamer states that have h-bonds with this rotamer state
		for ( int jj = 1; jj <= ig_->get_num_states_for_node(second_node_ind); ++jj ) {
			Real one_body_j = (ig_->get_one_body_energy_for_node_state( second_node_ind, jj ))/scmult_1b;
			// if symmetric interface design need to check one-body energy to make sure rotamer doesn't clash with itself accross the interface
			if ( symmetric_ && one_body_j > clash_threshold_ ) {
				continue;
			}
			auto * pdedge = static_cast< PDEdge * >(*edge_iter);
			Real twobody( (new_node_ind < second_node_ind) ? (pdedge->get_two_body_energy(newstate, jj))/scmult_2b : (pdedge->get_two_body_energy(jj, newstate))/scmult_2b );
			if ( twobody < 0.0 ) {
				twobody = this->scale_twobody_energy( twobody, rotamer_sets_->rotamer_set_for_moltenresidue(new_node_ind)->rotamer(newstate)->name1(), rotamer_sets_->rotamer_set_for_moltenresidue(second_node_ind)->rotamer(jj)->name1() );
			}
			if ( twobody < hb_threshold ) {
				init_sc +=  twobody;
				bool cycle = false;
				// if we get back to any of the starting residues, stop:
				//  it's more efficient to branch overlapping networks later than detect many duplicates and subnetworks through rec explosion
				//if ( start_res_vec_.find( res2 ) == start_res_vec_.end() ) {
				if ( start_res_vec_.find( res2 ) == start_res_vec_.end() || network_rec_count < 3 ) { //otherwise we could just cycle back and forth between 2 starting reisdues and never go deeper
					// If new res does not clash with residues already in the network...
					if ( !( check_clash( residues, second_ni, jj, res2, init_sc, cycle ) ) ) {
						if ( cycle ) { // If found a cycle, store network
							store_network( residues, init_sc, false, true, false, false );
						} else { // Else, add it to the network and keep going!  Recursive call to rec_trav
							//init_sc += one_body_j;
							recursive_traverse( second_node_ind, jj, res2, newres, residues, network_rec_count, init_sc, hb_threshold ); //recursive call
							instance_rec_call_cnt++;
						}
					}
				} else { // if we are back to original starting rotamer, don't need check_clash else check_clash
					auto const start_state = static_cast<int>(residues.front()->rot_index);
					if ( residues.front()->resnum == res2 && start_state == jj && residues.size() > 2 ) {
						store_network( residues, init_sc, true, true, false, false );
						//return; 16/03
					} else if ( !( check_clash( residues, second_ni, jj, res2, init_sc, cycle ) ) ) {
						// if we get back to any of the starting residues, stop:
						//  it's more efficient to branch overlapping networks later than detect many duplicates and subnetworks through rec explosion
						store_network( residues, init_sc, false, cycle, false, false );
					}
				}
			}
		}
	}
	// If this instance of rec_trav doesn't make any recursive calls or store a network, then store the current network.
	//   In other words, this node of this h-bond network isn't able to form any other possible h-bonds, so stop and store where we're at
	//if ( store_subnetworks_ || ( instance_rec_call_cnt == 0 && residues.size() > 2 ) ){
	if ( store_subnetworks_ || ( instance_rec_call_cnt == 0 && residues.size() > 2 ) ||
			rotamer_sets_->rotamer_set_for_moltenresidue(new_node_ind)->rotamer(newstate)->name1() == 'S' ||
			rotamer_sets_->rotamer_set_for_moltenresidue(new_node_ind)->rotamer(newstate)->name1() == 'T' ||
			rotamer_sets_->rotamer_set_for_moltenresidue(new_node_ind)->rotamer(newstate)->name1() == 'Y' ) {
		store_network( residues, init_sc, false, false, false, false ); // network complete, store it
	}
	if ( secondary_search_ && instance_rec_call_cnt == 0 && !second_search ) {
		recursive_traverse( new_node_ind, newstate, newres, prevres, residues, network_rec_count, init_sc, secondary_threshold_, true ); //call this same recursive instance again with lower threshold
		instance_rec_call_cnt++;
	}
}//rec_trav

void
HBNet::MC_traverse_IG( ){

	if ( ! hb_database_ ) {
		hb_database_ = core::scoring::hbonds::HBondDatabase::get_database();
	}

	if ( symmetric_ ) {
		static_cast< core::conformation::symmetry::SymmetricConformation & >( orig_pose_->conformation() ).recalculate_transforms();
		static_cast< core::pack::hbonds::MCHBNetInteractionGraph & >( * ig_ ).find_symmetric_hbonds( *symm_info_, *orig_pose_, onebody_hb_threshold_ );
	}

	//Remove edges with no incident edges. If they can not be part of a network then they are wasting our time.
	for ( utility::graph::EdgeListIterator it = hbond_graph_->edge_list_begin(), end = hbond_graph_->edge_list_end();
			it != end; ) {

		auto * edge = static_cast< AtomLevelHBondEdge * >( *it );
		if ( hbond_graph_->get_node( edge->get_first_node_ind() )->num_edges() == 1 &&
				hbond_graph_->get_node( edge->get_second_node_ind() )->num_edges() == 1
				) {
			//delete edges that have no neighboring (incident) edges and thus can not create a network
			++it;//very important to increment before deleting the edge that "it" is pointing to
			hbond_graph_->delete_edge( edge );
		} else if ( monte_carlo_seed_is_dead_end( edge ) ) {
			//delete edges that have neighbors but there does not exist an incident edge that this edge can coexist with
			++it;
			hbond_graph_->delete_edge( edge );
		} else {
			//fill in hbond info
			core::pack::hbonds::determine_atom_level_edge_info(
				* edge,
				* rotamer_sets_,
				* hb_database_,
				orig_pose_->energies().tenA_neighbor_graph(),
				* orig_pose_,
				secondary_threshold_,
				symm_info_
			);
			++it;
		}
	}

	core::Size const scmult_1b( symmetric_ ? symm_info_->score_multiply_factor() : 1 );

	//Do not consider hbonds where one side clashes with the background
	//Also, find satisfying interactions with background
	for ( core::Size node_id = 1; node_id <= hbond_graph_->num_nodes(); ++node_id ) {
		AtomLevelHBondNode * node = hbond_graph_->get_hbondnode( node_id );
		core::Real const one_body_1 = ig_->get_one_body_energy_for_node_state( node->moltenres(), node->local_rotamer_id() ) / scmult_1b;

		if ( one_body_1 > clash_threshold_ ) {
			hbond_graph_->drop_all_edges_for_node( node_id );
		} else if ( node->num_edges() ) { //do not waste time on nodes with no edges
			if ( node->num_edges() > 0 ) {
				core::pack::hbonds::determine_atom_level_node_info( * node, * rotamer_sets_, core_residues_ );
				core::pack::hbonds::find_satisfying_interactions_with_background( * node, * rotamer_sets_, * packer_neighbor_graph_, * ala_pose_, secondary_threshold_ );
			}
		}
	}


	////////////////
	//Look for seeds
	for ( utility::graph::EdgeListConstIterator it = hbond_graph_->const_edge_list_begin();
			it != hbond_graph_->const_edge_list_end(); ++it ) {
		auto const * edge = static_cast< AtomLevelHBondEdge const * >( *it );

		AtomLevelHBondNode const * node1 = hbond_graph_->get_hbondnode( edge->get_first_node_ind() );
		AtomLevelHBondNode const * node2 = hbond_graph_->get_hbondnode( edge->get_second_node_ind() );

		core::Size const resid1 = rotamer_sets_->moltenres_2_resid( node1->moltenres() );
		core::Size const resid2 = rotamer_sets_->moltenres_2_resid( node2->moltenres() );

		core::Size const rot1 = node1->local_rotamer_id();
		core::Size const rot2 = node2->local_rotamer_id();

		if ( edge_can_yield_monte_carlo_seed( resid1, resid2 )
				&& edge->energy() <= monte_carlo_seed_threshold_
				&& pair_meets_starting_criteria( resid1, rot1, resid2, rot2 )
				&& ! monte_carlo_seed_is_dead_end( edge ) ) {
			monte_carlo_seeds_.push_back( edge );
		}
	}

}


///@details traverse EnergyGraph of static pose to find all native networks
void
HBNet::traverse_native( Pose const & pose, Real const hb_threshold )
{
	for ( core::Size res : start_res_vec_ ) {
		utility::graph::Node const * node = pose.energies().energy_graph().get_node( res );

		for ( utility::graph::EdgeListConstIterator egraph_it = node->const_edge_list_begin();
				egraph_it != node->const_edge_list_end(); ++egraph_it ) {
			auto const * eedge = static_cast< core::scoring::EnergyEdge const * > (*egraph_it);
			Size other_res = eedge->get_other_node(node->get_node_index())->get_node_index();

			Size res1_ind(res),res2_ind(other_res);
			if ( symmetric_ ) {
				res1_ind = get_ind_res( *orig_pose_, res);
				res2_ind = get_ind_res( *orig_pose_, other_res);
			}

			if ( eedge->operator[](core::scoring::hbond_sc) < hb_threshold ) {
				// Initialize lists
				utility::vector1< HBondResStructCOP > residues;
				residues.push_back(
					pointer::make_shared< HBondResStruct >(
					res1_ind,
					0,
					pose.residue( res1_ind ).name1(),
					orig_pose_->pdb_info()->chain( res1_ind ),
					pose.residue( res1_ind ).is_protein(),
					false,
					pose.residue( res1_ind ).is_ligand()
					)
				);
				rec_trav_native( pose, res2_ind, res1_ind, residues, hydrogen_bond_threshold_ ); // recursive call
			}
		}
	}
}

void
HBNet::rec_trav_native( Pose const & pose, Size new_res, Size prev_res, utility::vector1< HBondResStructCOP > residues, Real const hb_threshold )
{
	bool found_new_res_in_net = false;
	for ( utility::vector1< HBondResStructCOP >::const_iterator i = residues.begin(); i != residues.end(); ++i ) {
		if ( (*i)->resnum == new_res ) {
			store_network(residues, 0.0, false, false, true );
			return;
		}
	}
	//residues.push_back( HBondResStructCOP( new HBondResStruct( new_res, 0, pose.residue(new_res).name1(), pose.pdb_info()->chain(new_res), pose.residue(new_res).is_protein(), pose.residue(new_res).is_water(), pose.residue(new_res).is_ligand() ) ) );

	residues.push_back(
		pointer::make_shared< HBondResStruct >(
		new_res,
		platform::uint( 0 ),
		pose.residue( new_res ).name1(),
		orig_pose_->pdb_info()->chain( new_res ),
		pose.residue( new_res ).is_protein(),
		pose.residue( new_res ).name1() == 'w',
		pose.residue( new_res ).is_ligand()
		)
	);

	utility::graph::Node const * node = pose.energies().energy_graph().get_node( new_res );
	for ( utility::graph::EdgeListConstIterator egraph_it = node->const_edge_list_begin();
			egraph_it != node->const_edge_list_end(); ++egraph_it ) {
		auto const * eedge = static_cast< core::scoring::EnergyEdge const * > (*egraph_it);
		Size other_res = eedge->get_other_node(node->get_node_index())->get_node_index();

		if ( symmetric_ ) {
			other_res = get_ind_res( *orig_pose_, other_res);
		}

		if ( eedge->operator[](core::scoring::hbond_sc) < hb_threshold ) {
			if ( other_res == prev_res ) {
				continue;
			} else if ( residues.front()->resnum == other_res ) {
				if ( store_subnetworks_ ) {
					store_network( residues, 0.0, true, false, true );
					found_new_res_in_net = true;
				} else {
					rec_trav_native( pose, other_res, new_res, residues, hydrogen_bond_threshold_ ); // recursive call
				}
			} else {
				rec_trav_native( pose, other_res, new_res, residues, hydrogen_bond_threshold_ );
				found_new_res_in_net = true; // recursive call
			}
		}
	}
	if ( !found_new_res_in_net ) {
		store_network(residues, 0.0, false, false, false, true );
	}
}

///@details Check if a new rotamer state clashes with any rotamer states already in a given h-bond network
///    return of true = it clashes, false = no clashes
bool
HBNet::check_clash(utility::vector1< HBondResStructCOP > const & residues, platform::uint new_node_ind, Size newstate, Size newres, Real & init_score, bool & cycle)
{
	cycle = false;
	bool clash = false;
	core::PackerEnergy twobody(0.0);

	for ( const auto & residue : residues ) {
		platform::uint otherstate( residue->rot_index );
		platform::uint other_node_ind(rotamer_sets_->resid_2_moltenres(residue->resnum));
		auto const new_ind = static_cast<int>(new_node_ind);
		auto const other_ind = static_cast<int>(rotamer_sets_->resid_2_moltenres(residue->resnum));
		auto const new_st = static_cast<int>(newstate);
		auto const other_st = static_cast<int>(otherstate);

		if ( newres == residue->resnum ) {
			if ( newstate == otherstate ) { // If rotamer state already in network, we found a cycle, no need to clash check
				cycle = true;
				clash = false;
				return clash;
			} else {
				// Else, same residue, different rotamer states, so return clash = true;
				clash = true;
				return clash;
			}
		} else if ( ig_->find_edge(new_ind, other_ind) != nullptr ) { // this function is order-independent; returns the same Edge regardless of order
			core::pack::interaction_graph::EdgeBase * temp_edge = ig_->find_edge(new_ind, other_ind);
			auto * temp_pdedge = static_cast< PDEdge * >(temp_edge);
			int first_node_ind = temp_pdedge->get_first_node_ind(); // need to get the first node of the Edge (lower-index node) (we're getting this from the IG)
			platform::uint first_ni = first_node_ind;
			Size res1 = rotamer_sets_->moltenres_2_resid(first_ni); // get residue number that corresopnds to the first node of Edge

			// temp_pdedge->get_two_body_energy() is order-DEPENDENT; lower index (first index) must always be first)
			//    so we need to determine does first_ni correspond to new_node_ind?
			if ( first_ni == new_node_ind && res1 == newres ) { // second half of if is sanity check; shouldn't be necessary but best to be safe
				twobody = temp_pdedge->get_two_body_energy( new_st, other_st );
			} else if ( first_ni == other_node_ind && res1 == residue->resnum ) {
				// or other_node_ind?
				twobody = temp_pdedge->get_two_body_energy( other_st, new_st );
			} else {
				runtime_assert( ( first_ni != new_node_ind && res1 == newres ) || ( first_ni != other_node_ind && res1 == residue->resnum ) );
			}

			// now that we have the correct two-body energy from the IG, check if there is a clash
			if ( symmetric_ ) {
				twobody = twobody/(symm_info_->score_multiply_factor());
			}
			if ( twobody >= clash_threshold_ ) {
				clash = true;
				return clash;
			} else if ( twobody >= hydrogen_bond_threshold_ ) { //include fa_rep values for residues that may contact; h-bonding residues ( < hydrogen_bond_threshold) already included
				init_score += twobody;
			}
		}
	}
	return clash;
}//check_clash

///@details Function to check whether two h-bond networks clash with eachother
///    return of true = they clash, false = no clashes (networks i and j are compatible)
bool
HBNet::net_clash( HBondNetStruct const & i, HBondNetStruct const & j )
{
	return net_clash( i.residues, j.residues );
}

///@details Function to check whether two h-bond networks clash with eachother
///    return of true = they clash, false = no clashes (networks i and j are compatible)
bool
HBNet::net_clash(utility::vector1< HBondResStructCOP > const & residues_i, utility::vector1< HBondResStructCOP > const & residues_j)
{
	if ( monte_carlo_branch_ ) return monte_carlo_net_clash( residues_i, residues_j ) ;

	core::PackerEnergy twobody(0.0);

	for ( const auto & res_i : residues_i ) {
		int const state_i = res_i->rot_index;
		int const node_i = rotamer_sets_->resid_2_moltenres(res_i->resnum);
		for ( const auto & res_j : residues_j ) {
			int const state_j = res_j->rot_index;
			int const node_j = rotamer_sets_->resid_2_moltenres(res_j->resnum);
			// If networks share a residue, make sure it is the same rotamer state; otherwise, return clash = true;
			//    except at *(start_res_vec_.begin()) in ligand()binding design case: check for compatible ligand states in merge_networks()
			//NEED TO CHECK LIGAND COMPATIBILITY HERE //NEED TO FIX
			if ( (this->ligand()) && ( res_i->resnum == (this->ligand()) || res_j->resnum == (this->ligand()) ) ) {
				continue;
			} else if ( res_i->resnum == res_j->resnum && ( state_i != state_j ) ) {
				if ( rotamer_sets_ != nullptr ) {
					if ( res_i->aa != res_j->aa ) {
						return true;
					} else {
						core::conformation::ResidueCOP r_i = rotamer_sets_->rotamer_set_for_residue(res_i->resnum)->rotamer(state_i);
						core::conformation::ResidueCOP r_j = rotamer_sets_->rotamer_set_for_residue(res_j->resnum)->rotamer(state_j);
						//if ( core::scoring::automorphic_rmsd(*r_i, *r_j, true) > SC_RMSD_CUTOFF ) {
						if ( core::scoring::residue_sc_rmsd_no_super( r_i, r_j, true ) > SC_RMSD_CUTOFF ) {
							return true;
						}
					}
				} else {
					return true;
				}
			}

			if ( ig_ != nullptr  ) { //if IG still exists, use 2-body for faster lookup
				if ( ig_->find_edge(node_i, node_j) != nullptr ) { // If edge exist s, look up 2-body energy for rotamer state_i and state_j
					core::pack::interaction_graph::EdgeBase * temp_edge = ig_->find_edge(node_i, node_j);
					auto * temp_pdedge = static_cast< PDEdge * >(temp_edge);
					int first_node_ind = temp_pdedge->get_first_node_ind();

					if ( first_node_ind == node_i ) {
						twobody = temp_pdedge->get_two_body_energy( state_i, state_j );
					} else if ( first_node_ind == node_j ) {
						twobody = temp_pdedge->get_two_body_energy( state_j, state_i );
					} else {
						runtime_assert( first_node_ind == node_i || first_node_ind == node_j );
					}
					if ( symmetric_ ) {
						twobody = twobody/(symm_info_->score_multiply(res_i->resnum,res_j->resnum));
					}
					if ( twobody >= clash_threshold_ ) {
						return true;
					}
				}
			} else { // safeguard: if IG has been cleared or reset to null, we can still check clash
				Pose ala_copy = *ala_pose_;
				core::conformation::ResidueCOP rot_i = rotamer_sets_->rotamer_set_for_moltenresidue(node_i)->rotamer(res_i->rot_index);
				core::conformation::ResidueCOP rot_j = rotamer_sets_->rotamer_set_for_moltenresidue(node_j)->rotamer(res_j->rot_index);

				ala_copy.replace_residue(res_i->resnum, *rot_i, false);
				ala_copy.replace_residue(res_j->resnum, *rot_j, false);
				(*scorefxn_)(ala_copy); //sfxn->eval methods ASSUME THE POSE HAS BEEN SCORED?
				core::scoring::EnergyMap new_emap;
				init_scorefxn_->eval_cd_2b_sc_sc(*rot_i, *rot_j, ala_copy, new_emap);
				Real farep_val = new_emap.get(core::scoring::fa_rep);
				if ( farep_val > clash_threshold_ ) {
					return true;
				}
			}
		}
	}
	return false;
}//net_clash

bool one_network_is_subset_of_other( utility::vector1< HBondResStructCOP > const & residues_i, utility::vector1< HBondResStructCOP > const & residues_j ){

	utility::vector1< HBondResStructCOP > const & larger =  ( residues_i.size() > residues_j.size() ? residues_i : residues_j );
	utility::vector1< HBondResStructCOP > const & smaller = ( residues_i.size() > residues_j.size() ? residues_j : residues_i );

	for ( HBondResStructCOP const & small_guy : smaller ) {
		bool exists_in_larger = false;
		for ( HBondResStructCOP const & big_guy : larger ) {
			if ( *small_guy == *big_guy ) {
				exists_in_larger = true;
				break;
			}
		}
		if ( ! exists_in_larger ) return false;
	}

	return true;
}

bool
HBNet::monte_carlo_net_clash( utility::vector1< HBondResStructCOP > const & residues_i, utility::vector1< HBondResStructCOP > const & residues_j ) const {

	if ( residues_i.size() == 0 || residues_j.size() == 0 ) return false;
	if ( one_network_is_subset_of_other( residues_i, residues_j) ) return true;

	utility::vector1< core::Size > global_rots1;
	global_rots1.reserve( residues_i.size() );
	//core::Size const offset1 = rotamer_sets_->nrotamer_offset_for_moltenres( rotamer_sets_->resid_2_moltenres( residues_i.front()->resnum ) );
	for ( HBondResStructCOP const & res : residues_i ) {
		core::Size const offset1 = rotamer_sets_->nrotamer_offset_for_moltenres( rotamer_sets_->resid_2_moltenres( res->resnum ) );
		global_rots1.push_back( offset1 + res->rot_index );
	}

	utility::vector1< core::Size > global_rots2;
	global_rots2.reserve( residues_j.size() );
	//core::Size const offset2 = rotamer_sets_->nrotamer_offset_for_moltenres( rotamer_sets_->resid_2_moltenres( residues_j.front()->resnum ) );
	for ( HBondResStructCOP const & res : residues_j ) {
		core::Size const offset2 = rotamer_sets_->nrotamer_offset_for_moltenres( rotamer_sets_->resid_2_moltenres( res->resnum ) );
		global_rots2.push_back( offset2 + res->rot_index );
	}

	for ( core::Size const rot1 : global_rots1 ) {
		AtomLevelHBondNode const * node1 = hbond_graph_->get_hbondnode( rot1 );
		for ( core::Size const rot2 : global_rots2 ) {
			if ( node1->clashes( rot2 ) ) return true;

			// check for networks that have common intersections that rotamers are compatible
			if ( rotamer_sets_->moltenres_for_rotamer(rot1) == rotamer_sets_->moltenres_for_rotamer(rot2) && ( rot1 != rot2 ) ) {
				if ( rotamer_sets_->rotamer(rot1)->name1() != rotamer_sets_->rotamer(rot2)->name1() ) {
					return true;
				} else {
					core::conformation::ResidueCOP r_1 = rotamer_sets_->rotamer(rot1);
					core::conformation::ResidueCOP r_2 = rotamer_sets_->rotamer(rot2);
					//if ( core::scoring::automorphic_rmsd(*r_i, *r_j, true) > SC_RMSD_CUTOFF ) {
					if ( core::scoring::residue_sc_rmsd_no_super( r_1, r_2, true ) > SC_RMSD_CUTOFF ) {
						return true;
					}
				}
			}
		}
	}

	return false;
}

bool
HBNet::network_already_stored( utility::vector1< HBondResStructCOP > & residues, utility::vector1< HBondResStructCOP > & i_residues )
{
	bool already_stored(false);
	if ( residues.size() == i_residues.size() ) {
		std::sort( residues.begin(), residues.end(), compare_hbond_residues() );
		std::sort( i_residues.begin(), i_residues.end(), compare_hbond_residues() );

		already_stored = true;
		std::vector< HBondResStructCOP >::const_iterator j = i_residues.begin();
		for ( std::vector< HBondResStructCOP >::const_iterator k = residues.begin(); k != residues.end(); ++k ) {
			if ( ( (*j)->resnum == (*k)->resnum ) && ( (*j)->aa == (*k)->aa ) ) {
				if ( rotamer_sets_ == nullptr ) continue;
				core::conformation::ResidueCOP rot1(rotamer_sets_->rotamer_set_for_moltenresidue(rotamer_sets_->resid_2_moltenres((*j)->resnum))->rotamer((*j)->rot_index));
				core::conformation::ResidueCOP rot2(rotamer_sets_->rotamer_set_for_moltenresidue(rotamer_sets_->resid_2_moltenres((*k)->resnum))->rotamer((*k)->rot_index));
				//if ( core::scoring::automorphic_rmsd(*rot1, *rot2, true /*superpose*/ ) > SC_RMSD_CUTOFF ) {
				if ( core::scoring::residue_sc_rmsd_no_super( rot1, rot2, true /*final group only*/ ) > SC_RMSD_CUTOFF ) {
					already_stored = false;
					break;
				}
			} else {
				already_stored = false;
				break;
			}
			++j;
		}
	}
	return already_stored;
}

void
HBNet::store_network(utility::vector1< HBondResStructCOP > residues, Real init_score, bool term_w_start, bool term_w_cycle, bool score_now, bool native )
{
	bool already_stored(false);
	for ( auto & i : network_vector_ ) {
		already_stored = network_already_stored( residues, i->residues );

		if ( already_stored ) {
			if ( init_score < i->score && !score_now ) {
				//replace i with new network with better score
				i->is_native = native;
				i->term_w_start = term_w_start;
				i->term_w_cycle = term_w_cycle;
				i->residues = residues;
				if ( (this->ligand()) ) {
					i->lig_state_list.push_back((*(find_HBondResStruct(i->residues, this->ligand())))->rot_index);
				}
				i->score = init_score; //don't need to adjust for symmetric_ because did that when tracking IG energies to generate init_score
			} else if ( this->ligand() ) {
				i->lig_state_list.push_back((*(find_HBondResStruct(i->residues, this->ligand())))->rot_index);
			}
			return;
		}
	}
	if ( !( already_stored ) ) {
		HBondNetStructOP new_net = pointer::make_shared< HBondNetStruct >();
		new_net->is_native = native;
		new_net->term_w_start = term_w_start;
		new_net->term_w_cycle = term_w_cycle;
		new_net->residues = residues;
		new_net->lig_state_list.clear();
		if ( (this->ligand()) ) {
			new_net->lig_state_list.push_back((*(find_HBondResStruct((new_net)->residues, this->ligand())))->rot_index);
		}
		if ( score_now ) {
			Pose copy = *ala_pose_;
		} else {
			new_net->score = init_score; //don't need to adjust for symmetric_ because did that when tracking IG energies to generate init_score
			network_vector_.push_back(new_net);
		}
	}
}//store_network

void
HBNet::minimize_network( Pose & pose, HBondNetStruct & network, bool residues_already_placed /* true */ )
{
	if ( !residues_already_placed ) {
		place_rots_on_pose( pose, network, network.is_native );
		pose.update_residue_neighbors();
	}
	core::kinematics::MoveMapOP mm = pointer::make_shared< core::kinematics::MoveMap >();
	mm->set_bb( false ); mm->set_chi( false ); mm->set_jump( false );
	for ( utility::vector1< HBondResStructCOP >::const_iterator rit=network.residues.begin(); rit != network.residues.end(); rit++ ) {
		if ( (this->ligand()) && (*rit)->resnum == (this->ligand()) ) {
			continue;
		}
		mm->set_chi( (*rit)->resnum, true );
	}
	// TODO NEED TO ADD CODE HERE FOR BRIDGING_WATER CASE; NEED TO CHECK RESIDUE NUMBERS OF NETWORK RESIDUES AFTER WATER
	minimization_packing::MinMoverOP min_mover;
	if ( symmetric_ ) {
		core::pose::symmetry::make_symmetric_movemap( pose, *mm );
	}
	min_mover = pointer::make_shared< minimization_packing::MinMover >( mm, scorefxn_, "dfpmin_armijo_nonmonotone", 0.0001, true );
	min_mover->apply(pose);
	pose.update_residue_neighbors();

	//update roatmers and store, so we don't have to do this again
	Size old_total( (symmetric_) ? symm_info_->num_independent_residues() : orig_pose_->total_residue() );
	network.rotamers.clear();
	for ( Size r = 1; r <= old_total; r++ ) {
		if ( find_HBondResStruct( network.residues, r ) != network.residues.end() ) {
			network.rotamers.push_back( pose.residue(r).clone() );
		}
	}
}


void
HBNet::score_network_on_pose( Pose & pose, HBondNetStruct & i )
{
	if ( i.scored ) return;

	//Pose ala_copy = *ala_pose_; // TODO does scoring get copied too? looks like it does, but NEED TO CHECK THIS
	//ala_copy.update_residue_neighbors();
	//(*scorefxn_)(ala_copy); //already scored

	//set baseline energies from background pose
	core::scoring::EnergyMap & baseline_emap = ala_pose_->energies().total_energies();
	core::Real total_baseline = baseline_emap.get(core::scoring::total_score);
	core::Real hbond_bb_sc_baseline = baseline_emap.get(core::scoring::hbond_bb_sc);

	(*scorefxn_)(pose);
	core::scoring::EnergyMap & new_emap = pose.energies().total_energies();
	Real pose_total = new_emap.get(core::scoring::total_score) - total_baseline;
	Real e( pose_total );
	if ( (new_emap.get( core::scoring::hbond_bb_sc) - hbond_bb_sc_baseline) <= MIN_HB_E_CUTOFF ) {
		i.term_w_bb = true;
	} else {
		i.term_w_bb = false;
	}
	//core::Real cycle_bonus = -2.0;
	//if ( i.term_w_cycle ) e = e + cycle_bonus;
	//if ( (this->ligand()) && i.term_w_start ) e = e + ret2lig_bonus_;
	if ( symmetric_ ) e = e/(symm_info_->num_bb_clones() + 1.0);
	//normalize
	e = e/((Real)i.residues.size());
	i.score = e;
	i.scored = true;
}

void
HBNet::merge_2_branched_networks(utility::vector1< HBondResStructCOP > const & residues1, utility::vector1< HBondResStructCOP > const & residues2, utility::vector1< HBondResStructCOP > & new_residues)
{
	for ( const auto & res1 : residues1 ) {
		bool found(false);
		for ( utility::vector1< HBondResStructCOP >::const_iterator n = new_residues.begin(); n != new_residues.end(); ++n ) {
			if ( res1->resnum == (*n)->resnum ) {
				found = true;
				break;
			}
		}
		if ( !found ) {
			new_residues.push_back(res1);
		}
	}
	for ( const auto & res2 : residues2 ) {
		bool found(false);
		for ( utility::vector1< HBondResStructCOP >::const_iterator n = new_residues.begin(); n != new_residues.end(); ++n ) {
			if ( res2->resnum == (*n)->resnum ) {
				found = true;
				break;
			}
		}
		if ( !found ) {
			new_residues.push_back(res2);
		}
	}
}

void
HBNet::merge_2_branched_networks(HBondNetStruct const & i, HBondNetStruct const & j, HBondNetStructOP new_network)
{
	if ( (this->ligand()) && (i.lig_state_list.empty() || j.lig_state_list.empty()) ) {
		return;
	}
	utility::vector1< HBondResStructCOP > new_residues(0);
	merge_2_branched_networks(i.residues, j.residues, new_residues);

	new_network->sort_first_by_tot_unsat = true;
	new_network->scored = false;
	new_network->residues = new_residues;
	new_network->is_native = ( i.is_native && j.is_native ) ? true : false;
	new_network->is_extended = ( ( i.is_native && !j.is_native ) || ( !i.is_native && j.is_native ) || ( i.is_extended || j.is_extended ) ) ? true : false;
	new_network->term_w_bb = ( i.term_w_bb || j.term_w_bb ) ? true : false;
	new_network->term_w_cycle = ( i.term_w_cycle || j.term_w_cycle ) ? true : false;
	new_network->term_w_start = ( i.term_w_start || j.term_w_start ) ? true : false;

	if ( this->ligand() ) {
		std::vector<platform::uint> new_lig_state_list(1000);
		std::vector<platform::uint> liglist1 = i.lig_state_list;
		std::vector<platform::uint> liglist2 = j.lig_state_list;
		std::sort( liglist1.begin(), liglist1.end() );
		std::sort( liglist2.begin(), liglist2.end() );
		std::vector<platform::uint>::iterator lit;
		lit=std::set_intersection(liglist1.begin(), liglist1.end(), liglist2.begin(), liglist2.end(), new_lig_state_list.begin());
		new_lig_state_list.resize(lit-new_lig_state_list.begin());
		if ( new_lig_state_list.empty() ) {
			new_lig_state_list.push_back((*(find_HBondResStruct(i.residues, this->ligand())))->rot_index);
			new_lig_state_list.push_back((*(find_HBondResStruct(j.residues, this->ligand())))->rot_index);
			new_network->lig_state_list = new_lig_state_list;
		} else {
			new_network->lig_state_list = new_lig_state_list;
		}
	}
	//score and num_unsat, etc. should be set by place_rotamers_and_score()
	new_network->lig_num_unsatisfied = 0;
	new_network->num_unsat_Hpol = 0;
	new_network->num_heavy_unsat = 0;
	new_network->score = 0.0;
}

//consider Ser and Thr idential for benchmarking purposes
bool
HBNet::networks_identical_aa_sequence( HBondNetStruct const & i, HBondNetStruct const & j )
{
	if ( i.residues.size() == j.residues.size() ) {
		std::vector< char > i_aa(0);
		for ( const auto & residue : i.residues ) {
			if ( residue->aa == 'T' || residue->aa == 'S' ) {
				i_aa.push_back( 'S' );
			} else {
				i_aa.push_back( residue->aa );
			}
		}
		std::vector< char > j_aa(0);
		for ( const auto & residue : j.residues ) {
			if ( residue->aa == 'T' || residue->aa == 'S' ) {
				j_aa.push_back( 'S' );
			} else {
				j_aa.push_back( residue->aa );
			}
		}
		std::sort( i_aa.begin(), i_aa.end() );
		std::sort( j_aa.begin(), j_aa.end() );
		std::vector<int> a(20);
		std::vector<int>::iterator ait;
		ait=std::set_symmetric_difference(i_aa.begin(), i_aa.end(), j_aa.begin(), j_aa.end(), a.begin());
		a.resize(ait-a.begin());

		if ( a.size() == 0 ) {
			return true;
		}
	}
	return false;
}

bool
HBNet::is_sub_residues( utility::vector1< HBondResStructCOP > & residues1, utility::vector1< HBondResStructCOP > & residues2 )
{
	bool branch(false);
	return is_sub_residues( residues1, residues2, branch );
}

///@details identical networks are considered subsets of each other here by default unless specified by true_if_identical=false
///  std::includes() returns false if 2 vectors are identical
bool
HBNet::is_sub_residues( utility::vector1< HBondResStructCOP > & residues1, utility::vector1< HBondResStructCOP > & residues2, bool & branch, bool true_if_identical /* true */ )
{
	if ( residues1.empty() || residues2.empty() ) {
		return false;
	}
	// If networks is not a subset and has > 1 symmetric difference
	std::sort( residues1.begin(), residues1.end(), compare_hbond_residues() );
	std::sort( residues2.begin(), residues2.end(), compare_hbond_residues() );
	std::vector<HBondResStructCOP> intersec;
	std::vector<HBondResStructCOP> symm_diff;
	std::set_intersection(residues1.begin(), residues1.end(), residues2.begin(), residues2.end(), std::back_inserter(intersec), compare_hbond_residues() );
	std::set_symmetric_difference(residues1.begin(), residues1.end(), residues2.begin(), residues2.end(), std::back_inserter(symm_diff), compare_hbond_residues() );

	// return true if networks identical
	if ( true_if_identical && ( residues1.size() == residues2.size() && symm_diff.size() == 0 ) ) {
		return true;
	}

	//to be a subnetwork, have to overlap at at least one position (2 for ligand-binding design because will overlap at ligand by default)
	//if ( ( (this->ligand()) && v.size() > 1 ) || ( !(this->ligand()) && v.size() > 0 )){
	if ( intersec.size() > 0 ) {
		if ( residues1.size() > residues2.size() && ( std::includes( residues1.begin(), residues1.end(), residues2.begin(), residues2.end(), compare_hbond_residues() ) ) ) {
			return true;
		} else if ( residues2.size() > residues1.size() && ( std::includes( residues2.begin(), residues2.end(), residues1.begin(), residues1.end(), compare_hbond_residues() ) ) ) {
			return true;
		}
		if ( symm_diff.size() >= 2 ) { //only branches if at least 1 diff is unique to each list, !subreslist (e.g. if list1 = {1,2} and list2 = {1,2,3,4}, v2 = 2 but list2 includes list1)
			branch = true;
		}
	}
	return false;
}

/////////////////////////////////////////////////////////////////////////
// NOTE: branch_overlapping_networks() is the bane of classic HBNet's existence
//  In some cases, it's helpful to be exhaustive; but if HBNet ever takes FOREVER to run, this funciton in the culprit
//     It's a catch-22 because this will only happen if many many h-bond networks are found;
//       but that doesn't do you any good if nothing is ever returned due to prohibitevly long runtimes
//  Jack Maguire's Mone Carlo (MC) approach bypasses the need for this function;
//     MC-HBNet builds networks stochastically to completion, meaning that:
//       1) networks can be returned as soon as they are found, so even short runtimes yield results (though not exhaustive)
//       2) the order of traversal/addition no longer matters, eliminating need for the combining/branching function
/////////////////////////////////////////////////////////////////////////

///@details Critical function to HBNet; because of depth-first traversal, need to combine networks for branched amino acids, e.g. ASN that can make multiple h-bonds
///   The reason we need this is because bread-first additions ("branches") are order-dependent, meaning we need to try ALL combinations
///   Branching is also a poor name as it is not "branching" in a computational or graph-theory sense
void
HBNet::branch_overlapping_networks()
{
	if ( network_vector_.empty() ) {
		if ( verbose_ && TR.visible() ) {
			TR << " NO NETWORKS TO BRANCH; RETURNING WITHOUT DOING ANYTHING." << std::endl;
		}
		return;
	}
	merged_vecs_.clear(); //unnecessary if function only used internally
	for ( auto & i : network_vector_ ) {
		i->net_indices.clear();
	}
	if ( TR.visible() ) TR << "recording net_indices " << std::endl;
	for ( auto i = network_vector_.begin(); i != network_vector_.end(); ++i ) {
		// Compare network i with all j, and store best-scoring j that doesn't clash with i
		for ( auto j = i; ++j != network_vector_.end(); ) {
			// If networks overlap by more than 2 residue (will overlap at ligand), continue.
			bool branch(false);
			if ( is_sub_residues( (*i)->residues, (*j)->residues, branch ) || !branch ) {
				continue;
			}

			//TODO need check for small molecule ligands here
			if ( only_native_ || !(net_clash( **i, **j )) ) {
				(*i)->net_indices.push_back( j - network_vector_.begin() );
				(*j)->net_indices.push_back( i - network_vector_.begin() );
			}
		}
	}
	if ( TR.visible() ) TR << "starting recursion " << std::endl;
	for ( auto i = network_vector_.begin(); i != network_vector_.end(); ++i ) {
		Size i_pos = i - network_vector_.begin();
		std::vector< Size > add_index_vec;
		add_index_vec.push_back( i_pos );

		for ( auto j = (*i)->net_indices.begin(); j != (*i)->net_indices.end(); ++j ) {
			if ( *j <= ( i_pos ) ) {
				continue;
			}
			add_index_vec.push_back( *j );

			if ( store_subnetworks_ ) {
				bool found( false );
				std::sort( add_index_vec.begin(), add_index_vec.end() );
				for ( auto & merged_vec : merged_vecs_ ) {
					//std::sort(merged_vec.begin(), merged_vec.end()); // SHOULD ALWAYS ALEADY BE SORTED!
					std::vector< int > v2(100);
					auto it2 = std::set_symmetric_difference(merged_vec.begin(), merged_vec.end(), add_index_vec.begin(), add_index_vec.end(), v2.begin());
					v2.resize(it2-v2.begin());

					if ( v2.size() == 0 ) {
						found = true;
						break;
					}
				}
				if ( !found ) {
					merged_vecs_.push_back(add_index_vec);
				}
			}
			std::sort( network_vector_[*j]->net_indices.begin(), network_vector_[*j]->net_indices.end());
			std::sort( (*i)->net_indices.begin(), (*i)->net_indices.end() );
			std::vector< Size > v(1000);
			std::vector< Size >::iterator it;
			it=std::set_intersection((*i)->net_indices.begin(), (*i)->net_indices.end(), (network_vector_[*j])->net_indices.begin(), (network_vector_[*j])->net_indices.end(), v.begin());
			v.resize(it-v.begin());

			if ( v.size() == 0 ) {
				if ( !store_subnetworks_ ) {
					bool found( false );
					std::sort( add_index_vec.begin(), add_index_vec.end() );
					for ( auto & merged_vec : merged_vecs_ ) {
						//std::sort(merged_vec.begin(), merged_vec.end()); // SHOULD ALWAYS ALEADY BE SORTED!
						std::vector< int > v2(100);
						auto it2 = std::set_symmetric_difference(merged_vec.begin(), merged_vec.end(), add_index_vec.begin(), add_index_vec.end(), v2.begin());
						v2.resize(it2-v2.begin());

						if ( v2.size() == 0 ) {
							found = true;
							break;
						}
					}
					if ( !found ) {
						merged_vecs_.push_back(add_index_vec);
					}
				}
				add_index_vec.pop_back();
				continue;
			} else if ( v.size() == 1 ) {
				add_index_vec.push_back( v.front() );
				bool found( false );
				std::sort( add_index_vec.begin(), add_index_vec.end() );
				for ( auto & merged_vec : merged_vecs_ ) {
					//std::sort(merged_vec.begin(), merged_vec.end()); // SHOULD ALWAYS ALEADY BE SORTED!
					std::vector< int > v2(100);
					auto it2 = std::set_symmetric_difference(merged_vec.begin(), merged_vec.end(), add_index_vec.begin(), add_index_vec.end(), v2.begin());
					v2.resize(it2-v2.begin());

					if ( v2.size() == 0 ) {
						found = true;
						break;
					}
				}
				if ( !found ) {
					merged_vecs_.push_back(add_index_vec);
				}
				add_index_vec.pop_back();
			} else {
				rec_set_intersection( add_index_vec, v, *j );
			}
			add_index_vec.pop_back();
		}
	}
	if ( TR.visible() ) TR << "finished recursion " << std::endl;
	if ( TR.visible() ) TR << "merged_vecs_.size() =  " << merged_vecs_.size() << std::endl;
	if ( TR.visible() ) TR << "finished temp_merge_vec, actually merging networks " << std::endl;
	Size mcount(1);
	std::vector< HBondNetStructOP > temp_net_vec(0);
	for ( auto & merged_vec : merged_vecs_ ) {
		if ( merged_vec.size()>1 ) {
			auto m1 = merged_vec.begin();
			HBondNetStructOP temp_HBondNetStruct = network_vector_[*m1];
			std::string nets_being_merged("");
			if ( verbose_ ) {
				nets_being_merged = "       "+utility::to_string(mcount)+". branching together networks = ";
				nets_being_merged = nets_being_merged+" ( #"+(utility::to_string(*m1+1))+": ";\
					std::string tempstr = print_list_to_string( *temp_HBondNetStruct, false, false, false, false);
				//std::string tempstr = print_list_to_string(temp_HBondNetStruct.rotlist, temp_HBondNetStruct.term_w_start,
				//                                           temp_HBondNetStruct.term_w_cycle, temp_HBondNetStruct.term_w_bb );
				nets_being_merged = nets_being_merged+tempstr+" ) + ";
			}
			for ( auto m2 = m1; ++m2 != merged_vec.end(); ) {
				if ( verbose_ ) {
					nets_being_merged = nets_being_merged+" ( #"+(utility::to_string((*m2+1)))+": ";
					std::string tempstr = print_list_to_string( *(network_vector_[*m2]), false, false, false, false);
					//std::string tempstr = print_list_to_string(network_vector_[*m2].rotlist, network_vector_[*m2].term_w_start,
					//                                           network_vector_[*m2].term_w_cycle, network_vector_[*m2].term_w_bb);
					nets_being_merged = nets_being_merged+tempstr+" ) + ";
				}
				HBondNetStructOP new_network( pointer::make_shared< HBondNetStruct >() );
				merge_2_branched_networks( *temp_HBondNetStruct, *(network_vector_[*m2]), new_network );
				temp_HBondNetStruct = new_network;
			}
			if ( verbose_ && TR.visible() ) {
				TR << nets_being_merged << std::endl;
			}

			temp_net_vec.push_back( temp_HBondNetStruct );
			mcount++;
		}
	}
	//push the merged networks to the back of network_vector_
	//network_vector_.clear(); DON'T WANT TO CLEAR HERE
	for ( auto & netit : temp_net_vec ) {
		network_vector_.push_back( netit );
	}
	// Sort all individual networks by best score, so start with highest scoring
	std::sort( network_vector_.begin(), network_vector_.end(), compare_net_vec() );
	merged_vecs_.clear();
}//branch_overlapping networks

// used by branch_overlapping() to efficiently search for all combinations of compatible networks that can be merged
void
HBNet::rec_set_intersection( std::vector< Size > add_index_vec, std::vector< Size > next_index_vec, Size pos )
{
	for ( auto i = next_index_vec.begin(); i != next_index_vec.end(); ++i ) {
		if ( *i <=  pos  ) {
			continue;
		}
		add_index_vec.push_back( *i );
		if ( store_subnetworks_ ) {
			bool found( false );
			std::sort( add_index_vec.begin(), add_index_vec.end() );
			for ( auto & merged_vec : merged_vecs_ ) {
				//std::sort(merged_vec.begin(), merged_vec.end()); // SHOULD ALWAYS ALEADY BE SORTED!
				std::vector< int > v2(100);
				auto it2 = std::set_symmetric_difference(merged_vec.begin(), merged_vec.end(), add_index_vec.begin(), add_index_vec.end(), v2.begin());
				v2.resize(it2-v2.begin());

				if ( v2.size() == 0 ) {
					found = true;
					break;
				}
			}
			if ( !found ) {
				merged_vecs_.push_back(add_index_vec);
			}
		}
		std::sort( next_index_vec.begin(),  next_index_vec.end());
		std::sort( (network_vector_[*i])->net_indices.begin(), (network_vector_[*i])->net_indices.end() );
		std::vector< Size > v(1000);
		std::vector< Size >::iterator it;
		it=std::set_intersection(next_index_vec.begin(), next_index_vec.end(), (network_vector_[*i])->net_indices.begin(),
			(network_vector_[*i])->net_indices.end(), v.begin());
		v.resize(it-v.begin());
		if ( v.size() == 0 ) {
			if ( !store_subnetworks_ ) {
				bool found( false );
				std::sort( add_index_vec.begin(), add_index_vec.end() );
				for ( auto & merged_vec : merged_vecs_ ) {
					//std::sort(merged_vec.begin(), merged_vec.end()); // SHOULD ALWAYS ALEADY BE SORTED!
					std::vector< int > v2(100);
					auto it2 = std::set_symmetric_difference(merged_vec.begin(), merged_vec.end(), add_index_vec.begin(), add_index_vec.end(), v2.begin());
					v2.resize(it2-v2.begin());

					if ( v2.size() == 0 ) {
						found = true;
						break;
					}
				}
				if ( !found ) {
					merged_vecs_.push_back(add_index_vec);
				}
			}
			add_index_vec.pop_back();
			continue;
		} else if ( v.size() == 1 ) {
			add_index_vec.push_back(v.front());
			bool found( false );
			std::sort( add_index_vec.begin(), add_index_vec.end() );
			for ( auto & merged_vec : merged_vecs_ ) {
				//std::sort(merged_vec.begin(), merged_vec.end()); // SHOULD ALWAYS ALEADY BE SORTED!
				std::vector< int > v2(100);
				auto it2 = std::set_symmetric_difference(merged_vec.begin(), merged_vec.end(), add_index_vec.begin(), add_index_vec.end(), v2.begin());
				v2.resize(it2-v2.begin());

				if ( v2.size() == 0 ) {
					found = true;
					break;
				}
			}
			if ( !found ) {
				merged_vecs_.push_back(add_index_vec);
			}
			add_index_vec.pop_back();
		} else {
			rec_set_intersection( add_index_vec, v, *i );
		}
		add_index_vec.pop_back();
	}
}

///@details heart of MC-HBNet: use the MC seeds to stochastically grow networks to completion
void HBNet::MC_branch_overlapping_networks(){

	TR << "starting monte carlo branching protocol" << std::endl;

	utility::vector0< std::list< utility::vector1< core::Size > > > hash_table;//each list holds a bunch of vectors of rotamer ids whose sum adds up to a number that defines the hash
	//the goal of this hash table is to prevent us form discovering the same network twice
	hash_table.resize( 839 ); //arbitrary large prime number
	std::list< NetworkState > designed_networks;

	//Data for the purpose of tracer output:
	core::Size next_milestone = total_num_mc_runs_ / 10;
	short int next_milestone_title = 1;

	core::Size const num_seeds = monte_carlo_seeds_.size();
	if ( num_seeds == 0 ) {
		TR << "MC HBNet could not find any seed hydrogen bonds to branch from." << std::endl;
		return;
	}

	for ( core::Size run_id = 1; run_id <= total_num_mc_runs_; ++run_id ) {
		core::Size const rand_index = 1 + numeric::random::rg().uniform() * num_seeds;
		AtomLevelHBondEdge const * const starting_state = monte_carlo_seeds_[ rand_index ];

		if ( ! --next_milestone ) {
			TR << "\t" << next_milestone_title << "0% done with branching" << std::endl;
			next_milestone = total_num_mc_runs_ / 10;
			++next_milestone_title;
		}

		NetworkState current_network_state( starting_state, hbond_graph_ );

		AtomLevelHBondNode * next_node = get_next_node( current_network_state );
		while ( next_node ) {//keep looping until there are no more nodes to add
			add_residue_to_network_state( current_network_state, next_node );

			if ( symmetric_ || network_state_is_satisfied( current_network_state ) ) {
				register_new_network( current_network_state, hash_table, designed_networks );
			}
			next_node = get_next_node( current_network_state );
		}//while true

		if ( symmetric_ || network_state_is_satisfied( current_network_state ) ) {
			//symmetric poses get a free pass from these filters because it is not obvious yet if the networks are satisfied
			register_new_network( current_network_state, hash_table, designed_networks );
		}

	}//run_id

	TR << "Monte carlo branching protocol found " << designed_networks.size() << " networks." << std::endl;

	utility::vector1< DecoratedNetworkState > designed_networks_for_sorting;
	designed_networks_for_sorting.reserve( designed_networks.size() );

	for ( NetworkState const & ns : designed_networks ) {
		core::Real const saturation_est = estimate_saturation( ns );

		unsigned int n_unsat_Hpols = 0;//assumes no heavy unsats
		for ( auto const & uint_vec_pair : ns.unsatisfied_sc_atoms_const() ) {
			n_unsat_Hpols += uint_vec_pair.second.size();
		}
		designed_networks_for_sorting.emplace_back( saturation_est, n_unsat_Hpols, & ns );
	}
	std::sort( designed_networks_for_sorting.begin(), designed_networks_for_sorting.end() );

	append_to_network_vector( designed_networks_for_sorting );
	//std::sort( network_vector_.begin(), network_vector_.end(), compare_net_vec() );
}

void add_network_resids_to_pose( core::pose::Pose & pose, utility::vector1< core::Size > resids, std::string name_of_subset="HBNet" )
{
	using namespace core::select::residue_selector;

	ResidueSubsetOP selected_residues = pointer::make_shared< utility::vector1< bool > > ( pose.size(), false );
	for ( core::Size resid : resids ) {
		assert( resid <= pose.size() );
		(*selected_residues)[ resid ] = true;
	}

	residue_selectors::StoreResidueSubsetMover store_subset_mover;
	store_subset_mover.set_subset( selected_residues, name_of_subset, true );
	store_subset_mover.apply( pose );
}


utility::vector1< core::Size >
HBNet::place_rots_on_pose( pose::Pose & pose, HBondNetStruct & i, bool use_pose )
{
	if ( !( pose.pdb_info() ) ) {
		pose.pdb_info( pointer::make_shared< core::pose::PDBInfo >( pose, true ) );
	}

	utility::vector1< core::Size > resids;
	resids.reserve( i.rotamers.size() );
	if ( !(i.rotamers.empty()) ) {
		for ( const auto & r : i.rotamers ) {
			resids.push_back( r->seqpos() );
			pose.replace_residue( r->seqpos(), *r, false );
			pose.pdb_info()->add_reslabel( r->seqpos(), "HBNet" );
		}
	} else if ( use_pose ) {
		for ( const auto & residue : i.residues ) {
			resids.push_back( residue->resnum );
			pose.replace_residue( residue->resnum, orig_pose_->residue( residue->resnum ), false );
			pose.pdb_info()->add_reslabel( residue->resnum, "HBNet" );
		}
	} else if ( rotamer_sets_ != nullptr ) {
		for ( const auto & residue : i.residues ) {
			resids.push_back( residue->resnum );
			ResidueCOP copy_rot(rotamer_sets_->rotamer_for_moltenres(rotamer_sets_->resid_2_moltenres((platform::uint)(residue->resnum)),residue->rot_index));
			pose.replace_residue( residue->resnum, *copy_rot, false );
			pose.pdb_info()->add_reslabel( residue->resnum, "HBNet" );
		}
	}
	add_network_resids_to_pose( pose, resids, "HBNet" );
	//    if ( place_waters ){
	//        place_bridging_waters_on_pose( pose, i );
	//        if ( pack_waters && !(i.bridging_water_oxygens.empty()) ){
	//            optimize_waters( pose, i );
	//        }
	//    }

	if ( store_network_scores_in_pose_ ) {
		core::pose::setPoseExtraScore( pose, "HBNet_NumUnsatHpol", core::Real( i.num_unsat_Hpol ) );
		core::pose::setPoseExtraScore( pose, "HBNet_Saturation", i.percent_hbond_capacity );
		core::pose::setPoseExtraScore( pose, "HBNet_Score", i.score );
	}

	pose.update_residue_neighbors();
	//pose.update_pose_chains_from_pdb_chains();

	return resids;
}

// to very quickly eliminate networks with 1 or more unsatisfied heavy-atom donors or acceptors
bool
HBNet::quick_and_dirty_network_has_heavy_atom_unsat( Pose const & pose, HBondNetStruct const & network )
{
	for ( const auto & res : network.residues ) {
		if ( res_is_core( res->resnum ) ) {
			for ( const auto & acc_pos : pose.residue(res->resnum).accpt_pos_sc() ) {
				//Size base( pose.residue(res->resnum) );
				if ( quick_and_dirty_heavy_atom_is_unsat( pose, id::AtomID( acc_pos, res->resnum ) ) ) {
					return true;
				}
			}
			//for ( const auto & chemical::AtomIndices : pose.residue(res->resnum).Hpos_polar_sc() ){
			Size prev_base(0);
			for ( const auto & hpol : pose.residue(res->resnum).Hpos_polar_sc() ) {
				Size base( pose.residue(res->resnum).atom_base(hpol) );
				if ( prev_base == base ) continue;
				if ( quick_and_dirty_heavy_atom_is_unsat( pose, id::AtomID( base, res->resnum ) ) ) {
					return true;
				}
				prev_base = base;
			}
		}
	}
	return false;
}

// to very quickly eliminate networks with 1 or more unsatisfied heavy-atom donors or acceptors
bool
HBNet::quick_and_dirty_heavy_atom_is_unsat( Pose const & pose, id::AtomID const at_id )
{
	Real const heavy_atom_don_acc_cutoff(3.5);
	Real const dist_2(heavy_atom_don_acc_cutoff*heavy_atom_don_acc_cutoff);
	for ( Size res = 1; res <= pose.total_residue(); ++res ) {
		if ( res == at_id.rsd() ) continue;
		if ( pose.xyz( at_id ).distance_squared( pose.residue(res).nbr_atom_xyz() ) < 100.0 ) { // 10^2
			if ( pose.residue( at_id.rsd() ).atom_type( at_id.atomno() ).is_acceptor() ) {
				Size prev_base(0);
				for ( const auto & hpol : pose.residue(res).Hpos_polar() ) {
					Size base( pose.residue(res).atom_base(hpol) );
					if ( prev_base == base ) continue;
					if ( pose.xyz( at_id ).distance_squared( pose.residue(res).xyz(base) ) < dist_2 ) {
						return false;
					}
					prev_base = base;
				}
			}
			if ( pose.residue( at_id.rsd() ).atom_type( at_id.atomno() ).is_donor() ) {
				for ( const auto & acc : pose.residue(res).accpt_pos() ) {
					if ( pose.xyz( at_id ).distance_squared( pose.residue(res).xyz(acc) ) < dist_2 ) {
						return false;
					}
				}
			}
		}
	}
	return true;
}

// hacky check based on distance for cases where we only have water O (not explicit H's)
bool
HBNet::atom_hbonds_to_bridging_water( Pose const & pose, id::AtomID const at_id )
{
	//TR << "checking atom_hbonds_to_bridging_water" << std::endl;
	Real const hpol_dist(2.2); //water O (acc) to Hpol dist: generous, error on side of satisfaction
	Real const hpol_dist_sq(hpol_dist*hpol_dist);
	Real const acc_dist(3.4); //water O (don) to Acc dist: generous, error on side of satisfaction
	Real const acc_dist_2(acc_dist*acc_dist);
	for ( Size res = 1; res <= pose.total_residue(); ++res ) {
		//if ( pose.residue(res).is_water() && pose.xyz( at_id ).distance_squared( pose.residue( res ).nbr_atom_xyz()) < 16.0 ) { // 4^2; all water types have neighbor radius < 4
		if ( pose.residue(res).name1() == 'w' && pose.xyz( at_id ).distance_squared( pose.residue( res ).nbr_atom_xyz()) < 16.0 ) { // 4^2; all water types have neighbor radius < 4
			//TR << "we found a water! " << std::endl;
			Vector waterO( pose.residue(res).xyz(1) );//water O should always be first;
			//Vector waterO( pose.residue(res).xyz("O") );
			if ( pose.residue( at_id.rsd() ).atom_type( at_id.atomno() ).is_acceptor() && pose.xyz( at_id ).distance_squared( waterO ) < acc_dist_2 ) {
				//TR << "it H-BONDS!" << std::endl;
				return true;
			} else if ( pose.residue( at_id.rsd() ).atom_is_polar_hydrogen( at_id.atomno() ) && pose.xyz( at_id ).distance_squared( waterO ) < hpol_dist_sq ) {
				//TR << "it H-BONDS!" << std::endl;
				return true;
			}
		}
	}
	//TR << "it DOES NOT h-bond" << std::endl;
	return false;
}

void
HBNet::update_core_and_boundary_residues( core::pose::Pose const & pose )
{
	if ( core_selector_ ) {
		core_residues_ = core_selector_->apply( pose );
	} else {
		core::select::residue_selector::LayerSelectorOP core_layer( pointer::make_shared< core::select::residue_selector::LayerSelector >() );
		core_layer->set_layers( true, false, false );
		//core_layer->set_use_sc_neighbors( true ); // now true by default
		core_layer->set_cutoffs( 4.4 /* core */, 2.0 /* surface */ ); //default core is 5.2, here we use 4.4 to be overly cautious with unsats
		runtime_assert( core_layer->use_sc_neighbors() );
		core_residues_ = core_layer->apply( pose );
	}
}

void
HBNet::find_unsats( Pose const & pose, HBondNetStruct & network )
{
	// THIS SHOULD BE UPDATED AND REPLACED WITH CLEANER CODE

	Size num_unsatisfied_Hpol(0);
	Size num_heavy_unsat(0);

	// update core residues
	//core::select::residue_selector::ResidueSubset temp_core_residues = get_core_residues();
	//core::select::residue_selector::ResidueSubset temp_boundary_residues = get_boundary_residues();
	//if ( bridging_waters_ ){
	//    update_core_and_boundary_residues( pose );
	//}

	runtime_assert(!(network.hbond_vec.empty()));

	if ( verbose_ && TR.visible() ) network.hbond_set->show(TR);

	//TODO NEED TO ADD option to use own bunsat calc here

	// used for percent_hbond_capacity
	Size total_polar_groups_that_could_hbond(0);
	Size polar_groups_making_hbonds(0);

	std::vector< Size > resnums(0);
	resnums.reserve( network.residues.size() );
	for ( utility::vector1< HBondResStructCOP >::const_iterator res = network.residues.begin(); res != network.residues.end(); ++res ) {
		//for ( std::set< Size >::const_iterator res = actual_hbond_residues.begin(); res != actual_hbond_residues.end(); ++res ){
		resnums.push_back( (*res)->resnum );
	}
	//    if ( bridging_waters_ && !(network.waterrots.empty()) ){
	//        Size anchor_res( (symmetric_) ? symm_info_->num_independent_residues() : orig_pose_->total_residue() );
	//        Size new_total( pose.total_residue() );
	//        if ( symmetric_ ){
	//            core::conformation::symmetry::SymmetricConformation const & newSymmConf(dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation()));
	//            new_total = newSymmConf.Symmetry_Info()->num_independent_residues();
	//        }
	//        for ( Size r = anchor_res; r <= new_total; r++ ){
	//            if ( pose.residue(r).is_water() && pose.pdb_info()->res_haslabel(r, "HBNet" ) ){
	//                resnums.push_back(r);
	//            }
	//        }
	//    }

	for ( const auto resnum : resnums ) {

		if ( verbose_ && TR.visible() ) network.hbond_set->show( pose, resnum, true, TR );

		bool skip_backbones(true);
		//bool skip_ligands(true);
		//        if ( !(pose.residue( resnum ).is_water()) ) { // need to handle waters separately
		//for ( Size a_index = 1; a_index <= pose.residue(resnum).natoms(); ++a_index ) {
		for ( Size a_index = 1; a_index <= pose.residue(resnum).nheavyatoms(); ++a_index ) {
			if ( ( skip_backbones && !( pose.residue( resnum ).atom_is_backbone( a_index ) ) ) || ((this->ligand()) && resnum == (this->ligand()) ) ) {
				if ( pose.residue( resnum ).atom_type( a_index ).is_donor() && pose.residue( resnum ).atomic_charge( a_index ) != 0.0 && pose.residue( resnum ).atom_type(a_index).name() != "OH" ) {
					Size h_count(0);
					Size h_unsat(0);
					for ( Size hatm = pose.residue(resnum).attached_H_begin(a_index); hatm <= pose.residue(resnum).attached_H_end(a_index); ++hatm ) {
						//TR.Debug << "hatm = " << pose.residue(resnum).atom_type(hatm).atom_type_name() << std::endl;
						h_count++;
						total_polar_groups_that_could_hbond++;
						id::AtomID const at_id( hatm, resnum);

						if ( !((network.hbond_set)->nhbonds( at_id, false /* include_only_allowed */ ) ) ) {
							//if ( !((network.hbond_set)->nhbonds( at_id, false /* include_only_allowed */ )) &&
							//( !bridging_waters_ || !network.has_bridging_wat || !(network.waterrots.empty()) || !atom_hbonds_to_bridging_water( pose, at_id ) ) ){ // if waterrots !empty, no need to check this because waters optimized and should be in hbond_vec
							//if ( ( !((network.hbond_set)->nhbonds( at_id, false /* include_only_allowed */ )) && ( !bridging_waters_ || !(network.has_bridging_wat) || !(network.waterrots.empty()) ) ) ||
							//( bridging_waters_ && network.has_bridging_wat && network.waterrots.empty() && !atom_hbonds_to_bridging_water( pose, at_id ) ) ){
							if ( res_is_core( resnum ) ) {
								num_unsatisfied_Hpol++;
								h_unsat++;
								network.unsat_Hpols.push_back( at_id );
								TR.Debug << " res " << resnum << " atom " << hatm << " " << pose.residue(resnum).atom_name(hatm) << " is unsat" << std::endl;
							}
						} else {
							polar_groups_making_hbonds++;
						}
					}
					if ( h_unsat == h_count ) {
						num_heavy_unsat++;
						if ( verbose_ && TR.visible() ) TR << "h_unsat == h_count; heavy unsat is residue " << resnum << pose.residue(resnum).name1() << " atom " << a_index << std::endl;
					}
				} else if ( pose.residue( resnum ).atom_type( a_index ).is_acceptor() && pose.residue( resnum ).atomic_charge( a_index ) != 0.0 ) { //important for ligand case
					total_polar_groups_that_could_hbond++;
					// if it's a carbonyl, total polar groups that could hbond is 2 (two lone pairs)
					// if hydroxyl, can both donate and accept, so total groups that could hbond is 2
					//    should be easier way to check if SP2, but atom_type does not have property accessor function like residue_type does
					utility::vector1< std::string > atom_properties = pose.residue(resnum).atom_type( a_index ).get_all_properties();
					if ( pose.residue(resnum).atom_type( a_index ).name() == "OH" || std::find( atom_properties.begin(), atom_properties.end(), "SP2_HYBRID" ) != atom_properties.end() ) {
						total_polar_groups_that_could_hbond++;
					}

					id::AtomID const at_id( a_index, resnum);
					// Careful!  For hydroxyls, this will not count any donor h-bonds; need to check Hpol for that!
					Size num_hbonds = (network.hbond_set)->nhbonds( at_id, false /* include_only_allowed */ );
					polar_groups_making_hbonds += num_hbonds;

					//TR << " " << resnum << pose.residue(resnum).name1() << " atom " << pose.residue(resnum).atom_type( a_index ).name() << "; num_hbonds = " << num_hbonds << std::endl;
					//TR << ": total_polar_groups_that_could_hbond = " << total_polar_groups_that_could_hbond << "; polar_groups_making_hbonds = " << polar_groups_making_hbonds << std::endl;

					if ( ! num_hbonds ) { // num_hbonds == 0
						//                        if ( !((network.hbond_set)->nhbonds( at_id, false /* include_only_allowed */ ))
						//                            && ( !bridging_waters_ || !network.has_bridging_wat || !(network.waterrots.empty()) || !atom_hbonds_to_bridging_water( pose, at_id ) ) ){ // if waterrots !empty, no need to check this because waters optimized and should be in hbond_vec
						//if ( ( !( (network.hbond_set)->nhbonds( at_id , false /* include_only_allowed */ )) && ( !bridging_waters_ || !(network.has_bridging_wat) || !(network.waterrots.empty()) ) ) ||
						//( bridging_waters_ && network.has_bridging_wat && network.waterrots.empty() && !atom_hbonds_to_bridging_water( pose, at_id ) ) ){
						// if res is buried (aka "core") then we penalize it; if bridging waters: expicits checked below, but if don't have them yet, check based on water O
						if ( res_is_core( resnum ) ) {
							if ( pose.residue(resnum).atom_type( a_index ).name() == "OH" ) { //do not want to double penalize hydroxyls
								core::Size hatm = pose.residue(resnum).attached_H_begin( a_index );
								if ( !((network.hbond_set)->nhbonds( id::AtomID( hatm, resnum ), false /* include_only_allowed */ )) ) {
									//                                    if ( !((network.hbond_set)->nhbonds( id::AtomID( hatm, resnum ), false /* include_only_allowed */ )) &&
									//                                        ( !bridging_waters_ || !network.has_bridging_wat || !(network.waterrots.empty()) || !atom_hbonds_to_bridging_water( pose, id::AtomID( hatm, resnum ) ) ) ){ // if waterrots !empty, no need to check this
									//if ( ( !((network.hbond_set)->nhbonds( id::AtomID( hatm, resnum ), false /* include_only_allowed */ )) &&
									//( !bridging_waters_ || !(network.has_bridging_wat) || !(network.waterrots.empty()) ) ) ||
									//( bridging_waters_ && network.has_bridging_wat && network.waterrots.empty() && !atom_hbonds_to_bridging_water( pose, id::AtomID( hatm, resnum ) ) ) ){
									num_unsatisfied_Hpol++;
									num_heavy_unsat++;
									if ( verbose_ && TR.visible() ) TR << "heavy unsat is residue " << at_id.rsd() << pose.residue(at_id.rsd()).name1() << " atom " << a_index << std::endl;
									network.unsat_Hpols.push_back( at_id );
								} else {
									polar_groups_making_hbonds++;
								}
							} else {
								num_unsatisfied_Hpol++;
								num_heavy_unsat++;
								if ( verbose_ && TR.visible() ) TR << "heavy unsat is residue " << at_id.rsd() << pose.residue(at_id.rsd()).name1() << " atom " << at_id.atomno() << std::endl;
								network.unsat_accs.push_back( at_id );
								TR.Debug << " res " << resnum << " atom " << a_index << " " << pose.residue(resnum).atom_name(a_index) << " is heavy unsat" << std::endl;
								TR.Debug << " res " << resnum << " atom " << a_index << " " << pose.residue(resnum).atom_name(a_index) << " is unsat" << std::endl;
								//TR << " res_is_core( resnum ) = " << res_is_core( resnum ) << " and atom_sasa[ at_id ] = " << atom_sasa[ at_id ] << std::endl;
							}
						}
					} else {
						if ( pose.residue(resnum).atom_type( a_index ).name() == "OH" ) {
							core::Size hatm = pose.residue(resnum).attached_H_begin( a_index );
							if ( !((network.hbond_set)->nhbonds( id::AtomID( hatm, resnum ), false /* include_only_allowed */ )) ) {
								//                                if ( !((network.hbond_set)->nhbonds( id::AtomID( hatm, resnum ), false /* include_only_allowed */ )) &&
								//                                    ( !bridging_waters_ || !network.has_bridging_wat || !(network.waterrots.empty()) || !atom_hbonds_to_bridging_water( pose, id::AtomID( hatm, resnum ) ) ) ){ // if waterrots !empty, no need to check this
								if ( res_is_core( resnum ) && ( hydroxyls_must_donate_ || ( pose.residue(resnum).name1() == 'Y' && tyr_hydroxyls_must_donate_ )) ) {
									num_unsatisfied_Hpol++;
									num_heavy_unsat++;
									if ( verbose_ && TR.visible() ) TR << "heavy unsat is residue " << at_id.rsd() << pose.residue(at_id.rsd()).name1() << " atom " << a_index << std::endl;
									network.unsat_Hpols.push_back( at_id );
									TR.Debug << "res " << resnum << " atom " << a_index << " " << pose.residue(resnum).atom_name( a_index ) << " is heavy unsat" << std::endl;
								}
							} else {
								polar_groups_making_hbonds++;
							}
						}
					}
				}
			}
		}
		//        }
		//        else if ( pose.residue(resnum).is_water() && !(network.waterrots.empty()) ) { // need a more clear check: !(network.waterrots.empty() means waters have been packed
		//            //water_count++;
		//            // for waters, require at least 2 h-bonds:
		//            //    this works for current water types, but need more general solution
		//            Size a_index( 1 );
		//            Size num_hbonds_for_water(0);
		//            polar_atom_count++;
		//            //TR << "water count = " << water_count << std::endl;
		//            //TR << "number water O h-bonds" << (network.hbond_set)->nhbonds( id::AtomID( a_index, resnum ), false ) << std::endl;
		//            if ( !((network.hbond_set)->nhbonds( id::AtomID( a_index, resnum ), false /* include_only_allowed */ )) ){
		//                num_unsatisfied_Hpol++;
		//            } else {
		//                num_hbonds_for_water++;
		//                polar_atms_making_hbonds++;
		//            }
		//            polar_atom_count++;
		//            core::Size hatm1 = pose.residue(resnum).attached_H_begin( a_index );
		//            //TR << "number hatm1 h-bonds" << (network.hbond_set)->nhbonds( id::AtomID( hatm1, resnum ), false ) << std::endl;
		//            //Size hatm1( 2 ); // this should work?  atoms numbered starting with 1; not sure if attached_H works for waters yet??
		//            if ( !((network.hbond_set)->nhbonds( id::AtomID( hatm1, resnum ), false /* include_only_allowed */ )) ){
		//                num_unsatisfied_Hpol++;
		//            } else {
		//                num_hbonds_for_water++;
		//                polar_atms_making_hbonds++;
		//            }
		//            polar_atom_count++;
		//            core::Size hatm2 = pose.residue(resnum).attached_H_end( a_index );
		//            //TR << "number hatm2 h-bonds" << (network.hbond_set)->nhbonds( id::AtomID( hatm2, resnum ), false ) << std::endl;
		//            //Size hatm2( 3 ); // this should work?  atoms numbered starting with 1; not sure if attached_H works for waters yet??
		//            if ( !((network.hbond_set)->nhbonds( id::AtomID( hatm2, resnum ), false /* include_only_allowed */ ))) {
		//                num_unsatisfied_Hpol++;
		//            } else {
		//                num_hbonds_for_water++;
		//                polar_atms_making_hbonds++;
		//            }
		//            if ( num_hbonds_for_water < 2 ) num_heavy_unsat++;
		//        }
	}
	network.num_unsat_Hpol = num_unsatisfied_Hpol;
	network.num_heavy_unsat = num_heavy_unsat;
	network.percent_hbond_capacity = (Real)polar_groups_making_hbonds/(Real)total_polar_groups_that_could_hbond; // NOW UPDATED FOR BRIDGING WATERS (CODE ADDED ABOVE)

	//if ( bridging_waters_ ){
	//    set_core_residues( temp_core_residues );
	//    set_boundary_residues( temp_boundary_residues );
	//}
} // find_unsats

//checks resnum and aa but not rotamer
//consider Ser and Thr identical
bool
HBNet::residues_identical( utility::vector1< HBondResStructCOP > & residues1, utility::vector1< HBondResStructCOP > & residues2 )
{
	if ( residues1.size() != residues2.size() ) {
		return false;
	}
	std::sort( residues1.begin(), residues1.end(), compare_hbond_residues() );
	std::sort( residues2.begin(), residues2.end(), compare_hbond_residues() );

	utility::vector1< HBondResStructCOP >::const_iterator res2 = residues2.begin();
	for ( utility::vector1< HBondResStructCOP >::const_iterator res1 = residues1.begin(); res1 != residues1.end(); ++res1 ) {
		if ( ( (*res1)->resnum == (*res2)->resnum ) &&
				( ( ( (*res1)->aa == (*res2)->aa ) ) || ( ( (*res1)->aa == 'S' || (*res1)->aa == 'T' ) && ( (*res2)->aa == 'S' || (*res2)->aa == 'T' ) ) ) ) {
			++res2;
			continue;
		} else {
			return false;
		}
	}
	return true;
}

bool
HBNet::residues_not_unique( utility::vector1< HBondResStructCOP > & residues1, utility::vector1< HBondResStructCOP > & residues2 )
{
	if ( residues1.size() != residues2.size() ) {
		return false;
	}
	std::sort( residues1.begin(), residues1.end(), compare_hbond_residues() );
	std::sort( residues2.begin(), residues2.end(), compare_hbond_residues() );

	utility::vector1< HBondResStructCOP >::const_iterator res2 = residues2.begin();
	for ( utility::vector1< HBondResStructCOP >::const_iterator res1 = residues1.begin(); res1 != residues1.end(); ++res1 ) {
		if ( ( (*res1)->resnum == (*res2)->resnum ) &&
				( ( ( (*res1)->aa == (*res2)->aa ) ) ||
				( ( (*res1)->aa == 'N' || (*res1)->aa == 'Q' ) && ( (*res2)->aa == 'N' || (*res2)->aa == 'Q' ) ) ||
				( ( (*res1)->aa == 'N' || (*res1)->aa == 'D' ) && ( (*res2)->aa == 'N' || (*res2)->aa == 'D' ) ) ||
				( ( (*res1)->aa == 'Q' || (*res1)->aa == 'E' ) && ( (*res2)->aa == 'Q' || (*res2)->aa == 'E' ) ) ||
				( ( (*res1)->aa == 'D' || (*res1)->aa == 'E' ) && ( (*res2)->aa == 'D' || (*res2)->aa == 'E' ) ) ||
				( ( (*res1)->aa == 'S' || (*res1)->aa == 'T' ) && ( (*res2)->aa == 'S' || (*res2)->aa == 'T' ) ) ) ) {
			++res2;
			continue;
		} else {
			return false;
		}
	}
	return true;
}

bool
HBNet::networks_unique( HBondNetStruct const & i, HBondNetStruct const & j, bool no_surface /* true */ )
{
	utility::vector1< HBondResStructCOP > residues_i(0);
	utility::vector1< HBondResStructCOP > residues_j(0);

	if ( no_surface ) {
		for ( const auto & residue : i.residues ) {
			if ( res_is_core(residue->resnum) || res_is_boundary(residue->resnum) ) {
				residues_i.push_back( residue );
			}
		}
		for ( const auto & residue : j.residues ) {
			if ( res_is_core(residue->resnum) || res_is_boundary(residue->resnum) ) {
				residues_j.push_back( residue );
			}
		}
	} else {
		residues_i = i.residues;
		residues_j = j.residues;
	}
	return !(residues_not_unique( residues_i, residues_j ));
}

// 02/28/15 changing behavior, will only keep 1 network with unique resnum/aa identity; networks with same resnum/aa seq but unique rotamers
// considered through final scoring
//   Not releveant for MC-HBNet
void
HBNet::remove_replicate_networks( Size same_max/*=1*/ )
{
	std::sort( network_vector_.begin(), network_vector_.end(), compare_net_vec() ); //sort all networks by score
	std::vector< HBondNetStructOP > temp_net_vec(0);
	//std::vector< Size > net_count(0);
	if ( verbose() && TR.visible() ) TR << "         DELETING REPLICATE NETWORKS THAT HAVE IDENTIAL RESNUM AND AA SEQ; ONLY KEEPING BEST SCORING ONE" << std::endl;
	for ( auto & i : network_vector_ ) {
		bool only_water(true);
		Size num_protein_res(0);
		for ( utility::vector1< HBondResStructCOP >::const_iterator k = i->residues.begin(); k != i->residues.end(); ++k ) {
			if ( (*k)->is_protein || (*k)->is_ligand ) {
				num_protein_res++;
				if ( num_protein_res>=2 ) {
					only_water = false;
					break;
				}
			}
		}
		if ( only_water ) {
			continue;
		}
		bool reached_limit = false;
		//Size temp_net_count = 0;
		Size temp_net_count = 1;
		for ( auto vit = temp_net_vec.begin(); vit != temp_net_vec.end(); ++vit ) {
			if ( residues_identical( i->residues, (*vit)->residues ) ) { //&& ( !bridging_waters_ || ( i->bridging_water_oxygens.size() == (*vit)->bridging_water_oxygens.size() ) ) ){
				//net_count[temp_net_count] = net_count[temp_net_count] + 1;
				temp_net_count++;
				//if ( net_count[temp_net_count] > same_max ) {
				if ( temp_net_count > same_max ) {
					reached_limit = true;
					// don't want break here because need to go through all
				}
			} else if ( only_native_ && is_sub_residues( i->residues, (*vit)->residues ) ) { // checks subsets both ways
				// for native case, ensure that we don't store subnetworks
				if ( i->residues.size() > (*vit)->residues.size() ) {
					**vit = *i;
					auto vit2 = vit+1;
					for ( ; vit2 != temp_net_vec.end(); ) {
						if ( is_sub_residues( i->residues, (*vit2)->residues ) ) {
							vit2 = temp_net_vec.erase(vit2);
						} else {
							++vit2;
						}
					}
				}
				reached_limit = true;
				break;
			}
		}
		if ( !reached_limit ) {
			temp_net_vec.push_back( i );
		} else if ( verbose_ && TR.visible() ) {
			TR << "reached limit, continuing: network " << print_list_to_string(*i) << std::endl;
		}
	}
	//clear network_vector_ then push back all of the ones that passed the cutoff (from temp_net_vec)
	network_vector_.clear();
	for ( auto & k : temp_net_vec ) {
		network_vector_.push_back( k );
	}
	if ( verbose_ && TR.visible() ) {
		TR << " AFTER REMOVING REPLICATES: " << network_vector_.size() << " NETWORKS" << std::endl;
		TR << " ========================================================" << std::endl;
	}
}//remove_replicate_networks

Size
HBNet::get_num_native_rot( Pose & pose, utility::vector1< HBondResStructCOP > const & residues, Real sc_rmsd_cut, bool super )
{
	Size num_native(0);
	for ( const auto & residue : residues ) {
		if ( ( (this->ligand()) && residue->resnum == (this->ligand()) ) || residue->resnum > orig_pose_->total_residue() ) {
			continue; //does not count ligand in order to be fair
		}
		if ( !( pose.residue(residue->resnum).name1() == orig_pose_->residue(residue->resnum).name1() ) ) {
			if ( (pose.residue(residue->resnum).name1() != 'S' && pose.residue(residue->resnum).name1() != 'T') || (orig_pose_->residue(residue->resnum).name1() != 'S' && orig_pose_->residue(residue->resnum).name1() != 'T') ) {
				continue;
			} else {
				num_native++; // for benchmarking purposes, treat S/T same
			}
			continue;
		}
		//Real sc_rmsd( core::scoring::automorphic_rmsd( pose.residue((*res)->resnum), orig_pose_->residue((*res)->resnum), super /*superpose*/ ) );
		Real sc_rmsd( core::scoring::residue_sc_rmsd_no_super( pose.residue(residue->resnum).clone(), orig_pose_->residue(residue->resnum).clone(), super /*final group only*/ ) );
		if ( sc_rmsd <= sc_rmsd_cut ) {
			num_native++;
		}
	}
	return num_native;
}


Size
HBNet::get_num_native_seq( core::pose::Pose & pose, utility::vector1< HBondResStructCOP > const & residues )
{
	Size num_native(0);
	for ( const auto & residue : residues ) {
		if ( ( (this->ligand()) && residue->resnum == (this->ligand()) ) || residue->resnum >= orig_pose_->total_residue() ) {
			continue; //does not count ligand in ligand() case
		}
		bool aa_same = false;
		if ( pose.residue(residue->resnum).name1() == orig_pose_->residue(residue->resnum).name1() ) {
			aa_same = true;
		} else if ( ( pose.residue(residue->resnum).name1() == 'S' || pose.residue(residue->resnum).name1() == 'T' )
				&& ( orig_pose_->residue(residue->resnum).name1() == 'S' || orig_pose_->residue(residue->resnum).name1() == 'T' ) ) {
			aa_same = true;
		}
		if ( aa_same ) {
			num_native++;
		}
	}
	return num_native;
}

void
HBNet::output_networks( bool finalize )
{
	if ( finalize ) {
		//output_net_vec_.clear(); //if going to finalize, clear output_net_vec
		output_vector_.clear();
	}
	if ( TR.visible() ) TR << " Designed these new networks that meet your criteria: " << std::endl << "\t" << print_headers() << print_additional_headers() << std::endl;
	std::sort( network_vector_.begin(), network_vector_.end(), compare_net_vec() ); //sort all new networks
	Size count(1);
	//for ( auto i = network_vector_.begin(); i != network_vector_.end(); ++i ) { // outer loop
	auto i = network_vector_.begin();
	for ( ; i != network_vector_.end(); ) { // outer loop
		// for this options case, need to check if existing (native) networks have been extended
		//  if they have, choose either native or extended, whichever is better
		if ( extend_existing_networks_ ) {
			bool erase_this_network = false;
			auto nit( native_networks_.begin() );
			for ( ; nit != native_networks_.end(); ) { // inner loop
				//bool dummy(false);
				// if the native network is a subset of network i, then i is extended;
				//if ( is_sub_residues( (*nit)->residues, (*i)->residues, dummy, true /* return true_if_identical (networks same size) */) ){
				if ( is_sub_residues( (*nit)->residues, (*i)->residues ) ) {
					// if extended is bigger (and meets criteria) or same size but few unsats, then erase native
					//   bool operator< overridden by HBondNetStruct to define best network
					Size i_size( (symmetric_ && !((*i)->asymm_residues.empty()) ) ? (*i)->asymm_residues.size() : (*i)->residues.size() );
					Size nit_size( (symmetric_ && !(*nit)->asymm_residues.empty()) ? (*nit)->asymm_residues.size() : (*nit)->residues.size() );
					if ( i_size > nit_size || ( i_size == nit_size && (*i)->num_unsat_Hpol < (*nit)->num_unsat_Hpol ) ) {
						(*i)->is_extended = true; // record that this is an extended network
						nit = native_networks_.erase( nit );
					} else {
						// keep native and don't keep new network; we won't add it to output_net_vec_
						// skip this network (it won't get returned)
						erase_this_network = true;
						break; //because we found a native that's better, exit inner loop
					}
				} else {
					++nit;
				}
			}
			//logic for determing if network is native or extended from staring HBNet PDBInfoLabels:
			if ( erase_this_network ) {
				i = network_vector_.erase(i);
				continue;
			} // move on to the next new network, continue outer loop
		}

		Pose ala_copy = *ala_pose_; // deep copy
		//place_rots_on_pose( ala_copy, **i, (*i)->is_native, bridging_waters_, (*i)->waterrots.empty() );
		place_rots_on_pose( ala_copy, **i, (*i)->is_native );
		if ( (*i)->hbond_vec.empty() ) {
			(*scorefxn_)(ala_copy); // NEED TO SCORE BEFORE get_hbond_atom_pairs()
			get_hbond_atom_pairs( **i, ala_copy );
			(*i)->total_hbonds = (*i)->hbond_vec.size();
		}
		(*i)->id = count;
		std::string outstring( ( pdb_numbering() ) ? print_network_w_pdb_numbering( ala_copy, **i, true ) : print_network( **i ) );
		outstring = outstring + this->print_additional_info_for_net( **i, ala_copy );
		(*i)->outstring = outstring;

		if ( write_network_pdbs_ ) write_network_pdb( *i );

		if ( TR.visible() ) TR << "\t" << outstring << std::endl;
		if ( finalize ) {
			//std::vector< Size > net_ids(0);
			//net_ids.push_back( (*i)->id );
			std::set< Size > net_ids;
			net_ids.insert( (*i)->id );

			if ( keep_existing_networks_ ) {
				if ( (*i)->is_extended ) {
					for ( auto j = i; ++j != network_vector_.end(); ) {
						if ( (*j)->is_extended ) {
							bool compatible( true );
							for ( auto & net_id : net_ids ) {
								if ( net_clash( *(get_network_by_id(net_id)), **j ) ) {
									compatible = false;
									break;
								}
							}
							if ( compatible ) net_ids.insert( (*j)->id );
						}
					}
				} else {
					for ( auto & k : network_vector_ ) {
						if ( k->is_extended ) {
							bool compatible( true );
							for ( const auto & net_id : net_ids ) {
								if ( net_clash( *(get_network_by_id(net_id)), *k ) ) {
									compatible = false;
									break;
								}
							}
							if ( compatible ) net_ids.insert( k->id );
						}
					}

				}
			}
			bool already_stored(false);
			for ( const auto & output_vec : output_vector_ ) {
				if ( output_vec == net_ids ) {
					already_stored = true;
				}
			}
			if ( !already_stored ) {
				output_vector_.push_back( net_ids );
			}
			//output_net_vec_.push_back( *i );
		}
		count++;
		i++;
	}
	if ( find_native_ || keep_existing_networks_ ) {
		if ( TR.visible() && native_networks_.size() > 0 ) {
			TR << std::endl << " Keeping these existing (native) networks that meet your criteria: "
				<< std::endl << "\t" << print_headers() << print_additional_headers() << std::endl;
		}
		std::sort( native_networks_.begin(), native_networks_.end(), compare_net_vec() ); //sort native networks
		Size count(1);
		std::string comment_str( print_headers() + this->print_additional_headers() );
		Pose & out_pose( ( output_poly_ala_background_ ) ? *ala_pose_ : *orig_pose_ );
		for ( auto & native_network : native_networks_ ) {
			native_network->id = count;

			// if these options are set, we want to retain the native networks and csts in each pose that HBNet passes back
			//   call private nonconst_get_orig_pose() and modify the original from which all other poses make deep copies before being modified and passed back
			if ( keep_existing_networks_ || only_native_ ) {
				//add_reslabels_to_pose( nonconst_get_orig_pose(), *native_network );
				//place_rots_on_pose( out_pose, *native_network, native_network->is_native, bridging_waters_, native_network->waterrots.empty() );
				place_rots_on_pose( out_pose, *native_network, native_network->is_native );
				if ( native_network->hbond_vec.empty() ) {
					(*scorefxn_)(out_pose); // NEED TO SCORE BEFORE get_hbond_atom_pairs()
					get_hbond_atom_pairs( *native_network, out_pose );
					find_unsats( out_pose, *native_network );
					native_network->total_hbonds = native_network->hbond_vec.size();
				}
				if ( !native_network->scored ) {
					score_network_on_pose( out_pose, *native_network  );
				}
				core::scoring::constraints::ConstraintSetOP cst_op( out_pose.constraint_set()->clone() ); // get existing csts and add to them
				//set_constraints( out_pose, *cst_op, native_network, write_cst_files_ && !bridging_waters_ );
				set_constraints( out_pose, *cst_op, * native_network, write_cst_files_ );
				out_pose.constraint_set(cst_op); // add constraints to the pose
			}
			std::string outstring( ( pdb_numbering() ) ? print_network_w_pdb_numbering( out_pose, *native_network, true ) : print_network( *native_network ) );
			outstring = outstring + this->print_additional_info_for_net( *native_network, out_pose );
			native_network->outstring = outstring;
			if ( TR.visible() ) TR << "\t" << outstring << std::endl;
			comment_str = comment_str + "\n" + native_network->outstring;

			if ( write_network_pdbs_ ) write_network_pdb( native_network );

			count++;
		}
		if ( keep_existing_networks_ || only_native_ ) {
			core::pose::add_comment( out_pose, "HBNet native networks: ", "\n" + comment_str + "\n" );
		}
		if ( TR.visible() ) {
			TR << std::endl;
		}
	}
}// output_networks

core::pack::task::PackerTaskOP
HBNet::create_ptask(core::pose::Pose const & pose, bool initialize_from_commandline/*=false*/)
{
	using namespace core::pack::task;
	TR<<" Creating packer task based on specified task operations..."<< std::endl;
	runtime_assert( task_factory_ != nullptr );
	if ( initialize_from_commandline ) {
		task_factory_->push_back( pointer::make_shared< operation::InitializeFromCommandline >() );
	}
	PackerTaskOP task = task_factory_->create_task_and_apply_taskoperations( pose );
	for ( core::Size resid = 1; resid <= pose.size(); ++resid ) {
		if ( task->being_designed( resid ) ) {
			task->nonconst_residue_task( resid ).allow_aa( core::chemical::aa_gly );
		}
	}
	return task;
}

bool
HBNet::task_is_valid(Pose const & pose) const
{
	if ( task_->total_residue() != pose.total_residue() ) return false;
	for ( Size i(1); i <= pose.total_residue(); ++i ) {
		chemical::ResidueTypeCOP r = pose.residue_type_ptr(i);
		if ( ! task_->residue_task(i).is_original_type( r ) ) return false;
	}
	return true;
}

void
HBNet::set_symmetry( Pose & pose )
{
	//add_background_energies_ = 0; //we do not want to store background and intrares energies at 1-body for symmetry
	//because 2-body energies of symm clones with themselves stored here

	init_scorefxn_ = core::scoring::symmetry::symmetrize_scorefunction( *init_scorefxn_ );
	scorefxn_ = core::scoring::symmetry::symmetrize_scorefunction( *scorefxn_ );

	if ( !core::pose::symmetry::is_symmetric(pose) ) {
		return;
	}
	if ( TR.visible() ) TR << "         POSE IS SYMMETRIC" << std::endl;
	symmetric_ = true;
	auto const & SymmConf(dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation()));
	symm_info_ = SymmConf.Symmetry_Info();
	Size num_components(symm_info_->get_num_components());
	multi_component_ = num_components >= 2;
	if ( TR.visible() ) {
		TR << "         SETTING UP SYMMETRY_INFO: # independent residues = " << symm_info_->num_independent_residues()
			<< "; # bb clones = " << symm_info_->num_bb_clones() << "; # total symm interfaces = " << symm_info_->num_interfaces() << std::endl;
	}
	if ( multi_component_ ) {
		if ( TR.visible() ) TR << "         DETECTED MULTICOMPONENT SYMMETRY: # sym_dof interfaces = " << core::pose::symmetry::sym_dof_names( pose ).size() << std::endl;
	}
	for ( Size ic = 1; ic <= pose.conformation().num_chains(); ++ic ) {
		Size ic_begin = pose.conformation().chain_begin(ic);
		Size ic_end = pose.conformation().chain_end(ic);
		char chain = pose.chain(ic_begin);
		chain_bounds_[chain].first = ic_begin;
		chain_bounds_[chain].second = ic_end;
	}
}

Size
HBNet::get_ind_res( Pose const & pose, Size const res_i)
{
	Size resi_ind(res_i);
	if ( symmetric_ && res_i > symm_info_->num_independent_residues() ) {
		char resi_chain = pose.chain(res_i);
		if ( multi_component_ ) {
			//symm_info_->component_lower_bound()
			std::map<char,std::pair<Size,Size> > const & component_bounds = symm_info_->get_component_bounds();
			char resi_comp = symm_info_->get_component_of_residue(res_i);
			resi_ind = res_i - chain_bounds_[ resi_chain ].first + component_bounds.at( resi_comp ).first;
		} else {
			resi_ind = res_i - chain_bounds_[resi_chain].first + 1;
		}
	}
	return resi_ind;
}

//// TODO Need to finish
//void
//HBNet::write_pymol_files( std::string name )
//{
//    utility::io::ozstream pml_output_stream;
//    std::string pml_fname( name );
//    std::string pdb_fname( pml_fname );
//    std::string str_cst = ".cst";
//    pdb_fname.replace(pdb_fname.find(str_cst),str_cst.length(),".pdb");
//    pml_output_stream.open(pml_fname, std::ios_base::out);
//    pml_output_stream << "# pml file for " << pml_fname;
//    //pml_output_stream << p->outstring << std::endl;
//    pml_output_stream << "load " << pdb_fname << std::endl;
//    pml_output_stream << "unset ignore_case," << std::endl;
//    pml_output_stream << "hide everything, all" << std::endl;
//    pml_output_stream << "bg_color white" << std::endl;
//    pml_output_stream << "space cmyk" << std::endl;
//    pml_output_stream << "color grey90, all" << std::endl;
//    pml_output_stream << "show cartoon, all" << std::endl;
//
//    if ( (this->ligand()) ) {
//        pml_output_stream << "select ligand, resi " << (this->ligand()) << std::endl;
//        pml_output_stream << "util.cbam ligand" << std::endl;
//    } else {
//        pml_output_stream << "select start_residues, resi ";
//        for ( std::set< core::Size >::iterator vit = start_res_vec_.begin(); vit != start_res_vec_.end(); ++vit ) {
//            pml_output_stream << *vit << "+";
//        }
//        pml_output_stream << std::endl;
//        //pml_output_stream << "color red, start_residues" << std::endl;
//    }
//    pml_output_stream << "select DesignResidues, resi ";
//    pml_output_stream << std::endl;
//    pml_output_stream << std::endl;
//    //std::string network_residues("select hbond_network, ");
//    //pml_output_stream << network_residues << std::endl;
//    pml_output_stream << "show sticks, hbond_network" << std::endl;
//    pml_output_stream << "util.cbao hbond_network" << std::endl; //color by atom element, with carbon bright orange
//    pml_output_stream << "#drawing dashed lines for each h-bond in network:" << std::endl;
//    //pml_output_stream << pymol_str_stream.str();
//    pml_output_stream.close();
//}

void
HBNet::write_network_pdb( HBondNetStructOP & p ) // better to pass object pointer or reference?
{
	if ( p->network_pdb_written ) return;

	runtime_assert( p->id != 0 );
	Pose ala_copy = *ala_pose_;
	//place_waters_on_pose( ala_copy, p );
	//place_rots_on_pose( ala_copy, *p, p->is_native, bridging_waters_, p->waterrots.empty() );
	place_rots_on_pose( ala_copy, *p, p->is_native );
	std::string current_out_tag = jd2::current_output_name();
	//from jd2 output_tag; 1st pose will be _0001_1.pdb, 0001.pdb if only 1 pose; additional poses with MPM are _0001_2.pdb, _0001_3.pdb...

	std::string net_prefix("");
	if ( p->is_native ) net_prefix = "native";
	else if ( p->is_extended ) net_prefix = "extended";
	std::string ext( ".pdb" );
	std::string pdb_tag( get_file_name( p->id, net_prefix, ext ) );
	utility::file::FileName pdb_out(pdb_tag);
	std::ostringstream oss;
	oss << basic::options::option[ basic::options::OptionKeys::out::prefix ]() << pdb_out.base();
	pdb_out.base( oss.str() );
	if ( basic::options::option[ basic::options::OptionKeys::out::path::all ].user() ) {
		pdb_out.path(basic::options::option[ basic::options::OptionKeys::out::path::all ]().path());
		pdb_out.vol(basic::options::option[ basic::options::OptionKeys::out::path::all ]().vol());
	} else {
		pdb_out.path("");
		pdb_out.vol("");
	}
	pdb_tag = pdb_out.name();

	core::pose::add_comment( ala_copy, "HBNet details: ", "\n" + print_headers() + this->print_additional_headers() + "\n" + p->outstring );

	//to print csts to pdb info
	core::scoring::constraints::ConstraintSetOP cst_op( ala_copy.constraint_set()->clone() ); // get existing csts and add to them
	set_constraints( ala_copy, *cst_op, * p, false );

	ala_copy.dump_pdb(pdb_tag);

	p->network_pdb_written = true;
}

// assumes the network rotamers are already placed onto the pose
void
HBNet::set_constraints( Pose & pose, core::scoring::constraints::ConstraintSet & constraints, HBondNetStruct & p, bool write_cst_file /* false */ )
{
	runtime_assert( !(p.hbond_vec.empty()) );

	utility::io::ozstream cst_output_stream;
	std::string cst_fname("");
	std::string comment_str( print_headers() + this->print_additional_headers() );
	comment_str = "#" + comment_str + "\n" + "#" + p.outstring + "\n";

	if ( write_cst_file  && !(p.cst_file_written ) ) {
		std::string prefix("");
		if ( p.is_native ) prefix = "native_";
		else if ( p.is_extended ) prefix = "extended_";
		cst_fname = get_file_name( p.id, prefix, ".cst" );
		core::pose::add_comment( pose, "HBNet cst filename for "+prefix+"network_"+utility::to_string(p.id), cst_fname );
		cst_output_stream.open(cst_fname, std::ios_base::out);
		cst_output_stream << "# " << cst_fname << std::endl;
		cst_output_stream << "# " << p.outstring << std::endl;
	}
	for ( utility::vector1< HBondCOP >::const_iterator i = p.hbond_vec.begin(); i != p.hbond_vec.end(); ++i ) {
		Size drsd((*i)->don_res());
		Size arsd((*i)->acc_res());
		Size don_hatm((*i)->don_hatm());
		Size datm(pose.residue(drsd).atom_base((int)(don_hatm)));
		Size aatm((*i)->acc_atm());
		// we only want to turn on constraints for h-bonds involving network residues and atoms that will not change/be designed away downstream!
		Real xtol(0.20);
		numeric::xyzVector<core::Real> don_coordinates = pose.residue(drsd).atom(datm).xyz();
		numeric::xyzVector<core::Real> acc_coordinates = pose.residue(arsd).atom(aatm).xyz();
		Real don_acc_dist = don_coordinates.distance(acc_coordinates);
		if ( don_acc_dist <= 3.0 && don_acc_dist >= 2.6 ) {
			don_acc_dist = 2.8; //idealized value
		} else if ( don_acc_dist > 3.0 ) {
			don_acc_dist = (don_acc_dist+2.6)/2.0;
			xtol = don_acc_dist - 2.6;
		}
		id::AtomID acc_id( aatm, arsd );
		id::AtomID don_id( datm, drsd );
		Real lb( don_acc_dist - xtol );
		Real ub( don_acc_dist + xtol );
		Real sd( 0.2 );
		core::scoring::func::FuncFactory func_fact;
		core::scoring::func::FuncOP const dist_func = func_fact.new_func( "BOUNDED" );
		std::ostringstream bounded_vals;
		std::string tag( "hbond" );
		bounded_vals << utility::to_string(lb) << " " << utility::to_string(ub) << " " << utility::to_string(sd) << " " << tag;
		std::istringstream in(bounded_vals.str());
		dist_func->read_data(in);
		//TODO need more control of how strongly this is penalized
		//NEED option for angle constraint
		core::scoring::constraints::AtomPairConstraintCOP hbond_cst( pointer::make_shared< core::scoring::constraints::AtomPairConstraint >( acc_id, don_id, dist_func ) );
		constraints.add_constraint( hbond_cst );
		//(*i)->show(in_pose, 1, TR);
		Real angle_AHD((*i)->get_AHDangle(pose));
		Real angle_BAH((*i)->get_BAHangle(pose));
		Real dist_AH((*i)->get_HAdist(pose));
		Real hb_energy((*i)->energy());//unweighted h-bond energy
		std::ostringstream cst_str_stream;
		cst_str_stream << "# angle_AHD = " << angle_AHD << ", angle_BAH = " << angle_BAH << ", dist_AH = " << dist_AH << angle_BAH << ", unweighted hb_energy = " << hb_energy << "\n";
		cst_str_stream << "#AtomPair " << pose.residue(arsd).atom_name((int)(aatm)) << " " << arsd << " " << pose.residue(drsd).atom_name((int)(datm)) << " " << drsd
			<< " BOUNDED " << lb << " " << ub << " " << sd << " " << tag << std::endl;
		comment_str = comment_str + cst_str_stream.str();
		if ( write_cst_file && !(p.cst_file_written ) ) {
			cst_output_stream << cst_str_stream.str();
		}
	}
	core::pose::add_comment( pose, "HBNet constraints for network_" + utility::to_string(p.id) + ": ", "\n" + comment_str + "\n" );
	if ( write_cst_file  && !(p.cst_file_written ) ) {
		cst_output_stream.close();
		p.cst_file_written = true;
	}
}

std::string
HBNet::get_file_name( Size id, std::string net_prefix, std::string extension )
{
	std::string current_out_tag = jd2::current_output_name();
	std::string fname = current_out_tag+(basic::options::option[ basic::options::OptionKeys::out::suffix ]())+"_"+net_prefix+"network_"+utility::to_string( id )+extension;
	utility::file::FileName outfile(fname);
	std::ostringstream oss;
	oss << basic::options::option[ basic::options::OptionKeys::out::prefix ]() << outfile.base();
	outfile.base( oss.str() );
	if ( basic::options::option[ basic::options::OptionKeys::out::path::all ].user() ) {
		outfile.path(basic::options::option[ basic::options::OptionKeys::out::path::all ]().path());
		outfile.vol(basic::options::option[ basic::options::OptionKeys::out::path::all ]().vol());
	} else {
		outfile.path("");
		outfile.vol("");
	}
	return outfile.name();
}

bool
HBNet::water_clashes(
	Pose const & pose,
	Vector const water_O
)
{
	Real const acc_cutoff( 2.0 );
	Real const acc_cutoff_sq( acc_cutoff * acc_cutoff );

	//NEED TO add check for symmetric waters that clash with themselves
	// hacky check to see if water O is on z-axis (in which case it will clash with symmetric copies of itself)
	if ( symmetric_ && ( (water_O.x()+water_O.y())*(water_O.x()+water_O.y()) < acc_cutoff_sq ) ) return true;
	for ( Size res = 1; res <= pose.total_residue(); ++res ) {
		if ( water_O.distance_squared( pose.residue( res ).nbr_atom_xyz()) > 144.0 ) { /* 12^2 */ // could get by with 10?
			continue;
		}
		if ( water_oxygen_clashes_with_residue( water_O, pose.residue(res) ) ) {
			return true;
		}
	}
	return false;
}

bool
HBNet::water_oxygen_clashes_with_residue( Vector const water_oxygen, Size const resnum, int const rot_state )
{
	//Should be able to look up atom coords in rotamer_sets_; rotamer_building_functions orient and transform coords onto pose backbone
	if ( rotamer_sets_ != nullptr && rotamer_sets_->resid_2_moltenres(resnum) != 0 ) { //returns false if res1 not molten (not design/repack-able)
		return water_oxygen_clashes_with_residue( water_oxygen, *(rotamer_sets_->rotamer_set_for_moltenresidue(rotamer_sets_->resid_2_moltenres(resnum))->rotamer(rot_state)) );
	}
	return false;
}

bool
HBNet::water_oxygen_clashes_with_residue( Vector const water_oxygen, Residue const & res )
{
	//should hard code these constants
	Real const hpol_cutoff( 1.0 );
	Real const hpol_cutoff_sq( hpol_cutoff * hpol_cutoff );
	Real const acc_cutoff( 2.0 );
	Real const acc_cutoff_sq( acc_cutoff * acc_cutoff );
	Real const clash_cutoff( 1.8 );
	Real const clash_cutoff_sq( clash_cutoff * clash_cutoff );

	for ( Size at2 = 1; at2 <= res.natoms(); ++at2 ) {
		//skip virtual atoms but not virtual waters!
		std::string resname( res.name() );
		//we want to compare against all waters, even if VRT (but skip other VRT atoms)
		if ( res.atom_type( at2 ).lj_wdepth() == 0.0 && !(resname == "HOH_V" || resname == "PWAT_V") ) continue; //fd: update following water simplification
		Real const dist2( water_oxygen.distance_squared( res.xyz( at2 ) ) );
		if ( !(res.atom_is_backbone( at2 )) && res.atom_is_polar_hydrogen( at2 ) && dist2 < hpol_cutoff_sq ) {
			return true;
			//        } else if ( ( res.atom_type(at2).is_acceptor() || res.is_water() ) && dist2 < acc_cutoff_sq ) {
		} else if ( ( res.atom_type(at2).is_acceptor() || res.name1() == 'w' ) && dist2 < acc_cutoff_sq ) {
			return true;
		} else if ( dist2 < clash_cutoff_sq ) {
			return true;
		}
	}
	return false;
}

///@details minimal h-bond network information is stored in memory and used here to place the rotamers onto a copy of input pose, returend on-the-fly as requested
/// returns networks in order of ranking and automatically turns on constraints (added to Pose that's returned)
pose::PoseOP
HBNet::get_additional_output()
{
	pose::PoseOP out_pose;
	if ( output_vector_.size() == 0 ) {
		return nullptr;
	} else if ( output_vector_.size() == 1 ) {
		// Last iteration, do not create a new pose copy
		out_pose = ( output_poly_ala_background_ ) ? ala_pose_ : orig_pose_;
	} else {
		out_pose = (
			output_poly_ala_background_ ?
			pointer::make_shared< pose::Pose >( * ala_pose_ ) :
			pointer::make_shared< pose::Pose >( * orig_pose_ ) );
	}
	std::set< Size > net_ids( output_vector_.back() ); // back returs ref to back
	output_vector_.pop_back(); // pop_back() deletes/removes last, does not return;

	std::string comment_str( print_headers() + this->print_additional_headers() );
	if ( net_ids.size() > 1 && TR.visible() ) {
		TR << std::endl << "combining networks ";
	}
	for ( const auto & net_id : net_ids ) {

		if ( net_ids.size() > 1 && TR.visible() ) {
			TR << net_id << "  ";
		}
		//TR << "net_id = " << net_id << std::endl;
		runtime_assert( get_network_by_id(net_id) != nullptr );
		//place_rots_on_pose( *out_pose, *(get_network_by_id( net_id )), (get_network_by_id( net_id ))->is_native, bridging_waters_, (get_network_by_id( net_id ))->waterrots.empty() );
		place_rots_on_pose( *out_pose, *(get_network_by_id( net_id )), (get_network_by_id( net_id ))->is_native );
		//add_reslabels_to_pose( *out_pose, *(get_network_by_id( net_id )) ); // needs to be done on the fly when adding waters
		comment_str = comment_str + "\n" + (get_network_by_id( net_id ))->outstring;
		//have to do this each time if adding residues (waters)
		core::scoring::constraints::ConstraintSetOP cst_op( out_pose->constraint_set()->clone() ); // get existing csts and add to them
		//        if ( bridging_waters_ ){
		//            //NEED TO SCORE BEFORE get_hbond_atom_pairs()
		//            ( *scorefxn_ )( *out_pose ); // TODO NEED BETTER WAY THAN SCORING EACH TIME AND OVERRIDING hbond_vec
		//            get_hbond_atom_pairs( *(get_network_by_id( net_id )), *out_pose ); // if residues (waters) have been added, need to update
		//        }
		set_constraints( *out_pose, *cst_op, * get_network_by_id( net_id ), write_cst_files_ );
		out_pose->constraint_set(cst_op); // add constraints to the pose
	}
	if ( TR.visible() ) TR << std::endl;
	core::pose::add_comment( *out_pose, "HBNet designed networks: ", "\n" + comment_str + "\n" );

	//(*out_pose).update_residue_neighbors(); // now taken care of in place_rots_on_pose()
	( *scorefxn_ )( *out_pose ); // always want to return a scored pose; and set_constraints() clears energies!! so need to re-score before calling get_hbond_atom_pairs()
	return out_pose;
}

void
HBNet::search_IG_for_networks( Pose & )
{
	// traverse IG and enumerate all possible h-bond networks given parameters
	if ( ig_ != nullptr ) {
		if ( TR.visible() ) {
			TR << " ==================================================================" << std::endl;
			TR << " ============     PERFORMING H-BOND NETWORK DESIGN     ============" << std::endl;
			TR << " ==================================================================" << std::endl;
		}

		if ( monte_carlo_branch_ ) {
			MC_traverse_IG( );
		} else {
			traverse_IG( hydrogen_bond_threshold_ );
		}

	}
}

void
HBNet::prepare_output()
{
	output_networks(true); //will add single networks to output vector
}

Size
HBNet::num_core_res( HBondNetStruct const & network ){
	Size num_core(0);
	for ( const auto & residue : network.residues ) {
		if ( res_is_core( residue->resnum ) ) ++num_core;
	}
	return num_core;
}

Size
HBNet::num_boundary_res( HBondNetStruct const & network ){
	Size num_boundary(0);
	for ( const auto & residue : network.residues ) {
		if ( res_is_boundary( residue->resnum ) ) ++num_boundary;
	}
	return num_boundary;
}

void
HBNet::select_best_networks()
{

	// TODO now that Jack's AtomLevelHBondGraph enables high-integrity unsat check on-the-fly,
	//   this function can be refactored to be more efficient, and written in a cleaner fashion

	//min_network_size, min_core_res, and min_boundary_res apply to ASU;
	//  HBNetStapleInterface contains additional options to check for symmetric cases that span entire symmetric interface
	auto i = network_vector_.begin();
	for ( ; i != network_vector_.end(); ) {
		//if ( (*i)->residues.size() < 3 && min_network_size_ > 2 ) {
		if ( (*i)->residues.size() < min_network_size_ || (*i)->residues.size() > max_network_size_ ) {
			//TR << "i = " << print_list_to_string( (*i)->residues ) << "ERASING! Size = " << (*i)->residues.size() << std::endl;
			i = network_vector_.erase( i );
		} else if ( min_core_res_ && num_core_res( **i ) < min_core_res_ ) {
			i = network_vector_.erase( i );
		} else if ( min_boundary_res_ && num_boundary_res( **i ) < min_boundary_res_ ) {
			i = network_vector_.erase( i );
		} else if ( !(network_meets_initial_criteria( **i )) ) {
			i = network_vector_.erase( i );
		} else {
			++i;
		}
	}
	//iterator over all the networks, erase those that do not meet specified criteria
	Pose ala_copy = *ala_pose_;
	utility::vector1< core::Size > previous_resids;

	i = network_vector_.begin();
	for ( ; i != network_vector_.end(); ) {

		//ADD FUNCTION HERE TO CHECK IF NETWORK CONTAINS SPECIFICIED AA TYPES


		//place network residues on background pose for scoring and evaluation
		//place_rots_on_pose( ala_copy, **i, (*i)->is_native, bridging_waters_, pack_waters_before_unsat_check_ );
		for ( auto /*integral*/ resid : previous_resids ) {
			ala_copy.replace_residue( resid, ala_pose_->residue( resid ), false );
		}
#ifndef NDEBUG
		for ( core::Size resid = 1; resid <= ala_pose_->size(); ++resid ) {
			debug_assert(
				ala_copy.residue( resid ).name1() == ala_pose_->residue( resid ).name1()
			);
		}
#endif
		previous_resids = place_rots_on_pose( ala_copy, **i, (*i)->is_native );
		//ala_copy.update_residue_neighbors();

		// very quickly identify networks that have heavy atom donors or acceptors unsatisfied and get rid of them!
		if ( no_heavy_unsats_allowed_ &&  quick_and_dirty_network_has_heavy_atom_unsat( ala_copy, **i ) ) {
			if ( verbose_ && TR.visible() ) TR << "i = " << print_list_to_string( **i ) << "ERASING! did not pas quick_and_dirty_unsat_check" << std::endl;
			i = network_vector_.erase( i );
		} else {
			if ( minimize_ ) {
				minimize_network( ala_copy, **i, true /* residues_already_placed */ );
			}
			// Before get_hbond_atom_pairs(), NEED TO BE SURE POSE IS SCORED and neighbors and energies are updated: residue_neighbors_updated
			if ( !ala_copy.energies().energies_updated() ) {
				(*scorefxn_)(ala_copy);
			}
			get_hbond_atom_pairs( **i, ala_copy ); //default sets bb_exclusion to false

			//if ( (*i)->hbond_vec.empty() || ( bridging_waters_ && (*i)->has_bridging_wat && (*i)->bridging_water_oxygens.empty() ) ) {
			if ( (*i)->hbond_vec.empty() ) {
				if ( verbose_ ) TR << "i = " << print_list_to_string( **i ) << "ERASING! h-bonds empty" << std::endl;
				i = network_vector_.erase( i );
			} else {
				//find all buried unsatisfied polars in the network
				//  num unsats also applies only to ASU
				find_unsats( ala_copy, **i ); //anything < 0.0 considered h-bond during satisfaction check

				//move some of the to meets_criteria NEED TO FIX
				if ( ( no_heavy_unsats_allowed_ && (*i)->num_heavy_unsat > 0 ) || (*i)->num_unsat_Hpol > max_unsat_Hpol_ ) {
					if ( verbose_ && TR.visible() ) {
						TR << "i = " << print_list_to_string( **i ) << "ERASING! num unsat = " << (*i)->num_unsat_Hpol << "; num_heavy_unsat = " << (*i)->num_heavy_unsat << std::endl;
					}
					i = network_vector_.erase( i );
				} else {
					//                    if ( bridging_waters_ && (*i)->has_bridging_wat && !((*i)->bridging_water_oxygens.empty()) && !pack_waters_before_unsat_check_ ){ // pack waters now then update networks
					//                        optimize_waters( ala_copy, **i ); // will score the pose, so get_hbond_atom_pairs() is good
					//                        get_hbond_atom_pairs( **i, ala_copy ); //default sets bb_exclusion to false
					//                        if ( (*i)->hbond_vec.empty() ){
					//                            i = network_vector_.erase( i );
					//                            continue;
					//                        }
					//                        find_unsats( ala_copy, **i );
					//                    }
					if ( (*i)->num_unsat_Hpol > max_unsat_Hpol_ || (*i)->percent_hbond_capacity < min_percent_hbond_capacity_ ||
							( this->ligand() && (*i)->lig_num_unsatisfied > max_lig_unsat_ ) || !(this->network_meets_final_criteria( ala_copy, **i )) ) {
						if ( verbose_ && TR.visible() ) TR << "i = " << print_list_to_string( **i ) << "ERASING! percent_hbond_capacity = " << (*i)->percent_hbond_capacity << "unsat = " << (*i)->num_unsat_Hpol << std::endl;
						i = network_vector_.erase( i );
					} else {
						score_network_on_pose( ala_copy, **i );
						if ( (*i)->score > upper_score_limit_ ) {
							i = network_vector_.erase( i );
						} else {
							++i;
						}
					}
				}
			}
		}
	}
	//score the networks
	if ( verbose_ && TR.visible() ) {
		TR << "NUMBER OF NETWORKS BEFORE SCORE = " << network_vector_.size() << std::endl;
		TR << " SCORING THE NETWORKS: " << std::endl;
	}
	//TODO NEED TO add charge-charge repulsion check
	//score_networks( minimize_ /* false */ ); // now scored above
	if ( min_unique_networks_ > 1 ) {
		Size num_unique(1);
		for ( auto i = network_vector_.begin(); i != network_vector_.end(); ++i ) {
			for ( auto j = i; ++j != network_vector_.end(); ) {
				if ( networks_unique( **i, **j ) ) {
					num_unique++;
				}
				if ( num_unique >= min_unique_networks_ ) {
					break;
				}
			}
			if ( num_unique >= min_unique_networks_ ) {
				break;
			}
		}
		if ( num_unique < min_unique_networks_ ) {
			network_vector_.clear();
		}
	}
} // select_best_networks()

void
HBNet::trim_rotamers( Pose & pose ){
	this->trim_additional_rotamers( pose );
}

void
HBNet::setup_packer_task_and_starting_residues( Pose const & pose )
{
	if ( task_factory_ == nullptr ) {
		task_factory_ = pointer::make_shared< task::TaskFactory >( );
	}
	task_ = create_ptask( pose ); //set task (which residues are designable/repackable
	if ( start_res_vec_.empty() ) {
		utility::vector1<bool> is_repack = task_->repacking_residues();
		runtime_assert(is_repack.size() == pose.total_residue());
		// TODO NEED TO make sure task_ is fully processed and trimmed to asymmetric unit for symmetric cases by this point
		for ( Size r = 1; r <= pose.total_residue(); ++r ) {
			if ( pose.residue(r).is_protein() ) {
				if ( task_->design_residue((int)r) || is_repack[r]==1 ) {
					start_res_vec_.insert(r);
				}
			}
		}
	}
}

void
HBNet::get_native_networks( Pose const & pose )
{
	utility::vector1<Size> residues_to_ala(0);
	for ( Size r = 1; r <= pose.total_residue(); ++r ) {
		if ( pose.residue(r).is_ligand() || pose.residue_type(r).name() == "VRT" ) {
			continue;
		}
		residues_to_ala.push_back(r);
	}
	// set ala_pose_
	ala_pose_ = pointer::make_shared< Pose >( pose );
	toolbox::pose_manipulation::construct_poly_XXX_pose("ALA", *ala_pose_, residues_to_ala, true, true, true);
	ala_pose_->update_residue_neighbors();
	(*scorefxn_)(*ala_pose_); // score now so we don't have to later; several functions require ala_pose_ to be scored

	bool orig_only_native = only_native_;
	only_native_ = true;
	//    bool orig_bridging_waters = bridging_waters_;
	//    if ( ( find_native_ || keep_existing_networks_ ) && !only_native_ ){
	//        bridging_waters_ = false;
	//    }

	traverse_native( pose, hydrogen_bond_threshold_ );
	branch_overlapping_networks();
	remove_replicate_networks();

	// if want to find natives to with criteria specified in HBNet mover defition (independent of input HBNet PDBInfoLabels)
	if ( find_native_ ) {
		select_best_networks();
	} else {
		for ( auto & native_network : native_networks_ ) {
			Pose ala_copy = *ala_pose_;
			//place_rots_on_pose( ala_copy, *native_network, native_network->is_native, bridging_waters_, native_network->waterrots.empty() );
			place_rots_on_pose( ala_copy, *native_network, native_network->is_native );
			if ( native_network->hbond_vec.empty() ) {
				(*scorefxn_)(ala_copy); // NEED TO SCORE BEFORE get_hbond_atom_pairs()
				get_hbond_atom_pairs( *native_network, ala_copy );
				find_unsats( ala_copy, *native_network );
				native_network->total_hbonds = native_network->hbond_vec.size();
			}
			if ( !native_network->scored ) {
				score_network_on_pose( ala_copy, *native_network  );
			}
		}
	}
	only_native_ = orig_only_native;
	//bridging_waters_ = orig_bridging_waters;
	native_networks_ = network_vector_;
	network_vector_.clear();
}

void
HBNet::setup( Pose const & pose )
{
	this->setup_packer_task_and_starting_residues( pose );

	if ( task_->rotamer_links_exist() ) {
		if ( TR.visible() ) TR << " WARNING!!!! ROTAMER LINKS DETECTED: HBNet can't currently handle rotamer links!" << std::endl;
		rotamer_links_ = task_->rotamer_links();
	}
	if ( find_native_ ) {
		if ( start_res_vec_.empty() ) {
			if ( TR.visible() ) TR << "No starting residues defined; traversing the entire static pose:" << std::endl;
			Size total( (symmetric_) ? symm_info_->num_independent_residues() : pose.total_residue() );
			for ( Size res = 1; res <= total; ++res ) {
				start_res_vec_.insert(res);
			}
		}
		get_native_networks( pose );

		if ( only_native_ ) {
			return;
		} else {
			// if we found native networks, add them to the "HBNet" info label residues in case keep_existing or extend_exisitng options set to true
			core::select::residue_selector::ResiduePDBInfoHasLabelSelectorOP hbnet_info_label( pointer::make_shared< core::select::residue_selector::ResiduePDBInfoHasLabelSelector >( "HBNet" ) );
			input_hbnet_info_residues_ = hbnet_info_label->apply( pose );
		}
	}
	if ( only_extend_existing_ ) {
		start_res_vec_.clear();
		for ( Size i = 1; i <= input_hbnet_info_residues_.size(); ++i ) {
			if ( input_hbnet_info_residues_[i] ) {
				start_res_vec_.insert( i );
			}
		}
	}

	core::pack::task::operation::RestrictAbsentCanonicalAASOP set_hbnet_des_res = pointer::make_shared< core::pack::task::operation::RestrictAbsentCanonicalAAS >();
	set_hbnet_des_res->keep_aas( des_residues_ );
	set_hbnet_des_res->include_residue(0);
	set_hbnet_des_res->apply( pose, *task_ );

	//HBNet requires PDInteractionGraph (or DensePDInteractionGraph).  Need to check to make sure NOT using lin_mem_ig
	task_->and_linmem_ig( false );
	if ( ( *task_ ).linmem_ig() && TR.visible() ) {
		TR << " ERROR: EXITING HBNet: You cannot use linmem_ig with HBNet; Please remove linmem_ig from HBNet task_operations." << std::endl;
		TR << " You can define a SetIGType task_operations in your XML and add it to the task_operations of other movers in your protocol, just not HBNet." << std::endl;
	}
	runtime_assert(!(( *task_ ).linmem_ig()));
	runtime_assert( task_is_valid( pose ) );

	if ( show_task_ && TR.visible() ) {
		//task_->show();
		task_->show(TR);
	}
	// Setup IG ig_:
	if ( symmetric_ ) {
		// need to truncate task to independent residues only
		utility::vector1<bool> allow_repacked( pose.total_residue(), false );
		for ( Size res=1; res <= pose.total_residue(); ++res ) {
			if ( pose.residue(res).aa() != core::chemical::aa_vrt && symm_info_->fa_is_independent(res) ) {
				allow_repacked.at(res) = true;
			}
		}
		task_->restrict_to_residues( allow_repacked ); //"vector boolean is based on residue position, disables packing at false positions does nothing to true positions.  Cannot turn on packing.
		rotamer_set::symmetry::SymmetricRotamerSetsOP sym_rotset_op( pointer::make_shared< rotamer_set::symmetry::SymmetricRotamerSets >() );
		rotamer_sets_ = sym_rotset_op;
	} else {
		rotamer_sets_ = RotamerSetsOP( pointer::make_shared< rotamer_set::RotamerSets >() );
	}
	// mechanism to keep/extend networks by HBNet InfoLabels in input pose
	if ( keep_existing_networks_ || extend_existing_networks_ ) {
		//for every position with HBNet info label in input pdb

		std::set< core::Size > temp_start_res_vec_ = start_res_vec_; // save starting residues
		if ( keep_existing_networks_ ) {
			start_res_vec_.clear();
		}
		for ( Size i = 1; i <= input_hbnet_info_residues_.size(); ++i ) {
			if ( input_hbnet_info_residues_[i] ) {
				if ( keep_existing_networks_ ) {
					start_res_vec_.insert(i); // assumes all HBNet InfoLabel residues all in ASU in input pose (always true unless pose manually modified in between HBNet runs)
					if ( !extend_existing_networks_ ) {
						task_->nonconst_residue_task( i ).prevent_repacking();
					}
				}
				if ( extend_existing_networks_ ) {
					task_->nonconst_residue_task( i ).reset(); // in case previously set to prevent_repacking()
					// turning on optH restricts rotamers to only protons flips during packing
					//   this simultaneously allows for
					task_->nonconst_residue_task( i ).or_optimize_h( true );
					task_->nonconst_residue_task( i ).or_include_current( true );
				}
			}
		}
		// if we've already found natives, skip this; else, detect the existing native networks as defined by HBNet InfoLabels and store them for later
		if ( keep_existing_networks_ && !find_native_ ) {
			get_native_networks( pose );
		}
		start_res_vec_ = temp_start_res_vec_;

		if ( find_native_ ) {
			//for each rotamer in native_networks_: set them to not be packable/designable in task_
			for ( auto & native_network : native_networks_ ) {
				for ( auto & residue : native_network->residues ) {
					if ( keep_existing_networks_ && !extend_existing_networks_ ) {
						task_->nonconst_residue_task( residue->resnum ).prevent_repacking();
					} else if ( extend_existing_networks_ ) {
						// do any of these flip it to packable if task_factory alraedy designates it as prevent_repacking()?
						task_->nonconst_residue_task( residue->resnum ).restrict_to_repacking(); //ensure packable
						task_->nonconst_residue_task( residue->resnum ).or_include_current( true );
						task_->nonconst_residue_task( residue->resnum ).sample_proton_chi(true);
					}
				}
			}
		}
	}
	// added with Benjamin Basanta; makes start selector positions fixed
	if ( keep_start_selector_rotamers_fixed_ ) {
		for ( const auto & start_res : start_res_vec_ ) {
			if ( !task_->being_packed( start_res ) ) {
				//TR << "Residue " << start_res << " was set to not repack,
				//only the input rotamer will be used."<< std::endl;
				const chemical::AA aa = task_->nonconst_residue_task( start_res ).get_original_residue();
				task_->nonconst_residue_task( start_res ).reset(); // in case previously set to prevent_repacking()
				task_->nonconst_residue_task( start_res ).allow_aa( aa );
				task_->nonconst_residue_task( start_res ).or_optimize_h( true );
				task_->nonconst_residue_task( start_res ).or_include_current( true );
			}
		}
	}
	utility::vector1<Size>  residues_to_ala(0);
	utility::vector1<bool> is_repack = task_->repacking_residues();
	runtime_assert(is_repack.size() == pose.total_residue());
	for ( Size r = 1; r <= pose.total_residue(); ++r ) {
		if ( pose.residue(r).is_protein() ) {
			if ( task_->design_residue((int)r) || is_repack[r]==1 ) {
				residues_to_ala.push_back(r);
			}
		}
	}
	TR.Debug << " Storing background version of pose; design/repack residues set to Poly-ALA (keeping Pro,Gly,Cys-disulfide):" << std::endl;
	ala_pose_ = pointer::make_shared< pose::Pose >( pose );
	toolbox::pose_manipulation::construct_poly_XXX_pose("ALA", *ala_pose_, residues_to_ala, true, true, true);
	ala_pose_->update_residue_neighbors();
	(*scorefxn_)(*ala_pose_); // score now so we don't have to later; several functions require pose to be scored
}

void
HBNet::run( Pose & pose )
{
	if ( only_native_ ) {
		return; // no need to do design
	}
	//pose.update_residue_neighbors();
	(*init_scorefxn_).setup_for_packing( pose, task_->repacking_residues(), task_->designing_residues() );
	packer_neighbor_graph_ = create_packer_graph( pose, *init_scorefxn_, task_ ); //even in symmetric cases, packer graph will have a node for every residue in the pose

	rotamer_sets_->set_task( task_ );
	rotamer_sets_->initialize_pose_for_rotsets_creation(pose);
	rotamer_sets_->build_rotamers( pose, *init_scorefxn_, packer_neighbor_graph_ );

	//remove unwanted rotamers from rotamer_sets_ before IG is populated and searched.
	trim_rotamers( pose );

	rotamer_sets_->prepare_sets_for_packing( pose, *init_scorefxn_ );
	if ( monte_carlo_branch_ ) {
		initialize_hbond_graph();
		ig_ = pointer::make_shared< MCHBNetInteractionGraph >( hbond_graph_, rotamer_sets_, hydrogen_bond_threshold_, clash_threshold_ );
	} else {
		ig_ = InteractionGraphFactory::create_interaction_graph( *task_, *rotamer_sets_, pose, *init_scorefxn_, *packer_neighbor_graph_ );
	}
	ig_->initialize( *rotamer_sets_ );

	PrecomputedPairEnergiesInteractionGraphOP pig( utility::pointer::dynamic_pointer_cast< PrecomputedPairEnergiesInteractionGraph > ( ig_ ) );
	runtime_assert(pig);

	if ( symmetric_ ) {
		for ( platform::uint ii = 1; ii <= rotamer_sets_->nmoltenres(); ++ii ) {
			utility::vector1< core::PackerEnergy > one_body_energies((rotamer_sets_->rotamer_set_for_moltenresidue(ii))->num_rotamers() );
			hbnet_symm_one_body_energies(pose, * rotamer_sets_->rotamer_set_for_moltenresidue(ii), *init_scorefxn_, *task_, *packer_neighbor_graph_, one_body_energies );
			pig->add_to_nodes_one_body_energy( ii, one_body_energies );
		}
	}
	// precompute_two_body_energies() now vitural and dervied; no longer need to cast to SymmetricRotamerSets first
	rotamer_sets_->precompute_two_body_energies( pose, *init_scorefxn_, packer_neighbor_graph_, pig, true );

	if ( TR.visible() ) {
		TR << " built " << rotamer_sets_->nrotamers() << " rotamers at "
			<< rotamer_sets_->nmoltenres() << " positions." << std::endl;
	}
	if ( TR.visible() ) TR << " IG: " << ig_->getTotalMemoryUsage() << " bytes" << std::endl;

	//Need to check that starting residues are in the design shell (IG); if not, remove them
	utility::vector1< platform::uint > const molten_resvec(rotamer_sets_->moltenres_2_resid_vector());
	auto resvecit = start_res_vec_.begin();
	for ( ;  resvecit != start_res_vec_.end(); ) {
		bool in_design_shell = false;
		for ( core::Size mrvit : molten_resvec ) {
			platform::uint resvecit_uint(*resvecit);
			if ( resvecit_uint == mrvit ) {
				in_design_shell = true;
				break;
			}
		}
		if ( !( in_design_shell ) ) {
			if ( TR.visible() && verbose_ ) TR << " WARNING: residue " << *resvecit << " is not in the IG; I'm removing it from start_resnums .  If you wish to consider this residue, change your design shell task operations" << std::endl;
			//does NOT return iterator referencing same location after removal like std::vector does!
			start_res_vec_.erase( resvecit++ );
		} else {
			++resvecit;
		}
	}
	runtime_assert( !(start_res_vec_.empty()));

	this->search_IG_for_networks( pose ); // can be overriden for special cases of handling and searching IG

	if ( ! monte_carlo_branch_ ) {
		if ( TR.visible() ) TR << " INITIAL NUMBER OF H-BOND NETWORKS FOUND: " << network_vector_.size() << std::endl;
		if ( max_replicates_before_branch_ > 0 ) {
			remove_replicate_networks( max_replicates_before_branch_ );
			if ( TR.visible() ) {
				TR << " NUMBER OF NETWORKS AFTER INITIAL REMOVE_REP = " << network_vector_.size() << std::endl;
			}
		}
		if ( TR.visible() ) {
			TR << " MERGING NETWORKS TOGETHER TO FORM COMPLETE NETWORKS" << std::endl;
		}
		branch_overlapping_networks();
	} else {
		MC_branch_overlapping_networks();
	}
	if ( TR.visible() ) TR << " NUMBER OF H-BOND NETWORKS AFTER BRANCH: " << network_vector_.size() << std::endl;
	remove_replicate_networks( max_replicates_before_unsat_check_ );
	if ( TR.visible() ) TR << " NUMBER OF NETWORKS AFTER REMOVE REPLICATES = " << network_vector_.size() << std::endl;

	select_best_networks();
	if ( TR.visible() ) TR << " NUMBER OF NETWORKS AFTER SELECTION = " << network_vector_.size() << std::endl;
	remove_replicate_networks( max_rep_ );
	if ( TR.visible() ) TR << " NUMBER OF NETWORKS AFTER REMOVE_REP = " << network_vector_.size() << std::endl;

}//run

void
HBNet::apply( Pose & pose )
{
	if ( !( pose.pdb_info() ) ) {
		pose.pdb_info( pointer::make_shared< core::pose::PDBInfo >( pose, true ) );
	}
	// Delete all HBNet comments in pose
	//   If keep_existing or extend_existing, the input in networks will be re-detected and these comments regenerated
	std::map< std::string, std::string > comments = core::pose::get_all_comments(pose);
	for ( std::map< std::string, std::string >::const_iterator com = comments.begin(); com != comments.end(); ++com ) {
		if ( com->first.substr(0, 5 ) != "HBNet" ) {
			core::pose::delete_comment( pose, com->first );
		}
	}
	//need to reset in case apply called multiple times from code with different poses (e.g. nstruct > 1)
	network_vector_.clear();
	rotamer_sets_.reset(); // resets OP to 0
	ig_.reset(); // resets OP to 0
	monte_carlo_seeds_.clear();

	//shouldn't ever be NULL
	if ( scorefxn_ == nullptr ) {
		scorefxn_ = core::scoring::get_score_function(); // set to default
	}
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		set_symmetry(pose);
	}

	if ( !only_native_ && task_factory_ != nullptr ) {
		task::PackerTaskOP raw_task = create_ptask(pose);
		if ( raw_task->rotamer_links_exist() ) {
			if ( TR.visible() ) TR << " WARNING!!!! ROTAMER LINKS DETECTED: HBNet can't currently handle rotamer links!" << std::endl;
			rotamer_links_ = raw_task->rotamer_links();
		}
	}

	pose.update_residue_neighbors();
	(*init_scorefxn_)(pose); // safeguard to ensure input pose is scored with same sfxn used to propogate IG / HBondGraph

	orig_pose_ = pointer::make_shared< Pose >( pose ); //pose::PoseOP pointing to original pose
	// HBNet returns mutliple poses, generated on-the-fly by places networks onto copies of this pose

	if ( verbose_ && TR.visible() ) {
		TR << " ========================================================================" << std::endl;
		if ( only_native_ ) {
			TR << " WILL ONLY CONSIDER NATIVE INPUT ROTAMERS IN NETWORK SEARCH" << std::endl;
			TR << " NO DESIGN OR DESIGN OPTIONS WILL BE CONSIDERED" << std::endl;
		} else {
			TR << " hydrogen bond threshold, hb_threshold = " << hydrogen_bond_threshold_ << std::endl;
			TR << " clash_threshold = " << clash_threshold_ << std::endl;
			TR << " Residue types allowed for h-bond network design = " << des_residues_ << std::endl;
		}
		TR << " store_subnetworks = " << store_subnetworks_ << std::endl;
		TR << " ========================================================================" << std::endl;
	}

	// put all residues in selector in start_res_vec_ (we check later to ensure that these positions are actually in IG)
	//    ResidueSubset is a bool vector of all positions in the pose; can use Vikram's ReferencePose if need to add/delete residues
	if ( start_selector_ ) {
		start_res_vec_ = core::select::get_residue_set_from_subset( start_selector_->apply( pose ) );
	}
	if ( core_selector_ ) {
		core_residues_ = core_selector_->apply( pose );
	} else {
		core::select::residue_selector::LayerSelectorOP core_layer( pointer::make_shared< core::select::residue_selector::LayerSelector >() );
		core_layer->set_layers( true, false, false );
		//core_layer->set_use_sc_neighbors( true ); // now true by default
		core_layer->set_cutoffs( 4.4 /* core */, 2.0 /* surface */ ); //default core is 5.2, here we use 4.4 to be overly cautious with unsats
		runtime_assert( core_layer->use_sc_neighbors() );
		core_residues_ = core_layer->apply( pose );
	}
	if ( boundary_selector_ ) {
		boundary_residues_ = boundary_selector_->apply( pose );
	} else {
		core::select::residue_selector::LayerSelectorOP boundary_layer( pointer::make_shared< core::select::residue_selector::LayerSelector>() );
		boundary_layer->set_layers( false, true, false );
		//boundary_layer->set_use_sc_neighbors( true );
		boundary_layer->set_cutoffs( 4.0 /* core */, 2.0 /* surface */ ); //default core is 5.2, here we use 4 to be overly cautious with unsats
		boundary_residues_ = boundary_layer->apply( pose );
	}
	if ( verbose_ && TR.visible() ) {
		for ( Size r = 1; r <= pose.total_residue(); ++r ) {
			if ( res_is_core(r) ) {
				TR << "res " << r << " is core" << std::endl;
			}
		}
	}
	if ( keep_existing_networks_ || extend_existing_networks_ ) {
		//Need to check: are PDBInfoLabels kept after pose made symmetric? what about after BGS?
		core::select::residue_selector::ResiduePDBInfoHasLabelSelectorOP hbnet_info_label( pointer::make_shared< core::select::residue_selector::ResiduePDBInfoHasLabelSelector >( "HBNet" ) );
		input_hbnet_info_residues_ = hbnet_info_label->apply( pose );
	}
	// setup method common to all HBNet
	setup( pose );
	// run method common to all HBNet
	run( pose );
	TR.flush();
	this->prepare_output();
	if ( only_native_ ) {
		// in only_native_ case, HBNet reslabels and comments get added to cached orig_pose_, but apply() is passed Pose & pose REFERENCE, so need update "pose"
		pose = nonconst_get_orig_pose();
		if ( native_networks_.empty() ) set_last_move_status( moves::FAIL_RETRY );
		else set_last_move_status( moves::MS_SUCCESS );
		return;
	}
	TR.Debug << "After prepare_output, network_vector_.size() = " << network_vector_.size() << std::endl;
	if ( allow_no_hbnets_ && output_vector_.size() == 0 ) {
		if ( TR.visible() ) TR << "DID NOT FIND SOLUTIONS THAT MEET YOUR CRITERIA! BUT continuing because allow_no_hbnets is set" << std::endl;
		return;
	}
	if ( output_vector_.size() == 0 )  {
		if ( TR.visible() ) TR << "DID NOT FIND SOLUTIONS THAT MEET YOUR CRITERIA! EXITING" << std::endl;
		set_last_move_status( moves::FAIL_RETRY );
		ig_.reset(); // resets OP to 0
		//lkb_.reset();
		packer_neighbor_graph_.reset();
		return;
	} else {
		set_last_move_status( moves::MS_SUCCESS );
	}
	// need to reverse because no pop_front for std::vector -- could use list but already set up for sorting functions with vectors of HBondNetStructOP's
	std::reverse( output_vector_.begin(), output_vector_.end());

	std::set< Size > net_ids( output_vector_.back() );
	output_vector_.pop_back(); // pop_back() removes without returning

	Pose ala_copy = *ala_pose_; // deep copy
	if ( output_poly_ala_background_ ) pose = ala_copy; // "pose" is & ref

	std::string comment_str( print_headers() + this->print_additional_headers() );
	if ( net_ids.size() > 1 && TR.visible() ) {
		TR << std::endl << "combining networks ";
	}
	for ( const auto & net_id : net_ids ) {

		if ( net_ids.size() > 1 && TR.visible() ) {
			TR << net_id << "  ";
		}
		//TR << "net_id = " << *net_id << std::endl;
		runtime_assert( get_network_by_id(net_id) != nullptr );
		//place_rots_on_pose( pose, *(get_network_by_id( net_id )), (get_network_by_id( net_id ))->is_native, bridging_waters_, (get_network_by_id( net_id ))->waterrots.empty() );
		place_rots_on_pose( pose, *(get_network_by_id( net_id )), (get_network_by_id( net_id ))->is_native );
		//add_reslabels_to_pose( pose, *(get_network_by_id( net_id )) ); // needs to be done on the fly when adding waters
		comment_str = comment_str + "\n" + (get_network_by_id( net_id ))->outstring;
		core::scoring::constraints::ConstraintSetOP cst_op( pose.constraint_set()->clone() ); // get existing csts and add to them

		//        if ( bridging_waters_ ){
		//            //NEED TO SCORE BEFORE get_hbond_atom_pairs()
		//            ( *scorefxn_ )( pose ); // TODO NEED BETTER WAY THAN SCORING EACH TIME AND OVERRIDING hbond_vec
		//            get_hbond_atom_pairs( *(get_network_by_id( net_id )), pose ); // if residues (waters) have been added, need to update
		//        }
		set_constraints( pose, *cst_op, *get_network_by_id( net_id ), write_cst_files_ );
		pose.constraint_set(cst_op); // add constraints to the pose
	}
	if ( TR.visible() ) TR << std::endl;
	core::pose::add_comment( pose, "HBNet Design details: ", "\n" + comment_str + "\n" );

	//pose.update_residue_neighbors(); //now taken care of in place_rots_on_pose()
	( *scorefxn_ )( pose );

	ig_.reset(); // resets OP to 0 / nullptr
	packer_neighbor_graph_.reset();
	if ( TR.visible() ) TR.flush();
	if ( TR.Warning.visible() ) TR.Warning.flush();
	if ( TR.Debug.visible() ) TR.Debug.flush();
}// apply

std::string HBNet::get_name() const {
	return mover_name();
}

std::string HBNet::mover_name() {
	return "HBNet";
}

void
HBNet::attributes_for_hbnet( utility::tag::AttributeList & attlist )
{
	using namespace utility::tag;
	attlist
		+ XMLSchemaAttribute::attribute_w_default(
		"hb_threshold", xsct_real,
		"2-body h-bond energy cutoff to define rotamer pairs that h-bond. "
		"I've found that -0.5 without ex1-ex2 is the best starting point. "
		"If using ex1-ex2, try -0.75. This parameter is the most important and requires some tuning; "
		"the tradeoff is that the more stringent (more negative), the faster it runs but you miss a lot of networks; "
		"too positive and it will run forever; using ex1-ex2 results in many redundant networks that end up being filtered out anyway.",
		"-0.5" )
		+ XMLSchemaAttribute::attribute_w_default(
		"onebody_hb_threshold", xsct_real,
		"h-bond energy cutoff with their symmetric clones to define rotamer pairs that h-bond. ",
		"-0.4" )
		+ XMLSchemaAttribute::attribute_w_default(
		"charge_charge_rep_cutoff", xsct_real,
		"new feature, not fully tested, use at your own risk!;",
		"1.0" )
		+ XMLSchemaAttribute::attribute_w_default(
		"clash_threshold", xsct_real,
		"fa_rep energy above which two rotamers are considered clashing; used in all HBNet clash checks",
		"1.0" )
		+ XMLSchemaAttribute::attribute_w_default(
		"upper_score_limit", xsct_real,
		"upper score limit for network against background pose to weed out anything terrible",
		"15.0" )
		+ XMLSchemaAttribute::attribute_w_default(
		"find_native_networks", xsct_rosetta_bool,
		"Will find and report native networks in input pose but will also do design; "
		"for keep_existing_networks=true or extend_existing_networks=true, "
		"in addition to starting from \"HBNet\" PDBInfoLabel tags, "
		"HBnet will find all native networks in your input pose that meet your criteria (specified by options).",
		"false")
		+ XMLSchemaAttribute::attribute_w_default(
		"store_network_scores_in_pose", xsct_rosetta_bool,
		"Will store internal evaluation metrics in the pose so they will be included in score.sc files"
		"and accessible via ReadPoseExtraScoreFilter.",
		"false")
		+ XMLSchemaAttribute::attribute_w_default(
		"output_poly_ala_background", xsct_rosetta_bool,
		"for returned poses, rather than place network onto input pose, place into poly-ALA background; "
		" where all positions that are designable or packable (according to user-defined task operations)"
		" are poly-alanine",
		"false")
		+ XMLSchemaAttribute::attribute_w_default(
		"find_only_native_networks", xsct_rosetta_bool,
		"HBnet will find only find native networks in your input pose that meet your criteria "
		"(specified by options) and return a single pose with csts for those networks. "
		"If this option is true, all options below are overridden and HBNet does not search for new networks",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"keep_existing_networks", xsct_rosetta_bool,
		"In addition to design, Keeps existing networks from the input pose for each pose returned by HBNet, "
		"and turn on csts for all; existing networks are identified by default by \"HBNet\" PDBInfoLabel "
		"tags in the input pose, but can also be detected anew by setting find_native_networks=1. "
		"If keep_existing_networks=true and extend_existing_networks=false, then HBNet internally turns off "
		"design at input network positions (prevent_repacking); new networks are searched for and designed at "
		"the other positions based on your task_ops. If keep_existing_networks=true and extend_existing_networks=true, "
		"then HBNet internally puts only the input rotamer of each network residue into the IG to try to extend; "
		"an extended network will replace its native network if it is better (otherwise native networks retained).",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"extend_existing_networks", xsct_rosetta_bool,
		"Detects existing networks and tries to extend them, and also will search for new networks at other positions "
		"based on your criteria. Existing networks identified by HBNet PDBInfoLabel tags in the input pose, "
		"but can also be detected anew by setting find_native_networks=1. For existing networks, HBNet internally "
		"puts only the input rotamer of each network residue into the IG to try to extend; "
		"an extended network will replace its native network if it is better (otherwise native networks retained).",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"only_extend_existing", xsct_rosetta_bool,
		"Will not look for new networks at other positions; will only try to extend and improve the existing networks.",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"minimize", xsct_rosetta_bool,
		"new feature, not fully tested, use at your own risk!; "
		"will minimize network on background pose before scoring/evaluating",
		"true" )
		+ XMLSchemaAttribute::attribute_w_default(
		"store_subnetworks", xsct_rosetta_bool,
		"store subnetworks as separate networks that can be passed back; not recommended",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"secondary_search", xsct_rosetta_bool,
		"if during IG traversal, a search trajectory terminates in rotamer than cannot make any h-bonds, "
		"search again from that rotamer using a lower hb_threshold (-0.25)",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"secondary_threshold", xsct_real,
		"for secondary search option, what is the second (lower) hydrogen_bond_threshold to use",
		"-0.25" )
		+ XMLSchemaAttribute::attribute_w_default(
		"write_network_pdbs", xsct_rosetta_bool,
		"writes out pdb files of only the network (in poly-Ala background where any designable/packable "
		"residue is Ala -- rest of pose untouched); "
		"this is simply for debugging and visualizing the network as detected by HBNet",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"write_cst_files", xsct_rosetta_bool,
		"writes out Rosetta .cst constraint files for each network",
		"true" )
		+ XMLSchemaAttribute::attribute_w_default(
		"max_replicates", xsct_non_negative_integer,
		"max replicate networks that are allowed (defined as networks with same residue positions and aa types), "
		"but different rotamers; deafult is that replicates are not allowed",
		"1" )
		+ XMLSchemaAttribute::attribute_w_default(
		"max_replicates_before_branch", xsct_non_negative_integer,
		"hack that improves runtime without significantly affecting results; use at your own risk! "
		"if set to be greater than 0, is the max replicate networks that are allowed before branching (defined as networks with same residue positions and aa types)",
		"0" )
		+ XMLSchemaAttribute::attribute_w_default(
		"max_replicates_before_unsat_check", xsct_non_negative_integer,
		"hack that improves runtime without significantly affecting results; use at your own risk! "
		"max replicate networks that are allowed before unsat check is performed (defined as networks with same residue positions and aa types)",
		"1" )
		+ XMLSchemaAttribute::attribute_w_default(
		"tyr_hydroxyls_must_donate", xsct_rosetta_bool,
		"new feature, not fully tested, use at your own risk!; "
		"require that TYR must donate to be considered satisfied",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"hydroxyls_must_donate", xsct_rosetta_bool,
		"new feature, not fully tested, use at your own risk!; "
		"require that hydroxyls must donate their H to be considered satisfied",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"use_pdb_numbering", xsct_rosetta_bool,
		"in tracer output, use pdb numbering and chains when printing out network info (instead of Rosetta sequential numbering)",
		"true" )
		+ XMLSchemaAttribute::attribute_w_default(
		"no_heavy_unsats_allowed", xsct_rosetta_bool,
		"see max_unsat_Hpol for details",
		"1" )
		+ XMLSchemaAttribute::attribute_w_default(
		"show_task", xsct_rosetta_bool,
		"will show packer task by printing to TR, aka what's designable and packable and fixed based on your taskops",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"min_network_size", xsct_non_negative_integer,
		"minimum number of residues required for each network. Note: in symmetric cases, "
		"this refers to the number of residues in the ASU!",
		"3" )
		+ XMLSchemaAttribute::attribute_w_default(
		"max_network_size", xsct_non_negative_integer,
		"maximum number of residues required for each network",
		"15" )
		+ XMLSchemaAttribute::attribute_w_default(
		"min_unique_networks", xsct_non_negative_integer,
		"minimum number of networks with unique aa composition / rotamers",
		"1" )
		+ XMLSchemaAttribute::attribute_w_default(
		"min_core_res", xsct_non_negative_integer,
		"minimum core residues each network must have (as defined by core selector)",
		"0" )
		+ XMLSchemaAttribute::attribute_w_default(
		"min_boundary_res", xsct_non_negative_integer,
		"minimum boundary residues a network must contain",
		"0" )
		+ XMLSchemaAttribute::attribute_w_default(
		"max_unsat_Hpol", xsct_non_negative_integer,
		"maximum number of buried unsatisfied polar Hpol atoms allowed in each network. "
		"Note that if there are heavy atom donors or acceptors that are buried and unsatisfied, "
		"those networks are thrown out; only unsatisfied Hpols where the heavy "
		"atom donor is making at least one hydrogen bond. This behavior can be overridden "
		"to allow heavy atom unsats by setting no_heavy_unsats_allowed=false.",
		"5" )
		+ XMLSchemaAttribute::attribute_w_default(
		"min_percent_hbond_capacity", xsct_real,
		"minimum percent_hbond_capacity a network must have, where 1.0 means the maximum number of h-bonds are made, "
		"considering all polar atoms in the rotamers that comprise the network; unlike satisfaction, this metric is independent of burial; "
		" in rare cases, percent_hbond_capacity can be greater than 1, for example a hydroxyl could participate "
		"in 3 h-bonds even though we define max capacity as 2.",
		"0.5" )
		+ XMLSchemaAttribute::attribute_w_default(
		"verbose", xsct_rosetta_bool,
		"print out all HBNet tracer statements; only useful for debugging",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"design_residues", xs_string,
		"string of one-letter AA codes; which polar residues types do you want to include "
		"in the search; the default is all AA's that can potentially make h-bonds, "
		"further restricted by the task_operations you pass.",
		"STRKHYWNQDE" )
		+ XMLSchemaAttribute::attribute_w_default(
		"use_only_input_rot_for_start_res", xsct_rosetta_bool,
		"keep the input rotamer for starting residues fixed; only sample proton chis during network search",
		"false")
		+ XMLSchemaAttribute(
		"start_resnums", xsct_residue_number_cslist,
		"comma delimited list of residue positions to start network search from "
		"(e.g. \"1,2,3,5\"); now is better to use start_selector residue selector")
		+ XMLSchemaAttribute(
		"start_selector", xs_string,
		"residue selector that tells HBNet which residues to start from (will only search for networks that include these resid")
		+ XMLSchemaAttribute(
		"core_selector", xs_string,
		"residue selector that defines what HBNet considers \"core\"; "
		"used in buriedness determination for unsats; "
		"default is layer selector default using sidechain neighbors(core=5.2)." )
		+ XMLSchemaAttribute(
		"boundary_selector", xs_string,
		"residue selector to define what HBNet considers boundary during comparisons/scoring of networks" )
		+ XMLSchemaAttribute(
		"allow_no_hbnets", xsct_rosetta_bool,
		"allows the program to continue when no hbnets are found" )
		//Monte carlo options:
		+ XMLSchemaAttribute::attribute_w_default(
		"monte_carlo", xsct_rosetta_bool,
		"Step right up and try your luck with this stochastic HBNet protocol",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"monte_carlo_branch", xsct_rosetta_bool,
		"Outadated equivalent to monte_carlo",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"monte_carlo_seed_must_be_buried", xsct_rosetta_bool,
		"(If monte_carlo_branch) only branch from hbonds where at least one residue is buried. Effectively, this results in only finding networks that have at least one buried residue.",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"monte_carlo_seed_must_be_fully_buried", xsct_rosetta_bool,
		"(If monte_carlo_branch) only branch from hbonds where both residues are buried. This results in only finding networks that have at least one buried hbond but this does not prevent having additional exposed hbonds.",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"max_mc_nets", xsct_non_negative_integer,
		"(if monte_carlo_branch) This is experimental and so far it does not look super useful. Maximum number of networks that the monte carlo protocol will store. Loose rule of thumb is to make this 10x the max number of nets you want to end up with. The reason you do not want this to be too large is that each of these networks goes through a relatively expensive quality check. A value of zero (which is default) represents having no limit.",
		"0" )
		+ XMLSchemaAttribute::attribute_w_default(
		"total_num_mc_runs", xsct_non_negative_integer,
		"(if monte_carlo_branch) number of monte carlo runs to be divided over all the seed hbonds. A single monte carlo run appears to take roughly 1 ms (very loose estimate).",
		"10000" )
		+ XMLSchemaAttribute::attribute_w_default(
		"seed_hbond_threshold", xsct_real,
		"(if monte_carlo_branch) Maybe you only want to branch from strong hbonds. If this value is -1.2, for example, then only hbonds with a strength of -1.2 or lower will be branched from.",
		"0" );

	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	rosetta_scripts::attributes_for_parse_score_function( attlist );

}
void HBNet::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attributes_for_hbnet( attlist );

	moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"HBNet is a method to explicitly detect and design hydrogen bond networks within Rosetta. "
		"It functions as a mover within the RosettaScripts framework and will exhaustively search "
		"for all networks within the design space that you define with TaskOperations",
		attlist );
}

void HBNet::initialize_hbond_graph(){
	debug_assert( rotamer_sets_ );
	core::Size const nrot = rotamer_sets_->nrotamers();
	hbond_graph_ = pointer::make_shared< AtomLevelHBondGraph >( nrot );
	core::pack::hbonds::init_node_info( * hbond_graph_, * rotamer_sets_ );
}

bool HBNet::edge_can_yield_monte_carlo_seed( core::Size resid1, core::Size resid2 ) const{

	if ( monte_carlo_seed_must_be_fully_buried_ ) {
		if ( ! ( core_residues_[ resid1 ] && core_residues_[ resid2 ] ) ) {
			return false;
		}
	}

	if ( monte_carlo_seed_must_be_buried_ ) {
		if ( !core_residues_[ resid1 ] && !core_residues_[ resid2 ] ) {
			return false;
		}
	}

	if ( start_res_vec_.size() > 0 ) {
		if ( std::find( start_res_vec_.begin(), start_res_vec_.end(), resid1 ) == start_res_vec_.end() &&
				std::find( start_res_vec_.begin(), start_res_vec_.end(), resid2 ) == start_res_vec_.end() ) {
			return false;
		}
	}

	return true;
}

bool HBNet::edge_can_yield_monte_carlo_seed( core::pack::interaction_graph::EdgeBase const * edge ) const{
	core::Size const resid1 = rotamer_sets_->moltenres_2_resid( edge->get_first_node_ind() );
	core::Size const resid2 = rotamer_sets_->moltenres_2_resid( edge->get_second_node_ind() );
	return edge_can_yield_monte_carlo_seed( resid1, resid2 );
}

AtomLevelHBondNode * HBNet::get_next_node( NetworkState & current_state ){

	std::list< std::pair< core::Size, AtomLevelHBondNode * > > num_nbrs_per_node;
	core::Size total_num_nbrs = 0;

	bool all_nodes_satisfied( true );
	for ( auto const & size_atominfo_pair : current_state.unsatisfied_sc_atoms_const() ) {
		if ( size_atominfo_pair.second.size() ) {
			if ( ! size_atominfo_pair.second[ 1 ].is_hydrogen() ) { //TODO?
				all_nodes_satisfied = false;
				break;
			}
		}
	}


	for ( AtomLevelHBondNode const * const current_node : current_state.nodes() ) {

		if ( ! all_nodes_satisfied ) {
			unsigned int const mres = current_node->moltenres();
			if ( ! current_state.mres_has_heavy_unsats( mres ) ) continue;
		}
		core::Size const current_node_ind = current_node->global_rotamer_id();

		for ( utility::graph::EdgeListConstIterator it = current_node->const_edge_list_begin(); it != current_node->const_edge_list_end(); ++it ) {
			AtomLevelHBondNode * const proposed_new_node = hbond_graph_->get_hbondnode( (*it)->get_other_ind( current_node_ind ) );
			if ( ! node_is_compatible( current_state, proposed_new_node ) ) {
				continue;
			}

			if ( ! all_nodes_satisfied ) { // check if this edge satisfies a current unsat atom in the network
				auto const * const al_edge = static_cast< AtomLevelHBondEdge const * >( *it );
				if ( ! edge_satisfies_heavy_unsat_for_node( current_state, current_node, al_edge ) ) continue;
			}
			// if everything is satisfied, we'll still try to grow (and network state will have already been stored)

			core::Size num_nbrs = 1;//add dummy of 1
			core::Size const proposed_new_node_ind = proposed_new_node->global_rotamer_id();
			for ( utility::graph::EdgeListConstIterator it2 = proposed_new_node->const_edge_list_begin(); it2 != proposed_new_node->const_edge_list_end(); ++it2 ) {
				AtomLevelHBondNode const * const proposed_nbr_node = hbond_graph_->get_hbondnode( ( *it2 )->get_other_ind( proposed_new_node_ind ) );
				if ( node_is_compatible( current_state, proposed_nbr_node ) ) {
					++num_nbrs;
				}
			}

			total_num_nbrs += num_nbrs;
			num_nbrs_per_node.emplace_back( num_nbrs, proposed_new_node );
		}//it

	}//current_node

	if ( total_num_nbrs == 0 ) {
		return nullptr;
	}

	core::Real const rg_uniform = numeric::random::rg().uniform();
	core::Size const random_nbr_count = 1 + rg_uniform * total_num_nbrs;
#ifndef NDEBUG
	if ( random_nbr_count > total_num_nbrs ) {
		utility_exit_with_message( "Error: random_nbr_count > total_num_nbrs" );
		return num_nbrs_per_node.back().second;
	}
#endif

	core::Size nbr_counter = 0;
	for ( std::pair< core::Size, AtomLevelHBondNode * > & candidate : num_nbrs_per_node ) {
		nbr_counter += candidate.first;
		if ( nbr_counter >= random_nbr_count ) {
			return candidate.second;
		}
	}

	utility_exit_with_message( "Dead code has been reached" );
	return nullptr;
}


bool HBNet::monte_carlo_seed_is_dead_end( AtomLevelHBondEdge const * monte_carlo_seed ){
	AtomLevelHBondNode const * const nodeA = hbond_graph_->get_hbondnode( monte_carlo_seed->get_first_node_ind()  );
	AtomLevelHBondNode const * const nodeB = hbond_graph_->get_hbondnode( monte_carlo_seed->get_second_node_ind() );

	core::Size const nodeA_ind = nodeA->global_rotamer_id();
	core::Size const nodeB_ind = nodeB->global_rotamer_id();

	if ( nodeA_ind == nodeB_ind ) return false;//Symmetry

	for ( utility::graph::EdgeListConstIterator it_A = nodeA->const_edge_list_begin();
			it_A != nodeA->const_edge_list_end(); ++it_A ) {
		AtomLevelHBondNode const * const temp_node = hbond_graph_->get_hbondnode( (*it_A)->get_other_ind( nodeA_ind ) );

		if ( temp_node->moltenres() == nodeB->moltenres() ) continue;
		if ( nodeB->clashes( temp_node->global_rotamer_id() ) ) continue;

		return false;
	}

	for ( utility::graph::EdgeListConstIterator it_B = nodeB->const_edge_list_begin();
			it_B != nodeB->const_edge_list_end(); ++it_B ) {
		AtomLevelHBondNode const * const temp_node = hbond_graph_->get_hbondnode( (*it_B)->get_other_ind( nodeB_ind ) );

		if ( temp_node->moltenres() == nodeA->moltenres() ) continue;
		if ( nodeA->clashes( temp_node->global_rotamer_id() ) ) continue;

		return false;
	}

	return true;
}

void HBNet::append_to_network_vector( utility::vector1< DecoratedNetworkState > const & designed_networks ){

	//using max_replicates_before_unsat_check_
	std::string const poly_A_sequence( task_->num_to_be_packed(), 'A' );
	std::map< std::string, unsigned int > num_replicates_for_sequence;

	std::vector< HBondNetStructOP > & network_vector = HBNet::get_net_vec();
	core::Size net_count = 0;

	for ( DecoratedNetworkState const & dec_state : designed_networks ) {
		NetworkState const & state = * dec_state.network_state;

		utility::vector1< HBondResStructCOP > HBondResStructs;
		HBondResStructs.reserve( state.nodes().size() );

		std::string sequence = poly_A_sequence;

		for ( AtomLevelHBondNode const * const hbond_node : state.nodes() ) {
			core::Size const mres = hbond_node->moltenres();
			core::Size const resid = rotamer_sets_->moltenres_2_resid( mres );
			core::Size const rot =  hbond_node->local_rotamer_id();
			auto res_ptr = rotamer_sets_->rotamer_set_for_moltenresidue( mres )->rotamer( rot );

			sequence[ mres - 1 ] = res_ptr->name1();

			HBondResStructCOP const new_HBondResStruct(
				pointer::make_shared< HBondResStruct >(
				resid,
				rot,
				res_ptr->name1(),
				chain_for_moltenres( mres ),
				res_ptr->is_protein(), //is_protein
				false, //is_solvent
				res_ptr->is_ligand()  //is_ligand
				)
			);
			HBondResStructs.push_back( new_HBondResStruct );
		}//for network_res


		if ( num_replicates_for_sequence.find( sequence ) == num_replicates_for_sequence.end() ) {
			num_replicates_for_sequence[ sequence ] = 1;
		} else {
			if ( ++num_replicates_for_sequence[ sequence ] > max_replicates_before_unsat_check_ ) {
				//skip this duplicate
				continue;
			}
		}

		HBondNetStructOP network_struct = pointer::make_shared< HBondNetStruct >();
		network_struct->is_native = false;
		network_struct->term_w_start = false;

		network_struct->term_w_cycle = state.edges().size() >= state.nodes().size();
		network_struct->residues = HBondResStructs;

		network_struct->lig_state_list.clear();
		network_struct->total_hbonds = state.edges().size();//only counts 1 hbond per edge
		network_struct->score = state.score();

		network_vector.push_back( network_struct );

		if ( ++net_count == max_mc_nets_ ) break;
	}

}

void HBNet::add_residue_to_network_state( NetworkState & current_network_state, AtomLevelHBondNode * node_being_added ) const {

	current_network_state.add_polar_atoms( node_being_added );

	for ( AtomLevelHBondNode const * const existing_node : current_network_state.nodes() ) {

		//check for hbond
		AtomLevelHBondEdge const * const hbond_edge =
			static_cast< AtomLevelHBondEdge * > ( node_being_added->find_edge( existing_node->get_node_index() ) );
		if ( hbond_edge ) {
			current_network_state.edges().push_back( hbond_edge );
			current_network_state.append_twobody_energy( hbond_edge->energy() );

			//Declare hbonding atoms to be satisfied!
			for ( HBondInfo const & hbond : hbond_edge->hbonds() ) {
				bool const first_node_is_donor = hbond.first_node_is_donor();
				unsigned short int const local_atom_id_A = hbond.local_atom_id_A();//Acceptor
				unsigned short int const local_atom_id_D = hbond.local_atom_id_D();//Donor
				unsigned short int const local_atom_id_H = hbond.local_atom_id_H();//Hydrogen

				unsigned int const existing_node_mres = existing_node->moltenres();
				unsigned int const new_node_mres = node_being_added->moltenres();

				bool const new_node_is_first_node = new_node_mres < existing_node_mres;
				if ( new_node_is_first_node == first_node_is_donor ) {
					//node_being_added is donor
					AtomLevelHBondNode::remove_atom_info_from_vec_stable( current_network_state.get_unsats_for_mres( existing_node_mres )->second, local_atom_id_A );

					utility::vector1< AtomInfo > & donor_atom_vec = current_network_state.get_unsats_for_mres( new_node_mres )->second;
					AtomLevelHBondNode::remove_atom_info_from_vec_stable( donor_atom_vec, local_atom_id_D );
					AtomLevelHBondNode::remove_atom_info_from_vec_stable( donor_atom_vec, local_atom_id_H );
				} else {
					//node_being_added is acc
					AtomLevelHBondNode::remove_atom_info_from_vec_stable( current_network_state.get_unsats_for_mres( new_node_mres )->second, local_atom_id_A );

					utility::vector1< AtomInfo > & donor_atom_vec = current_network_state.get_unsats_for_mres( existing_node_mres )->second;
					AtomLevelHBondNode::remove_atom_info_from_vec_stable( donor_atom_vec, local_atom_id_D );
					AtomLevelHBondNode::remove_atom_info_from_vec_stable( donor_atom_vec, local_atom_id_H );
				}

				//unsatisfied_sc_atoms
			}//hbond
		}//hbond_edge
	}//existing_node

	if ( symmetric_ ) { //TODO is symmtery the only case where this can happen?
		AtomLevelHBondEdge const * const hbond_edge =
			static_cast< AtomLevelHBondEdge * > ( node_being_added->find_edge( node_being_added->get_node_index() ) );
		if ( hbond_edge ) {
			current_network_state.edges().push_back( hbond_edge );
			current_network_state.append_twobody_energy( hbond_edge->energy() );

			//Declare hbonding atoms to be satisfied!
			unsigned int const new_node_mres = node_being_added->moltenres();

			for ( HBondInfo const & hbond : hbond_edge->hbonds() ) {
				AtomLevelHBondNode::remove_atom_info_from_vec_stable(
					current_network_state.get_unsats_for_mres( new_node_mres )->second,
					hbond.local_atom_id_A()
				);

				utility::vector1< AtomInfo > & donor_atom_vec = current_network_state.get_unsats_for_mres( new_node_mres )->second;
				AtomLevelHBondNode::remove_atom_info_from_vec_stable( donor_atom_vec, hbond.local_atom_id_D() );
				AtomLevelHBondNode::remove_atom_info_from_vec_stable( donor_atom_vec, hbond.local_atom_id_H() );
			}//hbond
		}//hbond_edge
	}//symm

	current_network_state.nodes().push_back( node_being_added );
	current_network_state.score( current_network_state.full_twobody_energy() / current_network_state.nodes().size() );
}

bool HBNet::network_state_is_satisfied( NetworkState & current_state ) const {

	for ( auto const & uint_vec_pair : current_state.unsatisfied_sc_atoms_const() ) {
		if ( uint_vec_pair.second.size() != 0 ) {
			if ( uint_vec_pair.second[ 1 ].is_hydrogen() ) {
				//has unsat hydrogens, no heavy unsats
				//TODO what do we do for hydrogens??? keep a running total?
			} else {
				//at least one heavy unsat
				return false;
			}
		}
	}
	return true;
}

void HBNet::score_network_state( NetworkState & state ) const {
	state.score( state.full_twobody_energy() / state.nodes().size() );
}

bool HBNet::register_new_network(
	NetworkState & current_network_state,
	utility::vector0< std::list< utility::vector1< core::Size > > > & hash_table,
	std::list< NetworkState > & designed_networks
) const {

	utility::vector1< core::Size > rotamer_ids;
	rotamer_ids.reserve( current_network_state.nodes().size() );

	core::Size sum = 0;
	for ( AtomLevelHBondNode const * const hbond_node : current_network_state.nodes() ) {
		core::Size const rotamer_id = hbond_node->global_rotamer_id();
		rotamer_ids.push_back( rotamer_id );
		sum += rotamer_id;
	}
	std::sort( rotamer_ids.begin(), rotamer_ids.end() );

	core::Size const hash = sum % hash_table.size();//hash_table is 0-indexed to make this easier
	for ( utility::vector1< core::Size > const & rot_vec : hash_table[ hash ] ) {
		if ( rot_vec.size() == rotamer_ids.size() ) {
			bool same = true;
			for ( core::Size ii=1; ii<=rot_vec.size(); ++ii ) {
				if ( rot_vec[ ii ] != rotamer_ids[ ii ] ) {
					same = false;
					break;
				}
			}
			if ( same ) {
				return false;
			}
		}
	}

	score_network_state( current_network_state );
	hash_table[ hash ].push_back( rotamer_ids );
	designed_networks.push_back( current_network_state );

	return true;
}

bool HBNet::node_is_compatible( NetworkState const & current_state, AtomLevelHBondNode const * node_being_added ){
	for ( AtomLevelHBondNode const * const existing_node : current_state.nodes() ) {
		if ( existing_node->moltenres() == node_being_added->moltenres() ) return false;
		if ( node_being_added->clashes( existing_node->global_rotamer_id() ) ) return false;
	}
	return true;
}

bool HBNet::network_state_is_done_growing( NetworkState const & current_state, AtomLevelHBondGraphCOP hbond_graph ){

	for ( AtomLevelHBondNode const * const current_node : current_state.nodes() ) {
		core::Size const current_node_ind = current_node->global_rotamer_id();

		for ( utility::graph::EdgeListConstIterator it = current_node->const_edge_list_begin(); it != current_node->const_edge_list_end(); ++it ) {
			core::Size const other_ind = (*it)->get_other_ind( current_node_ind );
			AtomLevelHBondNode const * const proposed_new_node = hbond_graph->get_hbondnode( other_ind );
			if ( node_is_compatible( current_state, proposed_new_node ) ) {
				return false;
			}
		}
	}

	return true;
}

core::Real HBNet::estimate_saturation( NetworkState const & network_state ) const {
	using namespace core::scoring::hbonds::graph;
	unsigned short int total_polar_groups_that_could_hbond( 0 );

	for ( AtomLevelHBondNode const * node : network_state.nodes() ) {
		core::conformation::Residue const & rotamer = * rotamer_sets_->rotamer( node->global_rotamer_id() );

		unsigned char const n_Hpol = rotamer.Hpol_index().size();
		total_polar_groups_that_could_hbond += n_Hpol;

		for ( unsigned char atomno = rotamer.first_sidechain_atom(),
				end = rotamer.nheavyatoms();
				atomno <= end; ++atomno ) {
			auto const & atom_type = rotamer.atom_type( atomno );
			if ( ! atom_type.is_donor() && ! atom_type.is_acceptor() ) continue;

			if ( rotamer.atom_type( atomno ).name() == "OH" ) {
				++total_polar_groups_that_could_hbond;
			} else {
				std::string const element = rotamer.atom_type( atomno ).element();
				unsigned char num_lone_pairs = 1;
				if ( element == "O" ) {
					num_lone_pairs = 2;
				} else if ( element == "N" ) {
					num_lone_pairs = 1;
				} else {
					TR << "MC HBNet is unprepared to handle this polar atom's element: " << element << std::endl;
				}
				//Let's not worry about charged atoms
				//auto const charge = rotamer.atomic_charge( atomno );
				//num_lone_pairs -= charge;
				total_polar_groups_that_could_hbond += num_lone_pairs;
			}
		}
	}

	//Is this fair - at least for a quick estimate?
	//Real const polar_groups_making_hbonds = 2 * network_state.edges().size();
	Real polar_groups_making_hbonds( 0 );
	for ( AtomLevelHBondEdge const * edge : network_state.edges() ) {
		polar_groups_making_hbonds += 2 * edge->hbonds().size();
	}

	return polar_groups_making_hbonds / total_polar_groups_that_could_hbond;
}

std::string HBNetCreator::keyname() const {
	return HBNet::mover_name();
}

moves::MoverOP
HBNetCreator::create_mover() const {
	return pointer::make_shared< HBNet >();
}

void HBNetCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HBNet::provide_xml_schema( xsd );
}

} //hbnet
} //protocols
