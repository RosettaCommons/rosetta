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
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/numbers.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Core headers
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/selection.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/io/Remarks.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSet_.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/InteractionGraphFactory.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh>
#include <core/pack/interaction_graph/PDInteractionGraph.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/chemical/AtomType.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/facts/FACTSEnergy.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
//#include <core/graph/Graph.hh>
//#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyGraph.hh>
//#include <core/scoring/methods/BridgingWaterEnergy.hh>
//#include <core/scoring/BridgingWaterPotential.hh>
//#include <core/scoring/methods/BridgingWaterInfo.hh>
#include <core/id/AtomID.hh>
#include <core/select/residue_selector/ResiduePDBInfoHasLabelSelector.hh>
#include <core/select/residue_selector/ResidueNameSelector.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/LayerSelector.hh>

// Protocols headers
#include <protocols/enzdes/AddorRemoveCsts.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
//#include <protocols/rigid/RB_geometry.hh>
//#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <protocols/scoring/Interface.hh>
//#include <protocols/toolbox/task_operations/PreventResiduesFromRepackingOperation.hh>
#include <protocols/simple_moves/MakePolyXMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>


using namespace core;
using namespace pose;
using namespace pack;
using namespace rotamer_set;
using namespace interaction_graph;
using namespace conformation;
using namespace scoring::hbonds;

namespace protocols {
namespace hbnet {

//static basic::Tracer TR( "protocols.hbnet.HBNet" );
static THREAD_LOCAL basic::Tracer TR( "protocols.hbnet.HBNet" );

HBNet::HBNet( ) :
	protocols::moves::Mover( "HBNet" ),
	//use_enzdes_cst_(0),
	benchmark_(0),
	write_network_pdbs_(0),
	write_cst_files_(1),
	find_native_(0),
    only_native_(0),
    keep_existing_networks_(0),
    extend_existing_networks_(0),
    only_extend_existing_(0),
	verbose_(0),
	symmetric_(0),
	multi_component_(0),
	show_task_(0),
	minimize_(1),
	//bridging_waters_(0),
	dump_resfile_(0),
	start_from_csts_(0),
	tyr_hydroxyls_must_donate_(0),
    hydroxyls_must_donate_(0),
	use_pdb_numbering_(1),
	no_heavy_unsats_allowed_(1),
	min_network_size_(3),
	max_network_size_(15),
	min_unique_networks_(1),
	min_core_res_(0),
	min_boundary_res_(0),
	max_unsat_(5),
	max_lig_unsat_(15),
	max_rep_(1),
	des_residues_("STRKHYWNQDE"),
	constraint_resfile_property_("NATRO"),
	start_res_vec_( /* NULL */ ),
	network_vector_( /* NULL */ ),
	native_networks_( /* NULL */ ),
	merged_vecs_( /* NULL */ ),
	//output_net_vec_(/* NULL */ ),
    output_vector_(/* NULL */ ),
	pore_radius_(2.5),
	atom_burial_cutoff_(0.5),
	hydrogen_bond_threshold_(-0.5),
	onebody_hb_threshold_(-0.4),
	charge_charge_rep_cutoff_(1.0),
	//bw_bb_sc_cutoff_(-1.0),
	//bw_cutoff_(-0.5),
	clash_threshold_(1.0),
	task_factory_( /* NULL */ ),
	init_scorefxn_(core::scoring::ScoreFunctionFactory::create_score_function( "HBNet" )),
	scorefxn_( /* NULL */ ),
	rotamer_sets_( /* NULL */ ),
	ig_( /* NULL */ ),
	//bw_ig_(0),
	rotamer_links_(0),
	store_subnetworks_(0), // should decide on default behavior and get rid of this
	secondary_search_(0),
	secondary_threshold_(-0.25),
	upweight_twobody_(1.0),
	start_selector_( /* NULL */ ),
	core_selector_( /* NULL */ ),
	boundary_selector_( /* NULL */ ),
	core_residues_(0),
	boundary_residues_(0),
    hbnet_info_residues_(0)
{
    //need to be able to use beta if called from code or another mover
    if ( basic::options::option[ basic::options::OptionKeys::corrections::beta ].value(true) ) {
        init_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "HBNet_beta" );
    }
}

HBNet::HBNet( std::string const name ) :
	protocols::moves::Mover( name ),
	//use_enzdes_cst_(0),
	benchmark_(0),
	write_network_pdbs_(0),
	write_cst_files_(1),
	find_native_(0),
    only_native_(0),
    keep_existing_networks_(0),
    extend_existing_networks_(0),
    only_extend_existing_(0),
	verbose_(0),
	symmetric_(0),
	multi_component_(0),
	show_task_(0),
	minimize_(1),
	//bridging_waters_(0),
	dump_resfile_(0),
	start_from_csts_(0),
	tyr_hydroxyls_must_donate_(0),
    hydroxyls_must_donate_(0),
	use_pdb_numbering_(1),
	no_heavy_unsats_allowed_(1),
	min_network_size_(3),
	max_network_size_(15),
	min_unique_networks_(1),
	min_core_res_(0),
	min_boundary_res_(0),
	max_unsat_(5),
	max_lig_unsat_(15),
	max_rep_(1),
	des_residues_("STRKHYWNQDE"),
	constraint_resfile_property_("NATRO"),
    start_res_vec_( /* NULL */ ),
    network_vector_( /* NULL */ ),
    native_networks_( /* NULL */ ),
    merged_vecs_( /* NULL */ ),
    //output_net_vec_(/* NULL */ ),
    output_vector_(/* NULL */ ),
    pore_radius_(2.5),
    atom_burial_cutoff_(0.5),
    hydrogen_bond_threshold_(-0.5),
    onebody_hb_threshold_(-0.4),
    charge_charge_rep_cutoff_(1.0),
    //bw_bb_sc_cutoff_(-1.0),
    //bw_cutoff_(-0.5),
    clash_threshold_(1.0),
    task_factory_( /* NULL */ ),
    init_scorefxn_(core::scoring::ScoreFunctionFactory::create_score_function( "HBNet" )),
    scorefxn_( /* NULL */ ),
    rotamer_sets_( /* NULL */ ),
    ig_( /* NULL */ ),
    //bw_ig_(0),
    rotamer_links_(0),
    store_subnetworks_(0), // should decide on default behavior and get rid of this
    secondary_search_(0),
    secondary_threshold_(-0.25),
    upweight_twobody_(1.0),
    start_selector_( /* NULL */ ),
    core_selector_( /* NULL */ ),
    boundary_selector_( /* NULL */ ),
    core_residues_(0),
    boundary_residues_(0),
    hbnet_info_residues_(0)
{
    //need to be able to use beta if called from code or another mover
    if ( basic::options::option[ basic::options::OptionKeys::corrections::beta ].value(true) ) {
        init_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "HBNet_beta" );
    }
}

//constructor to be called from code, NEED TO CLEAN THIS UP!
HBNet::HBNet( core::scoring::ScoreFunctionCOP scorefxn,
	Size max_unsat,
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
	//bool minimize /*true*/
) :
	protocols::moves::Mover( "HBNet" ),
	//use_enzdes_cst_(0),
	benchmark_(0),
	write_network_pdbs_(0),
	write_cst_files_(1),
	find_native_(find_native),
    only_native_(only_native),
    keep_existing_networks_(keep_existing),
    extend_existing_networks_(extend_existing),
    only_extend_existing_(only_extend),
	verbose_(0),
	symmetric_(0),
	multi_component_(0),
	show_task_(0),
	minimize_(1),
	//bridging_waters_(0),
	dump_resfile_(0),
	start_from_csts_(0),
	tyr_hydroxyls_must_donate_(0),
    hydroxyls_must_donate_(0),
	use_pdb_numbering_(1),
	no_heavy_unsats_allowed_(1),
	min_network_size_(min_network_size),
	max_network_size_(max_network_size),
	min_unique_networks_(1),
	min_core_res_(0),
	min_boundary_res_(0),
	max_unsat_(max_unsat),
	max_lig_unsat_(15),
	max_rep_(1),
	des_residues_(des_residues),
	constraint_resfile_property_("NATRO"),
    start_res_vec_( /* NULL */ ),
    network_vector_( /* NULL */ ),
    native_networks_( /* NULL */ ),
    merged_vecs_( /* NULL */ ),
    //output_net_vec_(/* NULL */ ),
    output_vector_(/* NULL */ ),
    pore_radius_(2.5),
    atom_burial_cutoff_(0.5),
    hydrogen_bond_threshold_(hb_threshold),
    onebody_hb_threshold_(-0.4),
    charge_charge_rep_cutoff_(1.0),
    //bw_bb_sc_cutoff_(-1.0),
    //bw_cutoff_(-0.5),
    clash_threshold_(1.0),
    task_factory_( /* NULL */ ),
    init_scorefxn_(core::scoring::ScoreFunctionFactory::create_score_function( "HBNet" )),
    scorefxn_(scorefxn->clone()),
    rotamer_sets_( /* NULL */ ),
    ig_( /* NULL */ ),
    //bw_ig_(0),
    rotamer_links_(0),
    store_subnetworks_(0), // should decide on default behavior and get rid of this
    secondary_search_(0),
    secondary_threshold_(-0.25),
    upweight_twobody_(1.0),
    start_selector_( /* NULL */ ),
    core_selector_( /* NULL */ ),
    boundary_selector_( /* NULL */ ),
    core_residues_(0),
    boundary_residues_(0),
    hbnet_info_residues_(0)
{
    //need to be able to use beta if called from code or another mover
    if ( basic::options::option[ basic::options::OptionKeys::corrections::beta ].value(true) ) {
        init_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "HBNet_beta" );
    }
    if ( only_native_ ) find_native_ = true;
    if ( only_extend_existing_ ) extend_existing_networks_ = true;
}

//// TODO NEED: Copy constructor only copy over key settings and parameters;


std::string
HBNetCreator::keyname() const {
	return HBNetCreator::mover_name();
}

protocols::moves::MoverOP
HBNetCreator::create_mover() const {
	return protocols::moves::MoverOP( new HBNet );
}

std::string
HBNetCreator::mover_name() {
	return "HBNet";
}

protocols::moves::MoverOP
HBNet::clone() const {

	return( protocols::moves::MoverOP( new HBNet( *this ) ) );
}

protocols::moves::MoverOP
HBNet::fresh_instance() const {

	return protocols::moves::MoverOP( new HBNet );
}

//destructor
HBNet::~HBNet(){}

core::pack::task::TaskFactoryOP
HBNet::task_factory() const
{
	return task_factory_;
}

void
HBNet::task_factory( core::pack::task::TaskFactoryOP task_factory )
{
	task_factory_ = task_factory;
}

void
HBNet::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose )
{
	hydrogen_bond_threshold_ = tag->getOption<core::Real>( "hb_threshold" ,-0.5);
	//bw_cutoff_ = tag->getOption<core::PackerEnergy>("bw_cutoff",-0.5);
	onebody_hb_threshold_ = tag->getOption<core::PackerEnergy>("onebody_hb_threshold",-0.4);
	charge_charge_rep_cutoff_ = tag->getOption<core::PackerEnergy>("charge_charge_rep_cutoff",1.0);
	clash_threshold_ = tag->getOption<core::PackerEnergy>("clash_threshold",1.0);
    find_native_ = tag->getOption<bool>( "find_native_networks", 0 );
    only_native_ = tag->getOption<bool>( "find_only_native_networks", 0 );
    keep_existing_networks_ = tag->getOption<bool>( "keep_existing_networks", 0 );
    extend_existing_networks_ = tag->getOption<bool>( "extend_existing_networks", 0 );
    only_extend_existing_ = tag->getOption<bool>( "only_extend_existing", 0 );
	minimize_ = tag->getOption< bool > ("minimize",1);
	store_subnetworks_ = tag->getOption< bool >("store_subnetworks",0);
	secondary_search_ = tag->getOption< bool >("secondary_search",0);
	secondary_threshold_ = tag->getOption< Real >("secondary_threshold",-0.25);
	write_network_pdbs_ = tag->getOption< bool >( "write_network_pdbs", 0 );
	write_cst_files_ = tag->getOption< bool >( "write_cst_files", 1 );
	max_rep_ = tag->getOption< Size >( "max_replicates", 1 );
	constraint_resfile_property_ = tag->getOption< bool >( "constraint_resfile_property", "NATRO");
	pore_radius_ = tag->getOption< Real >( "pore_radius", 2.5 );
	atom_burial_cutoff_ = tag->getOption< Real >( "atom_burial_cutoff", 0.5);
	//bridging_waters_ = tag->getOption<bool>( "bridging_waters", 0 );
	//bridging_water_rec_limit_= tag->getOption<Size>( "bridging_water_rec_limit", 2 );
	//bw_binary_cut_ = tag->getOption<Real>( "bw_binary_cut", -1.0 );
	tyr_hydroxyls_must_donate_ = tag->getOption<bool>( "tyr_hydroxyls_must_donate", 0 ); // only for unsat check -- should move to filter
    hydroxyls_must_donate_ = tag->getOption<bool>( "hydroxyls_must_donate", 0 ); // only for unsat check -- should move to filter
	use_pdb_numbering_ = tag->getOption<bool>( "use_pdb_numbering", 1 );
	no_heavy_unsats_allowed_ = tag->getOption<bool>( "no_heavy_unsats_allowed", 1 );
	show_task_ = tag->getOption<bool>( "show_task", 0 );
	min_network_size_ = tag->getOption<Size>( "min_network_size", 3 );
	max_network_size_ = tag->getOption<Size>( "max_network_size", 15 );
	min_unique_networks_ = tag->getOption<Size>( "min_unique_networks", 1 );
	min_core_res_ = tag->getOption<Size>( "min_core_res", 0 );
	min_boundary_res_ = tag->getOption<Size>( "min_boundary_res", 0 );
	max_unsat_ = tag->getOption<Size>( "max_unsat", 5 );
	dump_resfile_ = tag->getOption< bool >( "dump_resfile", 0 );
	start_from_csts_ = tag->getOption< bool >( "start_from_csts", 0);
	//use_enzdes_cst_ = tag->getOption< bool >( "use_enzdes_cst", 0 );

    if ( only_native_ ) find_native_ = true;
    if ( only_extend_existing_ ) extend_existing_networks_ = true;
//	if ( start_from_csts_ ) {
//		use_enzdes_cst_ = true;
//	}
//	if ( use_enzdes_cst_ ) {
//		basic::options::option[ basic::options::OptionKeys::run::preserve_header ].value(true);
//        write_cst_files_ = true;
//	}

	benchmark_ = tag->getOption<bool>( "benchmark", 0);
	verbose_ = tag->getOption<bool>( "verbose", 0);

	if ( tag->hasOption("design_residues") ) {
		des_residues_ = tag->getOption<std::string>("design_residues","STRKHYWNQDE");
	}

	// get task operations
	if ( tag->hasOption("task_operations") ) {
		task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	}

	if ( tag->hasOption("start_resnums") ) {
		std::string string_resnums( tag->getOption< std::string >( "start_resnums", "" ) );
		start_res_vec_ = pose::get_resnum_list( string_resnums, pose );
	}
	if ( tag->hasOption("start_selector") ) {
		if ( tag->hasOption("start_resnums" ) ) {
			if ( TR.visible() ) TR << "WARNING: cannot use both start_resnums and start_selector options; start_selector will be used" << std::endl;
			start_res_vec_.clear();
		}
		std::string const selector_name ( tag->getOption< std::string >( "start_selector" ) );
		if ( TR.visible() ) TR << "Set selector name to " << selector_name << "." << std::endl;
		core::select::residue_selector::ResidueSelectorCOP selector;
		try {
			selector = data.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from AddCompositionConstraintMover::parse_tag()\n" + e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_message );
		}
		runtime_assert( selector );
		start_selector_ = selector->clone();
	}
	if ( tag->hasOption("core_selector") ) {
		std::string const selector_name ( tag->getOption< std::string >( "core_selector" ) );
		if ( TR.visible() ) TR << "Set selector name to " << selector_name << "." << std::endl;
		core::select::residue_selector::ResidueSelectorCOP selector;
		try {
			selector = data.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from AddCompositionConstraintMover::parse_tag()\n" + e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_message );
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
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from AddCompositionConstraintMover::parse_tag()\n" + e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_message );
		}
		runtime_assert( selector );
		boundary_selector_ = selector->clone();
		
        //NEED to apply residue_selectros at apply time -- that is done now
	}
	if ( basic::options::option[ basic::options::OptionKeys::corrections::beta ].value(true) ) {
		init_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "HBNet_beta" );
	}

	core::scoring::ScoreFunctionOP new_score_function( protocols::rosetta_scripts::parse_score_function( tag, data ) );
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
	// Start the IG traversal at pre-deterimined starting residues (e.g. ligand or interface residues)
	for (unsigned long res1 : start_res_vec_) {

			platform::uint const first_ni = rotamer_sets_->resid_2_moltenres(res1);//returns 0 if res1 not molten (not design/repack-able)
		if ( first_ni == 0 ) continue;
		int const first_node_ind = (int)(first_ni);

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
			for ( int ii = 1; ii <= ig_->get_num_states_for_node(first_node_ind); ++ii ) {

				//if ( !(this->state_is_starting_aa_type( res1_ind, ii )) ) continue;

				//ResidueOP rot1(rotamer_sets_->rotamer_set_for_moltenresidue(first_ni)->nonconst_rotamer(ii)); //nonconst so we can add to its cache
				// if symmetric_ need to check one-body energy to make sure rotamer doesn't clash with itself accross the interface
				//    SymmetricRotamerSets compute_energies() adds 2-body e to the 1-body e for cases where symm residues that interact with own clones
				if ( symmetric_ ) {
					Real one_body_i = (ig_->get_one_body_energy_for_node_state( first_node_ind, ii ))/scmult_1b;
					if ( one_body_i > clash_threshold_ ) {
						continue;
					}
					// check for residues that hydrogen bond with their symmetric clones
					if ( one_body_i < onebody_hb_threshold_ ) {
						utility::vector1< HBondResStructCOP > residues(0);
						residues.push_back( HBondResStructCOP( new hbond_res_struct( res1_ind, (platform::uint)(ii), (rotamer_sets_->rotamer_set_for_moltenresidue(first_ni)->rotamer(ii)->name1()), orig_pose_->pdb_info()->chain(res1_ind), 1, 0, 0 ) ) );
						store_network( residues, one_body_i, true, true, false, false ); //network complete, store it
					}
				}
				for ( int jj = 1; jj <= ig_->get_num_states_for_node(second_node_ind); ++jj ) {
					//ResidueOP rot2(rotamer_sets_->rotamer_set_for_moltenresidue(second_ni)->nonconst_rotamer(jj)); //nonconst so we can add to its cache

					if ( !(this->pair_meets_starting_criteria( res1_ind, ii, res2_ind, jj )) ) continue;

					if ( symmetric_ ) {
						Real one_body_j = (ig_->get_one_body_energy_for_node_state( second_node_ind, jj ))/scmult_1b;
						if ( one_body_j > clash_threshold_ ) {
							continue;
						} else if ( one_body_j < onebody_hb_threshold_ ) {
							utility::vector1< HBondResStructCOP > residues(0);
							residues.push_back( HBondResStructCOP( new hbond_res_struct( res2_ind, (platform::uint)(jj), (rotamer_sets_->rotamer_set_for_moltenresidue(second_ni)->rotamer(jj)->name1()), orig_pose_->pdb_info()->chain(res2_ind), 1, 0, 0 ) ) );
							store_network( residues, one_body_j, true, true, false, false ); //network complete, store it
						}
					}
					PDEdge const & pdedge =  static_cast< PDEdge const & >(edge);
					Real twobody( (first_node_ind < second_node_ind) ? (pdedge.get_two_body_energy(ii, jj))/scmult_2b : (pdedge.get_two_body_energy(jj, ii))/scmult_2b );
					twobody = twobody * ( this->upweight_starting_twobody_energy() );
                    
					// if energy < threshold, we found an h-bond, if doesn't clash, add the rotamer to the current network
					if ( twobody < hb_threshold ) {
						core::Real init_sc = twobody;
						utility::vector1< HBondResStructCOP > residues(0);
						residues.push_back( HBondResStructCOP( new hbond_res_struct( res1_ind, (platform::uint)(ii), rotamer_sets_->rotamer_set_for_moltenresidue(first_ni)->nonconst_rotamer(ii)->name1(), orig_pose_->pdb_info()->chain(res1_ind), 1, 0, 0 ) ) );
						recursive_traverse( second_node_ind, jj, res2_ind, res1_ind, residues, 1, init_sc, hb_threshold ); //recursive traverse of IG
					}
				}
			}
		}
	}
}

void
HBNet::recursive_traverse( int const new_node_ind, int const newstate, Size const newres, Size const prevres,
	utility::vector1< HBondResStructCOP > residues, Size network_rec_count, Real init_sc, Real const hb_threshold, bool const second_search /* false */ )
{
	if ( !second_search ) {
		//ResidueOP new_rot( rotamer_sets_->rotamer_set_for_moltenresidue(new_node_ind)->nonconst_rotamer(newstate) ); //nonconst so we can add to its cache
		residues.push_back( HBondResStructCOP( new hbond_res_struct( newres, (platform::uint)(newstate), (rotamer_sets_->rotamer_set_for_moltenresidue(new_node_ind)->rotamer(newstate)->name1()), orig_pose_->pdb_info()->chain(newres), 1, 0, 0 ) ) );
	}
	Size prev_resnum = prevres;
	network_rec_count++; // numer of recursive calls; network_rec_count = # rotamers in the current h-bond network

	if ( network_rec_count >= max_network_size_ ) { //arbitrary limit for size of h-bond network to prevent explosion
		store_network( residues, init_sc, false, false, false, false );
		return;
	}
	Size instance_rec_call_cnt(0); // Keep track of how many times this particular instance makes new recursive calls

	// Traverse IG starting from current node/state
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
			//ResidueOP rot2(rotamer_sets_->rotamer_set_for_moltenresidue(second_ni)->nonconst_rotamer(jj)); //nonconst so we can add to its cache
			Real one_body_j = (ig_->get_one_body_energy_for_node_state( second_node_ind, jj ))/scmult_1b;
			// if symmetric interface design need to check one-body energy to make sure rotamer doesn't clash with itself accross the interface
			if ( symmetric_ && one_body_j > clash_threshold_ ) {
				continue;
			}
			PDEdge * pdedge = static_cast< PDEdge * >(*edge_iter);
			Real twobody( (new_node_ind < second_node_ind) ? (pdedge->get_two_body_energy(newstate, jj))/scmult_2b : (pdedge->get_two_body_energy(jj, newstate))/scmult_2b );

			if ( twobody < hb_threshold ) {
				init_sc +=  twobody;
				bool cycle = false;
				// if we get back to any of the starting residues, stop:
				//  it's more efficient to branch overlapping networks later than detect many duplicates and subnetworks through rec explosion
				if ( start_res_vec_.find( res2 ) == start_res_vec_.end() ) {
					// If new res does not clash with residues already in the network...
					if ( !( check_clash( residues, second_ni, jj, res2, init_sc, cycle ) ) ) {
						if ( cycle ) { // If found a cycle, store network
							store_network( residues, init_sc, false, true, false, false );
							//return; 16/03
						} else { // Else, add it to the network and keep going!  Recursive call to rec_trav
							//init_sc += one_body_j;
							recursive_traverse( second_node_ind, jj, res2, newres, residues, network_rec_count, init_sc, hb_threshold ); //recursive call
							instance_rec_call_cnt++;
						}
					}
				} else { // if we are back to original starting rotamer, don't need check_clash else check_clash
					int const start_state = static_cast<int>(residues.front()->rot_index);
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
	//if ( store_subnetworks_ || instance_rec_call_cnt == 0 )
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

///@details Recursivey traverse EnergyGraph of static pose to find all native networks
void
HBNet::traverse_native( Pose & pose, Real const hb_threshold )
{
	//if ( !( pose.energies().energies_updated() ) ){
	pose.update_residue_neighbors();
	( *scorefxn_ )( pose ); //score pose so we can get EnergyGraph (EG)

	for (unsigned long res : start_res_vec_) {
		utility::graph::Node * node = pose.energies().energy_graph().get_node( res );

		for ( utility::graph::EdgeListIterator egraph_it = node->edge_list_begin();
				egraph_it != node->edge_list_end(); ++egraph_it ) {
			core::scoring::EnergyEdge * eedge = static_cast< core::scoring::EnergyEdge * > (*egraph_it);
			Size other_res = eedge->get_other_node(node->get_node_index())->get_node_index();

			Size res1_ind(res),res2_ind(other_res);
			if ( symmetric_ ) {
				res1_ind = get_ind_res( *orig_pose_, res);
				res2_ind = get_ind_res( *orig_pose_, other_res);
			}

			if ( eedge->operator[](core::scoring::hbond_sc) < hb_threshold ) {
				// Initialize lists
				utility::vector1< HBondResStructCOP > residues(0);
				residues.push_back( HBondResStructCOP( new hbond_res_struct( res1_ind, 0, pose.residue(res1_ind).name1(), pose.pdb_info()->chain(res1_ind),
					pose.residue(res1_ind).is_protein(), 0, pose.residue(res1_ind).is_ligand() ) ) );

				rec_trav_native( pose, res2_ind, res1_ind, residues, hydrogen_bond_threshold_ ); // recursive call
			}
		}
	}
}

void
HBNet::rec_trav_native(Pose & pose, Size new_res, Size prev_res, utility::vector1< HBondResStructCOP > residues, Real const hb_threshold )
{
	bool found_new_res_in_net = false;
	for ( utility::vector1< HBondResStructCOP >::const_iterator i = residues.begin(); i != residues.end(); ++i ) {
		if ( (*i)->resnum == new_res ) {
			store_network(residues, 0.0, false, false, true );
			return;
		}
	}
	residues.push_back( HBondResStructCOP( new hbond_res_struct( new_res, 0, pose.residue(new_res).name1(), pose.pdb_info()->chain(new_res), pose.residue(new_res).is_protein(), 0, pose.residue(new_res).is_ligand() ) ) );

	utility::graph::Node * node = pose.energies().energy_graph().get_node( new_res );
	for ( utility::graph::EdgeListConstIterator egraph_it = node->const_edge_list_begin();
			egraph_it != node->const_edge_list_end(); ++egraph_it ) {
		core::scoring::EnergyEdge const * eedge = static_cast< core::scoring::EnergyEdge const * > (*egraph_it);
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

	for (const auto & residue : residues) {
		platform::uint otherstate( residue->rot_index );
		platform::uint other_node_ind(rotamer_sets_->resid_2_moltenres(residue->resnum));
		int const new_ind = static_cast<int>(new_node_ind);
		int const other_ind = static_cast<int>(rotamer_sets_->resid_2_moltenres(residue->resnum));
		int const new_st = static_cast<int>(newstate);
		int const other_st = static_cast<int>(otherstate);

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
		} else if ( ig_->find_edge(new_ind, other_ind) != 0 ) { // this function is order-independent; returns the same Edge regardless of order
			core::pack::interaction_graph::EdgeBase * temp_edge = ig_->find_edge(new_ind, other_ind);
			core::pack::interaction_graph::PDEdge * temp_pdedge = static_cast< PDEdge * >(temp_edge);
			int first_node_ind = temp_pdedge->get_first_node_ind(); // need to get the first node of the Edge (lower-index node) (we're getting this from the IG)
			platform::uint first_ni = first_node_ind;
			Size res1 = rotamer_sets_->moltenres_2_resid(first_ni); // get residue number that corresopnds to the first node of Edge

			//temp_pdedge->get_two_body_energy() is order-DEPENDENT; lower index (first index) must always be first)
			// so we need to determine does first_ni correspond to new_node_ind?
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
HBNet::net_clash(hbond_net_struct const & i, hbond_net_struct const & j)
{
	return net_clash(i.residues, j.residues);
}

///@details Function to check whether two h-bond networks clash with eachother
///    return of true = they clash, false = no clashes (networks i and j are compatible)
bool
HBNet::net_clash(utility::vector1< HBondResStructCOP > const & residues_i, utility::vector1< HBondResStructCOP > const & residues_j)
{
	core::PackerEnergy twobody(0.0);

	for (const auto & res_i : residues_i) {
		int const state_i = static_cast<int const>(res_i->rot_index);
		int const node_i = static_cast<int const>(rotamer_sets_->resid_2_moltenres(res_i->resnum));
		for (const auto & res_j : residues_j) {
			int const state_j = static_cast<int const>(res_j->rot_index);
			int const node_j = static_cast<int const>(rotamer_sets_->resid_2_moltenres(res_j->resnum));
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
						//Real sc_rmsd(core::scoring::automorphic_rmsd(*r_i, *r_j, true));
						//if ( core::scoring::automorphic_rmsd(*r_i, *r_j, true) > SC_RMSD_CUTOFF ) {
						// NEED TO COMPARE THESE 2
                        			if ( core::scoring::residue_sc_rmsd_no_super( r_i, r_j, true ) > SC_RMSD_CUTOFF ) {
							return true;
						}
					}
				} else {
					return true;
				}
			}

			if ( ig_ != nullptr  ) { //if IG still exists, use 2-body for faster lookup
				if ( ig_->find_edge(node_i, node_j) != 0 ) { // If edge exist s, look up 2-body energy for rotamer state_i and state_j
					core::pack::interaction_graph::EdgeBase * temp_edge = ig_->find_edge(node_i, node_j);
					core::pack::interaction_graph::PDEdge * temp_pdedge = static_cast< PDEdge * >(temp_edge);
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
	for (auto & i : network_vector_) {
		already_stored = network_already_stored( residues, i->residues );

		if ( already_stored ) {
			if ( init_score < i->score && !score_now ) {
				//replace i with new network with better score
                i->is_native = native;
				i->term_w_start = term_w_start;
				i->term_w_cycle = term_w_cycle;
				i->residues = residues;
				if ( (this->ligand()) ) {
					i->lig_state_list.push_back((*(find_hbond_res_struct(i->residues, this->ligand())))->rot_index);
				}
				i->score = init_score; //don't need to adjust for symmetric_ because did that when tracking IG energies to generate init_score
				//if (normalize_)
				//    (*i)->score = (*i)->score/(residues.size());
			} else if ( this->ligand() ) {
				i->lig_state_list.push_back((*(find_hbond_res_struct(i->residues, this->ligand())))->rot_index);
			}
			return;
		}
	}
	if ( !( already_stored ) ) {
		HBondNetStructOP new_net( new hbond_net_struct() );
        new_net->is_native = native;
		new_net->term_w_start = term_w_start;
		new_net->term_w_cycle = term_w_cycle;
		new_net->residues = residues;
		new_net->lig_state_list.clear();
		if ( (this->ligand()) ) {
			new_net->lig_state_list.push_back((*(find_hbond_res_struct((new_net)->residues, this->ligand())))->rot_index);
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
HBNet::score_networks( bool minimize )
{
	auto i = network_vector_.begin();

	for ( ; i != network_vector_.end(); ) {
		if ( !((*i)->scored) ) {
			Pose ala_copy = *ala_pose_;
			place_rotamers_and_score( ala_copy, **i ,minimize );
            //TODO NEED TO add charge-charge repulsion check
			++i;
		} else {
			++i;
		}
	}
}

//assumes residues are already on the pose and pose is scored
// RUNS BUT IS NOT USEFUL YET
bool
HBNet::has_charge_charge_repulsion( Pose & pose, utility::vector1< HBondResStructCOP > const & residues )
{
	for ( auto r1 = residues.begin(); r1 != residues.end(); ++r1 ) {
		//Real fa_elec(0.0);
		utility::graph::Node * node = pose.energies().energy_graph().get_node( (*r1)->resnum );
		for ( utility::graph::EdgeListConstIterator egraph_it = node->const_edge_list_begin();
				egraph_it != node->const_edge_list_end(); ++egraph_it ) {
			core::scoring::EnergyEdge const * eedge = static_cast< core::scoring::EnergyEdge const * > (*egraph_it);
			Size other_res = eedge->get_other_node(node->get_node_index())->get_node_index();

			bool not_network(true);
			for (const auto & residue : residues) {
				if ( other_res == residue->resnum ) not_network = false;
			}
			if ( not_network ) {
				continue;
			}

			Real fa_elec( eedge->operator[]( core::scoring::fa_elec ) );

			if ( fa_elec > charge_charge_rep_cutoff_ ) {
				//                TR << "fa_elec = " << fa_elec << "; res1 = " << (*r1)->resnum << "; other_res = " << other_res << std::endl;
				return true;
			}
		}
	}
	return false;
}

void
HBNet::merge_2_branched_networks(utility::vector1< HBondResStructCOP > const & residues1, utility::vector1< HBondResStructCOP > const & residues2, utility::vector1< HBondResStructCOP > & new_residues)
{
	for (const auto & res1 : residues1) {
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
	for (const auto & res2 : residues2) {
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
HBNet::merge_2_branched_networks(hbond_net_struct const & i, hbond_net_struct const & j, HBondNetStructOP new_network)
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
			new_lig_state_list.push_back((*(find_hbond_res_struct(i.residues, this->ligand())))->rot_index);
			new_lig_state_list.push_back((*(find_hbond_res_struct(j.residues, this->ligand())))->rot_index);
			new_network->lig_state_list = new_lig_state_list;
		} else {
			new_network->lig_state_list = new_lig_state_list;
		}
	}
	//score and num_unsat, etc. should be set by place_rotamers_and_score()
	new_network->lig_num_unsatisfied = 0;
	new_network->num_unsat = 0;
	new_network->num_heavy_unsat = 0;
	new_network->score = 0.0;
}

////TODO NEED TO remove this and come up with better solution for combining networks in HBNetStapleInterface
//void
//HBNet::merge_2_networks(hbond_net_struct const & i, hbond_net_struct const & j, HBondNetStructOP new_network)
//{
//	if ( (this->ligand()) && (i.lig_state_list.empty() || j.lig_state_list.empty()) ) {
//		return;
//	}
//
//	//for metrics used to sort the networks, compute avg (otherwise merged will always score worse than individual)
//    new_network->is_native = false;
//    new_network->is_extended = false;
//    new_network->term_w_bb = ( i.term_w_bb || j.term_w_bb ) ? true : false;
//	new_network->term_w_cycle = ( i.term_w_cycle || j.term_w_cycle ) ? true : false;
//	new_network->term_w_start = ( i.term_w_start || j.term_w_start ) ? true : false;
//	new_network->sort_first_by_tot_unsat = ( i.sort_first_by_tot_unsat || j.sort_first_by_tot_unsat );
//	new_network->scored = ( i.scored && j.scored );
//	new_network->outstring = i.outstring + "\n" + j.outstring;
//	new_network->num_unsat = (i.num_unsat + j.num_unsat)/2.0;
//	new_network->num_heavy_unsat = (i.num_heavy_unsat + j.num_heavy_unsat)/2.0;
//	new_network->lig_num_unsatisfied = (i.lig_num_unsatisfied + j.lig_num_unsatisfied)/2.0;//need better way for ligand
//	new_network->score = (i.score + j.score)/2.0;
//	new_network->residues = i.residues;
//	new_network->residues.insert( new_network->residues.end(), j.residues.begin(), j.residues.end() );
//	new_network->asymm_residues = i.asymm_residues;
//	new_network->asymm_residues.insert( new_network->asymm_residues.end(), j.asymm_residues.begin(), j.asymm_residues.end() );
//	new_network->unsat_Hpols = i.unsat_Hpols;
//	new_network->unsat_Hpols.insert( new_network->unsat_Hpols.end(), j.unsat_Hpols.begin(), j.unsat_Hpols.end() );
//	new_network->unsat_accs = i.unsat_accs;
//	new_network->unsat_accs.insert( new_network->unsat_accs.end(), j.unsat_accs.begin(), j.unsat_accs.end() );
//	new_network->hbond_vec = i.hbond_vec;
//	new_network->hbond_vec.insert( new_network->hbond_vec.end(), j.hbond_vec.begin(), j.hbond_vec.end() );
//	//new_network->rotamers = i.rotamers;
//	//new_network->rotamers.insert( new_network->rotamers.end(), j.rotamers.begin(), j.rotamers.end() );
//	//new_network->waterrots = i.waterrots;
//	//new_network->waterrots.insert( new_network->waterrots.end(), j.waterrots.begin(), j.waterrots.end() );
//
//	//NEED TO FIX
//	if ( this->ligand() ) {
//		std::vector<platform::uint> new_lig_state_list(1000);
//		std::vector<platform::uint> liglist1 = i.lig_state_list;
//		std::vector<platform::uint> liglist2 = j.lig_state_list;
//		std::sort( liglist1.begin(), liglist1.end() );
//		std::sort( liglist2.begin(), liglist2.end() );
//		std::vector<platform::uint>::iterator lit;
//		lit=std::set_intersection(liglist1.begin(), liglist1.end(), liglist2.begin(), liglist2.end(), new_lig_state_list.begin());
//		new_lig_state_list.resize(lit-new_lig_state_list.begin());
//		if ( new_lig_state_list.empty() ) {
//			new_lig_state_list.push_back((*(find_hbond_res_struct(i.residues, this->ligand())))->rot_index);
//			new_lig_state_list.push_back((*(find_hbond_res_struct(j.residues, this->ligand())))->rot_index);
//			new_network->lig_state_list = new_lig_state_list;
//		} else {
//			new_network->lig_state_list = new_lig_state_list;
//		}
//	}
//}

//consider Ser and Thr idential for benchmarking purposes
bool
HBNet::networks_identical_aa_sequence( hbond_net_struct const & i, hbond_net_struct const & j )
{
	if ( i.residues.size() == j.residues.size() ) {
		std::vector< char > i_aa(0);
		for (const auto & residue : i.residues) {
			if ( residue->aa == 'T' || residue->aa == 'S' ) {
				i_aa.push_back( 'S' );
			} else {
				i_aa.push_back( residue->aa );
			}
		}
		std::vector< char > j_aa(0);
		for (const auto & residue : j.residues) {
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

//identical networks are considered subsets of each other here by default unless specified by true_if_identical=false
//  std::includes() returns false if 2 vectors are identical
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

//find networks that can be branched together (e.g. branch off the same single residue and should really be the same network)
void
HBNet::branch_overlapping_networks()
{
	if ( network_vector_.empty() ) {
		if ( verbose_ && TR.visible() ) {
			TR << " NO NETWORKS TO BRANCH; RETURNING WITHOUT DOING ANYTHING." << std::endl;
		}
		return;
	}
	//Size icount(1);
	merged_vecs_.clear(); //unnecessary if function only used internally
	for (auto & i : network_vector_) {
		i->net_indices.clear();
	}
	for ( auto i = network_vector_.begin(); i != network_vector_.end(); ++i ) {
		//Size jcount(1);
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
			//jcount++;
		}
		//icount++;
	}
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
				merged_vecs_.push_back(add_index_vec);
			}

			std::sort( network_vector_[*j]->net_indices.begin(),  network_vector_[*j]->net_indices.end());
			std::sort( (*i)->net_indices.begin(), (*i)->net_indices.end() );
			std::vector< Size > v(1000);
			std::vector< Size >::iterator it;
			it=std::set_intersection((*i)->net_indices.begin(), (*i)->net_indices.end(), (network_vector_[*j])->net_indices.begin(), (network_vector_[*j])->net_indices.end(), v.begin());
			v.resize(it-v.begin());

			if ( v.size() == 0 ) {
				if ( !store_subnetworks_ ) {
					merged_vecs_.push_back(add_index_vec);
				}
				add_index_vec.pop_back();
				continue;
			} else if ( v.size() == 1 ) {
				add_index_vec.push_back( v.front() );
				merged_vecs_.push_back(add_index_vec);
				add_index_vec.pop_back();
			} else {
				rec_set_intersection( add_index_vec, v, *j );
			}
			add_index_vec.pop_back();
		}
	}
	std::vector< std::vector< Size > > temp_merge_vec(0);
	Size mcount(1);
	for (auto & merged_vec : merged_vecs_) {
		if ( merged_vec.empty() ) continue; //should never happen
		std::sort(merged_vec.begin(), merged_vec.end());
		bool found(false);
		for (auto & t : temp_merge_vec) {
			std::sort(t.begin(), t.end());
			std::vector< int > v2(50);
			auto it2 = std::set_symmetric_difference(merged_vec.begin(), merged_vec.end(), t.begin(), t.end(), v2.begin());
			v2.resize(it2-v2.begin());

			if ( v2.size() == 0 ) {
				found = true;
				break;
			}
		}
		if ( !found ) {
			temp_merge_vec.push_back(merged_vec);
		}
	}
	merged_vecs_.clear();
	for ( std::vector< std::vector< Size > >::const_iterator mit = temp_merge_vec.begin(); mit != temp_merge_vec.end(); ++mit ) {
		merged_vecs_.push_back( *mit );
	}
	std::vector< HBondNetStructOP > temp_net_vec(0);
	for (auto & merged_vec : merged_vecs_) {
		if ( merged_vec.size()>1 ) {
			auto m1 = merged_vec.begin();
			HBondNetStructOP temp_hbond_net_struct = network_vector_[*m1];
			std::string nets_being_merged("");
			if ( verbose_ ) {
				nets_being_merged = "       "+utility::to_string(mcount)+". branching together networks = ";
				nets_being_merged = nets_being_merged+" ( #"+(utility::to_string(*m1+1))+": ";\
					std::string tempstr = print_list_to_string(temp_hbond_net_struct->residues, false, false, false, false);
				//std::string tempstr = print_list_to_string(temp_hbond_net_struct.rotlist, temp_hbond_net_struct.term_w_start,
				//                                           temp_hbond_net_struct.term_w_cycle, temp_hbond_net_struct.term_w_bb );
				nets_being_merged = nets_being_merged+tempstr+" ) + ";
			}
			for ( auto m2 = m1; ++m2 != merged_vec.end(); ) {
				if ( verbose_ ) {
					nets_being_merged = nets_being_merged+" ( #"+(utility::to_string((*m2+1)))+": ";
					std::string tempstr = print_list_to_string(network_vector_[*m2]->residues, false, false, false, false);
					//std::string tempstr = print_list_to_string(network_vector_[*m2].rotlist, network_vector_[*m2].term_w_start,
					//                                           network_vector_[*m2].term_w_cycle, network_vector_[*m2].term_w_bb);
					nets_being_merged = nets_being_merged+tempstr+" ) + ";
				}
				HBondNetStructOP new_network( new hbond_net_struct() );
				merge_2_branched_networks( *temp_hbond_net_struct, *(network_vector_[*m2]), new_network );
				temp_hbond_net_struct = new_network;
			}
			if ( verbose_ && TR.visible() ) {
				TR << nets_being_merged << std::endl;
			}

			temp_net_vec.push_back( temp_hbond_net_struct );
			mcount++;
		}
	}
	//push the merged networks to the back of network_vector_
	//network_vector_.clear(); DON'T WANT TO CLEAR HERE
	for (auto & netit : temp_net_vec) {
		network_vector_.push_back( netit );
	}
	// Sort all individual networks by best score, so start with highest scoring
	std::sort( network_vector_.begin(), network_vector_.end(), compare_net_vec() );
	merged_vecs_.clear();
}//branch_overlapping networks

//used by branch_overlapping() to efficiently search for all combinations of compatible networks that can be merged
void
HBNet::rec_set_intersection( std::vector< Size > add_index_vec, std::vector< Size > next_index_vec, Size pos )
{
	for ( auto i = next_index_vec.begin(); i != next_index_vec.end(); ++i ) {
		if ( *i <=  pos  ) {
			continue;
		}
		add_index_vec.push_back( *i );
		if ( store_subnetworks_ ) {
			merged_vecs_.push_back(add_index_vec);
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
				merged_vecs_.push_back(add_index_vec);
			}
			add_index_vec.pop_back();
			continue;
		} else if ( v.size() == 1 ) {
			add_index_vec.push_back(v.front());
			merged_vecs_.push_back(add_index_vec);
			add_index_vec.pop_back();
		} else {
			rec_set_intersection( add_index_vec, v, *i );
		}
		add_index_vec.pop_back();
	}
}

void
HBNet::place_rots_on_pose( pose::Pose & pose, utility::vector1< HBondResStructCOP > const & residues, bool use_pose /* false */ )
{
	//need better solution here than re-copying residue datacache!
	if ( use_pose ) {
		for (const auto & residue : residues) {
			pose.replace_residue( residue->resnum, orig_pose_->residue( residue->resnum ), false );
		}
	} else {
		assert( rotamer_sets_ != nullptr );
		for (const auto & residue : residues) {
			ResidueCOP copy_rot(rotamer_sets_->rotamer_for_moltenres(rotamer_sets_->resid_2_moltenres((platform::uint)(residue->resnum)),residue->rot_index));
			pose.replace_residue( residue->resnum, *copy_rot, false );
		}
	}
}

void
HBNet::find_unsats( Pose const & pose, hbond_net_struct & network )
{
	Size num_unsatisfied(0);
	Size num_heavy_unsat(0);

	runtime_assert(!(network.hbond_vec.empty()));

	if ( verbose_ && TR.visible() ) network.hbond_set->show(TR);

	//TODO NEED TO ADD option to use own bunsat calc here
    
	Size polar_atom_count(0);
	Size polar_atms_making_hbonds(0);
	for ( utility::vector1< HBondResStructCOP >::const_iterator res = network.residues.begin(); res != network.residues.end(); ++res ) {
		//for ( std::set< Size >::const_iterator res = actual_hbond_residues.begin(); res != actual_hbond_residues.end(); ++res ){
		Size const resnum( (*res)->resnum );
		//Size const resnum( *res );
		//Size unsats_so_far(num_unsatisfied);

		if ( verbose_ && TR.visible() ) network.hbond_set->show( pose, resnum, true, TR );

		for ( Size a_index = 1; a_index <= pose.residue(resnum).natoms(); ++a_index ) {
			if ( !( pose.residue( resnum ).atom_is_backbone( a_index ) ) || ((this->ligand()) && resnum == (this->ligand()) ) ) {

				if ( pose.residue( resnum ).atom_type( a_index ).is_donor() && pose.residue( resnum ).atomic_charge( a_index ) != 0.0 ) {
					if ( !(pose.residue( resnum ).atom_type(a_index).name() == "OH") ) { //do not want to double penalize hydroxyls, checked both O and H above
						Size h_count(0);
						Size h_unsat(0);
						for ( Size hatm = pose.residue(resnum).attached_H_begin(a_index); hatm <= pose.residue(resnum).attached_H_end(a_index); ++hatm ) {
							//TR.Debug << "hatm = " << pose.residue(resnum).atom_type(hatm).atom_type_name() << std::endl;
							h_count++;
							polar_atom_count++;
							id::AtomID const at_id( hatm, resnum);

							if ( !( (network.hbond_set)->nhbonds( at_id, false /* include_only_allowed */ ) ) ) { // unsat
								//if ( res_is_core( resnum ) || atom_sasa[ at_id ] < atom_burial_cutoff_ ){
								if ( res_is_core( resnum ) ) {
									num_unsatisfied++;
									h_unsat++;
									network.unsat_Hpols.push_back( at_id );
									TR.Debug << " res " << resnum << " atom " << hatm << " " << pose.residue(resnum).atom_name(hatm) << " is unsat" << std::endl;
								}
							} else {
								polar_atms_making_hbonds++;
							}
						}
						if ( h_unsat == h_count ) {
							if ( pose.residue(resnum).atom_type( a_index ).name() == "OH" ) {
								if ( hydroxyls_must_donate_ || ( pose.residue(resnum).name1() == 'Y' && tyr_hydroxyls_must_donate_ ) ) {
									num_heavy_unsat++;
									TR.Debug << " res " << resnum << " atom " << a_index << " " << pose.residue(resnum).atom_name(a_index) << " is heavy unsat" << std::endl;
								}
							} else {
								num_heavy_unsat++;
								TR.Debug << " res " << resnum << " atom " << a_index << " " << pose.residue(resnum).atom_name(a_index) << " is heavy unsat" << std::endl;
							}
							num_unsatisfied++;
						}
					}
				}
				if ( pose.residue( resnum ).atom_type( a_index ).is_acceptor() && pose.residue( resnum ).atomic_charge( a_index ) != 0.0 ) { //important for ligand case
					polar_atom_count++;
					id::AtomID const at_id( a_index, resnum);

					if ( !( (network.hbond_set)->nhbonds( at_id , false /* include_only_allowed */ ) ) ) { // unsat
						if ( res_is_core( resnum ) ) {
							if ( pose.residue(resnum).atom_type( a_index ).name() == "OH" ) {
								core::Size hatm = pose.residue(resnum).attached_H_begin( a_index );
								if ( !( (network.hbond_set)->nhbonds( id::AtomID( hatm, resnum ), false /* include_only_allowed */ )) ) { // only count this atom as unsat if we haven't already counted the hydrogen
									num_unsatisfied++;
									num_heavy_unsat++;
									network.unsat_accs.push_back( at_id );
									TR.Debug << "res " << resnum << " atom " << a_index << " " << pose.residue(resnum).atom_name( a_index ) << " is heavy unsat" << std::endl;
								}
							} else {
								num_unsatisfied++;
								num_heavy_unsat++;
								network.unsat_accs.push_back( at_id );
								TR.Debug << " res " << resnum << " atom " << a_index << " " << pose.residue(resnum).atom_name(a_index) << " is heavy unsat" << std::endl;
								TR.Debug << " res " << resnum << " atom " << a_index << " " << pose.residue(resnum).atom_name(a_index) << " is unsat" << std::endl;
								//TR << " res_is_core( resnum ) = " << res_is_core( resnum ) << " and atom_sasa[ at_id ] = " << atom_sasa[ at_id ] << std::endl;
							}
						}
					} else {
						polar_atms_making_hbonds++;
					}
				}
			}
		}
	}
	network.num_unsat = num_unsatisfied;
	network.num_heavy_unsat = num_heavy_unsat;
	network.connectivity = (Real)polar_atms_making_hbonds/(Real)polar_atom_count;
	//TR << "polar_atms_making_hbonds = " << polar_atms_making_hbonds << "; polar_atom_count = " << polar_atom_count << "; connectivity = " << network.connectivity << std::endl;
}

core::PackerEnergy
HBNet::place_rotamers_and_score(pose::Pose & pose, hbond_net_struct & i, bool minimize/*=false*/)
{
	pose.update_residue_neighbors();
	//if (!(pose.energies().energies_updated()))
	(*scorefxn_)(pose);
	core::scoring::EnergyMap & baseline_emap = pose.energies().total_energies();
	core::PackerEnergy total_baseline = baseline_emap.get(core::scoring::total_score);
	core::PackerEnergy hbond_bb_sc_baseline = baseline_emap.get(core::scoring::hbond_bb_sc);

	place_rots_on_pose( pose, i.residues, i.is_native );
	pose.update_residue_neighbors();
	//how to handle h2o water here? need jump to be connected to correct atom?
    //TODO NEED TO ADD BETTER OPTIONS FOR MINIMIZATION
	if ( minimize ) { //minimize chi's of network residues with hack sfxn
		core::kinematics::MoveMapOP mm = core::kinematics::MoveMapOP( new core::kinematics::MoveMap );
		mm->set_bb( false ); mm->set_chi( false ); //mm->set_chi(repack);
		for ( utility::vector1< HBondResStructCOP >::const_iterator rit=i.residues.begin(); rit != i.residues.end(); rit++ ) {
			//            if ( (this->ligand()) && !minimize_lig_ && (*rit)->resnum == (this->ligand()) )
			//                continue;
			if ( (this->ligand()) && (*rit)->resnum == (this->ligand()) ) {
				continue;
			}
			mm->set_chi( (*rit)->resnum, true );
		}
		protocols::simple_moves::MinMoverOP min_mover;
		if ( symmetric_ ) {
			core::pose::symmetry::make_symmetric_movemap( pose, *mm );
			min_mover = protocols::simple_moves::MinMoverOP( new protocols::simple_moves::symmetry::SymMinMover( mm, scorefxn_, "dfpmin_armijo_nonmonotone", 0.0001, true ) );
		} else {
			min_mover = protocols::simple_moves::MinMoverOP( new protocols::simple_moves::MinMover( mm, scorefxn_, "dfpmin_armijo_nonmonotone", 0.0001, true ) );
		}
		min_mover->apply(pose);
	}
	pose.update_residue_neighbors(); // Score pose
	(*scorefxn_)(pose);
	core::scoring::EnergyMap & new_emap = pose.energies().total_energies();
	PackerEnergy pose_total = new_emap.get(core::scoring::total_score) - total_baseline;
	PackerEnergy e( pose_total );
	if ( (new_emap.get( core::scoring::hbond_bb_sc) - hbond_bb_sc_baseline) <= MIN_HB_E_CUTOFF ) {
		i.term_w_bb = true;
	} else {
		i.term_w_bb = false;
	}

	//NEED TO FIX
	core::Real cycle_bonus = 2.0;
	if ( i.term_w_cycle ) e = e + cycle_bonus;
	//if ( (this->ligand()) && i.term_w_start ) e = e + ret2lig_bonus_;
	if ( symmetric_ ) e = e/(symm_info_->num_bb_clones() + 1.0);
	i.score = e;
	i.scored = true;
	return e;
}

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
HBNet::networks_unique( hbond_net_struct const & i, hbond_net_struct const & j, bool no_surface /* true */ )
{
	utility::vector1< HBondResStructCOP > residues_i(0);
	utility::vector1< HBondResStructCOP > residues_j(0);

	if ( no_surface ) {
		for (const auto & residue : i.residues) {
			if ( res_is_core(residue->resnum) || res_is_boundary(residue->resnum) ) {
				residues_i.push_back( residue );
			}
		}
		for (const auto & residue : j.residues) {
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
void
HBNet::remove_replicate_networks( Size same_max/*=1*/ )
{
	std::sort( network_vector_.begin(), network_vector_.end(), compare_net_vec() ); //sort all networks by score
	std::vector< HBondNetStructOP > temp_net_vec(0);
	std::vector< Size > net_count(0);
	if ( verbose() && TR.visible() ) TR << "         DELETING REPLICATE NETWORKS THAT HAVE IDENTIAL RESNUM AND AA SEQ; ONLY KEEPING BEST SCORING ONE" << std::endl;
	for (auto & i : network_vector_) {
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
		if ( only_water ) { // NEED TO implement this better
			continue;
		}
		bool reached_limit = false;
		Size ind_count = 0;

		runtime_assert( !(i->residues.empty()) );

		for (auto & vit : temp_net_vec) {
			if ( residues_identical( i->residues, vit->residues ) ) {
				net_count[ind_count] = net_count[ind_count] + 1;
				if ( net_count[ind_count] > same_max ) {
					reached_limit = true;
				}
			} else {
				continue;
			}
			ind_count++;
		}
		if ( reached_limit ) {
			if ( verbose_ && TR.visible() ) TR << "reached limit, continuing: network " << print_list_to_string(i->residues) << std::endl;
			continue;
		}
		net_count.push_back(1);
		temp_net_vec.push_back( i );
	}
	//clear network_vector_ then push back all of the ones that passed the cutoff (from temp_net_vec)
	network_vector_.clear();
	for (auto & k : temp_net_vec) {
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
	for (const auto & residue : residues) {
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
	for (const auto & residue : residues) {
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
	for ( auto i = network_vector_.begin(); i != network_vector_.end(); ++i ) { // outer loop
        // for this options case, need to check if existing (native) networks have been extended
        //  if they have, choose either native or extended, whichever is better
        bool skip_this_network = false;
        if ( extend_existing_networks_ ){
            auto nit( native_networks_.begin() );
            for ( ; nit != native_networks_.end(); ) { // inner loop
                bool dummy(false);
                // if the native network is a subset
                if ( is_sub_residues( (*nit)->residues, (*i)->residues, dummy, true ) ){
                    // if extended is bigger (and meets criteria) or same size but few unsats, then erase native
                    //   bool operator< overridden by hbond_net_struct to define best network
                    Size i_size( (symmetric_ && !((*i)->asymm_residues.empty()) ) ? (*i)->asymm_residues.size() : (*i)->residues.size() );
                    Size nit_size( (symmetric_ && !(*nit)->asymm_residues.empty()) ? (*nit)->asymm_residues.size() : (*nit)->residues.size() );
                    if ( i_size > nit_size ||
                        ( i_size == nit_size && (*i)->num_unsat < (*nit)->num_unsat ) )
                    {
                        (*i)->is_extended = true; // record that this is an extended network
                        nit = native_networks_.erase( nit );
                    } else {
                        // keep native and don't keep new network; we won't add it to output_net_vec_
                        // skip this networks (it won't get returned)
                        skip_this_network = true;
                        break; //because we found a native that's better, exit inner loop
                        ++nit;
                    }
                }
                else
                    ++nit;
            }
            //logic for determing if network is native or extended from staring HBNet PDBInfoLabels:
            if ( skip_this_network ) continue; // move on to the next new network, exit outer loop
            Size hbnet_label_count(0);
            for (auto & residue : (*i)->residues) {
                if ( residue->resnum <= hbnet_info_residues_.size() && hbnet_info_residues_[residue->resnum] ){
                    hbnet_label_count++;
                }
            }
            if ( hbnet_label_count == (*i)->residues.size() )
                (*i)->is_native = true;
            else if ( hbnet_label_count > 2 )
                (*i)->is_extended = true;
        }
		(*i)->id = count;
        //TR << "(*i)->id =" << (*i)->id << " and count = " << count << std::endl;
        //std::string outstring( ((*i)->is_extended) ? "extended" : "" );
		std::string outstring( ( pdb_numbering() ) ? print_network_w_pdb_numbering( get_orig_pose(), **i, true ) : print_network( **i ) );
		outstring = outstring + this->print_additional_info_for_net( **i );
		(*i)->outstring = outstring;
        
        if ( write_network_pdbs_ ) write_network_pdb( *i );
        
        if ( TR.visible() ) TR << "\t" << outstring << std::endl;
		if ( finalize ) {
            std::vector< Size > net_ids(0);
            net_ids.push_back( (*i)->id );
            
            if ( keep_existing_networks_ ){
                if ( ( (*i)->is_native || (*i)->is_extended ) ){
                    for ( auto j = i; ++j != network_vector_.end(); ){
                        if ( (*j)->is_native || (*j)->is_extended ){
                            bool compatible( true );
                            for (unsigned long & net_id : net_ids){
                                if ( net_clash( *(get_network_by_id(net_id)), **j ) ){
                                    compatible = false;
                                    break;
                                }
                            }
                            if ( compatible ) net_ids.push_back( (*j)->id );
                        }
                    }
                }
                else {
                    for (auto & k : network_vector_) {
                        if ( k->is_native || k->is_extended ){
                            bool compatible( true );
                            for (unsigned long & net_id : net_ids){
                                if ( net_clash( *(get_network_by_id(net_id)), *k ) ){
                                    compatible = false;
                                    break;
                                }
                            }
                            if ( compatible ) net_ids.push_back( k->id );
                        }
                    }
                    
                }
            }
            output_vector_.push_back( net_ids );
			//output_net_vec_.push_back( *i );
		}
		count++;
    }
    if ( find_native_ ){
        if ( TR.visible() ) TR << std::endl << " Keeping the following existing (native) networks that meet your criteria: "
            << std::endl << "\t" << print_headers() << print_additional_headers() << std::endl;
        std::sort( native_networks_.begin(), native_networks_.end(), compare_net_vec() ); //sort native networks
        Size count(1);
        std::string comment_str( print_headers() + this->print_additional_headers() );
        for (auto & native_network : native_networks_) {
            native_network->id = count;
            std::string outstring( ( pdb_numbering() ) ? print_network_w_pdb_numbering( get_orig_pose(), *native_network, true ) : print_network( *native_network ) );
            outstring = outstring + this->print_additional_info_for_net( *native_network );
            native_network->outstring = outstring;
            if ( TR.visible() ) TR << "\t" << outstring << std::endl;
            
            // if these options are set, we want to retain the native networks and csts in each pose that HBNet passes back
            //   call private nonconst_get_orig_pose() and modify the original from which all other poses make deep copies before being modified and passed back
            if ( keep_existing_networks_ || only_native_ ){
                add_reslabels_to_pose( nonconst_get_orig_pose(), *native_network );
                core::scoring::constraints::ConstraintSetOP cst_op( nonconst_get_orig_pose().constraint_set()->clone() ); // get existing csts and add to them
                set_constraints( nonconst_get_orig_pose(), *cst_op, native_network, write_cst_files_ );
                nonconst_get_orig_pose().constraint_set(cst_op); // add constraints to the pose
                comment_str = comment_str + "\n" + native_network->outstring;
            }
            count++;
        }
        if ( keep_existing_networks_ || only_native_ )
            core::pose::add_comment( nonconst_get_orig_pose(), "HBNet Native details: ", "\n" + comment_str + "\n" );
        if ( TR.visible() )
            TR << std::endl;
    }
}// output_networks

core::pack::task::PackerTaskOP
HBNet::create_ptask(core::pose::Pose & pose, bool initialize_from_commandline/*=false*/)
{
	using namespace core::pack::task;
	TR<<" Creating packer task based on specified task operations..."<< std::endl;
	runtime_assert( task_factory_ != nullptr );
	if ( initialize_from_commandline ) {
		task_factory_->push_back( core::pack::task::operation::TaskOperationOP( new operation::InitializeFromCommandline ) );
	}
	PackerTaskOP task = task_factory_->create_task_and_apply_taskoperations( pose );
	return task;
}

//copied from PackRotamersMover
bool
HBNet::task_is_valid(Pose const & pose) const
{
	if ( task_->total_residue() != pose.total_residue() ) return false;
	for ( Size i(1); i <= pose.total_residue(); ++i ) {
		chemical::ResidueTypeCOP r = pose.residue_type(i).get_self_ptr();
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
	core::conformation::symmetry::SymmetricConformation const & SymmConf(dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation()));
	symm_info_ = SymmConf.Symmetry_Info();
	//max_network_size_ = max_network_size_*symm_info_->num_bb_clones();
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
HBNet::get_ind_res( Pose & pose, Size res_i)
{
	Size resi_ind(res_i);
	if ( symmetric_ && res_i > symm_info_->num_independent_residues() ) {
		char resi_chain = pose.chain(res_i);
		if ( multi_component_ ) {
			//symm_info_->component_lower_bound()
			std::map<char,std::pair<Size,Size> > component_bounds = symm_info_->get_component_bounds();
			char resi_comp = symm_info_->get_component_of_residue(res_i);
			resi_ind = res_i - chain_bounds_[resi_chain].first + component_bounds[resi_comp].first;
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
HBNet::write_network_pdb( HBondNetStructOP p ) // better to pass object pointer or reference?
{
    if ( p->network_pdb_written ) return;
    
    runtime_assert( p->id != 0 );
    Pose ala_copy = *ala_pose_;
    place_rots_on_pose( ala_copy, p->residues, p->is_native );
    //place_waters_on_pose( ala_copy, p );
    std::string current_out_tag = protocols::jd2::JobDistributor::get_instance()->current_output_name();
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
    add_reslabels_to_pose( ala_copy, *p );
    
    ala_copy.dump_pdb(pdb_tag);
    
    p->network_pdb_written = true;
}

// assumes the network rotamers are already placed onto the pose
void
HBNet::set_constraints( Pose & pose, core::scoring::constraints::ConstraintSet & constraints, HBondNetStructOP p, bool write_cst_file /* false */ )
{
    runtime_assert( !(p->hbond_vec.empty()) );
    
    utility::io::ozstream cst_output_stream;
    std::string cst_fname("");
    
    if ( write_cst_file ){
        std::string prefix("");
        if ( p->is_native ) prefix = "native_";
        else if ( p->is_extended ) prefix = "extended_";
        cst_fname = get_file_name( p->id, prefix, ".cst" );
        core::pose::add_comment( pose, "cst_filename for "+prefix+"network_"+utility::to_string(p->id), cst_fname );
        cst_output_stream.open(cst_fname, std::ios_base::out);
        cst_output_stream << "# " << cst_fname << std::endl;
        cst_output_stream << "# " << p->outstring << std::endl;
    }
    for ( utility::vector1< HBondCOP >::const_iterator i = p->hbond_vec.begin(); i != p->hbond_vec.end(); ++i ) {
        Real xtol(0.20);
        Size drsd((*i)->don_res());
        Size don_hatm((*i)->don_hatm());
        Size datm(pose.residue(drsd).atom_base((int)(don_hatm)));
        Size arsd((*i)->acc_res());
        Size aatm((*i)->acc_atm());
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
        core::scoring::constraints::AtomPairConstraintCOP hbond_cst( new core::scoring::constraints::AtomPairConstraint( acc_id, don_id, dist_func ) );
        constraints.add_constraint( hbond_cst );
        if ( write_cst_file && !(p->cst_file_written) ){
            //(*i)->show(in_pose, 1, TR);
            Real angle_AHD((*i)->get_AHDangle(pose));
            Real angle_BAH((*i)->get_BAHangle(pose));
            Real dist_AH((*i)->get_HAdist(pose));
            Real hb_energy((*i)->energy());//unweighted h-bond energy
            std::ostringstream cst_str_stream;
            cst_str_stream << "# angle_AHD = " << angle_AHD << ", angle_BAH = " << angle_BAH << ", dist_AH = " << dist_AH << angle_BAH << ", unweighted hb_energy = " << hb_energy << "\n";
            cst_str_stream << "AtomPair " << pose.residue(arsd).atom_name((int)(aatm)) << " " << arsd << " " << pose.residue(drsd).atom_name((int)(datm)) << " " << drsd
            << " BOUNDED " << lb << " " << ub << " " << sd << " " << tag << std::endl;
            cst_output_stream << cst_str_stream.str();
        }
    }
    if ( write_cst_file ){
        cst_output_stream.close();
        p->cst_file_written = true;
    }
}

//fast check if water ideal O position clashes with any existing atoms in the pose (templated from Chris King's SemiExplicitWaterBuriedUnsatCalc)
bool
HBNet::water_clashes(
	Pose const & pose,
    utility::graph::GraphOP packer_neighbor_graph,
	Size const anchor_i,
	Size const anchor_j,
	Vector const water_O
	//Real const clash_dist_cut /* 1.5 */
)
{
	Real const waterO_2_Hpol_cutoff( 1.2 );
	Real const waterO_2_Hpol_cutoff_sq( waterO_2_Hpol_cutoff * waterO_2_Hpol_cutoff );
	Real const waterO_cutoff( 2.4 );
	Real const waterO_cutoff_sq( waterO_cutoff * waterO_cutoff );

    for ( utility::graph::Graph::EdgeListConstIter
			ir  = packer_neighbor_graph->get_node( anchor_i )->const_edge_list_begin(),
			ire = packer_neighbor_graph->get_node( anchor_i )->const_edge_list_end();
			ir != ire; ++ir ) {
		Size const res2( (*ir)->get_other_ind( anchor_i ) );

		if ( water_O.distance_squared( pose.residue( res2 ).nbr_atom_xyz()) > 144.0 ) { /* 12^2 */
			continue;
		}
		for ( Size at2 = 1; at2 <= pose.residue( res2 ).natoms(); ++at2 ) {
			//skip virtual atoms!
			if ( pose.residue( res2 ).atom_type( at2 ).lj_wdepth() == 0.0 ) continue;
			id::AtomID atid2( at2, res2 );

			Real const dist2( water_O.distance_squared( pose.xyz( atid2 ) ) );
			if ( pose.residue( res2 ).atom_is_polar_hydrogen( at2 ) && dist2 < waterO_2_Hpol_cutoff_sq ) {
				return true;
			} else if ( dist2 < waterO_cutoff_sq ) {
				return true;
			}
		}
	}
    for ( utility::graph::Graph::EdgeListConstIter
			ir  = packer_neighbor_graph->get_node( anchor_j )->const_edge_list_begin(),
			ire = packer_neighbor_graph->get_node( anchor_j )->const_edge_list_end();
			ir != ire; ++ir ) {
		Size const res2( (*ir)->get_other_ind( anchor_j ) );

		if ( water_O.distance_squared( pose.residue( res2 ).nbr_atom_xyz()) > 144.0 ) { /* 12^2 */
			continue;
		}
		for ( Size at2 = 1; at2 <= pose.residue( res2 ).natoms(); ++at2 ) {
			//skip virtual atoms!
			if ( pose.residue( res2 ).atom_type( at2 ).lj_wdepth() == 0.0 ) continue;
			id::AtomID atid2( at2, res2 );

			Real const dist2( water_O.distance_squared( pose.xyz( atid2 ) ) );
			if ( pose.residue( res2 ).atom_is_polar_hydrogen( at2 ) && dist2 < waterO_2_Hpol_cutoff_sq ) {
				return true;
			} else if ( dist2 < waterO_cutoff_sq ) {
				return true;
			}
		}
	}
	return false;
}
    
std::string
HBNet::get_file_name( Size id, std::string net_prefix, std::string extension )
{
    std::string current_out_tag = protocols::jd2::JobDistributor::get_instance()->current_output_name();
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
	conformation::Residue const & water
)
{
	Real const hbond_cutoff( 1.2 );
	Real const hbond_cutoff_sq( hbond_cutoff * hbond_cutoff );
	Real const cutoff( 1.5 );
	Real const cutoff_sq( cutoff * cutoff );

	for ( Size iatid = 1; iatid <= water.natoms(); ++iatid ) {
		Vector const at1_xyz(water.xyz(iatid));
		for ( Size res2 = 1; res2 <= pose.total_residue(); ++res2 ) {
			if ( at1_xyz.distance_squared( pose.residue( res2 ).nbr_atom_xyz()) > 144.0 ) { /* 12^2 */
				continue;
			}
			for ( Size at2 = 1; at2 <= pose.residue( res2 ).natoms(); ++at2 ) {
				//skip virtual atoms!
				if ( pose.residue( res2 ).atom_type( at2 ).lj_wdepth() == 0.0 ) continue;
				id::AtomID atid2( at2, res2 );
				Real const dist2( at1_xyz.distance_squared( pose.xyz( atid2 ) ) );
				if ( pose.residue( res2 ).atom_is_polar_hydrogen( at2 ) && water.atom_name(iatid) == "O" && dist2 < hbond_cutoff_sq ) {
					return true;
				} else if ( pose.residue( res2 ).atom_type(at2).is_acceptor() && (water.atom_name(iatid) == "H1" || water.atom_name(iatid) == "H2") && dist2 < hbond_cutoff_sq ) {
					return true;
				} else if ( dist2 < cutoff_sq ) {
					return true;
				}

			}
		}
	}
	return false;
}

pose::PoseOP
HBNet::get_additional_output()
{
	pose::PoseOP out_pose;
    if ( output_vector_.size() == 0 ){
	//if ( output_net_vec_.size() == 0 ) {
		//set_last_move_status(protocols::moves::FAIL_RETRY); TODO NEED TO CHECK THIS IS CORRECT want this commented
		//        if (kill_if_no_more_networks_found()) // if true, will kill the rest of the nstruct runs once all h-bond networks are tried
		//            complete(true);
		return NULL;
    } else if ( output_vector_.size() == 1 ){
	//} else if ( output_net_vec_.size() == 1 ) { // Last iteration, do not create a new pose reference
		out_pose = orig_pose_;
	} else {
		out_pose = pose::PoseOP( new Pose( *orig_pose_ ) ); // May produce additional output poses, create copy of the target
        // TODO double check that this preserves pdb_info and res_labels
	}
	//HBondNetStructOP p(output_net_vec_.back());
	//output_net_vec_.pop_back();
//	place_rots_on_pose( *out_pose, p->residues, p->is_native );
//	//place_waters_on_pose( *out_pose, p );
//	(*out_pose).update_residue_neighbors();
//	
//    //TODO MOVE TO output_networks()
//    write_files( *out_pose, p );

    std::vector< Size > net_ids( output_vector_.back() ); // back returs ref to back
    output_vector_.pop_back(); // pop_back() deletes/removes last, does not return;
    
    std::string comment_str( print_headers() + this->print_additional_headers() );
    for (unsigned long & net_id : net_ids){
        //TR << "net_id = " << *net_id << std::endl;
        runtime_assert( get_network_by_id(net_id) != nullptr );
        place_rots_on_pose( *out_pose, (get_network_by_id( net_id ))->residues, (get_network_by_id( net_id ))->is_native );
        //place_waters_on_pose( *out_pose, p );
        //update neighbors
        //write files
        
        add_reslabels_to_pose( *out_pose, *(get_network_by_id( net_id )) );
        comment_str = comment_str + "\n" + (get_network_by_id( net_id ))->outstring;
        
        core::scoring::constraints::ConstraintSetOP cst_op( out_pose->constraint_set()->clone() ); // get existing csts and add to them
        set_constraints( *out_pose, *cst_op, get_network_by_id( net_id ), write_cst_files_ );
        out_pose->constraint_set(cst_op); // add constraints to the pose
        
        //write_files( *out_pose, get_network_by_id( *net_id ) );
    }
    core::pose::add_comment( *out_pose, "HBNet Design details: ", "\n" + comment_str + "\n" );
    
    (*out_pose).update_residue_neighbors();
	( *scorefxn_ )( *out_pose );
	return out_pose;
}

void
HBNet::setup( Pose & pose )
{
	if ( task_factory_ == nullptr ) {
		task_factory_ = task::TaskFactoryOP( new task::TaskFactory );
	}
	task_ = create_ptask( pose ); //set task (which residues are designable/repackable
	//dangerous, if empty default is to start at every designable/packable position in the pose
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
HBNet::search_IG_for_networks( Pose & )
{
    // traverse IG and enumerate all possible h-bond networks given parameters
    if ( ig_ != nullptr ){
        if ( TR.visible() ) {
            TR << " ==================================================================" << std::endl;
            TR << " ============     PERFORMING H-BOND NETWORK DESIGN     ============" << std::endl;
            TR << " ==================================================================" << std::endl;
            traverse_IG( hydrogen_bond_threshold_ );
        }
    }
}

void
HBNet::prepare_output()
{
	output_networks(true); //will add single networks to output vector
}

//NEED TO MOVE THIS TO IT'S OWN CLASS! //TODO
void
HBNet::benchmark_with_native( core::pose::Pose & pose ) // I changed this to & pose, may need to copy Pose
{
	if ( TR.visible() ) TR << "BENCHMARKING AGAINST NATIVE NETWORKS: " << std::endl;
	std::sort( network_vector_.begin(), network_vector_.end(), compare_net_vec() ); //sort all networks

	//for benchmarking/analysis:
	Size count = 1;
	Size ten = 10;
	Size twenty = 20;
	Size fifty = 50;
	Size total_rots = 0;
	Size total_correct_seq = 0;
	Size total_correct_rot = 0;
	Size native_seq = 0;
	Size native_rot = 0;
	Real seq_in_top_ten = 0;
	Real seq_rec_top_ten = 0;
	Real seq_rec_top_twenty = 0;
	Real seq_rec_top_fifty = 0;
	Real rot_in_top_ten = 0;
	Real rot_rec_top_ten = 0;
	Real rot_rec_top_twenty = 0;
	Real rot_rec_top_fifty = 0;
	Real seq_in_top_fiftyperc = 0.0;
	Real seq_in_top_twentyfiveperc = 0.0;
	Real rot_in_top_fiftyperc = 0.0;
	Real rot_in_top_twentyfiveperc = 0.0;
	Real fifty_perc = 0.5*(network_vector_.size());
	Real twentyfive_perc = 0.25*(network_vector_.size());
	Size int_fifty_perc = static_cast<Size>(utility::round(fifty_perc));
	Size int_twentyfive_perc = static_cast<Size>(utility::round(twentyfive_perc));

	for (auto & i : network_vector_) {
		//bool is_native_seq(0);
		//bool is_native_rot(0);
		if ( benchmark_ ) {
			Real  net_size = i->residues.size();
			if ( (this->ligand()) ) { //not fair to count ligand
				for ( utility::vector1< HBondResStructCOP >::const_iterator listit = i->residues.begin(); listit != i->residues.end(); ++listit ) {
					if ( (*listit)->resnum == (this->ligand()) ) {
						net_size = net_size -1;
					}
				}
			}
			pose::PoseOP ala_copy = pose::PoseOP( new pose::Pose( *ala_pose_) ); // make copy of polyalanine pose
			place_rots_on_pose( *ala_copy, i->residues, i->is_native );
			//place_waters_on_pose(*ala_copy, *i);
			ala_copy->update_residue_neighbors();

			Size num_native_seq(get_num_native_seq( *ala_copy, i->residues ) );
			Size num_native_rot(get_num_native_rot( *ala_copy, i->residues, SC_RMSD_CUTOFF, 1 ) );
			Real native_seq_rec = num_native_seq;
			native_seq_rec = native_seq_rec/net_size;
			Real native_rot_rec = num_native_rot;
			native_rot_rec = native_rot_rec/net_size;

			//            if (native_seq_rec==1)
			//                is_native_seq=1;
			//            if (native_rot_rec==1)
			//                is_native_rot=1;

			std::string network( ( pdb_numbering() ) ? ( print_list_to_string( get_orig_pose(), i->residues ) ) : (print_list_to_string( i->residues ) ) );
			std::string outstring = "NETWORK "+utility::to_string(count)+": " + network + " has "+utility::to_string(100*native_seq_rec)+"% native sequence, and has "+utility::to_string(100*native_rot_rec)+"% native rotamers\n                       ";
			if ( TR.visible() ) TR << outstring;

			if ( native_seq_rec == 1 ) {
				native_seq++;
			}
			if ( native_rot_rec == 1 ) {
				native_rot++;
			}
			total_rots = total_rots + net_size;
			total_correct_seq = total_correct_seq + num_native_seq;
			total_correct_rot = total_correct_rot + num_native_rot;
			if ( count == ten ) {
				seq_in_top_ten = 100.0*(native_seq/10.0);
				rot_in_top_ten = 100.0*(native_rot/10.0);
				seq_rec_top_ten = total_correct_seq;
				seq_rec_top_ten = 100.0*(seq_rec_top_ten/total_rots);
				rot_rec_top_ten = total_correct_rot;
				rot_rec_top_ten = 100.0*(rot_rec_top_ten/total_rots);
			} else if ( count == twenty ) {
				seq_rec_top_twenty = total_correct_seq;
				seq_rec_top_twenty = 100.0*(seq_rec_top_twenty/total_rots);
				rot_rec_top_twenty = total_correct_rot;
				rot_rec_top_twenty = 100.0*(rot_rec_top_twenty/total_rots);
			} else if ( count == fifty ) {
				seq_rec_top_fifty = total_correct_seq;
				seq_rec_top_fifty = 100.0*(seq_rec_top_fifty/total_rots);
				rot_rec_top_fifty = total_correct_rot;
				rot_rec_top_fifty = 100.0*(rot_rec_top_fifty/total_rots);
			}
			if ( count == int_fifty_perc ) {
				seq_in_top_fiftyperc = 100.0*(native_seq/(static_cast<Real>(int_fifty_perc)));
				rot_in_top_fiftyperc = 100.0*(native_rot/(static_cast<Real>(int_fifty_perc)));
			} else if ( count == int_twentyfive_perc ) {
				seq_in_top_twentyfiveperc = 100.0*(native_seq/(static_cast<Real>(int_twentyfive_perc)));
				rot_in_top_twentyfiveperc = 100.0*(native_rot/(static_cast<Real>(int_twentyfive_perc)));
			}
		}
		count++;
	}
	// For analysis, writes out to TR benchmark statistics and details:

	Real native_seq_perc = 100.0*((static_cast< core::Real >(native_seq))/network_vector_.size());
	Real native_rot_perc = 100.0*((static_cast< core::Real >(native_rot))/network_vector_.size());
	Real seq_rec_tot = total_correct_seq;
	seq_rec_tot = 100.0*(seq_rec_tot/total_rots);
	Real rot_rec_tot = total_correct_rot;
	rot_rec_tot = 100.0*(rot_rec_tot/total_rots);
	if ( network_vector_.size() < ten ) {
		seq_rec_top_ten = seq_rec_tot;
		rot_rec_top_ten = rot_rec_tot;
	}
	if ( network_vector_.size() < twenty ) {
		seq_rec_top_twenty = seq_rec_tot;
		rot_rec_top_twenty = rot_rec_tot;
	}
	if ( network_vector_.size() < fifty ) {
		seq_rec_top_fifty = seq_rec_tot;
		rot_rec_top_fifty = rot_rec_tot;
	}
	if ( TR.visible() ) {
		TR << " ===============================================================" << std::endl;
		TR << " Sequence recovery = " << seq_rec_top_ten << " % top 10, " << seq_rec_top_twenty << " % top 20, " << seq_rec_top_fifty << " % top 50, " << seq_rec_tot << " % TOTAL, " << std::endl;
		TR << " Rotamer recovery = " << rot_rec_top_ten << " % top 10, " << rot_rec_top_twenty << " % top 20, " << rot_rec_top_fifty << " % top 50, " << rot_rec_tot << " % TOTAL, " << std::endl;
		TR << " NUMBER OF NETWORKS THAT HAVE 100% NATIVE SEQUENCE = " << native_seq << "; " << native_seq_perc << " % total;  ";
		TR << seq_in_top_twentyfiveperc << " % of top 25%, " << seq_in_top_fiftyperc << " % of top 50%, " << seq_in_top_ten << " % of top 10 networks, " << std::endl;
		TR << " NUMBER OF NETWORKS THAT HAVE 100% NATIVE ROTAMERS = " << native_rot << "; " << native_rot_perc << " % total;  ";
		TR << rot_in_top_twentyfiveperc << " % of top 25%, " << rot_in_top_fiftyperc << " % of top 50%, " << rot_in_top_ten << " % top 10, " << std::endl;
	}
	//Real native_perc(0.0);
	native_seq_perc = 0.0;
	store_subnetworks_ = false;

	std::vector< HBondNetStructOP > orig_net_vec = network_vector_;
	network_vector_.clear();
	native_networks_.clear();
	//bool orig_design = design_;
	//bool orig_only_perfect_pairs = only_perfect_pairs_;
	Size orig_max_unsat = max_unsat_;
	Size orig_min_network_size = min_network_size_;
	Size orig_max_network_size = max_network_size_;
	//Real orig_upper_limit = upper_score_limit_;
	bool orig_min = minimize_;

	//ensure pose is scored
	pose.update_residue_neighbors();
	(*init_scorefxn_)(pose);

	//find_native_networks( pose );

	max_unsat_ = orig_max_unsat;
	min_network_size_ = orig_min_network_size;
	max_network_size_ = orig_max_network_size;
	//upper_score_limit_ = orig_upper_limit;
	minimize_ = orig_min;

	native_networks_ = network_vector_;
	network_vector_.clear();
	network_vector_ = orig_net_vec;
	Size lig((this->ligand()) ? 1 : 0);
	Size num_native_nets(native_networks_.size());
	Size found_native(0);
	Size found_native_sequence_50(0);
	Size found_native_sequence_66(0);
	Size found_native_sequence_75(0);
	Size found_native_w_subnetwork_50(0);
	Size found_native_w_subnetwork_66(0);
	Size found_native_w_subnetwork_75(0);
	for (auto & native_network : native_networks_) {
		bool found_nat_sub_50 = false;
		bool found_nat_seq_50 = false;
		bool found_nat_sub_66 = false;
		bool found_nat_seq_66 = false;
		bool found_nat_sub_75 = false;
		bool found_nat_seq_75 = false;
		for (auto & netvec_it : network_vector_) {
			//if ( networks_identical_aa_sequence( **native_it, **netvec_it, true ) ){
			if ( residues_identical( native_network->residues, netvec_it->residues ) ) {
				found_nat_seq_50 = true;
				found_nat_seq_66 = true;
				found_nat_seq_75 = true;
				Size num_native_rot = get_num_native_rot( pose, netvec_it->residues );
				if ( num_native_rot == (netvec_it->residues.size()) ) {
					found_native++;
					found_nat_sub_50 = true;
					found_nat_sub_66 = true;
					found_nat_sub_75 = true;
					break;
				} else if ( 1.0*num_native_rot/(netvec_it->residues.size() - lig) >= 0.5 ) {
					found_nat_sub_50 = true;
					if ( 1.0*num_native_rot/(netvec_it->residues.size() - lig) >= 0.66 ) {
						found_nat_sub_66 = true;
					}
					if ( 1.0*num_native_rot/(netvec_it->residues.size() - lig) >= 0.75 ) {
						found_nat_sub_75 = true;
					}
					break;
				}
			} else if ( is_sub_residues( native_network->residues, netvec_it->residues ) && ( netvec_it->residues.size() - lig) >= 0.5*(native_network->residues.size() - lig ) ) {
				found_nat_seq_50 = true;
				if ( 1.0*( netvec_it->residues.size() - lig) >= 0.66*(native_network->residues.size() - lig ) ) {
					found_nat_seq_66 = true;
				}
				if ( 1.0*( netvec_it->residues.size() - lig) >= 0.75*(native_network->residues.size() - lig ) ) {
					found_nat_seq_75 = true;
				}
				Size num_native_rot = get_num_native_rot( pose, netvec_it->residues );
				if ( 1.0*num_native_rot/(netvec_it->residues.size() - lig) >= 0.5 ) {
					found_nat_sub_50 = true;
				}
				if ( 1.0*num_native_rot/(netvec_it->residues.size() - lig) >= 0.66 ) {
					found_nat_sub_66 = true;
				}
				if ( 1.0*num_native_rot/(netvec_it->residues.size() - lig) >= 0.75 ) {
					found_nat_sub_75 = true;
				}
				break;
			}
		}
		if ( found_nat_sub_50 ) {
			found_native_w_subnetwork_50++;
		}
		if ( found_nat_seq_50 ) {
			found_native_sequence_50++;
		}
		if ( found_nat_sub_66 ) {
			found_native_w_subnetwork_66++;
		}
		if ( found_nat_seq_66 ) {
			found_native_sequence_66++;
		}
		if ( found_nat_sub_75 ) {
			found_native_w_subnetwork_75++;
		}
		if ( found_nat_seq_75 ) {
			found_native_sequence_75++;
		}
	}
	//native_perc = 1.0*found_native/num_native_nets;
	//native_seq_perc = 1.0*found_native_sequence/num_native_nets;
	if ( TR.visible() ) {
		TR << " ===============================================================" << std::endl;
		TR << " WHAT PERCENTAGE OF THE NATIVE NETWORKS WERE RECAPITULATED?" << std::endl;
		TR << " There are " << num_native_nets << " native networks in the pose." << std::endl;
		TR << " Recapitulated " << found_native << " of the complete native networks;" << std::endl;
		TR << "   " << found_native_w_subnetwork_50 << " networks that are a subset with >= 50% rotamers identical to the native network" << std::endl;
		TR << "   " << found_native_w_subnetwork_66 << " networks that are a subset with >= 66.6% rotamers identical to the native network" << std::endl;
		TR << "   " << found_native_w_subnetwork_75 << " networks that are a subset with >= 75% rotamers identical to the native network" << std::endl;
		TR << "   " << found_native_sequence_50 << " networks that are a subset with >= 50% sequence identical to the native network" << std::endl;
		TR << "   " << found_native_sequence_66 << " networks that are a subset with >= 66.6% sequence identical to the native network" << std::endl;
		TR << "   " << found_native_sequence_75 << " networks that are a subset with >= 75% sequence identical to the native network" << std::endl;
	}
}

Size
HBNet::num_core_res( hbond_net_struct const & network ){
	Size num_core(0);
	for (const auto & residue : network.residues) {
		if ( res_is_core( residue->resnum ) ) ++num_core;
	}
	return num_core;
}

Size
HBNet::num_boundary_res( hbond_net_struct const & network ){
	Size num_boundary(0);
	for (const auto & residue : network.residues) {
		if ( res_is_boundary( residue->resnum ) ) ++num_boundary;
	}
	return num_boundary;
}

void
HBNet::select_best_networks()
{
	auto i = network_vector_.begin();
	for ( ; i != network_vector_.end(); ) {
		if ( (*i)->residues.size() < 3 && min_network_size_ > 2 ) {
			//if (verbose_) TR << "i = " << print_list_to_string( (*i)->residues ) << "ERASING! Size = " << (*i)->residues.size() << std::endl;
			i = network_vector_.erase( i );
		} else if ( min_core_res_ && num_core_res( **i ) < min_core_res_ ) {
			i = network_vector_.erase( i );
		} else if ( min_boundary_res_ && num_boundary_res( **i ) < min_boundary_res_ ) {
			i = network_vector_.erase( i );
		} else {
			++i;
		}
	}

	//iterator over all the networks, erase those that do not meet specified criteria
	//std::vector< HBondNetStructOP >::iterator i = network_vector_.begin();
	i = network_vector_.begin();
	for ( ; i != network_vector_.end(); ) {
		//place network residues on background pose for scoring and evaluation
		Pose ala_copy = *ala_pose_;
		place_rots_on_pose( ala_copy, (*i)->residues, (*i)->is_native );
		ala_copy.update_residue_neighbors();
		(*init_scorefxn_)(ala_copy); //need this?

		//ADD FUNCTION HERE TO CHECK IF NETWORK CONTAINS SPECIFICIED AA TYPES
		//make an asymmetric representation of symmetric networks for proper scoring
		Size net_size = (*i)->residues.size();

		if ( net_size < min_network_size_ || net_size > max_network_size_ ) {
			if ( verbose_ && TR.visible() ) TR << "i = " << print_list_to_string( (*i)->residues ) << "ERASING! Size = " << net_size << std::endl;
			i = network_vector_.erase( i );
		} else {
			get_hbond_atom_pairs( **i, ala_copy );
			(*i)->total_hbonds = (*i)->hbond_vec.size();
			if ( (*i)->hbond_vec.empty() ) {
				if ( verbose_ ) TR << "i = " << print_list_to_string( (*i)->residues ) << "ERASING! h-bonds empty" << std::endl;
				i = network_vector_.erase( i );
			} else {
				//find all buried unsatisfied polars in the network
				//find_unsats( ala_copy, *i, new_packer_graph );

				find_unsats( ala_copy, **i );

				if ( symmetric_ && !multi_component_ && !( orig_pose_->total_residue() > 500 ) ) {
					bool network_spans_entire_symm_interface(false);
					for ( utility::vector1< HBondCOP >::const_iterator h = (*i)->hbond_vec.begin(); h != (*i)->hbond_vec.end(); ++h ) {
						//(*i)->show(in_pose, 1, TR);
						Size drsd((*h)->don_res());
						//Size don_hatm((*h)->don_hatm());
						Size arsd((*h)->acc_res());
						//Size aatm((*h)->acc_atm());
						if ( arsd == get_ind_res( *orig_pose_, drsd) || drsd == get_ind_res( *orig_pose_, arsd) ) {
							network_spans_entire_symm_interface = true;
							break;
						}
					}
					//utility::vector1< HBondResStructCOP > temp_residues( (*i)->residues );
					if ( network_spans_entire_symm_interface ) {
						for ( utility::vector1< HBondResStructCOP >::const_iterator ir = (*i)->residues.begin(); ir != (*i)->residues.end(); ++ir ) {
							if ( find_hbond_res_struct( (*i)->asymm_residues, (*ir)->resnum ) == (*i)->asymm_residues.end() ) {
								(*i)->asymm_residues.push_back( *ir );
							}
							utility::vector1<Size> resi_clones( get_symm_info()->bb_clones( (*ir)->resnum ) );
							for ( utility::vector1<Size>::const_iterator r = resi_clones.begin(); r != resi_clones.end(); ++r ) {
								if ( find_hbond_res_struct( (*i)->residues, *r ) == (*i)->residues.end() ) {
									(*i)->asymm_residues.push_back( HBondResStructCOP( new hbond_res_struct( *r, 0, (*ir)->aa, orig_pose_->pdb_info()->chain(*r), orig_pose_->residue(*r).is_protein(), 0, orig_pose_->residue(*r).is_ligand() ) ) );
								}
							}
						}
					}
				}
				//move some of the to meets_criteria NEED TO FIX
				if ( ( no_heavy_unsats_allowed_ && (*i)->num_heavy_unsat > 0 ) || (*i)->num_unsat > max_unsat_ || ( this->ligand() && (*i)->lig_num_unsatisfied > max_lig_unsat_ ) ||
						!(this->network_meets_criteria( ala_copy, **i )) ) {
					if ( verbose_ && TR.visible() ) TR << "i = " << print_list_to_string( (*i)->residues ) << "ERASING! num unsat = " << (*i)->num_unsat << "; num_heavy_unsat = " << (*i)->num_heavy_unsat << std::endl;
					i = network_vector_.erase( i );
				} else {
					++i;
				}
			}
		}
	}
	//score the networks
	if ( verbose_ && TR.visible() ) {
		TR << "NUMBER OF NETWORKS BEFORE SCORE = " << network_vector_.size() << std::endl;
		TR << " SCORING THE NETWORKS: " << std::endl;
	}
	score_networks( minimize_ );
    
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
    if ( extend_existing_networks_ ){
        //check that native_networks_ is not empty or nullptr
        //for each rotamer in native_networks_: remove them from rotamer_sets_
        // if keep_existing_networks_ true also, then will replace existing network with extended network IF EXTENDED IS BETTER
        if ( only_extend_existing_ ) start_res_vec_.clear();
        for ( Size i = 1; i <= hbnet_info_residues_.size(); ++i ) {
            if ( hbnet_info_residues_[i] && rotamer_sets_->rotamer_set_for_residue( i ) != nullptr ){
                //TR.Debug << "dropping HBNet rotamer, i = " << i << std::endl;
                // Removee all rotamers except input rotamer
                // if extend_existing_networks_ = true then input rotamer will be retained: or_include_current(true) added to RLT in setup
                utility::vector1<bool> rots_to_remove( rotamer_sets_->rotamer_set_for_residue( i )->num_rotamers(), true);
                rotamer_sets_->rotamer_set_for_residue( i )->drop_rotamers( rots_to_remove );
                start_res_vec_.insert( i );
            }
        }
        for (auto & native_network : native_networks_){
            for (auto & residue : native_network->residues){
                if ( rotamer_sets_->rotamer_set_for_residue( residue->resnum ) != nullptr ){
                    //TR.Debug << "dropping native_network rotamer, (*net_it)->resnum = " << (*net_it)->resnum << std::endl;
                    // Removee all rotamers except input rotamer
                    // if extend_existing_networks_ = true then input rotamer will be retained: or_include_current(true) added to RLT in setup
                    utility::vector1<bool> rots_to_remove( rotamer_sets_->rotamer_set_for_residue( residue->resnum )->num_rotamers(), true);
                    rotamer_sets_->rotamer_set_for_residue( residue->resnum )->drop_rotamers( rots_to_remove );
                    start_res_vec_.insert( residue->resnum );
                }
            }
        }
    }
    this->trim_additional_rotamers( pose );
}

void
HBNet::run( Pose & pose )
{
    if ( find_native_ ) {
        std::set< core::Size > temp_start_res_vec_ = start_res_vec_; // save starting residues
        if ( start_res_vec_.empty() ) {
            if ( TR.visible() ) TR << "No starting residues defined; traversing the entire static pose:" << std::endl;
            Size total( (symmetric_) ? symm_info_->num_independent_residues() : pose.total_residue() );
            for ( Size res = 1; res <= total; ++res ) {
                start_res_vec_.insert(res);
            }
        }
        utility::vector1<Size> residues_to_ala(0);
        for ( Size r = 1; r <= pose.total_residue(); ++r ) {
            if ( pose.residue(r).is_ligand() || pose.residue_type(r).name() == "VRT" ) {
                continue;
            }
            residues_to_ala.push_back(r);
        }
        ala_pose_ = PoseOP( new Pose(pose) );
        protocols::toolbox::pose_manipulation::construct_poly_XXX_pose("ALA", *ala_pose_, residues_to_ala, 1, 1, 1);
        ala_pose_->update_residue_neighbors();
        (*scorefxn_)(*ala_pose_); // score now so we don't have to later; several functions require pose to be scored
        
        traverse_native( pose, hydrogen_bond_threshold_ );
        
        bool orig_only_native = only_native_;
        only_native_ = true; // need better solution than this
        branch_overlapping_networks();
        remove_replicate_networks( max_rep_ );
        // if want to find natives to with criteria specified in HBNet mover defition (independent of input HBNet PDBInfoLabels)
        if ( find_native_ ) select_best_networks();
        only_native_ = orig_only_native;
        
        start_res_vec_ = temp_start_res_vec_;
        native_networks_ = network_vector_;
        network_vector_.clear();
        
        if ( only_native_ )
            return;
    }

    utility::vector1<Size> residues_to_ala(0);
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
    ala_pose_ = pose::PoseOP( new Pose(pose) );
    protocols::toolbox::pose_manipulation::construct_poly_XXX_pose("ALA", *ala_pose_, residues_to_ala, 1, 1, 1);
    ala_pose_->update_residue_neighbors();
    (*scorefxn_)(*ala_pose_); // score now so we don't have to later; several functions require pose to be scored
    
    core::pack::task::operation::RestrictAbsentCanonicalAASOP set_hbnet_des_res = core::pack::task::operation::RestrictAbsentCanonicalAASOP( new core::pack::task::operation::RestrictAbsentCanonicalAAS() );
    set_hbnet_des_res->keep_aas( des_residues_ );
    set_hbnet_des_res->include_residue(0);
    set_hbnet_des_res->apply( pose, *task_ );
    
    //HBNet requires PDInteractionGraph (or DensePDInteractionGraph).  Need to check to make sure NOT using lin_mem_ig
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
        rotamer_set::symmetry::SymmetricRotamerSetsOP sym_rotset_op( new rotamer_set::symmetry::SymmetricRotamerSets() );
        rotamer_sets_ = sym_rotset_op;
    } else {
        rotamer_sets_ = RotamerSetsOP(new rotamer_set::RotamerSets());
    }
    //need mechanism to keep/extend networks by HBNet InfoLabel; else detect networks based on criteria
    if ( keep_existing_networks_ || extend_existing_networks_ ){
        //for every position with HBNet info label in input pdb
        for ( Size i = 1; i <= hbnet_info_residues_.size(); ++i ) {
            if ( hbnet_info_residues_[i] ){
                if ( keep_existing_networks_ && !extend_existing_networks_ )
                    task_->nonconst_residue_task( i ).prevent_repacking();
                else if ( extend_existing_networks_ ){
                    // NEED TO CHECK: do any of these flip it to packable if task_factory alraedy designates it as prevent_repacking()?
                    task_->nonconst_residue_task( i ).restrict_to_repacking(); //ensure packable
                    //task_->nonconst_residue_task( i ).allow_aa(pose.residue(i).aa()); //ensure packable
                    task_->nonconst_residue_task( i ).or_include_current( true );
                    task_->nonconst_residue_task( i).sample_proton_chi(true);
                }
            }
        }
        //for each rotamer in native_networks_: set them to not be packable/designable in task_
        for (auto & native_network : native_networks_){
            for (auto & residue : native_network->residues){
                if ( keep_existing_networks_ && !extend_existing_networks_ )
                    task_->nonconst_residue_task( residue->resnum ).prevent_repacking();
                else if ( extend_existing_networks_ ){
                    // do any of these flip it to packable if task_factory alraedy designates it as prevent_repacking()?
                    task_->nonconst_residue_task( residue->resnum ).restrict_to_repacking(); //ensure packable
                    //task_->nonconst_residue_task( (*net_it)->resnum ).allow_aa(pose.residue((*net_it)->resnum).aa()); //ensure packable
                    task_->nonconst_residue_task( residue->resnum ).or_include_current( true );
                    task_->nonconst_residue_task( residue->resnum ).sample_proton_chi(true);
                }
            }
        }
    }
    pose.update_residue_neighbors();
    (*init_scorefxn_).setup_for_packing( pose, task_->repacking_residues(), task_->designing_residues() );
    packer_neighbor_graph_ = create_packer_graph( pose, *init_scorefxn_, task_ ); //even in symmetric cases, packer graph will have a node for every residue in the pose
    
    rotamer_sets_->set_task( task_ );
    rotamer_sets_->initialize_pose_for_rotsets_creation(pose);
    rotamer_sets_->build_rotamers( pose, *init_scorefxn_, packer_neighbor_graph_ );
    
    //remove unwanted rotamers from rotamer_sets_ before IG is populated and searched.
    trim_rotamers( pose );
    
    // NEED TO MOVE TO HBNetLigand
    //    if ( start_from_csts_ ){
    //        for (std::set< Size >::const_iterator resvecit = start_res_vec_.begin();  resvecit != start_res_vec_.end(); ++resvecit){
    //            //TR << "resvecit = " << *resvecit << "; id_for_current_rotamer = " << rotamer_sets_->rotamer_set_for_moltenresidue(rotamer_sets_->resid_2_moltenres((platform::uint)(*resvecit)))->id_for_current_rotamer() << std::endl;
    //            rotamer_sets_->clear_rotamer_set_for_residue( *resvecit, !(task_->include_current( *resvecit )));
    //        }
    //    }
    
    rotamer_sets_->prepare_sets_for_packing( pose, *init_scorefxn_ );
    ig_ = InteractionGraphFactory::create_interaction_graph( *task_, *rotamer_sets_, pose, *init_scorefxn_ );
    ig_->initialize( *rotamer_sets_ );
    
    PrecomputedPairEnergiesInteractionGraphOP pig( utility::pointer::dynamic_pointer_cast< PrecomputedPairEnergiesInteractionGraph > ( ig_ ) );
    runtime_assert(pig);
    
    if ( symmetric_ ) {
        for ( platform::uint ii = 1; ii <= rotamer_sets_->nmoltenres(); ++ii ) {
            utility::vector1< core::PackerEnergy > one_body_energies((rotamer_sets_->rotamer_set_for_moltenresidue(ii))->num_rotamers() );
            hbnet_symm_one_body_energies(pose, rotamer_sets_->rotamer_set_for_moltenresidue(ii), *init_scorefxn_, *task_, packer_neighbor_graph_, one_body_energies );
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
        for (unsigned long mrvit : molten_resvec) {
            platform::uint resvecit_uint(*resvecit);
            if ( resvecit_uint == mrvit ) {
                in_design_shell = true;
                break;
            }
        }
        if ( !( in_design_shell ) ) {
            if ( TR.visible() ) TR << " WARNING: residue " << *resvecit << " is not in the IG; I'm removing it from start_resnums .  If you wish to consider this residue, change your design shell task operations" << std::endl;
            //does NOT return iterator referencing same location after removal like std::vector does!
            start_res_vec_.erase( resvecit++ );
        } else {
            ++resvecit;
        }
    }
    runtime_assert( !(start_res_vec_.empty()));
    
    this->search_IG_for_networks( pose ); // can be overriden for special cases of handling and searching IG
    
    if ( TR.visible() ) {
        TR << " INITIAL NUMBER OF H-BOND NETWORKS FOUND: " << network_vector_.size() << std::endl;
        TR << " BRANCHING NETWORKS TOGETHER TO FORM COMPLETE NETWORKS: " << std::endl;
    }
    branch_overlapping_networks();
    if ( TR.visible() ) TR << " NUMBER OF H-BOND NETWORKS AFTER BRANCH: " << network_vector_.size() << std::endl;
    remove_replicate_networks( max_rep_ );
    if ( TR.visible() ) TR << "NUMBER OF NETWORKS AFTER REMOVE_REP = " << network_vector_.size() << std::endl;
    
    select_best_networks();
}//run
    
void
HBNet::apply( Pose & pose )
{
	if ( !( pose.pdb_info() ) ) {
		pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( pose, true ) ) );
	}
	//NEED TO FIX THIS, IS THIS NEEDED IF MOVER CALLED MULTIPLE TIMES?
	network_vector_.clear();
	rotamer_sets_.reset(); // resets OP to 0
	ig_.reset(); // resets OP to 0

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

	if ( scorefxn_ == nullptr ) {
		scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "talaris2013_cst" );
	}

	pose.update_residue_neighbors();
	(*init_scorefxn_)(pose); //need this?

	HBondSet bg_hbond_set;
	bg_hbond_set.setup_for_residue_pair_energies( pose, false/*calculate_derivative*/, true/*backbone_only*/ );
	//NEED TO FIX THIS!!!

	hb_database_ = HBondDatabase::get_database( bg_hbond_set.hbond_options().params_database_tag());

	orig_pose_ = PoseOP( new Pose(pose) ); //pose::PoseOP pointing to original pose
	//output_pose_ = pose.get_self_ptr(); //careful need to reset this pointer before apply() returns otherwise will create problems...

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
		core::select::residue_selector::ResidueSubset rs( start_selector_->apply( pose ) );
		start_res_vec_.clear();
		for ( Size r = 1; r <= rs.size(); ++r ) {
			if ( rs[ r ] ) {
				start_res_vec_.insert(r);
			}
		}
	}
	if ( core_selector_ ) {
		core_residues_ = core_selector_->apply( pose );
	} else {
		core::select::residue_selector::LayerSelectorOP core_layer( new core::select::residue_selector::LayerSelector() );
		core_layer->set_layers( true, false, false );
		//core_layer->set_use_sc_neighbors( true ); // now true by default
        // TODO need to set default to be more lenient, 3.6 or 4
		core_residues_ = core_layer->apply( pose );
	}
	if ( boundary_selector_ ) {
		boundary_residues_ = boundary_selector_->apply( pose );
	} else {
		core::select::residue_selector::LayerSelectorOP boundary_layer( new core::select::residue_selector::LayerSelector() );
		boundary_layer->set_layers( false, true, false );
		//boundary_layer->set_use_sc_neighbors( true );
		boundary_residues_ = boundary_layer->apply( pose );
	}
    if ( keep_existing_networks_ || extend_existing_networks_ ){
        //Need to check: are PDBInfoLabels kept after pose made symmetric? what about after BGS?
        core::select::residue_selector::ResiduePDBInfoHasLabelSelectorOP hbnet_info_label( new core::select::residue_selector::ResiduePDBInfoHasLabelSelector( "HBnet" ) );
        hbnet_info_residues_ = hbnet_info_label->apply( pose );
    }
	this->setup( pose );
    //run method common to all HBNet runs
    run( pose );
	TR.flush();
	this->prepare_output();
    if ( only_native_ ){
        if ( native_networks_.size() == 0 ) set_last_move_status( protocols::moves::FAIL_RETRY );
        else set_last_move_status( protocols::moves::MS_SUCCESS );
        return;
    }
	TR.Debug << "After prepare_output, network_vector_.size() = " << network_vector_.size() << std::endl;
	//TR.Debug << "output_net_vec_.size() = " << output_net_vec_.size() << std::endl;
	//if ( output_net_vec_.size() == 0 ) {
    if ( output_vector_.size() == 0 ) {
		if ( TR.visible() ) TR << "DID NOT FIND SOLUTIONS THAT MEET YOUR CRITERIA! EXITING" << std::endl;
		set_last_move_status( protocols::moves::FAIL_RETRY );
		ig_.reset(); // resets OP to 0
		packer_neighbor_graph_.reset();
		return;
	} else {
		set_last_move_status( protocols::moves::MS_SUCCESS );
	}
	if ( benchmark_ ) {
		benchmark_with_native( pose );
	}
    ////this shouldn't be necessary but is a safeguard
	//std::sort( output_net_vec_.begin(), output_net_vec_.end(), compare_net_vec() );
	//reverse the vector so when we pop_back() the best is at the back
	//std::reverse( output_net_vec_.begin(), output_net_vec_.end());
    
    // need to reverse because no pop_front for std::vector -- could use list but already set up for sorting functions with vectors of HBondNetStructOP's
    std::reverse( output_vector_.begin(), output_vector_.end());

    std::vector< Size > net_ids( output_vector_.back() );
    output_vector_.pop_back(); // pop_back() removes without returning
    
    std::string comment_str( print_headers() + this->print_additional_headers() );
    for (unsigned long & net_id : net_ids){
        //TR << "net_id = " << *net_id << std::endl;
        runtime_assert( get_network_by_id(net_id) != nullptr );
        place_rots_on_pose( pose, (get_network_by_id( net_id ))->residues, (get_network_by_id( net_id ))->is_native );
        //place_waters_on_pose( *out_pose, p );
        //update neighbors
        //write files
        
        add_reslabels_to_pose( pose, *(get_network_by_id( net_id )) );
        comment_str = comment_str + "\n" + (get_network_by_id( net_id ))->outstring;
        
        core::scoring::constraints::ConstraintSetOP cst_op( pose.constraint_set()->clone() ); // get existing csts and add to them
        set_constraints( pose, *cst_op, get_network_by_id( net_id ), write_cst_files_ );
        pose.constraint_set(cst_op); // add constraints to the pose
    }
    core::pose::add_comment( pose, "HBNet Design details: ", "\n" + comment_str + "\n" );
    
    pose.update_residue_neighbors();
    ( *scorefxn_ )( pose );
    
//	HBondNetStructOP p(output_net_vec_.back());
//	output_net_vec_.pop_back();
//	place_rots_on_pose( pose, p->residues, p->is_native );
//	//place_waters_on_pose( pose, p );
//	pose.update_residue_neighbors();
//    
//    write_files( pose, p );

	ig_.reset(); // resets OP to 0 / nullptr
	packer_neighbor_graph_.reset();
	if ( TR.visible() ) TR.flush();
	if ( TR.Warning.visible() ) TR.Warning.flush();
	if ( TR.Debug.visible() ) TR.Debug.flush();
	//output_pose_.reset(); //If make PoseOP that directly points to original Pose & pose of apply(), need to reset it!
}// apply

} //hbnet
} //protocols
