// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/SnugDockProtocol.cc
/// @brief Dock and antigen to an antibody while optimizing the rigid body orientation of the VH and VL chains and
/// performing CDR loop minimization.
/// @details
///
///
/// @author Jianqing Xu ( xubest@gmail.com )
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )
/// @author Jeliazko Jeliazkov ( jeliazkov@jhu.edu )

// Unit headers
#include <protocols/antibody/snugdock/SnugDockProtocol.hh>

// Package headers
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyNumberingConverterMover.hh>
#include <protocols/antibody/snugdock/SnugDock.hh>
#include <protocols/antibody/util.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loop_modeling/LoopProtocol.hh>
#include <protocols/loop_modeler/LoopModeler.hh>

// Project headers
#include <protocols/docking/DockingProtocol.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>

#include <core/select/util/interface_vector_calculate.hh>

// constraints

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

// Utility headers
#include <utility/tools/make_vector1.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <set>

#include <utility/stream_util.hh> // AUTO IWYU For operator<<

static basic::Tracer TR( "protocols.antibody.SnugDockProtocol" );

using namespace core;
using namespace core::scoring;

namespace protocols {
namespace antibody {
namespace snugdock {

///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief default constructor
SnugDockProtocol::SnugDockProtocol() : Mover() {
	init();
}

/// @brief copy constructor
SnugDockProtocol::SnugDockProtocol( SnugDockProtocol const & rhs ) : Mover(rhs) {
	init_for_equal_operator_and_copy_constructor( *this, rhs );
}

/// @brief assignment operator
SnugDockProtocol & SnugDockProtocol::operator=( SnugDockProtocol const & rhs ) {
	//abort self-assignment
	if ( this == &rhs ) return *this;
	Mover::operator=( rhs );
	init_for_equal_operator_and_copy_constructor( *this, rhs );
	return *this;
}

//destructor
SnugDockProtocol::~SnugDockProtocol() = default;

/// @brief Each derived class must specify its name.
std::string SnugDockProtocol::get_name() const {
	return type();
}

//@brief clone operator, calls the copy constructor
protocols::moves::MoverOP
SnugDockProtocol::clone() const {
	return utility::pointer::make_shared< SnugDockProtocol >( *this );
}

/// @brief fresh_instance returns a default-constructed object for JD2
protocols::moves::MoverOP
SnugDockProtocol::fresh_instance() const {
	return utility::pointer::make_shared< SnugDockProtocol >();
}

/// @brief This mover retains state such that a fresh version is needed if the input Pose is about to change
bool SnugDockProtocol::reinitialize_for_new_input() const {
	return true;
}

void SnugDockProtocol::register_options() {
	docking::DockingProtocol::register_options();
	SnugDock::register_options();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// END OF BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

void SnugDockProtocol::init() {
	Mover::type( "SnugDockProtocol" );

	set_default();
	init_from_options();

}

void SnugDockProtocol::init_for_equal_operator_and_copy_constructor( SnugDockProtocol & lhs, SnugDockProtocol const & rhs ) {
	// copy all data members from rhs to lhs
	lhs.antibody_info_ = rhs.antibody_info_;

	// Movers
	lhs.low_res_refine_cdr_h2_ = rhs.low_res_refine_cdr_h2_;
	lhs.low_res_refine_cdr_h3_ = rhs.low_res_refine_cdr_h3_;
	lhs.docking_ = rhs.docking_;

}

void SnugDockProtocol::set_default() {
	auto_generate_kink_constraint_ = true;
	high_res_kink_constraint_ = true;
	low_res_kink_constraint_ = true;
	constrain_vlvh_qq_ = true;
	ab_has_light_chain_ = false;
}

void SnugDockProtocol::init_from_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ OptionKeys::antibody::auto_generate_h3_kink_constraint ].user() ) {
		auto_generate_kink_constraint( option[ OptionKeys::antibody::auto_generate_h3_kink_constraint ]() );
	}
	if ( option[ OptionKeys::antibody::h3_loop_csts_hr ].user() ) {
		high_res_kink_constraint( option[ OptionKeys::antibody::h3_loop_csts_hr ]() );
	}
	if ( option[ OptionKeys::antibody::h3_loop_csts_lr ].user() ) {
		low_res_kink_constraint( option[ OptionKeys::antibody::h3_loop_csts_lr ]() );
	}
	if ( option[ OptionKeys::antibody::constrain_vlvh_qq ].user() ) {
		constrain_vlvh_qq( option[ OptionKeys::antibody::constrain_vlvh_qq ]() );
	}
}

void SnugDockProtocol::apply( Pose & pose ) {
	using namespace protocols::analysis;
	using namespace basic::options;

	TR << "Beginning apply function of " + get_name() + "." << std::endl;

	//if ( ! antibody_info_ ) setup_objects( pose );
	setup_objects( pose );

	show( TR );
	// Note: native flag cannot be used because VRTs are prepended!
	TR << "Using universal FoldTree." << std::endl;

	// set constraints, if specified
	if ( auto_generate_kink_constraint() ) {
		antibody::kink_constrain_antibody_H3( pose, CDR_loops_[h3].stop() - 2 );
	}

	TR << "Beginning application of " + docking()->get_name() + "." << std::endl;
	docking()->apply( pose );

	// before re-numbering or comparing or anything, revert fold tree here:
	// delete VRTs and correctly label chains in PDBInfo!
	// Note: AntibodyInfo is reliant on correctly labeled H & L chains.
	// take advantage of VRTs being located first in the "universal" pose/FT
	while ( pose.residue(1).is_virtual_residue() ) {
		pose.delete_residue_slow(1);
	}

	//JAB - Needed as delete_residue_slow clears the energies cache
	ScoreFunctionOP scorefxn = get_score_function();
	scorefxn->score(pose);

	if ( option[ OptionKeys::antibody::output_ab_scheme].user() ) {
		AntibodyNumberingConverterMover converter = AntibodyNumberingConverterMover();
		converter.apply( pose ); // might become broken?
	}

	// move Jared's code here
	if ( antibody_info_->antigen_present() ) {
		TR << "Running Interface AnalyzerMover" << std::endl;
		// use the first jump (should always be true with universal set up
		core::pose::DockingPartners partners = core::pose::DockingPartners::docking_partners_from_string( option[ OptionKeys::docking::partners ]() );
		InterfaceAnalyzerMover analyzer = InterfaceAnalyzerMover( partners /* get jump from partners flag */, false /* tracer */, nullptr /* let it construct sf*/, false /* compute_packstat */ , false /* pack_input */,  true /* pack_separated */) ;

		analyzer.apply(pose); //Adds PoseExtraScore_Float to be output into scorefile.
	}
}

void SnugDockProtocol::setup_objects( Pose & pose ) {
	TR << "Setting up data for " + get_name() + "." << std::endl;

	/// AntibodyInfo is used to store information about the Ab-Ag complex and to generate useful helper objects based on
	/// that information (e.g. the various FoldTrees that are needed for SnugDock).
	antibody_info_ = utility::pointer::make_shared< AntibodyInfo >( pose );

	// Convert FT to "Universal FT", so we don't have to repeatedly change back and forth.
	// This also should remove the necessity for LH_A chain order.
	setup_ab_ag_foldtree( pose, antibody_info_ );
	// from now on use a pose with leading VRTs and a single FT
	// ab/ag_residues are set above
	tf_ = repack_tf_from_residue_sets( pose, ab_residues_, ag_residues_ );

	setup_loop_refinement_movers();
	// re-initialize docking each time, since setup_objects must be called each time to configure the FT
	// this is probably sub-optimal
	docking_=nullptr;

	// do constraint setting things first, as this requires re-syncing flags
	if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_file ].user() ) {
		docking()->set_use_constraints( true );
	}
	if ( constrain_vlvh_qq() ) {
		// hard-coded weight for Q-Q cst... maybe a bad idea?
		docking()->set_cst_weight( 1.0 );
	}


	if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_file ].user() || constrain_vlvh_qq() ) {
		docking()->set_no_filters( true );
		docking()->sync_objects_with_flags(); // should fix lose of snugdock mover
	}

	// set tf/low-res stuff
	docking()->set_task_factory( tf_ );
	docking()->add_additional_low_resolution_step( low_res_refine_cdr_h2_ );
	docking()->add_additional_low_resolution_step( low_res_refine_cdr_h3_ );

	SnugDockOP high_resolution_phase( new SnugDock );
	// if debugging
	//high_resolution_phase->debug();
	// pass info on new FT here!
	high_resolution_phase->set_CDR_loops(CDR_loops_);
	high_resolution_phase->set_ab_has_light_chain(ab_has_light_chain_);
	high_resolution_phase->set_task_factory(tf_);
	if ( ab_has_light_chain_ ) { high_resolution_phase->set_vh_vl_jump( vh_vl_jump_ ); }
	// pass on kink constraint to high-res docking move
	if ( high_res_kink_constraint() ) {
		high_resolution_phase->high_res_kink_constraint( true );
	}
	docking()->set_docking_highres_mover( high_resolution_phase );

}

void SnugDockProtocol::setup_loop_refinement_movers() {
	using core::scoring::ScoreFunctionFactory;
	using core::scoring::ScoreFunctionOP;

	if ( ! antibody_info_ ) {
		using utility::excn::Exception;
		throw CREATE_EXCEPTION(Exception,  "A valid AntibodyInfo instance is required to setup " + get_name() + "'s centroid loop "
			+ "refinement movers." );
	}

	/// FIXME: The chain break weight configuration and constraint weight should be handled by RefineOneCDRLoop.
	ScoreFunctionOP low_res_loop_refinement_scorefxn = ScoreFunctionFactory::create_score_function("cen_std", "score4L");
	low_res_loop_refinement_scorefxn->set_weight( scoring::chainbreak, 1.0 );
	low_res_loop_refinement_scorefxn->set_weight( scoring::overlap_chainbreak, 10./3. );
	low_res_loop_refinement_scorefxn->set_weight( scoring::atom_pair_constraint, 100 );

	// update low-res sfxn with kink constraint, unlike H3 modeling, constriants are not enable in low-res by default
	// also weights are hard coded, so this should be refactored later
	if ( auto_generate_kink_constraint() and low_res_kink_constraint() ) {
		low_res_loop_refinement_scorefxn->set_weight( scoring::dihedral_constraint, 1.0 );
		low_res_loop_refinement_scorefxn->set_weight( scoring::angle_constraint, 1.0 );
	}

	// if csts were specified for DockingProtocol, then pose has csts
	// so just add atom_pair constraints to the sfxn to include them
	if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_file ].user() ) {
		// set weight to user specified or 1.0
		core::Real cst_weight = 1.0;
		if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_weight ].user() ) {
			cst_weight = basic::options::option[ basic::options::OptionKeys::constraints::cst_weight ];
		}
		low_res_loop_refinement_scorefxn->set_weight( scoring::atom_pair_constraint, cst_weight );
	}

	low_res_refine_cdr_h2_ = refine_loop();
	low_res_refine_cdr_h3_ = refine_loop();
	// pass Loop objects from above
	low_res_refine_cdr_h2_->set_loop( CDR_loops_[h2] );
	low_res_refine_cdr_h3_->set_loop( CDR_loops_[h3] );
	// pass scorefunction from above
	low_res_refine_cdr_h2_->set_cen_scorefxn( low_res_loop_refinement_scorefxn );
	low_res_refine_cdr_h3_->set_cen_scorefxn( low_res_loop_refinement_scorefxn );

}

// add function to set up the taskfactor with the new fold tree
core::pack::task::TaskFactoryOP SnugDockProtocol::repack_tf_from_residue_sets(Pose const & pose, std::set< core::Size > const & ab_residues, std::set< core::Size > const & ag_residues) {

	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	// new style creation of a task factory
	// set it based on distance of neighbors
	// use interface_vector_calculate then pass vector
	TaskFactoryOP tf( new TaskFactory() );
	// restrict to repacking & typical settings
	tf->push_back( utility::pointer::make_shared< RestrictToRepacking >() );
	tf->push_back( utility::pointer::make_shared< InitializeFromCommandline >() );
	tf->push_back( utility::pointer::make_shared< IncludeCurrent >() );
	tf->push_back( utility::pointer::make_shared< NoRepackDisulfides >() );
	// check if resfile option is given, then read from it
	if ( basic::options::option[ basic::options::OptionKeys::packing::resfile ].user() ) tf->push_back( utility::pointer::make_shared< core::pack::task::operation::ReadResfile >() );

	// identify interface residues hard coded numbers are taken from function defaults
	utility::vector1<bool> non_interacting_residues = core::select::util::calc_interacting_vector( pose, ab_residues, ag_residues, 10.0, 5.5, 75.0, 9.0 );
	// invert the vector to get non-interactive residues
	non_interacting_residues.flip();

	// Residue-Level TaskOperation to prevent repacking on non-interacting residues
	PreventRepackingRLTOP prevent_repacking = utility::pointer::make_shared< PreventRepackingRLT >();
	OperateOnResidueSubsetOP subset_op = utility::pointer::make_shared< OperateOnResidueSubset >( prevent_repacking, non_interacting_residues);
	tf->push_back(subset_op);

	return tf;

}
// function returns a virtual residue positioned at COM with correct XY orientation for docking
core::conformation::ResidueOP SnugDockProtocol::place_VRT_at_residue_COM( Pose const & pose, utility::vector1<bool> residues) {
	using namespace core::conformation;
	using namespace core::chemical;
	// Create VRT
	ResidueTypeCOP vrt_res_type = virtual_type_for_pose( pose );
	ResidueOP      vrt_res      = ResidueFactory::create_residue( *vrt_res_type );

	// compute residue nearest to COM that VRT will imitate
	core::Size com_resnum = core::pose::residue_center_of_mass(pose, residues);

	// set VRT coordinates to mirror residue coordinates
	// ORIG = CA
	// X = N
	// Y = C (of res - 1)
	vrt_res->set_xyz( "ORIG", pose.residue(com_resnum).xyz("CA") );
	vrt_res->set_xyz( "X", pose.residue(com_resnum).xyz("N") );
	vrt_res->set_xyz( "Y", pose.residue(com_resnum - 1).xyz("C"));

	return vrt_res;
}

// function returns a virtual residue positioned at COM with correct XY orientation for docking
core::conformation::ResidueOP SnugDockProtocol::place_VRT_at_residue_COM( Pose const & pose, core::Size start, core::Size stop) {
	using namespace core::conformation;
	using namespace core::chemical;
	// Create VRT
	ResidueTypeCOP vrt_res_type = virtual_type_for_pose( pose );
	ResidueOP      vrt_res      = ResidueFactory::create_residue( *vrt_res_type );

	// compute residue nearest to COM that VRT will imitate
	core::Size com_resnum = core::pose::residue_center_of_mass(pose, start, stop);

	// set VRT coordinates to mirror residue coordinates
	// ORIG = CA
	// X = N
	// Y = C (of res - 1)
	vrt_res->set_xyz( "ORIG", pose.residue(com_resnum).xyz("CA") );
	vrt_res->set_xyz( "X", pose.residue(com_resnum).xyz("N") );
	vrt_res->set_xyz( "Y", pose.residue(com_resnum - 1).xyz("C"));

	return vrt_res;
}


void SnugDockProtocol::setup_ab_ag_foldtree( Pose & pose, AntibodyInfoOP antibody_info ) {
	//
	// If original pose is unaltered, make const? **
	//
	// Returns a pose with the some added VRTs and a universal FT
	//
	// General idea is to place VRTs at the COMs for docking
	// and then connect those VRTs to the corresponding n-termini.
	// This should simplify loop modeling on individual chains.
	//
	// For example:
	//     VRT_AB_COM---VRT_VH_COM---VRT_VL_COM
	//      |                |           |
	//      |                |           |
	//      |              VH_Nter     VL_Nter
	//      |
	//     VRT_AG_COM---VRT_AG_A_COM---VRT_AG_B_COM---...
	//                       |           |             |
	//                       |           |             |
	//                     A_Nter       B_Nter        X_Nter

	using namespace core::kinematics;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::pose;

	// start with a vanilla FT and vanilla Pose
	PoseOP pose_vrt( new Pose );
	// initialize PDBInfo
	pose_vrt->pdb_info( utility::pointer::make_shared< PDBInfo >( *pose_vrt, true ) );

	// counter for number of VRTs
	core::Size n_vrts = 0;


	// extract VRTs and prepend to pose, depending on number of chains
	ResidueTypeCOP vrt_res_type = virtual_type_for_pose( pose );
	ResidueOP      vrt_res      = ResidueFactory::create_residue( *vrt_res_type );

	// set up pose with VRTs, we'll need 4+ for Ab-Ag docking: Ab COM, VH COM, (VL COM), Ag COM, A COM, (B COM), ...
	pose_vrt->append_residue_by_jump(*vrt_res, 1); // append one VRT for Ab COM (always res1)
	core::Size ab_com_resnum = 1;
	n_vrts++;

	pose_vrt->append_residue_by_jump(*vrt_res, 1); // append one VRT for Ag COM (always res2)
	core::Size ag_com_resnum = 2;
	n_vrts++;

	// create a map to track COM virtual residues
	// DOES NOT STORE antibody or antigen COM, only for chains
	std::map < std::string, core::Size > com_resnum_from_chain;

	// check for number of Ab chains, append one COM VRT for each chain
	// it will be important to restore these chains to the final pose

	// actually this function returns ["L", "H"] or ["H"], so we're better off with our own
	utility::vector1<std::string> antibody_chain_chars = antibody_info->get_antibody_chains();

	// add H chain (always present) and maybe L chain
	pose_vrt->append_residue_by_jump(*vrt_res, ab_com_resnum + n_vrts - 2);
	n_vrts++;
	com_resnum_from_chain["H"] = n_vrts;
	if ( antibody_chain_chars.index("L") != 0 ) {
		ab_has_light_chain_ = true;
		pose_vrt->append_residue_by_jump(*vrt_res, ab_com_resnum + n_vrts - 2);
		n_vrts++;
		com_resnum_from_chain["L"] = n_vrts;
	}

	// do the same for the Ag
	utility::vector1<std::string> antigen_chain_chars = antibody_info->get_antigen_chains();

	for ( auto chain : antigen_chain_chars ) {
		// append to ag_vrt + n_vrts - 3/4 (one for each COM VRT, depending on # Ab chains)
		pose_vrt->append_residue_by_jump(*vrt_res, ag_com_resnum + n_vrts - antibody_chain_chars.size() - 2);
		n_vrts++;
		com_resnum_from_chain[chain] = n_vrts;
	}

	// fold tree construction each VRT goes to another VRT + N-terminus
	FoldTreeOP ft( new FoldTree() );
	core::Size current_jump_num = 1;

	// connect VRTs
	// 1->2 is always VRT_AB_COM -- VRT_AG_COM
	TR.Debug << "Adding edge from Ab_COM (" << ab_com_resnum << ") to Ag_COM (" << ag_com_resnum << "), #" << current_jump_num << std::endl;
	ft->add_edge(ab_com_resnum, ag_com_resnum, current_jump_num);
	++current_jump_num; // increment jump number (jump 1, initially)
	// 2->3 will be Ab-VH
	TR.Debug << "Adding edge from Ab_COM (" << ab_com_resnum << ") to H (" << com_resnum_from_chain["H"] << "), #" << current_jump_num << std::endl;
	ft->add_edge(ab_com_resnum, com_resnum_from_chain["H"], current_jump_num);
	++current_jump_num;
	// Vh-Vl is 3->4, if present
	if ( ab_has_light_chain_ ) {
		TR.Debug << "Adding edge from H(" << com_resnum_from_chain["H"] << ") to L(" << com_resnum_from_chain["L"] << "), #" << current_jump_num << std::endl;
		ft->add_edge(com_resnum_from_chain["H"], com_resnum_from_chain["L"], current_jump_num);
		vh_vl_jump_ = current_jump_num;
		++current_jump_num;
	}

	// add edge for first antigen chain (ALWAYS present, special case jumps from COM)
	TR.Debug << "Adding edge from Ag_COM (" << ag_com_resnum << ") to Ag_X (" << com_resnum_from_chain[antigen_chain_chars.front()] << "), #" << current_jump_num << std::endl;
	ft->add_edge(ag_com_resnum, com_resnum_from_chain[antigen_chain_chars.front()], current_jump_num);
	++current_jump_num;
	// add edges for all other antigen chains, jump from each other
	// hence int iteration rather than auto
	if ( antigen_chain_chars.size() > 1 ) {
		for ( core::Size i = 2; i <= antigen_chain_chars.size(); ++i ) {
			TR.Debug << "Adding edge from Ag_X (" << com_resnum_from_chain[antigen_chain_chars[i-1]] << ") to Ag_Y (" << com_resnum_from_chain[antigen_chain_chars[i]] << "), #" << current_jump_num << std::endl;
			ft->add_edge(com_resnum_from_chain[antigen_chain_chars[i-1]], com_resnum_from_chain[antigen_chain_chars[i]], current_jump_num);
			++current_jump_num;
		}
	}

	// we should now have connected COM vrts, but where are the actual proteins?!
	// we'll add these below jumps should be appended automatically (?)

	// split the original pose by chain and append each chain to the correct VRT
	// also, while splitting compute individual chain residues closest to COM
	std::map < std::string, ResidueOP > com_VRT_from_chain;

	// heavy chain should be append first to match CDR loop placement later on
	// split particular chain
	PoseOP single_chain = pose.split_by_chain( get_chain_id_from_chain( "H", pose ) );
	// compute residue nearest to COM that VRT will imitate
	com_VRT_from_chain["H"] = place_VRT_at_residue_COM( *single_chain, 1, single_chain->size() );

	// append chain to new pose
	pose_vrt->append_pose_by_jump(*single_chain, com_resnum_from_chain["H"]);

	// light chain is next, if present
	if ( ab_has_light_chain_ ) {
		// split particular chain
		PoseOP single_chain = pose.split_by_chain( get_chain_id_from_chain( "L", pose ) );
		// store in order ORIG coords (equivalent to CA), X = N, Y = C (of -1 res)
		com_VRT_from_chain["L"] = place_VRT_at_residue_COM( *single_chain, 1, single_chain->size() );

		// append chain to new pose
		pose_vrt->append_pose_by_jump(*single_chain, com_resnum_from_chain["L"]);
	}

	// since new pose only contains antibody residues, we can get antibody COM here!
	TR.Debug << "Calculating Ab COM from residue " << n_vrts + 1 << " to " << pose_vrt->size() << std::endl;

	ResidueOP ab_com_VRT = place_VRT_at_residue_COM( *pose_vrt, n_vrts + 1, pose_vrt->size() );

	for ( auto chain : antigen_chain_chars ) {
		// split particular chain
		PoseOP single_chain = pose.split_by_chain( get_chain_id_from_chain( chain, pose ) );
		// compute chain COM
		com_VRT_from_chain[chain] = place_VRT_at_residue_COM(*single_chain, 1, single_chain->size());
		// append chain to new pose
		pose_vrt->append_pose_by_jump(*single_chain, com_resnum_from_chain[chain]);
	}

	// now we have to do a little work to get the antigen center of mass
	// we'll use the original pose for ease (no VRTs)
	// bool vectors to assign region
	utility::vector1<bool> res_in_antigen(pose.size(), false);

	// chain ids
	utility::vector1<core::Size> antigen_chains = antibody_info->get_antigen_chain_ids(pose);

	// iterate over all residues in input pose and assign region
	for ( core::Size i=1; i <= pose.size(); ++i ) {
		// check if current residue is in antigen chain
		if ( antigen_chains.has_value(pose.chain(i)) ) {
			res_in_antigen[i] = true;
		}
	}

	// get residue nearest to center of mass for VRT placement
	ResidueOP antigen_com_VRT = place_VRT_at_residue_COM( pose, res_in_antigen );

	// print and set COMs
	TR.Debug << "VRT residue nearest to antibody COM is " << ab_com_VRT << std::endl;

	// set ORIG to CA COM
	pose_vrt->set_xyz( core::id::AtomID( 1, ab_com_resnum ), ab_com_VRT->xyz("ORIG") );
	// set X to N COM (CA parent)
	pose_vrt->set_xyz( core::id::AtomID( 2, ab_com_resnum ), ab_com_VRT->xyz("X"));
	// set Y to C (of res - 1) COM (N parent)
	pose_vrt->set_xyz( core::id::AtomID( 3, ab_com_resnum ), ab_com_VRT->xyz("Y"));

	TR.Debug << "VRT residue nearest to antigen COM is " << antigen_com_VRT << std::endl;
	// set ORIG to CA COM
	pose_vrt->set_xyz( core::id::AtomID( 1, ag_com_resnum ), antigen_com_VRT->xyz("ORIG"));
	// set X to N COM
	pose_vrt->set_xyz( core::id::AtomID( 2, ag_com_resnum ), antigen_com_VRT->xyz("X"));
	// set Y to C COM
	pose_vrt->set_xyz( core::id::AtomID( 3, ag_com_resnum ), antigen_com_VRT->xyz("Y"));

	// COMs for individual chains
	for ( auto chain : antibody_chain_chars ) {
		TR.Debug << "VRT residue nearest to chain " << chain << " is " << com_VRT_from_chain[chain] << std::endl;
		// set ORIG to CA COM
		pose_vrt->set_xyz( core::id::AtomID( 1, com_resnum_from_chain[chain] ), com_VRT_from_chain[chain]->xyz("ORIG"));
		// set X to N COM
		pose_vrt->set_xyz( core::id::AtomID( 2, com_resnum_from_chain[chain] ), com_VRT_from_chain[chain]->xyz("X"));
		// set Y to C COM
		pose_vrt->set_xyz( core::id::AtomID( 3, com_resnum_from_chain[chain] ), com_VRT_from_chain[chain]->xyz("Y"));
	}
	for ( auto chain : antigen_chain_chars ) {
		TR.Debug << "VRT residue nearest to chain " << chain << " is " << com_VRT_from_chain[chain] << std::endl;
		// set ORIG to CA COM
		pose_vrt->set_xyz( core::id::AtomID( 1, com_resnum_from_chain[chain] ), com_VRT_from_chain[chain]->xyz("ORIG"));
		// set X to N COM
		pose_vrt->set_xyz( core::id::AtomID( 2, com_resnum_from_chain[chain] ), com_VRT_from_chain[chain]->xyz("X"));
		// set Y to C COM
		pose_vrt->set_xyz( core::id::AtomID( 3, com_resnum_from_chain[chain] ), com_VRT_from_chain[chain]->xyz("Y"));
	}

	// OK back to the FoldTree
	// make the assumption that chains are linearly continunous and insert peptide edges

	utility::vector1< core::Size > all_chain_starts; // n-termini for jumps
	utility::vector1< core::Size > all_chain_ends; // c-termini for peptide edges

	utility::vector1< std::pair< core::Size,core::Size > > continuous_chain_segments;

	core::Size i = n_vrts+1; // start at first real residue
	core::Size current_chain = pose_vrt->chain(i); // chains are numeric, the first is 1?
	core::Size seg_start = i; // initialize segment start
	all_chain_starts.push_back(i); // first chain start

	while ( i <= pose_vrt->size() ) {
		// end of continuous segment condition
		if ( current_chain != pose_vrt->chain(i) ) {
			continuous_chain_segments.push_back( {seg_start,i-1} );
			current_chain = pose_vrt->chain(i);
			seg_start = i;
			all_chain_starts.push_back(i);
			all_chain_ends.push_back(i-1);
		}
		++i;
	}

	all_chain_ends.push_back(pose_vrt->size()); // last chain end

	TR.Debug << "Chain starts: " << all_chain_starts << "." << std::endl;
	TR.Debug << "Chain ends: " << all_chain_ends << "." << std::endl;

	// sanity check
	if ( all_chain_ends.size() != all_chain_starts.size() ) {
		utility_exit_with_message("In SnugDock Protocol: Number of chain starts != number of chain ends!");
	}

	// update ft -- needs jumps from COM VRTs to N-termini
	// needs peptide edges from N-termini to C-termini
	for ( core::Size i = 1; i <= all_chain_starts.size(); ++i ) {

		core::Size nter_resn = all_chain_starts[i];
		core::Size cter_resn = all_chain_ends[i];

		// get com_vert resnum by matching chains, assumes correctly assigned chains
		TR.Debug << "Adding jump from COM VRT to N-term: " << com_resnum_from_chain[pose_vrt->pdb_info()->chain(nter_resn)] << " " << nter_resn << "." << std::endl;
		ft->add_edge(com_resnum_from_chain[pose_vrt->pdb_info()->chain(nter_resn)], nter_resn, current_jump_num);
		++current_jump_num;

		// peptide edge!
		ft->add_edge(nter_resn, cter_resn, -1);

		// store residue information here in sets
		// antibody chains come first (2 at most)
		for ( core::Size j = nter_resn; j <= cter_resn; ++j ) {
			if ( ab_has_light_chain_ ) {
				// H and L come first and have passed, so store antigen residues
				if ( i > 2 ) { ag_residues_.insert(j); }
				// otherwise, store antibody residues
				else { ab_residues_.insert(j); }
			} else {
				// H comes first and has passed, so store antigen residues
				if ( i > 1 ) { ag_residues_.insert(j); }
				// otherwise, store antibody residues
				else { ab_residues_.insert(j); }
			}
		}
	}

	// next, we need CDRs in our FT!

	// careful if flipping chain positions!
	// extract loops from individual chains (can we feed Ab info just H/L???) and place correctly!
	PoseOP ab_only( new Pose );
	if ( ab_has_light_chain_ ) {
		ab_only = pose.split_by_chain().at(core::pose::get_chain_id_from_chain("H", pose));
		ab_only->append_pose_by_jump( *pose.split_by_chain().at(core::pose::get_chain_id_from_chain("L", pose)), 1 );
	} else {
		ab_only = pose.split_by_chain().at(core::pose::get_chain_id_from_chain("H", pose));
	}

	// first, get old CDR loops
	// careful, this returns the default definition of the loops
	//loops::LoopsOP old_cdr_loops = antibody_info->get_CDR_loops( pose );
	AntibodyInfoOP HL_info( new AntibodyInfo(*ab_only, North) );
	loops::LoopsOP old_cdr_loops = HL_info->get_CDR_loops( *ab_only );
	// try me if above fails
	//loops::LoopsOp old_cdr_loops = loops::LoopsOp( new loops::Loops() );
	//antibody_info->get_AllCDRs_in_loopsop();


	// then, update for VRT offset
	// assuming only VRT resiudes precede antibody residues in new FT
	// pair loops with region so as to not lose track
	core::Size loop_enum_counter = 1;

	for ( auto loops_it = old_cdr_loops->begin(); loops_it != old_cdr_loops->end(); ++loops_it ) {

		loops::Loop loop = *loops_it;
		TR.Debug << "Adding loop: " << loop_enum_counter << std::endl;

		TR.Debug << "Old loop position: " << loop.start() << " " << loop.stop() << " " << loop.cut() << std::endl;
		TR.Debug << "Old loop residues: " << ab_only->residue(loop.start()).name() << " " << ab_only->residue(loop.stop()).name() << " " << ab_only->residue(loop.cut()).name() << std::endl;
		// create new loop
		// loops are ordered by enum h1 (1), h2 (2), h3 (3), l1 (4), l2 (5), l3 (6) ...
		// hence the counter
		CDR_loops_[static_cast<CDRNameEnum>(loop_enum_counter)] = loops::Loop(loop.start() + n_vrts, loop.stop() + n_vrts, loop.cut() + n_vrts);
		++loop_enum_counter;

		// update fold tree
		// do not add flanking residues, rely on North (structural) CDR definition
		TR.Debug << "New loop position: " << loop.start() + n_vrts << " " << loop.stop() + n_vrts << " " << loop.cut() + n_vrts << std::endl;
		TR.Debug << "New loop residues: " << pose_vrt->residue(loop.start() + n_vrts).name() << " " << pose_vrt->residue(loop.stop() + n_vrts).name() << " " << pose_vrt->residue(loop.cut() + n_vrts).name() << std::endl;

		ft->new_jump(loop.start() + n_vrts, loop.stop() + n_vrts, loop.cut() + n_vrts);
	}
	// easy!

	// check mah fold tree
	TR << "Is my new SnugDock fold tree ok? " << ft->check_fold_tree() << std::endl;

	// apply new ft to new pose
	pose_vrt->fold_tree(*ft);

	// update original pose
	pose = *pose_vrt;

	// for debugging purposes, convert VRTs to reals and dump
	//for ( i=1; i<=n_vrts; ++i ) {
	// pose.virtual_to_real(i);
	//}
	pose.pdb_info()->obsolete(false);
	//pose.dump_pdb("universal_foldtree_pose.pdb");
	// TODO: write function to remove VRTs and restore simple FT so -native flag doesn't break
}

protocols::loop_modeler::LoopModelerOP SnugDockProtocol::refine_loop() const {
	protocols::loop_modeler::LoopModelerOP refine_kic ( utility::pointer::make_shared<  protocols::loop_modeler::LoopModeler >() );

	// if fragments are given on the command-line, run fKIC, else just NGK (default)
	if ( basic::options::option[ basic::options::OptionKeys::loops::frag_files ].user() ) {
		refine_kic->setup_kic_with_fragments_config();
	}

	// do not build
	refine_kic->disable_build_stage();

	// pose is centroid, disable fullatom
	refine_kic->disable_fullatom_stage();

	// do not change FT during loop modeling
	refine_kic->trust_fold_tree();

	// If testing, run reduced cycles
	if ( basic::options::option[ basic::options::OptionKeys::run::test_cycles ] ) {
		refine_kic->centroid_stage()->mark_as_test_run();
	}

	return refine_kic;
}

docking::DockingProtocolOP SnugDockProtocol::docking() const {
	using namespace docking;
	using namespace core::scoring;
	if ( ! docking_ ) {
		/// The full DockingProtocol is used with a custom high resolution phase and post-low-resolution phase
		/// All FoldTrees will be setup through AntibodyInfo so DockingProtocol's autofoldtree setup is disabled.
		docking_ = utility::pointer::make_shared< docking::DockingProtocol >(utility::tools::make_vector1<core::SSize>(1), false, false, false, nullptr, nullptr);
	}
	return docking_;
}

void
SnugDockProtocol::show( std::ostream & out ) const {
	out << *this;
}

std::ostream & operator<<(std::ostream& out, SnugDockProtocol const & snugdockprotocol ) {
	if ( snugdockprotocol.antibody_info_ ) {
		out << snugdockprotocol.get_name() << " has been configured to operate on an Antibody-Antigen complex with the "
			<< "following information:" << std::endl;
		out << * snugdockprotocol.antibody_info_ << std::endl;
	} else {
		out << snugdockprotocol.get_name() << " has not been used yet.  " << snugdockprotocol.get_name()
			<<"'s data initialization occurs the first time its apply method is called." << std::endl;
	}
	return out;
}

} // namespace snugdock
} // namespace antibody
} // namespace protocols
