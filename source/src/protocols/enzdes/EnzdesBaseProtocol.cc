// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/enzdes/EnzdesBaseProtocol.cc
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu


//unit headers
#include <protocols/enzdes/EnzdesBaseProtocol.hh>

//package headers
#include <protocols/enzdes/AddorRemoveCsts.hh>
#include <protocols/enzdes/EnzdesTaskOperations.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesSeqRecoveryCache.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>
#include <protocols/toolbox/task_operations/LimitAromaChi2Operation.hh>
#include <protocols/enzdes/enzdes_util.hh>
#include <protocols/enzdes/ModifyStoredLigandRBConfsMovers.hh>
#include <protocols/ligand_docking/ligand_functions.hh> //for minimizing ligand torsions
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

//project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
// AUTO-REMOVED #include <core/io/pdb/file_data.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/options/option.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
 //needed for adding variant types in cst opt
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <basic/Tracer.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

#include <protocols/simple_moves/MinMover.hh>
// AUTO-REMOVED #include <utility/string_util.hh>

#include <utility/io/izstream.hh>


// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/id/SequenceMapping.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>



namespace protocols{
namespace enzdes{

static basic::Tracer tr("protocols.enzdes.EnzdesBaseProtocol");

EnzdesBaseProtocol::EnzdesBaseProtocol():
	LigandBaseProtocol(),
	bb_min_allowed_dev_(0.5),
  loop_bb_min_allowed_dev_(0.5),
	lig_superposition_file_read_(false),
	rb_min_(true),
  exclude_protein_protein_fa_elec_(false) //sboyken 01/22/14; changed true to false; with Talaris2013 don't want this on by default; only with old_estat
{
		rb_min_jumps_.clear();
		Mover::type( "EnzdesFixBBProtocol" );
		restype_set_ = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) ;
		//cst_io_ = new toolbox::match_enzdes_util::EnzConstraintIO(restype_set_);
		//catalytic_res_.clear();
		atoms_to_superimpose_on_.clear();
		bb_min_allowed_dev_ = basic::options::option[ basic::options::OptionKeys::enzdes::bb_min_allowed_dev];
		loop_bb_min_allowed_dev_ = basic::options::option[ basic::options::OptionKeys::enzdes::loop_bb_min_allowed_dev];
		include_all_design_targets_in_design_interface_ = basic::options::option[ basic::options::OptionKeys::enzdes::include_catres_in_interface_detection];

		//		std::string score_patch = basic::options::option[ basic::options::OptionKeys::score::patch ];

		//		if( score_patch == "" )	reduced_sfxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "enzdes_polyA_min" );
		//		else reduced_sfxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "enzdes_polyA_min", score_patch );
		using namespace basic::options;
		using namespace core::scoring;
		reduced_sfxn_ = ScoreFunctionFactory::create_score_function( "enzdes_polyA_min", option[ OptionKeys::score::patch ]() );

		if( basic::options::option[basic::options::OptionKeys::enzdes::chi_min].user() )chi_min_ = true;
		else chi_min_= false;

		if( basic::options::option[basic::options::OptionKeys::enzdes::bb_min].user() )bb_min_ = true;
		else bb_min_= false;

		bb_backrub_ = false; //not really used here so false, but can be set to true using the set_minimize_options function

		if( basic::options::option[basic::options::OptionKeys::enzdes::min_all_jumps].user() )min_all_jumps_ = true;
		else min_all_jumps_ = false;

		if( basic::options::option[basic::options::OptionKeys::enzdes::minimize_ligand_torsions].user() ){
			minimize_ligand_torsions_=true;
			minimize_all_ligand_torsions_= false;
			lig_min_stddev_= basic::options::option[basic::options::OptionKeys::enzdes::minimize_ligand_torsions].value();
		}
		else {
			if( basic::options::option[basic::options::OptionKeys::enzdes::minimize_all_ligand_torsions].user() ){
				minimize_ligand_torsions_=true;
				minimize_all_ligand_torsions_= true;
				lig_min_stddev_= basic::options::option[basic::options::OptionKeys::enzdes::minimize_ligand_torsions].value();
			}
			else {
				minimize_ligand_torsions_=false;
				minimize_all_ligand_torsions_= false;
				lig_min_stddev_= 0.0;
			}
		}
		if( basic::options::option[basic::options::OptionKeys::enzdes::fix_catalytic_aa].user() ) fix_catalytic_aa_ = true;
		else fix_catalytic_aa_= false;

		//if( basic::options::option[basic::options::OptionKeys::score::weights].user() ){
		//	std::string weights_tag = basic::options::option[basic::options::OptionKeys::score::weights];
		//	scorefxn_->initialize_from_file( basic::database::full_name( "scoring/weights/"+weights_tag+".wts" ) );
		//}
		if ( basic::options::option[ basic::options::OptionKeys::score::weights ].user() ) {
			scorefxn_ = core::scoring::getScoreFunction(); // This call handles the database vs working directory resolution -- DONT SUBVERT OR DUPLICATE IT
		}
		else {
			scorefxn_ = ScoreFunctionFactory::create_score_function( "talaris2013_cst", option[ OptionKeys::score::patch ]() ); //02/25/14 sboyken; changed default to talaris2013_cst
			/*			if( score_patch == "" ) scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function("enzdes");
			else scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function("enzdes", score_patch);
			*/

		}

		if (scorefxn_->has_zero_weight( core::scoring::coordinate_constraint ) ){
			constraint_weights_[core::scoring::coordinate_constraint] = 1.0;
		}
		else{
			constraint_weights_[core::scoring::coordinate_constraint] = scorefxn_->weights()[core::scoring::coordinate_constraint];
		}
		if (scorefxn_->has_zero_weight( core::scoring::atom_pair_constraint ) ){
			constraint_weights_[core::scoring::atom_pair_constraint] = 1.0;
		}
		else{
			constraint_weights_[core::scoring::atom_pair_constraint] = scorefxn_->weights()[core::scoring::atom_pair_constraint];
		}
		if (scorefxn_->has_zero_weight( core::scoring::angle_constraint ) ){
			constraint_weights_[core::scoring::angle_constraint] = 1.0;
		}
		else{
			constraint_weights_[core::scoring::angle_constraint] = scorefxn_->weights()[core::scoring::angle_constraint];
		}
		if (scorefxn_->has_zero_weight( core::scoring::dihedral_constraint ) ){
			constraint_weights_[core::scoring::dihedral_constraint] = 1.0;
		}
		else{
			constraint_weights_[core::scoring::dihedral_constraint] = scorefxn_->weights()[core::scoring::dihedral_constraint];
		}

		if( basic::options::option[basic::options::OptionKeys::enzdes::favor_native_res].user() || basic::options::option[ basic::options::OptionKeys::in::file::pssm ].user() ){
			if (scorefxn_->has_zero_weight( core::scoring::res_type_constraint ) ){
				constraint_weights_[core::scoring::res_type_constraint] = 1.0;
			}
			else{
				constraint_weights_[core::scoring::res_type_constraint] = scorefxn_->weights()[core::scoring::res_type_constraint];
			}
		}
		else constraint_weights_[core::scoring::res_type_constraint] = scorefxn_->weights()[core::scoring::res_type_constraint];

		enable_constraint_scoreterms();

		if( basic::options::option[ basic::options::OptionKeys::docking::ligand::old_estat ].user() ){
			exclude_protein_protein_fa_elec_ = basic::options::option[ basic::options::OptionKeys::docking::ligand::old_estat ];
		}
		if( exclude_protein_protein_fa_elec_ ){
			core::scoring::methods::EnergyMethodOptions options( scorefxn_->energy_method_options() );
			options.exclude_protein_protein_fa_elec( true );
			scorefxn_->set_energy_method_options( options );
		}


		//set the native pose if requested
		if( basic::options::option[basic::options::OptionKeys::in::file::native].user() ){
			core::pose::PoseOP natpose = new core::pose::Pose();
			core::import_pose::pose_from_pdb( *natpose, basic::options::option[basic::options::OptionKeys::in::file::native].value() );
			(*scorefxn_)( *natpose);
			this->set_native_pose( natpose );
		}

		//increase the chainbreak weight. 1.0 is apparently not enough for some constraints
		//scorefxn_->set_weight( core::scoring::chainbreak, 10.0 );

} //EnzdesBaseProtocol constructor


std::string
EnzdesBaseProtocol::get_name() const {
	return "EnzdesBaseProtocol";
}

utility::vector1< core::Size >
EnzdesBaseProtocol::catalytic_res( core::pose::Pose const & pose ) const
{
	using namespace core;
	utility::vector1< Size > to_return;
	protocols::toolbox::match_enzdes_util::EnzdesCacheableObserverCOP enz_obs( toolbox::match_enzdes_util::get_enzdes_observer( pose ) ); // toolbox::match_enzdes_util::get_enzdes_observer() for const pose can return NULL
	if( enz_obs ) {
		toolbox::match_enzdes_util::EnzdesCstCacheCOP cstcache (enz_obs->cst_cache() );
		if( cstcache ) to_return = cstcache->enzcst_io()->ordered_constrained_positions( pose );
	}

	if( to_return.size() == 0 ){
		for(core::Size i = 1, i_end = pose.total_residue(); i <= i_end; ++i){
			if( pose.residue_type( i ).is_ligand() ) to_return.push_back( i );
		}
	}
	return to_return;
}

std::set< core::Size > const &
EnzdesBaseProtocol::design_targets( core::pose::Pose const & pose ) const
{
	design_targets_.clear();

	for(core::Size i = 1; i <=  pose.total_residue(); ++i){
		if( enzutil::is_catalytic_seqpos( pose, i ) ) design_targets_.insert( i ); //Changed to is_catalytic_seqpos to avoid problems with running parser with no constraints -mdsmith
	}
		//if no positions are constrained, we'll put the ligands into the array
	if( design_targets_.size() == 0 ){
		for(core::Size i = 1; i <=  pose.total_residue(); ++i){
			if( pose.residue_type(i).is_ligand() ) design_targets_.insert( i );
		}
	}
	return design_targets_;
}


void
EnzdesBaseProtocol::register_options()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( in::file::native );
	option.add_relevant( in::file::s );
	option.add_relevant( in::file::l );
	option.add_relevant( OptionKeys::score::weights );
	option.add_relevant( OptionKeys::score::patch );
	option.add_relevant( OptionKeys::packing::soft_rep_design );

	protocols::enzdes::DetectProteinLigandInterface::register_options();
	protocols::enzdes::ProteinLigandInterfaceUpweighter::register_options();
	protocols::simple_moves::MinMover::register_options();
	protocols::simple_moves::PackRotamersMover::register_options();

	option.add_relevant( OptionKeys::enzdes::detect_design_interface );
	option.add_relevant( OptionKeys::enzdes::include_catres_in_interface_detection );
	option.add_relevant( OptionKeys::enzdes::fix_catalytic_aa );
	option.add_relevant( OptionKeys::enzdes::ex_catalytic_rot );
	option.add_relevant( OptionKeys::enzdes::cst_min );
	option.add_relevant( OptionKeys::enzdes::chi_min );
	option.add_relevant( OptionKeys::enzdes::bb_min );
	option.add_relevant( OptionKeys::enzdes::bb_min_allowed_dev );
	option.add_relevant( OptionKeys::enzdes::loop_bb_min_allowed_dev );
	option.add_relevant( OptionKeys::enzdes::enz_debug );
	option.add_relevant( OptionKeys::enzdes::no_packstat_calculation );
	option.add_relevant( OptionKeys::enzdes::compare_native );
	option.add_relevant( OptionKeys::enzdes::favor_native_res );
	option.add_relevant( OptionKeys::docking::ligand::old_estat );

}


core::pack::task::PackerTaskOP
EnzdesBaseProtocol::create_enzdes_pack_task(
	core::pose::Pose & pose,
	bool design
){

	using namespace core::pack::task;
	using namespace basic::options;

	//make sure the design targets are up to date
	design_targets( pose );

  DetectProteinLigandInterfaceOP detect_enzdes_interface = new DetectProteinLigandInterface();
  detect_enzdes_interface->set_design(design);
	if( include_all_design_targets_in_design_interface_ ){
		detect_enzdes_interface->set_design_target_res( design_targets_ );
	}
  TaskFactory taskfactory;
  taskfactory.push_back( new operation::InitializeFromCommandline );
  taskfactory.push_back( detect_enzdes_interface);
	if( design ) { // upweight ligand interactions and remove bad aromatic rotamers during design only
		taskfactory.push_back( new ProteinLigandInterfaceUpweighter() );
		taskfactory.push_back( new toolbox::task_operations::LimitAromaChi2Operation() );
	}
	if( toolbox::match_enzdes_util::get_enzdes_observer( pose ) ){
		taskfactory.push_back( new AddRigidBodyLigandConfs() );
	}
	if( basic::options::option[basic::options::OptionKeys::enzdes::detect_design_interface].user() ){
		SetCatalyticResPackBehaviorOP catpack = new SetCatalyticResPackBehavior();
		catpack->set_fix_catalytic_aa( this->fix_catalytic_aa_ );
		taskfactory.push_back( catpack );
	}
	if( basic::options::option[basic::options::OptionKeys::enzdes::run_ligand_motifs].user() ){
		taskfactory.push_back( new AddLigandMotifRotamers() );
	}

	PackerTaskOP task = taskfactory.create_task_and_apply_taskoperations( pose );
	task->append_rotamerset_operation( unboundrot_ );

	setup_sequence_recovery_cache( pose, *task );
	return task;

} //create_enzdes_pack_task

void
EnzdesBaseProtocol::setup_sequence_recovery_cache(
	core::pose::Pose & pose,
	core::pack::task::PackerTask const & task
) const
{
	//Initiate sequence recovery cache in the enzdes observer
	//Set wt sequence at the same time. That is the first and only time
	//that the wt sequence gets initiated
  if ( ! toolbox::match_enzdes_util::get_enzdes_observer( pose ) -> get_seq_recovery_cache() ){
    toolbox::match_enzdes_util::get_enzdes_observer( pose ) -> set_seq_recovery_cache( new toolbox::match_enzdes_util:: EnzdesSeqRecoveryCache );
  	toolbox::match_enzdes_util::get_enzdes_observer( pose ) -> get_seq_recovery_cache() -> set_sequence( pose );
  }

  //keep track of what residues we are designing
 	std::set < core::Size > designing_residues;
	for(core::Size jj=1; jj<=pose.total_residue(); ++jj) {
		if (pose.residue(jj).is_protein() && task.being_designed(jj) ) {
			designing_residues.insert( jj );
		}
	}
	//update or initiate what residues are designed in the EnzdesSeqRecoveryCache
	toolbox::match_enzdes_util::get_enzdes_observer( pose ) -> get_seq_recovery_cache()
		-> set_designable_residues( designing_residues );
}

///@details This function will modify the fold tree and add constraints to the pose if used with the bb_min and minimize_ligand_torsions options. Be WARNED!
core::kinematics::MoveMapOP
EnzdesBaseProtocol::create_enzdes_movemap(
  core::pose::Pose & pose,
  core::pack::task::PackerTaskCOP task,
  bool min_all_jumps
) const
{
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap();
	movemap->set_jump( false );
	movemap->set_chi( false );
	movemap->set_bb( false );
	core::Size jump_id (pose.num_jump()); //default jump-id last
	utility::vector1< bool > allow_move_bb(pose.total_residue(), false );

	for(core::Size i = 1, i_end = pose.total_residue(); i <= i_end; ++i) {
		if ( task->pack_residue(i) && pose.residue(i).is_polymer() ){
			if(chi_min_ ) movemap->set_chi(i, true);
			if(bb_min_ ) {
				allow_move_bb[i] = true;
			}
		}
	}

 if( rb_min_jumps().size() > 0 ){
	 tr<<"rb_min_jumps was set. setting the following rb dofs to true: ";
	 foreach( core::Size const rb, rb_min_jumps() ){
		 tr<<rb<<',';
		 movemap->set_jump( rb, true );
	 }
 }
 else if (rb_min_) {
    for (core::Size i = 1;i<=pose.num_jump();++i){
			if (min_all_jumps) movemap->set_jump( i, true );
			else {
      	core::Size upstream_jump_res, downstream_jump_res;
      	upstream_jump_res = pose.fold_tree().upstream_jump_residue( i );
      	downstream_jump_res = pose.fold_tree().downstream_jump_residue( i );
      	if ( pose.residue( upstream_jump_res ).is_ligand()  || pose.residue( downstream_jump_res ).is_ligand()  ) {
					movemap->set_jump( i, true );
					jump_id= i;
				}
			}
    }
  }

	//for the minimization to work properly, each residue that is minimized has to be part of a stretch of at least 4 other residues
	//so we need to allow a couple of more residues to move
	//core::Real initial_cbreak;
	if( bb_min_ || bb_backrub_ ){
		core::Size window = 4;

		enzutil::make_continuous_true_regions_in_bool_vector( allow_move_bb, window );
		tr.Info << "Doing a pose minimization... the backbone is allowed to move at positions: ";
		for(core::Size i = 1; i <= pose.total_residue(); ++i){
			if( allow_move_bb[i] && pose.residue(i).is_polymer() ){
				movemap->set_bb(i, true);
				tr.Info << i <<", ";
			}
		}
		tr.Info << std::endl;
		if (bb_min_) setup_bbmin_ft_and_csts( pose, allow_move_bb, jump_id ); //NOTE: this will modify fold tree and constraints set of the pose
	}

  //allow ligand minimization if requested
	if (minimize_ligand_torsions_) {
		core::Size const lig_id = get_ligand_id(pose, jump_id);
		using namespace protocols::ligand_docking;
		core::Real lig_min_stddev = lig_min_stddev_;
		bool minimize_all_ligand_torsions = minimize_all_ligand_torsions_;
		if (lig_min_stddev==0.0){
			tr.Info<< "**WARNING!! ligand minimization requested without specifying allowed deviation value!! Setting it to 10.0 degrees and allowing all torsions to minimize"<<std::endl;
			lig_min_stddev= 10.0;
			minimize_all_ligand_torsions = true;
		}
		constrain_ligand_torsions(pose,lig_min_stddev,minimize_all_ligand_torsions); //NOTE: this will modify the constraints set of the pose
		movemap->set_chi(lig_id, true);
	}
return movemap;
}

void
EnzdesBaseProtocol::setup_bbmin_ft_and_csts(
core::pose::Pose & pose,
utility::vector1< bool > allow_move_bb,
core::Size jump_id ) const
{
	core::Size const lig_id = jump_id !=0 ? get_ligand_id(pose, jump_id): 0;
  //restraining function for Calphas. should allow fairly liberal movement ~0.1A from the original position,
  //but severly limits movement beyond this
  core::scoring::func::FuncOP ss_ca_restr_func = new core::scoring::constraints::BoundFunc( 0, bb_min_allowed_dev_, 0.1, "CAdis");
  core::scoring::func::FuncOP loop_ca_restr_func = new core::scoring::constraints::BoundFunc( 0, loop_bb_min_allowed_dev_, 0.1, "CAdis");

	//note flo feb '11 reordering the foldtree will change residue types and can thus have an
	//effect on constraints, which is a nasty bug bc different constraints will be enforced.
	//to prevent this, we clone the constraint set from the old pose to the new.
 	//prolly not the cheapest way to do stuff, but wutevs..
	if( pose.constraint_set()->has_constraints() ){
		core::pose::Pose init_pose = pose;
  	reorder_foldtree_around_mobile_regions( pose, jump_id, allow_move_bb, lig_id );
 		pose.constraint_set( init_pose.constraint_set()->remapped_clone( init_pose, pose ) );
	}
	else reorder_foldtree_around_mobile_regions( pose, jump_id, allow_move_bb, lig_id );

  core::scoring::dssp::Dssp ss_pose(pose);
  utility::vector1< bool > allow_move_bb_loop(pose.total_residue(), false );
  utility::vector1< bool > allow_move_bb_ss(pose.total_residue(), false );
  for (core::Size i=1; i<=pose.total_residue(); ++i){
		if( !pose.residue(i).is_protein() ) continue;
    if(ss_pose.get_dssp_secstruct(i) == ' ' && allow_move_bb[i]){
      allow_move_bb_loop[i]=true;
      } else if(allow_move_bb[i]){
    	  allow_move_bb_ss[i]=true;
      }
    }
  restrain_protein_Calphas(pose, allow_move_bb_ss, ss_ca_restr_func );
  restrain_protein_Calphas(pose, allow_move_bb_loop, loop_ca_restr_func );
}

void
EnzdesBaseProtocol::enzdes_pack(
	core::pose::Pose & pose,
	core::pack::task::PackerTaskCOP task,
	core::scoring::ScoreFunctionCOP scorefxn,
	core::Size cycles,
	bool minimize_after_packing,
	bool pack_unconstrained,
	bool favor_native
) const
{

	if( pack_unconstrained ) remove_enzdes_constraints( pose, true );

	if( favor_native ){
		toolbox::match_enzdes_util::get_enzdes_observer( pose )->setup_favor_native_constraints( pose, task, *(this->get_native_pose()) );
	}

	core::pack::task::PackerTaskCOP usetask = task;
	core::scoring::ScoreFunctionCOP packsfxn;

	for( core::Size cycle = 1; cycle <= cycles; ++cycle){

		//soft_rep we'll only do this if the task is to be designed and not the last cycle,
		//because we really don't want any clashes
		bool soft_rep ( basic::options::option[basic::options::OptionKeys::packing::soft_rep_design] && usetask->design_any() && (cycle < cycles ) );

		if( soft_rep ) packsfxn = soft_scorefxn_;
		else packsfxn = scorefxn;
		protocols::simple_moves::PackRotamersMoverOP enzdes_pack = new protocols::simple_moves::PackRotamersMover(packsfxn, usetask);

		enzdes_pack->apply(pose);

		if(minimize_after_packing) cst_minimize(pose, task);

		usetask = enzutil::recreate_task( pose, *task );

	} //cycle loop


	if( basic::options::option[basic::options::OptionKeys::enzdes::favor_native_res].user() && favor_native){
		toolbox::match_enzdes_util::get_enzdes_observer( pose )->remove_favor_native_constraints( pose );
	}

	if( pack_unconstrained) add_pregenerated_enzdes_constraints( pose );

	(*scorefxn)( pose );

	if( task->design_any() ){
		toolbox::match_enzdes_util::EnzConstraintIOCOP cstio( enzutil::get_enzcst_io( pose ) );
		if( cstio )	cstio->update_pdb_remarks_for_backbone_params( pose );
	}

} //design function


void
EnzdesBaseProtocol::setup_enzdes_constraints(
	core::pose::Pose & pose,
	bool allow_missing_remark_blocks
) const {

	AddOrRemoveMatchCsts cstmover;
	cstmover.set_cst_action( ADD_NEW );
	cstmover.set_accept_blocks_missing_header( allow_missing_remark_blocks );
	cstmover.apply( pose );

	tr.Info << "Catalytic residues (pose numbering) are: ";
	for( core::Size i = 1; i <= pose.total_residue(); ++i){
		if( is_catalytic_position( pose, i ) ) tr.Info << i << " ";
	}
	tr.Info << std::endl;

} //setup enzdes constraints function


void
EnzdesBaseProtocol::remove_enzdes_constraints(
	core::pose::Pose & pose,
	bool keep_covalent
) const
{
	AddOrRemoveMatchCsts cstmover;
	cstmover.set_cst_action( REMOVE );
	cstmover.set_keep_covalent( keep_covalent );
	cstmover.apply( pose );
}

void
EnzdesBaseProtocol::add_pregenerated_enzdes_constraints(
	core::pose::Pose & pose
) const
{
	AddOrRemoveMatchCsts cstmover;
	cstmover.set_cst_action( ADD_PREGENERATED );
	cstmover.apply( pose );
}


void
EnzdesBaseProtocol::cst_minimize(
	core::pose::Pose & pose,
	core::pack::task::PackerTaskCOP task,
	bool cst_opt
) const
{


	core::pose::Pose old_Pose = pose; //copy old pose
	core::scoring::ScoreFunctionOP min_scorefxn;
	//core::scoring::constraints::ConstraintSetOP saved_constraints;

	if(cst_opt){

		min_scorefxn = reduced_sfxn_;


		utility::vector1< core::Size >positions_to_replace;

		for(core::Size i = 1, i_end = pose.total_residue(); i <= i_end; ++i) {

			if( task->pack_residue(i) && ( ! is_catalytic_position( pose, i ) ) ) positions_to_replace.push_back( i );
		}

		protocols::toolbox::pose_manipulation::construct_poly_ala_pose( pose, positions_to_replace, true, true, true );

	}

	else { min_scorefxn = scorefxn_; }

	if(basic::options::option[basic::options::OptionKeys::enzdes::enz_debug] ){
		//debug stage: only interested in constraints minimization for now
		//restrict move map to only move ligand along jump and the three
		//catalytic res
		min_scorefxn = reduced_sfxn_->clone();
		min_scorefxn->reset();
		min_scorefxn->set_weight(core::scoring::chainbreak, 1.0);
		min_scorefxn->set_weight(core::scoring::omega, 0.0 );
		min_scorefxn->set_weight(core::scoring::coordinate_constraint, 1.0);
		min_scorefxn->set_weight(core::scoring::atom_pair_constraint, 1.0);
		min_scorefxn->set_weight(core::scoring::angle_constraint, 1.0);
		min_scorefxn->set_weight(core::scoring::dihedral_constraint, 1.0);

	}

	//temporarily set chainbreak score very high, to prevent chainbreaks from opening
	//this has sometimes been observed in case of highly constrained systems
	core::Real orig_cbreak_weight = min_scorefxn->weights()[ core::scoring::chainbreak ];
	if( orig_cbreak_weight < 10.0 ) min_scorefxn->set_weight( core::scoring::chainbreak, 10.0);

	//create movemap
	 core::kinematics::MoveMapOP movemap = create_enzdes_movemap( pose, task, min_all_jumps_);
	//setting up move map done, now do minimization
	protocols::simple_moves::MinMoverOP dfpMinTightTol = new protocols::simple_moves::MinMover( movemap, min_scorefxn, "dfpmin_armijo_nonmonotone_atol", 0.02, true /*use_nblist*/ );
	dfpMinTightTol->apply(pose);

	min_scorefxn->set_weight( core::scoring::chainbreak, orig_cbreak_weight);
	(*min_scorefxn)( pose );

	MinimizeStoredRBConfs stored_rb_min( min_scorefxn );
	stored_rb_min.apply( pose );
	//core::Real post_cbreak = pose.energies().total_energies()[ core::scoring::chainbreak ];
	//std::cerr << "chainbreak score postmin " << post_cbreak << std::endl;


	//now we have to reinstate the orginial pose, so the task is still valid
	if(cst_opt){

		for(core::Size i = 1, i_end = pose.total_residue(); i <= i_end; ++i) {

			if ( task->pack_residue(i) && ( !is_catalytic_position( pose, i ) )  ){
				pose.replace_residue( i, old_Pose.residue(i), true );
			}
		}
	}

	if(bb_min_ || minimize_ligand_torsions_){

		pose.constraint_set( old_Pose.constraint_set()->clone() );
		pose.fold_tree( old_Pose.fold_tree() );
		if (bb_min_){
		//put back the right variants
			for(core::Size i = 1, i_end = pose.total_residue(); i <= i_end; ++i) {

				if( !pose.residue_type( i ).variants_match( old_Pose.residue_type( i ) ) ){

					utility::vector1< std::string > const new_var_types( pose.residue_type( i ).variant_types() );
					utility::vector1< std::string > const old_var_types( old_Pose.residue_type( i ).variant_types() );
					for( utility::vector1< std::string >::const_iterator newvars = new_var_types.begin(); newvars  != new_var_types.end(); ++newvars ){
						if( !old_Pose.residue_type( i ).has_variant_type( *newvars ) ) core::pose::remove_variant_type_from_pose_residue( pose, *newvars, i );
					}

					for( utility::vector1< std::string >::const_iterator oldvars = old_var_types.begin(); oldvars  != old_var_types.end(); ++oldvars ){
						if( !pose.residue_type( i ).has_variant_type( *oldvars ) ) core::pose::add_variant_type_to_pose_residue( pose, *oldvars, i );
					}


				} //if variants don't match
			}
		} //if bb_min
		(*min_scorefxn)( pose ); //just to be safe
	} //if bb_min || minimize_ligand
	if(basic::options::option[basic::options::OptionKeys::enzdes::enz_debug] ){
		pose.dump_scored_pdb("aftermin_pose.pdb", *min_scorefxn);
	}

}

bool
EnzdesBaseProtocol::is_catalytic_position( core::pose::Pose const & pose, core::Size const seqpos ) const
{
	return enzutil::is_catalytic_seqpos( pose, seqpos );
}

void
EnzdesBaseProtocol::enable_constraint_scoreterms(){

	scorefxn_->set_weight(core::scoring::coordinate_constraint, constraint_weights_[core::scoring::coordinate_constraint] );
	scorefxn_->set_weight(core::scoring::atom_pair_constraint, constraint_weights_[core::scoring::atom_pair_constraint] );
	scorefxn_->set_weight(core::scoring::angle_constraint, constraint_weights_[core::scoring::angle_constraint] );
	scorefxn_->set_weight(core::scoring::dihedral_constraint, constraint_weights_[core::scoring::dihedral_constraint] );
	scorefxn_->set_weight(core::scoring::res_type_constraint, constraint_weights_[core::scoring::res_type_constraint] );

}

void
EnzdesBaseProtocol::disable_constraint_scoreterms(){

	constraint_weights_.clear();
	constraint_weights_[core::scoring::coordinate_constraint] = scorefxn_->weights()[core::scoring::coordinate_constraint];
	constraint_weights_[core::scoring::atom_pair_constraint] = scorefxn_->weights()[core::scoring::atom_pair_constraint];
	constraint_weights_[core::scoring::angle_constraint] = scorefxn_->weights()[core::scoring::angle_constraint];
	constraint_weights_[core::scoring::dihedral_constraint] = scorefxn_->weights()[core::scoring::dihedral_constraint];
	constraint_weights_[core::scoring::res_type_constraint] = scorefxn_->weights()[core::scoring::res_type_constraint];

	scorefxn_->set_weight(core::scoring::coordinate_constraint, 0.0 );
	scorefxn_->set_weight(core::scoring::atom_pair_constraint, 0.0 );
	scorefxn_->set_weight(core::scoring::angle_constraint, 0.0 );
	scorefxn_->set_weight(core::scoring::dihedral_constraint, 0.0 );
	scorefxn_->set_weight(core::scoring::res_type_constraint, 0.0 );

}


bool
EnzdesBaseProtocol::exchange_ligands_in_pose(
	core::pose::Pose & pose,
	bool check_bb_clashes,
	core::scoring::ScoreFunctionCOP scofx
){

	if( !lig_superposition_file_read_ ) read_ligand_superposition_file(basic::options::option[basic::options::OptionKeys::enzdes::change_lig].value() );

	utility::vector1< core::Size > ligs_to_exchange;
	for( core::Size i = 1; i <= pose.total_residue(); ++i ){
		if( pose.residue(i).name3() == res_to_superimpose_.first ) ligs_to_exchange.push_back( i );
	}

	core::conformation::Residue new_res( restype_set_->name_map( res_to_superimpose_.second ), true );

	for( utility::vector1< core::Size >::const_iterator pos_it = ligs_to_exchange.begin();
			 pos_it != ligs_to_exchange.end(); ++pos_it )
		{

			pose.replace_residue( *pos_it, new_res, atoms_to_superimpose_on_ );

			//we should probably do this
			pose.update_residue_neighbors();

			// if we don't care for backbone clashes, we are already done
			if( check_bb_clashes ){

				//otherwise, we now have to check for bb_clashes with the new residue

				utility::vector1< core::conformation::ResidueCOP > accepted_rotamers;
				get_non_bb_clashing_rotamers( pose, *pos_it, scofx, accepted_rotamers );

				//DEBUG SHIT
				//core::pose::Pose debugpose = pose;
				//core::Size dei(0);
				//for( utility::vector1< core::conformation::ResidueOP >::iterator rot_it = accepted_rotamers.begin(); rot_it != accepted_rotamers.end(); ++rot_it ){
				//dei++;
				//debugpose.replace_residue( *pos_it, **rot_it, true );
				//debugpose.dump_pdb( "rodebug"+utility::to_string( dei ) );
				//}
				//DEBUG SHIT over

				if( accepted_rotamers.size() == 0 ) return false;

				tr << "There are " << accepted_rotamers.size() << " rotamers that don't clash with the backbone." << std::endl;

				//we should put one of the non_clashing rotamers into the pose, just to make sure...
				pose.replace_residue( *pos_it, *accepted_rotamers[1], true );
			}

			//now we also have to change the remarks
			core::pose::Remarks & remarks = pose.pdb_info()->remarks();
			for( std::vector< core::pose::RemarkInfo >::iterator rem_it = remarks.begin(); rem_it != remarks.end(); ++rem_it ){

				std::string chainA(""), resA(""),chainB(""),resB("");
				core::Size cst_block(0), exgeom_id(0);
				int seqposA(0), seqposB(0);
				if( toolbox::match_enzdes_util::split_up_remark_line( rem_it->value, chainA, resA, seqposA, chainB, resB, seqposB, cst_block, exgeom_id ) ){

					bool line_changed( false );

					if( resA == res_to_superimpose_.first ){
						if( seqposA == (int) *pos_it || seqposA == 0 ) {
							resA = res_to_superimpose_.second;
							line_changed=true;
						}
					}
					else if( resB == res_to_superimpose_.first ){
						if( seqposB == (int) *pos_it || seqposB == 0 ) {
							resB = res_to_superimpose_.second;
							line_changed=true;
						}
					}

					if( line_changed ){
						rem_it->value = toolbox::match_enzdes_util::assemble_remark_line( chainA, resA, seqposA, chainB, resB, seqposB, cst_block, exgeom_id );
					}
				}
			}
		} //lig_positions to exchange


	return true;

} //exchange_ligands_in_pose


core::scoring::ScoreFunctionCOP
EnzdesBaseProtocol::reduced_scorefxn() const {
	return reduced_sfxn_;
}

core::scoring::ScoreFunctionOP
EnzdesBaseProtocol::reduced_scorefxn(){
  return reduced_sfxn_;
}

core::Real
EnzdesBaseProtocol::design_targets_score(
	core::pose::Pose const & pose
) const
{
	core::Real return_val(0.0);
	using namespace core::scoring;

	for( std::set< core::Size >::const_iterator des_it = design_targets_.begin();
			 des_it != design_targets_.end(); ++des_it ){

		return_val += pose.energies().residue_total_energy( *des_it );
	}

	return return_val;
}//design targets score


void
EnzdesBaseProtocol::remap_resid(
	core::pose::Pose const & pose,
		core::id::SequenceMapping const & //smap
)
{
	//cst_io_->remap_resid( smap );
	design_targets( pose );
}



void
EnzdesBaseProtocol::generate_explicit_ligand_rotamer_poses(
	core::pose::Pose const & orig_pose,
	utility::vector1< core::pose::PoseOP > & ligrot_poses,
	core::scoring::ScoreFunctionCOP scofx
)
{

	ligrot_poses.clear();

	//std::cerr << "staring generate explictligrot function, looking for residues of type " << basic::options::option[ basic::options::OptionKeys::enzdes::process_ligrot_separately ].value() << std::endl;

	for( core::Size i = 1; i <= orig_pose.total_residue(); ++ i ){

		if( orig_pose.residue( i ).name3() != basic::options::option[ basic::options::OptionKeys::enzdes::process_ligrot_separately ].value() ) continue;

		utility::vector1< core::conformation::ResidueCOP > accepted_rotamers;
		get_non_bb_clashing_rotamers( orig_pose, i, scofx, accepted_rotamers );

		for( utility::vector1< core::conformation::ResidueCOP >::const_iterator rot_it = accepted_rotamers.begin();
				 rot_it != accepted_rotamers.end(); ++rot_it ){

			core::pose::PoseOP lig_pose = new core::pose::Pose( orig_pose );
			lig_pose->replace_residue( i, **rot_it, true );

			ligrot_poses.push_back( lig_pose );

		}
		//std::cerr << "There are " << accepted_rotamers.size() << " non_bb clashing residues at positon " << i << std::endl;
	}

} //generate_explicit_ligand_rotamer_poses

void
EnzdesBaseProtocol::read_ligand_superposition_file( std::string filename )
{

	lig_superposition_file_read_ = false;
	atoms_to_superimpose_on_.clear();

	bool switch_info_found( false );

	utility::io::izstream filedata( filename.c_str() );
	std::istringstream line_stream;
	std::string line("");

	if ( !filedata ) {
		std::cerr << "ERROR:: Unable to open ligand superposition info file: "
				<< filename << std::endl;
		std::exit( 1 );
	}

	while( !filedata.eof() ) {
		std::string key("");
		getline(filedata,line);
		line_stream.clear();
		line_stream.str(line);
		line_stream >> key;

		if( key == "SWITCH_NAME" ){

			line_stream >> res_to_superimpose_.first >> res_to_superimpose_.second;
			switch_info_found = true;

		}

		if( key == "SUPERIMPOSE" ){

			std::string buf1(""), buf2("");
			line_stream >> buf1 >> buf2;
			atoms_to_superimpose_on_.push_back( std::pair< std::string, std::string > ( buf1, buf2 ) );

		}
	}

	if( !switch_info_found ){
		utility_exit_with_message("Error: "+filename+" does not specify which residues to switch!");
	}

	if( atoms_to_superimpose_on_.size() < 3 ){
		utility_exit_with_message("Error: "+filename+" specifies less than 3 atoms to superimpose on, unambiguous superposition not possible");
	}

	tr << "read superposition info from " << filename << " ...";
	lig_superposition_file_read_ = true;

} //read_ligand_superposition_file

void
EnzdesBaseProtocol::set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn ){

	scorefxn_ = scorefxn->clone();

}

utility::vector1< core::Size >
EnzdesBaseProtocol::rb_min_jumps() const{ return rb_min_jumps_; }

void
EnzdesBaseProtocol::rb_min_jumps( utility::vector1< core::Size > const v ){
	rb_min_jumps_ = v; }

}//namespace enzdes
}//namespace protocols
