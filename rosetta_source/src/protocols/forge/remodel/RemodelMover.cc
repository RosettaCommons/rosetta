// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/remodel/RemodelMover.cc
/// @brief
/// @author Possu Huang (possu@u.washington.edu)
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/remodel/RemodelMover.hh>
#include <protocols/forge/remodel/RemodelLoopMover.hh>
#include <protocols/forge/remodel/RemodelEnzdesCstModule.fwd.hh>
#include <protocols/protein_interface_design/dock_design_filters.hh>

// package headers
#include <protocols/forge/build/BuildInstruction.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/SegmentInsert.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/methods/chainbreak_eval.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/forge/methods/util.hh>
#include <protocols/enzdes/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <core/pack/pack_rotamers.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/moves/symmetry/SetupForSymmetryMover.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

// project headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/io/pdb/file_data.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/Remarks.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh> // for pdbinfo
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <protocols/jumping/Dssp.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/LoopMover_CCD.hh>
#include <protocols/loops/LoopMover_KIC.hh>
#include <protocols/moves/PackRotamersMover.hh>
#include <protocols/toolbox/pose_manipulation.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.hh>
#include <protocols/toolbox/task_operations/PreventChainFromRepackingOperation.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>
#include <protocols/forge/remodel/RemodelDesignMover.hh>
#include <protocols/forge/remodel/RemodelAccumulator.hh>
#include <protocols/forge/remodel/RemodelEnzdesCstModule.hh>
#include <protocols/viewer/viewers.hh>

//test
#include <protocols/protein_interface_design/movers/DockAndRetrieveSidechains.hh>
#include <protocols/moves/ReturnSidechainMover.hh>
#include <protocols/moves/SwitchResidueTypeSetMover.hh>

// C++ headers
#include <utility>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>


#if (defined WIN32) && (defined WIN_PYROSETTA)
	#include <utility/tag/Tag.hh>
#endif

namespace protocols {
namespace forge {
namespace remodel {


static basic::Tracer TR( "protocols.forge.remodel.RemodelMover" );


/// @brief default constructor
RemodelMover::RemodelMover() :
	Super( "RemodelMover" ),
//	use_fullmer_( false ),
//	use_sequence_bias_( false ),
	max_linear_chainbreak_( 0.07 ),
	//centroid_loop_mover_str_( "quick_ccd" ),
	centroid_loop_mover_str_( "RemodelLoopMover" ),
	redesign_loop_neighborhood_( false ),
	dr_cycles_( basic::options::option[basic::options::OptionKeys::remodel::dr_cycles] ),
	centroid_sfx_( core::scoring::ScoreFunctionFactory::create_score_function( "remodel_cen" )),
	fullatom_sfx_( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::STANDARD_WTS, core::scoring::SCORE12_PATCH ) )
{
	register_user_options();
	if (basic::options::option[basic::options::OptionKeys::packing::soft_rep_design]){
		TR << "SWITCHING FULLATOM FUNCITON TO SOFT_REP_DESIGN" << std::endl;
		fullatom_sfx_ = core::scoring::ScoreFunctionFactory::create_score_function(core::scoring::SOFT_REP_DESIGN_WTS);
	}
}

void
RemodelMover::register_user_options()
{
	using namespace basic::options;
	using namespace core::scoring;

	//set optional weights
	if ( option[OptionKeys::remodel::hb_lrbb].user() ){
		centroid_sfx_-> set_weight(hbond_lr_bb, option[OptionKeys::remodel::hb_lrbb]);
	}

	if ( option[OptionKeys::remodel::hb_srbb].user() ){
		centroid_sfx_-> set_weight(hbond_sr_bb, option[OptionKeys::remodel::hb_srbb]);
	}
	if ( option[OptionKeys::remodel::rg].user() ){
		centroid_sfx_-> set_weight(rg, option[OptionKeys::remodel::rg]);
	}

	if ( option[OptionKeys::remodel::rsigma].user() ){
		centroid_sfx_-> set_weight(rsigma, option[OptionKeys::remodel::hb_lrbb]);
	}
	if ( option[OptionKeys::remodel::ss_pair].user() ){
		centroid_sfx_-> set_weight(ss_pair, option[OptionKeys::remodel::ss_pair]);
	}
}


/// @brief copy constructor
RemodelMover::RemodelMover( RemodelMover const & rval ) :
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	//manager_( rval.manager_ ),
	design_info_( rval.design_info_ ),
	//use_fullmer_( rval.use_fullmer_ ),
	//use_sequence_bias_( rval.use_sequence_bias_ ),
	max_linear_chainbreak_( rval.max_linear_chainbreak_ ),
	centroid_loop_mover_str_( rval.centroid_loop_mover_str_ ),
	redesign_loop_neighborhood_( rval.redesign_loop_neighborhood_ ),
	//resfile_( rval.resfile_ ),
	dr_cycles_( rval.dr_cycles_ ),
	centroid_sfx_( rval.centroid_sfx_->clone() ),
	fullatom_sfx_( rval.fullatom_sfx_->clone() )
{
	if ( rval.vlb_.get() ) {
		vlb_ = new VarLengthBuild( *rval.vlb_ );
	}
}


/// @brief default destructor
RemodelMover::~RemodelMover() {}


/// @brief clone this object
RemodelMover::MoverOP RemodelMover::clone() {
	return new RemodelMover( *this );
}


/// @brief create this type of object
RemodelMover::MoverOP RemodelMover::fresh_instance() {
	return new RemodelMover();
}


/// @brief the centroid level score function, default "remodel_cen"
RemodelMover::ScoreFunction const & RemodelMover::centroid_scorefunction() const {
	return *centroid_sfx_;
}


/// @brief the full-atom level score function, default score12
RemodelMover::ScoreFunction const & RemodelMover::fullatom_scorefunction() const {
	return *fullatom_sfx_;
}


/// @brief add instruction to the manager of this RemodelMover (no copy)
/// @param[in] bi BuildInstruction
/// @param[in] aa_during_design_refine The allowed amino acid sequence
///  during design.  Only applicable to BuildInstructions like
///  SegmentRebuild and SegmentInsert.  Make sure the length of this
///  string matches up properly.  Default empty string.
/*void RemodelMover::add_instruction(
	BuildInstructionOP bi,
	String const & aa_during_design_refine
)
{
	manager_.add( bi );
	if ( !aa_during_design_refine.empty() ) {
		design_info_.push_back( std::make_pair( bi->original_interval(), aa_during_design_refine ) );
	}

	// additional instruction means we'll need a new re-init the VLB, so
	// go ahead and drop the existing one
	vlb_ = 0;
}
*/

/// @brief create directed dependency between two instructions
/*
void RemodelMover::create_directed_dependency(
	BuildInstructionOP u,
	BuildInstructionOP v
)
{
	manager_.create_directed_dependency( u, v );
}
*/

/// @brief set the centroid level score function
void RemodelMover::centroid_scorefunction( ScoreFunction const & sfx ) {
	centroid_sfx_ = sfx.clone();
}


/// @brief set the centroid level score function
void RemodelMover::centroid_scorefunction( ScoreFunctionOP sfx ) {
	centroid_sfx_ = sfx->clone();
}


/// @brief set the full-atom level score function
void RemodelMover::fullatom_scorefunction( ScoreFunction const & sfx ) {
	fullatom_sfx_ = sfx.clone();
}


/// @brief set the full-atom level score function
void RemodelMover::fullatom_scorefunction( ScoreFunctionOP sfx ) {
	fullatom_sfx_ = sfx->clone();
}


/// @brief apply defined moves to given Pose
void RemodelMover::apply( Pose & pose ) {
	using core::pose::metrics::CalculatorFactory;
	using basic::MetricValue;
	using protocols::jumping::Dssp;
	using protocols::moves::MS_SUCCESS;
	using protocols::moves::FAIL_DO_NOT_RETRY;
	using protocols::moves::FAIL_RETRY;
	using protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator;
	using protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator;
	using namespace core::scoring;
	using namespace protocols::protein_interface_design;
	using namespace basic::options;

#if defined GL_GRAPHICS
  protocols::viewer::add_conformation_viewer( pose.conformation(), "Remodel" );
#endif

	//store the starting pose for KIC confirmation RMSD calculation
	native_pose_ = pose;

	// assign secondary structure
	Dssp dssp( pose );

	protocols::forge::remodel::RemodelData remodel_data;
	protocols::forge::remodel::WorkingRemodelSet working_model;
	// read blueprint
	remodel_data.getLoopsToBuildFromFile();

	TR << pose.total_residue() << std::endl;
	ObjexxFCL::FArray1D_char dsspSS( pose.total_residue() );
	dssp.dssp_reduced(dsspSS);
	TR << "input PDB dssp assignment:" << std::endl;

	for (Size i = 1; i<= pose.total_residue(); i++){
		TR << dsspSS(i);
	}
	TR << std::endl;

	remodel_data.updateWithDsspAssignment( dsspSS );
	dssp.insert_ss_into_pose( pose );


  // insertion
	if (option[OptionKeys::remodel::domainFusion::insert_segment_from_pdb].user()){
		TR << "INSERT SEGMENT FROM PDB" << std::endl;
		remodel_data.collectInsertionPose();
	}

	// score the pose first
	core::scoring::ScoreFunctionOP sfx = core::scoring::ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH );
	(*sfx)( pose );

	working_model.workingSetGen( pose, remodel_data );

	remodel_data_ = remodel_data; // will use the movemap for natro definition later
	working_model_=working_model;

/*  DEBUG
	std::set<core::Size> up = working_model.manager.undefined_positions();
  for ( std::set<core::Size>::iterator i = up.begin(); i!=up.end(); i++){
	  TR << *i << std::endl;
		}
     std::set<core::Size> uup = working_model.manager.union_of_intervals_containing_undefined_positions();
	for ( std::set<core::Size>::iterator i = uup.begin(); i!=uup.end(); i++){
		TR << *i <<  " UUP" <<  std::endl;
	}
*/
//	Pose testArc;
//	testArc = pose;
if (working_model.manager.size()!= 0){
  if (!basic::options::option[basic::options::OptionKeys::remodel::bypass_fragments]){
		working_model.manager.modify(pose);
		}else{
		working_model.manager.dummy_modify(pose.total_residue());
		}
//	protocols::forge::methods::restore_residues(working_model.manager.original2modified(), testArc, pose);
//	pose.dump_pdb("testArcRestore.pdb");
	//testArc=pose;
	manager_ = working_model.manager;
	core::pose::renumber_pdbinfo_based_on_conf_chains(
																												pose,
																												true ,  // fix chain
																												true, // start_from_existing_numbering
																												false, // keep_insertion_code
																												false // rotate_chain_id
	);

}
//	manager_.dummy_modify(testArc.n_residue());
//	core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID, true);
//	core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD, true);
//	protocols::forge::methods::restore_residues(manager_.original2modified(), testArc, pose);
//	pose.dump_pdb("testArcRestore2.pdb");
//	protocols::forge::methods::restore_residues(manager_.original2modified(), testArc, pose);
//	pose.update_residue_neighbors();
//	pose.dump_pdb("testArcRestore3.pdb");
	//testArc.dump_pdb("testArcRestoreSrc3.pdb");
	//protocols::moves::SwitchResidueTypeSetMover to_all_atom( core::chemical::FA_STANDARD);
	//protocols::moves::ReturnSidechainMover recover_sidechains( testArc);
	//to_all_atom.apply(pose);
	//recover_sidechains.apply(pose);
	//pose.dump_pdb("MoverREstore.pdb");

	//initialize symmetry
  if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
      protocols::moves::symmetry::SetupForSymmetryMover pre_mover;
      pre_mover.apply( pose );
			// Remodel assumes chain ID is ' '
			core::pose::PDBInfoOP pdb_info ( pose.pdb_info() );
			for ( Size i=1; i<= pdb_info->nres(); ++i ){
				pdb_info->chain(i,' ');
			}
			pose.pdb_info( pdb_info );
  }

/* moved to VLB
	if ( basic::options::option[basic::options::OptionKeys::enzdes::cstfile].user() ){

		RemodelEnzdesCstModuleOP cstOP = new RemodelEnzdesCstModule(remodel_data);

		//RemodelEnzdesCstModule cst(remodel_data);
		cstOP->use_backbone_only_blocks();
		cstOP->apply(pose);
		cstOP->enable_constraint_scoreterms(centroid_sfx_);
	}
*/
  Size i = basic::options::option[basic::options::OptionKeys::remodel::num_trajectory];
	Size num_traj = i; //need this for checkpointing math
	Size prev_checkpoint = 0;

	RemodelAccumulator accumulator(working_model);

	if (option[OptionKeys::remodel::checkpoint]){
		prev_checkpoint = accumulator.recover_checkpoint();
		if (prev_checkpoint >= i){
			i = 0;
		} else {
			i = i - prev_checkpoint;
		}
	}
 if (working_model.manager.size() != 0){
	// setup calculators
		CalculatorFactory::Instance().remove_calculator( neighborhood_calc_name() );
		CalculatorFactory::Instance().register_calculator(
      neighborhood_calc_name(),
      new NeighborhoodByDistanceCalculator( manager_.union_of_intervals_containing_undefined_positions() )
    );
}

	RemodelDesignMover designMover(remodel_data, working_model, fullatom_sfx_);

	while ( i > 0){
		// do centroid build
		TR << "BUILD CYCLE REMAINING " << i << std::endl;
		core::kinematics::FoldTree originalTree = pose.fold_tree();
		if (working_model.manager.size() != 0){
			if ( !centroid_build( pose, working_model.manager ) ) { // build failed
				set_last_move_status( FAIL_RETRY );
				continue;
			//	return;
			}
		}

		designMover.set_state("stage");

		//handle constraints as soon as centroid is done.  If applying sidechain
		//constraints, replace residue to the right AA right away.
		if ( basic::options::option[basic::options::OptionKeys::enzdes::cstfile].user() ){

			RemodelEnzdesCstModuleOP cstOP = new RemodelEnzdesCstModule(remodel_data);

			//RemodelEnzdesCstModule cst(remodel_data);
			//safety
			pose.remove_constraints();
			//wipe out cst_cache
			protocols::enzdes::get_enzdes_observer( pose ) -> set_cst_cache( NULL );
			//wipe out observer too
			pose.observer_cache().set( core::pose::datacache::CacheableObserverType::ENZDES_OBSERVER, NULL , false);

			//cstOP->remove_constraints_from_pose(pose,true /*keep covalent*/, false /*fail if missing*/);

			cstOP->use_all_blocks();
			cstOP->apply(pose);
			cstOP->enable_constraint_scoreterms(fullatom_sfx_);
			designMover.scorefunction(fullatom_sfx_);
		}

		if (basic::options::option[OptionKeys::remodel::build_disulf].user()){
			utility::vector1<std::pair <Size, Size> > disulf_partners;
			bool disulfPass = false;
			disulfPass = designMover.find_disulfides_in_the_neighborhood(pose, disulf_partners);
			if (disulfPass != true){
				continue;
			}

			for (utility::vector1<std::pair<Size,Size> >::iterator itr = disulf_partners.begin(); itr != disulf_partners.end(); itr++){
				Pose disulf_copy_pose = pose;
				utility::vector1<std::pair <Size, Size> > single_disulf;
				single_disulf.push_back(*itr);
				designMover.make_disulfide(disulf_copy_pose, single_disulf);
				designMover.apply(disulf_copy_pose);
				//for now, accept all disulf build, as it is hard enough to do
				//already.  Accept instead of cst filter?
				//accumulator.apply(disulf_copy_pose);
				if ( basic::options::option[basic::options::OptionKeys::enzdes::cstfile].user() ){
					ScoreTypeFilter const  pose_constraint( fullatom_sfx_, atom_pair_constraint, 10 );
					bool CScore(pose_constraint.apply( disulf_copy_pose ));
					if (!CScore){  // if didn't pass, rebuild
						continue;
					}
					else {
						accumulator.apply(disulf_copy_pose);
					}
				} else {
					accumulator.apply(disulf_copy_pose);
				}
			}
		}
		else {
			designMover.apply(pose);

			if ( basic::options::option[basic::options::OptionKeys::enzdes::cstfile].user() ){
						ScoreTypeFilter const  pose_constraint( fullatom_sfx_, atom_pair_constraint, 10 );
						bool CScore(pose_constraint.apply( pose ));
						if (!CScore){  // if didn't pass, rebuild
							continue;
						}
						else {
							accumulator.apply(pose);
						}
			} else {
				accumulator.apply(pose);
			}
		}

		if (option[OptionKeys::remodel::checkpoint]){
		// debug:
			TR << "writing chkpnt at step " << num_traj-i+prev_checkpoint << std::endl;
			accumulator.write_checkpoint(num_traj-i-prev_checkpoint);
		}
		//restore foldtree
		pose.fold_tree(originalTree);
		i--;
	}

	//take the lowest member and the cluster center
	//accumulator.shrink_cluster();
	std::vector<core::pose::PoseOP> results;
	if (accumulator.cluster_switch()){
		results = accumulator.clustered_best_poses();
		//results = accumulator.clustered_top_poses(option[OptionKeys::remodel::collect_clustered_top]);
	}
	else {
		results = accumulator.contents_in_pose_store();
	}

	//seriously refine the poses
	Size filecount = 1;

	TR << "clustered poses count: " << results.size() << std::endl;
	for(std::vector<core::pose::PoseOP>::iterator it = results.begin(), end= results.end(); it!= end; it++){
		bool bypass_refinement = basic::options::option[basic::options::OptionKeys::remodel::quick_and_dirty].user();
		if (working_model.manager.size() == 0 ) bypass_refinement = true;
		if ( !bypass_refinement ){

	//	std::stringstream SS1;
	//	SS1 << "pre-ref_" << filecount << ".pdb";
	//	(*(*it)).dump_scored_pdb(SS1.str(), *fullatom_sfx_);

			TR << "aggressively refine" << std::endl;
			if (basic::options::option[basic::options::OptionKeys::remodel::use_pose_relax]){
				if (!design_refine_seq_relax(*(*it), designMover)){
					TR << "WARNING: DESIGN REFINE SEQ RELAX FAILED!! (one should never see this)" << std::endl;
					continue;
				}
			}else {
				if (! design_refine(*(*it), designMover)){
					TR << "WARNING: DESIGN REFINE FAILED TO CLOSE STRUCTURE!!" << std::endl;
					continue;
				}
			}
		} else // simple design
		{
			designMover.set_state("finish");
			designMover.apply(*(*it));
		}

		if ( basic::options::option[basic::options::OptionKeys::remodel::run_confirmation].user()){
			if (!confirm_sequence(*(*it))){
				TR << "WARNING: STRUCTURE DID NOT PASS KIC CONFIRMATION!!" << std::endl;
				continue;
			}
		}

		std::stringstream SS;
		SS << filecount << ".pdb";

		//this is to make sure that the final scoring is done with SCORE12
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::STANDARD_WTS, core::scoring::SCORE12_PATCH );

		(*(*it)).dump_scored_pdb(SS.str(), *scorefxn);
		filecount++;
	}


		// setup calculators
		CalculatorFactory::Instance().remove_calculator( neighborhood_calc_name() );
		CalculatorFactory::Instance().register_calculator(
			neighborhood_calc_name(),
			new NeighborhoodByDistanceCalculator( manager_.union_of_intervals_containing_undefined_positions() )
		);


/*
	// do design-refine iteration
	if ( dr_cycles_ > 0 ) {

		if ( !design_refine( pose ) ) { // design-refine failed
			set_last_move_status( FAIL_RETRY );
			return;
		}

	}
*/
	// if we've gotten to this point, then the structure has been
	// built properly
	set_last_move_status( MS_SUCCESS );

	// setup the PoseMetricCalculators and add them to the evaluators in the
	// JobOutputter
	CalculatorFactory::Instance().remove_calculator( loops_buns_polar_calc_name() );
	CalculatorFactory::Instance().remove_calculator( neighborhood_buns_polar_calc_name() );

	CalculatorFactory::Instance().register_calculator(
		loops_buns_polar_calc_name(),
		new BuriedUnsatisfiedPolarsCalculator(
			"default",
			"default",
			manager_.union_of_intervals_containing_undefined_positions()
		)
	);

	MetricValue< std::set< Size > > loops_neighborhood;
	pose.metric( neighborhood_calc_name(), "neighbors", loops_neighborhood );
	CalculatorFactory::Instance().register_calculator(
		neighborhood_buns_polar_calc_name(),
		new BuriedUnsatisfiedPolarsCalculator(
			"default",
			"default",
			loops_neighborhood.value()
		)
	);
}

std::string
RemodelMover::get_name() const {
	return "RemodelMover";
}

bool RemodelMover::centroid_build(
	Pose & pose,
	protocols::forge::build::BuildManager & manager
) {
	manager_ = manager;
	if (basic::options::option[basic::options::OptionKeys::remodel::bypass_fragments]){
		TR << "-=BYPASSING FRAGMENT BUILD (REFINE ONLY) =-" << std::endl;
		return true;
	}
	if (centroid_build( pose )){
	//update external manager
		//manager = manager_;
		return true;
	}
	else {
		return false;
	}
}



/// @brief run the centroid level build stage
/// @return true if loop closed, false otherwise
bool RemodelMover::centroid_build(
	Pose & pose
) {
	using core::scoring::STANDARD_WTS;
	using core::scoring::SCORE12_PATCH;
	using core::scoring::ScoreFunctionOP;
	using core::scoring::ScoreFunctionFactory;
	using namespace basic::options;
	using protocols::moves::MS_SUCCESS;

	using core::util::switch_to_residue_type_set;
	using protocols::forge::methods::restore_residues;
	//using protocols::toolbox::pose_manipulation::construct_poly_uniq_restype_pose;
	using namespace protocols::forge::components;

	// safety, clear the energies object
	pose.energies().clear();

	// make backup Pose for transferring sidechains
	Pose archive_pose = pose;
	Pose modified_archive_pose = archive_pose;
	//manager_.modify( modified_archive_pose );

	// ensure modified_archive_pose is completely full-atom, otherwise mismatch
	// will occur when restoring sidechains at the end of the procedure
	bool mod_ap_is_full_atom = true;
	for ( Size i = 1, ie = modified_archive_pose.n_residue(); mod_ap_is_full_atom && i != ie; ++i ) {
		mod_ap_is_full_atom &= ( modified_archive_pose.residue( i ).residue_type_set().name() == core::chemical::FA_STANDARD );
	}

	if ( !mod_ap_is_full_atom ) {
		core::util::switch_to_residue_type_set( modified_archive_pose, core::chemical::FA_STANDARD );
	}
/*
	// flip to poly-ala-gly-pro-disulf pose, only in the rebuilt segment
	utility::vector1< Size > protein_residues;
	for (std::set<Size>::iterator it=rebuild.begin(), end=rebuild.end(); it != end; it++){
	//for ( Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
		if ( pose.residue( *it ).is_protein() ) {
			protein_residues.push_back( *it );
			TR<< "turning these to ala: " << *it << std::endl;
		}
	}
  TR << "default building restype: " << "ALA" << std::endl;
	construct_poly_uniq_restype_pose( pose, protein_residues, core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map("ALA"), true, true, true );
*/
	// run VLB to build the new section, if no segments have been added/deleted
	// we use the same VLB so that fragment caching works properly
	if ( !vlb_.get() ) {
		vlb_ = new VarLengthBuild( manager_ , remodel_data_ );
	}

	vlb_->scorefunction( centroid_sfx_ );
	vlb_->vall_memory_usage( VLB_VallMemoryUsage::CLEAR_IF_CACHING_FRAGMENTS );
	vlb_->use_fullmer( option[OptionKeys::remodel::use_same_length_fragments] );
//	vlb_->max_linear_chainbreak( max_linear_chainbreak_ );
	vlb_->loop_mover_str( centroid_loop_mover_str_ );
	vlb_->restart_mode(true);

	if ( option[OptionKeys::remodel::use_blueprint_sequence] ) {
		vlb_->new_sequence_override( remodel_data_.sequence );
	}


	vlb_->apply( pose );

	if ( vlb_->get_last_move_status() == MS_SUCCESS ) {

		// record the used manager w/ all mapping info
	//	manager_ = vlb_->manager();

		// safety, clear all the energies before restoring full-atom residues and
		// scoring
		pose.energies().clear();

		// Swap back original sidechains.  At the moment this is a two step process
		// in case any sidechains from SegmentInsert and the like that aren't in the
		// original archive pose need to be transferred.
		restore_residues( modified_archive_pose, pose );
		// since pose is setup modified in RemodelMover, only one step will do
		//restore_residues( manager_.original2modified(), archive_pose, pose );
		// go ahead and score w/ full-atom here; we do this in case there are no
		// design-refine cycles -- it's useful to have e.g. rama in the output
		(*fullatom_sfx_)( pose );

		return true; // loop closed

	} else {

		pose = archive_pose;

	}

	return false; // false if loop not closed
}

bool RemodelMover::design_refine_seq_relax(
	Pose & pose,
	RemodelDesignMover & designMover
)
{
	using core::kinematics::FoldTree;
	using core::pack::task::operation::RestrictResidueToRepacking;
	using core::pack::task::operation::RestrictResidueToRepackingOP;
	using core::pack::task::operation::RestrictToRepacking;
	using core::scoring::STANDARD_WTS;
	using core::scoring::SCORE12_PATCH;
	using core::scoring::ScoreFunctionOP;
	using core::scoring::ScoreFunctionFactory;
	using protocols::forge::build::SegmentInsert;
	using namespace protocols::loops;
	using protocols::loops::Loops;
	using protocols::loops::LoopMover_Refine_CCD;
	using protocols::moves::PackRotamersMover;
	using protocols::toolbox::task_operations::RestrictToNeighborhoodOperation;

	using core::pose::annotated_to_oneletter_sequence;
	using protocols::forge::methods::intervals_to_loops;
	using protocols::forge::methods::linear_chainbreak;
	using protocols::loops::remove_cutpoint_variants;

	typedef protocols::forge::build::BuildManager::Positions Positions;

	// collect new regions/positions
	std::set< Interval > loop_intervals = manager_.intervals_containing_undefined_positions();
	Original2Modified original2modified_interval_endpoints = manager_.original2modified_interval_endpoints();

	// collect loops
	Loops loops = intervals_to_loops( loop_intervals.begin(), loop_intervals.end() );

	protocols::forge::methods::fill_non_loop_cst_set(pose, loops);

	// safety, clear the energies object
	pose.energies().clear();
	ScoreFunctionOP sfx = fullatom_sfx_->clone();
//turning on weights
  sfx->set_weight(core::scoring::coordinate_constraint, 1.0 );
  sfx->set_weight(core::scoring::atom_pair_constraint, 1.0 );
  sfx->set_weight(core::scoring::angle_constraint, 1.0 );
  sfx->set_weight(core::scoring::dihedral_constraint, 1.0 );
  sfx->set_weight(core::scoring::res_type_constraint, 1.0);
	protocols::relax::FastRelax relaxMover(sfx);

	// run design-refine cycle
	for ( Size i = 0; i < dr_cycles_; ++i ) {

		designMover.set_state("finish");
		designMover.apply(pose);
		relaxMover.apply(pose);
	}

//turning off weights
  sfx->set_weight(core::scoring::coordinate_constraint, 0.0 );
  sfx->set_weight(core::scoring::atom_pair_constraint, 0.0 );
  sfx->set_weight(core::scoring::angle_constraint, 0.0 );
  sfx->set_weight(core::scoring::dihedral_constraint, 0.0 );
  sfx->set_weight(core::scoring::res_type_constraint, 0.0);

	(*sfx)( pose );

	return true;
}

/// @brief run the design-refine stage
/// @return currently always true
bool RemodelMover::design_refine(
	Pose & pose,
	RemodelDesignMover & designMover
)
{
	using core::kinematics::FoldTree;
	using core::pack::task::operation::RestrictResidueToRepacking;
	using core::pack::task::operation::RestrictResidueToRepackingOP;
	using core::pack::task::operation::RestrictToRepacking;
	using namespace core::pack::task::operation;
	using namespace protocols::toolbox::task_operations;
	using core::scoring::STANDARD_WTS;
	using core::scoring::SCORE12_PATCH;
	using core::scoring::ScoreFunctionOP;
	using core::scoring::ScoreFunctionFactory;
	using protocols::forge::build::SegmentInsert;
	using namespace protocols::loops;
	using protocols::loops::Loops;
	using protocols::loops::LoopMover_Refine_CCD;
	using protocols::moves::PackRotamersMover;
	using protocols::toolbox::task_operations::RestrictToNeighborhoodOperation;

	using core::pose::annotated_to_oneletter_sequence;
	using protocols::forge::methods::intervals_to_loops;
	using protocols::forge::methods::linear_chainbreak;
	using protocols::loops::remove_cutpoint_variants;

	typedef protocols::forge::build::BuildManager::Positions Positions;

	// collect new regions/positions
	std::set< Interval > loop_intervals = manager_.intervals_containing_undefined_positions();
	Original2Modified original2modified_interval_endpoints = manager_.original2modified_interval_endpoints();

	// collect loops
	Loops loops = intervals_to_loops( loop_intervals.begin(), loop_intervals.end() );

	// refine Mover used doesn't setup a fold tree, so do it here
	//FoldTree loop_ft = protocols::forge::methods::fold_tree_from_loops( pose, loops );
	FoldTree loop_ft;
	protocols::loops::fold_tree_from_loops( pose, loops, loop_ft, true /*term cut*/);

	// save original fold tree
	FoldTree original_ft = pose.fold_tree();

	// define the score function
	ScoreFunctionOP sfx = fullatom_sfx_->clone();

	// setup the refine TaskFactory
	TaskFactoryOP refine_tf = generic_taskfactory();
	refine_tf->push_back( new RestrictToNeighborhoodOperation( neighborhood_calc_name() ) );
	refine_tf->push_back( new RestrictToRepacking() );

	// safety, clear the energies object
	pose.energies().clear();

	// run design-refine cycle
	for ( Size i = 0; i < dr_cycles_; ++i ) {

		// design the new section
	//	PackRotamersMover design( sfx );
	//	design.task_factory( design_tf );
	//	design.apply( pose );
		designMover.set_state("finish");
		designMover.apply(pose);

		// set loop topology
		pose.fold_tree( loop_ft );


		if (!basic::options::option[basic::options::OptionKeys::remodel::swap_refine_confirm_protocols].user()){
			// refine the new section
			LoopMover_Refine_CCD refine( loops, sfx );
			core::kinematics::MoveMap combined_mm;

			combined_mm.import(remodel_data_.natro_movemap_);
			combined_mm.import( manager_.movemap() );

			//modify task to accept NATRO definition

			utility::vector1<core::Size> natroPositions;
			for (int i = 1; i<= pose.total_residue(); i++){
				if (remodel_data_.natro_movemap_.get_chi(i) == false){
					std::cout << "NATRO for position: " << i << std::endl;
					natroPositions.push_back(i);
				}
			}
			OperateOnCertainResiduesOP natroRes = new OperateOnCertainResidues;
			natroRes->residue_indices( natroPositions );
			natroRes->op( new PreventRepackingRLT );
			refine_tf->push_back( natroRes );

			refine.false_movemap( combined_mm );
			refine.set_task_factory( refine_tf );
			refine.apply( pose );
		}
		else {
			LoopMover_Refine_KIC KIC(loops);
			KIC.apply(pose);
		}

		// remove cutpoint variants -- shouldn't this happen at the end
		// of the refine Mover?
		remove_cutpoint_variants( pose );

		// set original topology
		pose.fold_tree( original_ft );
	//debug
//	std::stringstream SS;
//	SS << "RefineStage" << i << ".pdb";
//	pose.dump_pdb(SS.str());
	}

	// must score one last time since we've removed variants and set
	// new topology, otherwise component energies not correct for
	// e.g. structure output
	(*sfx)( pose );


	// evaluate all chainbreaks using linear chainbreak
	bool cbreaks_pass = true;
	for ( Loops::const_iterator l = loops.begin(), le = loops.end(); l != le && cbreaks_pass; ++l ) {
		if ( l->cut() > 0 ) {
			Real const c = linear_chainbreak( pose, l->cut() );
			TR << "design_refine: final chainbreak = " << c << std::endl;
			cbreaks_pass = c <= max_linear_chainbreak_;
		}
	}

	//return cbreaks_pass;
	return true; //FOR NOW!!  change me back!
}

bool RemodelMover::confirm_sequence(core::pose::Pose & pose ) {
	using core::kinematics::FoldTree;
	using namespace protocols::loops;
	using namespace protocols::forge::methods;
 	using protocols::forge::build::Interval;
	using protocols::toolbox::task_operations::RestrictToNeighborhoodOperation;
	using core::pack::task::operation::RestrictToRepacking;

	core::pose::Pose archive_pose = pose;  //for rmsd

	std::set< Interval > loop_intervals = manager_.intervals_containing_undefined_positions();
	//pose.dump_pdb("pre_KICpose.pdb");

  // collect loops
  Loops confirmation_loops = intervals_to_confirmation_loops( loop_intervals.begin(), loop_intervals.end(), pose.total_residue() );

 // refine Mover used doesn't setup a fold tree, so do it here
  FoldTree loop_ft;
	protocols::loops::fold_tree_from_loops( pose, confirmation_loops, loop_ft, true );
	TR << "confirmation loops tree" << loop_ft << std::endl;

  // save original fold tree
  FoldTree original_ft = pose.fold_tree();
	// switch over to new tree
	pose.fold_tree(loop_ft);

	//LoopMover_Refine_KIC KIC(confirmation_loops);

	TR << "fold tree entering confirmation: " << pose.fold_tree() << std::endl;

	//KIC.apply(pose);

		if (basic::options::option[basic::options::OptionKeys::remodel::swap_refine_confirm_protocols].user()){
		TR << "REFINE USING CCD" << std::endl;
			// refine the new section
			// setup the refine TaskFactory
			//
			//	protocols::forge::remodel::RemodelLoopMover scramble_mover(confirmation_loops);
			//	scramble_mover.randomize_stage(pose);

	TaskFactoryOP refine_tf = generic_taskfactory();
	refine_tf->push_back( new RestrictToNeighborhoodOperation( neighborhood_calc_name() ) );
	refine_tf->push_back( new RestrictToRepacking() );

	LoopMover_Refine_CCD refine( confirmation_loops, fullatom_sfx_ );
			core::kinematics::MoveMap combined_mm;

			combined_mm.import(remodel_data_.natro_movemap_);
			combined_mm.import( manager_.movemap() );

			refine.false_movemap( combined_mm );
			refine.set_task_factory( refine_tf );
			refine.apply( pose );
		}
		else {
		TR << "REFINE USING KIC" << std::endl;
			LoopMover_Refine_KIC KIC(confirmation_loops);
			KIC.apply(pose);
		}

	//reset to original foldtree
	pose.fold_tree(original_ft);

	//pose.dump_pdb("post_KICpose.pdb");

	//rmsd_calculation:

	core::Real sum_sd = 0;
	core::Real sum_sd_native = 0;
	core::Real sum_sd_archive2native=0;
	core::Size atom_count = 0;

	for ( Loops::iterator it = confirmation_loops.v_begin(), end = confirmation_loops.v_end(); it!=end; it++) {
					for (core::Size i = it->start(); i <= it->stop(); ++i){
						core::Real dist_squared = (pose.residue(i).xyz( "CA" ) - archive_pose.residue(i).xyz( "CA" ) ).length_squared();
						core::Real dist_squared_native = (pose.residue(i).xyz( "CA" ) - native_pose_.residue(i).xyz( "CA" ) ).length_squared();
						core::Real dist_squared_archive2native = (archive_pose.residue(i).xyz( "CA" ) - native_pose_.residue(i).xyz( "CA" ) ).length_squared();
						sum_sd = sum_sd + dist_squared;
						sum_sd_native = sum_sd_native + dist_squared_native;
						sum_sd_archive2native = sum_sd_archive2native + dist_squared_archive2native;
						atom_count++;
						if(false){
									std::cout << " Backbone atom= " <<  "CA res: " << i  << " dist_squared= " << dist_squared << std::endl;
									std::cout << " Backbone atom= " <<  "CA res: " << i  << " dist_squared_native= " << dist_squared_native << std::endl;
									std::cout << " Backbone atom= " <<  "CA res: " << i  << " dist_squared_archive2native= " << dist_squared_archive2native << std::endl;
							}
					}
	}
	sum_sd = sum_sd/atom_count;
	sum_sd_native = sum_sd_native/atom_count;
	sum_sd_archive2native = sum_sd_archive2native/atom_count;

	core::Real rmsd = sqrt(sum_sd);
	core::Real rmsd_native = sqrt(sum_sd_native);
	core::Real rmsd_archive2native = sqrt(sum_sd_archive2native);

	core::pose::PDBInfoOP temp_pdbinfo = pose.pdb_info();

	core::pose::RemarkInfo remark;
	remark.value = "KIC confirmation RMSD: " + utility::to_string(rmsd) + " to native RMSD: " + utility::to_string(rmsd_native) ;
	temp_pdbinfo->remarks().push_back( remark );

  remark.value = " ARCHIVE2NATIVE RMSD: " + utility::to_string(rmsd_archive2native);

	temp_pdbinfo->remarks().push_back( remark );

	pose.pdb_info(temp_pdbinfo);

	TR << "RMSD of KIC confirmation: " << rmsd << std::endl;
	TR << "RMSD of KIC confirmation to native: " << rmsd_native << std::endl;
	TR << "RMSD of ARCHIVE to NATIVE: " << rmsd_archive2native << std::endl;

	 //currently the confirmation stage is not setup as filter so always return
	 //true
	if (rmsd <= 1){
		return true;
	}
	else {
		return true; //for now CHANGE IT BACK!!
	}
}


/// @brief return a TaskFactory useable as a starting point for either
///  design or refinement
RemodelMover::TaskFactoryOP RemodelMover::generic_taskfactory() {
	using core::pack::task::operation::IncludeCurrent;
	using core::pack::task::operation::InitializeFromCommandline;
	using core::pack::task::operation::ReadResfile;
	using core::pack::task::operation::ReadResfileOP;
	using core::pack::task::TaskFactory;
	using core::pack::task::operation::NoRepackDisulfides;

	TaskFactoryOP tf = new TaskFactory();

	tf->push_back( new InitializeFromCommandline() ); // also inits -ex options
	tf->push_back( new IncludeCurrent() ); // enforce keeping of input sidechains
	tf->push_back( new NoRepackDisulfides() );

	// load resfile op only if requested
/*	if ( !resfile_.empty() ) {
		ReadResfileOP rrf = new ReadResfile();
		rrf->filename( resfile_ );
		tf->push_back( rrf );
	}
*/
	return tf;
}


/// @brief process a continuous design string, adding appropriate operations
///  to the TaskFactory
void RemodelMover::process_continuous_design_string(
	Interval const & original_interval,
	String const & design_str,
	Original2Modified const & original2modified_interval_endpoints,
	TaskFactoryOP design_tf
)
{
	using core::pack::task::operation::RestrictAbsentCanonicalAAS;

	using core::chemical::aa_from_oneletter_code;

	Size const offset = original2modified_interval_endpoints.find( original_interval.left )->second;
	for ( Size i = 0, ie = design_str.length(); i < ie; ++i ) {
		utility::vector1< bool > allowed_aa_types( 20, false );

		switch ( design_str.at( i ) ) {
			case 's': // surface case, no CFWY
				allowed_aa_types = allowed_surface_aa();
				break;
			case '.': // protocol default design
				continue;
			default: // regular case, single aa type
				allowed_aa_types[ aa_from_oneletter_code( design_str.at( i ) ) ] = true;
				break;
		}

		design_tf->push_back( new RestrictAbsentCanonicalAAS( i + offset, allowed_aa_types ) );
	}
}


/// @brief process a design string containing an insert, adding appropriate
///  operations to the TaskFactory
void RemodelMover::process_insert_design_string(
	Interval const & original_interval,
	String const & design_str,
	Original2Modified const & original2modified_interval_endpoints,
	TaskFactoryOP design_tf
)
{
	using core::pack::task::operation::RestrictAbsentCanonicalAAS;
	using core::pack::task::operation::RestrictResidueToRepacking;
	using core::pack::task::operation::RestrictResidueToRepackingOP;
	using protocols::forge::build::Interval;
	using protocols::forge::build::SegmentInsert;

	using core::chemical::aa_from_oneletter_code;

	char const insert_char = SegmentInsert::insertion_char();

	// Figure out the number of residues in each section.
	Interval const interval(
		original2modified_interval_endpoints.find( original_interval.left )->second,
		original2modified_interval_endpoints.find( original_interval.right )->second
	);

	Size const insert_char_idx = design_str.find( insert_char );
	Size const left_nres = insert_char_idx;
	Size const right_nres = design_str.size() - left_nres - 1;
	Size const insert_nres = interval.length() - left_nres - right_nres;

	// Make setup easy by building a new design string to expand the
	// insertion character into a series of the insertion character
	// the size of the insert.
	String aa = design_str;
	aa.replace( insert_char_idx, 1, insert_nres, insert_char );

	// setup TaskOperations
	RestrictResidueToRepackingOP repack_op = new RestrictResidueToRepacking();

	Size const left_offset = interval.left;
	for ( Size i = 0, ie = aa.size(); i < ie; ++i ) {
		utility::vector1< bool > allowed_aa_types( 20, false );

		if ( aa.at( i ) == insert_char ) { // repack only
			repack_op->include_residue( i + left_offset );
			continue;
		} else if ( aa.at( i ) == 's' ) { // surface case, no CFWY
			allowed_aa_types = allowed_surface_aa();
		} else if ( aa.at( i ) == '.' ) { // protocol default design
			continue;
		} else { // regular case, single aa type
			allowed_aa_types[ aa_from_oneletter_code( aa.at( i ) ) ] = true;
		}

		design_tf->push_back( new RestrictAbsentCanonicalAAS( i + left_offset, allowed_aa_types ) );
	}

	design_tf->push_back( repack_op );
}


/// @brief return a boolean vector specifying allowed a.a. when designing
///  on the surface
utility::vector1< bool > const & RemodelMover::allowed_surface_aa() {
	using core::chemical::aa_from_oneletter_code;

	static String surface_aa = "ADEGHIKLMNPQRSTV";
	static utility::vector1< bool > v( 20, false );

	for ( Size i = 0, ie = surface_aa.length(); i < ie; ++i ) {
		v[ aa_from_oneletter_code( surface_aa.at( i ) ) ] = true;
	}

	return v;
}


} // namespace remodel
} // namespace forge
} // namespace protocols
