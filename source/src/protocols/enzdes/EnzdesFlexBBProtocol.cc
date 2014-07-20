// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/enzdes/EnzdesFlexBBProtocol.cc
///
/// @brief
/// @author Florian Richter


#include <protocols/enzdes/EnzdesFlexBBProtocol.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesLoopsFile.hh>
// AUTO-REMOVED #include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/enzdes/enzdes_util.hh>
#include <protocols/flexpack/rotamer_set/FlexbbRotamerSets.hh>
#include <protocols/flexpack/FlexPacker.hh>

#include <core/fragment/BBTorsionAndAnglesSRFD.hh>
#include <core/fragment/FragData.hh>
#include <core/chemical/ResidueType.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/fragment/util.hh>
#include <core/graph/graph_util.hh> //for deleting edges from graph in calculating lig part sums
#include <basic/options/option.hh>
#include <core/pack/interaction_graph/PDInteractionGraph.hh> //for calculating lig part sums
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/IGEdgeReweightContainer.hh>
#include <core/pack/task/TaskFactory.hh> //task shit
#include <core/pack/pack_rotamers.hh>
#include <core/pack/packer_neighbors.hh> //for calculating lig part sums
#include <core/pack/rotamer_set/RotamerSets.hh> //for calculating lig part sums
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pose/PDB_Info.hh> //for getting pdb name
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/id/SequenceMapping.hh>
//#include <core/scoring/constraints/CoordinateConstraint.hh>
//#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <protocols/simple_filters/ScoreCutoffFilter.hh> // for filtering kinematic mover
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh> //input file reading
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/moves/MonteCarlo.hh>
// AUTO-REMOVED #include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <protocols/toolbox/IGEdgeReweighters.hh>
// AUTO-REMOVED #include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/TenANeighborGraph.hh>

#include <numeric/random/random.hh>

#include <basic/Tracer.hh>
#include <basic/prof.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/sort_predicates.hh>
#include <utility/file/file_sys_util.hh> //checking file existence for loop pdb reading

// Numeric headers
#include <numeric/statistics/functions.hh>

//C++ Headers
// AUTO-REMOVED #include <ctime>
//#include <math.h> //std::min ?
// option key includes

#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/util.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols{
namespace enzdes{

static basic::Tracer tr("protocols.enzdes.EnzdesFlexBBProtocol");

EnzdesFlexBBProtocol::~EnzdesFlexBBProtocol() {}

EnzdesFlexBBProtocol::EnzdesFlexBBProtocol()
	: EnzdesBaseProtocol(),
		enz_loops_file_(NULL),
		mc_kt_low_(basic::options::option[basic::options::OptionKeys::enzdes::mc_kt_low] ),
		mc_kt_high_(basic::options::option[basic::options::OptionKeys::enzdes::mc_kt_high] ),
		brub_min_atoms_( basic::options::option[basic::options::OptionKeys::backrub::min_atoms] ),
		brub_max_atoms_( basic::options::option[basic::options::OptionKeys::backrub::max_atoms] ),
		loop_ensemble_size_( basic::options::option[basic::options::OptionKeys::enzdes::single_loop_ensemble_size]),
		loopgen_trials_( basic::options::option[basic::options::OptionKeys::enzdes::loop_generator_trials] )
{
	flex_regions_.clear();

	//if( loop_ensemble_size_ % 2 != 0 ) loop_ensemble_size_++;
	//loop_ensemble_size_ /= 2; //divide by 2 because brub generates two reasonable structs for each run

	//brub_mover_ = new protocols::backrub::BackrubMover();
	//apparently it's better to initialize the backrub mover new for every structure
	brub_mover_ = NULL;
	brub_pivot_atoms_.push_back("CA");

	/*
	if( basic::options::option[ basic::options::OptionKeys::enzdes::enz_loops_file ].user() ){

		enz_loops_file_ = new toolbox::match_enzdes_util::EnzdesLoopsFile();

		if( !enz_loops_file_->read_loops_file( basic::options::option[ basic::options::OptionKeys::enzdes::enz_loops_file ] ) ){
			utility_exit_with_message("Reading enzdes loops file failed");
		}
	}
	*/
	if ( ! basic::options::option[ basic::options::OptionKeys::enzdes::kic_loop_sampling ] ) {
		if( scorefxn_->has_zero_weight( core::scoring::mm_bend ) ){ scorefxn_->set_weight( core::scoring::mm_bend, 1.0 ); }

		if( reduced_sfxn_->has_zero_weight( core::scoring::mm_bend ) ){ reduced_sfxn_->set_weight( core::scoring::mm_bend, 1.0 ); }
	}

	if( scorefxn_->has_zero_weight( core::scoring::rama ) ){

		if( reduced_sfxn_->has_zero_weight( core::scoring::rama ) ) {
			scorefxn_->set_weight( core::scoring::rama, 1.0 );
		}
		else scorefxn_->set_weight( core::scoring::rama, reduced_sfxn_->weights()[core::scoring::rama] );

	}

	if( reduced_sfxn_->has_zero_weight( core::scoring::rama ) ){ reduced_sfxn_->set_weight( core::scoring::rama, 1.0 ); }
}

void
EnzdesFlexBBProtocol::apply(
	core::pose::Pose & pose
){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::moves;
	using namespace core::pack::task;

	//if ( ! basic::options::option[ basic::options::OptionKeys::enzdes::kic_loop_sampling ] ) {
	//	brub_mover_ = new protocols::backrub::BackrubMover();
	//	brub_mover_->set_native_pose( & pose );
	//}
	//kinematic_mover_ = new protocols::loops::kinematic_closure::KinematicMover();

	pack_region_ala_pose_ = pose;

	flex_regions_.clear();
	fragment_counters_.clear();

	if( ! basic::options::option[basic::options::OptionKeys::in::file::native].user() ){

		core::pose::PoseOP natpose = new core::pose::Pose( pose );
		(*scorefxn_)( *natpose );
		this->set_native_pose( natpose );
	}


	// Scoring function already set up by superclass
	//tr.Info << "starting apply function..." << std::endl;

	//set up constraints (read cstfile, do mapping, etc, then add to pose)
	if( basic::options::option[basic::options::OptionKeys::enzdes::cstfile].user() ){
		enable_constraint_scoreterms();
		setup_enzdes_constraints( pose, false );
	}

	//create packer task (read resfile, etc)
	PackerTaskOP fixbb_pack_task;

	tr.Info << "Done setting up the task and constraints... " << std::endl;
	//score pose to make sure everything is initialised correctly
	(*scorefxn_)( pose );

	//cst opt stage, if demanded
	if(basic::options::option[basic::options::OptionKeys::enzdes::cst_opt]){
		fixbb_pack_task =  create_enzdes_pack_task( pose );
		tr.Info << "starting cst_opt minimization..." << std::endl;
		cst_minimize(pose, fixbb_pack_task, true);
		(*scorefxn_)( pose );
		tr.Info << "done cst_opt minimization." << std::endl;
	}


	if(basic::options::option[basic::options::OptionKeys::enzdes::cst_design]){

		fixbb_pack_task = create_enzdes_pack_task( pose ); //make a new task in case the ligand has moved a lot
		//setup_pose_metric_calculators( pose );

		determine_flexible_regions( pose, fixbb_pack_task );

		//we should also make those residues designable that have been deemed flexible
		PackerTaskOP design_pack_task_template = modified_task( pose, *fixbb_pack_task );

		//make a polyalanine copy of the pose
		utility::vector1< core::Size >positions_to_replace;
		utility::vector1< core::Size >all_pack_positions;
		for(core::Size i = 1, i_end = pose.total_residue(); i <= i_end; ++i) {
			if( design_pack_task_template->pack_residue(i) && !pose.residue( i ).is_ligand() ) {
				all_pack_positions.push_back( i );
				if( ! is_catalytic_position( pose, i ) ) positions_to_replace.push_back( i );
			}

		}
		core::pose::PoseOP poly_ala_pose = new core::pose::Pose(pose);

		protocols::toolbox::pose_manipulation::construct_poly_ala_pose( *poly_ala_pose, positions_to_replace, true, true, true );
		protocols::toolbox::pose_manipulation::construct_poly_ala_pose( pack_region_ala_pose_, all_pack_positions, true, true, true );

		core::scoring::EnergyMap const cur_emap = scorefxn_->weights();

		scorefxn_->set_weight( core::scoring::fa_sol, 0.0);


		if( !this->recover_loops_from_file( *poly_ala_pose ) ){

			for( core::Size regcount = 1; regcount <= flex_regions_.size(); ++regcount ){

				generate_ensemble_for_region( *poly_ala_pose, regcount );

					//flex_regions_[regcount]->sort_ensemble_by_designability( pose, scorefxn_, design_pack_task_template);

				flex_regions_[regcount]->hack_fillup_frag_designabilities();

				if( basic::options::option[basic::options::OptionKeys::enzdes::enz_debug]) break;
			} //loop over flexible regions
		}

		/// Quit if we're only trying to generate loop files for the flexible regions.
		if ( basic::options::option[ basic::options::OptionKeys::enzdes::dump_loop_samples ]() != "no"){

			std::string loops_pdb = pose.pdb_info()->name();
			utility::file::FileName fname( loops_pdb );
			loops_pdb = fname.base() + "_flex_loops.pdb";

			core::fragment::dump_frames_as_pdb( *poly_ala_pose, flex_regions_, loops_pdb, 2 );
		}

		if ( basic::options::option[ basic::options::OptionKeys::enzdes::dump_loop_samples ]() == "quit_afterwards" ) {
			utility_exit_with_message( "Exiting as requested by option enzdes::dump_loop_samples" );
		}

		scorefxn_->set_weight(core::scoring::fa_sol, cur_emap[ core::scoring::fa_sol ]);

		//ok, we have our ensembles sorted by designability, now assemble combinations
		//of ensemble members in order of best energies
		//note: we can't minimize bond angles yet, so set them to 0
		scorefxn_->set_weight( core::scoring::mm_bend, 0.0 );

		PackerTaskOP flex_pack_task = enzutil::recreate_task( pose, *design_pack_task_template );

		core::scoring::ScoreFunctionCOP flexpack_sfxn = scorefxn_;
		if( basic::options::option[basic::options::OptionKeys::packing::soft_rep_design] ) flexpack_sfxn = soft_scorefxn_;

		flexpack::FlexPackerOP flex_packer = new flexpack::FlexPacker( flex_pack_task, flex_regions_, flexpack_sfxn );

		time_t flex_start_time = time( NULL );
		flex_packer->apply( pose );
		time_t flex_end_time = time( NULL );

		tr << " flexpacker took " << long( flex_end_time - flex_start_time ) << " seconds. " << std::endl;

		if( option[OptionKeys::enzdes::cst_min] ) cst_minimize( pose, flex_pack_task );

		//should we do some fixbb design rounds afterward?
		if( option[OptionKeys::enzdes::design_min_cycles] > 1 ){

			core::Size fixbb_cycles = option[OptionKeys::enzdes::design_min_cycles] - 1;
			tr << "Doing " << fixbb_cycles << " rounds of fixbb design/min after flexpacking... " << std::endl;

			PackerTaskOP fix_pack_task = enzutil::recreate_task( pose, *design_pack_task_template );

			enzdes_pack( pose, fix_pack_task, scorefxn_, fixbb_cycles, basic::options::option[basic::options::OptionKeys::enzdes::cst_min], false, option[ OptionKeys::enzdes::favor_native_res].user() );


		}

		//do a repack without constraints
		if( ! basic::options::option[basic::options::OptionKeys::enzdes::no_unconstrained_repack]){
			PackerTaskOP repack_task = create_enzdes_pack_task( pose, false );
			enzdes_pack( pose, repack_task, scorefxn_, basic::options::option[basic::options::OptionKeys::enzdes::cst_min].user(), 1, true, false );
		}

	} //if cst_design

} //apply function


std::string
EnzdesFlexBBProtocol::get_name() const {
	return "EnzdesFlexBBProtocol";
}

void
EnzdesFlexBBProtocol::register_options()
{

	using namespace basic::options;

	option.add_relevant( OptionKeys::enzdes::dump_loop_samples );
	option.add_relevant( OptionKeys::loops::loop_file );
	option.add_relevant( OptionKeys::enzdes::enz_loops_file );
	option.add_relevant( OptionKeys::enzdes::cst_predock );
	option.add_relevant( OptionKeys::enzdes::kic_loop_sampling);
	option.add_relevant( OptionKeys::enzdes::mc_kt_low);
	option.add_relevant( OptionKeys::enzdes::mc_kt_high);
	option.add_relevant( OptionKeys::enzdes::single_loop_ensemble_size);
	option.add_relevant( OptionKeys::enzdes::loop_generator_trials );

	option.add_relevant( OptionKeys::enzdes::no_catres_min_in_loopgen );
	option.add_relevant( OptionKeys::enzdes::min_cacb_deviation );
	option.add_relevant( OptionKeys::enzdes::max_bb_deviation_from_startstruct );
	option.add_relevant( OptionKeys::enzdes::max_bb_deviation );
	option.add_relevant( OptionKeys::enzdes::checkpoint );

	protocols::backrub::BackrubMover::register_options();

}



bool
EnzdesFlexBBProtocol::is_flexible( core::Size seqpos ) const
{
	for( core::Size i = 1; i <= flex_regions_.size(); ++i ){
		if( ( seqpos >= flex_regions_[i]->start() ) && ( seqpos <= flex_regions_[i]->end() ) ) return true;
	}
	return false;
}

bool
EnzdesFlexBBProtocol::is_remodelable( core::Size seqpos ) const
{
	for( core::Size i = 1; i <= flex_regions_.size(); ++i ){
		if( !flex_regions_[i]->remodelable() ) continue;
		if( ( seqpos >= flex_regions_[i]->start() ) && ( seqpos <= flex_regions_[i]->end() ) ) return true;
	}

	return false;
}


/// As output, residue_position[i] is true for all neighbor residues including orginal input residues.
/// The function is used to find all neighboring residues of the loop residues in case they need to be
/// repacked or minimized in fullatom refinement.
void
EnzdesFlexBBProtocol::get_tenA_neighbor_residues(
  core::pose::Pose const & pose,
  utility::vector1<bool> & residue_positions
) const
{
  //make a local copy first because we will change content in residue_positions
	core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
	for ( Size i=1; i <= pose.total_residue(); ++i ) {
		if ( !is_remodelable(i) && !is_flexible(i) ) continue;
		core::graph::Node const * current_node( tenA_neighbor_graph.get_node(i)); // find neighbors for this node
		for ( core::graph::Node::EdgeListConstIter it = current_node->const_edge_list_begin();
			it != current_node->const_edge_list_end(); ++it ) {
			Size pos = (*it)->get_other_ind(i);
			if (pose.residue(pos).type().name() == "CYD") continue;
			residue_positions[ pos ] = true;
		//		tr << "residue pos TenAGraph " << pos << std::endl;
		}
	}
}


core::pack::task::PackerTaskOP
EnzdesFlexBBProtocol::modified_task(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask const & orig_task
){

	using namespace core::pack::task;

	tr << "Modifiyng task according to flexible regions. The following residues are additionally set to designing: ";

	PackerTaskOP mod_task = TaskFactory::create_packer_task( pose );
	mod_task->initialize_from_command_line();

	utility::vector1 < bool > remodel_loop_interface( pose.total_residue(), false );
	get_tenA_neighbor_residues(pose, remodel_loop_interface);

	for( core::Size i = 1; i <= pose.total_residue(); i++ ){

		//first, we need to copy the rotamer and rotamerset operations
		for( core::pack::rotamer_set::RotamerOperations::const_iterator rot_it = orig_task.residue_task(i).rotamer_operations().begin(); rot_it != orig_task.residue_task(i).rotamer_operations().end(); ++rot_it ){
			mod_task->nonconst_residue_task( i ).append_rotamer_operation( *rot_it );
		}
		for( core::pack::rotamer_set::RotSetOperationListIterator rotset_it = orig_task.residue_task(i).rotamer_set_operation_begin(); rotset_it != orig_task.residue_task(i).rotamer_set_operation_end(); ++rotset_it ){
			mod_task->nonconst_residue_task( i ).append_rotamerset_operation( *rotset_it );
		}

		if( is_catalytic_position( pose, i ) ) {

			if( !basic::options::option[basic::options::OptionKeys::enzdes::fix_catalytic_aa] ) mod_task->nonconst_residue_task(i).restrict_to_repacking();

			else mod_task->nonconst_residue_task(i).prevent_repacking();
		}

		/// APL TEMP HACK: Try to make the flexbb protocol as similar to the fixbb protocol as possible by not allowing design
		/// at extra positions.
		else if( orig_task.pack_residue( i ) && is_flexible( i ) && !is_remodelable( i ) && ( ! orig_task.design_residue( i )) ){
			mod_task->nonconst_residue_task(i).restrict_to_repacking();
		}

		else if( orig_task.pack_residue( i ) && !is_flexible( i ) && ( ! orig_task.design_residue( i )) ){
			mod_task->nonconst_residue_task(i).restrict_to_repacking();
		}

		else if (	remodel_loop_interface[i] && !orig_task.design_residue( i ) && !orig_task.pack_residue( i ) && !is_flexible( i )  && !is_remodelable( i ) ){
			mod_task->nonconst_residue_task(i).restrict_to_repacking();
			//tr << i << "+";
		}

		else if( !(orig_task.pack_residue( i ) || orig_task.design_residue( i ) ) && !is_flexible( i ) ){
			mod_task->nonconst_residue_task( i ).prevent_repacking() ;
		}

		//tr << std::endl;
		//decision done

		if( mod_task->design_residue(i) ) {

			//std::cerr << i << " is designing in mod task,  ";

			//ok, some restrictions apply
			if( pose.residue( i ).name3() == "CYS" && pose.residue( i ).has_variant_type( core::chemical::DISULFIDE ) ){
				mod_task->nonconst_residue_task( i ).restrict_to_repacking();
			}

			else{

				utility::vector1< bool > all_aas( core::chemical::num_canonical_aas, true );
				all_aas[ core::chemical::aa_cys ] = false;

				utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );

				for( ResidueLevelTask::ResidueTypeCOPListConstIter res_it = orig_task.residue_task( i ).allowed_residue_types_begin(); res_it != orig_task.residue_task( i ).allowed_residue_types_end(); ++res_it) {

					keep_aas[ (*res_it)->aa() ] = true;
				}


				if( orig_task.design_residue(i) ){
					mod_task->nonconst_residue_task(i).restrict_absent_canonical_aas( keep_aas );
				}

				else {
					tr << i << "+";
					mod_task->nonconst_residue_task(i).restrict_absent_canonical_aas( all_aas );
				}
			}
			//std::cerr << std::endl;
		}

	} //loop over all residues
	tr << std::endl;

	//don't forget to copy the upweighters
	if( orig_task.IGEdgeReweights() ) {
		for( utility::vector1< IGEdgeReweighterOP >::const_iterator it = orig_task.IGEdgeReweights()->reweighters_begin(); it != orig_task.IGEdgeReweights()->reweighters_end(); ++it){
			mod_task->set_IGEdgeReweights()->add_reweighter( *it );
		}
	}

	return mod_task;
}



void
EnzdesFlexBBProtocol::remap_resid(
	core::pose::Pose const & pose,
	core::id::SequenceMapping const & smap
){
	EnzdesBaseProtocol::remap_resid( pose, smap );
	for( core::Size regcount = 1; regcount <= flex_regions_.size(); ++regcount ){
		if( !flex_regions_[regcount]->remap_resid( pose, smap ) ) utility_exit_with_message("Failed to remap resid for flexible region");
	}
}

void
EnzdesFlexBBProtocol::add_flexible_region(
	core::Size start,
	core::Size end,
	core::pose::Pose const & pose,
	bool clear_existing
)
{
	if( clear_existing ) flex_regions_.clear();
	core::Size length( end - start +1 );
	flex_regions_.push_back( new EnzdesFlexibleRegion( flex_regions_.size() + 1, start, end, length, pose, this ) );
}


void
EnzdesFlexBBProtocol::determine_flexible_regions(
	core::pose::Pose const & pose,
	core::pack::task::PackerTaskCOP task
)
{
	//flex regions need to have a minimum length, since we don't want to be diversifying 2 or 3 residue stretches
	core::Size no_flex_regions(0);
	flex_regions_.clear();

	tr << "Determining regions to be treated as flexible... " << std::endl;

	//is there an enzdes loops file?
	if( toolbox::match_enzdes_util::get_enzdes_observer( pose ) && toolbox::match_enzdes_util::get_enzdes_observer( pose )->enzdes_loops_file() ) enz_loops_file_ = toolbox::match_enzdes_util::get_enzdes_observer( pose )->enzdes_loops_file();
	else if( !enz_loops_file_ && basic::options::option[ basic::options::OptionKeys::enzdes::enz_loops_file ].user() ){

		toolbox::match_enzdes_util::EnzdesLoopsFileOP loops_file = new toolbox::match_enzdes_util::EnzdesLoopsFile();

		if( !loops_file->read_loops_file( basic::options::option[ basic::options::OptionKeys::enzdes::enz_loops_file ] ) ){
			utility_exit_with_message("Reading enzdes loops file failed");
		}
		enz_loops_file_ = loops_file;
	}

	if( enz_loops_file_ ){

		//if the SpecialSegmentsObserver in the pose cache has been set,
		//loop start and end will be taken from it, in case another
		//loop has been previously remodeled with indels
		core::pose::datacache::SpecialSegmentsObserverCOP segob(NULL);
		if( pose.observer_cache().has( core::pose::datacache::CacheableObserverType::SPECIAL_SEGMENTS_OBSERVER )){
			core::pose::datacache::CacheableObserverCOP cache_ob = pose.observer_cache().get_const_ptr( core::pose::datacache::CacheableObserverType::SPECIAL_SEGMENTS_OBSERVER);
			segob = utility::pointer::static_pointer_cast< core::pose::datacache::SpecialSegmentsObserver const >( cache_ob );
			runtime_assert( segob->segments().size() == enz_loops_file_->num_loops() );
		}

		for( core::Size i(1); i <= enz_loops_file_->num_loops(); ++i){

			core::Size lstart(0), lstop(0);

			if( segob ){
				lstart = segob->segments()[i].first;
				lstop = segob->segments()[i].second - 1; //segment convention
			}
			else if ( enz_loops_file_->loop_info( i )->pose_numb() ) {

				lstart = enz_loops_file_->loop_info( i )->start();

				lstop = enz_loops_file_->loop_info( i )->stop();

			} else if ( enz_loops_file_->loop_info( i )->pdb_numb() ) {

				lstart =
					pose.pdb_info() -> pdb2pose(
						enz_loops_file_->loop_info( i )->start_pdb_chain(),
						enz_loops_file_->loop_info( i )->start_pdb() );

				lstop =
					pose.pdb_info() -> pdb2pose(
						enz_loops_file_->loop_info( i )->stop_pdb_chain(),
						enz_loops_file_->loop_info( i )->stop_pdb() );

			}

			core::Size length( lstop - lstart +1 );

			flex_regions_.push_back( new EnzdesFlexibleRegion( i, lstart, lstop, length, pose, this ) );
			tr << " " << lstart << "-" << lstop;

			core::Size min_length( enz_loops_file_->loop_info( i )->min_length() );
			core::Size max_length( enz_loops_file_->loop_info( i )->max_length() );

			//if( (min_length != length) || (max_length != length ) ){
				flex_regions_[i]->declare_remodelable( min_length, max_length );
				tr << " (remodelable with min_length " << min_length << " and max_length " << max_length << ")";
			//}
			tr << ", ";
		}
	}

	//is there a regular loops file?
	else if( basic::options::option[basic::options::OptionKeys::loops::loop_file].user() ){

        // we're hijacking the loop file reader
		loops::Loops loops_helper( true );

		tr << "reading information from loops file " << loops_helper.loop_file_name() << ": loops are " ;


		for( utility::vector1< loops::Loop >::const_iterator lit = loops_helper.v_begin(); lit != loops_helper.v_end(); ++lit){
			no_flex_regions++;
			core::Size lstart( lit->start() ), lstop( lit->stop() );
			flex_regions_.push_back( new EnzdesFlexibleRegion( no_flex_regions, lstart, lstop, lstop - lstart + 1, pose, this ) );
			tr << " " << lstart << "-" << lstop;

			if( lit->is_extended() ){

				core::Size min_length( lit->cut() ), max_length( (core::Size) lit->skip_rate() );

				flex_regions_[no_flex_regions]->declare_remodelable( min_length, max_length );

				tr << " (remodelable with min_length " << min_length << " and max_length " << max_length << ")";
			}
			tr << ", ";
		}
	} else {
		//if not, determine the flex regions automatically

		core::Size const min_flex_length = 6;
		utility::vector1< bool > flex_res( pose.total_residue(), false );

		for( core::Size i = 1; i<= pose.total_residue(); ++i){
			if( ( task->design_residue( i ) || is_catalytic_position( pose, i ) ) && pose.residue(i).is_polymer() ) flex_res[i] = true;
		}

		enzutil::make_continuous_true_regions_in_bool_vector( flex_res, min_flex_length );

		//make sure that the first residue isn't flexible and that no ligand was set to flexible
		flex_res[1] = false;
		for( core::Size i = 1; i<= pose.total_residue(); ++i){
			if( flex_res[i] && !pose.residue(i).is_polymer() ){

				core::Size lower( std::max( core::Size (1), i-1) );
				core::Size upper( std::min( i+1, pose.total_residue() ) );
				if( (flex_res[lower] && (lower != i ) ) && (flex_res[upper] && (upper != i ) ) ) utility_exit_with_message("Somehow a non polymer residue got into the middle of a flexible region.");

				flex_res[i] = false;
			}
		}

		for( core::Size i = 1; i<= pose.total_residue(); ++i){
			if( flex_res[i] ){
				core::Size j = i;


				while( flex_res[j] && (j <= pose.total_residue()) ) j++;

				no_flex_regions++;
				flex_regions_.push_back( new EnzdesFlexibleRegion( no_flex_regions, i, j - 1, (j - i), pose, this ) );
				//fragment_counters_.push_back( 1 );
				tr << "found " << i << "-" << j - 1 << ", ";
				i = j;
			}
		}
	} //automatic flex region determination

	for( core::Size i = 1; i <= flex_regions_.size(); ++i) fragment_counters_.push_back( 1 );
	tr << flex_regions_.size() << " flexible regions in total." << std::endl;

} //determine_flexible_regions function



void
EnzdesFlexBBProtocol::generate_ensemble_for_region(
	core::pose::Pose & pose,
	core::Size region
)
{
	time_t start_time = time(NULL);

	tr << "Starting to generate ensemble for region " << region << " ( aiming for " << loop_ensemble_size_ << " members, " << loopgen_trials_ << " brub trials for each ) ... " << std::endl;
	tr.flush();

	if ( ! basic::options::option[ basic::options::OptionKeys::enzdes::kic_loop_sampling ] ) {
		brub_mover_ = new protocols::backrub::BackrubMover();
		//brub_mover_->set_native_pose( & pose );
	}
	kinematic_mover_ = new protocols::loops::loop_closure::kinematic_closure::KinematicMover();

	(*reduced_scorefxn())( pose );

	minimize_cats_ = !( basic::options::option[basic::options::OptionKeys::enzdes::no_catres_min_in_loopgen].user()) && flex_regions_[region]->contains_catalytic_res();

	if( minimize_cats_ ) setup_catalytic_residue_minimization_for_region( pose, region );

	Size const region_size( flex_regions_[region]->positions().size() );

	Size const rbegin( flex_regions_[region]->positions()[ 1 ] );
	Size const rend( flex_regions_[region]->positions()[ flex_regions_[region]->positions().size()] );
	Size const region_middle( (region_size + 1) / 2 );
	Size rmid( flex_regions_[region]->positions()[ region_middle ] );

	/// NOTE: proline and pre-proline residues have a very sensitive Rama distribution; try to avoid
	/// them as pivot residues.
	if ( pose.residue_type( rmid ).aa() == core::chemical::aa_pro ) {
		// try to avoid prolines as the middle residue in the pivot.
		rmid += 1;
	} else if ( pose.residue_type( rmid + 1 ).aa() == core::chemical::aa_pro ) {
		// try also to avoid pre-pro residues.
		rmid -= 1;
	}

	core::kinematics::FoldTree ft_old = pose.fold_tree();
	core::kinematics::FoldTree ft_temp = pose.fold_tree();

	if ( rend <= pose.total_residue() - 3 && pose.chain( rbegin ) == pose.chain( rend + 2 ) ) {
		using namespace core::kinematics;
		//std::cout << "orig fold tree " << ft << std::endl;
		//ft_temp.new_jump( rbegin, rend+2, rend+1 );
		ft_temp.new_jump( rbegin, rend+3, rend+2 );
		//std::cout << "new fold tree? " << ft << std::endl;
		pose.fold_tree( ft_temp );
	}

	kinematic_mover_->set_pivots(rbegin, rmid, rend);

	protocols::loops::loop_closure::kinematic_closure::VicinitySamplingKinematicPerturberOP perturber = new protocols::loops::loop_closure::kinematic_closure::VicinitySamplingKinematicPerturber( &(*kinematic_mover_) );

	if ( basic::options::option[ basic::options::OptionKeys::enzdes::kic_loop_sampling ] ) {

		protocols::simple_filters::ScoreCutoffFilterOP bump_filter = new protocols::simple_filters::ScoreCutoffFilter();
		bump_filter->set_positions( flex_regions_[region]->positions() );
		bump_filter->set_score_type( core::scoring::fa_rep );
		bump_filter->set_cutoff( bump_filter->get_score( pose ) + region_size * 0.2 );

		kinematic_mover_->clear_filters();
		kinematic_mover_->set_rama_check( false );
		kinematic_mover_->set_hardsphere_bump_check( true );
		//kinematic_mover_->add_filter( bump_filter );
		kinematic_mover_->set_sfxn( reduced_scorefxn() );
		kinematic_mover_->set_do_sfxn_eval_every_iteration( false );


		core::Real sample_vicinity = 15.0; //( region_size > 10 ? 1.0 : 40.0 / region_size );
		//perturber->set_sample_vicinity( true );
		perturber->set_degree_vicinity( sample_vicinity );
		perturber->set_max_sample_iterations( 100 );

		kinematic_mover_->set_idealize_loop_first( false );
		kinematic_mover_->set_perturber( perturber );

		generate_alc_ensemble_for_region( pose, region );
		//if ( flex_regions_[ region ]->nr_frags() < loop_ensemble_size_ ) {
			/// try again choosing different pivot residues?
		//	generate_alc_ensemble_for_region( pose, region, true );
		//}
	} else {

		core::Real sample_vicinity = 2.0;
		//perturber->set_sample_vicinity( true );
		perturber->set_degree_vicinity( sample_vicinity );
		perturber->set_max_sample_iterations( 100 );

		kinematic_mover_->set_rama_check( false );
		kinematic_mover_->set_hardsphere_bump_check( false );
		kinematic_mover_->set_do_sfxn_eval_every_iteration( false );
		kinematic_mover_->set_idealize_loop_first( true );
		kinematic_mover_->set_perturber( perturber );

		generate_backrub_ensemble_for_region( pose, region );
	}

	minimize_cats_ = false;
	catmin_movemap_ = NULL;
	catmin_mover_ = NULL;
	brub_mover_ = NULL;
	pose.fold_tree( ft_old );

	time_t end_time = time(NULL);

	tr << " done generating ensemble for region " << region << " in " << long(end_time - start_time ) << " seconds. " << std::endl;
	tr.flush();

}

bool
EnzdesFlexBBProtocol::minimize_flexible_region(
	core::pose::Pose & pose,
		core::Size region,
		core::scoring::ScoreFunctionCOP scorefxn,
		std::set< core::Size > const & chi_to_move,
		bool const including_CA_angles,
		core::Real min_tolerance
)
{
	return flex_regions_[region]->minimize_region( pose, scorefxn, chi_to_move, including_CA_angles, min_tolerance  );
}

void
EnzdesFlexBBProtocol::generate_alc_ensemble_for_region(
	core::pose::Pose & pose,
	core::Size region
)
{
	using namespace core;
	using namespace core::fragment;
	using namespace core::chemical;

	flex_regions_[region]->set_target_proximity_to_starting_conformation(
		std::min( 0.2, basic::options::option[basic::options::OptionKeys::enzdes::max_bb_deviation_from_startstruct] / 2) );
	/// It's easy to generate tons of conformations; score and store them and then weed through them later.
	std::list< std::pair< Real, FragDataOP > > scored_confs;
	FragDataOP example_fragment = flex_regions_[region]->fragment_ptr( 1 )->clone();

	/// Faster kinematics and scoring if the right fold-tree topology is in place.
	core::pose::Pose local_pose( pose );
	(*reduced_scorefxn())( local_pose );

	//if rend is near the end of the pose, then it won't help to add a jump, so don't bother
	//if ( rend <= pose.total_residue() - 2 && pose.chain( rbegin ) == pose.chain( rend + 2 ) ) {
	//	using namespace core::kinematics;
	//	FoldTree ft = pose.fold_tree();
		//std::cout << "orig fold tree " << ft << std::endl;
	//	ft.new_jump( rbegin, rend+2, rend+1 );
		//std::cout << "new fold tree? " << ft << std::endl;
	//	local_pose.fold_tree( ft );
	//}

	//local_pose.dump_pdb( "alc_local_pose_" + utility::to_string( rbegin ) + "_" + utility::to_string( rend ) + ".pdb" );

	core::Size const rbegin( flex_regions_[region]->positions()[ 1 ] );
	core::Size const rend( flex_regions_[region]->positions()[ flex_regions_[region]->positions().size()] );
	core::Size const len = rend - rbegin + 1;

	utility::vector1< core::pose::PoseOP > loop_poses; loop_poses.reserve( loop_ensemble_size_ );
	core::pose::PoseOP native_loop_pose = new core::pose::Pose();
	utility::vector1< core::Real > rmsd_to_native; rmsd_to_native.reserve( loop_ensemble_size_ );
	flex_regions_[region]->fragment_as_pose( 1, *native_loop_pose, this->restype_set() );
	native_loop_pose->prepend_polymer_residue_before_seqpos( pose.residue( flex_regions_[region]->start() - 1) , 1, false );
	native_loop_pose->copy_segment( flex_regions_[region]->length() + 1, pose, 1, flex_regions_[region]->start() - 1 );
	loop_poses.push_back( native_loop_pose );
	core::pose::Pose loop_template_pose = *native_loop_pose;

	core::Real mc_decrement( (mc_kt_high_ - mc_kt_low_ ) / loopgen_trials_ );
	core::Real mc_temp( mc_kt_high_);

	protocols::moves::MonteCarlo mc( local_pose, *reduced_scorefxn(), mc_temp );


	static Size count_output = 1;
	core::Real kinmover_successrate(0.0);
	core::Size kinmover_successes(0), kintrials( loopgen_trials_ );

	//core::Size successive_failures(0);

	//core::Real const native_score = local_pose.energies().total_energies()[ core::scoring::total_score ];

	for( core::Size outerloop = 1; outerloop <= ( 5*loop_ensemble_size_) ; outerloop++ ){

		mc.reset( local_pose );
		mc_temp = mc_kt_high_;
		mc.set_temperature( mc_temp );

		kinmover_successes = 0;

		for ( core::Size kinits = 1; kinits <= kintrials; ++kinits ){

			//every 10th iteration reset the pose to the native to avoid drift
			//if( kinits % 10 == 0 ) example_fragment->apply( local_pose, flex_regions_[region]->start(), flex_regions_[region]->stop() );


			//if the kinematic mover has sampled itself to a region where it's
			//difficult to get closed solutions, put back one of the previous confs
			//if( successive_failures > 5 ){

			//	if( flex_regions_[region]->apply( numeric::random::random_range( 1, flex_regions_[region]->nr_frags() ), pose ) != flex_regions_[region]->length() ) utility_exit_with_message("unknown error when trying to apply a random fragment during ensemble generation.");
			//}

			//std::cerr << "outerloop " << outerloop << " Kinit " << kinits << ".... ";
			if( len > 3 ){
				core::Size p1( numeric::random::random_range( rbegin, rend - 3) );
				core::Size p2( numeric::random::random_range( p1+1, rend - 1) );
				core::Size p3( numeric::random::random_range( p2+1, rend) ) ;
				kinematic_mover_->set_pivots(p1, p2, p3);
			}

			kinematic_mover_->apply( local_pose );

			if( !kinematic_mover_->last_move_succeeded() ){
				//successive_failures++;
				continue;
			}
			//else successive_failures = 0;

			kinmover_successes++;

			(*reduced_scorefxn())( local_pose );

			//std::cerr << "Kinit " << kinits << " was successful.... " << std::endl;
			if( minimize_cats_ ) catmin_mover_->apply( local_pose );

			++count_output;
		//std::cout << "Found one sc: " << (*reduced_scorefxn())( local_pose ) << std::endl;
		///local_pose.dump_pdb( "test_sweep_" + utility::to_string( count_output ) + ".pdb" );
			Real score = local_pose.energies().total_energies()[ core::scoring::total_score ];

			FragDataOP newfrag = example_fragment->clone();
			newfrag->steal( local_pose, *flex_regions_[region] );
		//newfrag->apply( pose,  *flex_regions_[region] );
		//pose.dump_pdb( "test_sweep_after_apply_stolen_" + utility::to_string( count_output ) + ".pdb" );
			scored_confs.push_back( std::make_pair( score, newfrag ) );

			mc_temp = mc_temp - mc_decrement;
			mc.set_temperature( mc_temp );
			mc.boltzmann( local_pose );

		} //kinit iterations

		//on the first iteration, we determine how difficult it is for the kinematic mover
		//to generate confs for this particular loop
		//then we scale kintrials by the successrate
		if( kinmover_successes > 0  ){
			kinmover_successrate = ( (core::Real) kinmover_successes) / ( (core::Real) kintrials );
			kintrials = std::min( 5 * loopgen_trials_,  (core::Size)(loopgen_trials_ / kinmover_successrate) );

			//std::cerr << "setting kinmover succesrate to " << kinmover_successrate << " ( " << kinmover_successes << " successful moves ) and kintrials to " << kintrials << "    ";
		}

	//native_loop_pose->dump_pdb( "alc_native_loop_pose_" + utility::to_string( rbegin ) + "_" + utility::to_string( rend ) + ".pdb" );


		//std::sort( scored_confs.begin(), scored_confs.end(), utility::SortFirst< Real, FragDataOP >() );
		//scored_confs.sort( utility::SortFirst< Real, FragDataOP >() );

	/// check the fragments the first time.  we're only examining those that have a better score than the native
		/*
		for ( std::list< std::pair< Real, FragDataOP > >::iterator list_it = scored_confs.begin(),
						list_end = scored_confs.end(); list_it != list_end; ) {

			//std::list< std::pair< Real, FragDataOP > >::iterator iter_next = list_it;
			//iter_next++;

			if( list_it->first > native_score ){
				break;
			}

		//std::cout << "ScoredConfs sorted: " << ii << " " << scored_confs[ ii ].first << " accepted: " << flex_regions_[region]->nr_frags() << " nloop_poses: " << loop_poses.size()  << std::endl;
			list_it->second->apply( local_pose, flex_regions_[region]->start(), flex_regions_[region]->stop() );

			if(flex_regions_[region]->examine_new_loopconf(  local_pose, loop_template_pose, loop_poses, rmsd_to_native ) ){
				list_it = scored_confs.erase( list_it );
				if ( flex_regions_[region]->nr_frags() == loop_ensemble_size_ ) break;
			}
			else ++list_it;
		}
		*/

		//we remember both the lowest score pose as well as the last accepted one
		flex_regions_[region]->examine_new_loopconf(  mc.lowest_score_pose(), loop_template_pose, loop_poses, rmsd_to_native );
		if ( flex_regions_[region]->nr_frags() == loop_ensemble_size_ ) break;

		flex_regions_[region]->examine_new_loopconf(  mc.last_accepted_pose(), loop_template_pose, loop_poses, rmsd_to_native );
		if ( flex_regions_[region]->nr_frags() == loop_ensemble_size_ ) break;


		if ( flex_regions_[region]->nr_frags() == loop_ensemble_size_ ) break;

		//now put a random one of the previously generated regions into the pose, might help with generating diversity
		if( flex_regions_[region]->apply( numeric::random::random_range( 1, flex_regions_[region]->nr_frags() ), pose ) != flex_regions_[region]->length() ) utility_exit_with_message("unknown error when trying to apply a random fragment during ensemble generation.");

		//std::cerr << kinmover_successes << " kinsuccesses in outeriteration " << outerloop << ", num frags is " << flex_regions_[region]->nr_frags() << std::endl;

	} //outerloop

	scored_confs.sort( utility::SortFirst< Real, FragDataOP >() );

	mc.show_counters();

	Size count( 0 );
	Size const count_limit( basic::options::option[basic::options::OptionKeys::enzdes::max_bb_deviation ].user() ? 4: 2 );

	if( flex_regions_[region]->nr_frags() < loop_ensemble_size_ ){
		tr << "Not enough fragments after monte carlo (" << flex_regions_[region]->nr_frags() << " so far) , going through all stored configurations to find more." << std::endl;
	}

	while ( flex_regions_[region]->nr_frags() < loop_ensemble_size_ && count < count_limit ) {

		if (  count == 2 ) {
			/// Reduce the smoothness filter if we're not getting hits
			flex_regions_[region]->scale_target_proximity_to_other_conformations( 3 );
			flex_regions_[region]->set_target_proximity_to_starting_conformation(
				basic::options::option[basic::options::OptionKeys::enzdes::max_bb_deviation_from_startstruct]);
		}

		if ( count == 3 ) {
			/// Reduce the proximity-to-native filter if we're not getting hits
			flex_regions_[region]->scale_target_proximity_to_starting_conformation( 1.5 );
		}

		if ( count == 4 ) {
			/// Go further if we have to!
			flex_regions_[region]->scale_target_proximity_to_other_conformations( 2 );
		}

		for ( std::list< std::pair< Real, FragDataOP > >::iterator list_it = scored_confs.begin(),
						list_end = scored_confs.end(); list_it != list_end; ) {
			//std::cout << "ScoredConfs sorted: " << ii << " " << scored_confs[ ii ].first << " accepted: " << flex_regions_[region]->nr_frags() << " nloop_poses: " << loop_poses.size()  << std::endl;
			list_it->second->apply( local_pose, flex_regions_[region]->start(), flex_regions_[region]->stop() );

			if (flex_regions_[region]->examine_new_loopconf(  local_pose, loop_template_pose, loop_poses, rmsd_to_native ) ){
				list_it = scored_confs.erase( list_it );
				if ( flex_regions_[region]->nr_frags() == loop_ensemble_size_ ) break;
			}
			else ++list_it;
		}
		++count;
	}

	core::Real av_rmsd = numeric::statistics::mean( rmsd_to_native.begin(), rmsd_to_native.end(), 0.0 );
	core::Real std_dev_rmsd = numeric::statistics::std_dev_with_provided_mean( rmsd_to_native.begin(), rmsd_to_native.end(), av_rmsd );

	tr << " done generating ensemble for region " << region << ". " << flex_regions_[region]->nr_frags() - 1 << " new unique fragments were generated. Average rmsd/stdev to native is " << av_rmsd << " +- " << std_dev_rmsd << ". Kinematic Mover had a success rate of " << kinmover_successrate <<  std::endl;
}

void
EnzdesFlexBBProtocol::generate_backrub_ensemble_for_region(
	core::pose::Pose & pose,
	core::Size region
)
{
	//ObjexxFCL::FArray1D_bool flex_res( pose.total_residue(), false);
	//for( core::Size i = flex_regions_[region]->start() - 1; i <= flex_regions_[region]->end(); ++i) flex_res(i) = true;
	brub_mover_->set_native_pose( new core::pose::Pose( pose ) );
	brub_mover_->set_input_pose( brub_mover_->get_native_pose() );
	brub_mover_->clear_segments();
	brub_mover_->add_mainchain_segments( flex_regions_[region]->positions(), brub_pivot_atoms_, brub_min_atoms_, brub_max_atoms_ );
	brub_mover_->optimize_branch_angles( pose );

	utility::vector1< core::pose::PoseOP > loop_poses;
	core::pose::PoseOP native_loop_pose = new core::pose::Pose();

	utility::vector1< core::Real > rmsd_to_native;

	flex_regions_[region]->fragment_as_pose( 1, *native_loop_pose, this->restype_set() );

	native_loop_pose->prepend_polymer_residue_before_seqpos( pose.residue( flex_regions_[region]->start() - 1) , 1, false );


	native_loop_pose->copy_segment( flex_regions_[region]->length() + 1, pose, 1, flex_regions_[region]->start() - 1 );

	loop_poses.push_back( native_loop_pose );

	core::pose::Pose loop_template_pose = *native_loop_pose;

	//native_loop_pose->dump_pdb("natloop"+utility::to_string( region )+".pdb");

	core::Real mc_decrement( (mc_kt_high_ - mc_kt_low_ ) / loopgen_trials_ );
	core::Real mc_temp( mc_kt_high_);

	protocols::moves::MonteCarlo mc( pose, *reduced_scorefxn(), mc_temp );

	//core::Size kintrials(0), kinfails(0);

	for( core::Size i = 1; i <= ( 5 * loop_ensemble_size_) ; ++i ){

		mc.reset( pose );
		mc_temp = mc_kt_high_;
		mc.set_temperature( mc_temp );

		for( core::Size j = 1; j <= loopgen_trials_; ++j ){

			PROF_START( basic::BACKRUB_MOVER );
			brub_mover_->apply( pose );
			PROF_STOP( basic::BACKRUB_MOVER );

			//we have to idealize the bond angles
			//kinematic mover previously setup to do this
			//kinematic_mover_->apply( pose );
			//kintrials++;

			//if( !kinematic_mover_->last_move_succeeded() ){
				//std::cerr << "Kinematic mover idealize fail on iteration " << j << "  " ;
				//kinfails++;
				//continue;
			//}
			if( minimize_cats_ ){

				catmin_mover_->apply( pose );
				(*reduced_scorefxn() )( pose );

			}

			//we have to idealize the bond angles
			//kinematic mover previously setup to do this

			PROF_START( basic::MC_ACCEPT );
			mc_temp = mc_temp - mc_decrement;
			mc.set_temperature( mc_temp );
			mc.boltzmann( pose, brub_mover_->type() );
			PROF_STOP( basic::MC_ACCEPT );
		}

		//std::cerr << "rmsd lowest to native is " << core::scoring::rmsd_no_super_subset( *this->get_native_pose(), mc.lowest_score_pose(), flex_res, core::scoring::is_protein_CA ) << ", ";

		//std::cerr << "rmsd lastaccep to native is " << core::scoring::rmsd_no_super_subset( *this->get_native_pose(), mc.last_accepted_pose(), flex_res, core::scoring::is_protein_CA ) <<  std::endl;

		//if( i == 3){
		//	mc.lowest_score_pose().dump_pdb("reg"+utility::to_string( region )+"_ens5.pdb");
		//	mc.last_accepted_pose().dump_pdb("reg"+utility::to_string( region )+"_ens6.pdb");
		//}

		//we remember both the lowest score pose as well as the last accepted one
		flex_regions_[region]->examine_new_loopconf(  mc.lowest_score_pose(), loop_template_pose, loop_poses, rmsd_to_native );
		if ( flex_regions_[region]->nr_frags() == loop_ensemble_size_ ) break;

		flex_regions_[region]->examine_new_loopconf(  mc.last_accepted_pose(), loop_template_pose, loop_poses, rmsd_to_native );
		if ( flex_regions_[region]->nr_frags() == loop_ensemble_size_ ) break;

		//if( ! flex_regions_[ region ]->steal( mc.lowest_score_pose() ) ) utility_exit_with_message("Could not steal a just generated bbfragment from the pose");
		//if( ! flex_regions_[ region ]->steal( mc.last_accepted_pose() ) ) utility_exit_with_message("Could not steal a just generated bbfragment from the pose");


		//now put a random one of the previously generated regions into the pose, might help with generating diversity
		if( flex_regions_[region]->apply( numeric::random::random_range( 1, flex_regions_[region]->nr_frags() ), pose ) != flex_regions_[region]->length() ) utility_exit_with_message("unknown error when trying to apply a random fragment during ensemble generation.");

	} //outerloop

	//for ( Size ii = 2; ii <= flex_regions_[ region ]->nr_frags(); ++ii ) {
	//	core::pose::Pose copy_pose( pose );
	//	flex_regions_[region]->apply( ii, copy_pose );
	//	copy_pose.dump_pdb("fullpose_loopreg_" + utility::to_string( region ) + "_" + utility::to_string( ii ) + ".pdb" );
	//}
	//basic::prof_show();

	mc.show_counters();

	//put back the native conformation, just in case
	if( flex_regions_[region]->apply( 1, pose ) != flex_regions_[region]->length() ) utility_exit_with_message("unknown error when trying to reapply native fragment after generating ensemble.");


	core::Real av_rmsd = numeric::statistics::mean( rmsd_to_native.begin(), rmsd_to_native.end(), 0.0 );
	core::Real std_dev_rmsd = numeric::statistics::std_dev_with_provided_mean( rmsd_to_native.begin(), rmsd_to_native.end(), av_rmsd );

	tr  << flex_regions_[region]->nr_frags() - 1 << " new unique fragments were generated. Average rmsd/stdev to native is " << av_rmsd << " +- " << std_dev_rmsd << "." <<  std::endl;

	//core::Real kinfailrate = ( (core::Real) kinfails ) / ((core::Real) kintrials);
	//tr << "Kinematic mover had a fail rate of " << kinfailrate << std::endl;

	//some debug shit below
	/*
	flex_regions_[region]->apply( 6, pose ); //restore native pose
	pose.dump_pdb("regident5_"+utility::to_string( region )+".pdb" );

	flex_regions_[region]->apply( 1, pose ); //restore native pose
	pose.dump_pdb("reg"+utility::to_string( region )+"nat.pdb" );

	*/
} //generate_ensemble_for_region


/// @details figure out which combination of loop conformations is the next most promising one
bool
EnzdesFlexBBProtocol::assemble_next_best_loop_combination(
	core::pose::Pose & pose
)
{
	runtime_assert( flex_regions_.size() == fragment_counters_.size() );

	//2. figure out for which of the regions we still have frags to try
	utility::vector1< bool > valid_regions( flex_regions_.size(), false );

	for( core::Size i = 1; i <= flex_regions_.size(); ++i){
		if( fragment_counters_[i] < flex_regions_[i]->no_ranked_frags() ) valid_regions[i] = true;
	}

	core::Real lowest_deltaE(10000000 );
	core::Size best_region(0);

	//3. then go through the fragments and put in one where the gain in energy is best
	for( core::Size i = 1; i <= flex_regions_.size(); ++i){

		if( valid_regions[i] && ( flex_regions_[i]->deltaE_best( fragment_counters_[i] + 1 ) < lowest_deltaE ) ){
			lowest_deltaE = flex_regions_[i]->deltaE_best( fragment_counters_[i] + 1 );
			best_region = i;
		}
	}

	if( best_region == 0 ) utility_exit_with_message("Trying to assemble a new combination of loops even though all combinations have alrady been assembled.");

	tr << "Assembling next best loop conformation: fragment " << fragment_counters_[best_region] << " of region " << best_region << " with deltaE_best " << lowest_deltaE << "put into pose." << std::endl;
	fragment_counters_[best_region]++;

	for( core::Size i = 1; i <= flex_regions_.size(); ++i){
		flex_regions_[i]->apply_ranked_fragment( pose, fragment_counters_[i] );
	}

	return true;

} //assemble_next_best_loop_combination


/// @details returns false if the last combination is reached, true otherwise
bool
EnzdesFlexBBProtocol::hack_assemble_next_loop_combination(
	core::pose::Pose & pose
	)
{

		assert( flex_regions_.size() == fragment_counters_.size() );

		//std::cerr << "MEEP calling hack_assemble_next_loop_combination ";

		core::Size first_region_at_end(0);
		core::Size last_region_at_end(0);

		bool all_regions_at_end(true);

		for( core::Size j = flex_regions_.size(); j >= 1; --j){

			if( fragment_counters_[j] == flex_regions_[j]->no_ranked_frags() ){
				first_region_at_end = j;
				if( last_region_at_end == 0 ) last_region_at_end = j;
			}
			else all_regions_at_end = false;

		}

		if( first_region_at_end == 0 ){
			fragment_counters_[ flex_regions_.size() ]++;
			flex_regions_[ flex_regions_.size() ]->apply_ranked_fragment( pose, fragment_counters_[ flex_regions_.size() ] );
			//tr << "put in frag " << fragment_counters_[ flex_regions_.size() ] << " of region " << flex_regions_.size() << ", returning at first point " << std::endl;
			return true;
		}

		if( all_regions_at_end ) return false; //means we've exhausted every combination

		if( flex_regions_[ flex_regions_.size() ]->no_ranked_frags() != fragment_counters_[ flex_regions_.size() ] ){
			fragment_counters_[ flex_regions_.size() ]++;
			flex_regions_[ flex_regions_.size() ]->apply_ranked_fragment( pose, fragment_counters_[ flex_regions_.size() ] );
			//tr << "put in frag " << fragment_counters_[ flex_regions_.size() ] << " of region " << flex_regions_.size() << ", returning at second point " << std::endl;
			return true;
		}

		else{
			for( core::Size i = flex_regions_.size(); i >= first_region_at_end; --i ){

				if( flex_regions_[ i ]->no_ranked_frags() != fragment_counters_[ i ] ){

					fragment_counters_[ i ]++;
					flex_regions_[ i ]->apply_ranked_fragment( pose, fragment_counters_[ i ] );
					//tr << " put in frag " << fragment_counters_[ i ] << " of reg " << i << ", ";
					for( core::Size jj = i + 1; jj <= flex_regions_.size(); ++jj){
						fragment_counters_[ jj ] = 1;
						flex_regions_[ jj ]->apply_ranked_fragment( pose, fragment_counters_[ jj ] );
						//tr << " put in frag " << fragment_counters_[ jj ] << " reg " << jj << ", ";
					}
					//tr << "returning" << std::endl;
					return true;
				}
			}
		}

		//fragment_counters_[ last_region_at_end ] = 1;
		//fragment_counters_[ last_region_at_end - 1 ]
		return false;
}


bool
EnzdesFlexBBProtocol::recover_loops_from_file( core::pose::Pose const & pose)
{

	using namespace basic::options;
	//annoying: let's mute a couple of tracer channels that cause complaints when input files are missing
	//OXT and termini H
	//utility::vector1< std::string > muted_channels = option[OptionKeys::out::mute].value();
	//muted_channels.push_back("core.io.pdb.file_data");
	//muted_channels.push_back("core.conformation.Conformation");
	//option[OptionKeys::out::mute].value() = muted_channels;

	//option[OptionKeys::out::mute].value().push_back( std::string("core.io.pdb.file_data"));
	//option[OptionKeys::out::mute].value().push_back( std::string("core.conformation.Conformation"));

	std::string loops_pdb = option[OptionKeys::enzdes::checkpoint];

	if( loops_pdb == "" ){

		loops_pdb = pose.pdb_info()->name();
		utility::file::FileName fname( loops_pdb );

		loops_pdb = fname.base() + "_flex_loops.pdb";
	}

	if( option[OptionKeys::enzdes::checkpoint].user() ){

		if( utility::file::file_exists( loops_pdb  ) ){

			tr << "Recovering loop conformations from file " << loops_pdb << "... " << std::endl;
			if( core::fragment::fill_template_frames_from_pdb( pose, flex_regions_, loops_pdb ) ){
				tr << " done recovering loop conformations." << std::endl;
				return true;
			}

			else{
				utility_exit_with_message("Unknown error when trying to recover loop conformations from file "+loops_pdb+".");

			}
		}
	}

	return false;

} //recover_loops_from_file


void
EnzdesFlexBBProtocol::setup_catalytic_residue_minimization_for_region(
	core::pose::Pose const & pose,
	core::Size region )
{

	catmin_sfxn_ = new core::scoring::ScoreFunction();
	catmin_sfxn_->reset();
	catmin_sfxn_->set_weight( core::scoring::fa_rep, reduced_scorefxn()->get_weight( core::scoring::fa_rep ) );
	catmin_sfxn_->set_weight( core::scoring::fa_dun, reduced_scorefxn()->get_weight( core::scoring::fa_dun ) );
	catmin_sfxn_->set_weight( core::scoring::coordinate_constraint, reduced_scorefxn()->get_weight( core::scoring::coordinate_constraint ) );
	catmin_sfxn_->set_weight( core::scoring::atom_pair_constraint, reduced_scorefxn()->get_weight( core::scoring::atom_pair_constraint ) );
	catmin_sfxn_->set_weight( core::scoring::angle_constraint, reduced_scorefxn()->get_weight( core::scoring::angle_constraint ) );
	catmin_sfxn_->set_weight( core::scoring::dihedral_constraint, reduced_scorefxn()->get_weight( core::scoring::dihedral_constraint ) );

	//bool minimize_cats = !( basic::options::option[basic::options::OptionKeys::enzdes::no_catres_min_in_loopgen].user()) && flex_regions_[region]->contains_catalytic_res();

	//if( minimize_cats )
	catmin_movemap_ = new core::kinematics::MoveMap();
	catmin_movemap_->clear();

	tr << "Allowing minimization of the following catalytic residues during ensemble generation for region " << region << ": ";
	for( core::Size rescount = flex_regions_[region]->start(); rescount <= flex_regions_[region]->end(); ++rescount ){
		if( is_catalytic_position( pose, rescount ) ){
			catmin_movemap_->set_chi( rescount, true);
			tr << rescount << ", ";
		}
	}

	catmin_mover_ = new protocols::simple_moves::MinMover( catmin_movemap_, catmin_sfxn_, "linmin", 0.02, true /*use_nblist*/ );
	tr << std::endl;

} //setup_catalytic_residue_minimization_for_region

EnzdesFlexibleRegion::EnzdesFlexibleRegion(
	core::Size index_in,
	core::Size start,
	core::Size end,
	core::Size nr_res,
	core::pose::Pose const & pose,
	EnzdesFlexBBProtocolCAP enz_prot
) :
	core::fragment::Frame( start, end, nr_res ),
	index_(index_in),
	enzdes_protocol_( enz_prot ),
	design_targets_( enz_prot->design_targets( pose ) ),
	target_proximity_to_native_conformation_(
		basic::options::option[basic::options::OptionKeys::enzdes::max_bb_deviation_from_startstruct] ),
	target_proximity_to_other_conformations_(
		basic::options::option[basic::options::OptionKeys::enzdes::max_bb_deviation].user() ?
		basic::options::option[basic::options::OptionKeys::enzdes::max_bb_deviation] :
			0.0 ),
	remodelable_(false),
	remodel_min_length_(nr_res),
	remodel_max_length_(nr_res)
{

	native_conf_ = assemble_enzdes_fragdata( pose );

	core::Size addfrag_returnval = add_fragment( native_conf_ );
	if( addfrag_returnval != 1 ) utility_exit_with_message("Could not add the native conformation to the EnzdesFlexibleRegion, returnval is "+utility::to_string( addfrag_returnval ) + ".");

	positions_.reserve( end - start + 1 );
	for ( core::Size i = start; i <= end; ++i ) positions_.push_back( i );
	frag_designabilities_.clear();

} //FlexibleRegion constructor

EnzdesFlexibleRegion::~EnzdesFlexibleRegion(){}


bool
EnzdesFlexibleRegion::contains_catalytic_res() const
{

	for( std::set< core::Size >::const_iterator cat_it = design_targets_.begin();
			 cat_it != design_targets_.end(); ++cat_it ){
		if( this->contains_seqpos( *cat_it ) ) return true;
	}
	return false;
}

toolbox::match_enzdes_util::EnzdesLoopInfoCOP
EnzdesFlexibleRegion::enz_loop_info() const
{
	toolbox::match_enzdes_util::EnzdesLoopsFileCOP loop_file = enzdes_protocol_->enz_loops_file();

	if( !loop_file ){
		utility_exit_with_message("no enzdes loops file was read, but the info therein requested." );
	}

	assert( index_ <= loop_file->num_loops() );

	return loop_file->loop_info( index_ );

}

void
EnzdesFlexibleRegion::declare_remodelable(
	core::Size min_length,
	core::Size max_length
)
{
	remodelable_ = true;
	remodel_min_length_ = min_length;
	remodel_max_length_ = max_length;
}

core::Real
EnzdesFlexibleRegion::deltaE_best(
	core::Size const frag_rank
) const
{
	if( frag_rank > frag_designabilities_.size() ) utility_exit_with_message("Trying to add a fragment that is not ranked.");

	return frag_designabilities_[frag_rank].second - frag_designabilities_[1].second;
} //deltaE_best


core::fragment::FragDataOP
EnzdesFlexibleRegion::assemble_enzdes_fragdata(
	core::pose::Pose const & pose
)
{

	using namespace core::fragment;

	FragDataOP new_fragdata = new FragData();

// tex 8/16/08
// lines containing core::id::AtomID[3] aren't portable. /// WHY NOT?!!?!
#ifndef WIN32
	for( core::Size i = this->start(); i<= this->end(); ++i){

		SingleResidueFragDataOP srfd;
		//utility::vector1< std::pair< core::Size[3], core::Real > > angles;
		utility::vector1<std::pair<std::vector<core::Size>, core::Real> > angles;  // REQUIRED FOR WINDOWS

		//for now we only keep the Ca angle
		//core::id::AtomID atoms[3] = { core::id::AtomID (pose.residue_type( i ).atom_index( " N  "), i),
		//															core::id::AtomID (pose.residue_type( i ).atom_index( " CA "), i),
		//															core::id::AtomID (pose.residue_type( i ).atom_index( " C  "), i) };

		//fuckin' A, icc compiler doesn't like the following line that gcc is perfectly fine with:
		//angles.push_back( std::pair< core::id::AtomID[3], core::Real>( atoms, 0 ) );
		//so we have to create (and copy) everything explicityly :(
		//std::pair< core::Size[3], core::Real> newpair;
		std::pair<std::vector<core::Size>, core::Real> newpair;  // REQUIRED FOR WINDOWS
		//newpair.first[0] = pose.residue_type( i ).atom_index( " N  ");
		//newpair.first[1] = pose.residue_type( i ).atom_index( " CA ");
		//newpair.first[2] = pose.residue_type( i ).atom_index( " C  ");
		newpair.first.push_back(pose.residue_type( i ).atom_index( " N  "));
		newpair.first.push_back(pose.residue_type( i ).atom_index( " CA "));
		newpair.first.push_back(pose.residue_type( i ).atom_index( " C  "));
		newpair.second = 0;

		angles.push_back( newpair );

		if( enzdes_protocol_->is_catalytic_position( pose, i ) ){
			srfd = new BBTorsionAndAnglesSRFD( angles ); //temporary, we will do something different for the catalytic res
		}
		else { srfd = new BBTorsionAndAnglesSRFD( angles ); }

		new_fragdata->add_residue( srfd );
	}

#endif
	if( ! new_fragdata->steal( pose, this->start(), this->end() ) ) utility_exit_with_message("Some error occured when trying to steal enzdes_fragdata from the pose. fuckin' A ... ");

	return new_fragdata;

} //assemble_enzdes_fragdata

void
EnzdesFlexibleRegion::hack_fillup_frag_designabilities()
{

	for( core::Size i = 1; i<= this->nr_frags(); ++i )
	{

		core::Real fake_designability( -1 * i * this->index_ );
		frag_designabilities_.push_back( SizeRealPair( i, fake_designability ) );
	}

} //hack_fillup_frag_designabilities();

/// @details this function assumes that the pose is a poly ala pose of the residues to redesign
/// @details and that the backbone is in the native conformation
void
EnzdesFlexibleRegion::sort_ensemble_by_designability(
	core::pose::Pose const & ref_pose,
	core::scoring::ScoreFunctionOP scorefxn,
	core::pack::task::PackerTaskCOP task
)
{
	using namespace core::pack;
	using namespace core::pack::task;

	//make a modifiable copy
	core::pose::Pose pose = ref_pose;
	(*scorefxn)( pose );
	//pose.dump_pdb("sortstart_reg"+utility::to_string(index_)+".pdb");

	utility::vector1< utility::vector1< core::Real > > lig_part_sums;

	//ok, for each of the ensemble members, we design and record the energy with respect to the ligand

	//first, let's make the right task (and what residues to mutate to alanine in this particular context)
	utility::vector1< core::Size > other_design_res;
	std::set< core::Size > ten_A_neighbors = get_10A_neighbors( pose );

	PackerTaskOP looptask_template = task->clone();

	for( core::Size i = 1; i <= pose.total_residue(); ++i ){

		if( looptask_template->design_residue( i )
			&& !this->contains_seqpos( i )
			&& !enzdes_protocol_->is_catalytic_position(pose, i) )
			{
				looptask_template->nonconst_residue_task(i).prevent_repacking();
				other_design_res.push_back( i );
			}

		else if( looptask_template->pack_residue(i)
			&& (ten_A_neighbors.find( i ) == ten_A_neighbors.end() )
			&& !this->contains_seqpos( i ) )
			{
				looptask_template->nonconst_residue_task(i).prevent_repacking();
			}

		if( this->contains_seqpos( i ) && !looptask_template->design_residue(i) && !enzdes_protocol_->is_catalytic_position(pose, i) ){
			utility_exit_with_message("Task is corrupted when trying to screen ensemble for region "+utility::to_string(index_));
		}
	}

	protocols::toolbox::pose_manipulation::construct_poly_ala_pose( pose, other_design_res, true, true, true );
	//pose.dump_pdb("sortala_reg"+utility::to_string(index_)+".pdb");

	IGEdgeReweighterOP ig_up = new protocols::toolbox::ResidueGroupIGEdgeUpweighter( 0.5, this->positions(), other_design_res);
	looptask_template->set_IGEdgeReweights()->add_reweighter( ig_up );

	tr << "Beginning designability screen for " << this->nr_frags() << " ensemble members of flexible region " << index_ << "... " << std::endl;
	time_t start_time = time(NULL);

	//get the native designability score
	PackerTaskOP looptask = enzutil::recreate_task( pose, *looptask_template );

	pose.update_residue_neighbors();
	core::pack::pack_scorefxn_pose_handshake( pose, *scorefxn );
	scorefxn->setup_for_packing( pose, looptask->repacking_residues(), looptask->designing_residues() );

	utility::vector1< core::Real > native_lig_part_sum_compare;
	native_lig_part_sum_compare.push_back( calculate_rotamer_set_design_targets_partition_sum( pose, scorefxn, looptask ) );

	pack_rotamers_loop( pose, *scorefxn, looptask, 1 );
	(*scorefxn)( pose );

	core::Real native_backgroundE(0);
	core::Real native_designability = extract_lig_designability_score( pose, task, native_backgroundE );

	native_lig_part_sum_compare.push_back( native_designability );

	lig_part_sums.push_back( native_lig_part_sum_compare );

	std::list< SizeRealPair > frag_designability_list;
	//now get the designability for every of the fragment ensemble members
	//start at 2 because one is the native
	for( core::Size i = 2; i<= this->nr_frags(); ++i){

		if( this->apply( i, pose ) != this->length() ) utility_exit_with_message("unknown error when trying to apply a fragment for generating designability score.");

		pose.update_residue_neighbors();
		core::pack::pack_scorefxn_pose_handshake( pose, *scorefxn );
		scorefxn->setup_for_packing( pose, looptask->repacking_residues(), looptask->designing_residues() );

		PackerTaskOP looptask = enzutil::recreate_task( pose, *looptask_template );

		utility::vector1< core::Real > lig_part_sum_compare;
		lig_part_sum_compare.push_back( calculate_rotamer_set_design_targets_partition_sum( pose, scorefxn, looptask ) );

		pack_rotamers_loop( pose, *scorefxn, looptask, 1 );
		(*scorefxn)( pose );

		core::Real backgroundE(0);
		core::Real designability = extract_lig_designability_score( pose, task, backgroundE );

		lig_part_sum_compare.push_back( designability );

		if( (designability < native_designability) && (backgroundE < native_backgroundE) ) {
			frag_designability_list.push_back( SizeRealPair( i, designability ) );
			//pose.dump_pdb("loopdes_reg"+utility::to_string( index_ )+"_"+utility::to_string( frag_designability_list.size() )+".pdb");
		}

		lig_part_sums.push_back( lig_part_sum_compare );
	}
	//make sure the best fragment comes first
	frag_designability_list.sort( compare_SizeRealPairs );

	for( std::list< SizeRealPair >::const_iterator it = frag_designability_list.begin(); it != frag_designability_list.end(); ++it) frag_designabilities_.push_back( *it );

	time_t end_time = time(NULL);
	tr << "done with designability screen for region " << index_ << " in " << (long)(end_time - start_time ) << " seconds, " << frag_designabilities_.size() << " of " << this->nr_frags() << " ensemble conformations are potentially better than the native conformation." << std::endl;

	if( frag_designabilities_.size() > 0 ){
		tr << "Native designability score is " << native_designability << ", best designability score is " << frag_designabilities_[1].second << ", worst designability score is " << frag_designabilities_[ frag_designabilities_.size() ].second << "." << std::endl;
	//and then append the native as the last fragment
	}
	frag_designabilities_.push_back( SizeRealPair( 1, native_designability) );

	//pose.dump_pdb("sortend_reg"+utility::to_string(index_)+".pdb");
	//and lastly put back the native conformation
	if( this->apply( 1, pose ) != this->length() ) utility_exit_with_message("unknown error when trying to reapply the native conformation after getting designability scores for fragments.");


	//last thing, for now, write out the measured lig designabilitis and the number out of the partition sums
	for( utility::vector1< utility::vector1< core::Real > >::const_iterator lig_part_sum_it = lig_part_sums.begin();
			 lig_part_sum_it != lig_part_sums.end(); ++lig_part_sum_it ){

		tr << "LIGPARTCOMPARE ";
		for( utility::vector1< core::Real >::const_iterator comp_it = (*lig_part_sum_it).begin(); comp_it != (*lig_part_sum_it).end(); ++comp_it) {
			tr << *comp_it << " ";
		}
		tr << std::endl;
	}
} //sort_enseble_by_designability



core::Real
EnzdesFlexibleRegion::calculate_rotamer_set_design_targets_partition_sum(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunctionCOP scorefxn,
	core::pack::task::PackerTaskCOP task
) const
{

	core::Real boltzmann_fac( 1 / 0.6 );

	core::Real dtps(0);

	//first, build a rotamer set
	core::pack::rotamer_set::RotamerSetsOP rotsets( new core::pack::rotamer_set::RotamerSets() );
	core::graph::GraphOP packer_neighbor_graph = core::pack::create_packer_graph( pose, *scorefxn, task );

	rotsets->set_task( task );

	rotsets->build_rotamers( pose, *scorefxn, packer_neighbor_graph );

	rotsets->prepare_sets_for_packing( pose, *scorefxn );

	//ig = InteractionGraphFactory::create_interaction_graph( *task, *rotsets, pose, *scorefxn );
	core::pack::interaction_graph::PDInteractionGraphOP ig = new core::pack::interaction_graph::PDInteractionGraph( task->num_to_be_packed() );

	ig->initialize( *rotsets );

	rotsets->compute_one_body_energies( pose, *scorefxn, packer_neighbor_graph, ig );

	//now delete the unnecessary edges from the graph
	utility::vector1< core::Size > residue_groups( pose.total_residue(), 0 );
	for( std::set< core::Size >::const_iterator cat_it = design_targets_.begin(); cat_it != design_targets_.end(); ++cat_it){
		if( (*cat_it >= *(positions_.begin() ) ) && (*cat_it <= *(positions_.rbegin() ) ) ){
			residue_groups[ *cat_it ] = 2;
			//tr << "CATKEEPGRAPH: resi " << *cat_it << " is part of loop between " << *(positions_.begin() ) << " and " << *(positions_.rbegin() ) << std::endl;
		}
		else residue_groups[ *cat_it ] = 1;
	}
	core::graph::delete_all_intragroup_edges( *packer_neighbor_graph, residue_groups );

	//and then precompute the energies of (hopefully) only the positions/ligand interactions
	rotsets->precompute_two_body_energies(  pose, *scorefxn, packer_neighbor_graph, ig );

	for( utility::vector1< core::Size >::const_iterator pos_it = positions_.begin(); pos_it != positions_.end(); ++pos_it ){

		core::Size moltenid = rotsets->resid_2_moltenres( *pos_it);

		core::Real dtps_this_position(0.0);


		for( ig->reset_edge_list_iterator_for_node( moltenid ); !ig->edge_list_iterator_at_end(); ig->increment_edge_list_iterator() ){

			//core::pack::interaction_graph::PDEdge const & cur_edge = (core::pack::interaction_graph::PDEdge) ig->get_edge();
			core::pack::interaction_graph::PDEdge const & cur_edge =  static_cast< core::pack::interaction_graph::PDEdge const & > ( ig->get_edge() );

			core::Size targ_moltenid = cur_edge.get_other_ind( moltenid );

			core::Size lower_res = std::min( moltenid, targ_moltenid );
			core::Size upper_res = std::max( moltenid, targ_moltenid );

			for( int ii = 1; ii <= ig->get_num_states_for_node( lower_res ); ++ii){

				core::Real pos_1body_E = ig->get_one_body_energy_for_node_state( lower_res, ii );
				if( pos_1body_E < 0.0 ) pos_1body_E = 0.0;

				for( int jj = 1; jj <= ig->get_num_states_for_node( upper_res ); ++jj){

					core::Real E_this_state = pos_1body_E + ig->get_one_body_energy_for_node_state( upper_res, jj ) + cur_edge.get_two_body_energy( ii, jj );

					dtps_this_position += exp(- boltzmann_fac * E_this_state );

					//tr << "BOLTZ pos " << rotsets->moltenres_2_resid( lower_res ) << " state " << ii << " to pos " << rotsets->moltenres_2_resid( upper_res ) << " state " << jj << ": interaction E " << E_this_state << "; pure int E " <<  cur_edge.get_two_body_energy( ii, jj ) << "; boltzmann fac " << exp(- boltzmann_fac * E_this_state ) << std::endl;

				} // loop over num states for design target

			} // loop over num states of loop residue

		} //loop over design target edges for this position
		//tr << "DTPS pos " << *pos_it << " is " << dtps_this_position << std::endl;
		dtps += dtps_this_position;

	} //loop over positions of this loop

	return ( dtps / design_targets_.size() );

} //calculate_rotamer_set_design_target_partition_sum(


/// @details function under heavy development, will prolly change a lot in the coming weeks/months
/// @details main idea: look at how the conformation of this region in the input pose interacts with
/// @details the ligand as well as the protein background, then combine the two numbers in some way
core::Real
EnzdesFlexibleRegion::extract_lig_designability_score(
	core::pose::Pose const & pose,
	core::pack::task::PackerTaskCOP task,
	core::Real & backgroundE
)
{
	using namespace core::scoring;

	//first, figure out which residues we're designing against( ligand and catalytic positions,
	//and sum up the interaction with them
	core::Real per_res_design_target_interactionE(0);
	core::Real av_background_interactionE(0);

	EnergyMap const cur_weights = pose.energies().weights();

	for( std::set< Size >::const_iterator targ_it = design_targets_.begin(); targ_it != design_targets_.end(); ++targ_it ){

		for( core::graph::EdgeListConstIterator egraph_it = pose.energies().energy_graph().get_node( *targ_it )->const_edge_list_begin();
       egraph_it != pose.energies().energy_graph().get_node( *targ_it )->const_edge_list_end(); ++egraph_it){

			core::Size other_res = (*egraph_it)->get_other_ind( *targ_it );
			if( !this->contains_seqpos( other_res ) ) continue;

			EnergyEdge const * Eedge = static_cast< EnergyEdge const * > (*egraph_it);

			per_res_design_target_interactionE += Eedge->dot( cur_weights );

		}//loop over interactig partners of particular design target
	} //loop over design targets
	per_res_design_target_interactionE /= design_targets_.size();


	//then, get the interaction energy of this stretch with itself and neighboring packing residues
	std::set< Size > interacting_neighbors;

	for( core::Size i = this->start(); i <= this->end(); ++i ){

		for( core::graph::EdgeListConstIterator egraph_it = pose.energies().energy_graph().get_node( i )->const_edge_list_begin();
				 egraph_it != pose.energies().energy_graph().get_node( i )->const_edge_list_end(); ++egraph_it){

			core::Size other_res = (*egraph_it)->get_other_ind( i );
			if( task->design_residue( other_res )
				|| ( this->contains_seqpos( other_res ) && (other_res < i ) ) //no overcounting please
				|| ( design_targets_.find( other_res ) != design_targets_.end() ) ) continue;

			if( !this->contains_seqpos( other_res )
				&& ( interacting_neighbors.find( other_res) == interacting_neighbors.end() ) ) interacting_neighbors.insert( other_res );

			EnergyEdge const * Eedge = static_cast< EnergyEdge const * > (*egraph_it);

			av_background_interactionE += Eedge->dot( cur_weights );

		} //loop over interacting partners of particular residue
	} //loop over all residues of this region

	av_background_interactionE /= ( interacting_neighbors.size() + this->length() );

	backgroundE = av_background_interactionE;

	return per_res_design_target_interactionE;
	//return av_background_interactionE + (per_res_design_target_interactionE * basic::options::option[basic::options::OptionKeys::enzdes::lig_packer_weight] );

} //extract_region_designability_score

void
EnzdesFlexibleRegion::apply_ranked_fragment(
	core::pose::Pose & pose,
	core::Size frag_rank
)
{
	if( frag_rank > frag_designabilities_.size() ) utility_exit_with_message("Trying to access a fragment that is not ranked.");

	//std::list< SizeRealPair >::const_iterator list_it = frag_designabilities_.begin();
	//for( core::Size i = 2; i <= frag_rank; ++i) ++list_it;
	//tr << "putting ranked fragment " << frag_rank << " of region " << this->index_ << " into pose, translates to fragment " << frag_designabilities_[frag_rank].first << std::endl;

	this->apply( frag_designabilities_[frag_rank].first, pose );

} //apply_ranked_fragment

core::Real
EnzdesFlexibleRegion::get_region_mm_bend_score(
	core::pose::Pose const & pose ) const
{

	core::Real to_return(0.0);

	for( core::Size i = this->start(); i <= this->end(); ++i){
		to_return += pose.energies().residue_total_energies( i )[ core::scoring::mm_bend ];
	}

	return to_return;
} // get_region_mm_bend_score

bool
EnzdesFlexibleRegion::examine_new_loopconf(
	core::pose::Pose const & pose,
	core::pose::Pose & template_pose,
	utility::vector1< core::pose::PoseOP > & compare_poses,
	utility::vector1< core::Real > &  rmsd_to_native
)
{
	//static Size count_examined( 0 );
	//Size const sought_loop_id = 5;

	runtime_assert( compare_poses.size() > 0 );
	runtime_assert( template_pose.total_residue() == this->length() + 1);

	core::fragment::FragDataOP newfrag = this->fragment_ptr( 1 )->clone();
	if( !newfrag->steal( pose, *this ) ) utility_exit_with_message("unknown error when trying to steal fragment from pose for examination.");

	if( !newfrag->is_valid() ) utility_exit_with_message("unknown error when trying to steal fragment from pose for examination. fragment not valid");

	newfrag->apply( template_pose, 2, this->length() );
   /// FIX C and O on the last residue
	template_pose.set_phi( template_pose.total_residue(), pose.phi( positions_[ positions_.size() ] ) );
	template_pose.set_psi( template_pose.total_residue(), pose.psi( positions_[ positions_.size() ] ) );
	template_pose.set_omega( template_pose.total_residue(), pose.omega( positions_[ positions_.size() ] ) );

	//std::cout << "template: O " << template_pose.residue( template_pose.total_residue() ).xyz( "O" ).x();
	//std::cout << " " << template_pose.residue( template_pose.total_residue() ).xyz( "O" ).y();
	//std::cout << " " << template_pose.residue( template_pose.total_residue() ).xyz( "O" ).z() << std::endl;

	//std::cout << "regular: O " << pose.residue( positions_[ positions_.size() ] ).xyz( "O" ).x();
	//std::cout << " " << pose.residue( positions_[ positions_.size() ] ).xyz( "O" ).y();
	//std::cout << " " << pose.residue( positions_[ positions_.size() ] ).xyz( "O" ).z() << std::endl;

	template_pose.set_xyz(
		core::id::AtomID( template_pose.residue( template_pose.total_residue() ).atom_index( "O" ), template_pose.total_residue() ),
		pose.residue( positions_[ positions_.size() ] ).xyz( "O" ) );


	core::Real similarity_to_native = core::scoring::rmsd_no_super( *compare_poses[1], template_pose, core::scoring::is_protein_backbone );
		//if ( index_ == sought_loop_id ) {
		//	std::cout << "loop " << count_examined << " distance to native: " << similarity_to_native << std::endl;
		//}

	bool frag_close_to_native = ( similarity_to_native <= target_proximity_to_native_conformation_ );

	if( !frag_close_to_native ) return false;


	bool frag_unique(true);
	//std::cerr << "gotta compare against " << compare_poses.size() << "unique poses " << std::endl;
	Size count( 0 );
	for( utility::vector1< core::pose::PoseOP >::const_iterator pose_it = compare_poses.begin(); pose_it != compare_poses.end(); ++pose_it){

		count++;
		//std::cout << "rmsd with prepose " << count << " is " << core::scoring::rmsd_no_super( **pose_it, template_pose, core::scoring::is_protein_CA ) << ", biggest CA/CB dist is " << core::scoring::biggest_residue_deviation_no_super( **pose_it, template_pose, core::scoring::is_protein_CA_or_CB ) << std::endl;

		if( core::scoring::biggest_residue_deviation_no_super( **pose_it, template_pose, core::scoring::is_protein_CA_or_CB ) < basic::options::option[basic::options::OptionKeys::enzdes::min_cacb_deviation] ){
			frag_unique = false;
			return false;
		}
	}


	bool frag_close_to_another_fragment( true );
	if ( basic::options::option[basic::options::OptionKeys::enzdes::max_bb_deviation ].user() )
	{
		// Find the smallest maximum deviation the template pose has to all other accepted poses.
		core::Real smallest_max_dist = core::scoring::biggest_residue_deviation_no_super(
			*compare_poses[1], template_pose, core::scoring::is_protein_backbone );
		for ( Size ii = 2; ii <= compare_poses.size(); ++ii ) {
			core::Real ii_super = core::scoring::biggest_residue_deviation_no_super(
				*compare_poses[ii], template_pose, core::scoring::is_protein_backbone );
			if ( ii_super < smallest_max_dist ) {
				smallest_max_dist = ii_super;
			}
		}

		frag_close_to_another_fragment = smallest_max_dist < target_proximity_to_other_conformations_;

	}


	//if( compare_poses.size() == 1 ) { template_pose.dump_pdb("template_1.pdb"); std::cerr << " dumping " << std::endl; }
	if ( frag_unique && frag_close_to_native && frag_close_to_another_fragment ){
		this->add_fragment( newfrag );
		rmsd_to_native.push_back( core::scoring::rmsd_no_super( *compare_poses[1], template_pose, core::scoring::is_protein_backbone ) );
		compare_poses.push_back( new core::pose::Pose( template_pose ) );

		return true;

		//	if ( basic::options::option[ basic::options::OptionKeys::enzdes::dump_loop_samples ]() != "no" ) {
		//	template_pose.dump_pdb("loopreg_"+utility::to_string( index_ )+"_"+utility::to_string( compare_poses.size() )+".pdb");
		//	}
	}

	return false;

 /*else if ( ! pose_close_to_native ) {
		//std::cout << "rejected pose, unique, but not near enough to the native: " <<
		//	core::scoring::biggest_residue_deviation_no_super( *compare_poses[1], template_pose, core::scoring::is_protein_CA )
		//	<< std::endl;
		if ( index_ == sought_loop_id )  {
			std::cout << "dumping rejected_pose_not_close_to_native_" + utility::to_string( index_ ) + "_" + utility::to_string( ++count_examined ) + ".pdb" << std::endl;
			template_pose.dump_pdb("rejected_pose_not_close_to_native_" + utility::to_string( index_ ) + "_" + utility::to_string( count_examined ) + ".pdb" );
		}
	} else if ( ! pose_unique ) {
		if ( index_ == sought_loop_id ) {
			std::cout << "dumping rejected_pose_not_unique_" + utility::to_string( index_ ) + "_" + utility::to_string( ++count_examined ) + ".pdb" << std::endl;
			template_pose.dump_pdb("rejected_pose_not_unique_" + utility::to_string( index_ ) + "_" + utility::to_string( count_examined ) + ".pdb" );
		}
	} else {
		if ( index_ == sought_loop_id ) {
			std::cout << "dumping rejected_pose_not_near_anothre_frag_" + utility::to_string( index_ ) + "_" + utility::to_string( ++count_examined ) + ".pdb" << std::endl;
			template_pose.dump_pdb("rejected_pose_not_near_anothre_frag_" + utility::to_string( index_ ) + "_" + utility::to_string( count_examined ) + ".pdb" );
		}
	}*/

	/*
	if(this->nr_frags() == 6 ){
		template_pose.dump_pdb("looppose_"+utility::to_string( index_ )+"_5.pdb");
		std::cerr << "dumped looppose_5.pdb" << std::endl;
	}
	if(this->nr_frags() == 7 ){
		template_pose.dump_pdb("looppose_"+utility::to_string( index_ )+"_6.pdb");
		std::cerr << "dumped looppose_6.pdb" << std::endl;
	}
	*/

} //examine_new_loopconf



/// @details minimize the backbone of this pose over the fragment residues, including the
/// @details bond angles around Calpha if desired. NOTE: CA ANGLE MINIMIZATION UNTESTED
/// a chainbreak at the end of the region will be introduced. chainbreak weight will be
/// set to 100 to make sure that it stays closed. in case the chainbreak score is higher
/// after the min then before, function returns false. this might be the case when
/// the input conformation has bad clashes, which is the case e.g. for some conformations
/// that come out of remodel
bool
EnzdesFlexibleRegion::minimize_region(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scorefxn,
	std::set< core::Size > const & chi_to_move,
	bool const including_CA_angles,
	core::Real min_tolerance
)
{

	//1. we need to setup a correct foldtree, so no downstream residues get affected
	core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, this->end() );
	core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, this->end() + 1 );
	//core::Real mm_bend_bef = get_region_mm_bend_score( pose );

	core::kinematics::FoldTree old_fold_tree = pose.fold_tree();
	core::Size old_njump = old_fold_tree.num_jump();
	core::kinematics::FoldTree temp_fold_tree;
	core::kinematics::FoldTree const & f_const = old_fold_tree;
	//tr << "regmindebug start foldtree " << old_fold_tree << std::endl;
	core::scoring::ScoreFunctionOP min_scorefxn = scorefxn->clone();
	min_scorefxn->set_weight( core::scoring::chainbreak, 100.0 );

	//fold tree approach
	//the segment will be covered by an edge from start -1 to end
	//non-jump edges spanning the segment will be replaced by two edges
	//jump edges spanning the segment will be left untouched
	//non-jump edges going into the segment will go either until start or end
	//jump ed
	bool backward_into_seg(false), backward_outof_seg(false);
	utility::vector1< std::pair<core::Size,int> > jump_origins, jump_destinations;

	//find the edge that spans the end of this segment
	//tr << "regmindebug setting foldtree for region from " << this->start() << " to " << this->end() << std::endl;
	for( core::kinematics::FoldTree::const_iterator e = f_const.begin(); e != f_const.end(); ++e ){
		bool is_jump( e->is_jump() ), backward( e->start() > e->stop() ), start_in_seg( e->start() >= (int)this->start() && e->start() <= (int)this->end() ), stop_in_seg( e->stop() >= (int)this->start() && e->stop() <= (int)this->end() );
		bool span( backward ? ( (e->start() > (int)this->end()) && (e->stop() < (int)this->start()) ) : ( (e->start() < (int)this->start()) && (e->stop() > (int)this->end()) ) );
		//tr << "regmindebug dealing with edge from " << e->start() << " to " << e->stop() << " with label " << e->label() << " backward is " << backward << ", span is " << span << ", start_in_seg is " << start_in_seg << ", stop in seg is " << stop_in_seg << ", is_jump is " << is_jump << std::endl;

		//edges only in the segment will be discarded
		if( start_in_seg && stop_in_seg ){
			if( is_jump ) old_njump--;
			continue;
		}

		if( is_jump && start_in_seg ) jump_destinations.push_back( std::pair<core::Size, int>(e->stop(), e->label()) );
		else if( is_jump && stop_in_seg ) jump_origins.push_back( std::pair< core::Size, int>(e->start(), e->label() ) );
		else if( !is_jump && start_in_seg ){
			if( backward ){
				temp_fold_tree.add_edge( this->start(), e->stop(), e->label() );
				backward_outof_seg = true;
			}
			else temp_fold_tree.add_edge( this->end() + 1, e->stop(), e->label() );
		}
		else if( !is_jump && stop_in_seg ){
			if( backward ){
				temp_fold_tree.add_edge( e->start(), this->end() + 1, e->label() );
				backward_into_seg = true;
			}
			else temp_fold_tree.add_edge( e->start(), this->start(), e->label() );
		}
		else if (!is_jump && span ){
			if(backward ){
				//tr << "regmindebug dealing with backward span " << std::endl;
				temp_fold_tree.add_edge( e->start(), this->end() + 1, e->label() );
				temp_fold_tree.add_edge( this->start(), e->stop(), e->label() );
				backward_outof_seg = true;
				backward_into_seg = true;
			}
			else{
				//tr << "regmindebug dealing with forward span " << std::endl;
				temp_fold_tree.add_edge( e->start(), this->start(), e->label() );
				temp_fold_tree.add_edge( this->end() + 1, e->stop(), e->label() );
			}
		}

		else{
			temp_fold_tree.add_edge( *e );
		}
	} // loop over edges
	temp_fold_tree.add_edge( this->start(), this->end(), -1);  // add edge for segment
	core::Size jump_focus( this->start() );
	if( backward_into_seg && backward_outof_seg ){
		jump_focus = this->end() + 1;
		temp_fold_tree.add_edge(  core::kinematics::Edge(this->end() +1, this->start(), old_njump + 1 ,"CA", "CA", false ) );
	}
	else{
		//tr << "regmindebug adding jump across segment " << std::endl;
		temp_fold_tree.add_edge(  core::kinematics::Edge(this->start(), this->end() +1, old_njump + 1 ,"CA", "CA", false ) );
	}

	for( utility::vector1< std::pair<core::Size, int> >::const_iterator jump_o_it = jump_origins.begin(); jump_o_it != jump_origins.end(); ++jump_o_it )  temp_fold_tree.add_edge(  core::kinematics::Edge(jump_o_it->first, jump_focus, jump_o_it->second ,"CA", "CA", false ) );

	for( utility::vector1< std::pair< core::Size,int> >::const_iterator jump_d_it = jump_destinations.begin(); jump_d_it != jump_destinations.end(); ++jump_d_it )  temp_fold_tree.add_edge(  core::kinematics::Edge( jump_focus, jump_d_it->first, jump_d_it->second ,"CA", "CA", false ) );


	temp_fold_tree.delete_extra_vertices();
	//tr << "regmindebug new foldtree " << temp_fold_tree << std::endl;

	if( !temp_fold_tree.check_fold_tree() ) {
		utility_exit_with_message("Invalid fold tree after trying to set up for flexbb ca angle min");
	}

	//std::cerr << "frigging new fold tree at end has " << temp_fold_tree.num_cutpoint() << " cutpoints." << std::endl << temp_fold_tree << std::endl;
	pose.fold_tree( temp_fold_tree );

	//get the cur chainbreak
	core::Real totE_start = (*scorefxn)( pose );
	core::Real cbE_start = pose.energies().total_energies()[ core::scoring::chainbreak ];

	//2. now set up the correct movemap, including the CA bond angles
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap();
	movemap->clear();

	for( core::Size i = this->start(); i <= this->end(); ++i){

		movemap->set_chi( i, true );
		movemap->set_bb( i, true );

		if( including_CA_angles ){
			core::id::AtomID c_id ( pose.residue_type(i).atom_index("C"), i );
			core::id::DOF_ID ca_dof( c_id, core::id::PHI );
			movemap->set( ca_dof, true );
		}

	}

	for( std::set< core::Size >::const_iterator chi_it = chi_to_move.begin(); chi_it != chi_to_move.end(); ++chi_it ){
		movemap->set_chi( *chi_it, true );
	}

	//core::scoring::ScoreFunctionOP trial_score = new core::scoring::ScoreFunction( *scorefxn );
	//trial_score->set_weight( core::scoring::mm_bend, 0.0 );
	//(*trial_score)(pose);
	(*min_scorefxn)(pose);
	protocols::simple_moves::MinMoverOP dfpMinTightTol = new protocols::simple_moves::MinMover( movemap, min_scorefxn, "dfpmin_armijo_nonmonotone_atol", min_tolerance, true  );
	dfpMinTightTol->apply(pose);
	core::Real totE_end = (*scorefxn)(pose);
	core::Real cbE_end = pose.energies().total_energies()[ core::scoring::chainbreak ];

	//finally put back the old fold tree
	pose.fold_tree( old_fold_tree );
	core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_LOWER, this->end() );
	core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_UPPER, this->end() + 1 );

	(*scorefxn)(pose);
	//core::Real mm_bend_aft = get_region_mm_bend_score( pose );

	if( cbE_end > ( cbE_start + 0.1) ){
		tr << "FlexRegion minimization failed: cbE_start " << cbE_start << ", cbE_end " << cbE_end << "." << std::endl;
		return false;
	}

	if( totE_end > totE_start ){
		tr << "FlexRegion minimization failed: totE_start " << totE_start << ", totE_end " << totE_end << "." << std::endl;
		return false;
	}
	//lastly, we should clear the frags of this region
	this->clear();
	if( !this->steal( pose ) ) utility_exit_with_message("Unknown error when trying to save fragdata after minimizing region.");
	runtime_assert( this->nr_frags() == 1 );
	return true;

	//tr << "MMBEND before minimization, mm_bend score was " << mm_bend_bef << ", while after it is " << mm_bend_aft << std::endl;
} //minimize_region


core::Real
EnzdesFlexibleRegion::get_region_total_score( core::pose::Pose const & pose ) const
{

	using namespace core::scoring;

	core::Real to_return(0.0);

	for( core::Size i = this->start(); i <= this->end(); ++i ){
		to_return += pose.energies().residue_total_energy( i );
	}

	return to_return;
}

bool
EnzdesFlexibleRegion::remap_resid(
	core::pose::Pose const & pose,
	core::id::SequenceMapping const& smap
)
{
	core::Size newstart( smap[ this->start()] ), newend( smap[ this->end()] );
	core::Size newlength( newend - newstart + 1 );

	//only shifting?
	if( newlength == this->length() ){
		bool to_return = Super::align( smap );

		//super class handles checks for 0 in sequence mapping
		if( to_return ){

			for( core::Size i = 1; i <= positions_.size(); ++i){
				positions_[i] =  smap[ positions_[i] ];
			}
		}
		return to_return;
	}

	//if length of this region changed, need to proceed differently
	if( (newstart == 0) || (newend == 0 ) || !is_continuous() ) return false;
	positions_.clear();

	this->init_length( newstart, newend, newlength );
	for ( core::Size i = newstart; i <= newend; ++i ) positions_.push_back( i );

	this->clear();
	native_conf_ = assemble_enzdes_fragdata( pose );
	if( add_fragment( native_conf_ ) != 1 ) return false;
	frag_designabilities_.clear();
	return true;
}

/// @brief requires that the pose was scored
std::set< core::Size >
EnzdesFlexibleRegion::get_10A_neighbors(
	core::pose::Pose const & pose
) const {

	std::set< core::Size > ten_A_neighbors;
	ten_A_neighbors.clear();
	//shiat, seems to be complicated to have the scorefxn build the 12a neighbor graph,
	//so we'll take the ten a one for now...
	core::scoring::TenANeighborGraph const & cur_graph = pose.energies().tenA_neighbor_graph();

	for( core::Size i = this->start(); i <= this->end(); ++i ){

		for( core::graph::EdgeListConstIterator graph_it = cur_graph.get_node( i )->const_edge_list_begin();
				 graph_it != cur_graph.get_node( i )->const_edge_list_end(); ++graph_it){

			core::Size other_res = (*graph_it)->get_other_ind( i );
			if( !this->contains_seqpos( other_res )
				&& ( ten_A_neighbors.find( other_res) == ten_A_neighbors.end() ) )
				{
					ten_A_neighbors.insert( other_res );
				}
		}
	}
	return ten_A_neighbors;
} //determine_10A_neighbors(




void
EnzdesFlexBBProtocol::test_flexbb_rotamer_sets(
	core::pose::Pose const & pose,
	core::pack::task::PackerTaskCOP task
){

	using namespace flexpack::rotamer_set;

	FlexbbRotamerSetsOP flexset = new FlexbbRotamerSets( task );

	std::cerr << "FlexbbRotamerSet test: done initialising set." << std::endl;

	flexset->set_frames( pose, flex_regions_ );

	std::cerr << "FlexbbRotamerSet test: done setting frames." << std::endl;

	core::graph::GraphOP flex_graph = flexset->flexpack_neighbor_graph( pose, *scorefxn_, task );

	std::cerr << "FlexbbRotamerSet test: done setting up flexpack neighbor graph." << std::endl;

	flexset->build_rotamers( pose, *scorefxn_, *flex_graph );

	std::cerr << "FlexbbRotamerSet test: done building rotamers." << std::endl;

	flexset->dump_pdbs( pose, "flextest" );

} // test_flexbb_rotamer_sets


} //namespace enzdes
} //namespace protocols
