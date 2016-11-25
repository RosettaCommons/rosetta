// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/enzdes/EnzdesRemodelProtocol.cc
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, march 2009


#include <devel/enzdes/EnzdesRemodelProtocol.hh>
#include <devel/enzdes/EnzdesRemodelMoverCreator.hh>

#include <protocols/enzdes/AddorRemoveCsts.hh>
#include <protocols/enzdes/enzdes_util.hh>
//#include <protocols/enzdes/SecondaryMatchProtocol.hh> //for secmatch
#include <protocols/toolbox/match_enzdes_util/EnzConstraintParameters.hh> //for secmatch
#include <protocols/toolbox/match_enzdes_util/EnzdesLoopsFile.hh>
#include <protocols/toolbox/match_enzdes_util/EnzCstTemplateRes.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
//#include <protocols/enzdes/EnzConstraintIO.hh>

//stuff for rebuilding
#include <protocols/ligand_docking/ligand_functions.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/Interval.hh>
#include <protocols/forge/build/SegmentRebuild.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/remodel/ResidueVicinityRCG.hh>
#include <protocols/forge/constraints/InverseRotamersRCG.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh> //unnecessary include
#include <core/chemical/ResidueConnection.hh>

//#include <protocols/enzdes/EnzdesBaseProtocol.hh>
//#include <protocols/enzdes/EnzConstraintIO.hh>
//#include <protocols/enzdes/EnzdesMovers.hh>
//#include <core/chemical/ResidueType.hh>
//#include <basic/options/option.hh>
#include <core/pack/task/PackerTask_.hh>
//#include <core/pose/Pose.hh>
//#include <core/scoring/ScoreFunction.hh>
//#include <protocols/ligand_docking/LigandDockProtocol.hh>
//#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#if defined GL_GRAPHICS
#include <protocols/viewer/viewers.hh>
#endif
#include <protocols/match/MatcherMover.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/fragment/util.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/id/SequenceMapping.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/BoundConstraint.hh> //need function in this file
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>

#include <protocols/simple_filters/ScoreCutoffFilter.hh>
#include <protocols/simple_filters/PackerNeighborGraphFilter.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/toolbox/pose_metric_calculators/NonlocalContactsCalculator.hh>

// option key includes

#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// numeric includes
#include <numeric/random/random.hh>

//utility includes
#include <utility/tag/Tag.hh>

#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <platform/types.hh>
#include <utility/graph/Graph.hh>
#include <core/kinematics/FoldTree.hh>
#include <utility/string_util.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace devel {
namespace enzdes {

static THREAD_LOCAL basic::Tracer tr( "devel.enzdes.EnzdesRemodelProtocol" );

EnzdesRemodelProtocol::EnzdesRemodelProtocol() : EnzdesFlexBBProtocol() {
	this->reduced_sfxn_->set_weight(core::scoring::backbone_stub_constraint, 1.0 );
}

EnzdesRemodelProtocol::~EnzdesRemodelProtocol() = default;

void
EnzdesRemodelProtocol::apply(
	core::pose::Pose & pose
){

	using namespace protocols::moves;
	using namespace core::pack::task;

	// Scoring function already set up by superclass
	//tr.Info << "starting apply function..." << std::endl;

#if defined GL_GRAPHICS
	protocols::viewer::add_conformation_viewer( pose.conformation(), "EnzdesRemodel" );
#endif

	//score pose to make sure everything is initialised correctly
	(*scorefxn_)( pose );

	//set the native pose if requested
	if ( ! basic::options::option[basic::options::OptionKeys::in::file::native].user() ) {
		core::pose::PoseOP natpose( new core::pose::Pose( pose ) );
		(*scorefxn_)( *natpose );
		this->set_native_pose( natpose );
	}

	//set up constraints (read cstfile, do mapping, etc, then add to pose)
	if ( basic::options::option[basic::options::OptionKeys::enzdes::cstfile].user() ) {
		enable_constraint_scoreterms();
		setup_enzdes_constraints( pose, basic::options::option[basic::options::OptionKeys::enzdes::remodel_secmatch]  );
		(*scorefxn_)( pose );
	}

	//create packer task (read resfile, etc)
	PackerTaskOP orig_task = create_enzdes_pack_task( pose );

	//cst opt stage, if demanded
	if ( basic::options::option[basic::options::OptionKeys::enzdes::cst_opt] ) {

		tr.Info << "starting cst_opt minimization..." << std::endl;
		cst_minimize(pose, orig_task, true);
		(*scorefxn_)( pose );
		tr.Info << "done cst_opt minimization." << std::endl;
		orig_task =  create_enzdes_pack_task( pose );

	}

	determine_flexible_regions( pose, orig_task );

	core::id::SequenceMappingCOP seq_mapping;
	PackerTaskOP mod_task = modified_task( pose, *orig_task );

	for ( core::Size regcount = 1; regcount <= flex_regions_.size(); ++regcount ) {
		if ( !flex_regions_[ regcount ]->remodelable() ) continue;

		EnzdesRemodelMover remodel_mover( utility::pointer::dynamic_pointer_cast< EnzdesFlexBBProtocol >(get_self_ptr()), mod_task, flex_regions_[ regcount ] );

		remodel_mover.apply( pose );

		if ( remodel_mover.get_last_move_status() != protocols::moves::MS_SUCCESS ) {
			tr << "WARNING: remodeling of region " << regcount << " failed." << std::endl;
		}
		seq_mapping = remodel_mover.get_seq_mapping();
	}

	(*scorefxn_)(pose);

	//remap task according to indels after remodelling
	tr << "task remapping " << std::endl;
	mod_task -> remap_residue_level_tasks( seq_mapping, pose  );

	if ( basic::options::option[basic::options::OptionKeys::enzdes::cst_design] ) {

		core::Size design_min_cycles = basic::options::option[basic::options::OptionKeys::enzdes::design_min_cycles];

		tr.Info << "starting cst_design, " << design_min_cycles << " cycles of design/minimization ... " << std::endl;

		bool favor_native_res(false);
		if ( basic::options::option[basic::options::OptionKeys::enzdes::favor_native_res].user() ) favor_native_res = true;

		tr.Info << "starting design" << std::endl;
		enzdes_pack( pose, mod_task, scorefxn_, design_min_cycles, basic::options::option[basic::options::OptionKeys::enzdes::cst_min], false, favor_native_res );

	} //if cst_design

	(*scorefxn_)(pose);

	//remap PDB remark info for enzyme design constraints
	protocols::toolbox::match_enzdes_util::EnzConstraintIOCOP cstio( protocols::enzdes::enzutil::get_enzcst_io( pose ) );
	if ( cstio ) {
		for ( core::Size ii = 1; ii <=  cstio->num_mcfi_lists(); ++ ii ) {
			cstio-> enz_cst_params(ii) -> update_pdb_remarks( pose );
		}
	}
} //apply function

std::string
EnzdesRemodelProtocol::get_name() const {
	return "EnzdesRemodelProtocol";
}

// XRW TEMP std::string
// XRW TEMP EnzdesRemodelMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return EnzdesRemodelMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP EnzdesRemodelMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new EnzdesRemodelMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP EnzdesRemodelMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "EnzdesRemodelMover";
// XRW TEMP }


EnzdesRemodelMover::EnzdesRemodelMover()
:
	protocols::moves::Mover(),
	enz_prot_(/* NULL */),
	flex_region_( /* NULL */ ),
	remodel_trials_( basic::options::option[basic::options::OptionKeys::enzdes::remodel_trials] ),
	remodel_secmatch_(basic::options::option[basic::options::OptionKeys::enzdes::remodel_secmatch] ),
	reinstate_initial_foldtree_(false),
	keep_existing_aa_identities_(false),
	region_to_remodel_(1),
	start_to_current_smap_(/* NULL */),
	include_existing_conf_as_invrot_target_(false),
	ss_similarity_probability_( 1.0 - basic::options::option[basic::options::OptionKeys::enzdes::remodel_aggressiveness] ),
	max_allowed_score_increase_(0.0)
{
	predesign_filters_ = protocols::filters::FilterCollectionOP( new protocols::filters::FilterCollection() );
	postdesign_filters_ = protocols::filters::FilterCollectionOP( new protocols::filters::FilterCollection() );

	vlb_sfxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "remodel_cen" );
	//turn on the constraint weights in the remodel scorefunction for our purposes
	vlb_sfxn_->set_weight(core::scoring::coordinate_constraint, 1.0 );
	vlb_sfxn_->set_weight(core::scoring::atom_pair_constraint, 1.0);
	vlb_sfxn_->set_weight(core::scoring::angle_constraint, 1.0 );
	vlb_sfxn_->set_weight(core::scoring::dihedral_constraint, 1.0 );
	vlb_sfxn_->set_weight(core::scoring::backbone_stub_constraint, 1.0 );
	non_remodel_match_pos_.clear();
	init_aa_.clear();
	user_provided_ss_.clear();
}

EnzdesRemodelMover::EnzdesRemodelMover( EnzdesRemodelMover const & /*other*/ ) = default;

EnzdesRemodelMover::EnzdesRemodelMover(
	protocols::enzdes::EnzdesFlexBBProtocolOP enz_prot,
	core::pack::task::PackerTaskCOP task,
	protocols::enzdes::EnzdesFlexibleRegionCOP flex_region
):
	protocols::moves::Mover(),
	enz_prot_( enz_prot ),
	flex_region_( flex_region ),
	remodel_trials_( basic::options::option[basic::options::OptionKeys::enzdes::remodel_trials] ),
	remodel_secmatch_( basic::options::option[basic::options::OptionKeys::enzdes::remodel_secmatch] ),
	reinstate_initial_foldtree_(false),
	keep_existing_aa_identities_(false),
	region_to_remodel_(1),
	start_to_current_smap_(/* NULL */),
	include_existing_conf_as_invrot_target_(false),
	ss_similarity_probability_( 1.0 - basic::options::option[basic::options::OptionKeys::enzdes::remodel_aggressiveness] )
{

	using namespace basic::options;

	//if( flex_region->contains_catalytic_res() ) {
	// utility_exit_with_message("Requesting remodel of a region containing catalytic residues. This is not possible yet.");
	//}

	runtime_assert( flex_region == enz_prot->enz_flexible_region( flex_region->index() ) );

	predesign_filters_ = protocols::filters::FilterCollectionOP( new protocols::filters::FilterCollection() );
	postdesign_filters_ = protocols::filters::FilterCollectionOP( new protocols::filters::FilterCollection() );

	set_task( task );

	for ( core::Size i = flex_region_->start(); i <= flex_region_->stop(); ++i ) remodel_positions_.insert( i );

	tr << "ss_sim got set to " << ss_similarity_probability_ << std::endl;

	vlb_sfxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "remodel_cen" );

	//turn on the constraint weights in the remodel scorefunction for our purposes
	vlb_sfxn_->set_weight(core::scoring::coordinate_constraint, 1.0 );
	vlb_sfxn_->set_weight(core::scoring::atom_pair_constraint, 1.0);
	vlb_sfxn_->set_weight(core::scoring::angle_constraint, 1.0 );
	vlb_sfxn_->set_weight(core::scoring::dihedral_constraint, 1.0 );
	vlb_sfxn_->set_weight(core::scoring::backbone_stub_constraint, 1.0 );

	non_remodel_match_pos_.clear();
	user_provided_ss_.clear();
	init_aa_.clear();

} //enzdes remodel mover constructor

EnzdesRemodelMover::~EnzdesRemodelMover()= default;

protocols::moves::MoverOP
EnzdesRemodelMover::clone() const{
	return protocols::moves::MoverOP( new EnzdesRemodelMover( *this ) );
}


void
EnzdesRemodelMover::set_task( core::pack::task::PackerTaskCOP task )
{

	orig_task_ = task;

	other_design_positions_.clear();

	for ( core::Size i = 1; i <= task->total_residue(); ++i ) {

		if ( flex_region_->contains_seqpos( i ) ) continue;

		if ( orig_task_->residue_task( i ).being_designed() ) {
			other_design_positions_.insert( i );
		} else if ( orig_task_->residue_task( i ).being_packed() ) {
			other_repack_positions_.insert( i );
		}

	}

} //set task

void EnzdesRemodelMover::set_user_provided_ss( utility::vector1< std::string > const & user_ss )
{
	user_provided_ss_ = user_ss;
}


void
EnzdesRemodelMover::apply(
	core::pose::Pose & pose
)
{

	//in case necessary stuff hasn't been set yet
	if ( !enz_prot_ || !orig_task_ ) {
		this->initialize( pose );
	}
	//(*enz_prot_->scorefxn() )( pose );

	this->set_last_move_status( protocols::moves::FAIL_RETRY );

	this->set_native_pose( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( pose ) ) ) );

	//note: the pose gets reduced to poly-ala at all design positions here
	examine_initial_conformation( pose );

	setup_cached_observers( pose );
	//static core::Size num_nstruct(0);
	//num_nstruct++;

	for ( core::Size i = 1; i <= remodel_trials_; ++i ) {

		//first, get the new secondary structure of the loop
		std::string secstruct = generate_secstruct_string( pose );
		tr << "Score bef remodel is " << (*(enz_prot_->reduced_scorefxn()))(pose) << std::endl;

		//might be helpful to change the ligand around a bit
		//apply_random_lowE_ligconf( pose );
		//get new loop
		if ( !remodel_pose( pose, secstruct ) ) {
			tr << "Region " << region_to_remodel_ << " remodel trial " << i << ": fail remodel." << std::endl;
			continue;
		}
		//pose.dump_pdb("remdebug_aftremodel_"+utility::to_string( num_nstruct )+".pdb" );
		tr << "Score after remodel and before refine is " << (*(enz_prot_->reduced_scorefxn()))(pose) << std::endl;
		//some of the loops have nasty  clashes, so let's do some fa refinement
		if ( !refine_pose( pose ) ) {
			tr << "Region " << region_to_remodel_ << " remodel trial " << i << ": fail remodel." << std::endl;
			continue;
		}
		//pose.dump_pdb("remdebug_aftrefine_"+utility::to_string( num_nstruct )+".pdb" );
		tr << "Score after refine is " << (*(enz_prot_->reduced_scorefxn()))(pose) << std::endl;

		//apply predesign filters
		if ( !predesign_filters_->apply( pose ) ) {
			tr << "Region " << region_to_remodel_ << " remodel trial " << i << ": fail predesign filters." << std::endl;
			continue;
		}
		//if there are missing constraints in the pose, let's start the secondary matcher
		if ( !secmatch_after_remodel( pose ) ) {
			tr << "Region " << region_to_remodel_ << " remodel trial " << i << ": fail secmatch." << std::endl;
			continue;
		}
		//for now we will break as soon as we have the first new loop
		this->set_last_move_status( protocols::moves::MS_SUCCESS );
		tr << "Remodel mover for region " << region_to_remodel_ << " success with secstruct string " << secstruct << " in trial " << i <<"." << std::endl;
		break;

		//do design+min + unconstrained repack
		//apply post design filters
	}

	remove_cached_observers( pose );

	//maybe better to get rid of the task, to prepare for next call
	orig_task_ = nullptr;
}

// XRW TEMP std::string
// XRW TEMP EnzdesRemodelMover::get_name() const {
// XRW TEMP  return "EnzdesRemodelMover";
// XRW TEMP }

void
EnzdesRemodelMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,// datamap,
	Filters_map const &, //filters,
	protocols::moves::Movers_map const &, //movers,
	Pose const & //pose
)
{
	if ( tag->getName() != "EnzdesRemodelMover" ) {
		tr << "EnzdesRemodelMover received incompatible Tag " << tag << std::endl;
		runtime_assert(false);
		return;
	}

	if ( tag->hasOption("remodel_secmatch") ) {
		remodel_secmatch_ = tag->getOption<bool>( "remodel_secmatch", 1 );
	}
	if ( tag->hasOption("reinstate_foldtree") ) {
		reinstate_initial_foldtree_ = tag->getOption<bool>( "reinstate_foldtree", 1 );
	}
	if ( tag->hasOption("remodel_region") ) {
		region_to_remodel_ = tag->getOption<core::Size>( "remodel_region", 1 );
	}
	if ( tag->hasOption("include_existing_conf_as_invrot_target") ) {
		include_existing_conf_as_invrot_target_ = tag->getOption<bool>("include_existing_conf_as_invrot_target",1);
	}
}

void
EnzdesRemodelMover::initialize(
	core::pose::Pose & pose )
{

	if ( !enz_prot_ ) {
		enz_prot_ = protocols::enzdes::EnzdesFlexBBProtocolOP( new protocols::enzdes::EnzdesFlexBBProtocol() );
		if ( enz_prot_->reduced_scorefxn()->has_zero_weight( core::scoring::backbone_stub_constraint ) ) {
			enz_prot_->reduced_scorefxn()->set_weight( core::scoring::backbone_stub_constraint, 1.0 );
		}
	}
	if ( !orig_task_ ) orig_task_ = enz_prot_->create_enzdes_pack_task( pose );

	enz_prot_->determine_flexible_regions( pose, orig_task_ );
	//temporary hardcoded: only one region to be remodeled
	flex_region_ = enz_prot_->enz_flexible_region( region_to_remodel_ );

	pose.update_residue_neighbors(); // Neighbors needed for modified_task
	this->set_task( enz_prot_->modified_task( pose, *orig_task_ ) );
}

bool
EnzdesRemodelMover::remodel_pose(
	core::pose::Pose & pose,
	std::string secstruct
)
{
	using namespace protocols::forge;
	using namespace protocols::forge::build;
	using namespace protocols::forge::components;
	using protocols::moves::MS_SUCCESS;

	//debug
	//return true;

	enz_prot_->remove_enzdes_constraints( pose, false );

	Interval interval( flex_region_->start(), flex_region_->stop() );

	SegmentRebuildOP build_instr( new SegmentRebuild( interval, secstruct ) );

	BuildManager bmanager;

	bmanager.add( build_instr );

	VarLengthBuildOP vlb_op( new VarLengthBuild(bmanager) );
	VarLengthBuild & vlb = *vlb_op;

	//set the scorefunction that contains our constrained weights
	vlb.scorefunction( *vlb_sfxn_ );
	vlb.ignore_cmdline_enzdes_cstfile(true);

	setup_rcgs( pose, vlb );

	vlb.apply( pose );

	//fix the pdbinfo
	core::pose::renumber_pdbinfo_based_on_conf_chains( pose );
	pose.pdb_info()->obsolete( false );


	//score the pose in our relevant scorefunction
	(*enz_prot_->reduced_scorefxn())( pose );

	if ( vlb.get_last_move_status() != MS_SUCCESS ) {

		enz_prot_->add_pregenerated_enzdes_constraints( pose );
		return false;

	}

	if ( secstruct.size() != flex_region_->length() ) process_length_change( pose, vlb.manager().sequence_mapping() );

	if ( keep_existing_aa_identities_ ) {

		bool oddlength ( init_aa_.size() % 2 == 0 ? false : true );
		Size num_aa_each_side_to_replace( init_aa_.size()/2 );
		if ( flex_region_->length() < init_aa_.size() ) {
			num_aa_each_side_to_replace = flex_region_->length() / 2;
			if ( flex_region_->length() % 2 == 0 ) oddlength = false;
			else oddlength = true;
		}
		std::set< Size > res_to_repack;
		for ( Size i = 0; i < num_aa_each_side_to_replace; ++i ) {
			pose.replace_residue( flex_region_->start() + i, core::conformation::Residue( *(init_aa_[i+1].lock()), true), true);

			pose.replace_residue( flex_region_->stop() - i, core::conformation::Residue( *(init_aa_[init_aa_.size()-i].lock()), true), true);

			res_to_repack.insert( flex_region_->start() + i ); res_to_repack.insert( flex_region_->stop() - i );
		}
		if ( oddlength ) {
			res_to_repack.insert( flex_region_->start() + num_aa_each_side_to_replace );
			pose.replace_residue( flex_region_->start() + num_aa_each_side_to_replace, core::conformation::Residue( *(init_aa_[num_aa_each_side_to_replace+1].lock()), true), true);
		}
		core::pack::task::PackerTaskOP task( new core::pack::task::PackerTask_( pose ) );
		task->initialize_from_command_line();
		for ( Size i = 1; i <= task->total_residue(); ++i ) {
			if ( res_to_repack.find( i ) != res_to_repack.end() ) task->nonconst_residue_task( i ).restrict_to_repacking();
			else task->nonconst_residue_task( i ).prevent_repacking();
		}
		protocols::simple_moves::PackRotamersMoverOP packrot( new protocols::simple_moves::PackRotamersMover( enz_prot_->get_scorefxn(), task ) );
		packrot->apply( pose );
		(*enz_prot_->reduced_scorefxn())( pose );
	} // if(keep_existing_aa_identities)

	enz_prot_->add_pregenerated_enzdes_constraints( pose );

	return vlb.get_last_move_status() == MS_SUCCESS;


} //remodel_pose

utility::vector1< core::Size >
EnzdesRemodelMover::get_flex_region( ) const
{
	utility::vector1< core::Size> flex_pos;
	for ( core::Size i( flex_region_->start() );
			i<=flex_region_->stop(); ++i ) {
		flex_pos.push_back(i);
	}
	return flex_pos;
}


void
EnzdesRemodelMover::set_target_inverse_rotamers(utility::vector1< std::list < core::conformation::ResidueCOP > > & inv_rot)
{
	target_inverse_rotamers_.clear();
	target_inverse_rotamers_=inv_rot;
	tr << "New target inverse rotamers have been set via \"set_target_inverse_rotamers\"." << std::endl;
}


void
EnzdesRemodelMover::set_max_allowed_score_increase( core::Real sc_increase )
{
	max_allowed_score_increase_ = sc_increase;
}


bool
EnzdesRemodelMover::refine_pose(
	core::pose::Pose & pose
)
{

	//pose.dump_pdb("pose_bef_fa_refine.pdb");
	core::pose::Pose save_pose = pose;
	for ( auto & rcg : rcgs_ ) {

		rcg->add_remodel_constraints_to_pose( pose );

	}
	(*(enz_prot_->reduced_scorefxn()))(pose);
	//std::cerr << " done creating fa remcsts." << std::endl;
	tr << "backbone_stub constraint score before refine is " << pose.energies().total_energies()[ core::scoring::backbone_stub_constraint ] << std::endl;
	std::set< core::Size > dummy_repack_pos;
	if ( !enz_prot_->minimize_flexible_region( pose, flex_region_->index(), enz_prot_->reduced_scorefxn(), dummy_repack_pos, false, 0.01  ) ) {
		pose = save_pose;
		//pose.dump_pdb("pose_aftr_aftmin.pdb");
		return false;
	}
	//pose.dump_pdb("pose_aftr_aftmin.pdb");
	//std::cerr << " done refine minimizaing." << std::endl;
	enz_prot_->generate_ensemble_for_region( pose, flex_region_->index() );
	//pose.dump_pdb("pose_aftr_aftensemble.pdb");
	//std::cerr << " done creating refine ensemble." << std::endl;

	core::fragment::apply_best_scoring_fragdata( pose, *flex_region_, *(enz_prot_->reduced_scorefxn()) );

	//std::cerr << " done applying best frag." << std::endl;
	(*(enz_prot_->reduced_scorefxn()))(pose);
	//std::cerr << " done creating fa remcsts." << std::endl;
	tr << "backbone_stub constraint score after refine is " << pose.energies().total_energies()[ core::scoring::backbone_stub_constraint ] << std::endl;

	for ( auto & rcg : rcgs_ ) {
		rcg->remove_remodel_constraints_from_pose( pose );
	}

	//std::cerr << " done removing fa remcsts." << std::endl;
	(*(enz_prot_->reduced_scorefxn()))(pose);
	//pose.dump_pdb("pose_aft_fa_refine.pdb");
	return true;
} //refine pose


/// @detail  NOTE: the pose is reduced to poly a representation in this function
/// @details further this function sets up the following filters:
/// @details 1. PackerNeighborGraphFilter
/// @details 2. energy filter for reduced scorefunction not much higher than native conf
void
EnzdesRemodelMover::examine_initial_conformation(
	core::pose::Pose & pose
)
{
	predesign_filters_->clear();

	setup_postdesign_filters( pose );

	protocols::filters::FilterCOP png_filter = setup_packer_neighbor_graph_filter( pose );

	//create inverse rotamers if necessary
	if ( remodel_secmatch_ ) {
		create_target_inverse_rotamers( pose );
	}

	if ( keep_existing_aa_identities_ ) {
		init_aa_.clear();
		for ( Size i = flex_region_->start(); i <= flex_region_->stop(); ++i ) init_aa_.push_back( pose.residue_type( i ).get_self_weak_ptr() );
	}


	//reduce the pose to poly ala
	utility::vector1< core::Size > ala_pos;
	for ( core::Size i = 1; i <= pose.size(); ++i ) {

		if ( !pose.residue(i).is_protein() ) continue;

		if ( flex_region_->contains_seqpos( i ) && ( pose.residue_type( i ).aa() != core::chemical::aa_gly )  ) {
			ala_pos.push_back( i );
		} else if ( orig_task_->design_residue( i ) && ( pose.residue_type( i ).aa() != core::chemical::aa_gly  ) &&
				( pose.residue_type( i ).aa() != core::chemical::aa_pro ) ) {
			ala_pos.push_back( i );
		}
	}

	protocols::toolbox::pose_manipulation::construct_poly_ala_pose( pose, ala_pos, false, false, true );

	//score pose in reduced scorefxn to set up score filter correctly
	(*enz_prot_->reduced_scorefxn())( pose );

	protocols::simple_filters::ScoreCutoffFilterOP score_filter( new protocols::simple_filters::ScoreCutoffFilter() );


	//setting the predesign score filter to a pretty generous cutoff.
	//we only want this to prevent ridicoulous clashes to be rejected at this stage,
	//structures that are slightly worse can be improved in design
	core::Real predes_cutoff = pose.energies().total_energies()[ core::scoring::total_score ] + ( 10.0 * flex_region_->remodel_max_length() );
	if ( max_allowed_score_increase_ > predes_cutoff ) predes_cutoff = max_allowed_score_increase_;
	score_filter->set_cutoff( predes_cutoff );
	tr << "Cutoff for remodel score filter set to " << score_filter->cutoff() << "." << std::endl;

	//predesign_filters_.push_back( png_filter  );
	predesign_filters_->add_filter( score_filter );

	//last but not least, we have to initialise the start_to_cur sequence mapping
	start_to_current_smap_ = core::id::SequenceMappingOP( new core::id::SequenceMapping() );
	start_to_current_smap_->clear();

	for ( core::Size i = 1; i <= pose.size(); ++i ) start_to_current_smap_->push_back( i );

} //examine initial conformation


/// @details PackerNeighborGraphFilter, somewhat unusual/complicated:
/// @details every generated conformation needs to have a similar contact profile in a packer neighbor graph
/// @details as the original conformation in the nonlocal contacts graph, meaning:
/// @details a. there has to be the same number of contacts between the loop region and other designable regions
/// @details b. the packer neighbor graph of the loop region needs to have edges to undesignable buried nodes that
/// @details    have edges from the loop region in the nonlocal contacts graph of the original conformation
protocols::filters::FilterCOP
EnzdesRemodelMover::setup_packer_neighbor_graph_filter( core::pose::Pose const & pose) const
{

	using namespace utility::graph;

	protocols::toolbox::pose_metric_calculators::NonlocalContactsCalculator nl_calc( remodel_positions_, other_design_positions_ );
	basic::MetricValue< core::Size > mval_size;
	basic::MetricValue< utility::graph::GraphOP > mval_graph;

	nl_calc.get( (std::string) "region1_region2_nlcontacts_", mval_size, pose );
	nl_calc.get( (std::string) "nlcontacts_graph", mval_graph, pose );


	protocols::simple_filters::PackerNeighborGraphFilterOP png_filter( new protocols::simple_filters::PackerNeighborGraphFilter( orig_task_, enz_prot_->get_scorefxn() ) );

	png_filter->add_required_connections_between_regions( remodel_positions_, other_design_positions_, mval_size.value() );


	utility::graph::GraphCOP nl_graph = mval_graph.value();

	//now iterate through graph to get the residue required connections
	for ( Node::EdgeListConstIter edge_it( nl_graph->const_edge_list_begin() ), edge_end( nl_graph->const_edge_list_end() );
			edge_it != edge_end; ++edge_it ) {

		core::Size res1 = (*edge_it)->get_first_node_ind();
		core::Size res2 = (*edge_it)->get_second_node_ind();

		if ( remodel_positions_.find( res1 ) != remodel_positions_.end() ) {
			png_filter->add_required_connection_for_residue( res2 );
		} else if ( remodel_positions_.find( res2 ) != remodel_positions_.end() ) {
			png_filter->add_required_connection_for_residue( res1 );
		}

	}

	return png_filter;

}


std::string
EnzdesRemodelMover::generate_secstruct_string( core::pose::Pose & pose ) const {

	//A. if secondary structure strings in this mover have been set, we'll use one of these
	if ( user_provided_ss_.size() != 0 ) {
		std::string user_ss( user_provided_ss_[ numeric::random::rg().random_range(1, user_provided_ss_.size() ) ] );
		tr << "picking user-provided ss_string " << user_ss << "... ";
		return user_ss;
	}

	//B. if the user has specified secondary structure strings, let's use one of those
	if ( basic::options::option[ basic::options::OptionKeys::enzdes::enz_loops_file ].user()  &&
			flex_region_->enz_loop_info()->ss_strings_specified() ) {

		utility::vector1< std::string > const & possible_ss( flex_region_->enz_loop_info()->ss_strings() );
		std::string new_ss( possible_ss[ numeric::random::rg().random_range(1, possible_ss.size() ) ] );
		tr << "New secstruct (length " << new_ss.size() << ") for region " << region_to_remodel_ << " is " << new_ss << std::endl;
		return new_ss;
	}

	core::scoring::dssp::Dssp pose_ss( pose );
	pose_ss.insert_ss_into_pose( pose );

	//C. otherwise, we'll try to generate a secondary structure string of the new length
	// that matches the original string as closely as possible
	std::string orig_secstruct;

	for ( core::Size i = flex_region_->start(); i <= flex_region_->stop(); ++i ) {
		orig_secstruct.push_back( pose.secstruct( i ) );
	}

	core::Size low_length( flex_region_->remodel_min_length()), high_length( flex_region_->remodel_max_length());

	//std::cerr << "rm min length is " << low_length << ", rm max length is " << high_length << std::endl;

	core::Size new_length = numeric::random::rg().random_range( low_length, high_length );
	core::Size old_length = flex_region_->length();

	std::cout << "old length is " << old_length << ", new length is " << new_length << std::endl;

	//strategy to generate the new secstruct string:
	//1. create a fully degenerate string
	//2. set it to the values that are in the original string,
	//   this might leave some 'holes' if the new string is longer,
	//   or some not nice breaks if the new string is shorter
	//3. pass through again to try to fill the holes and introduce some randomness

	std::string new_secstruct;

	//1
	for ( core::Size i = 1; i <= new_length; ++i ) new_secstruct.push_back('D');

	//2a
	core::Size middle( new_length > old_length ? ( old_length / 2) : (new_length / 2 ) );

	for ( core::Size i = 0; i <= middle; ++i ) {

		new_secstruct[i] = orig_secstruct[i];

		new_secstruct[ new_length - i] = orig_secstruct[ old_length - i];

	}

	//2b. if the new length is bigger than the old length, we have to fill the hole
	if ( new_length > old_length ) {

		core::Size len_diff = (new_length - old_length);
		core::Size half_len_diff = (len_diff / 2) + 1;
		core::Size middle_plus_len_diff = middle + len_diff;

		for ( core::Size i = 1; i<= half_len_diff; ++ i ) {

			new_secstruct[ middle + i ] = new_secstruct[ middle ];
			new_secstruct[ middle_plus_len_diff - i ] = new_secstruct[ middle_plus_len_diff ];
		}

	}
	//std::cout << "new secstruct at second is " << new_secstruct << std::endl;

	//3
	core::Real prob = numeric::random::rg().uniform();

	if ( prob > ss_similarity_probability_ ) new_secstruct[0] = 'L';


	core::Size temp_stop( new_length - 1 );
	core::Real half_ss_sim( ss_similarity_probability_ / 2 );

	for ( core::Size i = 1; i < temp_stop; ++i ) {

		prob = numeric::random::rg().uniform();
		//std::cout << "step " << i << ", prob is " << prob << ",   ";

		if ( prob > ss_similarity_probability_ ) {

			//std::cout << "wtf, should set new_secstruct " << i << " to D, before change is " << new_secstruct[i];
			new_secstruct[ i ] = 'D';
			//std::cout << ", after change is " << new_secstruct[ i ] << std::endl;
		} else if ( prob > half_ss_sim ) new_secstruct[ i ] = new_secstruct[ i + 1 ];

		else new_secstruct[ i ] = new_secstruct[ i - 1 ];

	}

	prob = numeric::random::rg().uniform();

	if ( prob > ss_similarity_probability_ ) new_secstruct[ new_length - 1] = 'L';


	tr << "orig secstruct is " << orig_secstruct  << ", new secstruct is " << new_secstruct << std::endl;

	//done
	return new_secstruct;

}


void
EnzdesRemodelMover::apply_random_lowE_ligconf(
	core::pose::Pose & pose ) const
{
	core::scoring::ScoreFunctionCOP scofx(enz_prot_->reduced_scorefxn() );
	utility::vector1< core::Size > all_ligands( protocols::ligand_docking::get_ligand_seqpos( pose ) );

	for ( core::Size lig_num = 1; lig_num <= all_ligands.size(); ++lig_num ) {
		core::conformation::ResidueCOP cur_lig( pose.residue( all_ligands[lig_num] ).get_self_ptr() );
		utility::vector1< core::conformation::ResidueCOP > ligrots;
		enz_prot_->get_non_bb_clashing_rotamers( pose, all_ligands[lig_num], scofx, ligrots );
		ligrots.push_back( cur_lig ); //save current lig at the last position
		core::Size num_rots( ligrots.size() );
		utility::vector1< core::Real > remregion_scores( num_rots );

		for ( core::Size i =1; i <= num_rots; ++i ) {
			core::scoring::EnergyMap emap;
			for ( core::Size res(flex_region_->start()), resEnd( flex_region_->stop() ); res != resEnd; ++res ) {
				scofx->bump_check_backbone( *(ligrots[ i ]), pose.residue( res ), pose, emap );
			}
			remregion_scores[ i ] = scofx->weights().dot( emap );
		}// i loop over all conformers for this ligand

		utility::vector1< core::Size > lowE_ligs;
		for ( core::Size sco = 1; sco < num_rots; ++sco ) { //don't do last
			if ( remregion_scores[ sco ] <= remregion_scores[ num_rots ] ) lowE_ligs.push_back( sco );
		}
		lowE_ligs.push_back( num_rots );

		//now put a random one of the low scoring ligands into the pose
		core::Size ranconf( lowE_ligs[ numeric::random::rg().random_range(1, lowE_ligs.size() )] );
		tr << "Putting a ligand with remodel region score " << remregion_scores[ ranconf ] << " (orig ligand had " << remregion_scores[ num_rots ] << ") into pose." << std::endl;
		pose.replace_residue( all_ligands[lig_num ], *ligrots[ranconf], false );
	} //loop over all ligands in the pose
}


void
EnzdesRemodelMover::translate_atomnames_to_restype_set_atomids(
	core::pose::Pose const & pose,
	core::chemical::ResidueTypeSetCAP restype_set_cap,
	core::Size const seqpos,
	utility::vector1< std::string > const & atom_names,
	utility::vector1< core::Size > & atom_ids )
{
	core::chemical::ResidueTypeSetCOP restype_set( restype_set_cap );
	atom_ids.clear();

	for ( auto cur_at_name : atom_names ) {

		//lil tricky: the pose will be in centroid when the constraints are being generated,
		//so we have to make sure to put in the right atom ids
		std::string cur_res_name = pose.residue_type( seqpos ).name();

		//core::Size centroid_id = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID )->name_map( cur_res_name ).atom_index( cur_at_name );

		core::Size cur_at_id = restype_set->name_map( cur_res_name ).atom_index( cur_at_name );
		//at_ids.push_back( pose.residue_type( pos_it->first ).atom_index( *at_it ) );
		atom_ids.push_back( cur_at_id );

		//debug
		//std::cerr << "processing file info for cst pos " << seqpos << ", atom read was " << *at_it << ", translates to id " << pose.residue_type( seqpos ).atom_index( *at_it ) << " and this fa id translates to " << pose.residue_type(seqpos ).atom_name( pose.residue_type( seqpos ).atom_index( *at_it ) ) << ", translates to centroid id " << centroid_id << std::endl;
	}

}

protocols::forge::remodel::ResidueVicinityInfoOP
EnzdesRemodelMover::translate_res_interactions_to_rvinfos(
	core::pose::Pose const & pose,
	core::Size const targ_pos,
	core::Size const example_loop_seqpos,
	core::chemical::ResidueTypeSetCAP restype_set_cap,
	protocols::toolbox::match_enzdes_util::ResInteractions const & res_ints
)
{
	core::chemical::ResidueTypeSetCOP restype_set( restype_set_cap );
	core::Real dummy_val(0.0);
	utility::vector1< core::Size > targ_atom_ids, loopres_atom_ids;
	translate_atomnames_to_restype_set_atomids( pose, restype_set, targ_pos, res_ints.targ_atom_names(), targ_atom_ids );
	translate_atomnames_to_restype_set_atomids( pose, restype_set, example_loop_seqpos, res_ints.loopres_atom_names(), loopres_atom_ids );

	utility::vector1< core::Size > loopres_base_atom_ids, loopres_base2_atom_ids;
	utility::vector1< core::Size > targ_base_atom_ids, targ_base2_atom_ids;

	if ( res_ints.loopres_base_atom_names().size() != 0 ) translate_atomnames_to_restype_set_atomids( pose, restype_set, example_loop_seqpos, res_ints.loopres_base_atom_names(), loopres_base_atom_ids );
	if ( res_ints.loopres_base2_atom_names().size() != 0 ) translate_atomnames_to_restype_set_atomids( pose, restype_set, example_loop_seqpos, res_ints.loopres_base2_atom_names(), loopres_base2_atom_ids );
	if ( res_ints.targ_base_atom_names().size() != 0 ) translate_atomnames_to_restype_set_atomids( pose, restype_set, targ_pos, res_ints.targ_base_atom_names(), targ_base_atom_ids );
	if ( res_ints.targ_base2_atom_names().size() != 0 ) translate_atomnames_to_restype_set_atomids( pose, restype_set, targ_pos, res_ints.targ_base2_atom_names(), targ_base2_atom_ids );

	protocols::forge::remodel::ResidueVicinityInfoOP to_return( new protocols::forge::remodel::ResidueVicinityInfo( targ_pos, targ_atom_ids, loopres_atom_ids, res_ints.num_interactions() ) );

	//have to translate distance information in resints into proper function explicityly
	core::Real min_dis = std::max(0.0, res_ints.dis()->ideal_val() - res_ints.dis()->tolerance() );
	core::Real max_dis = res_ints.dis()->ideal_val() + res_ints.dis()->tolerance();
	core::Real force_k_dis = res_ints.dis()->force_const();
	to_return->set_dis( core::scoring::func::FuncOP( new core::scoring::constraints::BoundFunc(
		min_dis, max_dis, sqrt(1/ force_k_dis), "dis") ) );

	if ( (loopres_base_atom_ids.size() != 0 ) && res_ints.loop_ang() ) {
		runtime_assert( loopres_base_atom_ids.size() == loopres_atom_ids.size() );
		to_return->set_loopres_base_atoms( loopres_base_atom_ids );
		to_return->set_loop_ang( protocols::toolbox::match_enzdes_util::EnzConstraintParameters::convert_GeomSampleInfo_to_FuncOP( res_ints.loop_ang(), dummy_val ) );
	}

	if ( (targ_base_atom_ids.size() != 0 ) && res_ints.targ_ang() ) {
		runtime_assert( targ_base_atom_ids.size() == targ_atom_ids.size() );
		to_return->set_residue_base_atoms( targ_base_atom_ids );
		to_return->set_targ_ang( protocols::toolbox::match_enzdes_util::EnzConstraintParameters::convert_GeomSampleInfo_to_FuncOP( res_ints.targ_ang(), dummy_val ) );
	}

	if ( (loopres_base2_atom_ids.size() != 0 ) && res_ints.loop_dih() ) {
		runtime_assert( loopres_base2_atom_ids.size() == loopres_atom_ids.size() );
		runtime_assert( loopres_base_atom_ids.size() == loopres_base2_atom_ids.size() );
		to_return->set_loopres_base2_atoms( loopres_base2_atom_ids );
		to_return->set_loopres_base_atoms( loopres_base_atom_ids );
		to_return->set_loop_dih( protocols::toolbox::match_enzdes_util::EnzConstraintParameters::convert_GeomSampleInfo_to_FuncOP( res_ints.loop_dih(), dummy_val ) );
	}

	if ( (targ_base2_atom_ids.size() != 0 ) && res_ints.targ_dih() ) {
		runtime_assert( targ_base2_atom_ids.size() == targ_atom_ids.size() );
		runtime_assert( targ_base_atom_ids.size() == targ_base2_atom_ids.size() );
		to_return->set_residue_base2_atoms( targ_base2_atom_ids );
		to_return->set_targ_dih( protocols::toolbox::match_enzdes_util::EnzConstraintParameters::convert_GeomSampleInfo_to_FuncOP( res_ints.targ_dih(), dummy_val ) );
	}

	if ( (targ_base_atom_ids.size() != 0 ) && res_ints.lt_dih() ) {
		runtime_assert( targ_base_atom_ids.size() == targ_atom_ids.size() );
		runtime_assert( loopres_base_atom_ids.size() == loopres_atom_ids.size() );
		to_return->set_residue_base_atoms( targ_base_atom_ids );
		to_return->set_loopres_base_atoms( loopres_base_atom_ids );
		to_return->set_lt_dih( protocols::toolbox::match_enzdes_util::EnzConstraintParameters::convert_GeomSampleInfo_to_FuncOP( res_ints.lt_dih(), dummy_val ) );
	}

	return to_return;
}

/// @details putting a LengthEventCollector into the pose if not yet present
void
EnzdesRemodelMover::setup_cached_observers(
	core::pose::Pose & pose
){

	if ( !pose.observer_cache().has( core::pose::datacache::CacheableObserverType::LENGTH_EVENT_COLLECTOR ) ) {
		core::pose::datacache::LengthEventCollectorOP lencollect( new core::pose::datacache::LengthEventCollector() );

		pose.observer_cache().set( core::pose::datacache::CacheableObserverType::LENGTH_EVENT_COLLECTOR, lencollect );
	}
}

void
EnzdesRemodelMover::remove_cached_observers(
	core::pose::Pose & pose
){

	//nothing yet, the length event collector needs to stay because
	//stuff downstream needs it

	//fixing foldtree that somehow gets fucked up
	//by vlb stuff if we have several chains
	core::kinematics::FoldTree const & old_fold_tree( this->get_native_pose()->fold_tree() );
	core::kinematics::FoldTree new_fold_tree;
	for ( auto const & e : old_fold_tree ) {
		core::Size estart = (*start_to_current_smap_)[e.start() ];
		core::Size estop = (*start_to_current_smap_)[e.stop() ];
		new_fold_tree.add_edge( estart, estop, e.label() );
	}
	// Debugging information, as this foldtree reset has had issues in the past
	tr.Debug << "Resetting pose FoldTree." << std::endl;
	tr.Debug << "* Current Pose of size " << pose.size() << ": " << pose.annotated_sequence() << std::endl;
	tr.Debug << "Native Pose of size " << this->get_native_pose()->size() << ": " << this->get_native_pose()->annotated_sequence() << std::endl;
	tr.Debug << "Current FoldTree, size " << pose.fold_tree().nres() << ": " << pose.fold_tree() << std::endl;
	tr.Debug << "Starting (native) FoldTree, size " << old_fold_tree.nres() << ": " << old_fold_tree << std::endl;
	tr.Debug << "* Remodeled (new) FoldTree, size " << new_fold_tree.nres() << ": " << new_fold_tree << std::endl;
	pose.fold_tree( new_fold_tree );
}

void
EnzdesRemodelMover::process_length_change(
	core::pose::Pose & pose,
	core::id::SequenceMappingCOP smap
){
	enz_prot_->remap_resid( pose, *smap );
	core::id::combine_sequence_mappings( *start_to_current_smap_, *smap );
	protocols::enzdes::enzutil::create_remark_headers_from_cstcache( pose );

	for ( auto & rcg : rcgs_ ) {
		rcg->set_seqmap( smap );
	}

	for ( core::Size & non_remodel_match_po : non_remodel_match_pos_ ) {
		non_remodel_match_po = (*smap)[non_remodel_match_po];
	}
}


void
EnzdesRemodelMover::setup_postdesign_filters(
	core::pose::Pose const & pose
){

	if ( pose.size() > 1 ) {} //hack to prevent compiler warning for now

	postdesign_filters_->clear();


} //setup_postdesign_filters


bool
EnzdesRemodelMover::secmatch_after_remodel(
	core::pose::Pose & pose
)
{

	//if this option is off, return right away
	if ( ! remodel_secmatch_ ) return true;

	//if no parameters are missing, return right away
	protocols::toolbox::match_enzdes_util::EnzConstraintIOCOP cstio( protocols::enzdes::enzutil::get_enzcst_io( pose ) );
	if ( (!cstio) || cstio->enz_cst_params_missing_in_pose( pose ).size() == 0 ) return true;

	tr << "There are catalytic interactions missing in the pose, starting secondary matcher..." << std::endl;
	cstio->remove_constraints_from_pose( pose, false, true );

	core::conformation::ResidueCOP ligres;
	for ( core::Size i = 1; i <= cstio->num_mcfi_lists(); ++i ) {
		protocols::toolbox::match_enzdes_util::EnzConstraintParametersCOP cst_params( cstio->enz_cst_params( i ) );

		if ( cst_params->missing_in_pose( pose ) ) {
			core::Size present_template( cst_params->get_missing_template_other_res( pose )->param_index() );
			core::Size present_seqpos( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( i )->template_res_cache( present_template )->seqpos_map_begin()->first );
			if ( !ligres ) {
				ligres = pose.residue( present_seqpos ).get_self_ptr();
				break;
			}
		}
	}
	if ( !ligres || !(ligres->type().is_ligand()) ) ligres = core::conformation::ResidueCOP( core::conformation::ResidueOP( new core::conformation::Residue( *(cstio->mcfi_list( 1 )->mcfi( 1 )->allowed_restypes( 1 )[1]), true ) ) );

	protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->set_cst_cache( nullptr ); //wipe this out for now, matcher will overwrite
	protocols::match::MatcherMoverOP matcher_mover( new protocols::match::MatcherMover() );

	//generate the positions
	utility::vector1< core::Size > match_positions;
	for ( core::Size i = flex_region_->start(); i <= flex_region_->stop(); ++i ) match_positions.push_back( i );
	for ( utility::vector1< core::Size >::const_iterator other_match_pos_it = non_remodel_match_pos_.begin(); other_match_pos_it != non_remodel_match_pos_.end(); ++other_match_pos_it ) match_positions.push_back( *other_match_pos_it );
	matcher_mover->set_return_single_random_match( true );
	matcher_mover->set_match_positions( match_positions );
	matcher_mover->set_ligres( ligres );
	matcher_mover->apply( pose );

	protocols::enzdes::AddOrRemoveMatchCsts cstmover;
	cstmover.set_cst_action( protocols::enzdes::ADD_NEW );
	if ( matcher_mover->get_last_move_status() == protocols::moves::MS_SUCCESS ) {
		//cstio->add_constraints_to_pose( pose, enz_prot_->get_scorefxn(), false );
		cstmover.apply( pose );
		return true;
	}
	cstmover.set_accept_blocks_missing_header(true );
	cstmover.apply( pose );
	return false;

}  //secmatch_after_remodel


/// @brief setup up remodel constraints according to the info specified
/// in the enzdes loops file
void
EnzdesRemodelMover::setup_rcgs(
	core::pose::Pose & pose,
	protocols::forge::components::VarLengthBuild & vlb )
{

	vlb.clear_rcgs();
	rcgs_.clear();

	//rcgs from inverse rotamers
	if ( remodel_secmatch_ ) {
		setup_rcgs_from_inverse_rotamers( vlb );
	}

	//this function only does stuff if a loop file was read in,
	//bc otherwise no rcg info specified
	if ( ! basic::options::option[ basic::options::OptionKeys::enzdes::enz_loops_file ].user() ) return;

	using namespace protocols::enzdes;
	using namespace protocols::toolbox::match_enzdes_util;

	core::chemical::ResidueTypeSetCAP centroid_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );

	core::chemical::ResidueTypeSetCAP fa_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	EnzdesLoopInfoCOP loopinf = flex_region_->enz_loop_info();

	utility::vector1< protocols::forge::remodel::ResidueVicinityInfoOP > cen_rv_infos;

	utility::vector1< protocols::forge::remodel::ResidueVicinityInfoOP > fa_rv_infos;

	//loop over all the cst targets in this loop info, i.e. all cst target records in the loop file
	for ( auto const & cst_targ : loopinf->cst_interactions() ) {

		protocols::toolbox::match_enzdes_util::EnzCstTemplateResCOP targ_template;
		protocols::toolbox::match_enzdes_util::EnzCstTemplateResCacheCOP targ_template_cache;
		protocols::toolbox::match_enzdes_util::EnzConstraintIOCOP cstio( enzutil::get_enzcst_io( pose ) );
		if ( !cstio ) break;

		//which residue of the desired block are we interested in?
		if ( cst_targ.resA() == true ) {
			targ_template = cstio->enz_cst_params( cst_targ.cst_block() )->resA();
			targ_template_cache = protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( cst_targ.cst_block() )->template_res_cache( 1 );
		} else {
			targ_template = cstio->enz_cst_params( cst_targ.cst_block() )->resB();
			targ_template_cache = protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( cst_targ.cst_block() )->template_res_cache( 1 );
		}

		//loop over all positions that the targ residue is at
		//note: usually this will be one
		//old
		//for( std::map< Size, EnzCstTemplateResAtomsOP >::const_iterator pos_it = targ_template->respos_map_begin();
		//  pos_it != targ_template->respos_map_end(); ++pos_it ){
		for ( auto pos_it = targ_template_cache->seqpos_map_begin();
				pos_it != targ_template_cache->seqpos_map_end(); ++pos_it ) {

			cen_rv_infos.push_back( translate_res_interactions_to_rvinfos( pose, pos_it->first, flex_region_->start(), centroid_set, cst_targ ) );
			fa_rv_infos.push_back( translate_res_interactions_to_rvinfos( pose, pos_it->first, flex_region_->start(), fa_set, cst_targ ) );

		} // loop over all positions of cst_target
	} //loop over all cst_targets


	//now loop over all the other desired interactions
	for ( auto const & res_int : loopinf->res_interactions() ) {

		//target residue of this rcgs: we have to correct for eventual length changes that have happened already
		core::Size targ_res = (*start_to_current_smap_)[ res_int.targ_res() ];

		cen_rv_infos.push_back( translate_res_interactions_to_rvinfos( pose, targ_res, flex_region_->start(), centroid_set, res_int ) );
		fa_rv_infos.push_back( translate_res_interactions_to_rvinfos( pose, targ_res, flex_region_->start(), fa_set, res_int ) );

	}

	if ( cen_rv_infos.size() != 0 ) {
		using namespace protocols::forge::remodel;
		vlb.add_rcg( RemodelConstraintGeneratorOP( new ResidueVicinityRCG( flex_region_->start(), flex_region_->stop(), cen_rv_infos ) ) );
		rcgs_.push_back( protocols::forge::remodel::RemodelConstraintGeneratorOP( new ResidueVicinityRCG( flex_region_->start(), flex_region_->stop(), fa_rv_infos ) ) );
	}
}


void
EnzdesRemodelMover::setup_rcgs_from_inverse_rotamers(
	protocols::forge::components::VarLengthBuild & vlb
)
{
	for ( core::Size i =1; i <= target_inverse_rotamers_.size(); ++i ) {
		using namespace protocols::forge::constraints;
		using namespace protocols::forge::remodel;
		RemodelConstraintGeneratorOP newrcg( new InverseRotamersRCG( flex_region_->start(), flex_region_->stop(),target_inverse_rotamers_[i] ) );
		newrcg->set_id( "invrot1" );
		vlb.add_rcg( newrcg );

		newrcg = RemodelConstraintGeneratorOP( new InverseRotamersRCG( flex_region_->start(), flex_region_->stop(), target_inverse_rotamers_[i] ) );
		newrcg->set_id( "invrot2" );
		rcgs_.push_back( newrcg );
	}
}

void
EnzdesRemodelMover::create_target_inverse_rotamers(
	core::pose::Pose & pose )
{

	target_inverse_rotamers_.clear();
	non_remodel_match_pos_.clear();
	tr << "Beginning to build inverse rotamers...    " << std::endl;
	//protocols::toolbox::match_enzdes_util::EnzConstraintIOCOP enzcst_io (enz_prot_->cst_io() );
	protocols::toolbox::match_enzdes_util::EnzConstraintIOCOP enzcst_io( protocols::enzdes::enzutil::get_enzcst_io( pose ) );
	if ( !enzcst_io ) return;

	for ( core::Size i = enzcst_io->num_mcfi_lists(); i >= 1; --i ) {
		protocols::toolbox::match_enzdes_util::EnzConstraintParametersCOP cst_params( enzcst_io->enz_cst_params( i ) );
		bool invrots_build(false);

		//if the interaction isn't present yet, build inverse rotamers based on cst file
		if ( cst_params->missing_in_pose( pose ) ) {
			tr << "Building inverse rotamers for missing MatchConstraint " << i << "...   " << std::endl;
			core::Size param_index( cst_params->get_missing_template_other_res( pose )->param_index() );
			protocols::toolbox::match_enzdes_util::EnzCstTemplateResCacheCOP temp_res_cache( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( i )->template_res_cache( param_index ) );
			if ( temp_res_cache->seqpos_map_size() != 1 ) {
				utility_exit_with_message("When trying to build inverse rotamers in enzdes remodel, !=1 target residues were found to be saved in the cst cache");
			}
			core::Size targ_seqpos( temp_res_cache->seqpos_map_begin()->first );
			target_inverse_rotamers_.push_back( enzcst_io->mcfi_list( i )->inverse_rotamers_against_residue( param_index, pose.residue( targ_seqpos ).get_self_ptr() ) );
			invrots_build = true;

			//note: if the target residue is also in a region being remodelled,
			//we have to remove it's position from the cst cache
			if ( flex_region_->contains_seqpos( targ_seqpos ) ) {
				enzcst_io->remove_position_from_template_res_for_block( pose, targ_seqpos, i );
			}
		} else { //if(missing in pose )
			//otherwise we have to check whether a catalytic position happens to be in the flex_region
			//i.e. a preexisting catalytic interaction is being remodeled
			protocols::toolbox::match_enzdes_util::EnzdesCstParamCacheOP param_cache( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( i ) );
			runtime_assert( param_cache->template_res_cache_size() == 2 );

			for ( core::Size template_res = 1; template_res <= param_cache->template_res_cache_size(); ++template_res ) {

				for ( auto seqpos_it = param_cache->template_res_cache( template_res )->seqpos_map_begin(), seqpos_end = param_cache->template_res_cache( template_res )->seqpos_map_end();
						seqpos_it != seqpos_end; ) {

					if ( flex_region_->contains_seqpos( seqpos_it->first ) ) {
						core::Size seqpos( seqpos_it->first );

						core::Size other_res( template_res == 1 ? 2 : 1 );
						runtime_assert( param_cache->template_res_cache( other_res )->seqpos_map_size() == 1 );
						core::Size other_seqpos( param_cache->template_res_cache( other_res )->seqpos_map_begin()->first );

						if ( pose.residue_type( other_seqpos ).is_protein() ) non_remodel_match_pos_.push_back( other_seqpos );
						//complication: if this is a residue that plays a role in several cst blocks,
						//we only build inverse rotamers according to the geometry specified in the
						//earliest cst block that this residue appears in
						core::Size corresponding_res_block( template_res == 1 ? cst_params->resA()->corresponding_res_block() : cst_params->resB()->corresponding_res_block() );
						if ( (corresponding_res_block == 0 ) || ( corresponding_res_block > i ) ) {
							tr << "Catalytic residue for MatchConstraint " << i << " at position " << seqpos << " is in remodeled region, building inverse rotamers... " << std::endl;
							target_inverse_rotamers_.push_back( std::list<core::conformation::ResidueCOP> () );
							if ( include_existing_conf_as_invrot_target_ ) target_inverse_rotamers_[ target_inverse_rotamers_.size() ].push_back( core::conformation::ResidueCOP( core::conformation::ResidueOP( new core::conformation::Residue( pose.residue( seqpos ) ) ) ) );
							invrots_build = true;
							std::list< core::conformation::ResidueCOP > cur_inv_rots( enzcst_io->mcfi_list( i )->inverse_rotamers_against_residue( other_res, pose.residue( other_seqpos ).get_self_ptr() ) );
							if ( !cur_inv_rots.empty() /*size() != 0*/ ) target_inverse_rotamers_[ target_inverse_rotamers_.size() ].splice( target_inverse_rotamers_[ target_inverse_rotamers_.size() ].end(), cur_inv_rots );
						} else {
							tr << "Catalytic residue for MatchConstraint " << i << " at position " << seqpos << " is in remodeled region, but we're not building inverse rotamers because the same residue is also specified in block " << corresponding_res_block << " of the cstfile. Inverse rotamers will be built according to the geometry specified in that block." << std::endl;
						}
						//need to remove the position from the constraints,
						++seqpos_it; //crucial to avoid iterator becoming invalidated
						enzcst_io->remove_constraints_from_pose_for_block( pose, i, false );
						enzcst_io->remove_position_from_template_res_for_block( pose, seqpos, i );
						enzcst_io->clear_active_pose_constraints_for_block( pose, i );
						enz_prot_->design_targets( pose );
						protocols::enzdes::enzutil::remove_remark_header_for_geomcst( pose, i );
					} else ++seqpos_it;
				} //loop over seqpos of template res
			} //loop over template res
		} //if( missing in pose ){} else

		//if invrots were build in either of the cases, make the rcgs and possibly dump them
		if ( invrots_build ) {
			core::Size invrot_index( target_inverse_rotamers_.size() );
			tr << target_inverse_rotamers_[ invrot_index ].size() << " inverse rotamers were built for MatchConstraint " << i << "." << std::endl;

			if ( basic::options::option[basic::options::OptionKeys::enzdes::dump_inverse_rotamers] ) {
				core::Size modelcounter(1), atomcounter(1);
				std::string filename = "Inverse_rotamers_geomcst_"+utility::to_string( i )+".pdb";
				std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);

				for ( std::list< core::conformation::ResidueCOP >::const_iterator rot_it( target_inverse_rotamers_[ invrot_index ].begin() ), rot_end( target_inverse_rotamers_[ invrot_index ].end() ); rot_it != rot_end; ++rot_it ) {
					file << "MODEL  " << modelcounter << " \n";
					core::io::pdb::dump_pdb_residue( **rot_it, atomcounter, file );
					file << "ENDMDL   \n";
					modelcounter++;
					atomcounter = 1;
				}
				file.close();
			}
		} //if( invrots_build)
	} //loop over all mcfi lists
}


void
EnzdesRemodelMover::set_seq_mapping(
	core::id::SequenceMappingOP seq_map_in
) {
	start_to_current_smap_ = seq_map_in;
}

void
EnzdesRemodelMover::set_keep_existing_aa_identities(
	bool setting
)
{
	keep_existing_aa_identities_ = setting;
}

core::id::SequenceMappingCOP
EnzdesRemodelMover::get_seq_mapping() const {
	return start_to_current_smap_;
}

std::string EnzdesRemodelMover::get_name() const {
	return mover_name();
}

std::string EnzdesRemodelMover::mover_name() {
	return "EnzdesRemodelMover";
}

void EnzdesRemodelMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "remodel_secmatch", xsct_rosetta_bool, "XRW TO DO" )
		+ XMLSchemaAttribute( "reinstate_foldtree", xsct_rosetta_bool, "XRW TO DO" )
		+ XMLSchemaAttribute( "remodel_region", xsct_non_negative_integer, "XRW TO DO" )
		+ XMLSchemaAttribute( "include_existing_conf_as_invrot_target", xsct_rosetta_bool, "XRW TO DO" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string EnzdesRemodelMoverCreator::keyname() const {
	return EnzdesRemodelMover::mover_name();
}

protocols::moves::MoverOP
EnzdesRemodelMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new EnzdesRemodelMover );
}

void EnzdesRemodelMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	EnzdesRemodelMover::provide_xml_schema( xsd );
}



} //namespace enzdes
} //namespace protocols
