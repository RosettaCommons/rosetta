// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/FastDesign.cc
/// @brief The FastDesign
/// @detailed
/// @author Tom Linsky


//Unit Headers
#include <devel/denovo_design/FastDesign.hh>
#include <devel/denovo_design/FastDesignCreator.hh>

//Project Headers
//#include <protocols/denovo_design/RestrictWorstRegion.hh>
//#include <protocols/flxbb/filters/HighestEnergyRegion.hh>
//#include <protocols/flxbb/filters/DesignBySecondaryStructure.hh>
//#include <protocols/flxbb/filters/PsiPredFilter.hh>
#include <protocols/relax/Ramady.hh>
#include <protocols/relax/util.hh>

//Core Headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/util/SwitchResidueTypeSet.hh>

//Protocol Headers
#include <protocols/filters/Filter.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <utility/excn/Exceptions.hh>
#include <protocols/elscripts/util.hh>
#include <protocols/toolbox/task_operations/LimitAromaChi2Operation.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>


//Basic Headers
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

//Utility Headers
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/random/random.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

//C++ Headers
#include <fstream>



#ifdef GL_GRAPHICS
#include <protocols/viewer/viewers.hh>
#endif

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("devel.denovo_design.FastDesign");

using namespace core;
using namespace core::io::silent ;
using io::pdb::dump_pdb;
////////////////////////////////////////////////////////////////////////////////////////////////////

namespace devel {
namespace denovo_design {
////////////////////////////////////////////////////////////////////////////////////////////////////


using namespace ObjexxFCL;

std::string
FastDesignCreator::keyname() const
{
  return FastDesignCreator::mover_name();
}

protocols::moves::MoverOP
FastDesignCreator::create_mover() const {
  return new FastDesign();
}

std::string
FastDesignCreator::mover_name()
{
  return "FastDesign";
}

///  ---------------------------------------------------------------------------------
///  FastDesign main code:
///  ---------------------------------------------------------------------------------

/// @brief default constructor
FastDesign::FastDesign() :
	RelaxProtocolBase( "FastDesign" ),
	default_repeats_( basic::options::option[ basic::options::OptionKeys::relax::default_repeats ]() ),
  ramady_( basic::options::option[ basic::options::OptionKeys::relax::ramady ]() ),
	ramady_num_rebuild_( basic::options::option[ basic::options::OptionKeys::relax::ramady_max_rebuild ]() ),
	ramady_rms_limit_( basic::options::option[ basic::options::OptionKeys::relax::ramady_rms_limit ]() ),
	ramady_cutoff_( basic::options::option[ basic::options::OptionKeys::relax::ramady_cutoff ]() ),
	ramady_force_( basic::options::option[ basic::options::OptionKeys::relax::ramady_force ]() ),
	repack_( basic::options::option[ basic::options::OptionKeys::relax::chi_move]() ),
	test_cycles_( basic::options::option[ basic::options::OptionKeys::run::test_cycles ]() ),
	dumpall_( false ),
	script_max_accept_( basic::options::option[ basic::options::OptionKeys::relax::script_max_accept ]() ),
	symmetric_rmsd_( basic::options::option[ basic::options::OptionKeys::evaluation::symmetric_rmsd ]() ),
	allow_design_( true ),
	design_worst_( false ),
	design_by_psipred_( false ),
	design_by_frag_qual_( false ),
	only_design_changes_( false ),
	blueprint_file_( "" ),
	ramp_design_constraints_( false ),
	clear_designable_residues_( false ),
	rank_scorefxn_( NULL ),
	max_redesigns_( 5 ),
	resfile_( "" ),
	checkpoints_( "FastDesign" ),
	movemap_tag_( NULL ),
	dump_counter_( 0 ),
  //frag_qual_op_( NULL ),
	regions_to_design_( 1 ),
	run_count_( 0 ),
	cached_sequence_( "" )
  //worst_region_mover_( NULL )
{
	script_.clear();
	filters_.clear();
	allowed_aas_.clear();
	num_redesigns_.clear();
	read_script_file( "", default_repeats_ );
}

/// @brief copy constructor
FastDesign::FastDesign( FastDesign const & rval ) :
	RelaxProtocolBase( rval ),
	default_repeats_( rval.default_repeats_ ),
  ramady_( rval.ramady_ ),
	ramady_num_rebuild_( rval.ramady_num_rebuild_ ),
	ramady_rms_limit_( rval.ramady_rms_limit_ ),
	ramady_cutoff_( rval.ramady_cutoff_ ),
	ramady_force_( rval.ramady_force_ ),
	repack_( rval.repack_ ),
	test_cycles_( rval.test_cycles_ ),
	dumpall_( rval.dumpall_ ),
	script_max_accept_( rval.script_max_accept_ ),
	symmetric_rmsd_( rval.symmetric_rmsd_ ),
	allow_design_( rval.allow_design_ ),
	design_worst_( rval.design_worst_ ),
	design_by_psipred_( rval.design_by_psipred_ ),
	design_by_frag_qual_( rval.design_by_frag_qual_ ),
	only_design_changes_( rval.only_design_changes_ ),
	blueprint_file_( rval.blueprint_file_ ),
	ramp_design_constraints_( rval.ramp_design_constraints_ ),
	clear_designable_residues_( rval.clear_designable_residues_ ),
	rank_scorefxn_( rval.rank_scorefxn_ ),
	max_redesigns_( rval.max_redesigns_ ),
	resfile_( rval.resfile_ ),
	checkpoints_( rval.checkpoints_ ),
	script_( rval.script_ ),
	movemap_tag_( rval.movemap_tag_ ),
	dump_counter_( rval.dump_counter_ ),
	filters_( rval.filters_ ),
	allowed_aas_( rval.allowed_aas_ ),
	num_redesigns_( rval.num_redesigns_ ),
//frag_qual_op_( rval.frag_qual_op_ ),
	regions_to_design_( rval.regions_to_design_ ),
	run_count_( rval.run_count_ ),
	cached_sequence_( rval.cached_sequence_ )
//worst_region_mover_( rval.worst_region_mover_ )
{}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
FastDesign::~FastDesign() {}


/// Return a copy of ourselves
protocols::moves::MoverOP
FastDesign::clone() const {
	return new FastDesign(*this);
}

void
FastDesign::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose
) {
	set_scorefxn( protocols::rosetta_scripts::parse_score_function(tag, data)->clone() );
	core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
	mm->set_chi( true );
	mm->set_bb( true );
	mm->set_jump( true );

	set_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );

	// initially, all backbone torsions are movable
	protocols::rosetta_scripts::parse_movemap( tag, pose, mm, data, false);
	set_movemap(mm);

	default_repeats_ = tag->getOption< int >( "repeats", 8 );
	std::string script_file = tag->getOption< std::string >("relaxscript", "" );

	cartesian (tag->getOption< bool >( "cartesian", false ) );

	allow_design_ = tag->getOption< bool >( "allow_design", allow_design_ );
	design_worst_ = tag->getOption< bool >( "only_design_worst_region", design_worst_ );
	design_by_psipred_ = tag->getOption< bool >( "design_by_psipred", design_by_psipred_ );
	design_by_frag_qual_ = tag->getOption< bool >( "design_by_frag_qual", design_by_frag_qual_ );
	only_design_changes_ = tag->getOption< bool >( "design_changes", only_design_changes_ );
	if ( design_worst_ || design_by_psipred_ || design_by_frag_qual_ ) {
		allow_design_ = true;
	}
	blueprint_file_ = tag->getOption< std::string >( "blueprint", blueprint_file_ );
	// blueprint file must be specified if design_by_psipred is set
	if ( blueprint_file_ == "" && design_by_psipred_ ) {
		utility_exit_with_message( "Blueprint file must be specified if design_by_psipred is set." );
	}
	clear_designable_residues_ = tag->getOption< bool >( "clear_designable_residues", clear_designable_residues_ );
	ramp_design_constraints_ = tag->getOption< bool >( "ramp_design_constraints", ramp_design_constraints_ );
	dumpall_ = tag->getOption< bool >( "dumpall", dumpall_ );
	std::string const rank_scorefxn_name( tag->getOption< std::string >( "rank_scorefxn", "" ) );
	if ( rank_scorefxn_name == "" ) {
		rank_scorefxn_ = NULL;
	} else {
		rank_scorefxn_ = data.get< core::scoring::ScoreFunction * >( "scorefxns", rank_scorefxn_name )->clone();
	}

	max_redesigns_ = tag->getOption< core::Size >( "max_redesigns", max_redesigns_ );
	/*std::string const filters_str( tag->getOption< std::string >( "filters", "," ) );
	utility::vector1< std::string > const filter_list( utility::string_split( filters_str, ',' ) );
	for ( core::Size i=1; i<= filter_list.size(); ++i ) {
		if ( filter_list[i] != "" ) {
			TR << "Setting up filter " << filter_list[i] << std::endl;
			filters_.push_back( protocols::rosetta_scripts::parse_filter( filter_list[i], filters ) );
		}
		}*/
	/*
	std::string const design_worst_mover_str( tag->getOption< std::string >( "restrict_worst_mover", "" ) );
	if ( design_worst_mover_str != "" ) {
		moves::MoverOP mover = protocols::rosetta_scripts::parse_mover( design_worst_mover_str, movers );
		assert( utility::pointer::dynamic_pointer_cast< RestrictWorstRegion >( mover ) );
		worst_region_mover_ = utility::pointer::static_pointer_cast< RestrictWorstRegion >( mover );
		}*/

	utility::vector1< utility::tag::TagPtr > const & filter_tags( tag->getTags() );
	for ( core::Size i=1; i<=filter_tags.size(); ++i ) {
		utility::tag::TagPtr const layer_tag( filter_tags[i] );
		std::string op = layer_tag->getName(); // Should be "AND" just like in Generic MC
		std::string const filter_name( layer_tag->getOption< std::string >( "filter", "" ) );
		// tolerance means if the current value is within X of the starting value, it's OK
		// default = not used
		bool const use_tolerance( layer_tag->hasOption( "tolerance" ) );
		core::Real const tolerance( layer_tag->getOption< core::Real >( "tolerance", 0.0 ) );
		// threshold means if the current value is less than X, it's OK
		// default = not used
		bool const use_threshold( layer_tag->hasOption( "threshold" ) );
		core::Real const threshold( layer_tag->getOption< core::Real >( "threshold", 0.0 ) );
		// Tells if a high value is good or bad
		// default = low
		std::string const sample_type( layer_tag->getOption< std::string >( "sample_type", "low" ) );
		if ( filter_name != "" ) {
			// if the filter is named, add it to the list
			FilterParams fp;
			fp.filter = protocols::rosetta_scripts::parse_filter( filter_name, filters );
			if ( sample_type == "low" ) {
				fp.sample_low = true;
			} else if ( sample_type == "high" ) {
				fp.sample_low = false;
			} else {
				utility_exit_with_message( "Invalid sample type specified: " + sample_type );
			}
			fp.use_tolerance = use_tolerance;
			fp.tolerance = tolerance;
			fp.use_threshold = use_threshold;
			if ( fp.sample_low ) {
				fp.threshold = threshold;
			} else {
				fp.threshold = -threshold;
			}
			TR << " Read filter name=" << filter_name << " tolerance=" << tolerance << " threshold=" << threshold << " sample_type=" << sample_type << std::endl;
			filters_.push_back( fp );
		}
	}

	if( script_file == "" ) {
		read_script_file( "", default_repeats_ );
	} else {
		read_script_file( script_file );
	}
	initialize_aa_cache( pose );
}

void FastDesign::parse_def( utility::lua::LuaObject const & def,
	utility::lua::LuaObject const & score_fxns,
	utility::lua::LuaObject const & tasks,
	protocols::moves::MoverCacheSP ) {
	if( def["scorefxn"] ) {
		set_scorefxn( protocols::elscripts::parse_scoredef( def["scorefxn"], score_fxns ) );
	} else {
		set_scorefxn( score_fxns["score12"].to<core::scoring::ScoreFunctionSP>()->clone()  );
	}

	core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
	mm->set_chi( true );
	mm->set_bb( true );
	mm->set_jump( true );
	// initially, all backbone torsions are movable
	if( def["movemap"] )
		protocols::elscripts::parse_movemapdef( def["movemap"], mm );
	set_movemap(mm);

	if( def["tasks"] ) {
		core::pack::task::TaskFactoryOP new_task_factory( protocols::elscripts::parse_taskdef( def["tasks"], tasks ));
		if ( new_task_factory == 0) return;
		set_task_factory(new_task_factory);
	}

	default_repeats_ = def["repeats"] ? def["repeats"].to<int>() : 8;
	std::string script_file = def["relaxscript"] ? def["relaxscript"].to<std::string>() : "";

	cartesian (def["cartesian"] ? def["cartesian"].to<bool>() : false );

	if( script_file == "" ) {
		read_script_file( "", default_repeats_ );
	} else {
		read_script_file( script_file );
	}
}

/// @brief sets the number of repeats
void
FastDesign::set_repeats( core::Size const repeats )
{
	default_repeats_ = repeats;
}

// grab the score and remember the pose if the score is better then ever before.
void FastDesign::cmd_accept_to_best(
  const core::scoring::ScoreFunctionOP local_scorefxn,
  core::pose::Pose &pose,
  core::pose::Pose &best_pose,
  const core::pose::Pose &start_pose,
  core::Real       &best_score,
  core::Size       &accept_count
)
{
	using namespace core::scoring;
	using namespace core::conformation;
  core::Real score = (*local_scorefxn)( pose );
  if( ( score < best_score) || (accept_count == 0) ){
    best_score = score;
    best_pose = pose;
  }
  #ifdef BOINC_GRAPHICS
		using namespace protocols;
    boinc::Boinc::update_graphics_low_energy( best_pose, best_score  );
    boinc::Boinc::update_graphics_last_accepted( pose, score );
    //boinc::Boinc::update_mc_trial_info( total_count , "FastDesign" );  // total_count not defined
  #endif
  core::Real rms = 0, irms = 0;
  if ( core::pose::symmetry::is_symmetric( pose ) && symmetric_rmsd_ ) {
    rms = CA_rmsd_symmetric( *get_native_pose() , best_pose );
    irms = CA_rmsd_symmetric( start_pose , best_pose );
  } else {
    rms = native_CA_rmsd( *get_native_pose() , best_pose );
    irms = native_CA_rmsd( start_pose , best_pose );
  }
  TR << "MRP: " << accept_count << "  " << score << "  " << best_score << "  " << rms << "  " << irms << "  " << std::endl;
}


void FastDesign::do_minimize(
  core::pose::Pose &pose,
  core::Real tolerance,
  core::kinematics::MoveMapOP local_movemap,
	core::scoring::ScoreFunctionOP local_scorefxn
){
	using namespace core::scoring;
	using namespace core::conformation;

  protocols::simple_moves::MinMoverOP min_mover;
  if ( core::pose::symmetry::is_symmetric( pose ) )  {
    min_mover = new protocols::simple_moves::symmetry::SymMinMover( local_movemap, local_scorefxn, min_type_, tolerance, true );
  } else {
    min_mover = new protocols::simple_moves::MinMover( local_movemap, local_scorefxn, min_type_, tolerance, true );
  }
	min_mover->cartesian( cartesian_ );
  min_mover->apply( pose );
}


////////////////////////////////////////////////////////////////////////////////////////////////////
void FastDesign::apply( core::pose::Pose & pose ){
	using namespace core::scoring;
	using namespace core::conformation;
	using namespace core::pack;
	using namespace core::pack::task;
	using namespace core::kinematics;
	using namespace protocols;

	// increment run counter
	++run_count_;

	TR.Debug   << "================== FastDesign: " << script_.size() << " ===== Run: " << run_count_ << " ===================" << std::endl;
 	protocols::rosetta_scripts::parse_movemap( movemap_tag_, pose, get_movemap() );
	reset_status();

#if defined GL_GRAPHICS
    protocols::viewer::add_conformation_viewer( pose.conformation(), "TESTING");
#endif
	// One out of 10 times dont bother doing Ramady Relax. You may wonder why - the reason is that occasionally
	// there are phi-pso pairs outside of the commonly allowed ranges which are *correct* for some sort of quirk
	// of nature. In that case using ramady relax screws things up, so to allow for that chance only execture ramady relax
	// 90% of the time.
	bool do_rama_repair = ramady_;
	if( !ramady_force_ && numeric::random::uniform() <= 0.1 ) do_rama_repair = false;


	// Relax is a fullatom protocol so switch the residue set here. The rotamers
	// wont be properly packed at this stage but it doesnt really matter.
	// They'll get repacked shortly.
	if( !pose.is_fullatom() ){
		core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD);
	}

	// Make a local copy of movemap, called local movemap. The reason is that
	// the "constrain to coordinates" codes will wanna mess witht the movemap..
	core::kinematics::MoveMapOP local_movemap = get_movemap()->clone();
	initialize_movemap( pose, *local_movemap );
	relax::make_dna_rigid( pose, *local_movemap );
	set_movemap( local_movemap );

	// Deal with constraint options and add coodrinate constraints for all or parts of the structure.
	set_up_constraints( pose, *local_movemap );

	// Make sure we only allow symmetrical degrees of freedom to move and convert the local_movemap
	// to a local movemap
  if ( core::pose::symmetry::is_symmetric( pose )  )  {
    core::pose::symmetry::make_symmetric_movemap( pose, *local_movemap );
  }

	// start/best values
	pose::Pose start_pose( pose );
	pose::Pose best_pose( pose );
	utility::vector1< core::Real > best_score( filters_.size()+1, 1000000.0 );
	utility::vector1< core::Real > start_score;
	start_score.push_back( (*get_scorefxn())(start_pose) );
	for ( core::Size i=1; i<=filters_.size(); ++i ) {
		core::Real value( filters_[i].filter->report_sm( start_pose ) );
		if ( filters_[i].sample_low ) {
			start_score.push_back( value );
		} else {
			start_score.push_back( -value );
		}
	}

	// if necessary, create task op to find/design the worst region
	// save old task factory so it can be re-established later
	TaskFactoryOP old_task_factory( get_task_factory() );
	// create new task factory for modification
	TaskFactoryOP new_task_factory( new TaskFactory( *old_task_factory ) );

	// create a designaround operation to set all residues within 6 a of the residue in question as designable, everything else repackable.
	protocols::toolbox::task_operations::DesignAroundOperationOP design_around( create_worst_region_operation( pose ) );

	// check to see if we are disallowing a residue;
	// if we are, we should exit with FAIL_DO_NOT_RETRY status since changing a residue can affect montecarlo badly
	// monte carlo movers that use this mover should be enclosed in loopover movers which can be used to limit the number of disallowed positions and restart the monte carlos
	if ( get_last_move_status() == moves::FAIL_DO_NOT_RETRY ) {
		return;
	}

	utility::vector1< core::Size > not_allowed;
	//TR << "Num_redesigns=" << num_redesigns_ << std::endl;
	/*
	// sometimes, a taskop like layerdesign will restrict the identity current amino acid
	// There is then a chance that the designaround taskop added in the apply() method will restrict that residue to repacking
	// We need to look for amino acids set to repack-only which have a disallowed amino acid, and set them to designable
	core::pack::task::TaskFactory tmp_task_factory( *old_task_factory );

	if ( design_around ) {
		tmp_task_factory.push_back( design_around );
	}
	core::pack::task::PackerTaskOP tmp_task( tmp_task_factory.create_task_and_apply_taskoperations( pose ) );
	protocols::toolbox::task_operations::DesignAroundOperationOP des( new protocols::toolbox::task_operations::DesignAroundOperation() );
	des->design_shell( 0 );
	des->repack_shell( 8 );
	core::Size des_count( 0 );
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		if ( ! tmp_task->being_designed( i ) && ! tmp_task->bein ) {
			bool found( false );
			core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter res_iter;
			for ( res_iter = tmp_task->residue_task( i ).allowed_residue_types_begin();
						res_iter != tmp_task->residue_task( i ).allowed_residue_types_end();
						++res_iter ) {
				if ( pose.residue_type( i ).aa() == (*res_iter)->aa() ) {
					found = true;
					break;
				}
			}
			if ( !found ) {
				TR << "Designing/repacking residue " << i << " because it has a different type than any of the allowed types." << std::endl;
				//++des_count;
				des->include_residue( i );
				not_allowed.push_back( i );
			}
		}
		}
	if ( des_count > 0 ) {
		core::pack::task::TaskFactory des_task_factory( *old_task_factory );
		des_task_factory.push_back( des );
		core::pack::task::PackerTaskOP des_task( des_task_factory.create_task_and_apply_taskoperations( pose ) );
		protocols::simple_moves::PackRotamersMover designer( get_scorefxn(), des_task );
		designer.apply( pose );
		do_minimize( pose, 0.0001, local_movemap, get_scorefxn() );
		}*/

	if ( design_around ) {
		new_task_factory->push_back( design_around );
	}
	set_task_factory( new_task_factory );


	/*	utility::vector1< utility::vector1< bool > > mutations_by_pos;
 	// look for mutations that hurt psipred score and disallow them
	if ( true ) {
		// create psipred filter
		flxbb::filters::PsiPredFilterOP psipred( new flxbb::filters::PsiPredFilter( 0.0,
																																								"/work/javierbq/local/software/psipred/runpsipred_single",
																																								blueprint_file_,
																																								false ) );
		mutations_by_pos = check_and_disallow_mutations_by_filter( start_pose, get_task_factory()->create_task_and_apply_taskoperations( start_pose ), psipred );
	}
	*/
	// if requested, clear all residues marked as designable
	if ( clear_designable_residues_ ) {
		TR << "Clearing designable residues...";
		PackerTaskOP clear_task = get_task_factory()->create_task_and_apply_taskoperations( pose );
		for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
			if ( pose.residue( i ).is_protein() && clear_task->being_designed(i) ) {
				// skip glycine, proline
				if ( pose.residue( i ).name1() != 'G' || pose.residue( i ).name1() != 'P' ) {
					TR << i << "...";
					// mutate to alanine
					protocols::simple_moves::MutateResidue mut_res( i, 'A' );
					mut_res.apply( pose );
				} else {
					TR << "<" << i << ">...";
				}
			}
		}
		TR << std::endl;
		if ( dumpall_ ) {
			dump_pdb( pose, "clear_" );
		}
	}

	TR << "number of filters=" << filters_.size() << std::endl;

	// make a copy of the energy function too. Since we're going to be ramping around with weights,
	// we dont want to modify the existing scorefunction
	ScoreFunctionOP local_scorefxn( get_scorefxn()->clone() );

	// Remember the original weights - we're gonna be changing these during the ramp ups/downs
	core::scoring::EnergyMap full_weights = local_scorefxn()->weights();

	// Set up the packer task for packing.
	core::pack::task::PackerTaskOP task( create_packer_task( pose, local_movemap, not_allowed ) );

  protocols::simple_moves::PackRotamersMoverOP pack_full_repack_ = new protocols::simple_moves::PackRotamersMover( local_scorefxn, task );

	// If symmmetric pose then create a symmeteric rotamers mover
	if ( core::pose::symmetry::is_symmetric( pose ) )  {
    pack_full_repack_ = new protocols::simple_moves::symmetry::SymPackRotamersMover( local_scorefxn, task );
	}
	(*local_scorefxn)( pose );

	TR << *task << std::endl;

	/// prints something like this ***1***C***1*********2***C********3****C****2********3*****
	///                            **********xxxxxxxxxxxxx************************************
	if ( basic::options::option[ basic::options::OptionKeys::run::debug ]() ) {
		kinematics::simple_visualize_fold_tree_and_movemap_bb_chi( pose.fold_tree(),  *local_movemap, TR );
	}


	core::Size accept_count = 0;
	core::Size chk_counter = 0;

	int repeat_step=0;
	int repeat_count=-1;
	int total_repeat_count = 0;


	// Optain the native pose
	if( !get_native_pose() ) set_native_pose( new Pose( start_pose ) );

	// Statistic information
	std::vector< utility::vector1< core::Real > > best_score_log;
	std::vector< core::Real > curr_score_log;

	// for dry runs, literally just score the pose and quit.
	if( dry_run() ){
		(*local_scorefxn)( pose );
		return;
	}

	// Deal with disulphides - i have no idea what this does, ask Rob VErnon, he put this here.
	apply_disulfides(pose);

	int total_count=0;
	for ( core::Size ncmd = 0; ncmd < script_.size(); ncmd ++ ){
		total_count++;
		if ( basic::options::option[ basic::options::OptionKeys::run::debug ]() ) local_scorefxn->show( TR, pose );
		if ( basic::options::option[ basic::options::OptionKeys::run::debug ]() ) pose.constraint_set()->show_numbers( TR.Debug );

		// No MC is used, so update graphics manually
		#ifdef BOINC_GRAPHICS
			boinc::Boinc::update_graphics_current( pose );
		#endif

		RelaxScriptCommand cmd = script_[ncmd];

    TR.Debug << "Command: " << cmd.command << std::endl;

    if( cmd.command == "repeat" ){
			if( cmd.nparams < 1 ){ utility_exit_with_message( "ERROR: Syntax: " + cmd.command + "<number_of_repeats> " ); }
			repeat_count = (int) cmd.param1;
			repeat_step = ncmd;
		} else

		if( cmd.command == "endrepeat" ){
			TR.Debug << "CMD:  Repeat: " << repeat_count << std::endl;
			repeat_count -- ;
			total_repeat_count ++ ;
			if( repeat_count <= 0 ){}
			else{
				ncmd = repeat_step;
			}
		} else
		if( cmd.command == "dump" ){
			if( cmd.nparams < 1 ){ utility_exit_with_message( "ERROR: Syntax: " + cmd.command + "<number> " ); }
			pose.dump_pdb( "dump_" + right_string_of( (int) cmd.param1, 4, '0' ) + ".pdb" );
		}	else

		if( cmd.command.substr(0,7) == "dumpall" ){
				if( cmd.command.substr(8) == "true" )
					dumpall_ = true;
				if( cmd.command.substr(8) == "false" )
					dumpall_ = false;
		}	else

		if( cmd.command == "repack" ){
			//if( cmd.nparams < 0 ){ utility_exit_with_message( "More parameters expected after : " + cmd.command  ); }
			chk_counter++;
			std::string checkpoint_id = "chk" + string_of( chk_counter );
			if (!checkpoints_.recover_checkpoint( pose, get_current_tag(), checkpoint_id, true, true )){
				pack_full_repack_->apply( pose );
				checkpoints_.checkpoint( pose, get_current_tag(), checkpoint_id,  true );
			}
		}	else

		if( cmd.command == "min" ){
			if( cmd.nparams < 1 ){ utility_exit_with_message( "ERROR: Syntax " + cmd.command + " <min_tolerance>  " ); }

			chk_counter++;
			std::string checkpoint_id = "chk" + string_of( chk_counter );
			if (!checkpoints_.recover_checkpoint( pose, get_current_tag(), checkpoint_id, true, true )){
        do_minimize( pose, cmd.param1, local_movemap, local_scorefxn );
				checkpoints_.checkpoint( pose, get_current_tag(), checkpoint_id,  true );
			}
		}	else

		if( cmd.command.substr(0,5) == "scale" ){
			// no input validation as of now, relax will just die
			scoring::ScoreType scale_param = scoring::score_type_from_name(cmd.command.substr(6));
			local_scorefxn->set_weight( scale_param, full_weights[ scale_param ] * cmd.param1 );
		}   else

		if( cmd.command.substr(0,6) == "rscale" ){
				// no input validation as of now, relax will just die
				scoring::ScoreType scale_param = scoring::score_type_from_name(cmd.command.substr(7));
				local_scorefxn->set_weight( scale_param, full_weights[ scale_param ] * ((cmd.param2 - cmd.param1 ) * numeric::random::uniform() + cmd.param1 ));
		}   else

		if( cmd.command.substr(0,6) == "switch" ){
				// no input validation as of now, relax will just die
				if( cmd.command.substr(7) == "torsion" ) {
					TR << "Using AtomTreeMinimizer with dfp"  << std::endl;
					cartesian( false );
				} else if( cmd.command.substr(7) == "cartesian" ) {
					TR << "Using CartesianMinizer with lbfgs"  << std::endl;
					cartesian( true );
				}
		}   else

		if( cmd.command.substr(0,6) == "weight" ){
			// no input validation as of now, relax will just die
			scoring::ScoreType scale_param = scoring::score_type_from_name(cmd.command.substr(7));
			local_scorefxn->set_weight( scale_param, cmd.param1 );
			// I'm not too sure if the changing the default weight makes sense
			full_weights[ scale_param ] = cmd.param1;
		}   else
		if( cmd.command == "batch_shave" ){  // grab the score and remember the pose if the score is better then ever before.
		} else
		if( cmd.command == "show_weights" ){
			local_scorefxn->show(TR, pose);
		}   else

		if( cmd.command == "ramp_repack_min" ){
			if( cmd.nparams < 2 ){ utility_exit_with_message( "More parameters expected after : " + cmd.command  ); }

			// The first parameter is the relate repulsive weight
			local_scorefxn->set_weight( scoring::fa_rep, full_weights[ scoring::fa_rep ] * cmd.param1 );

			// The third paramter is the coordinate constraint weight
			if( constrain_coords() && (cmd.nparams >= 3) ){
				 local_scorefxn->set_weight( scoring::coordinate_constraint,  full_weights[ scoring::coordinate_constraint ] * cmd.param3 );
			}

			// we can also use the third parameter to ramp down design constraints
			if ( ramp_design_constraints_ && (cmd.nparams >= 3) ) {
				// TL: Added ramping of other constraints too -- do it the same way as the repulsive weights
				//core::Real const factor( 1 - cmd.param1 );
				core::Real const factor( cmd.param3 );
				local_scorefxn->set_weight( scoring::coordinate_constraint, full_weights[ scoring::coordinate_constraint ] * factor );
				local_scorefxn->set_weight( scoring::atom_pair_constraint, full_weights[ scoring::atom_pair_constraint ] * factor );
				local_scorefxn->set_weight( scoring::angle_constraint, full_weights[ scoring::angle_constraint ] * factor );
				local_scorefxn->set_weight( scoring::dihedral_constraint, full_weights[ scoring::dihedral_constraint ] * factor );
				TR << "[coordinate:atom_pair:angle:dihedral] = " << local_scorefxn->get_weight( scoring::coordinate_constraint ) << " : " << local_scorefxn->get_weight( scoring::atom_pair_constraint ) << " : " << local_scorefxn->get_weight( scoring::angle_constraint ) << " : " << local_scorefxn->get_weight( scoring::dihedral_constraint ) << std::endl;
			}

			// decide when to call ramady repair code
			if( total_repeat_count > 1 && repeat_count > 2 ){
				if( cmd.param1 < 0.2 ){
					if( do_rama_repair ){
						relax::fix_worst_bad_ramas( pose, ramady_num_rebuild_, 0.0, ramady_cutoff_, ramady_rms_limit_);
					}
				}
			}

			chk_counter++;
			std::string checkpoint_id = "chk" + string_of( chk_counter );
			if (!checkpoints_.recover_checkpoint( pose, get_current_tag(), checkpoint_id, true, true )){
				pack_full_repack_->apply( pose );
        do_minimize( pose, cmd.param2, local_movemap, local_scorefxn );
				checkpoints_.checkpoint( pose, get_current_tag(), checkpoint_id,  true );
			}


			// print some debug info
			if ( basic::options::option[ basic::options::OptionKeys::run::debug ]() ) {
				core::Real imedscore = (*local_scorefxn)( pose );
				core::pose::setPoseExtraScores( pose, "R" + right_string_of( total_count ,5,'0'), imedscore );
			}
		}	else


		if( cmd.command == "accept_to_best" ){
			// grab the score and remember the pose if the score is better then ever before.
			core::Real score = (*local_scorefxn)( pose );
			if ( rank_scorefxn_ ) {
				score = (*rank_scorefxn_)( pose );
			}

			bool improved( true );
			utility::vector1< core::Real > filt_scores;
			// first element of the vector is always the score
			filt_scores.push_back( score );

			// if there aren't filters, and the score is worse, this is not an improvement
			if ( filters_.size() == 0 && score > best_score[1] ) {
				improved = false;
			}
			// if any filters score worse, this is not an improvement
			for ( core::Size i=1; i<=filters_.size(); ++i ) {
				// apply each filter and compute score -- lower is better for all, so remember that
				core::Real value( filters_[i].filter->report_sm( pose ) );
				// if we want a higher score, make the value negative
				if ( ! filters_[i].sample_low ) {
					value = -value;
				}

				filt_scores.push_back( value );
				// if tolerance is being used, add/subtract the tol value for the final calculation
				if ( filters_[i].use_tolerance ) {
					if ( value > start_score[i+1] + filters_[i].tolerance ) {
						if ( ! filters_[i].use_threshold ||
								 ( filters_[i].use_threshold && value > filters_[i].threshold ) ) {
							improved = false;
							TR << "value for filter " << i << " is not improved; orig=" << best_score[i+1] << "; tol=" << filters_[i].tolerance << "; threshold=" << filters_[i].threshold << "; start=" << start_score[i+1] << std::endl;
						}
					}
				} else {
					if ( value > best_score[i+1] ) {
						if ( ! filters_[i].use_threshold ||
								 ( filters_[i].use_threshold && value > filters_[i].threshold ) ) {
							improved = false;
							TR << "value for filter " << i << " is not improved; orig=" << best_score[i+1] << "; threshold=" << filters_[i].threshold << std::endl;
						}
					}
				}
			}
			TR << filt_scores << std::endl;
			// if it's the first try, or if the result is improved, accept it
			if ( improved || ( accept_count == 0 ) ) {
				for ( core::Size i=1; i<=filt_scores.size(); ++i ) {
					best_score[i] = filt_scores[i];
				}
				best_pose = pose;
			}

			#ifdef BOINC_GRAPHICS
					boinc::Boinc::update_graphics_low_energy( best_pose, score  );
					boinc::Boinc::update_graphics_last_accepted( pose, score );
					boinc::Boinc::update_mc_trial_info( total_count , "FastDesign" );
			#endif
			core::Real rms = 0, irms = 0;
			if ( core::pose::symmetry::is_symmetric( pose ) && symmetric_rmsd_ ) {
				rms = CA_rmsd_symmetric( *get_native_pose() , best_pose );
				irms = CA_rmsd_symmetric( start_pose , best_pose );
			} else {
				rms = native_CA_rmsd( *get_native_pose() , best_pose );
				irms = native_CA_rmsd( start_pose , best_pose );
			}
			TR << "MRP: " << accept_count << "  " << score << "  " << best_score << "  " << rms << "  " << irms << "  " << std::endl;
			best_score_log.push_back( best_score );
			curr_score_log.push_back( score );

			accept_count++;
			if ( accept_count > script_max_accept_ ) break;
			if( test_cycles_ || dry_run() ) break;

		}	else


		if( cmd.command == "load_best" ){
			pose = best_pose;
		}	else
		if( cmd.command == "load_start" ){
			pose = start_pose;
		}	else

		if( cmd.command == "exit" ){
			utility_exit_with_message( "EXIT INVOKED" );
		}	else

		{
			utility_exit_with_message( "Unknown command: " + cmd.command );
		}


		if ( dumpall_ &&
				 ( cmd.command != "repeat" ) &&
				 ( cmd.command != "endrepeat" ) ) {
			dump_pdb( pose, cmd.command + "_" );
		}

		core::Real rms( -1.0 );
		if ( get_native_pose() ) {
			rms =  native_CA_rmsd( *get_native_pose() , pose );
			if ( core::pose::symmetry::is_symmetric( pose ) && symmetric_rmsd_ ) {
				rms = CA_rmsd_symmetric( *get_native_pose() , best_pose );
			}
		}
		core::Real irms =  native_CA_rmsd( start_pose , pose );
		if ( core::pose::symmetry::is_symmetric( pose ) && symmetric_rmsd_ ) {
			irms =  CA_rmsd_symmetric( start_pose , best_pose ); //need to make a symmetrical verision?
		}

		TR << "CMD: " <<  cmd.command << "  "
			 << (*local_scorefxn)( pose ) << "  ";
		if ( rank_scorefxn_ ) {
  	  TR << (*rank_scorefxn_)( pose ) << "  ";
		}
		TR << rms << "  "
			 << irms << "  "
			 << local_scorefxn->get_weight( scoring::fa_rep	)
			 << std::endl;

		// if we are allowing design, re-create the packer task
		if ( allow_design_ ) {
			// Set up the packer task for packing.
			task = create_packer_task( pose, local_movemap, not_allowed );

			// If symmmetric pose then create a symmeteric rotamers mover
			if ( core::pose::symmetry::is_symmetric( pose ) )  {
				pack_full_repack_ = new protocols::simple_moves::symmetry::SymPackRotamersMover( local_scorefxn, task );
			} else {
				pack_full_repack_ = new protocols::simple_moves::PackRotamersMover( local_scorefxn, task );
			}
		} // if allow_design
	}

	pose = best_pose;
	(*local_scorefxn)( pose );
	if ( rank_scorefxn_ ) {
		rank_scorefxn_->show( TR, pose );
	} else {
		local_scorefxn->show( TR, pose );
	}
	TR.flush();

	if ( basic::options::option[ basic::options::OptionKeys::run::debug ]() ) {
		for( Size j = 0; j < best_score_log.size(); j++ )
			for ( core::Size k=1; k<=best_score_log[j].size(); ++k )
				core::pose::setPoseExtraScores( pose, right_string_of(k,2,'0') + "B" + right_string_of(j,3,'0'), best_score_log[j][k] );

		for( Size j = 0; j < curr_score_log.size(); j ++ )
			core::pose::setPoseExtraScores( pose, "S" + right_string_of(j,3,'0'), curr_score_log[j] );
	}

	// restore task factory if we added anything
	if ( /*worst_region_mover_ ||*/ design_worst_ || design_by_psipred_ ) {
		set_task_factory( old_task_factory );
	}

	// compare to starting values.
	TR << "Starting values: " << start_score << std::endl;
	TR << "Final values   : " << best_score << std::endl;

	checkpoints_.clear_checkpoints();
}

std::string
FastDesign::get_name() const {
	return "FastDesign";
}


void FastDesign::read_script_file( const std::string &script_file, core::Size standard_repeats ){
	using namespace ObjexxFCL;
	script_.clear();
	std::vector< std::string > filelines;
	std::string line;

	runtime_assert( standard_repeats > 0 );
	if( script_file == "" ){
		TR << "================== Using default script ==================" << std::endl;
		filelines.push_back( "repeat " + string_of( standard_repeats )  );
		filelines.push_back( "ramp_repack_min 0.02  0.01     1.0"      );
		filelines.push_back( "ramp_repack_min 0.250 0.01     0.5"      );
		filelines.push_back( "ramp_repack_min 0.550 0.01     0.0"      );
		filelines.push_back( "ramp_repack_min 1     0.00001  0.0"      );
		filelines.push_back( "accept_to_best"                  );
		filelines.push_back( "endrepeat "                      );
	}else if (script_file == "NO CST RAMPING"){
		TR << "================== Using default script ==================" << std::endl;
		filelines.push_back( "repeat " + string_of( standard_repeats )  );
		filelines.push_back( "ramp_repack_min 0.02  0.01     1.0"      );
		filelines.push_back( "ramp_repack_min 0.250 0.01     1.0"      );
		filelines.push_back( "ramp_repack_min 0.550 0.01     1.0"      );
		filelines.push_back( "ramp_repack_min 1     0.00001  1.0"      );
		filelines.push_back( "accept_to_best"                  );
		filelines.push_back( "endrepeat "                      );
	}
	else{
		std::ifstream infile( script_file.c_str() );
		TR.Debug << "================== Reading script file: ==================" << std::endl;
		if (!infile.good()) {
			utility_exit_with_message( "[ERROR] Error opening script file '" + script_file + "'" );
		}
		while( getline(infile,line) ) {
			filelines.push_back( line );
		}
		infile.close();
	}

	int linecount=0;

	script_.clear();

	core::Size i;
	for( i =0; i< filelines.size(); i++ ){
		line = filelines[i];
		TR.Debug << line << std::endl;
		linecount++;
		utility::vector1< std::string > tokens ( utility::split( line ) );

		if ( tokens.size() > 0 ) {
			RelaxScriptCommand newcmd;
			newcmd.command = tokens[1];

			if (tokens.size() > 1) {newcmd.param1 = atof(tokens[2].c_str()); newcmd.nparams = 1;}
			if (tokens.size() > 2) {newcmd.param2 = atof(tokens[3].c_str()); newcmd.nparams = 2;}
			if (tokens.size() > 3) {newcmd.param3 = atof(tokens[4].c_str()); newcmd.nparams = 3;}
			if (tokens.size() > 4) {newcmd.param4 = atof(tokens[5].c_str()); newcmd.nparams = 4;}

			script_.push_back( newcmd );
		}
	}


}

void
FastDesign::makeDnaRigid( core::pose::Pose & pose, core::kinematics::MoveMapOP mm ){
	using namespace core::conformation;
	//if DNA present  so it doesn't move
	for ( Size i=1; i<=pose.total_residue() ; ++i )      {
		if( pose.residue(i).is_DNA()){
			TR << "turning off DNA bb and chi move" << std::endl;
			mm->set_bb( i, false );
			mm->set_chi( i, false );
		}
	}
	//end DNA rigidity
}

/// @brief outputs a pdb file with sequential numbering, with the specified prefix
void
FastDesign::dump_pdb( core::pose::Pose const & pose, std::string const & prefix ) {
	++dump_counter_;
	pose.dump_pdb( prefix + right_string_of( dump_counter_, 4, '0' ) + ".pdb" );
}

/// @brief initializes cache of allowed amino acids
void
FastDesign::initialize_aa_cache( core::pose::Pose const & pose )
{
	allowed_aas_.clear();
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		utility::vector1< bool > allowed_aas( core::chemical::num_canonical_aas, true );
		allowed_aas_.push_back( allowed_aas );
	}

	// initialize counts of redesigns
	num_redesigns_.clear();
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		utility::vector1< core::Size > counts( core::chemical::num_canonical_aas, 0 );
		num_redesigns_.push_back( counts );
	}

}

/// @brief initializes and creates the packer task
core::pack::task::PackerTaskOP
FastDesign::create_packer_task( core::pose::Pose const & pose,
																core::kinematics::MoveMapCOP /*local_movemap*/,
																utility::vector1< core::Size > const & /*allow_repacking*/ ) const
{
	core::pack::task::PackerTaskOP task;
	if ( get_task_factory() != 0 ) {
		task = get_task_factory()->create_task_and_apply_taskoperations( pose );
	} else {
		task = core::pack::task::TaskFactory::create_packer_task( pose ); //!!!!!!!!!!! << who tf put this comment here and why ?
	}

	/*utility::vector1<bool> allow_repack( pose.total_residue(), repack_ );
	if ( !basic::options::option[ basic::options::OptionKeys::relax::chi_move].user() ) {
		for ( Size pos = 1; pos<=pose.total_residue(); pos++ ) {
			allow_repack[ pos ] = local_movemap->get_chi( pos );
		}
	}

	if ( allow_design_ ) {
		task->initialize_from_command_line().restrict_to_residues( allow_repack );
	} else {
		task->initialize_from_command_line().restrict_to_repacking().restrict_to_residues( allow_repack );
		}*/
	task->initialize_from_command_line();
	task->or_include_current( true );

	// restrict design of disallowed amino acids
	if ( resfile_ != "" ) {
		core::pack::task::operation::ReadResfile resfile( resfile_ );
		resfile.apply( pose, *task );
		TR << "Read resfile and applied it." << std::endl;
	} else {
		for ( core::Size i=1; i<=allowed_aas_.size(); ++i ) {
			if ( ( max_redesigns_ > 0 ) && task->being_designed( i ) ) {
				//TR << "Restricting pos " << i << " to " << allowed_aas_[i] << std::endl;
				// count number of possible residues
				core::Size num_possible_mutations( task->nonconst_residue_task( i ).allowed_residue_types().size() );
				if ( num_possible_mutations <= 0 ) {
					TR << "WARNING: 0 mutations possible at residue " << i << std::endl;
				}
			}
			task->nonconst_residue_task( i ).restrict_absent_canonical_aas( allowed_aas_[i] );
		}
	}

	/*// if there are residues that need to be allowed to repack, do it
	if ( allow_repacking.size() ) {
		core::pack::task::operation::OperateOnCertainResidues oocr_repacking;
		for ( core::Size i=1; i<=allow_repacking.size(); ++i ) {
			TR << " Allowing repacking for residue " << allow_repacking[i] << std::endl;
			// allow native amino acid
			task->nonconst_residue_task( allow_repacking[i] ).allow_aa( pose.residue( allow_repacking[i] ).aa() );
		}
		oocr_repacking.op( new core::pack::task::operation::RestrictToRepackingRLT );
		oocr_repacking.residue_indices( allow_repacking );
		oocr_repacking.apply( pose, *task );
		}*/

	if( limit_aroma_chi2() ) {
		protocols::toolbox::task_operations::LimitAromaChi2Operation lp_op;
		lp_op.apply( pose, *task );
	}

	return task;
}

/// @brief checks each mutation at each designable position, and disallows those mutations that hurt the filter score
utility::vector1< utility::vector1< bool > >
FastDesign::check_and_disallow_mutations_by_filter( core::pose::Pose const & pose, core::pack::task::PackerTaskOP task, protocols::filters::FilterOP filter ) const
{
	TR << *task << std::endl;
	utility::vector1< utility::vector1< bool > > mutations_by_pos;
	core::Real const start_score( filter->report_sm( pose ) );
	// iterate through all positions
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		utility::vector1< bool > allowed_types( core::chemical::num_canonical_aas, false );
		// look through residue types allowed at this position
		core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter res_iter( task->residue_task( i ).allowed_residue_types_begin() );
		while ( res_iter != task->residue_task( i ).allowed_residue_types_end() ) {
			TR << "Testing " << (*res_iter)->name3() << " at pos " << i << std::endl;
			protocols::simple_moves::MutateResidue mut_res( i, (*res_iter)->name() );
			core::pose::Pose tmp_pose( pose );
			mut_res.apply( tmp_pose );
			core::Real value( filter->report_sm( tmp_pose ) );
			if ( value >= start_score ) {
				TR << "Allowing; start=" << start_score << "; this=" << value << std::endl;
				allowed_types[ core::chemical::aa_from_oneletter_code( (*res_iter)->name1() ) ] = true;
			} else {
				TR << "Rejecting; start=" << start_score << "; this=" << value << std::endl;
			}
			++res_iter;
		}

		mutations_by_pos.push_back( allowed_types );
		task->nonconst_residue_task( i ).restrict_absent_canonical_aas( allowed_types );
	}
	return mutations_by_pos;
}


/// @brief increments the number of times a particular residue (and aa) have been designed on
void
FastDesign::increment_num_redesigns( core::pose::Pose const & pose, core::Size const seqpos )
{
	++num_redesigns_[ seqpos ][ pose.residue( seqpos ).aa() ];
}


/// @brief checks how many times the residue at seqpos has been redesigned. If the value is larger than than the maximum allowed number of redesigns, return false, otherwise return true
bool
FastDesign::check_num_redesigns( core::pose::Pose const & pose, core::Size const seqpos ) const
{
	core::Size const value( num_redesigns_[ seqpos ][ pose.residue( seqpos ).aa() ] );
	return value >= max_redesigns_;
}

/// @brief Creates and returns a pointer to a taskoperation that designs the worst regions of a protein only. Returns NULL if we don't need the taskop
protocols::toolbox::task_operations::DesignAroundOperationOP
FastDesign::create_worst_region_operation( core::pose::Pose const & /*pose*/ )
{
	/*core::Real const shell_size( 6.0 );
	utility::vector1< flxbb::filters::HighestEnergyRegionOperationOP > worst_res_ops;
	if ( design_worst_ ) {
		flxbb::filters::HighestEnergyRegionOperationOP high_energy( new flxbb::filters::HighestEnergyRegionOperation() );
		high_energy->set_scorefxn( get_scorefxn() );
		high_energy->set_regions_to_design( pose.total_residue() );
		high_energy->set_region_shell( shell_size ); // to make it the same as the designaround operation
		worst_res_ops.push_back( high_energy );
	}
	if ( design_by_psipred_ ) {
		worst_res_ops.push_back( new flxbb::filters::DesignBySecondaryStructureOperation(
																	 blueprint_file_,
																	 "/work/javierbq/local/software/psipred/runpsipred_single",
																	 shell_size,
																	 1000.0,
																	 pose.total_residue(),
																	 true,
																	 false ) );
	}
	if ( design_by_frag_qual_ ) {
		// if we haven't created the frag qual taskop yet, create it
		if ( ! frag_qual_op_ ) {
			frag_qual_op_ = new flxbb::filters::DesignByFragmentQualityOperation(
											 "/work/javierbq/local/software/psipred/runpsipred_single" );
			frag_qual_op_->set_region_shell( shell_size );
		}
		worst_res_ops.push_back( frag_qual_op_ );
		}

	// if no task ops are created to select residues, we don't need this taskop and should return null
	if ( worst_res_ops.size() == 0 && ! worst_region_mover_ && ! only_design_changes_ ) {
		return NULL;
	}

	// now we can just apply a DesignAround operation using the residues that don't match
	protocols::toolbox::task_operations::DesignAroundOperationOP design_around;
	design_around = new protocols::toolbox::task_operations::DesignAroundOperation();
	design_around->design_shell( shell_size + 2 );
	design_around->repack_shell( 1000.0 );

	for ( core::Size i=1; i<=worst_res_ops.size(); ++i ) {
		// get list of residues that are worst
		utility::vector1< core::Size > sorted_resids( worst_res_ops[i]->get_residues_to_design( pose ) );

		// check for residues that have been designed too much
		core::Size j( 1 );
		utility::vector1< core::Size > resids_to_design;
		while ( ( j <= sorted_resids.size() ) && ( resids_to_design.size() < regions_to_design_ ) ) {
			resids_to_design.push_back( sorted_resids[j] );
		}

		// make a note of which residues we are designing (and their identities)
		for ( core::Size j=1; j<=resids_to_design.size(); ++j ) {
			increment_num_redesigns( pose, resids_to_design[j] );

			// if this is our first time designing this residue, disallow whatever the bad mutation is
			core::Size count( 0 );
			for ( core::Size aa=1; aa<=num_redesigns_[ resids_to_design[j] ].size(); ++aa ) {
				count += num_redesigns_[ resids_to_design[j] ][ aa ];
			}
			TR << " Total redesigns=" << count << std::endl;
			if ( count <= 1 ) {
				TR << "Disallowing since this is the first time redesigning at residue " << resids_to_design[j] << std::endl;
				num_redesigns_[ resids_to_design[j] ][ pose.residue( resids_to_design[j] ).aa() ] = max_redesigns_;
			}
			if ( ( max_redesigns_ > 0 ) && check_num_redesigns( pose, resids_to_design[j] ) ) {
				TR << "Disallowing " << pose.residue( resids_to_design[j] ).name1() << " at position " << resids_to_design[j] << std::endl;
				allowed_aas_[ resids_to_design[j] ][ pose.residue( resids_to_design[j] ).aa() ] = false;
				//restrict_aa( pose, resids_to_design[j] );
				//set_last_move_status( moves::FAIL_DO_NOT_RETRY );
			}
		}

		TR << "Residues to design are: ";
		for ( core::Size j=1; j<=resids_to_design.size(); ++j ) {
				TR << resids_to_design[j] << " ";
				design_around->include_residue( resids_to_design[j] );
		}
		TR << std::endl;
	}

	if ( worst_region_mover_ ) {
		utility::vector1< core::Size > const & resi( worst_region_mover_->last_residues_restricted() );
		TR << "Restriction mover reports residue=" << resi << std::endl;
		if ( resi.size() == 0 ) {
			//TODO: fix this stupid testing default
			design_around->include_residue( 2 );
		} else {
			for ( core::Size res=1; res<=resi.size(); ++res ) {
				design_around->include_residue( resi[res] );
			}
		}
	}
	return design_around;
	*/
	return NULL;
}


} // namespace denovo_design
} // namespace devel
