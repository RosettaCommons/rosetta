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
/// @details
/// @author Tom Linsky
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- Added support for D-amino acids.

//Unit Headers
#include <devel/denovo_design/FastDesign.hh>
#include <devel/denovo_design/FastDesignCreator.hh>

//Project Headers

//Core Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

//Protocol Headers
#include <basic/datacache/DataMap.hh>
#include <protocols/constraint_generator/AddConstraints.hh>
#include <protocols/constraint_generator/ConstraintGenerator.hh>
#include <protocols/constraint_generator/RemoveConstraints.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/toolbox/task_operations/LimitAromaChi2Operation.hh>

//Basic Headers
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

//Utility Headers
#include <utility/tag/Tag.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

//C++ Headers

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

static THREAD_LOCAL basic::Tracer TR( "devel.denovo_design.FastDesign" );

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
	return protocols::moves::MoverOP( new FastDesign() );
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
	FastRelax(),
	design_worst_( false ),
	design_by_psipred_( false ),
	design_by_frag_qual_( false ),
	only_design_changes_( false ),
	blueprint_file_( "" ),
	clear_designable_residues_( false ),
	rank_scorefxn_( /* NULL */ ),
	max_redesigns_( 5 ),
	dump_counter_( 0 ),
	//frag_qual_op_( NULL ),
	regions_to_design_( 1 ),
	run_count_( 0 ),
	cached_sequence_( "" ),
	cgs_()
	//worst_region_mover_( NULL )
{
	filters_.clear();
	allowed_aas_.clear();
	num_redesigns_.clear();
	set_enable_design( true );
	//read_script_file( "", default_repeats_ );
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
FastDesign::~FastDesign() {}


/// Return a copy of ourselves
protocols::moves::MoverOP
FastDesign::clone() const {
	return protocols::moves::MoverOP( new FastDesign(*this) );
}

void
FastDesign::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
) {
	using core::pack::task::operation::TaskOperationCOP;
	// make sure we create a task factory before parsing FastRelax::parse_my_tag
	// otherwise no design will occur
	core::pack::task::TaskFactoryOP local_tf( new core::pack::task::TaskFactory() );
	if ( get_task_factory() ) {
		local_tf = get_task_factory()->clone();
	} else {
		local_tf->push_back(TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline() ));
		if ( basic::options::option[ basic::options::OptionKeys::relax::respect_resfile]() &&
				basic::options::option[ basic::options::OptionKeys::packing::resfile].user() ) {
			local_tf->push_back(TaskOperationCOP( new core::pack::task::operation::ReadResfile() ));
			TR << "Using Resfile for packing step. " <<std::endl;
		} else {
			core::pack::task::operation::PreventRepackingOP turn_off_packing( new core::pack::task::operation::PreventRepacking() );
			for ( Size pos = 1; pos <= pose.total_residue(); ++pos ) {
				if ( ! get_movemap()->get_chi(pos) ) {
					turn_off_packing->include_residue(pos);
				}
			}
			local_tf->push_back(turn_off_packing);
		}
	}
	//Include current rotamer by default - as before.
	local_tf->push_back(TaskOperationCOP( new core::pack::task::operation::IncludeCurrent() ));

	if ( limit_aroma_chi2() ) {
		local_tf->push_back(TaskOperationCOP( new protocols::toolbox::task_operations::LimitAromaChi2Operation() ));
	}
	set_task_factory( local_tf );
	FastRelax::parse_my_tag( tag, data, filters, movers, pose );

	design_worst_ = tag->getOption< bool >( "only_design_worst_region", design_worst_ );
	design_by_psipred_ = tag->getOption< bool >( "design_by_psipred", design_by_psipred_ );
	design_by_frag_qual_ = tag->getOption< bool >( "design_by_frag_qual", design_by_frag_qual_ );
	only_design_changes_ = tag->getOption< bool >( "design_changes", only_design_changes_ );
	blueprint_file_ = tag->getOption< std::string >( "blueprint", blueprint_file_ );
	// blueprint file must be specified if design_by_psipred is set
	if ( blueprint_file_ == "" && design_by_psipred_ ) {
		utility_exit_with_message( "Blueprint file must be specified if design_by_psipred is set." );
	}
	clear_designable_residues_ = tag->getOption< bool >( "clear_designable_residues", clear_designable_residues_ );
	std::string const rank_scorefxn_name( tag->getOption< std::string >( "rank_scorefxn", "" ) );
	if ( rank_scorefxn_name == "" ) {
		rank_scorefxn_ = NULL;
	} else {
		rank_scorefxn_ = data.get< core::scoring::ScoreFunction * >( "scorefxns", rank_scorefxn_name )->clone();
	}

	max_redesigns_ = tag->getOption< core::Size >( "max_redesigns", max_redesigns_ );

	utility::vector1< std::string > cgs = utility::string_split( tag->getOption< std::string >( "cgs", "" ), ',' );
	for ( utility::vector1< std::string >::const_iterator cg=cgs.begin(); cg!=cgs.end(); ++cg ) {
		if ( cg->empty() ) continue;
		protocols::constraint_generator::ConstraintGeneratorCOP new_cg = data.get_ptr< protocols::constraint_generator::ConstraintGenerator const >( "CONSTRAINT_GENERATORS", *cg );
		if ( !new_cg ) {
			std::stringstream msg;
			msg << "FastDesign: Could not find a constraint generator named " << *cg << " in the data map.  Ensure it has been defined in an AddConstraints mover before being referenced by FastDesign."
				<< std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
		}
		cgs_.push_back( new_cg );
	}

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

	utility::vector1< utility::tag::TagCOP > const & filter_tags( tag->getTags() );
	for ( core::Size i=1; i<=filter_tags.size(); ++i ) {
		utility::tag::TagCOP const layer_tag( filter_tags[i] );
		// AMW: cppcheck notes this is never used
		//std::string op = layer_tag->getName(); // Should be "AND" just like in Generic MC
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

	initialize_aa_cache( pose );
}

/// @brief helper function to print out a vector1 of reals
void
print_list( basic::Tracer & tracer, utility::vector1< core::Real > const & list )
{
	tracer << "[ ";
	for ( core::Size i=1; i<=list.size(); ++i ) {
		tracer << list[i] << " ";
	}
	tracer << "]" << std::endl;
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

	TR.Debug   << "========================== FastDesign Run: " << run_count_ << " ===================" << std::endl;
	reset_status();

	// start/best values
	core::pose::Pose start_pose( pose );
	core::pose::Pose best_pose( pose );
	utility::vector1< core::Real > best_score;
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

	// check to see if we are disallowing a residue;
	// if we are, we should exit with FAIL_DO_NOT_RETRY status since changing a residue can affect montecarlo badly
	// monte carlo movers that use this mover should be enclosed in loopover movers which can be used to limit the number of disallowed positions and restart the monte carlos
	if ( get_last_move_status() == moves::FAIL_DO_NOT_RETRY ) {
		return;
	}

	utility::vector1< core::Size > not_allowed;
	//TR << "Num_redesigns=" << num_redesigns_ << std::endl;

	// create and print task
	core::pack::task::PackerTaskOP clear_task = get_task_factory()->create_task_and_apply_taskoperations( pose );
	clear_task->show( TR );

	// if requested, clear all residues marked as designable
	if ( clear_designable_residues_ ) {
		TR << "Clearing designable residues...";
		for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
			if ( pose.residue( i ).is_protein() && clear_task->being_designed(i) ) {
				if ( pose.residue(i).type().is_alpha_aa() ) {
					if ( pose.residue( i ).type().is_l_aa() ) { //Note that glycine is skipped, since it's not an L-amino acid.
						// skip proline
						if ( pose.residue( i ).type().aa() != core::chemical::aa_pro ) {
							TR << i << "...";
							// mutate to alanine
							protocols::simple_moves::MutateResidue mut_res( i, core::chemical::aa_ala );
							mut_res.apply( pose );
						} else {
							TR << "<" << i << ">...";
						}
					} else if ( pose.residue( i ).type().is_d_aa() ) {
						// skip D-proline
						if ( pose.residue( i ).type().aa() != core::chemical::aa_dpr ) {
							TR << i << "...";
							// mutate to D-alanine
							protocols::simple_moves::MutateResidue mut_res( i, core::chemical::aa_dal );
							mut_res.apply( pose );
						} else {
							TR << "<" << i << ">...";
						}
					} else { //Achiral case
						TR << "<" << i << ">...";
					}
				}
			}
		}
		TR << std::endl;
	}

	TR.Debug << "number of filters=" << filters_.size() << std::endl;

	// support for ramping reference weights
	modify_scripts_for_alternative_scorefunctions();

	FastRelax::apply( pose );

	best_score.push_back( (*get_scorefxn())(pose) );
	for ( core::Size i=1; i<=filters_.size(); ++i ) {
		core::Real value( filters_[i].filter->report_sm( start_pose ) );
		if ( filters_[i].sample_low ) {
			best_score.push_back( value );
		} else {
			best_score.push_back( -value );
		}
	}
	if ( rank_scorefxn_ ) {
		rank_scorefxn_->show( TR, pose );
	} else {
		get_scorefxn()->show( TR, pose );
	}
	TR.flush();

	// compare to starting values.
	TR.Debug << "Starting values: ";
	print_list( TR, start_score );
	TR.Debug << "Ending values  : ";
	print_list( TR, best_score );
}

std::string
FastDesign::get_name() const {
	return "FastDesign";
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

/// @brief sets constraint weights -- used with constraint ramping
void
FastDesign::set_constraint_weight(
	core::scoring::ScoreFunctionOP local_scorefxn,
	core::scoring::EnergyMap const & full_weights,
	core::Real const weight,
	core::pose::Pose & pose ) const
{
	runtime_assert( local_scorefxn != 0 );
	if ( cgs_.size() ) {
		protocols::constraint_generator::RemoveConstraints rm_csts( cgs_ );

		try {
			rm_csts.apply( pose );
		} catch ( protocols::constraint_generator::EXCN_RemoveCstsFailed const & e ) {
			// if removing constraints fails, we don't really care
			// they are only being removed to clean up before re-adding them below
		}
		if ( weight > 1e-6 ) {
			protocols::constraint_generator::AddConstraints( cgs_ ).apply( pose );
		}

	} else {
		local_scorefxn->set_weight( core::scoring::coordinate_constraint, full_weights[ core::scoring::coordinate_constraint ] * weight );
		local_scorefxn->set_weight( core::scoring::atom_pair_constraint, full_weights[ core::scoring::atom_pair_constraint ] * weight );
		local_scorefxn->set_weight( core::scoring::angle_constraint, full_weights[ core::scoring::angle_constraint ] * weight );
		local_scorefxn->set_weight( core::scoring::dihedral_constraint, full_weights[ core::scoring::dihedral_constraint ] * weight );
		TR << "[coordinate:atom_pair:angle:dihedral] = " << local_scorefxn->get_weight( core::scoring::coordinate_constraint ) << " : " << local_scorefxn->get_weight( core::scoring::atom_pair_constraint ) << " : " << local_scorefxn->get_weight( core::scoring::angle_constraint ) << " : " << local_scorefxn->get_weight( core::scoring::dihedral_constraint ) << std::endl;
	}
}

void
FastDesign::modify_scripts_for_alternative_scorefunctions()
{
	using namespace basic::options;

	// will attempt to modify the relax script for beta_nov15 only
	// options for any other score functions could be also added below
	std::vector< std::string > filelines;

	if ( !(FastRelax::script_file_specified_) &&
			(option[ OptionKeys::corrections::beta_nov15 ]() || option[ OptionKeys::corrections::beta_nov15_cart ] ) ) {

		// hard-coded reference weights for now...
		TR << "Calling correction for beta_nov15..." << std::endl;
		filelines.push_back( "repeat 4"                                );
		filelines.push_back( "reference 0.3     3.0     -2.2     -2.7     5.0     -0.2    0.5     5.2     -0.5     4.5     5.0     -1.8     -0.3     -1.5     0.0     -0.8     1.0     4.0     9.5     3.5" );
		filelines.push_back( "ramp_repack_min 0.02  0.01     1.0"      );

		filelines.push_back( "reference 1.32468 3.74979 -2.16074 -2.40953 1.99829 0.41816 0.27935 3.02374 -0.54958 2.20647 2.50235 -1.34026 -1.17821 -1.31095 0.08526 -0.27469 1.31675 3.12269 3.24099 1.18223" );
		filelines.push_back( "ramp_repack_min 0.250 0.01     0.5"      );

		filelines.push_back( "reference 1.32468 3.5979  -2.16074 -2.50953 1.79829 0.51816 0.17935 2.82374 -0.54958 2.00647 2.30235 -1.34026 -1.27821 -1.31095 0.08526 -0.27469 1.31675 3.02269 2.94099 1.00223" );
		filelines.push_back( "ramp_repack_min 0.550 0.01     0.0"      );

		filelines.push_back( "reference 1.32468 3.25479 -2.14574 -2.72453 1.21829 0.79816 -0.30065 2.30374 -0.71458 1.66147 1.65735 -1.34026 -1.64321 -1.45095 -0.09474 -0.28969 1.15175 2.64269 2.26099 0.58223" );
		filelines.push_back( "ramp_repack_min 1     0.00001  0.0"      );

		filelines.push_back( "accept_to_best"                  );
		filelines.push_back( "endrepeat "                      );
	}

	if ( filelines.size() > 0 ) {
		FastRelax::set_script_from_lines( filelines );
	}
}

} // namespace denovo_design
} // namespace devel
