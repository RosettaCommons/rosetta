// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/switches/GraftSwitchMover.cc
/// @brief Mover to graft a fucntional sequence optimally on a helical LOCKR type latch
/// @author Bobby Langan (robert.langan@gmail.com)

// Unit headers
#include <protocols/switches/GraftSwitchMover.hh>
#include <protocols/switches/GraftSwitchMoverCreator.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>

// core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <core/select/util.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/MoveMap.hh>

// unit headers
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/simple_moves/SimpleThreadingMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <protocols/simple_filters/BuriedUnsatHbondFilter.hh>
#include <protocols/simple_filters/ShapeComplementarityFilter.hh>

// numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <boost/filesystem/path.hpp>
#include <boost/iterator/indirect_iterator.hpp>

// Utility Headers
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/memory.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>

#define PACKER_ITER 1

namespace protocols {
namespace switches {

using basic::Warning;
using basic::t_warning;
static basic::Tracer TR( "protocols.switches.GraftSwitchMover" );

// GraftSwitchMoverCreator

std::string
GraftSwitchMoverCreator::keyname() const {
	return GraftSwitchMover::mover_name();
}

moves::MoverOP
GraftSwitchMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new GraftSwitchMover );
}

void GraftSwitchMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	GraftSwitchMover::provide_xml_schema( xsd );
}

void GraftSwitchMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "start", xsct_non_negative_integer, "Residue number defining the start of the threadable residue range.  Must be defined with 'end'. Overrides default behavior.", "0")
		+ XMLSchemaAttribute::attribute_w_default( "end", xsct_non_negative_integer, "Residue number to end the range of residues to allow threading.", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "n_term", xsct_rosetta_bool, "Default behaivor without setting start/end or selector is to find the c-terminal helix and set all positions designable.  Set to true if latch is on the n-terminus instead.", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "any_order", xsct_rosetta_bool, "With multiple sequences, set to false if you want to maintain the order of sequences you've defined.  Useful if you need one seq closer to one terminus than the other.", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "pack_neighbors", xsct_rosetta_bool, "Set to false if you do not want neighboring residues to repack during threading", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "graft_on_latch_loop", xsct_rosetta_bool, "Set to false if you do not want the loop between the cage/latch to be threaded", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "sequence", xs_string, "Comma separated list of sequences to thread. A '-' indicates a gap where NATAA will be default", "" )
		+ XMLSchemaAttribute::attribute_w_default( "important_residues", xs_string, "A comma-separated list of residue numbers (1 corresponding to first residue in 'sequence') that must be buried/caged. Separate important residues for separate sequences with a '/' in order listed for the 'sequences' tag", "" )
		+ XMLSchemaAttribute::attribute_w_default( "burial_cutoff", xsct_real, "Number of neighboring residues that determines extent of burial", "6")
		+ XMLSchemaAttribute::attribute_w_default( "selector", xs_string, "Residue Selector defining residues that can be threaded onto.  Overrides default behavior and start/end if defined.", "")
		+ XMLSchemaAttribute::attribute_w_default( "pack_min", xsct_rosetta_bool, "Do a final minimization and pack neighbors after threading", "true");
	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	rosetta_scripts::attributes_for_parse_score_function( attlist );

	moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Grafts a functional sequence onto the latch of a LOCKR type protein. Returns all threads where average degree of important residues is greater than the cutoff (defualt=6).  By default will find the c-terminal helix and set all positions as threadable.  Set n_term to true to find n-terminal helix instead.  Set start/end to define a custom range to begin threads.  Set selector to use a residue selector to define threadable residues and Rosetta will determine which positions can fit your sequence.", attlist );
}

// GraftSwitchMover

GraftSwitchMover::GraftSwitchMover() :
	protocols::moves::Mover("GraftSwitchMover"),
	start_graft_(0),
	end_graft_(0),
	n_terminus_(false),
	pack_neighbors_(true),
	graft_loop_(true),
	find_helix_(false),
	seq_in_any_order_(true),
	pack_min_(true),
	burial_cutoff_(6),
	output_list_()
{
}

GraftSwitchMover::~GraftSwitchMover() = default;

GraftSwitchMover::GraftSwitchMover( GraftSwitchMover const & other ) :
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( other )
{
	scorefxn_ = other.scorefxn_;
	task_factory_ = other.task_factory_;
	start_graft_ = other.start_graft_;
	end_graft_ = other.end_graft_;
	n_terminus_ = other.n_terminus_;
	pack_neighbors_ = other.pack_neighbors_;
	seq_in_any_order_ = other.seq_in_any_order_;
	pack_min_ = other.pack_min_;
	sequences_ = other.sequences_;
	important_residues_ = other.important_residues_;
	graft_loop_ = other.graft_loop_;
	burial_cutoff_ = other.burial_cutoff_;
	find_helix_ = other.find_helix_;
	threadable_residues_ = other.threadable_residues_;
	threadable_residues_selector_ = other.threadable_residues_selector_;
}

void
GraftSwitchMover::apply( Pose & pose )
{
	orig_pose_ = utility::pointer::make_shared< core::pose::Pose >( pose );

	if ( sequences_.size() == 0 ) {
		utility_exit_with_message("Sequence not set for threading.  Cannot continue.");
	}

	if ( !scorefxn_ ) {
		TR << "Using default score function" << std::endl;
		scorefxn_ = core::scoring::get_score_function();
	}

	//Get set of threadable residues from ResidueSelector
	if ( threadable_residues_selector_ ) {
		threadable_residues_ = core::select::get_residue_set_from_subset( threadable_residues_selector_->apply( pose ) );
	} else {
		//If no selector defined, and parse_my_tag determines it needs to find the terminal helix, do so
		if ( find_helix_ ) {
			if ( !n_terminus_ ) {
				TR << "Finding C-terminal Helix" << std::endl;
				start_graft_ = get_c_terminal_helix(orig_pose_);
				end_graft_ = orig_pose_->size();
			} else {
				TR << "Finding N-terminal Helix" << std::endl;
				start_graft_ = 1;
				end_graft_ = get_n_terminal_helix(orig_pose_);
			}
		}
		//If the mover doesn't find the helix itself, use the start and end definitions
		for ( core::Size i = start_graft_; i <= end_graft_; i++ ) {
			threadable_residues_.insert(i);
		}
	}

	//Series of tests to determine if threadable residues is viable for sequences defined
	if ( threadable_residues_.empty() ) utility_exit_with_message("No threadable residues found.");
	TR << "Residues available for threading are: ";
	for ( std::set<core::Size>::iterator i = threadable_residues_.begin(); i != threadable_residues_.end(); i++ ) {
		TR << *i << ", ";
	}
	TR << std::endl;
	core::Size number_residues_to_thread = 0;
	TR << "You are threading sequences: ";
	for ( core::Size i=1; i <= sequences_.size(); i++ ) {
		TR << sequences_[i] << ", ";
		number_residues_to_thread += sequences_[i].size();
	}
	TR << std::endl;
	TR << (seq_in_any_order_ ? "Threading sequences in any order." : "Threading sequences in order specified.") << std::endl;
	TR << "Threading " << number_residues_to_thread << " residues from " << sequences_.size() << " sequences onto " << threadable_residues_.size() << " designable residues." << std::endl;
	if ( threadable_residues_.size() < number_residues_to_thread ) utility_exit_with_message("The set of threadable residues is smaller than the length of the sequence(s) you want to thread.");

	//Find all possible starting positions for threading (throw error if none possible)
	std::list< utility::vector1< core::Size > > starting_positions = all_threading_start_combinations();
	if ( starting_positions.empty() ) utility_exit_with_message("The set of threadable residues is to stringent to fit the sequence(s) you want to thread. Could not find positions to start threading.");

	//make a copy of the input scorefunction b/c we're messing with the weights for initial soft scoring
	ScoreFunctionOP sfxn_soft( scorefxn_->clone() );
	TR << "Setting up soft scorefxn" << std::endl;
	sfxn_soft->set_weight( core::scoring::ScoreType::fa_rep, 1.0 );
	sfxn_soft->set_weight( core::scoring::ScoreType::fa_intra_rep_xover4, 1.0 ); //soft weights based off beta_nov16_soft compared to beta_nov16

	//Set up Threading mover and taskops
	TR << "Setting up SimpleThreadingMover" << std::endl;
	protocols::simple_moves::SimpleThreadingMoverOP threader = protocols::simple_moves::SimpleThreadingMoverOP( new protocols::simple_moves::SimpleThreadingMover() );
	threader->set_pack_neighbors(pack_neighbors_);
	threader->set_scorefxn( sfxn_soft );

	//do the threading for each possible combination of threading positions
	core::pose::PoseOP mypose = utility::pointer::make_shared< core::pose::Pose >( *orig_pose_ );
	std::list< utility::vector1< core::Size > >::iterator i = starting_positions.begin();
	while (  i != starting_positions.end() ) {
		utility::vector1<core::Size> thread_pos = *i;
		bool all_good = true;
		TR << "Trying set: " << thread_pos << std::endl;
		for ( core::Size j=1; j <= thread_pos.size(); j++ ) {
			std::string seq = sequences_[j];
			core::Size position = thread_pos[j];

			TR << "Starting thread of seq " << seq << " at residue " << position << std::endl;

			threader->set_sequence(seq, position);
			threader->apply( *mypose );
			if ( !init_burial_filter( *mypose, j, position ) ) {
				all_good = false;
				TR << "Burial insufficient for this sequence, removing redundancy and moving on. Erasing: " << *i << ", ";
				i = starting_positions.erase(i);
				std::list< utility::vector1< core::Size > >::iterator k = i;
				while ( k != starting_positions.end() ) {
					if ( (*k)[j] == position && *k == *i ) {
						TR << *k << ", ";
						i = k = starting_positions.erase(k);
					} else if ( (*k)[j] == position ) {
						TR << *k << ", ";
						k = starting_positions.erase(k);
					} else k++;
				}
				TR << std::endl;
				break;
			}
		}
		if ( all_good ) output_list_.push_back( *(i++) );
		mypose = utility::pointer::make_shared< core::pose::Pose >( *orig_pose_ );
	}

	if ( output_list_.size() == 0 ) {
		TR << "No positions found with adequite burial. Change neighboring residues burial cutoff or expand the range of threadable residues." << std::endl;
		set_last_move_status( moves::FAIL_RETRY );
		return;
	} else {
		set_last_move_status( moves::MS_SUCCESS );
	}

	utility::vector1<core::Size> thread_position( output_list_.front() );
	output_list_.pop_front();

	//Output first good thread (here we do the computationally hard minimization/packing)
	TR << "Outputting thread with ";
	for ( core::Size i=1; i <= thread_position.size(); i++ ) {
		TR << sequences_[i] << " at " << thread_position[i] << ", ";
	}
	TR << std::endl;

	for ( core::Size i=1; i <= thread_position.size(); i++ ) {
		threader->set_sequence(sequences_[i], thread_position[i]);
		threader->apply( pose );
	}

	if ( pack_min_ ) {
		//TaskOps
		core::pack::task::operation::PreventRepackingOP turn_off_packing = core::pack::task::operation::PreventRepackingOP( new core::pack::task::operation::PreventRepacking() );
		core::pack::task::operation::RestrictToRepackingOP turn_off_design = core::pack::task::operation::RestrictToRepackingOP( new core::pack::task::operation::RestrictToRepacking() );

		//Set up Minimizer
		core::kinematics::MoveMapOP mm = core::kinematics::MoveMapOP( new core::kinematics::MoveMap() );
		mm->set_bb( true );
		mm->set_chi( true );
		protocols::minimization_packing::MinMoverOP min = protocols::minimization_packing::MinMoverOP( new protocols::minimization_packing::MinMover( mm, scorefxn_, "lbfgs_armijo_nonmonotone", 0.0001, true ) );

		//Create vector of residues to repack (the threaded residues and neighbors)
		(*scorefxn_)(pose); //cause scoring got messed up somehow
		utility::vector1< bool > to_repack( pose.size(), false );
		for ( core::Size i = 1; i <= thread_position.size(); i++ ) {
			for ( core::Size j = thread_position[i]; j < thread_position[i] + sequences_[i].size(); j++ ) {
				to_repack[j] = true;
			}
		}
		if ( pack_neighbors_ ) {
			core::select::fill_neighbor_residues(pose, to_repack, 6.0);
		}
		TR << "Created vector of repackable residues" << std::endl;
		//Set up tasks to turn off design, and repack threaded+neighboring residues
		core::pack::task::PackerTaskOP task = task_factory_->create_packer_task( pose ); //want to create new task object each iteration and each packer run
		for ( core::Size resnum = 1; resnum <= pose.size(); resnum++ ) {
			if ( !to_repack[resnum] ) {
				turn_off_packing->include_residue( resnum );
			}
		}
		turn_off_packing->apply(pose, *task);
		turn_off_design->apply(pose, *task);

		//Set up packer
		protocols::minimization_packing::PackRotamersMoverOP packer = protocols::minimization_packing::PackRotamersMoverOP( new protocols::minimization_packing::PackRotamersMover(scorefxn_, task, PACKER_ITER) );

		TR << "Threaded and tasks set up" << std::endl;

		//Hardmin, Hardpack
		min->apply( pose );
		TR << "Done minimizing" << std::endl;
		packer->apply( pose );
		TR << "Done hardpacking" << std::endl;
	}
	( *scorefxn_ )( pose );
}

core::pose::PoseOP
GraftSwitchMover::get_additional_output() {
	core::pose::PoseOP out_pose;
	if ( output_list_.size() == 0 ) {
		return nullptr;
	} else if ( output_list_.size() == 1 ) {
		out_pose = orig_pose_;
	} else {
		out_pose = utility::pointer::make_shared< core::pose::Pose >( *orig_pose_ );
	}
	utility::vector1<core::Size> thread_position = output_list_.front();
	output_list_.pop_front();

	TR << "Outputting thread with ";
	for ( core::Size i=1; i <= thread_position.size(); i++ ) {
		TR << sequences_[i] << " at " << thread_position[i] << ", ";
	}
	TR << std::endl;

	//Thread
	protocols::simple_moves::SimpleThreadingMoverOP threader = protocols::simple_moves::SimpleThreadingMoverOP( new protocols::simple_moves::SimpleThreadingMover() );
	threader->set_pack_neighbors(pack_neighbors_);
	for ( core::Size i=1; i <= thread_position.size(); i++ ) {
		threader->set_sequence(sequences_[i], thread_position[i]);
		threader->apply( *out_pose );
	}

	if ( pack_min_ ) {
		//TaskOps
		core::pack::task::operation::PreventRepackingOP turn_off_packing = core::pack::task::operation::PreventRepackingOP( new core::pack::task::operation::PreventRepacking() );
		core::pack::task::operation::RestrictToRepackingOP turn_off_design = core::pack::task::operation::RestrictToRepackingOP( new core::pack::task::operation::RestrictToRepacking() );

		//Set up Minimizer
		core::kinematics::MoveMapOP mm = core::kinematics::MoveMapOP( new core::kinematics::MoveMap() );
		mm->set_bb( true );
		mm->set_chi( true );
		protocols::minimization_packing::MinMoverOP min = protocols::minimization_packing::MinMoverOP( new protocols::minimization_packing::MinMover( mm, scorefxn_, "lbfgs_armijo_nonmonotone", 0.0001, true ) );

		//Create vector of residues to repack (the threaded residues and neighbors
		(*scorefxn_)(*out_pose); //cause scoring got messed up somehow
		utility::vector1< bool > to_repack( out_pose->size(), false );
		for ( core::Size i = 1; i <= thread_position.size(); i++ ) {
			for ( core::Size j = thread_position[i]; j < thread_position[i] + sequences_[i].size(); j++ ) {
				to_repack[j] = true;
			}
		}
		if ( pack_neighbors_ ) {
			core::select::fill_neighbor_residues(*out_pose, to_repack, 6.0);
		}

		//Set up tasks to turn off design, and repack threaded+neighboring residues
		core::pack::task::PackerTaskOP task = task_factory_->create_packer_task( *out_pose ); //want to create new task object each iteration and each packer run
		for ( core::Size resnum = 1; resnum <= out_pose->size(); resnum++ ) {
			if ( !to_repack[resnum] ) {
				turn_off_packing->include_residue( resnum );
			}
		}
		turn_off_packing->apply(*out_pose, *task);
		turn_off_design->apply(*out_pose, *task);

		//Set up packer
		protocols::minimization_packing::PackRotamersMoverOP packer = protocols::minimization_packing::PackRotamersMoverOP( new protocols::minimization_packing::PackRotamersMover(scorefxn_, task, PACKER_ITER) );

		TR << "Threaded and tasks set up" << std::endl;

		//Hardmin, Hardpack
		min->apply( *out_pose );
		TR << "Done minimizing" << std::endl;
		packer->apply( *out_pose );
		TR << "Done hardpacking" << std::endl;
	}
	( *scorefxn_ )( *out_pose );

	return out_pose;
}

std::string
GraftSwitchMover::get_name() const { return mover_name(); }

std::string
GraftSwitchMover::mover_name() {
	return "GraftSwitchMover";
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
GraftSwitchMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	Pose const & pose
)
{
	if ( tag->hasOption("selector") ) {
		if ( tag->hasOption("start" ) || tag->hasOption("n_term") ) {
			TR << "Cannot use both start/end or n_term and selector options; start_selector will be used" << std::endl;
		}
		std::string const selector_name ( tag->getOption< std::string >( "selector" ) );
		if ( TR.visible() ) TR << "Set selector name to " << selector_name << "." << std::endl;
		core::select::residue_selector::ResidueSelectorCOP selector;
		try {
			selector = datamap.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name );
		} catch ( utility::excn::Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from AddCompositionConstraintMover::parse_tag()\n" + e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_message );
		}
		runtime_assert( selector );
		threadable_residues_selector_ = selector->clone();
	} else if ( tag->hasOption("start") && tag->hasOption("end") ) {
		if ( tag->hasOption("n_term") ) {
			TR << "Cannot use both start/end and n_term; start/end residues will be used" << std::endl;
		}
		start_graft_ = tag->getOption<core::Size>("start", 0);
		end_graft_ = tag->getOption<core::Size>("end", 0);
	} else {
		n_terminus_ = tag->getOption< bool >("n_term", false);
		graft_loop_ = tag->getOption< bool >("graft_on_latch_loop", true);
		find_helix_ = true;
	}

	pack_neighbors_ = tag->getOption< bool >("pack_neighbors", true);
	seq_in_any_order_ = tag->getOption< bool >("any_order", true);
	pack_min_ = tag->getOption< bool >("pack_min", true);
	burial_cutoff_ = tag->getOption< core::Size >("burial_cutoff", 6);
	sequences_ = utility::string_split( tag->getOption<std::string>("sequence", ""), ',' );

	std::string important_residues_string( tag->getOption<std::string>("important_residues", "") );
	utility::vector1< std::string > important_residues_string_vector = utility::string_split( important_residues_string, '/' );
	for ( core::Size i = 1; i <= important_residues_string_vector.size(); i++ ) {
		utility::vector1< std::string > subsplit = utility::string_split( important_residues_string_vector[i], ',' );
		utility::vector1< core::Size > subsplit_Size;
		for ( core::Size j = 1; j <= subsplit.size(); j++ ) {
			subsplit_Size.push_back( utility::string2Size( subsplit[j] ) );
		}
		important_residues_.push_back( subsplit_Size );
	}

	parse_score_function( tag, datamap, filters, movers, pose );
	parse_task_operations( tag, datamap, filters, movers, pose );
}

/// @brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
void
GraftSwitchMover::parse_score_function(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	ScoreFunctionOP new_score_function( protocols::rosetta_scripts::parse_score_function( tag, datamap ) );
	if ( new_score_function == nullptr ) return;
	score_function( new_score_function );
}

/// @brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
void
GraftSwitchMover::parse_task_operations(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == nullptr ) return;
	task_factory( new_task_factory );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
GraftSwitchMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new GraftSwitchMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
GraftSwitchMover::clone() const
{
	return protocols::moves::MoverOP( new protocols::switches::GraftSwitchMover( *this ) );
}

// setters

void GraftSwitchMover::score_function( ScoreFunctionCOP sf )
{
	runtime_assert( sf != nullptr );
	scorefxn_ = sf;
}

void GraftSwitchMover::task_factory( TaskFactoryCOP tf )
{
	runtime_assert( tf != nullptr );
	task_factory_ = tf;
}


/// @brief finds the residue number that starts the C-terminal helix
/// @details returns a single residue number as an integer
core::Size
GraftSwitchMover::get_c_terminal_helix(
	core::pose::PoseOP const& pose
){
	core::scoring::dssp::Dssp dssp_struct = core::scoring::dssp::Dssp( *pose );
	std::string secstruct = dssp_struct.get_dssp_secstruct();
	TR << "Secondary Structure: " << secstruct << std::endl;

	//TODO: Add option to include the loop in the graftable region (by default it is)
	bool last_loop = false;
	bool seen_helix = false;
	for ( core::Size i = secstruct.length()-1; i != 0; i-- ) {
		if ( last_loop && secstruct[i-1] == 'H' ) {
			return i;
		}

		char structure = secstruct[i];
		TR << "residue number " << i+1 << " with secondary structure " << structure << std::endl;
		if ( structure == 'H' ) {
			seen_helix = true;
			continue;
		} else if ( !seen_helix && !last_loop ) {
			continue;
		} else if ( seen_helix && structure == 'L' ) {
			last_loop = true;
			if ( !graft_loop_ ) {
				return i;
			}
		}
	}
	return 0;
}

/// @brief finds the residue number that starts the C-terminal helix
/// @details returns a single residue number as an integer
core::Size
GraftSwitchMover::get_n_terminal_helix(
	core::pose::PoseOP const& pose
){
	core::scoring::dssp::Dssp dssp_struct = core::scoring::dssp::Dssp( *pose );
	std::string secstruct = dssp_struct.get_dssp_secstruct();
	TR << "Secondary Structure: " << secstruct << std::endl;

	//TODO: Add option to include the loop in the graftable region (by default it is)
	bool first_loop = false;
	bool seen_helix = false;
	for ( core::Size i = 1; i < pose->size(); i++ ) {
		if ( first_loop && secstruct[i+1] == 'H' ) {
			return i;
		}
		char structure = secstruct[i];
		TR << "residue number " << i << " with secondary structure " << structure << std::endl;
		if ( structure == 'H' ) {
			seen_helix = true;
			continue;
		} else if ( !seen_helix && !first_loop ) {
			continue;
		} else if ( seen_helix && structure == 'L' ) {
			first_loop = true;
			if ( !graft_loop_ ) {
				return i;
			}
		}
	}
	return 0;
}

//initial filter - if key residues aren't buried at all don't bother repacking
/// @brief inital filter to determine if the key residues that need to be caged, are.
/// @details uses NeighborsByDistanceCalculator as metric for burial
bool
GraftSwitchMover::init_burial_filter(
	core::pose::Pose const& pose,
	core::Size index,
	core::Size thread_position
) {
	if ( important_residues_.size() < 1 ) return true;

	core::Real running_sum = 0;
	core::Size num_residues = 0;
	basic::MetricValue<core::Size> num_neighbors;
	for ( core::Size i=1; i <= important_residues_[index].size(); i++ ) {
		core::Size residue = important_residues_[index][i] + thread_position - 1;
		if ( residue <= pose.size() ) {
			protocols::pose_metric_calculators::NeighborsByDistanceCalculator count_neighbors(residue, 6.0);
			count_neighbors.get( "num_neighbors", num_neighbors, pose );
			TR << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << std::endl;
			TR << "Residue " << residue << " has " << num_neighbors.print() << " neighbors" << std::endl;
			TR << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << std::endl;
			running_sum = running_sum + core::Real( num_neighbors.value() );
			num_residues++;
		}
	}
	core::Real avg_degree = running_sum / num_residues;
	TR << "INIT avg degree is: " << avg_degree << std::endl;
	return avg_degree >= burial_cutoff_;
}

/// @brief returns all tuples of size n representing starting positions for each sequences
/// @details Creates a list of vectors (using a list for speed in deletion phase) where each vector is of length=number of sequences.  Each number in the vector corresponds to a starting position for threading that will fit the corresponding sequence.
std::list< utility::vector1< core::Size > >
GraftSwitchMover::all_threading_start_combinations() {
	std::list< utility::vector1< core::Size > > final_starting_positions;

	//Setup list of positions the first sequence can be in
	utility::vector1< core::Size > list_of_positions;
	for ( std::set<core::Size>::iterator j = threadable_residues_.begin(); j != threadable_residues_.end(); j++ ) {
		bool add = true;
		for ( core::Size k = *j; k < *j + sequences_[1].size(); k++ ) {
			if ( threadable_residues_.find( k ) == threadable_residues_.end()  ) add = false;
		}
		if ( add ) {
			list_of_positions.push_back( *j );
			final_starting_positions.push_back( list_of_positions );
			list_of_positions.clear();
		}
	}

	//iterate through the other sequences, delete positions where they can't fit, multiply items in the list for each position each sequence fits in
	//this won't fire if >1 sequence
	std::set< core::Size > working_threadable_residues = threadable_residues_;
	for ( core::Size i=2; i <= sequences_.size(); i++ ) {
		std::list< utility::vector1< core::Size > > working_starting_positions;
		for ( std::list< utility::vector1< core::Size > >::iterator j = final_starting_positions.begin(); j != final_starting_positions.end(); j++ ) {
			//go through all previously looked at seq at the positions dictated in *j and remove from working_threadable_residues
			//because given those starting positions, the residues that the sequences make up are taken and not designable
			for ( core::Size k=1; k < i; k++ ) {
				core::Size start_pos = (*j)[k];
				core::Size len_seq = sequences_[k].size();
				for ( core::Size l = start_pos; l < start_pos + len_seq; l++ ) {
					working_threadable_residues.erase(l);
				}
			}

			//looking at what's left, iterate through working_threadable_residues and add to copy of j if position works.  push j to end of working_starting_pos
			std::set<core::Size>::iterator k = working_threadable_residues.begin();
			if ( !seq_in_any_order_ ) k = working_threadable_residues.upper_bound( j->back() + sequences_[i-1].size() - 1 );
			while ( k != working_threadable_residues.end() ) {
				bool add = true;
				for ( core::Size l = *k; l < *k + sequences_[i].length(); l++ ) {
					if ( working_threadable_residues.find( l ) == working_threadable_residues.end()  ) add = false;
				}
				if ( add ) {
					utility::vector1< core::Size > copy_j = *j;
					copy_j.push_back(*k);
					working_starting_positions.push_back( copy_j );
				}
				k++;
			}
			working_threadable_residues = threadable_residues_;
		}
		final_starting_positions = working_starting_positions;
	}

	return final_starting_positions;
}

} // moves
} // protocols

