// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/SimpleThreadingMover.cc
/// @brief Very Simple class for threading a regional sequence onto a structure
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @modified NCAA support added by Vikram K. Mulligan (vmulligan@flatironinstitute.org

#include <protocols/simple_moves/SimpleThreadingMoverCreator.hh>
#include <protocols/simple_moves/SimpleThreadingMover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueType.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/util.hh>
#include <core/simple_metrics/metrics/SequenceMetric.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>

#include <basic/Tracer.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/tag/util.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <core/pack/task/PackerTask.hh> // AUTO IWYU For PackerTask
#include <core/pack/task/ResidueLevelTask.hh> // AUTO IWYU For ResidueLevelTask

static basic::Tracer TR("protocols.simple_moves.SimpleThreadingMover");

namespace protocols {
namespace simple_moves {
using namespace core::select;

/// @brief Initialization constructor.
SimpleThreadingMover::SimpleThreadingMover(std::string thread_sequence, core::Size start_position):
	protocols::moves::Mover("SimpleThreadingMover"),
	start_position_(start_position),
	thread_sequence_( thread_sequence )
{}

/// @brief Destructor.
SimpleThreadingMover::~SimpleThreadingMover()= default;


void
SimpleThreadingMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
)
{
	pack_neighbors_ = tag->getOption< bool >("pack_neighbors", pack_neighbors_);
	neighbor_dis_ = tag->getOption< core::Real >("neighbor_dis", neighbor_dis_);

	parsed_position_ = tag->getOption< std::string >("start_position");

	thread_sequence_ = tag->getOption< std::string >("thread_sequence", thread_sequence_);
	if ( tag->hasOption("scorefxn") ) {
		scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );
	}

	skip_unknown_mutant_ = tag->getOption< bool >("skip_unknown_mutant", skip_unknown_mutant_);
	pack_rounds_ = tag->getOption< core::Size >("pack_rounds", pack_rounds_);

	if ( tag->hasOption( "sequence_mode" ) ) {
		set_sequence_mode( tag->getOption<std::string>( "sequence_mode" ) );
	}
}


protocols::moves::MoverOP
SimpleThreadingMover::clone() const{
	return utility::pointer::make_shared< SimpleThreadingMover >(*this);

}

//SimpleThreadingMover & operator=( SimpleThreadingMover const & src){
// return SimpleThreadingMover(src);

//}

moves::MoverOP
SimpleThreadingMover::fresh_instance() const {
	return utility::pointer::make_shared< SimpleThreadingMover >();

}

/// @brief Provide the citation.
void
SimpleThreadingMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	using namespace basic::citation_manager;

	citations.add(
		utility::pointer::make_shared< UnpublishedModuleInfo >(
		mover_name(),
		CitedModuleType::Mover,
		"Jared Adolf-Bryfogle",
		"The Scripps Research Institute, La Jolla, CA",
		"jadolfbr@gmail.com"
		)
	);
}

void
SimpleThreadingMover::set_sequence(std::string thread_sequence, core::Size start_position){
	start_position_ = start_position;
	thread_sequence_ = thread_sequence;

}

void
SimpleThreadingMover::set_pack_neighbors(bool pack_neighbors){
	pack_neighbors_ = pack_neighbors;
}

void
SimpleThreadingMover::set_pack_rounds(core::Size pack_rounds){
	pack_rounds_ = pack_rounds;
}

void
SimpleThreadingMover::set_neighbor_distance(core::Real neighbor_dis){
	neighbor_dis_ = neighbor_dis;
}

/// @brief Set the sequence mode, by string.
/// @details This determines how the sequence is interpreted (one-letter codes, three-letter codes, etc.).
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
SimpleThreadingMover::set_sequence_mode (
	std::string const & mode_string_in
) {
	using namespace core::simple_metrics::metrics;
	SequenceMetricMode const mode( SequenceMetric::mode_enum_from_name( mode_string_in ) );
	runtime_assert_string_msg( mode != SequenceMetricMode::INVALID_MODE, "Error in SimpleThreadingMover::set_sequence_mode(): Could not parse \"" + mode_string_in + "\" as a valid sequence mode.  Valid modes are: " + SequenceMetric::allowed_output_modes() + "." );
	set_sequence_mode( mode );
}

/// @brief Set the sequence mode, by enum.
/// @details This determines how the sequence is interpreted (one-letter codes, three-letter codes, etc.).
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
SimpleThreadingMover::set_sequence_mode (
	core::simple_metrics::metrics::SequenceMetricMode const mode_in
) {
	runtime_assert_string_msg(
		static_cast< core::Size >(mode_in) > 0 && mode_in < core::simple_metrics::metrics::SequenceMetricMode::END_OF_LIST,
		"Error in SimpleThreadingMover::set_sequence_mode(): Mode not recognized."
	);
	sequence_mode_ = mode_in;
}

void
SimpleThreadingMover::set_scorefxn(core::scoring::ScoreFunctionCOP scorefxn){
	scorefxn_ = scorefxn;
}

bool
SimpleThreadingMover::get_pack_neighbors() const{
	return pack_neighbors_;

}

core::Real
SimpleThreadingMover::get_neighbor_distance() const {
	return neighbor_dis_;

}

void
SimpleThreadingMover::apply(
	core::pose::Pose& pose
) {
	using namespace core::pack::task;
	using namespace core::scoring;
	using namespace core::pack::task::operation;
	using namespace protocols::minimization_packing;
	using namespace core::simple_metrics::metrics;
	using namespace core::select::residue_selector;

	debug_assert( static_cast< core::Size >( sequence_mode_ ) > 0 && sequence_mode_ < SequenceMetricMode::END_OF_LIST ); //Should be true.

	if ( scorefxn_ == nullptr ) {
		scorefxn_ = core::scoring::get_score_function();
	}

	//This could have just as easily have been a task op.

	if ( thread_sequence_ == "0" ) {
		utility_exit_with_message("Sequence not set for threading.  Cannot continue.");
	}

	// Changed by VKM on 17 Sept 2019: Since we're using the MutateResidueMover under the hood, we can skip packing if necessary.
	if ( pack_rounds_ < 1 ) {
		TR.Warning << "Warning: pack_rounds is set to 0.  Mutations will be made, but sidechain geometry will likely be poor.  Subsequent repacking is highly recommended." << std::endl;
	}

	if ( parsed_position_ !=  "NA" ) {
		start_position_ = core::pose::parse_resnum(parsed_position_, pose);
	}

	TR << "Threading Sequence :" << thread_sequence_ << ":" << std::endl;

	//Thread the sequence:
	runtime_assert_string_msg( start_position_ > 0 && start_position_ <= pose.total_residue(), "Error in SimpleThreadingMover::apply(): The start position, " + std::to_string(start_position_) + ", is outside of the range of the " + std::to_string( pose.total_residue() ) + "-residue pose." );
	std::map< core::Size, std::string > const mutations = determine_mutations( start_position_,  utility::strip_whitespace( thread_sequence_ ), sequence_mode_, pose );
	ResidueIndexSelectorOP select_mutated_residues( utility::pointer::make_shared< ResidueIndexSelector >() );

	for ( std::map< core::Size, std::string >::const_iterator it( mutations.begin() ); it!=mutations.end(); ++it ) {
		if ( it->first > pose.total_residue() || it->first < 1 ) {
			TR.Warning << "Position " << it->first << " is not in the " << pose.total_residue() << "-residue pose.  Cannot mutate to " << it->second << "." << std::endl;
			continue;
		}
		TR << "Mutating position " << it->first << " to " << it->second << "." << std::endl;
		MutateResidue mutres( it->first, it->second );
		mutres.apply( pose );

		//Store this residue for later repacking.
		select_mutated_residues->append_index( it->first );
	}

	//Repack without design if we're repacking:
	if ( pack_rounds_ > 0 ) {
		TaskFactoryOP tf = utility::pointer::make_shared< TaskFactory >();
		tf->push_back( utility::pointer::make_shared< InitializeFromCommandline >() );
		tf->push_back( utility::pointer::make_shared<core::pack::task::operation::RestrictToRepacking>() );
		if ( pack_neighbors_ ) {
			NeighborhoodResidueSelectorOP select_neighborhood( utility::pointer::make_shared< NeighborhoodResidueSelector >( select_mutated_residues, neighbor_dis_ ) );
			NotResidueSelectorOP select_not_neighborhood( utility::pointer::make_shared< NotResidueSelector >() );
			select_not_neighborhood->set_residue_selector( select_neighborhood );
			OperateOnResidueSubsetOP prevent_repacking_not_neighborhood( utility::pointer::make_shared< OperateOnResidueSubset >( utility::pointer::make_shared< PreventRepackingRLT >(), select_not_neighborhood ) );
			tf->push_back( prevent_repacking_not_neighborhood );
		} else {
			NotResidueSelectorOP select_not_mutated( utility::pointer::make_shared< NotResidueSelector >() );
			select_not_mutated->set_residue_selector( select_mutated_residues );
			OperateOnResidueSubsetOP prevent_repacking_not_mutated( utility::pointer::make_shared< OperateOnResidueSubset >( utility::pointer::make_shared< PreventRepackingRLT >(), select_not_mutated ) );
			tf->push_back( prevent_repacking_not_mutated );
		}
		PackerTaskOP task = tf->create_task_and_apply_taskoperations(pose);

		scorefxn_->score(pose); //Segfault Protection.

		protocols::minimization_packing::PackRotamersMoverOP packer( utility::pointer::make_shared< protocols::minimization_packing::PackRotamersMover >( scorefxn_, task, pack_rounds_ ) );

		packer->apply( pose );
	}

	TR << "Complete" <<std::endl;

}


////////////// Creator /////////

std::string SimpleThreadingMover::get_name() const {
	return mover_name();
}

std::string SimpleThreadingMover::mover_name() {
	return "SimpleThreadingMover";
}

void SimpleThreadingMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	using namespace core::simple_metrics::metrics;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default(
		"pack_neighbors", xsct_rosetta_bool,
		"Option to pack neighbors while threading.  By default, only the mutated residues, and not the neighbors, are repacked.", "false");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"neighbor_dis", xsct_real,
		"Distance to repack neighbor side chains. Repack shell distance for each threaded residue.  Default 6.0 Angstroms.",
		"6.0");
	attlist + XMLSchemaAttribute::required_attribute(
		"start_position", xs_string,
		"Position to start thread. PDB numbering (like 30L) or Rosetta pose numbering. "
		"PDB numbering parsed at apply time to allow for pose-length changes prior to apply of this mover");
	attlist + XMLSchemaAttribute::required_attribute(
		"thread_sequence", xs_string,
		"The residue sequence that we will be grafting.  This can be provided as one-letter codes (e.g. \"RSTX[DASP]LNE\", comma-separated three-letter codes (e.g. \"ARG,SER,THR,DAS,LEU,ASN,GLU\"), base-names (e.g. \"ARG,SER,THR,DASP,LEU,ASN,GLU\"), or full names (e.g. \"ARG,SER:N_Methylation,THR,DASP,LEU,ASN,GLU\"), depending on the setting for sequence_mode." );
	attlist + XMLSchemaAttribute(
		"scorefxn", xs_string,
		"Optional Scorefunction name passed - setup in score function block");
	attlist + XMLSchemaAttribute(
		"skip_unknown_mutant", xsct_rosetta_bool,
		"Skip unknown amino acids in thread_sequence string instead of throwing an exception." );
	attlist + XMLSchemaAttribute::attribute_w_default(
		"pack_rounds", xsct_positive_integer,
		"Number of packing rounds for threading.  Set this to 0 to skip all packing, in which case the new side-chains will likely be in very poor conformations.  Defaults to 5.", "5");

	utility::vector1< std::string > const allowed_modes_vec( SequenceMetric::allowed_output_modes_as_vector() );
	utility::tag::add_schema_restrictions_for_strings( xsd, "SimpleThreadingMover_input_modes", allowed_modes_vec );
	attlist + XMLSchemaAttribute::attribute_w_default(
		"sequence_mode", "SequenceMetric_output_modes", "The format for the input sequence.  Allowed output formats are: " + SequenceMetric::allowed_output_modes() + ".", SequenceMetric::mode_name_from_enum( SequenceMetricMode::ONELETTER_CODE ) );

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)\n"
		"Modified: Vikram K. Mulligan (vmulligan@flatironinstitute.org) to add support for noncanonicals.\n"
		"This mover functions to thread the sequence of a region onto the given pose. "
		"Nothing fancy here. Useful when combined with -parser:string_vars option "
		"to replace strings within the RosettaScript. "
		"For more a more fancy comparative modeling protocol, please see the lovely RosettaCM",
		attlist );
}

/// @brief Given the sequence, the interpretation mode, and the start position, fill a map of position->mutation name.
std::map< core::Size, std::string >
SimpleThreadingMover::determine_mutations(
	core::Size const start_position,
	std::string const & sequence,
	core::simple_metrics::metrics::SequenceMetricMode const mode,
	core::pose::Pose const & pose
) const {
	using namespace core::simple_metrics::metrics;

	switch( mode ) {
	case SequenceMetricMode::ONELETTER_CODE :
		return determine_mutations_oneletter( start_position, sequence );
	case SequenceMetricMode::THREELETTER_CODE:
	case SequenceMetricMode::BASE_NAME :
	case SequenceMetricMode::FULL_NAME :
		return determine_mutations_comma_separated( start_position, sequence, mode, pose );
	default :
		utility_exit_with_message( "Program error in SimpleThreadingMover::determine_mutations(): The sequence input mode was improperly set!  Please consult a developer." );
	}
	// To keep compiler happy:
	return std::map< core::Size, std::string >();
}

/// @brief Given the as one-letter codes and the start position, fill a map of position->mutation name.
/// @details Supports possibility of positions of the form X[NCAA_NAME].
std::map< core::Size, std::string >
SimpleThreadingMover::determine_mutations_oneletter(
	core::Size const start_position,
	std::string const & sequence
) const {
	std::map< core::Size, std::string > outmap;
	core::Size curpos( start_position );
	for ( core::Size i(0), imax(sequence.length()); i<imax; ++i ) { //Loop through the string; strings are zero-based
		if ( core::chemical::oneletter_code_specifies_aa( sequence[i] ) && sequence[i] != 'X' && sequence[i] != 'Z'  && sequence[i] != 'z' && sequence[i] != 'w' ) {
			debug_assert( outmap.count(curpos) == 0 );
			outmap[curpos] = core::chemical::name_from_aa( core::chemical::aa_from_oneletter_code( sequence[i] ) );
			++curpos;
		} else if ( sequence[i] == 'X' ) {
			runtime_assert_string_msg( i < imax-1 && sequence[i+1] == '[', "Error in SimpleThreadingMover::determine_mutations_oneletter():  The character \"X\" must be followed by a square bracket (\"[\")." );
			core::Size j( i+1 );
			bool found(false);
			bool bracket_contents(false);
			for ( ; j<imax; ++j ) {
				if ( sequence[j] == ']' ) {
					found=true;
					break;
				} else if ( sequence[j] != '[' ) {
					bracket_contents = true;
				}
			}
			runtime_assert_string_msg(found, "Error in SimpleThreadingMover::determine_mutations_oneletter():  An opening bracket was found in the sequence with no closing bracket.");
			runtime_assert_string_msg(bracket_contents, "Error in SimpleThreadingMover::determine_mutations_oneletter():  Brackets were found in the sequence that enclose nothing.");
			debug_assert( outmap.count(curpos) == 0 );
			outmap[curpos] = sequence.substr( i+2, j-i-2 );
			++curpos;
			i=j;
		} else if ( sequence[i] == '-' ) {
			//Skip dashes in the sequence.
			++curpos;
		} else {
			if ( skip_unknown_mutant_ ) {
				TR.Warning << "Could not parse one-letter code \"" << sequence.substr(i, 1) << "\".  Skipping." << std::endl;
				++curpos;
			} else {
				utility_exit_with_message( "Error in SimpleThreadingMover::determine_mutations_oneletter():  Could not parse character \"" + sequence.substr(i,1) + "\" in sequence \"" + sequence + "\"." );
			}
		}
	}
	return outmap;
}

/// @brief Given the sequence as a comma-separated list of either three-letter codes, base names, or full names, plus the start position,
/// fill a map of position->mutation name.
std::map< core::Size, std::string >
SimpleThreadingMover::determine_mutations_comma_separated(
	core::Size const start_position,
	std::string const & sequence,
	core::simple_metrics::metrics::SequenceMetricMode const mode,
	core::pose::Pose const & pose
) const {
	static const std::string errmsg( "Error in SimpleThreadingMover::determine_mutations_comma_separated():  " );

	//Split the sequence by commas:
	utility::vector1< std::string > const splitseq( utility::string_split_multi_delim( sequence, " \t\n," ) );
	runtime_assert_string_msg( !splitseq.empty(), errmsg + "The sequence is empty, or could otherwise not be parsed." );

	std::map< core::Size, std::string > outmap;
	core::Size curpos( start_position );

	for ( core::Size i(1), imax(splitseq.size()); i<=imax; ++i ) { //Loop through the split sequence
		if ( splitseq[i][0] == '-' ) {
			//If a position is a dash or starts with a dash, skip it.
			++curpos;
			continue;
		}
		if ( mode == core::simple_metrics::metrics::SequenceMetricMode::THREELETTER_CODE ) {
			core::chemical::ResidueTypeFinder finder( *(pose.residue_type_set_for_pose()) );
			core::chemical::ResidueTypeCOP restype( finder.name3( splitseq[i] ).get_representative_type() );
			if ( skip_unknown_mutant_ ) {
				if ( restype == nullptr ) {
					TR.Warning << "Could not parse amino acid with three-letter code \"" + splitseq[i] + "\".  Skipping and continuing on." << std::endl;
					++curpos;
					continue;
				}
			} else {
				runtime_assert_string_msg( restype != nullptr, errmsg + "Could not find a suitable residue type for three-letter code \"" + splitseq[i] + "\"." );
			}
			debug_assert( outmap.count(curpos) == 0 );
			outmap[curpos] = restype->base_name();
			++curpos;
		} else { //Base name or full name
			debug_assert( outmap.count(curpos) == 0 );

			core::chemical::ResidueTypeCOP restype;
			if ( mode == core::simple_metrics::metrics::SequenceMetricMode::BASE_NAME ) {
				restype = core::chemical::ResidueTypeFinder( *pose.residue_type_set_for_pose() ).residue_base_name( splitseq[i] ).get_representative_type();
			} else {
				core::chemical::ResidueTypeFinder finder ( *pose.residue_type_set_for_pose() );
				utility::vector1< std::string > namesplit( utility::string_split( splitseq[i], ':' ) );
				core::Size const nentries( namesplit.size() );
				debug_assert( nentries > 0 ); // Should be true
				finder.residue_base_name( namesplit[1] );
				if ( nentries > 1 ) {
					namesplit.erase( namesplit.begin() ); //Delete the first entry.
					finder.patch_names( namesplit );
				}
				restype = finder.get_representative_type();
			}

			if ( skip_unknown_mutant_ ) {
				if ( restype == nullptr ) {
					TR.Warning << "Could not find suitable residue type for name \"" << splitseq[i] << "\".  Skipping and continuing on." << std::endl;
					++curpos;
					continue;
				}
			} else {
				runtime_assert_string_msg( restype != nullptr, errmsg + "Could not find a suitable residue type for name \"" + splitseq[i] + "\"." );
			}

			outmap[curpos] = splitseq[i];
			++curpos;
		}
	}

	return outmap;
}

std::string SimpleThreadingMoverCreator::keyname() const {
	return SimpleThreadingMover::mover_name();
}

protocols::moves::MoverOP
SimpleThreadingMoverCreator::create_mover() const {
	return utility::pointer::make_shared< SimpleThreadingMover >();
}

void SimpleThreadingMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SimpleThreadingMover::provide_xml_schema( xsd );
}




}//simple_moves
}//protocols
