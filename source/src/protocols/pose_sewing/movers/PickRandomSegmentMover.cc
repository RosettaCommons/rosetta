// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/movers/PickRandomSegmentMover.cc
/// @brief replaces pose with a random starting segment from a segment file
/// @author frankdt (frankdt@email.unc.edu)

// Unit headers
#include <protocols/pose_sewing/movers/PickRandomSegmentMover.hh>
#include <protocols/pose_sewing/movers/PickRandomSegmentMoverCreator.hh>
#include <protocols/pose_sewing/util.hh>
#include <protocols/moves/mover_schemas.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/filters/FilterFactory.hh>
#include <protocols/rosetta_scripts/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <cmath>

static basic::Tracer TR( "protocols.pose_sewing.movers.PickRandomSegmentMover" );

namespace protocols {
namespace pose_sewing {
namespace movers {

using namespace protocols::pose_sewing::data_storage;
using namespace protocols::filters;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
PickRandomSegmentMover::PickRandomSegmentMover():
	protocols::moves::Mover( PickRandomSegmentMover::mover_name() )
{
	pose_vector_ = data_storage::TerminalDSSPSortedPoseVectorOP(new data_storage::TerminalDSSPSortedPoseVector());

	set_ss_defaults();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
PickRandomSegmentMover::PickRandomSegmentMover( PickRandomSegmentMover const &  ) = default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PickRandomSegmentMover::~PickRandomSegmentMover(){}

std::map<char, SeedSettings >
PickRandomSegmentMover::get_ss_settings() const {
	return ss_settings_;
}

SeedSettings
PickRandomSegmentMover::create_ss_defaults() const {
	SeedSettings def_settings;

	def_settings.min_terminal_length_ = 1;
	def_settings.max_terminal_length_ = 100;

	return def_settings;
}

void
PickRandomSegmentMover::set_ss_defaults(){
	ss_settings_.clear();
	SeedSettings def_settings = create_ss_defaults();
	ss_settings_['A'] = def_settings;
}

SeedSettings
PickRandomSegmentMover::get_ss_defaults() const {
	if ( ss_settings_.count('A') ) {
		return ss_settings_.at('A');
	} else {
		return create_ss_defaults();
	}
}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
PickRandomSegmentMover::apply( core::pose::Pose& pose){

	utility::vector1<PoseWithTerminalSegmentsOfKnownDSSPOP> seed_vector;

	core::Size max_attempts = 1000;//Here to fail more quickly if it is to happen.

	if ( even_sampling_ ) {
		core::Size index = numeric::random::rg().random_range(1, segment_file_paths_.size());
		seed_vector = pose_vector_->populate_from_segment_file_and_get_random_vector_set(segment_file_paths_[index], "EHL","EHL",max_attempts);
	} else {
		seed_vector = pose_vector_->populate_from_segment_file_and_get_random_vector_set(segment_file_paths_, "EHL","EHL",max_attempts);
	}


	core::Size attempts = 0;
	PoseWithTerminalSegmentsOfKnownDSSPOP selected_pose = nullptr;
	core::pose::PoseOP seed_pose = nullptr;
	for ( core::Size i = 1; i <= seed_vector.size(); ++i ) {
		attempts+=1;
		bool passed_filters = true;

		char const c_term_ss = seed_vector[i]->get_C_term_DSSP();
		char const n_term_ss = seed_vector[i]->get_N_term_DSSP();
		char c_term_ss_set = 'A';
		char n_term_ss_set = 'A';
		if ( ss_settings_.count( c_term_ss) ) {
			c_term_ss_set = c_term_ss;
		}
		if ( ss_settings_.count( n_term_ss) ) {
			n_term_ss_set = n_term_ss;
		}

		core::Size const c_term_len = seed_vector[i]->get_C_term_length();
		core::Size const n_term_len = seed_vector[i]->get_N_term_length();

		TR << "Attempts: " << attempts << std::endl;
		if ( c_term_len < ss_settings_.at(c_term_ss_set).min_terminal_length_ ||
				c_term_len > ss_settings_.at(c_term_ss_set).max_terminal_length_ ||
				n_term_len < ss_settings_.at(n_term_ss_set).min_terminal_length_ ||
				n_term_len > ss_settings_.at(n_term_ss_set).max_terminal_length_ ) {
			TR << "Seed did not match terminal lengths." << std::endl;
			continue;
		} else if ( attempts > max_attempts ) {
			TR <<"Max attempts hit" << std::endl;
			break;
		} else {
			if ( c_term_ss == 'E' && n_term_ss == 'E' && match_terminal_sheet_lengths_ ) {
				int const diff = c_term_len - n_term_len;
				core::Size const offset = std::abs(diff);
				if ( offset <= terminal_sheet_length_max_delta_ ) {
					selected_pose = seed_vector[i];
				}
			} else {
				selected_pose = seed_vector[i];
			}

			if ( selected_pose && (filters_.size() || seg_filters_.count(selected_pose->get_segfile_path())) ) {
				seed_pose = selected_pose->get_source_pose_op(false /*copy pose*/);

				//Create MT, update from elements if there.

				//Doing this in two steps to avoid copying the filters into a full vector.
				for ( FilterCOP filter : filters_ ) {
					TR << "Running filter "<< filter->name() << std::endl;
					if ( ! filter->apply(*seed_pose) ) {
						TR << "Segment did not pass filter: " <<filter->name() << " continueing..." << std::endl;
						passed_filters = false;
						break;
					}
				}
				if ( passed_filters && seg_filters_.count(selected_pose->get_segfile_path()) ) {
					for ( FilterCOP filter : seg_filters_.at(selected_pose->get_segfile_path()) ) {
						TR << "Running Segfile-specific filter "<< filter->name() << std::endl;
						if ( ! filter->apply(*seed_pose) ) {
							TR << "Segment did not pass filter: " <<filter->name() << " continueing..." << std::endl;
							passed_filters = false;
							break;
						}
					}
				}
				if ( ! passed_filters ) {
					//TR << "Incrementing total attempts" << std::endl;
					continue; //Go to next segment!
				}
				if ( selected_pose == nullptr ) {
					continue; //Indicates EE, w/o terminal matches.
				}
			}

			//Without filtering, this is set to true.
			if ( passed_filters ) {
				break;
			};
		}

	}

	//Make sure we fail if for some reason nothing is found and our settings are too stringent.
	if ( selected_pose == nullptr ) {
		utility_exit_with_message("PickRandomSegmentMover: Could not find any seeds that match!");
	} else if ( seed_pose == nullptr ) {
		seed_pose = selected_pose->get_source_pose_op(false /*copy pose*/);
	}
	pose = *seed_pose;
	//Attach MT
	assert(! all_L_dssp(pose));


	TR << "pose secstruct: " << pose.secstruct() << std::endl;
	TR << "pose chain numbers: " << pose.chain(1) << " " << pose.chain(pose.size()) << std::endl;
	if ( pose.pdb_info()==nullptr ) {
		core::pose::PDBInfoOP new_info = core::pose::PDBInfoOP(new core::pose::PDBInfo(pose,true));
		pose.pdb_info(new_info);
	}
	for ( core::Size current_residue = 1; current_residue <= pose.size(); current_residue++ ) {
		pose.pdb_info()->add_reslabel(current_residue,selected_pose->get_filename());
	}

}

utility::vector1< std::string >
PickRandomSegmentMover::get_segment_file_paths() const{
	return segment_file_paths_;
}
data_storage::TerminalDSSPSortedPoseVectorOP
PickRandomSegmentMover::get_pose_vector() const{
	return pose_vector_;
}

void
PickRandomSegmentMover::set_segment_file_paths(utility::vector1<std::string> segment_file_paths){

	segment_file_paths_ = segment_file_paths;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
PickRandomSegmentMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
PickRandomSegmentMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data_map
)
{


	set_ss_defaults();
	filters_.clear();
	pose_vector_->clear_all_vectors();

	if ( tag->hasOption("segment_file_path") ) {
		std::string segment_file_path = tag->getOption<std::string>("segment_file_path");
		segment_file_paths_ = utility::string_split(segment_file_path, ',');
	}

	match_terminal_sheet_lengths_ = tag->getOption< bool >("match_terminal_sheet_lengths",
		match_terminal_sheet_lengths_);

	terminal_sheet_length_max_delta_ = tag->getOption< core::Size >("terminal_sheet_length_max_delta",
		terminal_sheet_length_max_delta_);


	//Add filters pre-defined or given as subelement
	if ( tag->hasOption("filters") ) {
		utility::vector1< std::string > filter_names = utility::string_split(tag->getOption< std::string >("filters"), ',');
		for ( std::string const & fname : filter_names ) {
			filters_.push_back(protocols::rosetta_scripts::parse_filter( fname, data_map ) );
			TR<<"Added filter "<<fname << std::endl;
		}
	}

	//Subelement filters
	utility::vector1< std::string > ss_elements;
	ss_elements.push_back("E");
	ss_elements.push_back("H");
	ss_elements.push_back("L");

	even_sampling_ = tag->getOption< bool >("even_sampling", even_sampling_);

	//Now we set per-SS and seg settings
	utility::vector1< std::string > ss_types;
	ss_types.push_back("E");
	ss_types.push_back("L");
	ss_types.push_back("H");

	for ( TagCOP inner_tag : tag->getTags() ) {

		std::string inner = inner_tag->getName();
		if ( inner == "Segment" ) {
			std::string seg_path = inner_tag->getOption< std::string >("segment_file_path");
			if ( ! segment_file_paths_.has_value(seg_path) ) {
				segment_file_paths_.push_back(seg_path);
			}

			if ( inner_tag->hasOption("filters") ) {
				utility::vector1< std::string > filter_names = utility::string_split(inner_tag->getOption< std::string >("filters"), ',');
				for ( std::string const & fname : filter_names ) {
					seg_filters_[seg_path].push_back(protocols::rosetta_scripts::parse_filter( fname, data_map ) );
					TR<<"Added filter "<< fname << " for seg " << seg_path << std::endl;
				}
			}
		} else if ( ss_types.has_value(inner) ) {
			char ss_char = inner[0];

			//Grab options
			core::Size max_term_size = inner_tag->getOption< core::Size >("max_terminal_length");
			core::Size min_term_size = inner_tag->getOption< core::Size >("min_terminal_length");


			//Setup the Settings struct.
			SeedSettings ss_setting;
			ss_setting.min_terminal_length_ = min_term_size;
			ss_setting.max_terminal_length_ = max_term_size;


			ss_settings_[ss_char] = ss_setting;
		} else {
			continue;
		}

	}

}
void PickRandomSegmentMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	using namespace protocols::moves;

	XMLSchemaAttribute min_term_size =  XMLSchemaAttribute::attribute_w_default( "min_terminal_length", xsct_positive_integer, "Minimum number of residues in each terminal SS block","1" );

	XMLSchemaAttribute max_term_size = XMLSchemaAttribute::attribute_w_default( "max_terminal_length", xsct_positive_integer, "Maximum number of residues in each terminal SS block","100" );


	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "segment_file_path", xs_string, "Path to directory containing segment file. Can be multiple paths, comma-delimited if you are feeling daring today. Can also be given in Segment subblock. ","NO_SEGMENT_FILE" )
		+ XMLSchemaAttribute::attribute_w_default("even_sampling", xsct_rosetta_bool, "Should we perform even sampling for multiple segment files? If false, those with higher numbers will dominate", "true")
		+ min_term_size
		+ max_term_size
		+ XMLSchemaAttribute::attribute_w_default( "match_terminal_sheet_lengths", xsct_rosetta_bool, "Match length of terminal SS blocks?  Uses terminal_length_delta_max to set allowable length mismatch.","false" )
		+ XMLSchemaAttribute( "filters", xs_string, "Comma-delimited list of previously defined filter(s) that will be run on each chosen segment.")
		+ XMLSchemaAttribute::attribute_w_default( "terminal_sheet_length_delta_max", xsct_positive_integer, "Allowable length mismatch if match_terminal_lengths is true","1" );

	//FilterFactory::get_instance()->define_filter_xml_schema( xsd );
	//subelements.add_group_subelement( & FilterFactory::get_instance()->filter_xml_schema_group_name );

	AttributeList ss_subelement_attlist;

	//Copy attributes and make them required.
	//Copy may not be nessessary, but this is to make sure both instances are different.
	XMLSchemaAttribute min_term_size_req = min_term_size; min_term_size_req.is_required(true);
	XMLSchemaAttribute max_term_size_req = max_term_size; max_term_size_req.is_required(true);

	//Blow away the default value or the schema check will complain.
	min_term_size_req.default_value("");
	max_term_size_req.default_value("");

	ss_subelement_attlist
		+ min_term_size_req
		+ max_term_size_req;

	AttributeList seg_subelement_attlist;
	seg_subelement_attlist
		+ XMLSchemaAttribute::required_attribute( "segment_file_path", xs_string, "Path to directory containing segment file. Can be multiple paths, comma-delimited.")
		+ XMLSchemaAttribute( "filters", xs_string, "Comma-delimited list of previously defined filter(s) that will be run on each chosen segment");

	XMLSchemaSimpleSubelementList subelements;
	subelements
		.add_simple_subelement("Segment", seg_subelement_attlist, "Add a segment file here, optionally control specific filters applied to this segment (in addition to globally set filters");

	subelements
		.add_simple_subelement( "E", ss_subelement_attlist, "Specify a few settings for specific secondary structure.  If given here, will override the default, or any settings given in the primary PickRandomSegmentMover tag" );

	subelements
		.add_simple_subelement( "L", ss_subelement_attlist, "Specify a few settings for specific secondary structure.  If given here, will override the default, or any settings given in the primary PickRandomSegmentMover tag" );

	subelements
		.add_simple_subelement( "H", ss_subelement_attlist, "Specify a few settings for specific secondary structure.  If given here, will override the default, or any settings given in the primary PickRandomSegmentMover tag" );


	//// Create the XSD ////

	//Since we have two types of subelements, we need to do this manually.
	std::string description = "Select a random segment from a segment file or set of segment files.  Replace the pose with the segment.  Used to generate seeds for SewAnything";


	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), description, attlist, subelements );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PickRandomSegmentMover::fresh_instance() const
{
	return utility::pointer::make_shared< PickRandomSegmentMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PickRandomSegmentMover::clone() const
{
	return utility::pointer::make_shared< PickRandomSegmentMover >( *this );
}

std::string PickRandomSegmentMover::get_name() const {
	return mover_name();
}

std::string PickRandomSegmentMover::mover_name() {
	return "PickRandomSegmentMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
PickRandomSegmentMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< PickRandomSegmentMover >();
}

std::string
PickRandomSegmentMoverCreator::keyname() const
{
	return PickRandomSegmentMover::mover_name();
}

void PickRandomSegmentMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PickRandomSegmentMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, PickRandomSegmentMover const & mover )
{
	mover.show(os);
	return os;
}

} //movers
} //pose_sewing
} //protocols
