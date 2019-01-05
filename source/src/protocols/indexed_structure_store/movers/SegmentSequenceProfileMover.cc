// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/indexed_structure_store/movers/SegmentSequenceProfileMover.cc

#include <algorithm>

#include "boost/format.hpp"
#include "boost/range/algorithm/copy.hpp"

#include <utility/io/ozstream.hh>

#include <core/sequence/SequenceProfile.hh>

#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

#include <protocols/indexed_structure_store/StructureStore.hh>
#include <protocols/indexed_structure_store/StructureStoreManager.hh>
#include <protocols/indexed_structure_store/search/QueryDatabase.hh>
#include <protocols/indexed_structure_store/search/QueryDatabase.json.hh>
#include <protocols/indexed_structure_store/vector_tools.hh>
#include <protocols/indexed_structure_store/utility.hh>
#include <protocols/indexed_structure_store/pose_utility.hh>

#include <protocols/indexed_structure_store/SegmentSequenceProfile.hh>
#include <protocols/indexed_structure_store/movers/SegmentSequenceProfileMover.hh>
#include <protocols/indexed_structure_store/movers/SegmentSequenceProfileMoverCreator.hh>

#include <protocols/simple_moves/FavorSequenceProfile.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.indexed_structure_store.movers.SegmentSequenceProfileMover" );

namespace protocols { namespace indexed_structure_store { namespace  movers {

std::string
SegmentSequenceProfileMover::pssm_output_order = {
'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};


moves::MoverOP
SegmentSequenceProfileMover::clone() const {
	return utility::pointer::make_shared< SegmentSequenceProfileMover >( *this );
}

void
SegmentSequenceProfileMover::apply(core::pose::Pose & pose) {

	using ArrayXaa = SegmentSequenceProfileResult::ArrayXaa;

	// Resolve residue selector into a contiguous residue span
	if ( !selector_ ) {
		utility_exit_with_message("SegmentSequenceProfileMover requires target residue selector.");
	}
	core::select::residue_selector::ResidueSubset subset = selector_->apply( pose );

	if ( std::none_of(subset.begin(), subset.end(),[](bool s){return s;}) ) {
		utility_exit_with_message("SegmentSequenceProfileMover residue selector did not result in selection.");
	}

	auto selection_begin = std::find(subset.begin(), subset.end(), true);
	auto selection_end = std::find(selection_begin, subset.end(), false);

	if ( std::any_of(selection_end, subset.end(), [](bool s){return s;}) ) {
		utility_exit_with_message("SegmentSequenceProfileMover residue selector did not result in contiguous selection.");
	}

	core::Size span_start = selection_begin - subset.begin();
	core::Size span_end = selection_end - subset.begin();

	TR.Info << "Resolved segment range: [" << span_start << ", " << span_end << "]" << std::endl;

	SegmentSequenceProfile profiler(profile_config_);
	SegmentSequenceProfileResult profile = profiler.segment_profile(
		*structure_store(), *structure_database(),
		pose, span_start + 1, span_end + 1
	);

	ArrayXaa full_profile = ArrayXaa::Zero(pose.total_residue(), core::chemical::num_canonical_aas);
	full_profile.block(
		span_start, 0, span_end - span_start, core::chemical::num_canonical_aas) = profile.log_odds;

	core::sequence::SequenceProfile seq_profile = to_sequence_profile(full_profile);

	if ( apply_profile_ ) {
		TR << "Applying profile via FavorSequenceProfile." << std::endl;

		protocols::simple_moves::FavorSequenceProfile fsp_mover;
		fsp_mover.set_scaling(scaling_);
		fsp_mover.set_weight(weight_);
		fsp_mover.set_profile(seq_profile);

		fsp_mover.apply(pose);
	}

	if ( !output_lookup_count_.empty() ) {
		TR.Info << boost::format("Storing lookup hit count %s : %i") % output_lookup_count_ % profile.query_results.size() << std::endl;
		core::pose::setPoseExtraScore(pose, output_lookup_count_, profile.query_results.size());
	}

	if ( !output_pssm_filename_.empty() ) {
		TR.Info << "Writing output pssm:"  << output_pssm_filename_ << std::endl;
		utility::io::ozstream outfile(output_pssm_filename_.c_str(), std::ios::out | std::ios::binary);
		if ( !outfile ) {
			utility_exit_with_message("Unable to open file for write: " + output_pssm_filename_);
		}

		outfile << protocols::jd2::current_output_name() << std::endl;
		outfile << "SegementSequenceProfile results"
			<< " hits: " << profile.query_results.size()
			<< " rmsd_tolerance: " << profiler.config.rmsd_tolerance
			<< " pseudocount: " << profiler.config.pseudocount << std::endl;

		write_profile(seq_profile, outfile);

		outfile.close();
	}

	if ( !output_pssm_inline_.empty() ) {
		if ( !protocols::jd2::jd2_used() ) {
			utility_exit_with_message("SegmentSequenceProfileMover inline pssm output not supported outside of jd2.");
		}

		std::ostringstream oss;
		std::string header_line = "#BEGIN_SEGMENT_SEQUENCE_PROFILE_PSSM " + output_pssm_inline_;
		std::string footer_line = "#END_SEGMENT_SEQUENCE_PROFILE_PSSM " + output_pssm_inline_;

		auto & output_strings = protocols::jd2::get_current_job()->get_strings();
		output_strings.remove_if(
			[header_line](std::string s){
				return s.substr(0, header_line.size()) == header_line;}
		);

		oss << header_line << std::endl;
		oss << "SegementSequenceProfile results"
			<< " hits: " << profile.query_results.size()
			<< " rmsd_tolerance: " << profiler.config.rmsd_tolerance
			<< " pseudocount: " << profiler.config.pseudocount
			<< std::endl;

		write_profile(seq_profile, oss);
		oss << footer_line << std::endl;

		output_strings.push_back(oss.str());
	}
}

core::sequence::SequenceProfile
SegmentSequenceProfileMover::to_sequence_profile(SegmentSequenceProfileResult::ArrayXaa profile) {

	utility::vector1<std::string> result_alphabet;

	for ( char aa : pssm_output_order ) {
		result_alphabet.push_back(std::string(1, aa));
	}

	utility::vector1<utility::vector1<core::Real>> result_profile;
	for ( core::SSize i = 0; i < profile.rows(); ++i ) {
		utility::vector1<core::Real> prow;

		for ( char aa : pssm_output_order ) {
			prow.push_back(
				profile(i, core::chemical::aa_from_oneletter_code(aa) - core::chemical::first_l_aa));
		}

		result_profile.push_back(prow);
	}

	core::sequence::SequenceProfile result;
	result.alphabet(result_alphabet);
	result.profile(result_profile);
	result.sequence(std::string(profile.rows(), 'X'));
	return result;
}

void
SegmentSequenceProfileMover::write_profile(core::sequence::SequenceProfile profile, std::ostream & out) {
	out << "          ";
	boost::copy(profile.alphabet(), std::ostream_iterator<std::string>(out, "  "));
	out << std::endl;

	boost::format line_format(
		"%5i %1s  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f");

	for ( core::Size i = 0; i < profile.length(); ++i ) {
		line_format.clear_binds();
		line_format.bind_arg(1, i+1);
		line_format.bind_arg(2, profile.sequence()[i]);
		for ( core::Size j = 0; j < profile.width(); ++j ) {
			line_format.bind_arg(
				line_format.cur_arg(),
				profile.prof_row(i+1)[j+1]
			);
		}
		out << line_format << std::endl;
	}
}

void SegmentSequenceProfileMover::structure_store_path(std::string path) {
	if ( path != structure_store_path_ ) {
		structure_store_.reset();
		structure_database_.reset();
		structure_store_path_ = path;
	}
}

std::string SegmentSequenceProfileMover::structure_store_path() {
	return structure_store_path_;
}

StructureStoreOP SegmentSequenceProfileMover::structure_store() {
	init_structure_store();
	return structure_store_;
}

search::StructureDatabaseOP SegmentSequenceProfileMover::structure_database() {
	init_structure_store();
	return structure_database_;
}


void
SegmentSequenceProfileMover::init_structure_store() {
	if ( !structure_store_ ) {
		runtime_assert(!structure_database_);
		structure_store_ = StructureStoreManager::get_instance()->load_structure_store(structure_store_path());
		structure_database_ = utility::pointer::make_shared< search::StructureDatabase >();
		structure_database_->initialize(structure_store_->residue_entries);
	}
}

void
SegmentSequenceProfileMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {
	structure_store_path( tag->getOption< std::string >( "structure_store" ) );
	if ( tag->hasOption("residue_selector") ) {
		selector_ = protocols::rosetta_scripts::parse_residue_selector(tag, data);
	}
	profile_config_.rmsd_tolerance = tag->getOption< core::Real >( "rmsd_tolerance", .5);
	profile_config_.pseudocount = tag->getOption< core::Real >( "pseudocount", 1);

	output_lookup_count_ = tag->getOption< std::string >( "output_lookup_count", "");
	output_pssm_filename_ = tag->getOption< std::string >( "output_pssm_filename", "");
	output_pssm_inline_ = tag->getOption< std::string >( "output_pssm_inline", "");

	apply_profile_ = tag->getOption<bool>("apply_profile", true);
	scaling_ = tag->getOption<std::string>("scaling", "prob");
	weight_ = tag->getOption<core::Real>("weight", 1.0);
}

void
SegmentSequenceProfileMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;

	XMLSchemaRestriction scaling_type;
	scaling_type.name( "favor_seqprof_scaling_type" );
	scaling_type.add_restriction( xsr_enumeration, "prob" );
	scaling_type.add_restriction( xsr_enumeration, "none" );
	scaling_type.add_restriction( xsr_enumeration, "global" );
	xsd.add_top_level_element( scaling_type );

	attlist + XMLSchemaAttribute::required_attribute("structure_store" , xs_string , "Path of target struture store.");
	core::select::residue_selector::attributes_for_parse_residue_selector_when_required(
		attlist, "residue_selector",
		"Residue selector spanning a contiguous segment region to profile." );
	attlist + XMLSchemaAttribute::attribute_w_default("rmsd_tolerance" , xsct_real, "Lookup endpoint rmsd tolerance.", ".5");
	attlist + XMLSchemaAttribute::attribute_w_default("pseudocount" , xsct_real, "Sequence profile pseudocount.", "1");
	attlist + XMLSchemaAttribute::attribute_w_default("output_lookup_count" , xs_string, "Output count of lookup hits as score entry.", "");
	attlist + XMLSchemaAttribute::attribute_w_default("output_pssm_filename" , xs_string, "Output generated log-odds pssm with given output filename.", "");
	attlist + XMLSchemaAttribute::attribute_w_default("output_pssm_inline" , xs_string, "Output generated log-odds pssm inline with pdb under given header name.", "");

	attlist + XMLSchemaAttribute::attribute_w_default("apply_profile", xsct_rosetta_bool, "Apply profile via FavorSequenceProfile.", "true");

	attlist + XMLSchemaAttribute::attribute_w_default( "scaling", "favor_seqprof_scaling_type",
		"Set how to scale the given values, see FavorSequenceProfile."
		"\"prob\"=Boltzmann-weighted probability based on the profile score"
		"\"global\"= global linear fixed-zero rescaling"
		"such that all (pre-weighted) values fall in the range of -1.0 to 1.0"
		"\"none\" does no adjustment of values.", "prob" );
	attlist + XMLSchemaAttribute::attribute_w_default( "weight", xsct_real, "Adjust the post-scaling strength of the constraints, see FavorSequenceProfile.", "1" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Perform alignment-based sequence profile generation. TODO fordas", attlist );
}

protocols::moves::MoverOP
SegmentSequenceProfileMoverCreator::create_mover() const {
	return utility::pointer::make_shared< SegmentSequenceProfileMover >();
}

std::string
SegmentSequenceProfileMover::get_name() const {
	return mover_name();
}

std::string
SegmentSequenceProfileMover::mover_name() {
	return "SegmentSequenceProfileMover";
}

std::string
SegmentSequenceProfileMoverCreator::keyname() const {
	return SegmentSequenceProfileMover::mover_name();
}

void
SegmentSequenceProfileMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SegmentSequenceProfileMover::provide_xml_schema(xsd);
}

} } }
