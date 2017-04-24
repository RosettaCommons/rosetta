// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/LongestContinuousPolarSegmentFilter.cc
/// @brief This filter computes the longest continuous stretch of polar residues within a pose or selection.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#include <protocols/simple_filters/LongestContinuousPolarSegmentFilter.hh>
#include <protocols/simple_filters/LongestContinuousPolarSegmentFilterCreator.hh>

//Core includes
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AA.hh>

//Protocols includes
#include <protocols/rosetta_scripts/util.hh>

//XSD includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.LongestContinuousPolarSegmentFilter" );

namespace protocols {
namespace simple_filters {

/// @brief Constructor.
LongestContinuousPolarSegmentFilter::LongestContinuousPolarSegmentFilter():
	protocols::filters::Filter( "LongestContinuousPolarSegmentFilter" ),
	exclude_chain_termini_(true),
	count_gly_as_polar_(true),
	filter_out_high_(true),
	cutoff_(5),
	residue_selector_(nullptr)
{}

/// @brief Destructor.
LongestContinuousPolarSegmentFilter::~LongestContinuousPolarSegmentFilter()
{}

/// @brief Parse tag to allow RosettaScripts XML to call this mover.
void
LongestContinuousPolarSegmentFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	set_exclude_chain_termini( tag->getOption<bool>( "exclude_chain_termini", exclude_chain_termini() ) );
	set_count_gly_as_polar( tag->getOption<bool>("count_gly_as_polar", count_gly_as_polar()) );
	set_filter_out_high( tag->getOption<bool>( "filter_out_high", filter_out_high() ) );
	set_cutoff( tag->getOption<core::Size>( "cutoff", cutoff() ) );

	core::select::residue_selector::ResidueSelectorCOP selector( protocols::rosetta_scripts::parse_residue_selector( tag, data ) );
	if ( selector != nullptr ) { set_residue_selector( selector ); }
}

protocols::filters::FilterOP
LongestContinuousPolarSegmentFilter::clone() const
{
	return protocols::filters::FilterOP( new LongestContinuousPolarSegmentFilter( *this ) );
}


protocols::filters::FilterOP
LongestContinuousPolarSegmentFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new LongestContinuousPolarSegmentFilter );
}

/// @brief Returns true if the structure passes the filter, false otherwise.
bool
LongestContinuousPolarSegmentFilter::apply( core::pose::Pose const &pose ) const
{
	core::Size longest_count, longest_start, longest_end;
	compute( pose, residue_selector_, exclude_chain_termini_, longest_count, longest_start, longest_end );

	if ( TR.visible() && longest_count != 0 ) {
		TR << "In the current pose, the longest stretch of polar residues is " << longest_count << " residues.  (Residue " << pose.residue(longest_start).name3() << longest_start << " through " << pose.residue(longest_end).name3() << longest_end << ".)" << std::endl;
	}

	if ( filter_out_high_ ) {
		if ( longest_count > cutoff_ ) {
			TR << "The length of the longest stretch of polar residues (" << longest_count << ") is greater than the maximum cutoff threshold (" << cutoff_ << ").  Filter fails." << std::endl;
			return false;
		} else {
			TR << "The length of the longest stretch of polar residues (" << longest_count << ") is less than or equal to than the maximum cutoff threshold (" << cutoff_ << ").  Filter passes." << std::endl;
			return true;
		}
	}

	if ( longest_count < cutoff_ ) {
		TR << "The length of the longest stretch of polar residues (" << longest_count << ") is less than the minimum cutoff threshold (" << cutoff_ << ").  Filter fails." << std::endl;
		return false;
	} else {
		TR << "The length of the longest stretch of polar residues (" << longest_count << ") is greater than or equal to than the minimum cutoff threshold (" << cutoff_ << ").  Filter passes." << std::endl;
		return true;
	}

	return true; //To make compiler happy
}

/// @brief Required for reporting score values.
core::Real
LongestContinuousPolarSegmentFilter::report_sm( core::pose::Pose const &pose ) const
{
	core::Size longest_count, longest_start, longest_end;
	compute( pose, residue_selector_, exclude_chain_termini_, longest_count, longest_start, longest_end );

	return static_cast< core::Real >( longest_count );
}

/// @brief Allows printing data to a stream.
void
LongestContinuousPolarSegmentFilter::report( std::ostream &os, core::pose::Pose const &pose ) const
{
	core::Size longest_count, longest_start, longest_end;
	compute( pose, residue_selector_, exclude_chain_termini_, longest_count, longest_start, longest_end );

	os << "The LongestContinuousPolarSegmentFilter reports that the longest stretch of polar residues in the pose, " << ( exclude_chain_termini_ ? "excluding chain terminal segments," : "including chain terminal segments," ) << " is " << longest_count << " residues long.";
	if ( longest_count > 0 ) {
		os << "  It runs from " << pose.residue(longest_start).name3() << longest_start << " through " << pose.residue(longest_end).name3() << longest_end << ".";
	}
	os << std::endl;
}

/// @brief Set the residue selector.
void
LongestContinuousPolarSegmentFilter::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	runtime_assert_string_msg( selector_in != nullptr, "Error in protocols::simple_moves::LongestContinuousPolarSegmentFilter::set_residue_selector(): The pointer passed to this function is null." );
	residue_selector_ = selector_in;
}


std::string LongestContinuousPolarSegmentFilter::name() const {
	return class_name();
}

std::string LongestContinuousPolarSegmentFilter::class_name() {
	return "LongestContinuousPolarSegment";
}

void LongestContinuousPolarSegmentFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default("exclude_chain_termini", xsct_rosetta_bool, "If true, polar stretches at the ends of chains are not counted.  If false, they are.  True by default.  Note that if this is set to true, an entirely polar chain will not be counted.", "true")
		+ XMLSchemaAttribute::attribute_w_default("count_gly_as_polar", xsct_rosetta_bool, "If true, glycine is considered \"polar\" for purposes of this filter.  True by default.", "true")
		+ XMLSchemaAttribute::attribute_w_default("filter_out_high", xsct_rosetta_bool, "If true, poses with more than the cutoff number of residues in the longest polar stretch will be rejected.  If false, poses with fewer than the cutoff number of residues in the longest polar stretch will be rejected.  True by default.", "true")
		+ XMLSchemaAttribute::attribute_w_default( "cutoff", xsct_non_negative_integer, "The maximum (or minimum, if \"filter_out_high\" is set to \"false\") number of residues in the longest polar stretch that will still allow the pose to pass this filter.  Default 5.", "5" )
		;

	protocols::rosetta_scripts::attributes_for_parse_residue_selector( attlist, "An optional, previously-defined residue selector.  If provided, the filter will only consider stretches of polar residues that have at least one residue in the selection.  Not used if not specified." );

	protocols::filters::xsd_type_definition_w_attributes(
		xsd,
		class_name(),
		"The LongestContinuousPolarSegment filter counts the number of amino acid residues in the longest continuous stretch of polar amino acids in any chain of your pose (or, optionally, in a residue selection).  Optionally, it can ignore the polar stretches at the N- and C-termini of chains.",
		attlist );
}


/// @brief Given a pose, actually compute the number of residues in the longest stretch, and return this value,
/// along with the start and end indices (in Rosetta numbering).
/// @details Returns 0, 0, 0 if no polar stretch could be found.
/// @param[in] pose The pose to analyse.
/// @param[in] selector An optional const-owning pointer to a ResidueSelector. If provided, only those stretches that have at least one residue selected by the ResidueSelector
/// will be counted.
/// @param[in] ignore_termini If true, stretches at the N- and C-termini of chains are ignored.  If false, they're counted.
/// @param[out] longest_stretch_res_count The number of polar residues in the longest polar stretch.  Set to 0 if no polar stretch is found.
/// @param[out] longest_stretch_start The index, in Rosetta numbering, of the first residue of the longest polar stretch.  Set to 0 if no polar stretch is found.
/// @param[out] longest_stretch_end The index, in Rosetta numbering, of the last residue of the longest polar stretch.  Set to 0 if no polar stretch is found.
void
LongestContinuousPolarSegmentFilter::compute(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSelectorCOP selector,
	bool const ignore_termini,
	core::Size &longest_stretch_res_count,
	core::Size &longest_stretch_start,
	core::Size &longest_stretch_end
) const {

	longest_stretch_res_count = 0;
	longest_stretch_start = 0;
	longest_stretch_end = 0;

	//Set up the residue selection:
	core::select::residue_selector::ResidueSubset selected_res( pose.total_residue(), true );
	if ( selector != nullptr ) {
		selected_res = selector->apply( pose );
	}

	//Counts, start, and end for a candidate stretch.
	core::Size count(0), start(0);
	bool is_Nterminal_stretch(false), is_polar_stretch(false);

	for ( core::Size ir(1), irmax(pose.total_residue()); ir<=irmax; ++ir ) { //Loop through all residues in the pose
		if ( !pose.residue_type(ir).is_protein() ) continue; //Only count protein residues.

		if ( !is_polar_stretch ) { //If we're not in a polar stretch.
			if ( is_counted( pose.residue_type(ir) ) ) { //If this is the first residue of a polar stretch
				is_Nterminal_stretch = (pose.residue(ir).connected_residue_at_lower() == 0); // Check whether we're at an N-terminus, and set is_Nterminal_stretch appropriately.
				is_polar_stretch = true; //We're now in a polar stretch.
				count = 1; //Start counting.
				start=ir;
			}
		} else { //If we ARE already in a polar stretch.
			if ( is_counted( pose.residue_type(ir) ) ) {
				++count; //Increment the count if the next residue is polar.
			}
			if ( !is_counted( pose.residue_type(ir) ) || pose.residue(ir).connected_residue_at_upper() == 0 ) { //If we have reached the end of a polar stretch, either by hitting a hydrophobic or by reaching the chain end.
				core::Size const end( is_counted( pose.residue_type(ir) ) ? ir : ir-1 ); //The last residue is either the previous one (if we're ending because this is apolar) or the current one (if we're ending because this is a C-teminus and polar).
				if ( //If we're either counting terminal stretches, or (this is not an N-terminal streatch and this is not a C-terminal stretch):
						!ignore_termini ||
						( !is_Nterminal_stretch && !(pose.residue(ir).connected_residue_at_upper() == 0) )
						) {
					if ( count > longest_stretch_res_count ) {
						//Check whether the current stretch has at least one selected residue:
						bool stretch_is_selected( selector == nullptr );
						if ( selector!=nullptr ) {
							for ( core::Size ir2(start); ir2<end; ++ir2 ) {
								if ( selected_res[ir2] == true ) {
									stretch_is_selected=true;
									break;
								}
							}
						}

						if ( stretch_is_selected ) {
							longest_stretch_res_count = count;
							longest_stretch_start = start;
							longest_stretch_end = end;
						}
					}
				}
				is_polar_stretch=false; //No longer in a polar stretch.
				count=0; //Resetting for safe measure (though not really necessary).
				start=0; //Resetting for safe measure (though not really necessary).
			}
		}
	}

	if ( longest_stretch_res_count == 0 ) {
		TR.Warning << "No polar stretch was found in this pose!" << std::endl;
	}
}

/// @brief Given a residue type, determine whether it's one of the types that this filter
/// should count.
/// @details Based on whether the residue type has the POLAR property.  Special-case exception
/// is made for glycine, depending on whether count_gly_as_polar_ is true.
bool
LongestContinuousPolarSegmentFilter::is_counted(
	core::chemical::ResidueType const &restype
) const {
	if ( restype.is_polar() ) return true;
	if ( count_gly_as_polar() && restype.aa() == core::chemical::aa_gly ) return true;
	return false;
}


/////////////// Creator ///////////////

protocols::filters::FilterOP
LongestContinuousPolarSegmentFilterCreator::create_filter() const
{
	return protocols::filters::FilterOP( new LongestContinuousPolarSegmentFilter );
}

std::string
LongestContinuousPolarSegmentFilterCreator::keyname() const
{
	return LongestContinuousPolarSegmentFilter::class_name();
}

void LongestContinuousPolarSegmentFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LongestContinuousPolarSegmentFilter::provide_xml_schema( xsd );
}

} //protocols
} //simple_filters
