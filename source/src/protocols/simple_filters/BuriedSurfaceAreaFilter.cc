// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/BuriedSurfaceAreaFilter.cc
/// @brief Calculates buried surface area (exposed surface area minus total surface area, on a per-residue basis).  Accepts
/// a residue selector to allow buried subsets to be considered.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#include <protocols/simple_filters/BuriedSurfaceAreaFilter.hh>
#include <protocols/simple_filters/BuriedSurfaceAreaFilterCreator.hh>

//Core includes
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/scoring/sasa.hh>
#include <core/id/AtomID_Map.hh>
#include <core/select/residue_selector/util.hh>

//XSD includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.BuriedSurfaceAreaFilter" );

namespace protocols {
namespace simple_filters {

BuriedSurfaceAreaFilter::BuriedSurfaceAreaFilter():
	protocols::filters::Filter( "BuriedSurfaceAreaFilter" ),
	residue_selector_(nullptr),
	select_only_FAMILYVW_(false),
	cutoff_area_(500),
	filter_out_low_(true),
	probe_radius_(2.0)
{

}

BuriedSurfaceAreaFilter::~BuriedSurfaceAreaFilter()
{}

void
BuriedSurfaceAreaFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector( protocols::rosetta_scripts::parse_residue_selector( tag, data ) );
	}
	set_select_only_FAMILYVW( tag->getOption<bool>( "select_only_FAMILYVW", select_only_FAMILYVW() ) );
	set_cutoff_buried_surface_area( tag->getOption<core::Real>( "cutoff_buried_surface_area", cutoff_buried_surface_area() ) );
	set_filter_out_low( tag->getOption<bool>( "filter_out_low", filter_out_low() ) );
	set_probe_radius( tag->getOption<core::Real>( "probe_radius", probe_radius() ) );
}

protocols::filters::FilterOP
BuriedSurfaceAreaFilter::clone() const
{
	return protocols::filters::FilterOP( new BuriedSurfaceAreaFilter( *this ) );
}


protocols::filters::FilterOP
BuriedSurfaceAreaFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new BuriedSurfaceAreaFilter );
}

bool
BuriedSurfaceAreaFilter::apply( core::pose::Pose const &pose ) const
{
	core::Real returnval;
	return compute( TR, pose, returnval );
}

core::Real
BuriedSurfaceAreaFilter::report_sm( core::pose::Pose const &pose ) const
{
	core::Real returnval;
	compute( TR, pose, returnval );
	return returnval;
}

void
BuriedSurfaceAreaFilter::report( std::ostream &os, core::pose::Pose const &pose ) const
{
	core::Real returnval;
	compute(os, pose, returnval);
}

std::string BuriedSurfaceAreaFilter::name() const {
	return class_name();
}

std::string BuriedSurfaceAreaFilter::class_name() {
	return "BuriedSurfaceArea";
}

void BuriedSurfaceAreaFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "select_only_FAMILYVW", xsct_rosetta_bool, "If true, then only hydrophobic residues and alanine (phe, ala, met, ile, leu, tyr, val, or trp) are counted.  If false (the default), then all residues are counted.  Note that this is combined with AND logic to create the final residue selection if a residue selector is also provided (that is, residues that are selected by the selector AND which are in the set FAMILYVW are used to compute the buried surface area).", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "filter_out_low", xsct_rosetta_bool, "If true (the default), then poses with buried surface area below the cutoff are rejected.  If false, then poses with buried surface area above the cutoff are rejected.", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "cutoff_buried_surface_area", xsct_real, "The buried surface area below which (or above which, if \"filter_out_low\" is false) a pose is rejected.  Defaults to 500 square Angstroms, an arbitrarily-chosen value.", "500" )
		+ XMLSchemaAttribute::attribute_w_default( "probe_radius", xsct_real, "The radius for the probe used in the rolling-ball algorithm for determining solvent-accessible surface area.  The buried surface area is the total minus the solvent accessible.  Defaults to 2.0 Angstroms.", "2.0" );
	core::select::residue_selector::attributes_for_parse_residue_selector_default_option_name( attlist, "An optional, previously-defined residue selector.  If provided, then only the selected residues are used in computing buried surface area.  If not provided, then all residues are used." );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "This filter computes the total buried surface area for a pose or a subset of a pose selected by a residue selector.  It discards poses with buried surface area below (or above) a user-defined threshold value.  Note that this filter operates ONLY on canonical L-amino acids, their D-amino acid counterparts, and glycine; it will compute a buried surface area of zero for all other types.", attlist );
}

/////////////// Setters ///////////////

/// @brief Sets the residue selector to use to select a subset of residues for which to calculate
/// buried surface area.
/// @details Copies the input owning pointer; does not clone.  This means that residue selectors could be
/// shared with other Rosetta modules.
void
BuriedSurfaceAreaFilter::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	runtime_assert_string_msg( selector_in != nullptr, "Error in protocols::simple_filters::BuriedSurfaceAreaFilter::set_residue_selector(): The pointer passed to this function was null!" );
	residue_selector_ = selector_in;
}

/// @brief Set whether only hydrophobic residues and alanine (FAMILYVW) are considered.  False by default.
/// @details The selection FAMILYVW is combined with the residue selector, if specified, using AND logic.
void
BuriedSurfaceAreaFilter::set_select_only_FAMILYVW(
	bool const setting
) {
	select_only_FAMILYVW_ = setting;
}

/// @brief Set the cutoff buried surface area below which (or above which, if filter_out_low_ is false) structures
/// are discarded.
void
BuriedSurfaceAreaFilter::set_cutoff_buried_surface_area(
	core::Real const &setting
) {
	runtime_assert_string_msg( setting >= 0.0, "Error in protocols::simple_filters::BuriedSurfaceAreaFilter::set_cutoff_buried_surface_area(): The cutoff area must be greater than or equal to zero." );
	cutoff_area_ = setting;
}

/// @brief Set whether structures with less than the cutoff buried area or greater than the cutoff
/// buried area are discarded.
/// @details If true, structures with less than the cutoff buried area are discrded.
void
BuriedSurfaceAreaFilter::set_filter_out_low(
	bool const setting
) {
	filter_out_low_ = setting;
}

/// @brief Set the probe radius for calculating solvent-accessible surface area.
///
void
BuriedSurfaceAreaFilter::set_probe_radius(
	core::Real const &setting
) {
	runtime_assert_string_msg( setting > 0.0, "Error in protocols::simple_filters::BuriedSurfaceAreaFilter::set_probe_radius(): The radius must be greater than 0 Angstroms." );
	probe_radius_ = setting;
}

/////// Private methods ///////////////

/// @brief Common function called by apply(), report(), and report_sm().
/// @details Does the actual computation.
bool
BuriedSurfaceAreaFilter::compute(
	std::ostream &os,
	core::pose::Pose const &pose,
	core::Real &buried_surf_area
) const {

	core::select::residue_selector::ResidueSubset selection( pose.total_residue(), true );
	compute_residue_selection( pose, selection ); //Uses the residue selector, the FAMILYVW selection option, and the requirement that we only select canonical L- and D-amino acids (and glycine).

	core::id::AtomID_Map< core::Real > atom_sasa; //Not used, but required by calc_per_atom_sasa function
	utility::vector1< core::Real > rsd_sasa; //Will store vector of residue SASAs.
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius(), false );

	utility::vector1< core::Real > rsd_buried_area( pose.total_residue(), 0.0 ); //Initialize a buried surface area vector (with entries for each residue) to zeroes.
	buried_surf_area = 0.0; //Used to calculate total.

	os << "RES\tBURIED AREA" << std::endl;

	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) { //Loop through all residues
		if ( !selection[i] ) {
			os << pose.residue_type(i).name3() << i << "\t0.0" << std::endl;
			continue; //Do nothing if a position is not selected.
		}

		char const oneletter_code( pose.residue_type(i).name1() );
		debug_assert( is_allowed_type( oneletter_code ) ); //Should be true.

		rsd_buried_area[i] = core::scoring::normalizing_area_total( oneletter_code ) - rsd_sasa[i];
		buried_surf_area += rsd_buried_area[i];
		os << pose.residue_type(i).name3() << i << "\t" << rsd_buried_area[i] << std::endl;
	}

	os << "TOTAL BURIED AREA = " << buried_surf_area << std::endl;
	bool const filter_passes( filter_out_low() ? buried_surf_area >= cutoff_buried_surface_area() : buried_surf_area <= cutoff_buried_surface_area() );

	os << "The BuriedSurfaceArea filter reports that this pose " << (filter_passes ? "passes" : "fails") << "." << std::endl;
	os.flush();

	return filter_passes;
}


/// @brief Is this a residue type for which we can calculate a total SASA value?
///
bool
BuriedSurfaceAreaFilter::is_allowed_type(
	char const oneletter_code
) const {
	return ( oneletter_code >= 'A' &&
		oneletter_code <= 'Y' &&
		oneletter_code != 'B' &&
		oneletter_code != 'J' &&
		oneletter_code != 'O' &&
		oneletter_code != 'U' &&
		oneletter_code != 'X') ;
}

/// @brief Is this a hydrophobic residue type (FAMILYVW)?
///
bool
BuriedSurfaceAreaFilter::is_FAMILYVW(
	char const oneletter_code
) const {
	return ( oneletter_code == 'F' ||
		oneletter_code == 'A' ||
		oneletter_code == 'M' ||
		oneletter_code == 'I' ||
		oneletter_code == 'L' ||
		oneletter_code == 'Y' ||
		oneletter_code == 'V' ||
		oneletter_code == 'W' ) ;
}

/// @brief Given a pose (input) and a ResidueSubset (output), compute the residues that should be
/// operated on.
void
BuriedSurfaceAreaFilter::compute_residue_selection(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset &selection
) const {

	debug_assert( pose.total_residue() == selection.size() ); //Should be true

	//Apply the residue selector, if there is one:
	if ( residue_selector() != nullptr ) {
		selection = residue_selector()->apply(pose);
	}

	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) { //Mask by allowed types:
		if ( !is_allowed_type(pose.residue_type(i).name1()) ) {
			selection[i] = false;
			continue;
		}
		if ( select_only_FAMILYVW() && !is_FAMILYVW( pose.residue_type(i).name1() ) ) selection[i] = false;
	}

}


/////////////// Creator ///////////////

protocols::filters::FilterOP
BuriedSurfaceAreaFilterCreator::create_filter() const
{
	return protocols::filters::FilterOP( new BuriedSurfaceAreaFilter );
}

std::string
BuriedSurfaceAreaFilterCreator::keyname() const
{
	return BuriedSurfaceAreaFilter::class_name();
}

void BuriedSurfaceAreaFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BuriedSurfaceAreaFilter::provide_xml_schema( xsd );
}

} //protocols
} //simple_filters
