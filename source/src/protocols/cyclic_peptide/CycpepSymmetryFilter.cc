// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/CycpepSymmetryFilter.cc
/// @brief A filter that examines a cyclic peptide's structure and returns TRUE if and only if it has a desired backbone symmetry.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#include <protocols/cyclic_peptide/CycpepSymmetryFilter.hh>
#include <protocols/cyclic_peptide/CycpepSymmetryFilterCreator.hh>

//Core headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

//Protocols headers
#include <protocols/rosetta_scripts/util.hh>

//XSD includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

//Basic and utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>

//Numeric headers
#include <numeric/angle.functions.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.cyclic_peptide.CycpepSymmetryFilter" );

namespace protocols {
namespace cyclic_peptide {

CycpepSymmetryFilter::CycpepSymmetryFilter():
	protocols::filters::Filter( "CycpepSymmetryFilter" ),
	symm_repeats_(2),
	mirror_symm_(false),
	selector_(),
	angle_threshold_(10.0)
	//TODO -- initialize variables here.
{

}

CycpepSymmetryFilter::~CycpepSymmetryFilter()
{}

void
CycpepSymmetryFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	using namespace core::select::residue_selector;
	ResidueSelectorCOP selector( protocols::rosetta_scripts::parse_residue_selector( tag, data ) );
	if ( selector != nullptr ) { set_selector(selector); }
	set_symm_repeats( tag->getOption<core::Size>("symmetry_repeats", symm_repeats() ) );
	set_mirror_symm( tag->getOption<bool>("mirror_symmetry", mirror_symm()) );
	set_angle_threshold( tag->getOption<core::Real>("angle_threshold", angle_threshold() ) );
}

/// @brief Sets the repeats in the symmetry that we're looking for
/// (e.g. 2 for c2 or s2 symmetry, 3 for c3, etc.).
void
CycpepSymmetryFilter::set_symm_repeats(
	core::Size const repeats_in
) {
	runtime_assert_string_msg( repeats_in > 1, "Error in protocols::cyclic_peptide::CycpepSymmetryFilter::set_symm_repeats(): The number of repeats must be greater than one.");
	symm_repeats_ = repeats_in;
}

/// @brief Sets the residue selector.
/// @details Note that the selector is NOT cloned.
void
CycpepSymmetryFilter::set_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	runtime_assert_string_msg( selector_in != nullptr, "Error in protocols::cyclic_peptide::CycpepSymmetryFilter::set_selector(): A null pointer was passed to this function!" );
	selector_ = selector_in;
}

/// @brief Set the cutoff, in degrees, that two mainchain dihedral values must lie within in order for two residues to
/// be considered to have the "same" value for that mainchain degree of freedom.
void
CycpepSymmetryFilter::set_angle_threshold(
	core::Real const &setting
) {
	runtime_assert_string_msg( setting > 0.0, "Error in protocols::cyclic_peptide::CycpepSymmetryFilter::set_angle_threshold(): A value greater than 0 is required." );
	angle_threshold_ = setting;
}


protocols::filters::FilterOP
CycpepSymmetryFilter::clone() const
{
	return protocols::filters::FilterOP( new CycpepSymmetryFilter( *this ) );
}


protocols::filters::FilterOP
CycpepSymmetryFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new CycpepSymmetryFilter );
}

bool
CycpepSymmetryFilter::apply( core::pose::Pose const &pose ) const
{
	if ( mirror_symm() ) {
		runtime_assert_string_msg( symm_repeats() % 2 == 0, "Error in protocols::cyclic_peptide::CycpepSymmetryFilter::apply(): Mirror symmetry was set, but the number of symmetry repeats is not even." );
	}

	utility::vector1< core::Size > reslist;
	get_residues( pose, selector(), reslist ); //Get the residues of the peptide.
	runtime_assert_string_msg( is_cyclic_peptide( pose, reslist ), "Error in protocols::cyclic_peptide::CycpepSymmetryFilter::apply(): The pose or selected residues do NOT represent a cyclic peptide." ); //Check that these residues represent a cyclic peptide.

	if ( reslist.size() % symm_repeats() != 0 ) {
		TR.Warning << "The CycpepSymmetryFilter was applied with " << (mirror_symm() ? "s" : "c")
			<< symm_repeats() << " symmetry to a " << (selector() ? "selection" : "pose" )
			<< " with " << reslist.size() << " residues.  The number of residues is not an integer multiple of the symmetry repeats.  Pose MUST be asymmetric.  Returning FALSE."
			<< std::endl;
		return false;
	}

	bool symm(true);

	for ( core::Size i=1, imax=( reslist.size() / symm_repeats() /*Already checked above that this is evenly divisible.*/ ); i<=imax; ++i ) {
		core::Size const res1( reslist[i] );

		for ( core::Size j=2; j<=symm_repeats(); ++j ) {
			core::Size const res2( reslist[(j-1)*imax+i] );
			bool const flip( mirror_symm() && j % 2 == 0 ); //Is this a mirror position?

			if (
					mainchain_torsions_differ( pose, res1, res2, flip )
					) {
				symm=false;
				break;
			}
		}
		if ( !symm ) break;
	}

	return symm;
}

core::Real
CycpepSymmetryFilter::report_sm( core::pose::Pose const &pose ) const
{
	return (apply(pose) ? 1.0 : 0.0);
}

void
CycpepSymmetryFilter::report( std::ostream &outstream, core::pose::Pose const &pose ) const
{
	outstream << "This peptide ";
	if ( apply(pose) ) {
		outstream << "is";
	} else {
		outstream << "is not";
	}
	outstream << (mirror_symm() ? " s" : " c") << symm_repeats() << " symmetric." << std::endl;
}

std::string CycpepSymmetryFilter::name() const {
	return class_name();
}

std::string CycpepSymmetryFilter::class_name() {
	return "CycpepSymmetryFilter";
}

void CycpepSymmetryFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute( "symmetry_repeats", xsct_non_negative_integer, "The number of repeats in this type of symmetry.  For example, for c3 symmetry, one would provide \"3\" as input.  Defaults to 2 (for c2 symmetry)." )
		+ XMLSchemaAttribute( "mirror_symmetry", xsct_rosetta_bool, "Is this a type of cyclic symmetry with mirror operations, such as s2, s4, or s6 symmetry?  If true, symmetry_repeats must be a multiple of 2.  Defaults to false (for c2 symmetry -- i.e. not s2 symmetry)." )
		+ XMLSchemaAttribute( "angle_threshold", xsct_real, "The cutoff, in degrees, for the difference between two dihedral angles before they are considered \"different\" angles.  This is used when comparing mainchain torsion values of differet residues.  Defaults to 10.0 degrees." )
		;

	core::select::residue_selector::attributes_for_parse_residue_selector_default_option_name( attlist, "An optional residue selector set to select the cyclic peptide.  If not provided, the whole pose is used." );

	protocols::filters::xsd_type_definition_w_attributes(
		xsd,
		class_name(),
		"The CycpepSymmetryFilter examines a cyclic peptide for internal backbone symmetry matching a particular cyclic symmetry type (e.g. c2, c3, c4, etc.) or mirror symmetry type (e.g. s2, s4, etc.), and returns true if and only if the backbone is symmetric.",
		attlist );
}

/****************************************
PRIVATE FUNCTIONS
*****************************************/

/// @brief Given a pose and an optional ResidueSelector, get a list of residues.
void
CycpepSymmetryFilter::get_residues(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSelectorCOP selector,
	utility::vector1<core::Size> &reslist_out
) const {
	reslist_out.clear();
	if ( selector == nullptr ) {
		reslist_out.resize(pose.total_residue());
		for ( core::Size i=1, imax=pose.total_residue(); i<=imax; ++i ) {
			reslist_out[i] = i;
		}
	} else {
		core::select::residue_selector::ResidueSubset const subset( selector->apply(pose) );
		for ( core::Size i=1, imax=subset.size(); i<=imax; ++i ) {
			if ( subset[i] ) reslist_out.push_back(i);
		}
	}
}


/// @brief Checks that the geometry to which we're applying this filter is actually a
/// cyclic peptide.
bool
CycpepSymmetryFilter::is_cyclic_peptide(
	core::pose::Pose const &pose,
	utility::vector1< core::Size > const &residues
) const {

	for ( core::Size i=1, imax=residues.size(); i<=imax; ++i ) {
		core::Size const res( residues[i] );
		core::conformation::Residue const &curres( pose.residue(res) );
		if ( !curres.has_lower_connect() || !curres.has_upper_connect() ) return false;
		core::chemical::ResidueType const &restype( pose.residue_type(res) );
		if ( !restype.is_alpha_aa() && !restype.is_beta_aa() && !restype.is_gamma_aa() && !restype.is_peptoid() && !restype.is_oligourea() ) return false;

		if ( i == 1 ) {
			if ( curres.connected_residue_at_resconn( restype.lower_connect_id() ) != residues[imax] ) return false;
		} else { //if( i > 1 )
			if ( curres.connected_residue_at_resconn( restype.lower_connect_id() ) != residues[i-1] ) return false;
		}

		if ( i == imax ) {
			if ( curres.connected_residue_at_resconn( restype.upper_connect_id() ) != residues[1] ) return false;
		} else { //if( i < imax )
			if ( curres.connected_residue_at_resconn( restype.upper_connect_id() ) != residues[i+1] ) return false;
		}
	}

	return true;
}

/// @brief Returns TRUE if the number of mainchain torsions differs, or if the values differ; FALSE otherwise.
/// @brief If flip is true, then torsion values are inverted prior to comparison.
bool
CycpepSymmetryFilter::mainchain_torsions_differ(
	core::pose::Pose const &pose,
	core::Size const res1,
	core::Size const res2,
	bool const flip
) const {
	using namespace numeric;

	core::Size const res1_ntors(pose.residue(res1).mainchain_torsions().size());
	core::Size const res2_ntors(pose.residue(res2).mainchain_torsions().size());
	if ( res1_ntors != res2_ntors ) return true;

	core::Real const flip_multiplier( flip ? -1.0 : 1.0 );

	for ( core::Size i=1, imax=res1_ntors; i<=imax; ++i ) {
		core::Real const val1( nonnegative_principal_angle_degrees ( pose.residue(res1).mainchain_torsions()[i] ) );
		core::Real const val2( nonnegative_principal_angle_degrees ( flip_multiplier * pose.residue(res2).mainchain_torsions()[i] ) );

		core::Real diff( std::abs(val1-val2) );
		if ( diff > 180.0 ) diff = 360.0 - diff;

		if ( diff > angle_threshold() ) return true;
	}

	return false;
}


/////////////// Creator ///////////////

protocols::filters::FilterOP
CycpepSymmetryFilterCreator::create_filter() const
{
	return protocols::filters::FilterOP( new CycpepSymmetryFilter );
}

std::string
CycpepSymmetryFilterCreator::keyname() const
{
	return CycpepSymmetryFilter::class_name();
}

void CycpepSymmetryFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CycpepSymmetryFilter::provide_xml_schema( xsd );
}

} //protocols
} //cyclic_peptide
