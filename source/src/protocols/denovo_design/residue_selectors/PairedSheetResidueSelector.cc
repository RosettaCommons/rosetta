// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/PairedSheetResidueSelector.hh
/// @brief  Selects residues that are involved in strand-strand pairings
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/residue_selectors/PairedSheetResidueSelector.hh>
#include <protocols/denovo_design/residue_selectors/PairedSheetResidueSelectorCreator.hh>

#include <core/select/residue_selector/util.hh>

// Protocol headers
#include <protocols/denovo_design/components/SegmentPairing.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>

// Core headers
#include <core/scoring/dssp/Dssp.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions
#include <core/pose/Pose.hh>

// Basic Headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.residue_selectors.PairedSheetResidueSelector" );

namespace protocols {
namespace denovo_design {
namespace residue_selectors {

/// @brief Constructor.
///
PairedSheetResidueSelector::PairedSheetResidueSelector():
	ResidueSelector(),
	secstruct_( "" ),
	sheet_topology_( "" ),
	use_dssp_( true )
{
}

/// @brief Destructor.
///
PairedSheetResidueSelector::~PairedSheetResidueSelector() {}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
PairedSheetResidueSelector::ResidueSelectorOP
PairedSheetResidueSelector::clone() const
{
	return ResidueSelectorOP( new PairedSheetResidueSelector(*this) );
}

PairedSheetResidueSelector::ResidueSelectorOP
PairedSheetResidueSelectorCreator::create_residue_selector() const
{
	return PairedSheetResidueSelector::ResidueSelectorOP( new PairedSheetResidueSelector );
}

std::string
PairedSheetResidueSelector::get_name() const
{
	return PairedSheetResidueSelector::class_name();
}

std::string
PairedSheetResidueSelector::class_name()
{
	return "PairedSheetResidueSelector";
}

std::string
PairedSheetResidueSelectorCreator::keyname() const
{
	return PairedSheetResidueSelector::class_name();
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
PairedSheetResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & )
{
	secstruct_ = tag->getOption< std::string >( "secstruct", secstruct_ );
	sheet_topology_ = tag->getOption< std::string >( "sheet_topology", sheet_topology_ );
	use_dssp_ = tag->getOption< bool >( "use_dssp", use_dssp_ );
}

void
PairedSheetResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	//Syntax Example:
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute(  "secstruct",      xs_string, "Secondary structure, e.g. \"EEEELLEEEE\"."
		" If undefined, sec struct will be chosen based on the value of the 'use_dssp' option" )
		+ XMLSchemaAttribute(  "sheet_topology", xs_string, "String describing sheet topology, of the format A-B.P.R, "
		"where A is the strand number of the first strand in primary space, "
		"B is the strand number of the second strand in primary space, "
		"P\t is 'P' for parallel and 'A' for antiparallel, and R is the register shift. "
		"E.g. \"1-2.A.-1\"." )
		+ XMLSchemaAttribute::attribute_w_default(  "use_dssp",       xsct_rosetta_bool, "Use dssp to aito-detect secondary structure.",  "true"  );
	core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(), "Selects residues that are involved in strand-strand pairings.", attributes );
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
///
void
PairedSheetResidueSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	PairedSheetResidueSelector::provide_xml_schema( xsd );
}

/// @brief Sets sheet topology to be filtered
void
PairedSheetResidueSelector::set_sheet_topology( std::string const & sheet_topology )
{
	sheet_topology_ = sheet_topology;
}

/// @brief Sets secondary structure to be used in calculations
void
PairedSheetResidueSelector::set_secstruct( std::string const & secstruct )
{
	secstruct_ = secstruct;
}

/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
PairedSheetResidueSelector::ResidueSubset
PairedSheetResidueSelector::apply( core::pose::Pose const & pose ) const
{
	using protocols::fldsgn::topology::StrandPairingSet;
	using protocols::fldsgn::topology::StrandPairings;
	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;

	std::string const secstruct = get_secstruct( pose );
	std::string const sheet_topology = get_sheet_topology( pose );

	ResidueSubset subset( pose.size(), false );

	// if no pairings are specified, assume all E residues are paired
	if ( sheet_topology.empty() ) {
		core::Size resid = 1;
		for ( auto ss : secstruct ) {
			if ( ss == 'E' ) subset[ resid ] = true;
			++resid;
		}
		return subset;
	}

	TR.Debug << "Using sheet topology " << sheet_topology << std::endl;

	// now we check the strand pairing definitions in the BluePrint -- if some of them are impossible, DSSP will never return 'E',
	// so we will want these residues to be L/E (which is specified as 'h').
	SS_Info2_OP ss_info( new SS_Info2( pose, secstruct ) );
	StrandPairingSet const spairs( sheet_topology, ss_info );
	StrandPairings const pairs( spairs.strand_pairings() );

	// determine paired residues
	for ( StrandPairings::const_iterator pair=pairs.begin(); pair!=pairs.end(); ++pair ) {
		for ( core::Size s1_resid=(*pair)->begin1(); s1_resid<=(*pair)->end1(); ++s1_resid ) {
			if ( !(*pair)->has_paired_residue( s1_resid ) ) continue;
			core::Size const s2_resid = (*pair)->residue_pair( s1_resid );
			TR.Debug << "Pairing found: " << std::make_pair( s1_resid, s2_resid ) << std::endl;
			subset[s1_resid] = true;
			subset[s2_resid] = true;
		}
	}

	return subset;
}

/// @brief Gets secondary structure string to be used in computation
/// @param[in] pose input pose
/// @returns String representing secondary structure (e.g. "LEEEEELLEEEEEL")
/// @details If secstruct_ is given, that is returned.
///          If use_dssp_ is true, DSSP secstruct is returned
///          If use_dssp_ is false and secstruct_ is not given, return pose's secondary structure
std::string
PairedSheetResidueSelector::get_secstruct( core::pose::Pose const & pose ) const
{
	if ( !secstruct_.empty() ) {
		TR.Debug << "Using user-provided secstruct: " << secstruct_ << std::endl;
		return secstruct_;
	}
	if ( use_dssp_ ) {
		core::scoring::dssp::Dssp dssp( pose );
		TR.Debug << "Using DSSP secstruct: " << dssp.get_dssp_secstruct() << std::endl;
		return dssp.get_dssp_secstruct();
	}
	TR.Debug << "Using pose secstruct: " << pose.secstruct() << std::endl;
	return pose.secstruct();
}

/// @brief Gets sheet topology string to be used in computation
/// @param[in] pose input pose
/// @returns String representing sheet topology (e.g. "1-2.A.0;2-3.P.1")
/// @details If sheet_topology_ is given, that is returned.
///          If sheet_topology_ is not given, and the pose has a cached StructureData,
///          the topology string is derived from that.
///          If a sheet topology cannot be found by the above methods, an error is
///          thrown
std::string
PairedSheetResidueSelector::get_sheet_topology( core::pose::Pose const & pose ) const
{
	if ( !sheet_topology_.empty() ) return sheet_topology_;
	components::StructureDataFactory const & factory = *components::StructureDataFactory::get_instance();
	if ( factory.has_cached_data( pose ) ) {
		return components::SegmentPairing::get_strand_pairings( factory.get_from_const_pose( pose ) );
	}

	TR << "Could not determine strand pairings! "
		<< "You must specify them using the \"sheet_topology\" option or attach a StructureData object to the pose. "
		<< "No residues will be selected."
		<< std::endl;
	return "";
}

} //protocols
} //denovo_design
} //residue_selectors
