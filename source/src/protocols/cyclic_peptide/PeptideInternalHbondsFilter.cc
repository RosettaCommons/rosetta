// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/PeptideInternalHbondsFilter.cc
/// @brief A filter that thinly wraps the PeptideInternalHbondsMetric.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#include <protocols/cyclic_peptide/PeptideInternalHbondsFilter.hh>
#include <protocols/cyclic_peptide/PeptideInternalHbondsFilterCreator.hh>

//Protocols inccludes
#include <protocols/filters/filter_schemas.hh>
#include <protocols/cyclic_peptide/PeptideInternalHbondsMetric.hh>

//Core includes
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/xml_util.hh>

//Basic includes
#include <basic/Tracer.hh>
#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

//Utility includes
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.cyclic_peptide.PeptideInternalHbondsFilter" );

namespace protocols {
namespace cyclic_peptide {

/// @brief Default constructor.
PeptideInternalHbondsFilter::PeptideInternalHbondsFilter() :
	protocols::filters::Filter( "PeptideInternalHbondsFilter" ),
	hbond_metric_(utility::pointer::make_shared< PeptideInternalHbondsMetric >())
{}

/// @brief Required in the context of the parser/scripting scheme.
protocols::filters::FilterOP
PeptideInternalHbondsFilter::fresh_instance() const
{
	return utility::pointer::make_shared< PeptideInternalHbondsFilter >();
}

/// @brief Make a copy of this filter and return an owning pointer to the copy.
protocols::filters::FilterOP
PeptideInternalHbondsFilter::clone() const
{
	return utility::pointer::make_shared< PeptideInternalHbondsFilter >( *this );
}

/// @brief returns true if the structure has internal hbonds equal to or greater than
/// the cutoff, false otherwise.
bool
PeptideInternalHbondsFilter::apply( core::pose::Pose const &pose ) const
{
	debug_assert( hbond_metric_!=nullptr );
	core::Real const hbond_count( hbond_metric_->calculate(pose) );
	bool const passed( hbond_count >= hbond_cutoff_ );
	TR << "Measured " << hbond_count << " hydrogen bonds; cutoff is " << hbond_cutoff_ << ".  Filter " << (passed ? "PASSED." : "FAILED.") << std::endl;
	return passed;
}

/// @brief Required for reporting score values.
core::Real
PeptideInternalHbondsFilter::report_sm( core::pose::Pose const &pose ) const
{
	debug_assert( hbond_metric_!=nullptr );
	return hbond_metric_->calculate( pose );
}

/// @brief allows printing data to a stream.
void
PeptideInternalHbondsFilter::report(
	std::ostream &outstream,
	core::pose::Pose const &pose
) const {
	hbond_metric_->report( outstream, pose );
}

std::string PeptideInternalHbondsFilter::name() const {
	return class_name();
}

std::string PeptideInternalHbondsFilter::class_name() {
	return "PeptideInternalHbondsFilter";
}

/// @brief Parse XML tag (to use this Filter in Rosetta Scripts).
void
PeptideInternalHbondsFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &datamap
) {
	debug_assert( hbond_metric_!=nullptr );
	hbond_metric_->configure_from_filter_tag( tag, datamap );

	//Add any filter-specific options here:
	set_hbond_cutoff( tag->getOption<core::Size>( "hbond_cutoff", hbond_cutoff() ) );
}

/// @brief Describe the XML interface in a machine-readable way.
void
PeptideInternalHbondsFilter::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) {
	using namespace utility::tag;
	AttributeList attlist;

	PeptideInternalHbondsMetric::provide_shared_xml_schema_elements( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default( "hbond_cutoff", xsct_non_negative_integer, "The threshold number of internal hydrogen bonds, below which the filter fails.  Defaults to 1.", "1");
	protocols::filters::xsd_type_definition_w_attributes(
		xsd,
		class_name(),
		"This filter filters based on a user-specified minimum number of internal backbone hydrogen bonds within the selection or pose.  Note that this filter is actually a thin wrapper for the PeptideInternalHbondsMetric.",
		attlist );
}

//////// Setters ////////

/// @brief Set whether we're counting backbone-backbone, backbone-sidechain, and sidechain-sidechain hydrogen bonds.
/// @details Throws if all are set to false.  Defaults to only backbone-backbone being true.
void
PeptideInternalHbondsFilter::set_hbond_types(
	bool const backbone_backbone_setting,
	bool const backbone_sidechain_setting,
	bool const sidechain_sidechain_setting
) {
	debug_assert( hbond_metric_ != nullptr );
	hbond_metric_->set_hbond_types( backbone_backbone_setting, backbone_sidechain_setting, sidechain_sidechain_setting );
}

/// @brief Set the number of residues apart in terms of covalent connectivity that two residues have to be in order for their
/// hydrogen bonds to be counted.
/// @details Defaults to 1.
void
PeptideInternalHbondsFilter::set_exclusion_distance(
	core::Size const setting
) {
	debug_assert( hbond_metric_ != nullptr );
	hbond_metric_->set_exclusion_distance( setting );
}

/// @brief Set the energy cutoff for counting a hydrogen bond.
/// @details Defaults to -0.25.
void
PeptideInternalHbondsFilter::set_hbond_energy_cutoff(
	core::Real const setting
) {
	debug_assert( hbond_metric_ != nullptr );
	hbond_metric_->set_hbond_energy_cutoff( setting );
}

/// @brief Set the scorefunction for identifying hydrogen bonds.  If this is not called, the default scorefunction is used.
void
PeptideInternalHbondsFilter::set_scorefxn(
	core::scoring::ScoreFunctionCOP const & sfxn_in
) {
	debug_assert( hbond_metric_ != nullptr );
	hbond_metric_->set_scorefxn( sfxn_in );
}

/// @brief Set the scorefunction for selecting the peptide.  If this is not called, whole pose is counted.
void
PeptideInternalHbondsFilter::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP const & res_selector_in
) {
	debug_assert( hbond_metric_ != nullptr );
	hbond_metric_->set_residue_selector( res_selector_in );
}

//////// Getters ////////

/// @brief Get whether we're counting backbone-backbone hbonds:
bool
PeptideInternalHbondsFilter::backbone_backbone() const {
	debug_assert( hbond_metric_ != nullptr );
	return hbond_metric_->backbone_backbone();
}

/// @brief Get whether we're counting backbone-sidechain hbonds:
bool
PeptideInternalHbondsFilter::backbone_sidechain() const {
	debug_assert( hbond_metric_ != nullptr );
	return hbond_metric_->backbone_sidechain();
}

/// @brief Get whether we're counting sidechain-sidechain hbonds:
bool
PeptideInternalHbondsFilter::sidechain_sidechain() const {
	debug_assert( hbond_metric_ != nullptr );
	return hbond_metric_->sidechain_sidechain();
}

/// @brief Get the number of residues apart in terms of covalent connectivity that two residues have to be in order for their
/// hydrogen bonds to be counted.
/// @details Defaults to 1.
core::Size
PeptideInternalHbondsFilter::exclusion_distance() const {
	debug_assert( hbond_metric_ != nullptr );
	return hbond_metric_->exclusion_distance();
}

/// @brief Get the energy cutoff for counting a hydrogen bond.
/// @details Defaults to -0.25.
core::Real
PeptideInternalHbondsFilter::hbond_energy_cutoff() const {
	debug_assert( hbond_metric_ != nullptr );
	return hbond_metric_->hbond_energy_cutoff();
}

/// @brief Get the scorefunction for identifying hydrogen bonds.  Note that if this is nullptr, the
/// default scorefunction is used.
core::scoring::ScoreFunctionCOP
PeptideInternalHbondsFilter::scorefxn() const {
	debug_assert( hbond_metric_ != nullptr );
	return hbond_metric_->scorefxn();
}

/// @brief Get the scorefunction for selecting the peptide.  If this nullptr, whole pose is counted.
core::select::residue_selector::ResidueSelectorCOP
PeptideInternalHbondsFilter::residue_selector() const {
	debug_assert( hbond_metric_ != nullptr );
	return hbond_metric_->residue_selector();
}

////////// CITATION MANAGER FUNCTIONS //////////

/// @brief Does this filter provide information about how to cite it?
/// @details Returns false.
bool
PeptideInternalHbondsFilter::filter_provides_citation_info() const {
	return false;
}

/// @brief Provide the citation.
/// @returns An empty vector for this filter, since it's unpublished.  Provides citations for the scorefunction
/// and residue selector, if available, however.
utility::vector1< basic::citation_manager::CitationCollectionCOP >
PeptideInternalHbondsFilter::provide_citation_info() const {
	utility::vector1< basic::citation_manager::CitationCollectionCOP > returnvec;
	debug_assert( hbond_metric_ != nullptr );
	core::select::residue_selector::ResidueSelectorCOP selector( hbond_metric_->residue_selector() );
	core::scoring::ScoreFunctionCOP sfxn( hbond_metric_->scorefxn() );

	if ( selector != nullptr ) {
		basic::citation_manager::merge_into_citation_collection_vector( selector->provide_citation_info(), returnvec );
	}
	if ( sfxn != nullptr ) {
		basic::citation_manager::merge_into_citation_collection_vector( sfxn->provide_citation_info(), returnvec );
	}
	return returnvec;
}

/// @brief Does this filter indicate that it is unpublished (and, by extension, that the author should be
/// included in publications resulting from it)?
/// @returns True, since this is unpublished.
bool
PeptideInternalHbondsFilter::filter_is_unpublished() const {
	return true;
}

/// @brief Provide a list of authors and their e-mail addresses, as strings.
/// @details Also provides unpublished author information for residue selectors and scorefunctions used.
utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >
PeptideInternalHbondsFilter::provide_authorship_info_for_unpublished() const {
	utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP > returnvec{
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		class_name(),
		basic::citation_manager::CitedModuleType::Filter,
		"Vikram K. Mulligan",
		"Systems Biology, Center for Computational Biology, Flatiron Institute",
		"vmulligan@flatironinstitute.org"
		)
		};

	debug_assert( hbond_metric_ != nullptr );
	core::select::residue_selector::ResidueSelectorCOP selector( hbond_metric_->residue_selector() );
	core::scoring::ScoreFunctionCOP sfxn( hbond_metric_->scorefxn() );

	if ( selector != nullptr ) {
		basic::citation_manager::merge_into_unpublished_collection_vector( selector->provide_authorship_info_for_unpublished(), returnvec );
	}
	if ( sfxn != nullptr ) {
		basic::citation_manager::merge_into_unpublished_collection_vector( sfxn->provide_authorship_info_for_unpublished(), returnvec );
	}
	return returnvec;
}

/////////////// Creator ///////////////

protocols::filters::FilterOP
PeptideInternalHbondsFilterCreator::create_filter() const
{
	return utility::pointer::make_shared< PeptideInternalHbondsFilter >( );
}

std::string
PeptideInternalHbondsFilterCreator::keyname() const
{
	return PeptideInternalHbondsFilter::class_name();
}

void PeptideInternalHbondsFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PeptideInternalHbondsFilter::provide_xml_schema( xsd );
}

} //cyclic_peptide
} //protocols
