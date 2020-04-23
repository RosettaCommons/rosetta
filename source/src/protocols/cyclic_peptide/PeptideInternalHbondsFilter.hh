// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/PeptideInternalHbondsFilter.hh
/// @brief A filter that thinly wraps the PeptideInternalHbondsMetric.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_cyclic_peptide_PeptideInternalHbondsFilter_hh
#define INCLUDED_protocols_cyclic_peptide_PeptideInternalHbondsFilter_hh

// Unit headers
#include <protocols/cyclic_peptide/PeptideInternalHbondsFilter.fwd.hh>
#include <protocols/cyclic_peptide/PeptideInternalHbondsMetric.fwd.hh>
#include <protocols/filters/Filter.hh>

// Protocols headers

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace cyclic_peptide {

///@brief A filter that thinly wraps the PeptideInternalHbondsMetric.
class PeptideInternalHbondsFilter : public protocols::filters::Filter {

public:

	/// @brief Default constructor.
	PeptideInternalHbondsFilter();

	/// @brief Default destructor.
	~PeptideInternalHbondsFilter() override = default;

	/// @brief Required in the context of the parser/scripting scheme.
	protocols::filters::FilterOP
	fresh_instance() const override;

	/// @brief Make a copy of this filter and return an owning pointer to the copy.
	protocols::filters::FilterOP
	clone() const override;

public:

	/// @brief returns true if the structure has internal hbonds equal to or greater than
	/// the cutoff, false otherwise.
	bool
	apply( core::pose::Pose const & pose ) const override;

	/// @brief Required for reporting score values.
	core::Real
	report_sm( core::pose::Pose const & pose ) const override;

	/// @brief allows printing data to a stream.
	void
	report( std::ostream & os, core::pose::Pose const & pose ) const override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	/// @brief Parse XML tag (to use this Filter in Rosetta Scripts).
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	/// @brief Describe the XML interface in a machine-readable way.
	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public: //Setters:

	/// @brief Set whether we're counting backbone-backbone, backbone-sidechain, and sidechain-sidechain hydrogen bonds.
	/// @details Throws if all are set to false.  Defaults to only backbone-backbone being true.
	void
	set_hbond_types(
		bool const backbone_backbone_setting,
		bool const backbone_sidechain_setting,
		bool const sidechain_sidechain_setting
	);

	/// @brief Set the number of residues apart in terms of covalent connectivity that two residues have to be in order for their
	/// hydrogen bonds to be counted.
	/// @details Defaults to 1.
	void
	set_exclusion_distance(
		core::Size const setting
	);

	/// @brief Set the energy cutoff for counting a hydrogen bond.
	/// @details Defaults to -0.25.
	void
	set_hbond_energy_cutoff(
		core::Real const setting
	);

	/// @brief Set the scorefunction for identifying hydrogen bonds.  If this is not called, the default scorefunction is used.
	void
	set_scorefxn(
		core::scoring::ScoreFunctionCOP const & sfxn_in
	);

	/// @brief Set the scorefunction for selecting the peptide.  If this is not called, whole pose is counted.
	void
	set_residue_selector(
		core::select::residue_selector::ResidueSelectorCOP const & res_selector_in
	);

	/// @brief Set the threshold number of internal hydrogen bonds, below which the filter fails.
	/// @details Defaults to 1.
	inline void set_hbond_cutoff( core::Size const setting ) { hbond_cutoff_ = setting; }

public: //Getters:

	/// @brief Get whether we're counting backbone-backbone hbonds:
	bool backbone_backbone() const;

	/// @brief Get whether we're counting backbone-sidechain hbonds:
	bool backbone_sidechain() const;

	/// @brief Get whether we're counting sidechain-sidechain hbonds:
	bool sidechain_sidechain() const;

	/// @brief Get the number of residues apart in terms of covalent connectivity that two residues have to be in order for their
	/// hydrogen bonds to be counted.
	/// @details Defaults to 1.
	core::Size exclusion_distance() const;

	/// @brief Get the energy cutoff for counting a hydrogen bond.
	/// @details Defaults to -0.25.
	core::Real hbond_energy_cutoff() const;

	/// @brief Get the scorefunction for identifying hydrogen bonds.  Note that if this is nullptr, the
	/// default scorefunction is used.
	core::scoring::ScoreFunctionCOP scorefxn() const;

	/// @brief Get the scorefunction for selecting the peptide.  If this nullptr, whole pose is counted.
	core::select::residue_selector::ResidueSelectorCOP residue_selector() const;

	/// @brief Get the threshold number of internal hydrogen bonds, below which the filter fails.
	/// @details Defaults to 1.
	inline core::Size hbond_cutoff() const { return hbond_cutoff_; }

public: //Accessors

	/// @brief Access the internal PeptideInternalHbondsMetric used by this filter (const access).
	inline PeptideInternalHbondsMetric const & hbond_metric() const { return *hbond_metric_; }

	/// @brief Access the internal PeptideInternalHbondsMetric used by this filter (non-const access).
	inline PeptideInternalHbondsMetric & hbond_metric() { return *hbond_metric_; }

public: //Functions needed for the citation manager

	/// @brief Does this filter provide information about how to cite it?
	/// @details Returns false.
	bool filter_provides_citation_info() const override;

	/// @brief Provide the citation.
	/// @returns An empty vector for this filter, since it's unpublished.  Provides citations for the scorefunction
	/// and residue selector, if available, however.
	utility::vector1< basic::citation_manager::CitationCollectionCOP > provide_citation_info() const override;

	/// @brief Does this filter indicate that it is unpublished (and, by extension, that the author should be
	/// included in publications resulting from it)?
	/// @returns True, since this is unpublished.
	bool filter_is_unpublished() const override;

	/// @brief Provide a list of authors and their e-mail addresses, as strings.
	/// @details Also provides unpublished author information for residue selectors and scorefunctions used.
	utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP > provide_authorship_info_for_unpublished() const override;

private: // Data:

	/// @brief The simple metric that we wrap.
	PeptideInternalHbondsMetricOP hbond_metric_;

	/// @brief The threshold number of internal hydrogen bonds, below which the filter fails.
	/// @details Defaults to 1.
	core::Size hbond_cutoff_ = 1;

};

} //cyclic_peptide
} //protocols

#endif //INCLUDED_protocols_cyclic_peptide_PeptideInternalHbondsFilter_hh
