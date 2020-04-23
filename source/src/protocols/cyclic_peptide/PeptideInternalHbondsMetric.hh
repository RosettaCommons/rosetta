// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/PeptideInternalHbondsMetric.hh
/// @brief A simple metric for measuring the hydrogen bonds within a peptide backbone, ignoring hydrogen
/// bonds to other parts of a pose.  Note that this metric does not count intra-residue hydrogen bonds.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_cyclic_peptide_PeptideInternalHbondsMetric_HH
#define INCLUDED_protocols_cyclic_peptide_PeptideInternalHbondsMetric_HH

#include <protocols/cyclic_peptide/PeptideInternalHbondsMetric.fwd.hh>
#include <protocols/cyclic_peptide/PeptideInternalHbondsMetricTests.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace cyclic_peptide {

/// @brief A simple metric for measuring the backbone hydrogen bonds within a peptide backbone,
/// ignoring hydrogen bonds to other parts of a pose.  Note that this metric does not count
/// intra-residue hydrogen bonds.
class PeptideInternalHbondsMetric : public core::simple_metrics::RealMetric{

	// Allow friendship with test classes to allow private functions to be tested.
	friend class ::PeptideInternalHbondsMetricTests;

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	PeptideInternalHbondsMetric() = default;

	/// @brief Copy constructor (not needed unless you need deep copies)
	PeptideInternalHbondsMetric( PeptideInternalHbondsMetric const & ) = default;

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~PeptideInternalHbondsMetric() override = default;

	/// @brief Duplicate this object, and return a smart pointer to the copy.
	core::simple_metrics::SimpleMetricOP clone() const override;

public:

	/////////////////////
	/// Metric Methods ///
	/////////////////////

	///Defined in RealMetric:
	///
	/// @brief Calculate the metric and add it to the pose as a score.
	///           labeled as prefix+metric+suffix.
	///
	/// @details Score is added through setExtraScorePose and is output
	///            into the score tables/file at pose output.
	//void
	//apply( core::pose::Pose & pose, prefix="", suffix="" ) override;

	///@brief Calculate the metric.
	core::Real
	calculate( core::pose::Pose const & pose ) const override;

public:

	///@brief Name of the class
	std::string
	name() const override;

	///@brief Name of the class for creator.
	static
	std::string
	name_static();

	/// @brief Return text indicating what we're measuring.
	std::string
	metric() const override;

public:

	/// @brief called by parse_my_tag -- should not be used directly
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	/// @brief Configure this metric from the tag passed to the companion filter (the
	/// PeptideInternalHbondsFilter).
	void
	configure_from_filter_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	);

	/// @brief Provide a description of the XML interface.
	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Provide descriptions of XML interface elements shared between this metric and the companion filter.
	static
	void
	provide_shared_xml_schema_elements(
		utility::tag::AttributeList & attlist
	);

	/// @brief allows printing data to a stream.
	void
	report(
		std::ostream &outstream,
		core::pose::Pose const &pose
	) const;

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

public: //Getters:

	/// @brief Get whether we're counting backbone-backbone hbonds:
	inline bool backbone_backbone() const { return backbone_backbone_; }

	/// @brief Get whether we're counting backbone-sidechain hbonds:
	inline bool backbone_sidechain() const { return backbone_sidechain_; }

	/// @brief Get whether we're counting sidechain-sidechain hbonds:
	inline bool sidechain_sidechain() const { return sidechain_sidechain_; }

	/// @brief Get the number of residues apart in terms of covalent connectivity that two residues have to be in order for their
	/// hydrogen bonds to be counted.
	/// @details Defaults to 1.
	inline core::Size exclusion_distance() const { return exclusion_distance_; }

	/// @brief Get the energy cutoff for counting a hydrogen bond.
	/// @details Defaults to -0.25.
	inline core::Real hbond_energy_cutoff() const { return hbond_energy_cutoff_; }

	/// @brief Get the scorefunction for identifying hydrogen bonds.  Note that if this is nullptr, the
	/// default scorefunction is used.
	inline core::scoring::ScoreFunctionCOP scorefxn() const { return scorefxn_; }

	/// @brief Get the scorefunction for selecting the peptide.  If this nullptr, whole pose is counted.
	inline core::select::residue_selector::ResidueSelectorCOP residue_selector() const { return residue_selector_; }

public: //Functions needed for the citation manager

	/// @brief Does this simple metric provide information about how to cite it?
	/// @details Returns false.
	bool simple_metric_provides_citation_info() const override;

	/// @brief Provide the citation.
	/// @returns An empty vector for this simple metric, since it's unpublished.  Provides citations for the scorefunction
	/// and residue selector, if available, however.
	utility::vector1< basic::citation_manager::CitationCollectionCOP > provide_citation_info() const override;

	/// @brief Does this simple metric indicate that it is unpublished (and, by extension, that the author should be
	/// included in publications resulting from it)?
	/// @returns True, since this is unpublished.
	bool simple_metric_is_unpublished() const override;

	/// @brief Provide a list of authors and their e-mail addresses, as strings.
	/// @details Also provides unpublished author information for residue selectors and scorefunctions used.
	utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP > provide_authorship_info_for_unpublished() const override;

private: //Methods:

	/// @brief Given a pose and a selection, count the internal hydrogen bonds.
	core::Size
	count_hbonds(
		core::pose::Pose const & original_pose,
		core::select::residue_selector::ResidueSubset const & selection
	) const;

	/// @brief Given a vector of selected resiudes and a pose, construct a vector of lists of allowed residue partners.
	/// @details This is an O(N^2) operation, unfortunately.  However, the selected_residues list should be small, and N
	/// is its length.
	utility::vector1< std::list < core::Size > >
	generate_allowed_partners(
		core::pose::Pose const & pose,
		utility::vector1< core::Size > const & selected_residues
	) const;

private: //Data:

	/// @brief Are we counting backbone-backbone hydrogen bonds?  Default true.
	bool backbone_backbone_ = true;

	/// @brief Are we counting backbone-sidechain hydrogen bonds?  Default false.
	bool backbone_sidechain_ = false;

	/// @brief Are we counting sidechain-sidechain hydrogen bonds?  Default false.
	bool sidechain_sidechain_ = false;

	/// @brief The number of residues apart in terms of covalent connectivity that two residues have to be in order for their
	/// hydrogen bonds to be counted.  Defaults to 1.
	core::Size exclusion_distance_ = 1;

	/// @brief The energy cutoff for counting a hydrogen bond.  Defaults to -0.25.
	core::Real hbond_energy_cutoff_ = -0.25;

	/// @brief The scorefunction to use for evaluating hydrogen bonds.  If none is provided, we'll grab the
	/// default scorefunction and cache it.
	mutable core::scoring::ScoreFunctionCOP scorefxn_;

	/// @brief The residue selector.  If none is provided, we evaluate the whole pose.
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //cyclic_peptide
} //protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_cyclic_peptide_PeptideInternalHbondsMetric )
#endif // SERIALIZATION

#endif //protocols_cyclic_peptide_PeptideInternalHbondsMetric_HH





