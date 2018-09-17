// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/simple_metrics/SequenceRecoveryMetric.hh
/// @brief Calculate sequence recovery statistics on a protein, relative to a reference.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_analysis_simple_metrics_SequenceRecoveryMetric_HH
#define INCLUDED_protocols_analysis_simple_metrics_SequenceRecoveryMetric_HH

#include <protocols/analysis/simple_metrics/SequenceRecoveryMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/sequence/SequenceProfile.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace analysis {
namespace simple_metrics {

///@brief Calculate sequence recovery statistics on a protein, relative to a reference.
///
/// There's several options for how the sequence recovery is calculated, depending on what parameters are set.
/// Each metric is only calculated over the set of residues specified by the residue selector.
///
/// * Standard - Used when PSSM isn't set but reference pose is
///      This is a strict match/no-match fraction.
/// * Pass/Fail PSSM - Used when PSSM is set and use_ave_pssm is not (does not use reference pose)
///      This is the PSSM recovery metric from DeLuca, Dorr & Meiler 2011 Biochem 50(40):8521
///      Residue identities with positive (or zero) values in the PSSM count as a match, those with negative vales as no-match.
/// * Ave PSSM - Used when PSSM is set, use_ave_pssm is true, and no reference PDB is given
///      This value is the average of the values in the PSSM matrix for the residue identities.
/// * Delta PSSM - Used when PSSM is set, use_ave_pssm is true, and a reference PDB is provided.
///      This value is the average of the change in value of the PSSM matrix (mut - ref)
///
/// For PSSM metrics, it's assumed that the Pose numbering of both the main and reference structure
/// matches the numbering of the PSSM.

class SequenceRecoveryMetric : public core::simple_metrics::RealMetric {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	SequenceRecoveryMetric();

	/// @brief Copy constructor (not needed unless you need deep copies)
	SequenceRecoveryMetric( SequenceRecoveryMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~SequenceRecoveryMetric() override;

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

	///@brief Set a reference pose
	void
	set_comparison_pose( core::pose::PoseCOP pose);

	///@brief Set a residue selector to calculate the recovery from a subset of residues
	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP residue_selector );

	///@brief Set a residue selector (for the reference) to calculate the recovery from a subset of residues
	/// This needs to select the same number of residues as the set_residue_selector
	/// If not set (or nullptr), will use the same list as the selector passed to set_residue_selector()
	void
	set_residue_selector_ref( core::select::residue_selector::ResidueSelectorCOP residue_selector );

	///@brief Load the PSSM for reference from the given file.
	void
	load_pssm( std::string const & pssm_filename );

	///@brief Set if we want to use a average PSSM value metric, or the pass/fail setting.
	void
	set_use_ave_pssm( bool setting ) { use_ave_pssm_ = setting; }

	///@brief Name of the class
	std::string
	name() const override;

	///@brief Name of the class for creator.
	static
	std::string
	name_static();

	///@brief Name of the metric
	std::string
	metric() const override;

public:

	/// @brief called by parse_my_tag -- should not be used directly
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	core::simple_metrics::SimpleMetricOP
	clone() const override;

private:

	/// @brief The subsection of residues to use for sequence recovery
	core::select::residue_selector::ResidueSelectorCOP res_select_;

	/// @brief The reference pose to use for sequence recovery purposes
	core::pose::PoseCOP ref_pose_;
	/// @brief The subsection of residues to use for sequence recovery on the reference side
	/// If not set, use the same residues as in the passed pose.
	core::select::residue_selector::ResidueSelectorCOP res_select_ref_;
	/// @brief If set, do pssm recovery,
	core::sequence::SequenceProfileOP ref_profile_;
	/// @brief If true, do an averaged PSSM value, rather than a pass/fail one.
	bool use_ave_pssm_ = false;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //protocols
} //analysis
} //simple_metrics

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_analysis_simple_metrics_SequenceRecoveryMetric )
#endif // SERIALIZATION

#endif //protocols_analysis_simple_metrics_SequenceRecoveryMetric_HH





