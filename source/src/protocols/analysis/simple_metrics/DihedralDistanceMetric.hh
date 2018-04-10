// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/DihedralDistanceMetric.hh
/// @brief A metric to calculate the dihedral distance between two poses or the input and the set cmd-line native.  Can set a subset of residues to calculate via ResidueSelector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_analysis_simple_metrics_DihedralDistanceMetric_HH
#define INCLUDED_protocols_analysis_simple_metrics_DihedralDistanceMetric_HH

#include <protocols/analysis/simple_metrics/DihedralDistanceMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <map>

namespace protocols {
namespace analysis {
namespace simple_metrics {

///@brief A metric to calculate the dihedral distance between two poses or the input and the set cmd-line native.  Can set a subset of residues to calculate via ResidueSelector.
///
///  This is the normalized metric in degrees
///
///@details
/// Metric described in:
///   North B, Lehmann A, Dunbrack RL. A new clustering of antibody CDR loop conformations. J Mol Biol 2011; 406:228-256.
///
/// NOTE: Only works for Protein and Carbohydrate residues as dihedral definitions in Rosetta are still aweful.
///        Metric is computed on the protein and carbohydrate backbone.
///
class DihedralDistanceMetric : public core::simple_metrics::RealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	DihedralDistanceMetric();

	DihedralDistanceMetric( core::select::residue_selector::ResidueSelectorCOP selector );

	/// @brief Copy constructor (not needed unless you need deep copies)
	DihedralDistanceMetric( DihedralDistanceMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~DihedralDistanceMetric() override;

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
	//apply( pose::Pose & pose, prefix="", suffix="" ) override;

	///@brief Calculate DihedralDistance metric.
	///  This is the normalized metric in degrees
	///
	///@details
	/// Metric described in:
	///North B, Lehmann A, Dunbrack RL. A new clustering of antibody CDR loop conformations. J Mol Biol 2011; 406:228-256.
	///
	core::Real
	calculate( core::pose::Pose const & pose ) const override;

public:

	///@brief Set a reference pose to calculate rmsd
	///
	void
	set_comparison_pose( core::pose::Pose const & pose );

	///@brief Set a residue selector to calculate total energy of a subset of residues.
	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP residue_selector );

	///@brief Set a residue selector for the comparison pose.
	void
	set_residue_selector_reference( core::select::residue_selector::ResidueSelectorCOP residue_selector_ref);

	///@brief Set a map to compute the Dihedral distance on input->reference residue numbers.
	void
	set_residue_mapping( std::map< core::Size, core::Size> const & res_map );

public:

	///@brief Set the boolean to include protein omega.
	/// Default false.
	///
	void
	set_include_protein_omega( bool include_omega );

	///@brief Load any set in:file:native as the reference pose.
	///  This is opt-in to save on loading time.
	void
	load_native_pose_as_reference();

public:

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

	core::select::residue_selector::ResidueSelectorCOP residue_selector_ = nullptr;
	core::select::residue_selector::ResidueSelectorCOP residue_selector_ref_ = nullptr;

	core::pose::PoseOP ref_pose_ = nullptr;
	bool include_protein_omega_ = false;

	std::map< core::Size, core::Size > res_map_;

};

} //core
} //simple_metrics
} //metrics



#endif //protocols_analysis_simple_metrics_DihedralDistanceMetric_HH





