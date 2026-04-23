// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/RMSDMetric.hh
/// @brief A metric to calculate the RMSD between two poses or the input.  Can set a subset of residues to calculate via ResidueSelector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_metrics_RMSDMetric_HH
#define INCLUDED_core_simple_metrics_metrics_RMSDMetric_HH

#include <core/simple_metrics/metrics/RMSDMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>
#include <core/scoring/rms_enum.hh>

// Core headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <map>

namespace core {
namespace simple_metrics {
namespace metrics {

///@brief A metric to calculate the RMSD between two poses.
/// Please note that by default it doesn't align poses. You need to set set_run_superimpose() to True for performing alignment.
/// You Can set a subset of residues to calculate via ResidueSelector.
///
/// By default RMSD is calcualted between all atoms other than hydrogens (all_heavy). This can be changed via **set_rmsd_type**
///
/// @details We match all corresponding atoms and do not fail if a residue does not match up.
///
class RMSDMetric : public core::simple_metrics::RealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	RMSDMetric();

	RMSDMetric( core::pose::PoseCOP ref_pose);

	RMSDMetric( core::pose::PoseCOP ref_pose, core::select::residue_selector::ResidueSelectorCOP selector );

	/// @brief Copy constructor (not needed unless you need deep copies)
	RMSDMetric( RMSDMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~RMSDMetric() override;

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

	///@brief Calculate the metric. This is the RMSD between the input and the set comparison pose.
	/// Deafult is to calculate all_heavy atoms - but this can be set.
	///
	///@details Make sure that reference pose and set pose are the same length or set a map to compare specific residues.
	///           We match all corresponding atoms for each residue to match.
	///
	core::Real
	calculate( core::pose::Pose const & pose ) const override;

public: //Functions needed for the citation manager

	/// @brief Provide the citation.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

public:

	///@brief Set a reference pose to calculate rmsd
	///
	void
	set_comparison_pose( core::pose::PoseCOP pose);


	///@brief Set a residue selector to calculate total energy of a subset of residues.
	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP residue_selector );

	///@breif Set a reference residue selector.  Both selectors should return the same number of residues.
	/// If not set, will use the main residue selector.
	void
	set_residue_selector_reference( core::select::residue_selector::ResidueSelectorCOP residue_selector );


	///@brief Optional - Set a residue selector specifically for superposition.
	void
	set_residue_selector_super( core::select::residue_selector::ResidueSelectorCOP residue_selector );

	///@breif Set a reference residue selector for superposition.  Both selectors should return the same number of residues.
	/// If not set, but super selector is set, will use that.
	void
	set_residue_selector_super_reference( core::select::residue_selector::ResidueSelectorCOP residue_selector );


	///@brief Set a map to compute the RMSD on input->reference residue numbers.
	void
	set_residue_mapping( std::map< core::Size, core::Size> const & rmsd_map );

	///@brief Set maps of cyclic pose to compute the RMSDs on input->reference residue numbers.
	void
	set_cyclic_residue_mappings( std::vector< std::map< core::Size, core::Size> > const & rmsd_maps );

	///@brief Run a superimpose on the residues selected in the residue selector (or all)
	/// default False.
	void
	set_run_superimpose( bool super );

	///@brief Set what we will be calculating the RMSD on.
	///
	/// By default RMSD is calcualted between all atoms other than hydrogens (all_heavy).
	///
	/// You can look up core/rms/rms_enum.hh for available options. Right now the available options are:
	/// rmsd_protein_bb_heavy: backbone heavy atoms N,CA,C
	/// rmsd_protein_bb_heavy_including_O: backbone heavy atoms N,CA,C,O
	/// rmsd_protein_bb_ca: CA (C-alpha)
	/// rmsd_sc_heavy: side chain heavy atoms
	/// rmsd_sc: all side chain atoms (including H)
	/// rmsd_all_heavy: all heavy atoms
	/// rmsd_all: all atoms
	///
	/// Example:
	/// metric = rosetta.core.simple_metrics.metrics.RMSDMetric(pose_1)
	/// metric.set_rmsd_type(rosetta.core.scoring.rmsd_protein_bb_heavy)
	/// metric.set_run_superimpose(True)
	/// rmsd = metric.calculate(pose_2)
	///
	void
	set_rmsd_type( scoring::rmsd_atoms rmsd_type );

	///@brief Set whether we are robust to atom mismatches for selected residues.
	//  By default we only match atoms that are corresponding. (True).
	//
	/// Set this to false to fail instead.
	///
	void
	set_corresponding_atoms_robust( bool robust );

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

	///@brief Setup Str to rmsd_atom map
	void
	setup_name_mapping();

private:

	core::select::residue_selector::ResidueSelectorCOP residue_selector_ = nullptr;
	core::select::residue_selector::ResidueSelectorCOP residue_selector_ref_ = nullptr;

	core::select::residue_selector::ResidueSelectorCOP residue_selector_super_ = nullptr;
	core::select::residue_selector::ResidueSelectorCOP residue_selector_super_ref_ = nullptr;

	std::map< core::Size, core::Size > rmsd_map_;
	std::vector< std::map< core::Size, core::Size > > rmsd_maps_;

	core::pose::PoseCOP ref_pose_ = nullptr;
	scoring::rmsd_atoms rmsd_type_ = scoring::rmsd_all_heavy;
	utility::vector1< std::string> override_atom_names_; //List of atom names to match '  CA ' for example.
	bool robust_ = true;

	std::map< std::string, scoring::rmsd_atoms > name_mapping_;
	bool superimpose_=false;
	bool cyclic_pose_ = false;

};


} // metrics
} // simple_metrics
} // core



#endif //protocols_analysis_simple_metrics_RMSDMetric_HH





