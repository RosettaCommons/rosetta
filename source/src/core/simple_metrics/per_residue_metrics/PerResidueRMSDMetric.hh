// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/PerResidueRMSDMetric.hh
/// @brief A per-residue metric thtat will calculate the RMSD for each residue given in a residue selector to a reference pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_per_residue_metrics_PerResidueRMSDMetric_HH
#define INCLUDED_core_simple_metrics_per_residue_metrics_PerResidueRMSDMetric_HH

#include <core/simple_metrics/per_residue_metrics/PerResidueRMSDMetric.fwd.hh>
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/scoring/rms_enum.hh>

// Core headers
#include <core/types.hh>
#include <core/id/AtomID.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <map>

namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

///@brief A per-residue metric thtat will calculate the RMSD for each residue given in a residue selector to a reference pose.
class PerResidueRMSDMetric : public core::simple_metrics::PerResidueRealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	PerResidueRMSDMetric();

	PerResidueRMSDMetric( core::pose::PoseCOP ref_pose );

	/// @brief Copy constructor (not needed unless you need deep copies)
	PerResidueRMSDMetric( PerResidueRMSDMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~PerResidueRMSDMetric() override;

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


	///Defined in PerResidueRealMetric
	//void
	//set_residue_selector( select::residue_selector::ResidueSelectorCOP selector );

	///@brief Calculate the metric.
	/// This map is Rosetta Resnum->value and includes only those residues selected.
	///
	///@details
	/// Return by value as this function can not STORE the result, it only calculates.
	/// Store the result in the pose by using the apply method, which calls this method and stores the result
	/// in the pose as ExtraScoreValues.
	std::map< core::Size, core::Real >
	calculate( core::pose::Pose const & pose ) const override;

public:

	///@brief Set a reference pose to calculate rmsd
	///
	void
	set_comparison_pose( core::pose::PoseCOP ref_pose );

	///@breif Set a reference residue selector.  Both selectors should return the same number of residues.
	void
	set_residue_selector_reference( core::select::residue_selector::ResidueSelectorCOP residue_selector );

	///@brief Set a map to compute the RMSD on input->reference residue numbers.
	void
	set_residue_mapping( std::map< core::Size, core::Size> const & rmsd_map );



	///@brief Set what we will be calculating the RMSD on.
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

	///@brief Create an AtomID map according to options set in this class.
	std::map< core::id::AtomID, core::id::AtomID >
	create_atom_id_map(core::pose::Pose const & pose) const;

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

	core::select::residue_selector::ResidueSelectorCOP residue_selector_ref_ = nullptr;
	std::map< core::Size, core::Size > rmsd_map_;

	core::pose::PoseCOP ref_pose_ = nullptr;
	scoring::rmsd_atoms rmsd_type_ = scoring::rmsd_all_heavy;
	utility::vector1< std::string> override_atom_names_; //List of atom names to match '  CA ' for example.
	std::map< scoring::rmsd_atoms, utility::vector1< std::string > > rmsd_atom_names_;
	bool robust_ = true;

	std::map< std::string, scoring::rmsd_atoms > name_mapping_;

};

} //core
} //simple_metrics
} //per_residue_metrics



#endif //core_simple_metrics_per_residue_metrics_PerResidueRMSDMetric_HH





