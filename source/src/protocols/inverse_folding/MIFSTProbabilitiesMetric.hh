// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/inverse_folding/MIFSTProbabilitiesMetric.hh
/// @brief A PerResidueProbabilitiesMetric that stores amino acid probabilities predicted by the MIF-ST model.
/// @author Moritz Ertelt (moritz.ertelt@gmail.com)

#ifndef INCLUDED_protocols_inverse_folding_MIFSTProbabilitiesMetric_HH
#define INCLUDED_protocols_inverse_folding_MIFSTProbabilitiesMetric_HH

#include <protocols/inverse_folding/MIFSTProbabilitiesMetric.fwd.hh>
#include <core/simple_metrics/PerResidueProbabilitiesMetric.hh>

// protocol headers
#include <protocols/inverse_folding/MIFST.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace inverse_folding {

///@brief A PerResidueProbabilitiesMetric that stores amino acid probabilities predicted by the MIF-ST model.
class MIFSTProbabilitiesMetric : public core::simple_metrics::PerResidueProbabilitiesMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	MIFSTProbabilitiesMetric();

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~MIFSTProbabilitiesMetric() override;


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
	//apply( pose::Pose & pose, prefix="", suffix="" ) override;


	/// @brief Set the residue selector that we'll be using.
	/// @details Passing nullptr results in no residue selector being used.
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );

	/// @brief Get the residue selector.
	/// @details If this returns nullptr, it means that no residue selector is being used.
	core::select::residue_selector::ResidueSelectorCOP residue_selector() const;

	/// @brief A second optional residue selector for feature selection.
	void set_feature_selector(core::select::residue_selector::ResidueSelectorCOP selector);

	///@brief set the multirun option
	void set_multirun( bool multirun );

	///@brief set the use_gpu option
	void set_use_gpu( bool use_gpu );

	///@brief Calculate the metric.
	/// This map is Rosetta Resnum->(AA, value) and includes only those residues selected.
	///
	///@details
	/// Return by value as this function can not STORE the result, it only calculates.
	/// Store the result in the pose by using the apply method, which calls this method and stores the result
	/// in the pose as ExtraScoreValues.
	std::map< core::Size, std::map< core::chemical::AA, core::Real > >
	calculate( core::pose::Pose const & pose ) const override;

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

	/// @brief This simple metric is unpublished.  It returns Moritz Ertelt as its author.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

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


private: //Data

	/// @brief An optional residue selector.
	core::select::residue_selector::ResidueSelectorCOP residue_selector_ = nullptr;

	/// @brief Residue selector for attention masking.
	core::select::residue_selector::ResidueSelectorCOP selector_two_ = nullptr;

	///@brief whether to run inference on all residues at once, or one-by-one
	bool multirun_ = true;

	///@brief whether to run inference on the GPU (if one is available)
	bool use_gpu_ = false;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //inverse_folding
} //protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_inverse_folding_MIFSTProbabilitiesMetric )
#endif // SERIALIZATION

#endif //protocols_inverse_folding_MIFSTProbabilitiesMetric_HH
