// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/DensityFitResidueSelector.hh
/// @brief  Select residues that have a bad fit to provided electron denisty map.
///
/// - Original Code and Logic - Reference (eLife 2016, Dimaio)
/// @author Ray Wang (wangyr@uw.edu)
/// @author Frank Dimaio (fdimaio@gmail.com)
///
///  - Logic into Residue Selector
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Sebastian Rämisch (raemisch@scripps.edu)

#ifndef INCLUDED_core_select_residue_selector_DensityFitResidueSelector_HH
#define INCLUDED_core_select_residue_selector_DensityFitResidueSelector_HH

// Unit headers
#include <core/select/residue_selector/DensityFitResidueSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/simple_metrics/per_residue_metrics/PerResidueDensityFitMetric.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

/// @brief Selects residues that are considered 'error' region with respect to the density resolution.
///
/// @details Uses weighted sum of density, density-compared-to-neighbors, rama (where applicable) and cart_bonded to compute
///
class DensityFitResidueSelector : public core::select::residue_selector::ResidueSelector {
public:
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSubset ResidueSubset;

public:

	/// @brief Constructor.
	DensityFitResidueSelector();

	DensityFitResidueSelector( simple_metrics::per_residue_metrics::PerResidueDensityFitMetricCOP den_fit_metric );

	/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
	DensityFitResidueSelector(DensityFitResidueSelector const & src);

public:

	/// @brief Destructor.
	~DensityFitResidueSelector() override;

	/// @brief Clone operator.
	/// @details Copy the current object (creating the copy on the heap) and return an owning pointer
	/// to the copy.  All ResidueSelectors must implement this.
	ResidueSelectorOP clone() const override;

	/// @brief "Apply" function.
	/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
	/// indicating whether each residue is selected ("true") or not ("false").
	ResidueSubset apply( core::pose::Pose const & input_pose ) const override;

	/// @brief XML parse.
	/// @details Parse RosettaScripts tags and set up this mover.
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) override;

	/// @brief Get the mover class name.
	std::string
	get_name() const override;

	/// @brief Get the mover class name.
	static std::string
	class_name();

	/// @brief Provide XSD information, enabling mechanical validation of input XML.
	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public:

	///@brief Set the density fit metric that will be used to do the calculation in this class.
	/// Clones into this class.
	void
	set_den_fit_metric( simple_metrics::per_residue_metrics::PerResidueDensityFitMetricCOP den_fit_metric);

	///@brief Edit the density fit metric used to do the calculation in this class.
	simple_metrics::per_residue_metrics::PerResidueDensityFitMetricOP
	get_den_fit_metric();


	///@brief Set the scorecut to use.  From -3 up.  -.5 is default.  Anything higher, is deemed 'good' for zscore.
	///  This is also the correlation for match_mode.  Correlation of 1.0 fits density completely.  .7 is a good cutoff.
	void
	set_score_cut( Real score_cut );

	///@brief Set to invert - IE select BAD density instead of GOOD.
	void
	set_invert( bool invert );

	///@brief Set up options to use any data already stored in the SM cache through its apply method.
	void
	set_cache_options( bool use_sm_cache, std::string const & prefix="", std::string const & suffix = "", bool fail_on_missing_cache=false);


private:

	core::Real score_cut_ = -.5; //Number usually goes from -3 up.
	bool invert_ = false; // Invert to select bad fit.

	core::simple_metrics::per_residue_metrics::PerResidueDensityFitMetricOP den_fit_metric_;

	bool use_cache_ = false;
	std::string prefix_ = "";
	std::string suffix_ = "";
	bool fail_on_missing_cache_ = false;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION


};


} //core
} //select
} //residue_selector

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_select_residue_selector_DensityFitResidueSelector )
#endif // SERIALIZATION

#endif //INCLUDEDcore_select_residue_selector_DensityFitResidueSelector_HH
