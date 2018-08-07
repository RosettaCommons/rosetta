// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/PerturbBundle.hh
/// @brief  Headers for PerturbBundle.cc.  Perturbs a helical bundle by altering the Crick parameters.
/// @details The bundle is centred on the origin, with the outer helix axis pointing along the
/// global z-axis.  This mover calls the PerturbBundleHelix mover.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_helical_bundle_PerturbBundle_hh
#define INCLUDED_protocols_helical_bundle_PerturbBundle_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/helical_bundle/PerturbBundle.fwd.hh>
#include <protocols/helical_bundle/PerturbBundleHelix.fwd.hh>
#include <protocols/helical_bundle/PerturbBundleHelix.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.fwd.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.fwd.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.hh>
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.fwd.hh>
#include <core/conformation/parametric/ParametersSet.hh>
#include <protocols/helical_bundle/BundleParametrizationCalculator.fwd.hh>

// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/ContingentFilter.fwd.hh>
#include <protocols/filters/ContingentFilter.hh>

// Project Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <numeric/constants.hh>

#include <set>

#include <core/grid/CartGrid.fwd.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace helical_bundle {

class PerturbBundle : public protocols::moves::Mover
{
public: //Typedefs

	typedef core::conformation::parametric::Parameters Parameters;
	typedef core::conformation::parametric::ParametersOP ParametersOP;
	typedef core::conformation::parametric::ParametersSet ParametersSet;
	typedef core::conformation::parametric::ParametersSetOP ParametersSetOP;

	typedef protocols::helical_bundle::parameters::BundleParameters BundleParameters;
	typedef protocols::helical_bundle::parameters::BundleParametersOP BundleParametersOP;
	typedef protocols::helical_bundle::parameters::BundleParametersCOP BundleParametersCOP;
	typedef protocols::helical_bundle::parameters::BundleParametersSet BundleParametersSet;
	typedef protocols::helical_bundle::parameters::BundleParametersSetOP BundleParametersSetOP;
	typedef protocols::helical_bundle::parameters::BundleParametersSetCOP BundleParametersSetCOP;

public:
	PerturbBundle();
	PerturbBundle( PerturbBundle const &src );
	~PerturbBundle() override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;


	/// @brief Actually apply the mover to the pose.
	void apply(core::pose::Pose & pose) override;


	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const &
	) override;

public:
	////////////////////////////////////////////////////////////////////////////////
	//          PUBLIC FUNCTIONS                                                  //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Set which bundle parameters set will be used, if more than one is defined in the pose's Conformation object.
	/// @details  A value of n indicates that the nth bundle paramets set encountered will be perturbed.
	void set_bundleparametersset_index(core::Size const val) { bundleparametersset_index_=val; return; }

	/// @brief Returns a value indicating which bundle parameters set will be used, if more than one is defined in the pose's Conformation object.
	/// @details  A value of n indicates that the nth bundle paramets set encountered will be perturbed.
	core::Size bundleparametersset_index() const { return bundleparametersset_index_; }

	/// @brief Clear the list of helices.
	/// @details Clears the individual_helix_calculators_ list.
	void reset_helices();

	/// @brief Add options for a new helix.
	/// @details Return value is a smart pointer to the BundleParametrizationCalculator for the new
	/// helix, cloned from the default_calculator_.
	BundleParametrizationCalculatorOP add_helix( core::Size const helix_index );

	/// @brief Access the default calculator (const access).
	inline BundleParametrizationCalculatorCOP default_calculator_cop() const { return default_calculator_; }

	/// @brief Access the default calculator (nonconst access).
	inline BundleParametrizationCalculatorOP default_calculator() { return default_calculator_; }

	/// @brief Access the calculator for a given helix (const access).
	BundleParametrizationCalculatorCOP individual_helix_calculator_cop( core::Size const helix_calculator_index ) const;

	/// @brief Access the calculator for a given helix (nonconst access).
	BundleParametrizationCalculatorOP individual_helix_calculator( core::Size const helix_calculator_index );

	/// @brief Get the number of individual helices that have been configured (which is *not* necessarily the number
	/// of helices in the pose, since we might only be perturbing a subset).
	inline core::Size n_helices() const { return individual_helix_calculators_.size(); }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE DATA                                                      //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief The BundleParametrizationCalculator object that will be used for setting up default parameters.
	BundleParametrizationCalculatorOP default_calculator_;

	/// @brief The BundleParametrizationCalculator objects for individual helices.  Cloned from the default calculator.
	utility::vector1< std::pair< core::Size, BundleParametrizationCalculatorOP > > individual_helix_calculators_;

	/// @brief Which set of bundle parameters (if there exists more than one) should the mover alter?
	/// @details Defaults to 1.  Higher values indicate the nth set encountered in the ParametersSet list.
	core::Size bundleparametersset_index_;


private:
	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE FUNCTIONS                                                 //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Confirms that a helix has not yet been defined.  Returns "true" if the helix
	/// has NOT been defined, false otherwise.
	bool helix_not_defined( core::Size const helix_index ) const;

	/// @brief Get a calculator for a particular helix.
	/// @details Returns nullptr if no calculator for this helix has been defined.
	BundleParametrizationCalculatorCOP get_calculator_for_helix( core::Size const helix_index ) const;

	/// @brief Perturb the helical bundle parameter values in the pose, subject to the options already set.
	/// @details Called by the apply() function.  Returns true for success, false for failure.
	bool perturb_values( BundleParametersSetOP params_set) const;

	/// @brief Rebuild the helical bundle conformation using the bundle parameter values in the pose.
	///
	void rebuild_conformation(
		core::pose::Pose &pose,
		BundleParametersSetOP params_set,
		core::Size const params_set_index,
		bool &failed
	) const;

	/// @brief Write out the perturbed Crick parameters.
	/// @details The "before" parameter determines whether this is a pre-perturbation
	/// report or a post-perturbation report.
	void write_report(
		BundleParametersSetOP params_set,
		bool const before
	) const;


}; //PerturbBundle class

} //namespace helical_bundle
} //namespace protocols

#endif //INCLUDED_protocols_helical_bundle_PerturbBundle_hh
