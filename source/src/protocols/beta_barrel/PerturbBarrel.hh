// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/PerturbBarrel.hh
/// @brief  Headers for PerturbBarrel.cc.  Perturbs a beta-barrel by altering the barrel parameters.
/// @details This mover calls the PerturbBarrelStrand mover.
/// @author Andy Watkins

#ifndef INCLUDED_protocols_beta_barrel_PerturbBarrel_hh
#define INCLUDED_protocols_beta_barrel_PerturbBarrel_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/beta_barrel/PerturbBarrel.fwd.hh>
#include <protocols/beta_barrel/parameters/BarrelParameters.fwd.hh>
#include <protocols/beta_barrel/parameters/BarrelParameters.hh>
#include <protocols/beta_barrel/parameters/BarrelParametersSet.fwd.hh>
#include <protocols/beta_barrel/parameters/BarrelParametersSet.hh>
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.fwd.hh>
#include <core/conformation/parametric/ParametersSet.hh>
#include <protocols/beta_barrel/BarrelParametrizationCalculator.fwd.hh>

// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// Project Headers
#include <utility/vector1.hh>



///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace beta_barrel {

class PerturbBarrel : public protocols::moves::Mover
{
public: //Typedefs

	typedef core::conformation::parametric::Parameters Parameters;
	typedef core::conformation::parametric::ParametersOP ParametersOP;
	typedef core::conformation::parametric::ParametersSet ParametersSet;
	typedef core::conformation::parametric::ParametersSetOP ParametersSetOP;

	typedef protocols::beta_barrel::parameters::BarrelParameters BarrelParameters;
	typedef protocols::beta_barrel::parameters::BarrelParametersOP BarrelParametersOP;
	typedef protocols::beta_barrel::parameters::BarrelParametersCOP BarrelParametersCOP;
	typedef protocols::beta_barrel::parameters::BarrelParametersSet BarrelParametersSet;
	typedef protocols::beta_barrel::parameters::BarrelParametersSetOP BarrelParametersSetOP;
	typedef protocols::beta_barrel::parameters::BarrelParametersSetCOP BarrelParametersSetCOP;

public:
	PerturbBarrel();
	PerturbBarrel( PerturbBarrel const &src );
	~PerturbBarrel() override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;


	/// @brief Actually apply the mover to the pose.
	void apply(core::pose::Pose & pose) override;


	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

public:
	////////////////////////////////////////////////////////////////////////////////
	//          PUBLIC FUNCTIONS                                                  //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Set which barrel parameters set will be used, if more than one is defined in the pose's Conformation object.
	/// @details  A value of n indicates that the nth barrel parameters set encountered will be perturbed.
	void set_barrelparametersset_index(core::Size const val) { barrelparametersset_index_=val; return; }

	/// @brief Returns a value indicating which barrel parameters set will be used, if more than one is defined in the pose's Conformation object.
	/// @details  A value of n indicates that the nth barrel parameters set encountered will be perturbed.
	core::Size barrelparametersset_index() const { return barrelparametersset_index_; }

	/// @brief Clear the list of strands.
	/// @details Clears the individual_strand_calculators_ list.
	void reset_strands();

	/// @brief Add options for a new strand.
	/// @details Return value is a smart pointer to the BarrelParametrizationCalculator for the new
	/// strand, cloned from the default_calculator_.
	BarrelParametrizationCalculatorOP add_strand( core::Size const strand_index );

	/// @brief Access the default calculator (const access).
	inline BarrelParametrizationCalculatorCOP default_calculator_cop() const { return default_calculator_; }

	/// @brief Access the default calculator (nonconst access).
	inline BarrelParametrizationCalculatorOP default_calculator() { return default_calculator_; }

	/// @brief Access the calculator for a given strand (const access).
	BarrelParametrizationCalculatorCOP individual_strand_calculator_cop( core::Size const strand_calculator_index ) const;

	/// @brief Access the calculator for a given strand (nonconst access).
	BarrelParametrizationCalculatorOP individual_strand_calculator( core::Size const strand_calculator_index );

	/// @brief Get the number of individual strands that have been configured (which is *not* necessarily the number
	/// of strands in the pose, since we might only be perturbing a subset).
	inline core::Size n_strands() const { return individual_strand_calculators_.size(); }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public: //Function overrides needed for the citation manager:

	/// @brief Provide the citation.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

private:
	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE DATA                                                      //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief The BarrelParametrizationCalculator object that will be used for setting up default parameters.
	BarrelParametrizationCalculatorOP default_calculator_;

	/// @brief The BarrelParametrizationCalculator objects for individual strands.  Cloned from the default calculator.
	utility::vector1< std::pair< core::Size, BarrelParametrizationCalculatorOP > > individual_strand_calculators_;

	/// @brief Which set of barrel parameters (if there exists more than one) should the mover alter?
	/// @details Defaults to 1.  Higher values indicate the nth set encountered in the ParametersSet list.
	core::Size barrelparametersset_index_;


private:
	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE FUNCTIONS                                                 //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Confirms that a strand has not yet been defined.  Returns "true" if the strand
	/// has NOT been defined, false otherwise.
	bool strand_not_defined( core::Size const strand_index ) const;

	/// @brief Get a calculator for a particular strand.
	/// @details Returns nullptr if no calculator for this strand has been defined.
	BarrelParametrizationCalculatorCOP get_calculator_for_strand( core::Size const strand_index ) const;

	/// @brief Perturb the beta-barrel parameter values in the pose, subject to the options already set.
	/// @details Called by the apply() function.  Returns true for success, false for failure.
	bool perturb_values( BarrelParametersSetOP params_set) const;

	/// @brief Rebuild the beta-barrel conformation using the barrel parameter values in the pose.
	///
	void rebuild_conformation(
		core::pose::Pose &pose,
		BarrelParametersSetOP params_set,
		core::Size const params_set_index,
		bool &failed
	) const;

	/// @brief Write out the perturbed barrel parameters.
	/// @details The "before" parameter determines whether this is a pre-perturbation
	/// report or a post-perturbation report.
	void write_report(
		BarrelParametersSetOP params_set,
		bool const before
	) const;


}; //PerturbBarrel class

} //namespace beta_barrel
} //namespace protocols

#endif //INCLUDED_protocols_beta_barrel_PerturbBarrel_hh
