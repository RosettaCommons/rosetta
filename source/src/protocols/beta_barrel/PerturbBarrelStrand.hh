// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/PerturbBarrelStrand.hh
/// @brief  Headers for PerturbBarrelStrand.cc.  Reads the barrel parameters for a single
/// strand in a pose and rebuilds the strand geometry accordingly.  The parameters are
/// presumed to have been perturbed by another mover.  This mover is intended to be called
/// by the PerturbBarrel mover, which handles the perturbation of the barrel parameters.
/// @author Andy Watkins

#ifndef INCLUDED_protocols_beta_barrel_PerturbBarrelStrand_hh
#define INCLUDED_protocols_beta_barrel_PerturbBarrelStrand_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/beta_barrel/PerturbBarrelStrand.fwd.hh>
#include <protocols/beta_barrel/parameters/BarrelParameters.fwd.hh>
#include <protocols/beta_barrel/parameters/BarrelParameters.hh>
#include <protocols/beta_barrel/parameters/BarrelParametersSet.fwd.hh>
#include <protocols/beta_barrel/parameters/BarrelParametersSet.hh>
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.fwd.hh>
#include <core/conformation/parametric/ParametersSet.hh>

// Scripter Headers
#include <protocols/moves/Mover.fwd.hh>

// Project Headers


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace beta_barrel {

class PerturbBarrelStrand : public protocols::moves::Mover
{
public:
	//Typedefs:
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
	PerturbBarrelStrand();
	~PerturbBarrelStrand() override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;


	/// @brief Actually apply the mover to the pose.
	void apply(core::pose::Pose & pose) override;

public: //Setters:

	/// @brief Sets the index of the ParametersSet object in the pose that's being perturbed.
	/// @details This ParametersSet contains the Parameters object that defines the strand that's being perturbed.
	void set_parameters_set_index( core::Size const index ) { parameters_set_index_=index; return; }

	/// @brief Sets the index of the Parameters object in the ParametersSet object in the pose that's being perturbed.
	/// @details This Parameters object contains the barrel parameters that define the strand that's being perturbed.
	void set_parameters_index( core::Size const index ) { parameters_index_=index; return; }


public: //Getters:

	/// @brief Returns the index of the ParametersSet object in the pose that's being perturbed.
	/// @details This ParametersSet contains the Parameters object that defines the strand that's being perturbed.
	core::Size parameters_set_index() const { return parameters_set_index_; }

	/// @brief Returns the index of the Parameters object in the ParametersSet object in the pose that's being perturbed.
	/// @details This Parameters object contains the barrel parameters that define the strand that's being perturbed.
	core::Size parameters_index() const { return parameters_index_; }

	bool last_apply_failed() const { return last_apply_failed_; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE DATA                                                      //
	////////////////////////////////////////////////////////////////////////////////

private:

	/// @brief The index of the ParametersSet object in the pose that's being perturbed.
	/// @details This ParametersSet contains the Parameters object that defines the strand that's being perturbed.
	core::Size parameters_set_index_;

	/// @brief The index of the Parameters object in the ParametersSet object in the pose that's being perturbed.
	/// @details This Parameters object contains the barrel parameters that define the strand that's being perturbed.
	core::Size parameters_index_;

	/// @brief Stores whether the last apply attempt failed.
	///
	bool last_apply_failed_;

	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE FUNCTIONS                                                 //
	////////////////////////////////////////////////////////////////////////////////

private:

	void set_last_apply_failed(bool const status) { last_apply_failed_=status; return;}

};

} //namespace beta_barrel
} //namespace protocols

#endif //INCLUDED_protocols_beta_barrel_PerturbBarrelStrand_hh
