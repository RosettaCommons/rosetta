// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/PerturbBarrelStrand.cc
/// @brief  Reads the barrel parameters for a single strand in a pose and rebuilds the strand
/// geometry accordingly.  The parameters are presumed to have been perturbed by another mover.
/// This mover is intended to be called by the PerturbBarrel mover, which handles the
/// perturbation of the barrel parameters.
/// @author Andy Watkins

// Unit Headers
#include <protocols/beta_barrel/PerturbBarrelStrand.hh>
#include <protocols/beta_barrel/PerturbBarrelStrandCreator.hh>
#include <protocols/beta_barrel/BarrelParametrizationCalculator.hh>

#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>

//Auto Headers
#include <core/pose/Pose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <protocols/moves/mover_schemas.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/conformation/Conformation.hh>


namespace protocols {
namespace beta_barrel {

static basic::Tracer TR("protocols.beta_barrel.PerturbBarrelStrand");


/// @brief Constructor for PerturbBarrelStrand mover.
PerturbBarrelStrand::PerturbBarrelStrand():
	Mover("PerturbBarrelStrand"),
	parameters_set_index_(0),
	parameters_index_(0),
	last_apply_failed_(false)
{
}


/// @brief Destructor for PerturbBarrelStrand mover.
PerturbBarrelStrand::~PerturbBarrelStrand() = default;


/// @brief Clone operator to create a pointer to a fresh PerturbBarrelStrand object that copies this one.
protocols::moves::MoverOP PerturbBarrelStrand::clone() const {
	return utility::pointer::make_shared< PerturbBarrelStrand >( *this );
}


/// @brief Fresh_instance operator to create a pointer to a fresh PerturbBarrelStrand object that does NOT copy this one.
protocols::moves::MoverOP PerturbBarrelStrand::fresh_instance() const {
	return utility::pointer::make_shared< PerturbBarrelStrand >();
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
void PerturbBarrelStrand::apply (core::pose::Pose & pose)
{
	if ( TR.visible() ) TR << "Perturbing a strand in a beta-barrel using the barrel parametrization equations." << std::endl;

	runtime_assert_string_msg( parameters_set_index()>0  && parameters_set_index()<=pose.conformation().n_parameters_sets(),
		"In protocols::beta_barrel::PerturbBarrelStrand::apply() : The index of the ParametersSet object is not an index that exists in the pose." );
	runtime_assert_string_msg( parameters_index()>0  && parameters_index()<=pose.conformation().parameters_set( parameters_set_index() )->n_parameters(),
		"In protocols::beta_barrel::PerturbBarrelStrand::apply() : The index of the Parameters object is not an index that exists in the ParametersSet object." );

	bool failed(false);

	BarrelParametersSetOP paramset = utility::pointer::dynamic_pointer_cast<BarrelParametersSet>( pose.conformation().parameters_set( parameters_set_index() ) );
	runtime_assert_string_msg( paramset, "In protocols::beta_barrel::PerturbBarrelStrand::apply() : The ParametersSet object is not a BarrelParametersSet." );
	BarrelParametersOP params = utility::pointer::dynamic_pointer_cast<BarrelParameters>( paramset->parameters( parameters_index() ) );
	runtime_assert_string_msg( params, "In protocols::beta_barrel::PerturbBarrelStrand::apply() : The Parameters object is not a BarrelParameters object." );

	BarrelParametrizationCalculator temp_calculator( false, params );

	failed = temp_calculator.build_strand( pose, params->first_residue()->seqpos(), params->last_residue()->seqpos() );

	set_last_apply_failed(failed);
	if ( failed ) {
		if ( TR.visible() ) TR << "Mover failed.  The barrel parameters do not generate sensible strand geometry.  Returning input pose." << std::endl;
		return; //At this point, the input pose has not been modified.
	}

	if ( TR.Debug.visible() ) TR.Debug << "Finished apply function." << std::endl;

	return;
}

////////////////////////////////////////////////////////////////////////////////


/// @brief Returns the name of this mover ("PerturbBarrelStrand").
std::string PerturbBarrelStrand::get_name() const {
	return mover_name();
}

std::string PerturbBarrelStrand::mover_name() {
	return "PerturbBarrelStrand";
}

void
PerturbBarrelStrand::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "This is a helper mover called by the PerturbBarrel mover.  It is NOT intended to be invoked directly from RosettaScripts.  As such, it has no configurable settings.", attlist );
}

std::string PerturbBarrelStrandCreator::keyname() const {
	return PerturbBarrelStrand::mover_name();
}

protocols::moves::MoverOP
PerturbBarrelStrandCreator::create_mover() const {
	return utility::pointer::make_shared< PerturbBarrelStrand >();
}

void PerturbBarrelStrandCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PerturbBarrelStrand::provide_xml_schema( xsd );
}


} //namespace beta_barrel
} //namespace protocols
