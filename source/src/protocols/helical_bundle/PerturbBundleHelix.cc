// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/PerturbBundleHelix.cc
/// @brief  Reads the Crick parameters for a piece of a pose from the input pose and sets the
/// mainchain torsions accordingly.  The parameters are presumed to have been perturbed by
/// another mover.  This mover is intended to be called by the PerturbBundle mover, which handles
/// the perturbation of the Crick parameters.
/// @details The bundle is centred on the origin, with the outer helix axis pointing along the
/// global z-axis.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit Headers
#include <protocols/helical_bundle/PerturbBundleHelix.hh>
#include <protocols/helical_bundle/PerturbBundleHelixCreator.hh>
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <protocols/helical_bundle/BundleParametrizationCalculator.hh>
#include <numeric/crick_equations/BundleParams.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/NamedAtomID.hh>
#include <core/scoring/rms_util.hh>
#include <core/pose/util.tmpl.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using basic::Error;
using basic::Warning;

//static numeric::random::RandomGenerator RG(741701);  // <- Magic number, do not change it!

namespace protocols {
namespace helical_bundle {

static basic::Tracer TR("protocols.helical_bundle.PerturbBundleHelix");

// XRW TEMP std::string
// XRW TEMP PerturbBundleHelixCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return PerturbBundleHelix::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP PerturbBundleHelixCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new PerturbBundleHelix );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP PerturbBundleHelix::mover_name()
// XRW TEMP {
// XRW TEMP  return "PerturbBundleHelix";
// XRW TEMP }


/// @brief Creator for PerturbBundleHelix mover.
PerturbBundleHelix::PerturbBundleHelix():
	Mover("PerturbBundleHelix"),
	parameters_set_index_(0),
	parameters_index_(0),
	last_apply_failed_(false)
{
}


/// @brief Destructor for PerturbBundleHelix mover.
PerturbBundleHelix::~PerturbBundleHelix() = default;


/// @brief Clone operator to create a pointer to a fresh PerturbBundleHelix object that copies this one.
protocols::moves::MoverOP PerturbBundleHelix::clone() const {
	return protocols::moves::MoverOP( new PerturbBundleHelix( *this ) );
}


/// @brief Fresh_instance operator to create a pointer to a fresh PerturbBundleHelix object that does NOT copy this one.
protocols::moves::MoverOP PerturbBundleHelix::fresh_instance() const {
	return protocols::moves::MoverOP( new PerturbBundleHelix );
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
void PerturbBundleHelix::apply (core::pose::Pose & pose)
{
	if ( TR.visible() ) TR << "Perturbing a helix in a helical bundle using the Crick equations." << std::endl;

	runtime_assert_string_msg( parameters_set_index()>0  && parameters_set_index()<=pose.conformation().n_parameters_sets(),
		"In protocols::helical_bundle::PerturbBundleHelix::apply() : The index of the ParametersSet object is not an index that exists in the pose." );
	runtime_assert_string_msg( parameters_index()>0  && parameters_index()<=pose.conformation().parameters_set( parameters_set_index() )->n_parameters(),
		"In protocols::helical_bundle::PerturbBundleHelix::apply() : The index of the Parameters object is not an index that exists in the ParametersSet object." );

	bool failed(false);

	BundleParametersSetOP paramset = utility::pointer::dynamic_pointer_cast<BundleParametersSet>( pose.conformation().parameters_set( parameters_set_index() ) );
	runtime_assert_string_msg( paramset, "In protocols::helical_bundle::PerturbBundleHelix::apply() : The ParametersSet object is not a BundleParametersSet." );
	BundleParametersOP params = utility::pointer::dynamic_pointer_cast<BundleParameters>( paramset->parameters( parameters_index() ) );
	runtime_assert_string_msg( params, "In protocols::helical_bundle::PerturbBundleHelix::apply() : The Parameters object is not a BundleParameters object." );

	BundleParametrizationCalculator temp_calculator( false, params );

	failed = temp_calculator.build_helix( pose, params->first_residue()->seqpos(), params->last_residue()->seqpos() );

	set_last_apply_failed(failed);
	if ( failed ) {
		if ( TR.visible() ) TR << "Mover failed.  The Crick parameters do not generate a sensible helix.  Returning input pose." << std::endl;
		return; //At this point, the input pose has not been modified.
	}

	if ( TR.Debug.visible() ) TR.Debug << "Finished apply function." << std::endl;

	return;
}

////////////////////////////////////////////////////////////////////////////////


/// @brief Returns the name of this mover ("PerturbBundleHelix").
std::string PerturbBundleHelix::get_name() const {
	return mover_name();
}

std::string PerturbBundleHelix::mover_name() {
	return "PerturbBundleHelix";
}

void
PerturbBundleHelix::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "This is a helper mover called by the PerturbBundle mover.  It is NOT intended to be invoked directly from RosettaScripts.  As such, it has no configurable settings.", attlist );
}

std::string PerturbBundleHelixCreator::keyname() const {
	return PerturbBundleHelix::mover_name();
}

protocols::moves::MoverOP
PerturbBundleHelixCreator::create_mover() const {
	return protocols::moves::MoverOP( new PerturbBundleHelix );
}

void PerturbBundleHelixCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PerturbBundleHelix::provide_xml_schema( xsd );
}


////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
//          (DELETED)
////////////////////////////////////////////////////////////////////////////////
// @brief parse XML (specifically in the context of the parser/Rosetta_scripting scheme)
//
/*void
PerturbBundleHelix::parse_my_tag(
utility::tag::TagCOP tag,
basic::datacache::DataMap & data_map,
protocols::filters::Filters_map const &filters,
protocols::moves::Movers_map const &movers,
core::pose::Pose const & pose
) {


return;
}*/


////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////


} //namespace helical_bundle
} //namespace protocols
