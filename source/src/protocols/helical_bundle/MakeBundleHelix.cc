// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/MakeBundleHelix.cc
/// @brief  Builds a single helix as part of a helical bundle.
/// @details The bundle is centred on the origin, with the outer helix axis pointing along the
/// global z-axis.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit Headers
#include <protocols/helical_bundle/MakeBundleHelix.hh>
#include <protocols/helical_bundle/MakeBundleHelixCreator.hh>
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <protocols/helical_bundle/BundleParametrizationCalculator.hh>
#include <numeric/crick_equations/BundleParams.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/conformation/parametric/RealVectorValuedParameter.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>
#include <core/conformation/parametric/SizeVectorValuedParameter.hh>
#include <core/conformation/parametric/SizeValuedParameter.hh>

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
#include <utility/pointer/memory.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using basic::Error;
using basic::Warning;



namespace protocols {
namespace helical_bundle {

static basic::Tracer TR("protocols.helical_bundle.MakeBundleHelix");

/// @brief Constructor for MakeBundleHelix mover.
MakeBundleHelix::MakeBundleHelix():
	Mover("MakeBundleHelix"),
	reset_pose_(true),
	helix_length_(10),
	residue_name_(),
	last_apply_failed_(false),
	calculator_( new BundleParametrizationCalculator )
{
	set_minor_helix_params_from_file("alpha_helix"); //By default, set the minor helix parameters to those of an alpha helix (read in from the database).
	residue_name_.push_back("ALA");
}


/// @brief Copy constructor for MakeBundleHelix mover.
MakeBundleHelix::MakeBundleHelix( MakeBundleHelix const &src ):
	protocols::moves::Mover( src ),
	reset_pose_(src.reset_pose_),
	helix_length_(src.helix_length_),
	residue_name_(src.residue_name_),
	last_apply_failed_(src.last_apply_failed_),
	calculator_( utility::pointer::static_pointer_cast< BundleParametrizationCalculator >( src.calculator_->clone() ) )
{
}

/// @brief Initialization constructor: initializes this MakeBundleHelix mover with a BundleParametrizationCalculator.
/// @details Input calculator is cloned.
MakeBundleHelix::MakeBundleHelix( BundleParametrizationCalculatorCOP input_calculator ) :
	Mover("MakeBundleHelix"),
	reset_pose_(true),
	helix_length_(10),
	residue_name_(),
	last_apply_failed_(false),
	calculator_( utility::pointer::static_pointer_cast< BundleParametrizationCalculator >( input_calculator->clone() ) )
{
}

/// @brief Destructor for MakeBundleHelix mover.
MakeBundleHelix::~MakeBundleHelix() {}


/// @brief Clone operator to create a pointer to a fresh MakeBundleHelix object that copies this one.
protocols::moves::MoverOP MakeBundleHelix::clone() const {
	return protocols::moves::MoverOP( new MakeBundleHelix ( *this ) );
}


/// @brief Fresh_instance operator to create a pointer to a fresh MakeBundleHelix object that does NOT copy this one.
protocols::moves::MoverOP MakeBundleHelix::fresh_instance() const {
	return protocols::moves::MoverOP( new MakeBundleHelix );
}

/// @brief Copy the parameter values for parameters that have not been set from the global parameters.
/// @details This function should be called before apply().
void
MakeBundleHelix::copy_unset_params_from_globals(
	BundleParametrizationCalculatorCOP global_calculator
) {
	calculator_->copy_unset_params_from_globals( global_calculator );
}

/// @brief Copy the parameter values for parameters that copy values from previous helices, from the previous helices.
/// @details This function should be called before apply().
/// @returns Returns true for failure, false for success.
bool
MakeBundleHelix::copy_params_from_previous_helices(
	core::pose::Pose const & prev_helices_pose
) {
	return calculator_->copy_params_from_previous_helices_makebundlehelix_style(prev_helices_pose);
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
void
MakeBundleHelix::apply (
	core::pose::Pose & pose
) {
	if ( TR.visible() ) TR << "Building a helix in a helical bundle using the Crick equations." << std::endl;

	//Initial checks:
	runtime_assert_string_msg(
		calculator_->size_parameter_cop( BPC_residues_per_repeat )->value() == residue_name_.size(),
		"In protocols::helical_bundle::MakeBundleHelix::apply(): The number of residues per repeat does not match the size of the list of residue types."
	);

	//

	//Create the pose object:
	core::pose::Pose helixpose;

	//Build the pose:
	if ( TR.Debug.visible() ) TR.Debug << "Doing initial build." << std::endl;
	protocols::cyclic_peptide::PeptideStubMover stubmover;
	stubmover.set_reset_mode( true );
	stubmover.reset_mover_data();
	core::Size repeat_index( calculator_->size_parameter_cop( BPC_repeating_unit_offset )->value() ); //Index in the repeating unit making up the minor helix
	core::Size const residues_per_repeat( calculator_->size_parameter_cop( BPC_residues_per_repeat )->value() );
	for ( core::Size i=1, imax=helix_length(); i<=imax; ++i ) {
		++repeat_index;
		if ( repeat_index > residues_per_repeat ) repeat_index=1;
		stubmover.add_residue ("Append", residue_name(repeat_index), 0, (i==1 ? true : false), "", 1, 0, "");
	}
	stubmover.apply(helixpose);

	set_last_apply_failed( calculator_->build_helix( helixpose ) );

	if ( last_apply_failed() ) {
		if ( TR.visible() ) TR << "Mover failed.  The Crick parameters do not generate a sensible helix.  Returning input pose." << std::endl;
		return; //At this point, the input pose has not been modified.
	}

	//Either reset the pose and replace it with the helix pose, or append the helix pose to the current pose, depending on the reset mode:
	if ( reset_pose() ) {
		if ( TR.Debug.visible() ) TR.Debug << "Clearing pose and adding helix to pose." << std::endl;
		pose.clear();
		pose=helixpose;
	} else {
		if ( TR.Debug.visible() ) TR.Debug << "Appending helix to pose." << std::endl;
		pose.append_pose_by_jump(helixpose, 1);
	}

	if ( TR.Debug.visible() ) TR.Debug << "Adding Crick parameter data to Conformation oject." << std::endl;
	add_parameter_info_to_pose( pose );

	if ( TR.Debug.visible() ) TR.Debug << "Finished apply function." << std::endl;

	return;
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
////////////////////////////////////////////////////////////////////////////////
/// @brief parse XML (specifically in the context of the parser/Rosetta_scripting scheme)
///
/*void
MakeBundleHelix::parse_my_tag(
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

/// @brief Add Crick parameter information to the Conformation object within the pose.
/// @details This function updates the calculator's Parameters object's links to residues within the pose,
/// and then clones the Parameters object into the pose's Conformation object.  The new Parameters
/// object will be added to a new ParameterSet object in the Conformation object.
void MakeBundleHelix::add_parameter_info_to_pose( core::pose::Pose &pose )
{
	core::Size const first_res( pose.size() - helix_length() + 1 );
	core::Size const last_res( pose.size() );

#ifdef NDEBUG
	parameters::BundleParametersOP output_parameters( utility::pointer::static_pointer_cast< parameters::BundleParameters >( calculator_->parameters_cop()->clone() ) );
#else
	parameters::BundleParametersOP output_parameters( utility::pointer::dynamic_pointer_cast< parameters::BundleParameters >( calculator_->parameters_cop()->clone() ) );
	debug_assert( output_parameters != nullptr );
#endif

	output_parameters->reset_residue_list();
	output_parameters->reset_sampling_and_perturbing_info();

	for ( core::Size ir=first_res; ir<=last_res; ++ir ) { //Loop through all of the helix residues.
		output_parameters->add_residue( pose.conformation().residue_cop( ir ) ); //Add owning pointers to the residue objects.
	}

	//Create a new ParametersSet object:
	parameters::BundleParametersSetOP output_parameters_set( utility::pointer::make_shared< BundleParametersSet >() );
	output_parameters_set->add_parameters( output_parameters );

	//Create a new ParametersSet in the Conformation object:
	pose.conformation().add_parameters_set( output_parameters_set );
}

/// @brief Set the minor helix parameters by reading them in from a file.
///
void
MakeBundleHelix::set_minor_helix_params_from_file (
	std::string const &filename
) {
	calculator_->init_from_file(filename);
}

std::string MakeBundleHelix::get_name() const {
	return mover_name();
}

std::string MakeBundleHelix::mover_name() {
	return "MakeBundleHelix";
}

void MakeBundleHelix::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "This is a helper mover called by the MakeBundle and BundleGridSampler movers.  It is NOT intended to be invoked directly from RosettaScripts.  As such, it has no configurable settings.", attlist );
}

std::string MakeBundleHelixCreator::keyname() const {
	return MakeBundleHelix::mover_name();
}

protocols::moves::MoverOP
MakeBundleHelixCreator::create_mover() const {
	return protocols::moves::MoverOP( new MakeBundleHelix );
}

void MakeBundleHelixCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MakeBundleHelix::provide_xml_schema( xsd );
}


} //namespace helical_bundle
} //namespace protocols
