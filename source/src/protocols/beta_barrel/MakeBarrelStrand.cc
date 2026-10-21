// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/MakeBarrelStrand.cc
/// @brief  Builds a single strand as part of a beta-barrel.
/// @details The barrel is centred on the origin, with the barrel axis pointing along the
/// global z-axis.
/// @author Andy Watkins

// Unit Headers
#include <protocols/beta_barrel/MakeBarrelStrand.hh>
#include <protocols/beta_barrel/MakeBarrelStrandCreator.hh>
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <protocols/beta_barrel/BarrelParametrizationCalculator.hh>
#include <core/conformation/parametric/SizeValuedParameter.hh>

#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/pointer/memory.hh>

//Auto Headers
#include <core/pose/Pose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <protocols/moves/mover_schemas.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/conformation/Conformation.hh>




namespace protocols {
namespace beta_barrel {

static basic::Tracer TR("protocols.beta_barrel.MakeBarrelStrand");

/// @brief Constructor for MakeBarrelStrand mover.
MakeBarrelStrand::MakeBarrelStrand():
	Mover("MakeBarrelStrand"),
	reset_pose_(true),
	strand_length_(7),
	residue_name_(),
	last_apply_failed_(false),
	calculator_( utility::pointer::make_shared< BarrelParametrizationCalculator >() )
{
	set_minor_helix_params_from_file("beta_strand"); //By default, set the minor helix parameters to those of a beta strand (read in from the database).
	residue_name_.push_back("ALA");
}


/// @brief Copy constructor for MakeBarrelStrand mover.
MakeBarrelStrand::MakeBarrelStrand( MakeBarrelStrand const &src ):
	protocols::moves::Mover( src ),
	reset_pose_(src.reset_pose_),
	strand_length_(src.strand_length_),
	residue_name_(src.residue_name_),
	last_apply_failed_(src.last_apply_failed_),
	calculator_( utility::pointer::static_pointer_cast< BarrelParametrizationCalculator >( src.calculator_->clone() ) )
{
}

/// @brief Initialization constructor: initializes this MakeBarrelStrand mover with a BarrelParametrizationCalculator.
/// @details Input calculator is cloned.
MakeBarrelStrand::MakeBarrelStrand( BarrelParametrizationCalculatorCOP input_calculator ) :
	Mover("MakeBarrelStrand"),
	reset_pose_(true),
	strand_length_(7),
	residue_name_(),
	last_apply_failed_(false),
	calculator_( utility::pointer::static_pointer_cast< BarrelParametrizationCalculator >( input_calculator->clone() ) )
{
}

/// @brief Destructor for MakeBarrelStrand mover.
MakeBarrelStrand::~MakeBarrelStrand() = default;


/// @brief Clone operator to create a pointer to a fresh MakeBarrelStrand object that copies this one.
protocols::moves::MoverOP MakeBarrelStrand::clone() const {
	return utility::pointer::make_shared< MakeBarrelStrand > ( *this );
}


/// @brief Fresh_instance operator to create a pointer to a fresh MakeBarrelStrand object that does NOT copy this one.
protocols::moves::MoverOP MakeBarrelStrand::fresh_instance() const {
	return utility::pointer::make_shared< MakeBarrelStrand >();
}

/// @brief Copy the parameter values for parameters that have not been set from the global parameters.
/// @details This function should be called before apply().
void
MakeBarrelStrand::copy_unset_params_from_globals(
	BarrelParametrizationCalculatorCOP global_calculator
) {
	calculator_->copy_unset_params_from_globals( global_calculator );
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
void
MakeBarrelStrand::apply (
	core::pose::Pose & pose
) {
	if ( TR.visible() ) TR << "Building a strand in a beta-barrel using the Crick equations." << std::endl;

	//Initial checks:
	runtime_assert_string_msg(
		calculator_->residues_per_repeat() == residue_name_.size(),
		"In protocols::beta_barrel::MakeBarrelStrand::apply(): The number of residues per repeat does not match the size of the list of residue types."
	);

	//Create the pose object:
	core::pose::Pose strandpose;

	//Build the pose:
	if ( TR.Debug.visible() ) TR.Debug << "Doing initial build." << std::endl;
	protocols::cyclic_peptide::PeptideStubMover stubmover;
	stubmover.set_reset_mode( true );
	stubmover.reset_mover_data();
	core::Size repeat_index( 0 ); //Index in the repeating unit making up the minor helix
	core::Size const residues_per_repeat( calculator_->residues_per_repeat() );
	for ( core::Size i=1, imax=strand_length(); i<=imax; ++i ) {
		++repeat_index;
		if ( repeat_index > residues_per_repeat ) repeat_index=1;
		stubmover.add_residue ("Append", residue_name(repeat_index), 0, (i==1 ? true : false), "", 1, 0, nullptr, "");
	}
	stubmover.apply(strandpose);

	set_last_apply_failed( calculator_->build_strand( strandpose ) );

	if ( last_apply_failed() ) {
		if ( TR.visible() ) TR << "Mover failed.  The Crick parameters do not generate a sensible strand.  Returning input pose." << std::endl;
		return; //At this point, the input pose has not been modified.
	}

	//Either reset the pose and replace it with the strand pose, or append the strand pose to the current pose, depending on the reset mode:
	if ( reset_pose() ) {
		if ( TR.Debug.visible() ) TR.Debug << "Clearing pose and adding strand to pose." << std::endl;
		pose.clear();
		pose=strandpose;
	} else {
		if ( TR.Debug.visible() ) TR.Debug << "Appending strand to pose." << std::endl;
		pose.append_pose_by_jump(strandpose, 1);
	}

	if ( TR.Debug.visible() ) TR.Debug << "Adding Crick parameter data to Conformation object." << std::endl;
	add_parameter_info_to_pose( pose );

	if ( TR.Debug.visible() ) TR.Debug << "Finished apply function." << std::endl;

	return;
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////

/// @brief Add Crick parameter information to the Conformation object within the pose.
/// @details This function updates the calculator's Parameters object's links to residues within the pose,
/// and then clones the Parameters object into the pose's Conformation object.  The new Parameters
/// object will be added to a new ParameterSet object in the Conformation object.
void MakeBarrelStrand::add_parameter_info_to_pose( core::pose::Pose &pose )
{
	core::Size const first_res( pose.size() - strand_length() + 1 );
	core::Size const last_res( pose.size() );

#ifdef NDEBUG
	parameters::BarrelParametersOP output_parameters( utility::pointer::static_pointer_cast< parameters::BarrelParameters >( calculator_->parameters_cop()->clone() ) );
#else
	parameters::BarrelParametersOP output_parameters( utility::pointer::dynamic_pointer_cast< parameters::BarrelParameters >( calculator_->parameters_cop()->clone() ) );
	debug_assert( output_parameters != nullptr );
#endif

	output_parameters->reset_residue_list();

	for ( core::Size ir=first_res; ir<=last_res; ++ir ) { //Loop through all of the strand residues.
		output_parameters->add_residue( pose.conformation().residue_cop( ir ) ); //Add owning pointers to the residue objects.
	}

	//Create a new ParametersSet object:
	parameters::BarrelParametersSetOP output_parameters_set( utility::pointer::make_shared< BarrelParametersSet >() );
	output_parameters_set->add_parameters( output_parameters );

	//Create a new ParametersSet in the Conformation object:
	pose.conformation().add_parameters_set( output_parameters_set );
}

/// @brief Set the minor helix parameters by reading them in from a file.
///
void
MakeBarrelStrand::set_minor_helix_params_from_file (
	std::string const &filename
) {
	calculator_->init_from_file(filename);
}

std::string MakeBarrelStrand::get_name() const {
	return mover_name();
}

std::string MakeBarrelStrand::mover_name() {
	return "MakeBarrelStrand";
}

void MakeBarrelStrand::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "This is a helper mover called by the MakeBarrel mover.  It is NOT intended to be invoked directly from RosettaScripts.  As such, it has no configurable settings.", attlist );
}

std::string MakeBarrelStrandCreator::keyname() const {
	return MakeBarrelStrand::mover_name();
}

protocols::moves::MoverOP
MakeBarrelStrandCreator::create_mover() const {
	return utility::pointer::make_shared< MakeBarrelStrand >();
}

void MakeBarrelStrandCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MakeBarrelStrand::provide_xml_schema( xsd );
}


} //namespace beta_barrel
} //namespace protocols
