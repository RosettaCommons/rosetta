// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/conformation/parametric/BundleParametrizationCalculator.cc
/// @brief  Function definitions for the BundleParametrizationCalculator class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit header
#include <protocols/helical_bundle/BundleParametrizationCalculator.hh>
#include <core/conformation/parametric/ParametrizationCalculator.hh>
#include <protocols/helical_bundle/util.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.hh>
#include <protocols/helical_bundle/parameters/OmegaBundleParameter.hh>

// Core headers
#include <core/conformation/parametric/Parameter.hh>
#include <core/conformation/parametric/BooleanValuedParameter.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>
#include <core/conformation/parametric/RealVectorValuedParameter.hh>
#include <core/conformation/parametric/SizeValuedParameter.hh>
#include <core/conformation/parametric/SizeVectorValuedParameter.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers

// Utility Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

// C++ headers
#include <tuple>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace helical_bundle {

static basic::Tracer TR( "protocols.helical_bundle.BundleParametrizationCalculator" );

/// @brief Constructor.
///
BundleParametrizationCalculator::BundleParametrizationCalculator( bool const use_degrees/*=false*/ ) :
	core::conformation::parametric::ParametrizationCalculator( utility::pointer::make_shared< protocols::helical_bundle::parameters::BundleParameters >() ),
	use_degrees_(use_degrees)
	//TODO -- initialize variables here.
{
	using namespace core::conformation::parametric;
	for ( core::Size i(1); i<static_cast<core::Size>( BPC_end_of_list ); ++i ) {
		if ( static_cast<BPC_Parameters>(i) == BPC_omega0 ) {
			add_custom_parameter(
				parameter_name_from_enum( static_cast< BPC_Parameters>(i) ),
				parameter_type_from_enum( static_cast< BPC_Parameters >(i) ),
				parameter_description_from_enum( static_cast< BPC_Parameters >(i) ),
				short_parameter_description_from_enum( static_cast< BPC_Parameters >(i) ),
				parameter_units_from_enum( static_cast< BPC_Parameters >(i) ),
				parameter_properties_from_enum( static_cast< BPC_Parameters>(i) ),
				utility::pointer::make_shared< parameters::OmegaBundleParameter >()
			);
		} else {
			add_parameter(
				parameter_name_from_enum( static_cast< BPC_Parameters>(i) ),
				parameter_type_from_enum( static_cast< BPC_Parameters >(i) ),
				parameter_description_from_enum( static_cast< BPC_Parameters >(i) ),
				short_parameter_description_from_enum( static_cast< BPC_Parameters >(i) ),
				parameter_units_from_enum( static_cast< BPC_Parameters >(i) ),
				parameter_properties_from_enum( static_cast< BPC_Parameters>(i) )
			);
		}
		//Set suffix for copy option:
		parameter(i)->set_copy_suffix("copies_helix");
	}

	//Set defaults for special parameters:
	set_use_degrees_for_parameters();
	real_parameter( BPC_epsilon )->set_default_value( 1.0 );
	boolean_parameter( BPC_set_dihedrals )->set_default_value( true );
	boolean_parameter( BPC_set_bondangles )->set_default_value( true );
	boolean_parameter( BPC_set_bondlengths )->set_default_value( true );
}

/// @brief Params constructor.
BundleParametrizationCalculator::BundleParametrizationCalculator(
	bool const use_degrees,
	protocols::helical_bundle::parameters::BundleParametersCOP params_in
) :
	core::conformation::parametric::ParametrizationCalculator( params_in->clone() ),
	use_degrees_( use_degrees )
{}

/// @brief Copy constructor.
/// @details Deep-copies the stored parameters (by calling the copy constructor of the base class).
BundleParametrizationCalculator::BundleParametrizationCalculator( BundleParametrizationCalculator const & src ) :
	core::conformation::parametric::ParametrizationCalculator(src),
	use_degrees_(src.use_degrees_)
{
}

/// @brief Destructor.
BundleParametrizationCalculator::~BundleParametrizationCalculator() {}

/// @brief Copy this object ( allocate actual memory for it )
core::conformation::parametric::ParametrizationCalculatorOP
BundleParametrizationCalculator::clone() const {
	return core::conformation::parametric::ParametrizationCalculatorOP( BundleParametrizationCalculatorOP( new BundleParametrizationCalculator(*this) ) );
}

/// @brief Read a Crick params file and initialize this calculator.
/// @details Triggers a read from disk!
void
BundleParametrizationCalculator::init_from_file(
	std::string const &filename
) {
	utility::vector1 <core::Real> r1;
	core::Real omega1(0.0);
	core::Real z1(0.0);
	utility::vector1 <core::Real> delta_omega1;
	utility::vector1 <core::Real> delta_z1;
	core::Size residues_per_repeat(1);
	utility::vector1 <core::Size> atoms_per_residue;
	read_minor_helix_params ( filename, r1, omega1, z1, delta_omega1, delta_z1, residues_per_repeat, atoms_per_residue );
	realvector_parameter( BPC_r1_peratom )->set_value( r1 );
	real_parameter( BPC_omega1 )->set_value( omega1 );
	real_parameter( BPC_z1 )->set_value( z1 );
	realvector_parameter( BPC_delta_omega1_peratom )->set_value( delta_omega1 );
	realvector_parameter( BPC_delta_z1_peratom )->set_value( delta_z1 );
	size_parameter( BPC_residues_per_repeat )->set_value( residues_per_repeat );
	sizevector_parameter( BPC_atoms_per_residue )->set_value( atoms_per_residue );
}

/// @brief Copy the parameter values for parameters that have not been set from the global parameters.
/// @details This function should be called before build_helix().
void
BundleParametrizationCalculator::copy_unset_params_from_globals(
	BundleParametrizationCalculatorCOP global_calculator
) {
	for ( core::Size i(1); i<static_cast<core::Size>(BPC_end_of_list); ++i ) {
		core::conformation::parametric::ParameterOP curparam( parameter(i) );
		if ( !curparam->value_was_set() ) {
			curparam->copy_value_from_parameter( global_calculator->parameter_cop(i), global_calculator->parameters_cop(), this->parameters_cop() );
		}
	}
}

/// @brief Copy the parameter values for parameters that copy values from previous helices, from the previous helices.
/// @details This function should be called before build_helix().
/// @returns Returns true for failure, false for success.
bool
BundleParametrizationCalculator::copy_params_from_previous_helices_makebundlehelix_style(
	core::pose::Pose const & prev_helices_pose
) {
	bool failed(false);
	for ( core::Size i(1); i<static_cast<core::Size>(BPC_end_of_list); ++i ) {
		core::conformation::parametric::ParameterOP curparam( parameter(i) );
		core::Size const param_index( curparam->copy_from_parameters_index() ); //Index of parameters object to copy *from*.
		if ( param_index != 0 ) {
			debug_assert( prev_helices_pose.conformation().n_parameters_sets() >= param_index );
			//The relevant object will be the 1st Parameters object in the Nth ParametersSet, given the current state of the pose.
			core::conformation::parametric::ParametersSetCOP params_set( prev_helices_pose.conformation().parameters_set( param_index ) );
#ifndef NDEBUG
			debug_assert( utility::pointer::dynamic_pointer_cast< protocols::helical_bundle::parameters::BundleParametersSet const >(params_set) != nullptr );
#endif
			debug_assert( params_set->n_parameters() == 1 );
			core::conformation::parametric::ParametersCOP params( params_set->parameters(1) );
#ifndef NDEBUG
			debug_assert( utility::pointer::dynamic_pointer_cast< protocols::helical_bundle::parameters::BundleParameters const >(params) != nullptr );
#endif
			failed = curparam->copy_value_from_parameter( params->parameter_cop(i), params, this->parameters_cop() );
			if ( failed ) return failed;
		}
	}
	return failed;
}

/// @brief Copy the parameter values for parameters that copy values from previous helices, from the previous helices.
/// @details This function should be called before build_helix().
/// @note This function assumes that the input pose is the whole shebang passed to the PerturbBundle mover.  Unlike
/// copy_params_from_previous_helices_makebundlehelix_style(), this functon does *not* modify the
/// BundleParametrizationCalculator.  Instead it copies the parameter value for the given parameter from the
/// indicated previous helix to the indicated parameter object.
/// @returns Returns true for failure, false for success.
/*static*/ bool
BundleParametrizationCalculator::copy_params_from_previous_helices_perturbbundle_style(
	parameters::BundleParametersSetCOP paramset,
	core::Size const copy_from_helix_index,
	BPC_Parameters const parameter_type,
	core::conformation::parametric::RealValuedParameterOP parameter_to_copy_to,
	protocols::helical_bundle::parameters::BundleParametersCOP current_helix_params
) {
	debug_assert( parameter_to_copy_to != nullptr );
	debug_assert( paramset != nullptr);
	core::Size helix_index(0);\
		bool failed(false);
	for ( core::Size iparams(1), iparamsmax(paramset->n_parameters()); iparams<=iparamsmax; ++iparams ) {
#ifdef NDEBUG
		//Release mode: static cast
		protocols::helical_bundle::parameters::BundleParametersCOP params( utility::pointer::static_pointer_cast< protocols::helical_bundle::parameters::BundleParameters const >( paramset->parameters( iparams ) ) );
#else
		//Debug mode: dynamic cast and check
		protocols::helical_bundle::parameters::BundleParametersCOP params( utility::pointer::dynamic_pointer_cast< protocols::helical_bundle::parameters::BundleParameters const >( paramset->parameters( iparams ) ) );
		debug_assert( params!=nullptr );
#endif
		++helix_index;
		if ( helix_index != copy_from_helix_index ) continue;

#ifdef NDEBUG
		core::conformation::parametric::RealValuedParameterCOP curparam( utility::pointer::static_pointer_cast< core::conformation::parametric::RealValuedParameter const>( params->parameter_cop( parameter_type ) ) );
#else
		core::conformation::parametric::RealValuedParameterCOP curparam( utility::pointer::dynamic_pointer_cast< core::conformation::parametric::RealValuedParameter const>( params->parameter_cop( parameter_type ) ) );
		debug_assert( curparam != nullptr );
#endif
		failed = parameter_to_copy_to->copy_value_from_parameter( curparam, params, current_helix_params );
		if ( failed ) return failed;
	}
	return failed;
}

/// @brief Attempts to build a helix based on the current Crick parameters.
/// @details Returns true for FAILURE and false for SUCCESS.  The object "pose" is replaced with a
/// new helix if this operation succeeds.
bool
BundleParametrizationCalculator::build_helix(
	core::pose::Pose &pose
) const {
	//Variables for the start residue and end residue of the helix in pose:
	core::Size const helix_start(1);
	core::Size const helix_end( pose.total_residue() );

	return( build_helix(pose, helix_start, helix_end) );
}

/// @brief Attempts to build a helix based on the current Crick parameters.
/// @details Returns true for FAILURE and false for SUCCESS.  The object "pose" is replaced with a
/// new helix if this operation succeeds.
/// @note This version takes parameters for the start and end of the helix.
bool
BundleParametrizationCalculator::build_helix(
	core::pose::Pose &pose,
	core::Size const helix_start,
	core::Size const helix_end
) const {
	static const std::string errmsg( "Error in protocols::helical_bundle::BundleParametrizationCalculator::build_helix(): " );

	runtime_assert_string_msg( helix_start < helix_end, errmsg + "The first helix residue index must be less than the last helix residue index." );
	runtime_assert_string_msg( helix_start > 0, errmsg + "The first residue of the helix must be greater than zero." );
	runtime_assert_string_msg( helix_start <= pose.total_residue(), errmsg + "The first residue of the helix must be less than or equal to the number of residues in the pose." );
	runtime_assert_string_msg( helix_end > 0, errmsg + "The last residue of the helix must be greater than zero." );
	runtime_assert_string_msg( helix_end <= pose.total_residue(), errmsg + "The last residue of the helix must be less than or equal to the number of residues in the pose." );

	//Generate a vector of vectors of x,y,z coordinates of mainchain atoms.
	//Note that this vector extends one extra residue in each direction.
	utility::vector1 < utility::vector1 < numeric::xyzVector< core::Real > > > atom_positions;
	if ( TR.Debug.visible() ) TR.Debug << "Generating atom positions." << std::endl;
	bool failed(false);
	generate_atom_positions(atom_positions, pose, helix_start, helix_end, real_parameter_cop( BPC_r0 )->value(),
		real_parameter_cop(BPC_omega0)->value(), real_parameter_cop(BPC_delta_omega0)->value(), real_parameter_cop(BPC_delta_t)->value(),
		real_parameter_cop(BPC_z1_offset)->value(), real_parameter_cop(BPC_z0_offset)->value(), real_parameter_cop(BPC_epsilon)->value(),
		boolean_parameter_cop(BPC_invert_helix)->value(), realvector_parameter_cop(BPC_r1_peratom)->value(), real_parameter_cop(BPC_omega1)->value(),
		real_parameter_cop(BPC_z1)->value(), realvector_parameter_cop( BPC_delta_omega1_peratom )->value(), real_parameter_cop(BPC_delta_omega1)->value(),
		realvector_parameter_cop(BPC_delta_z1_peratom)->value(), size_parameter_cop(BPC_residues_per_repeat)->value(),
		sizevector_parameter_cop(BPC_atoms_per_residue)->value(), size_parameter_cop(BPC_repeating_unit_offset)->value(),
		failed );

	if ( failed ) return failed; //Give up now if the parameters don't generate a sensible helix.

	// Make a temporary pose copy and place the atoms using the Crick equations.
	if ( TR.Debug.visible() ) TR.Debug << "Placing atoms." << std::endl;
	core::pose::Pose pose_copy( pose );
	place_atom_positions(pose_copy, atom_positions, helix_start, helix_end);

	//Copy the bond length values from the pose in which we set the x,y,z coordinates of mainchain atoms based on the Crick equations
	//to the ideal pose.
	if ( boolean_parameter_cop( BPC_set_bondlengths )->value() ) {
		if ( TR.Debug.visible() ) TR.Debug << "Copying bond length values." << std::endl;
		copy_helix_bondlengths(pose, pose_copy, helix_start, helix_end);
	}

	//Copy the bond angle values from the pose in which we set the x,y,z coordinates of mainchain atoms based on the Crick equations
	//to the ideal pose.
	if ( boolean_parameter_cop( BPC_set_bondangles )->value() ) {
		if ( TR.Debug.visible() ) TR.Debug << "Copying bond angle values." << std::endl;
		copy_helix_bondangles(pose, pose_copy, helix_start, helix_end);
	}

	//Copy the dihedral values from the pose in which we set the x,y,z coordinates of mainchain atoms based on the Crick equations
	//to the ideal pose.
	if ( boolean_parameter_cop( BPC_set_dihedrals )->value() ) {
		if ( TR.Debug.visible() ) TR.Debug << "Copying dihedral values." << std::endl;
		copy_helix_dihedrals(pose, pose_copy, helix_start, helix_end);
	}

	//Align the ideal pose to the pose in which we set the x,y,z coordinates of mainchain atoms.
	if ( TR.Debug.visible() ) TR.Debug << "Aligning to ideal helix." << std::endl;
	align_mainchain_atoms_of_residue_range(pose, pose_copy, helix_start, helix_end);

	return failed;
}


/// @brief Get the parameter type from the BPC_Parameters enum.
/// @details Returns PT_invalid_type if invalid.
core::conformation::parametric::ParameterType
BundleParametrizationCalculator::parameter_type_from_enum(
	BPC_Parameters param_enum
) {
	using namespace core::conformation::parametric;
	switch( param_enum ) {
	case BPC_r0 :
		return PT_generic_nonnegative_valued_real;
	case BPC_omega0 :
		return PT_angle;
	case BPC_delta_omega0 :
		return PT_angle;
	case BPC_delta_omega1 :
		return PT_angle;
	case BPC_delta_t :
		return PT_generic_real;
	case BPC_z0_offset :
		return PT_generic_real;
	case BPC_z1_offset :
		return PT_generic_real;
	case BPC_epsilon :
		return PT_generic_nonnegative_valued_real;
	case BPC_residues_per_repeat :
		return PT_generic_natural_number;
	case BPC_repeating_unit_offset :
		return PT_generic_whole_number;
	case BPC_atoms_per_residue :
		return PT_generic_natural_number_vector;
	case BPC_r1_peratom :
		return PT_generic_nonnegative_valued_real_vector;
	case BPC_omega1 :
		return PT_angle;
	case BPC_z1 :
		return PT_generic_real;
	case BPC_delta_omega1_peratom :
		return PT_angle_vector;
	case BPC_delta_z1_peratom :
		return PT_generic_real_vector;
	case BPC_invert_helix :
		return PT_boolean;
	case BPC_set_dihedrals :
		return PT_boolean;
	case BPC_set_bondangles :
		return PT_boolean;
	case BPC_set_bondlengths :
		return PT_boolean;
	case BPC_unknown_parameter :
		return PT_invalid_type;
	}
	return PT_invalid_type;
}

/// @brief Given a BPC_Parameters enum, get a short lay-language description (used for annotating output) of the parameter.
/// @details Returns "INVALID!!!" if invalid.
std::string const &
BundleParametrizationCalculator::parameter_description_from_enum(
	BPC_Parameters param_enum
) {
	static const utility::vector1< std::string > descriptions { //Must match the BPC_Parameters enum ordering:
		"Major helix radius, in Angstroms.",
		"Major helix twist per residue, stored in radians.",
		"Rotation of a helix about the z-axis, stored in radians.",
		"Rotation of a helix about its own axis, stored in radians.",
		"Offset along the polypeptide backbone, in residues.",
		"Offset along the global z-axis, in Angstroms.",
		"Offset along the superhelical path through space, in Angstroms.",
		"Lateral squash parameter/eccentricity of the cross-section of a bundle or barrel.",
		"Number of residues per repeating unit in a helix.",
		"Shift, in residues, of the repeating unit of a helix.",
		"Number of mainchain atoms per residue in the repeating unit of a helix -- a vector of integers.",
		"Minor helix radius -- a vector of real numbers in Angstroms, with one per atom in the repeating unit of a helix.  Read from Crick params file, and not normally set by hand.",
		"Minor helix twist per residue, stored in radians.  Read from Crick params file, and not normally set by hand, sampled, or perturbed.",
		"Minor helix rise per residue along the helix axis, in Angstroms.  Read from Crick params file, and not normally set by hand, sampled, or perturbed.",
		"Minor helix angular offsets of each mainchain atom in the repeating unit, in radians.  Read from Crick params file, and not normally set by hand.",
		"Minor helix axial offsets of each mainchain atom in the repeating unit, in Angstroms.  Read from Crick params file, and not normally set by hand.",
		"Inversion state of this helix -- true for inverted.",
		"True indicates that the parametric machinery will set mainchain torsion values.",
		"True indicates that the parametric machinery will allow mainchain bond angle values to deviate from ideality.",
		"True indicates that the parametric machinery will allow mainchain bond length values to deviate from ideality."
		};
	runtime_assert( param_enum < BPC_end_of_list && param_enum > 0);
	return descriptions[param_enum];
}

/// @brief Given a BPC_Parameters enum, get a one- or two-word lay-language description (used for annotating output) of the parameter.
/// @details Returns "INVALID!!!" if invalid.
std::string const &
BundleParametrizationCalculator::short_parameter_description_from_enum(
	BPC_Parameters param_enum
) {
	static const utility::vector1< std::string > descriptions { //Must match the BPC_Parameters enum ordering:
		"Major radius",
		"Major twist",
		"Rotation about z-axis",
		"Roll about helix axis",
		"Registry shift",
		"Major axial offset",
		"Minor axial offset",
		"Lateral squash",
		"Residues/repeat",
		"Repeat unit offset",
		"Atoms/residue",
		"Minor radius",
		"Minor twist",
		"Minor rise",
		"Atom angular offset",
		"Atom axial offset",
		"Invert helix",
		"Dihedrals set",
		"Nonideal bond angles",
		"Nonideal bond lengths"
		};
	runtime_assert( param_enum < BPC_end_of_list && param_enum > 0);
	return descriptions[param_enum];
}

/// @brief Given a BPC_Parameters enum, units for the parameter.
std::string const &
BundleParametrizationCalculator::parameter_units_from_enum(
	BPC_Parameters param_enum
) {
	static const utility::vector1< std::string > descriptions { //Must match the BPC_Parameters enum ordering:
		"Angstroms",
		"radians/residue",
		"radians",
		"radians",
		"residues",
		"Angstroms",
		"vert. Angstroms",
		"dimensionless",
		"dimensionless",
		"residues",
		"dimensionless",
		"Angstroms",
		"radians/residue",
		"Angstroms/residue",
		"radians",
		"Angstroms",
		"Boolean",
		"Boolean",
		"Boolean",
		"Boolean"
		};
	runtime_assert( param_enum < BPC_end_of_list && param_enum > 0);
	return descriptions[param_enum];
}

/// @brief Given a BPC_Parameters enum, get whether this parameter is one that the user can set, whether it can be sampled, whether it can be perturbed, and whether it's global
/// for a ParametersSet (helical bundle).
core::conformation::parametric::ParameterizationCalculatorProperties const &
BundleParametrizationCalculator::parameter_properties_from_enum(
	BPC_Parameters param_enum
) {
	static const utility::vector1< ParameterizationCalculatorProperties > paramproperties { //Must match the BPC_Parameters enum ordering:
		ParameterizationCalculatorProperties( true, true, true, true, false ), //r0
		ParameterizationCalculatorProperties( true, true, true, true, false ), //omega0
		ParameterizationCalculatorProperties( true, true, true, true, false ), //delta_omega0
		ParameterizationCalculatorProperties( true, true, true, true, false ),  //delta_omega1
		ParameterizationCalculatorProperties( true, true, true, true, false ), //delta_t
		ParameterizationCalculatorProperties( true, true, true, true, false ), //z0_offset
		ParameterizationCalculatorProperties( true, true, true, true, false ), //z1_offset
		ParameterizationCalculatorProperties( true, true, true, true, false ), //epsilon
		ParameterizationCalculatorProperties( false, false, false, false, false ), //residues_per_repeat
		ParameterizationCalculatorProperties( true, false, false, false, false ), //repeating_unit_offset
		ParameterizationCalculatorProperties( false, false, false, false, false ), //atoms_per_residue
		ParameterizationCalculatorProperties( true, false, false, false, false ), //r1_peratom
		ParameterizationCalculatorProperties( true, true, true, true, false ), //omega1
		ParameterizationCalculatorProperties( true, true, true, true, false ), //z1
		ParameterizationCalculatorProperties( true, false, false, false, false ), //delta_omega1_peratom
		ParameterizationCalculatorProperties( true, false, false, false, false ), //delta_z1_peratom
		ParameterizationCalculatorProperties( true, false, false, false, false ), //invert_helix
		ParameterizationCalculatorProperties( true, false, false, false, false ), //set_dihedrals
		ParameterizationCalculatorProperties( true, false, false, false, false ), //set_bondangles
		ParameterizationCalculatorProperties( true, false, false, false, false ), //set_bondlengths
		};
	debug_assert( param_enum < BPC_end_of_list && param_enum > 0 );
	return paramproperties[param_enum];
}


/// @brief Get the parameter name from the BPC_Parameters enum.
///
std::string const &
BundleParametrizationCalculator::parameter_name_from_enum(
	BPC_Parameters param_type
) {
	static const utility::vector1< std::string > paramnames { //Must match the BPC_Parameters enum ordering:
		"r0",
		"omega0",
		"delta_omega0",
		"delta_omega1",
		"delta_t",
		"z0_offset",
		"z1_offset",
		"epsilon",
		"residues_per_repeat",
		"repeating_unit_offset",
		"atoms_per_residue",
		"r1_peratom",
		"omega1",
		"z1",
		"delta_omega1_peratom",
		"delta_z1_peratom",
		"invert",
		"set_dihedrals",
		"set_bondangles",
		"set_bondlengths"
		};
	debug_assert( param_type < BPC_end_of_list && param_type > 0 );
	return paramnames[param_type];
}

/// @brief Get the BPC_Parameters enum from the parameter name.
/// @details Returns BPC_unknown_parameter if this can't be parsed.
BPC_Parameters
BundleParametrizationCalculator::parameter_enum_from_name(
	std::string const &name
) {
	for ( core::Size i(1); i < static_cast<core::Size>(BPC_end_of_list); ++i ) {
		if ( !name.compare( parameter_name_from_enum( static_cast<BPC_Parameters>(i) ) ) ) return static_cast<BPC_Parameters>(i);
	}
	return BPC_unknown_parameter;
}

/// @brief Set whether this calculator uses degrees or radians for the user-set angles.
/// @details Defaults to radians (false).
void
BundleParametrizationCalculator::set_use_degrees(
	bool const setting
) {
	use_degrees_ = setting;
	set_use_degrees_for_parameters();
}

/// @brief Set the perturbation type for all parameters.
void
BundleParametrizationCalculator::set_perturbation_type_globally(
	std::string const & perturbation_type_in
) {
	for ( core::Size i(1); i<=static_cast<core::Size>(BPC_last_parameter_to_be_sampled); ++i ) {
		core::conformation::parametric::RealValuedParameterOP curparam( real_parameter( i ) );
		if ( curparam == nullptr ) continue; //Not a real-valued parameter.
		core::conformation::parametric::RealPerturbationType const perttype_enum( core::conformation::parametric::RealValuedParameter::perturbation_type_enum_from_string( perturbation_type_in ) );
		runtime_assert_string_msg( perttype_enum > 0 && perttype_enum < core::conformation::parametric::RPT_unknown_type, "Error in BundleParametrizationCalculator::set_perturbation_type_globally(): Could not interpret perturbation type " + perturbation_type_in + "." );
		curparam->set_perturbation_type( perttype_enum );
	}
}

/// @brief Update whether certain parameters expect their inputs in degrees.
void
BundleParametrizationCalculator::set_use_degrees_for_parameters() {
	real_parameter( BPC_omega0 )->set_input_is_angle_in_degrees( use_degrees_ );
	real_parameter( BPC_delta_omega0 )->set_input_is_angle_in_degrees( use_degrees_ );
	real_parameter( BPC_delta_omega1 )->set_input_is_angle_in_degrees( use_degrees_ );
}

} // namespace helical_bundle
} // namespace protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::helical_bundle::BundleParametrizationCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::conformation::parametric::ParametrizationCalculator >( this ) );
	arc( CEREAL_NVP( use_degrees_ ) ); //bool
	//TODO -- archive private member variables here.
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::helical_bundle::BundleParametrizationCalculator::load( Archive & arc ) {
	arc( cereal::base_class< core::conformation::parametric::ParametrizationCalculator >( this ) );
	arc( use_degrees_ ); //bool
	//TODO -- de-archive private member variables here.
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::helical_bundle::BundleParametrizationCalculator );
CEREAL_REGISTER_TYPE( protocols::helical_bundle::BundleParametrizationCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_helical_bundle_BundleParametrizationCalculator )
#endif // SERIALIZATION
