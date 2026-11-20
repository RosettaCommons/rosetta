// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/BarrelParametrizationCalculator.cc
/// @brief  Parametrization calculator for beta-barrel backbone generation.
/// @author Andy Watkins

#include <protocols/beta_barrel/BarrelParametrizationCalculator.hh>
#include <protocols/beta_barrel/parameters/BarrelParameters.hh>

#include <core/conformation/parametric/RealValuedParameter.hh>
#include <core/conformation/parametric/RealVectorValuedParameter.hh>
#include <core/conformation/parametric/SizeValuedParameter.hh>
#include <core/conformation/parametric/SizeVectorValuedParameter.hh>
#include <core/conformation/parametric/BooleanValuedParameter.hh>

#include <protocols/helical_bundle/util.hh>

#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <utility/pointer/memory.hh>
#include <utility/exit.hh>

#include <numeric/constants.hh>
#include <cmath>

#ifdef    SERIALIZATION
#include <utility/serialization/serialization.hh>
#include <cereal/types/polymorphic.hpp>
#endif

namespace protocols {
namespace beta_barrel {

static basic::Tracer TR( "protocols.beta_barrel.BarrelParametrizationCalculator" );

BarrelParametrizationCalculator::BarrelParametrizationCalculator( bool const use_degrees ) :
	core::conformation::parametric::ParametrizationCalculator(
		utility::pointer::make_shared< parameters::BarrelParameters >() ),
	use_degrees_( use_degrees )
{
	using namespace core::conformation::parametric;
	for ( core::Size i(1); i < static_cast<core::Size>(BBPC_end_of_list); ++i ) {
		add_parameter(
			parameter_name_from_enum( static_cast<BBPC_Parameters>(i) ),
			parameter_type_from_enum( static_cast<BBPC_Parameters>(i) ),
			parameter_description_from_enum( static_cast<BBPC_Parameters>(i) ),
			short_parameter_description_from_enum( static_cast<BBPC_Parameters>(i) ),
			parameter_units_from_enum( static_cast<BBPC_Parameters>(i) ),
			parameter_properties_from_enum( static_cast<BBPC_Parameters>(i) )
		);
		parameter(i)->set_copy_suffix("copies_strand");
	}
	set_use_degrees_for_parameters();
	real_parameter( BBPC_epsilon )->set_default_value( 1.0 );
	boolean_parameter( BBPC_set_dihedrals )->set_default_value( true );
	boolean_parameter( BBPC_set_bondangles )->set_default_value( true );
	boolean_parameter( BBPC_set_bondlengths )->set_default_value( true );
}

BarrelParametrizationCalculator::BarrelParametrizationCalculator(
	bool const use_degrees,
	parameters::BarrelParametersCOP params_in
) :
	core::conformation::parametric::ParametrizationCalculator( params_in->clone() ),
	use_degrees_( use_degrees )
{}

BarrelParametrizationCalculator::BarrelParametrizationCalculator(
	BarrelParametrizationCalculator const & src
) :
	core::conformation::parametric::ParametrizationCalculator( src ),
	use_degrees_( src.use_degrees_ )
{}

BarrelParametrizationCalculator::~BarrelParametrizationCalculator() = default;

core::conformation::parametric::ParametrizationCalculatorOP
BarrelParametrizationCalculator::clone() const {
	return utility::pointer::make_shared< BarrelParametrizationCalculator >( *this );
}

void
BarrelParametrizationCalculator::init_from_file(
	std::string const & filename
) {
	utility::vector1< core::Real > r1;
	core::Real omega1( 0.0 );
	core::Real z1( 0.0 );
	utility::vector1< core::Real > delta_omega1;
	utility::vector1< core::Real > delta_z1;
	core::Size residues_per_repeat( 1 );
	utility::vector1< core::Size > atoms_per_residue;
	protocols::helical_bundle::read_minor_helix_params(
		filename, r1, omega1, z1, delta_omega1, delta_z1, residues_per_repeat, atoms_per_residue );
	realvector_parameter( BBPC_r1_peratom )->set_value( r1 );
	real_parameter( BBPC_omega1 )->set_value( omega1 );
	real_parameter( BBPC_z1 )->set_value( z1 );
	realvector_parameter( BBPC_delta_omega1_peratom )->set_value( delta_omega1 );
	realvector_parameter( BBPC_delta_z1_peratom )->set_value( delta_z1 );
	size_parameter( BBPC_residues_per_repeat )->set_value( residues_per_repeat );
	sizevector_parameter( BBPC_atoms_per_residue )->set_value( atoms_per_residue );
}

void
BarrelParametrizationCalculator::copy_unset_params_from_globals(
	BarrelParametrizationCalculatorCOP global_calculator
) {
	for ( core::Size i(1); i < static_cast<core::Size>(BBPC_end_of_list); ++i ) {
		core::conformation::parametric::ParameterOP curparam( parameter(i) );
		if ( !curparam->value_was_set() ) {
			curparam->copy_value_from_parameter(
				global_calculator->parameter_cop(i),
				global_calculator->parameters_cop(),
				this->parameters_cop() );
		}
	}
}

bool
BarrelParametrizationCalculator::build_strand(
	core::pose::Pose & pose
) const {
	core::Size const strand_start( 1 );
	core::Size const strand_end( pose.total_residue() );
	return build_strand( pose, strand_start, strand_end );
}

bool
BarrelParametrizationCalculator::build_strand(
	core::pose::Pose & pose,
	core::Size const strand_start,
	core::Size const strand_end
) const {
	static const std::string errmsg( "Error in protocols::beta_barrel::BarrelParametrizationCalculator::build_strand(): " );

	runtime_assert_string_msg( strand_start < strand_end, errmsg + "The first strand residue index must be less than the last." );
	runtime_assert_string_msg( strand_start > 0 && strand_start <= pose.total_residue(), errmsg + "Invalid strand start index." );
	runtime_assert_string_msg( strand_end > 0 && strand_end <= pose.total_residue(), errmsg + "Invalid strand end index." );

	utility::vector1< utility::vector1< numeric::xyzVector< core::Real > > > atom_positions;
	if ( TR.Debug.visible() ) TR.Debug << "Generating strand atom positions." << std::endl;
	bool failed( false );

	// omega0=0 for a closed barrel (strands parallel to barrel axis).
	// omega0!=0 for a solenoid/beta-helix (strands spiral around the axis).
	// epsilon=1.0 for circular cross-section; !=1 for elliptical.
	core::Real const z1_offset = 0.0;
	core::Size const repeating_unit_offset = 0;

	protocols::helical_bundle::generate_atom_positions(
		atom_positions, pose, strand_start, strand_end,
		real_parameter_cop( BBPC_r0 )->value(),
		real_parameter_cop( BBPC_omega0 )->value(),
		real_parameter_cop( BBPC_delta_omega0 )->value(),
		real_parameter_cop( BBPC_delta_t )->value(),
		z1_offset,
		real_parameter_cop( BBPC_delta_z0 )->value(),
		real_parameter_cop( BBPC_epsilon )->value(),
		boolean_parameter_cop( BBPC_invert_strand )->value(),
		realvector_parameter_cop( BBPC_r1_peratom )->value(),
		real_parameter_cop( BBPC_omega1 )->value(),
		real_parameter_cop( BBPC_z1 )->value(),
		realvector_parameter_cop( BBPC_delta_omega1_peratom )->value(),
		real_parameter_cop( BBPC_delta_omega1 )->value(),
		realvector_parameter_cop( BBPC_delta_z1_peratom )->value(),
		size_parameter_cop( BBPC_residues_per_repeat )->value(),
		sizevector_parameter_cop( BBPC_atoms_per_residue )->value(),
		repeating_unit_offset,
		failed );

	if ( failed ) return failed;

	if ( TR.Debug.visible() ) TR.Debug << "Placing strand atoms." << std::endl;
	core::pose::Pose pose_copy( pose );
	protocols::helical_bundle::place_atom_positions( pose_copy, atom_positions, strand_start, strand_end );

	if ( boolean_parameter_cop( BBPC_set_bondlengths )->value() ) {
		protocols::helical_bundle::copy_helix_bondlengths( pose, pose_copy, strand_start, strand_end );
	}
	if ( boolean_parameter_cop( BBPC_set_bondangles )->value() ) {
		protocols::helical_bundle::copy_helix_bondangles( pose, pose_copy, strand_start, strand_end );
	}
	if ( boolean_parameter_cop( BBPC_set_dihedrals )->value() ) {
		protocols::helical_bundle::copy_helix_dihedrals( pose, pose_copy, strand_start, strand_end );
	}
	protocols::helical_bundle::align_mainchain_atoms_of_residue_range( pose, pose_copy, strand_start, strand_end );

	return failed;
}


////////////////////////////////////////////////////////////////////////////////
//          STATIC METADATA FUNCTIONS                                        //
////////////////////////////////////////////////////////////////////////////////

core::conformation::parametric::ParameterType
BarrelParametrizationCalculator::parameter_type_from_enum(
	BBPC_Parameters param_enum
) {
	using namespace core::conformation::parametric;
	switch( param_enum ) {
	case BBPC_r0:                    return PT_generic_nonnegative_valued_real;
	case BBPC_omega0:                return PT_angle;
	case BBPC_delta_omega0:          return PT_angle;
	case BBPC_delta_z0:              return PT_generic_real;
	case BBPC_delta_omega1:          return PT_angle;
	case BBPC_delta_t:               return PT_generic_real;
	case BBPC_epsilon:               return PT_generic_nonnegative_valued_real;
	case BBPC_residues_per_repeat:   return PT_generic_natural_number;
	case BBPC_atoms_per_residue:     return PT_generic_natural_number_vector;
	case BBPC_r1_peratom:            return PT_generic_nonnegative_valued_real_vector;
	case BBPC_omega1:                return PT_angle;
	case BBPC_z1:                    return PT_generic_real;
	case BBPC_delta_omega1_peratom:  return PT_angle_vector;
	case BBPC_delta_z1_peratom:      return PT_generic_real_vector;
	case BBPC_invert_strand:         return PT_boolean;
	case BBPC_set_dihedrals:         return PT_boolean;
	case BBPC_set_bondangles:        return PT_boolean;
	case BBPC_set_bondlengths:       return PT_boolean;
	case BBPC_unknown_parameter:     return PT_invalid_type;
	}
	return PT_invalid_type;
}

std::string const &
BarrelParametrizationCalculator::parameter_description_from_enum(
	BBPC_Parameters param_enum
) {
	static const utility::vector1< std::string > descriptions {
		"Barrel/solenoid radius in Angstroms.",
		"Superhelical twist per residue, stored in radians. Zero for closed barrels; nonzero for solenoids/beta-helices.",
		"Azimuthal position of this strand around the barrel axis, stored in radians.",
		"Axial offset of this strand along the barrel axis, in Angstroms.",
		"Rotational offset of this strand about its own axis, stored in radians.",
		"Offset along the polypeptide backbone, in residues.",
		"Lateral squash parameter/eccentricity of the barrel cross-section.",
		"Number of residues per repeating unit in a strand.",
		"Number of mainchain atoms per residue in the repeating unit -- a vector of integers.",
		"Minor helix radius per atom, in Angstroms. Read from Crick params file.",
		"Minor helix twist per residue, stored in radians. Read from Crick params file.",
		"Minor helix rise per residue along the strand axis, in Angstroms. Read from Crick params file.",
		"Per-atom angular offsets in the repeating unit, in radians. Read from Crick params file.",
		"Per-atom axial offsets in the repeating unit, in Angstroms. Read from Crick params file.",
		"Inversion state of this strand -- true for inverted (antiparallel).",
		"True indicates that the parametric machinery will set mainchain torsion values.",
		"True indicates that the parametric machinery will allow mainchain bond angle values to deviate from ideality.",
		"True indicates that the parametric machinery will allow mainchain bond length values to deviate from ideality."
	};
	runtime_assert( param_enum < BBPC_end_of_list && param_enum > 0 );
	return descriptions[param_enum];
}

std::string const &
BarrelParametrizationCalculator::short_parameter_description_from_enum(
	BBPC_Parameters param_enum
) {
	static const utility::vector1< std::string > descriptions {
		"Barrel radius",
		"Superhelical twist",
		"Azimuthal position",
		"Axial offset",
		"Roll about strand axis",
		"Registry shift",
		"Lateral squash",
		"Residues/repeat",
		"Atoms/residue",
		"Minor radius",
		"Minor twist",
		"Minor rise",
		"Atom angular offset",
		"Atom axial offset",
		"Invert strand",
		"Dihedrals set",
		"Nonideal bond angles",
		"Nonideal bond lengths"
	};
	runtime_assert( param_enum < BBPC_end_of_list && param_enum > 0 );
	return descriptions[param_enum];
}

std::string const &
BarrelParametrizationCalculator::parameter_units_from_enum(
	BBPC_Parameters param_enum
) {
	static const utility::vector1< std::string > units {
		"Angstroms",
		"radians/residue",
		"radians",
		"Angstroms",
		"radians",
		"residues",
		"dimensionless",
		"dimensionless",
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
	runtime_assert( param_enum < BBPC_end_of_list && param_enum > 0 );
	return units[param_enum];
}

core::conformation::parametric::ParameterizationCalculatorProperties const &
BarrelParametrizationCalculator::parameter_properties_from_enum(
	BBPC_Parameters param_enum
) {
	static const utility::vector1< ParameterizationCalculatorProperties > props {
		//                                                        can_set, can_copy, can_sample, can_perturb, global
		ParameterizationCalculatorProperties( true,  true,  true,  true,  true  ), // r0 (global for all strands)
		ParameterizationCalculatorProperties( true,  true,  true,  true,  true  ), // omega0 (global; 0 for barrel, nonzero for solenoid)
		ParameterizationCalculatorProperties( true,  true,  true,  true,  false ), // delta_omega0
		ParameterizationCalculatorProperties( true,  true,  true,  true,  false ), // delta_z0
		ParameterizationCalculatorProperties( true,  true,  true,  true,  false ), // delta_omega1
		ParameterizationCalculatorProperties( true,  true,  true,  true,  false ), // delta_t
		ParameterizationCalculatorProperties( true,  true,  true,  true,  true  ), // epsilon (global)
		ParameterizationCalculatorProperties( false, false, false, false, false ), // residues_per_repeat
		ParameterizationCalculatorProperties( false, false, false, false, false ), // atoms_per_residue
		ParameterizationCalculatorProperties( true,  false, false, false, false ), // r1_peratom
		ParameterizationCalculatorProperties( true,  false, false, false, false ), // omega1
		ParameterizationCalculatorProperties( true,  false, false, false, false ), // z1
		ParameterizationCalculatorProperties( true,  false, false, false, false ), // delta_omega1_peratom
		ParameterizationCalculatorProperties( true,  false, false, false, false ), // delta_z1_peratom
		ParameterizationCalculatorProperties( true,  false, false, false, false ), // invert_strand
		ParameterizationCalculatorProperties( true,  false, false, false, false ), // set_dihedrals
		ParameterizationCalculatorProperties( true,  false, false, false, false ), // set_bondangles
		ParameterizationCalculatorProperties( true,  false, false, false, false ), // set_bondlengths
	};
	debug_assert( param_enum < BBPC_end_of_list && param_enum > 0 );
	return props[param_enum];
}

std::string const &
BarrelParametrizationCalculator::parameter_name_from_enum(
	BBPC_Parameters param_enum
) {
	static const utility::vector1< std::string > names {
		"r0",
		"omega0",
		"delta_omega0",
		"delta_z0",
		"delta_omega1",
		"delta_t",
		"epsilon",
		"residues_per_repeat",
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
	debug_assert( param_enum < BBPC_end_of_list && param_enum > 0 );
	return names[param_enum];
}

BBPC_Parameters
BarrelParametrizationCalculator::parameter_enum_from_name(
	std::string const & name
) {
	for ( core::Size i(1); i < static_cast<core::Size>(BBPC_end_of_list); ++i ) {
		if ( name == parameter_name_from_enum( static_cast<BBPC_Parameters>(i) ) ) return static_cast<BBPC_Parameters>(i);
	}
	return BBPC_unknown_parameter;
}

void
BarrelParametrizationCalculator::set_use_degrees(
	bool const setting
) {
	use_degrees_ = setting;
	set_use_degrees_for_parameters();
}

void
BarrelParametrizationCalculator::set_perturbation_type_globally(
	std::string const & perturbation_type_in
) {
	for ( core::Size i(1); i <= static_cast<core::Size>(BBPC_last_parameter_to_be_sampled); ++i ) {
		core::conformation::parametric::RealValuedParameterOP curparam( real_parameter(i) );
		if ( curparam == nullptr ) continue;
		core::conformation::parametric::RealPerturbationType const perttype_enum(
			core::conformation::parametric::RealValuedParameter::perturbation_type_enum_from_string( perturbation_type_in ) );
		runtime_assert_string_msg( perttype_enum > 0 && perttype_enum < core::conformation::parametric::RPT_unknown_type,
			"Error in BarrelParametrizationCalculator::set_perturbation_type_globally(): Could not interpret perturbation type " + perturbation_type_in + "." );
		curparam->set_perturbation_type( perttype_enum );
	}
}

core::Size
BarrelParametrizationCalculator::residues_per_repeat() const {
	return size_parameter_cop( BBPC_residues_per_repeat )->value();
}

void
BarrelParametrizationCalculator::set_use_degrees_for_parameters() {
	real_parameter( BBPC_omega0 )->set_input_is_angle_in_degrees( use_degrees_ );
	real_parameter( BBPC_delta_omega0 )->set_input_is_angle_in_degrees( use_degrees_ );
	real_parameter( BBPC_delta_omega1 )->set_input_is_angle_in_degrees( use_degrees_ );
}

} // namespace beta_barrel
} // namespace protocols

#ifdef    SERIALIZATION

template< class Archive >
void
protocols::beta_barrel::BarrelParametrizationCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::conformation::parametric::ParametrizationCalculator >( this ) );
	arc( CEREAL_NVP( use_degrees_ ) );
}

template< class Archive >
void
protocols::beta_barrel::BarrelParametrizationCalculator::load( Archive & arc ) {
	arc( cereal::base_class< core::conformation::parametric::ParametrizationCalculator >( this ) );
	arc( use_degrees_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::beta_barrel::BarrelParametrizationCalculator );
CEREAL_REGISTER_TYPE( protocols::beta_barrel::BarrelParametrizationCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_beta_barrel_BarrelParametrizationCalculator )
#endif
