// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/BarrelParametrizationCalculator.hh
/// @brief  Parametrization calculator for beta-barrel backbone generation.
/// @author Andy Watkins

#ifndef INCLUDED_protocols_beta_barrel_BarrelParametrizationCalculator_hh
#define INCLUDED_protocols_beta_barrel_BarrelParametrizationCalculator_hh

#include <protocols/beta_barrel/BarrelParametrizationCalculator.fwd.hh>
#include <core/conformation/parametric/ParametrizationCalculator.hh>

#include <protocols/beta_barrel/parameters/BarrelParametersSet.fwd.hh>
#include <protocols/beta_barrel/parameters/BarrelParameters.fwd.hh>
#include <core/conformation/parametric/RealValuedParameter.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#ifdef    SERIALIZATION
#include <cereal/types/polymorphic.fwd.hpp>
#endif

namespace protocols {
namespace beta_barrel {

/// @brief Parameters for a single strand in a beta-barrel.
enum BBPC_Parameters {
	// Sampled/perturbable parameters:
	BBPC_r0 = 1,           // Barrel radius (Angstroms), effectively global
	BBPC_delta_omega0,      // Azimuthal position of this strand (radians)
	BBPC_delta_z0,          // Axial offset of this strand (Angstroms)
	BBPC_delta_omega1,      // Rotational offset about the strand's own axis (radians)
	BBPC_last_parameter_to_be_sampled = BBPC_delta_omega1,

	// Invariant minor helix parameters (from .crick_params):
	BBPC_residues_per_repeat,
	BBPC_atoms_per_residue,
	BBPC_r1_peratom,
	BBPC_omega1,
	BBPC_z1,
	BBPC_delta_omega1_peratom,
	BBPC_delta_z1_peratom,

	// Control flags:
	BBPC_invert_strand,
	BBPC_set_dihedrals,
	BBPC_set_bondangles,
	BBPC_set_bondlengths,

	BBPC_unknown_parameter,
	BBPC_end_of_list = BBPC_unknown_parameter
};

class BarrelParametrizationCalculator : public core::conformation::parametric::ParametrizationCalculator
{
public:

	typedef core::conformation::parametric::ParameterizationCalculatorProperties ParameterizationCalculatorProperties;

public:

	BarrelParametrizationCalculator( bool const use_degrees = false );

	BarrelParametrizationCalculator( bool const use_degrees, parameters::BarrelParametersCOP params_in );

	BarrelParametrizationCalculator( BarrelParametrizationCalculator const & src );

	~BarrelParametrizationCalculator() override;

	core::conformation::parametric::ParametrizationCalculatorOP clone() const override;

public:

	void init_from_file( std::string const & filename );

	void copy_unset_params_from_globals( BarrelParametrizationCalculatorCOP global_calculator );

	/// @brief Build a single strand based on current parameters.
	/// @returns true for FAILURE, false for SUCCESS.
	bool build_strand( core::pose::Pose & pose ) const;

	/// @brief Build a single strand given explicit residue range.
	/// @returns true for FAILURE, false for SUCCESS.
	bool build_strand( core::pose::Pose & pose, core::Size strand_start, core::Size strand_end ) const;

	// Static metadata functions:
	static core::conformation::parametric::ParameterType parameter_type_from_enum( BBPC_Parameters param_enum );
	static std::string const & parameter_description_from_enum( BBPC_Parameters param_enum );
	static std::string const & short_parameter_description_from_enum( BBPC_Parameters param_enum );
	static std::string const & parameter_units_from_enum( BBPC_Parameters param_enum );
	static ParameterizationCalculatorProperties const & parameter_properties_from_enum( BBPC_Parameters param_enum );
	static std::string const & parameter_name_from_enum( BBPC_Parameters param_enum );
	static BBPC_Parameters parameter_enum_from_name( std::string const & name );

	void set_use_degrees( bool const setting );
	inline bool use_degrees() const { return use_degrees_; }

	void set_perturbation_type_globally( std::string const & perturbation_type_in );

	core::Size residues_per_repeat() const;

private:

	void set_use_degrees_for_parameters();

private:

	bool use_degrees_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif

}; // class BarrelParametrizationCalculator

} // namespace beta_barrel
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_beta_barrel_BarrelParametrizationCalculator )
#endif

#endif
