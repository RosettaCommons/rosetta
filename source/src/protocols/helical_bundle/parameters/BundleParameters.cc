// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/helical_bundle/parameters/BundleParameters.cc
/// @brief  A class for holding parameters for parametric helical bundle generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit header
#include <protocols/helical_bundle/parameters/BundleParameters.hh>

// Package headers
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.hh>

// Project headers

// Basic headers
#include <basic/basic.hh>
#include <basic/Tracer.hh>

// Numeric headers

// Utility Headers
#include <utility/assert.hh>
#include <utility/py/PyAssert.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <iostream>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace helical_bundle {
namespace parameters {

static basic::Tracer TR( "protocols.helical_bundle.parameters.BundleParameters" );

/// @brief Constructor.
///
BundleParameters::BundleParameters() :
	r0_(0.0),
	omega0_(0.0),
	delta_omega0_(0.0),
	residues_per_repeat_(1),
	repeating_unit_offset_(0),
	atoms_per_residue_(),
	r1_(),
	omega1_(0.0),
	delta_omega1_all_(0.0),
	z1_(0.0),
	delta_omega1_(),
	delta_z1_(),
	z1_offset_(0.0),
	z0_offset_(0.0),
	epsilon_(1.0),
	invert_helix_(false),
	delta_t_(0.0),
	allow_dihedrals_(true),
	allow_bondangles_(true),
	allow_bondlengths_(true)
{
}

BundleParameters::BundleParameters( BundleParameters const & src ) :
	core::conformation::parametric::Parameters( src ),
	r0_(src.r0()),
	omega0_(src.omega0()),
	delta_omega0_(src.delta_omega0()),
	residues_per_repeat_(src.residues_per_repeat()),
	repeating_unit_offset_( src.repeating_unit_offset() ),
	atoms_per_residue_(src.atoms_per_residue_),
	r1_(src.r1_),
	omega1_(src.omega1()),
	delta_omega1_all_(src.delta_omega1_all()),
	z1_(src.z1()),
	delta_omega1_(src.delta_omega1_),
	delta_z1_(src.delta_z1_),
	z1_offset_(src.z1_offset_),
	z0_offset_(src.z0_offset_),
	epsilon_(src.epsilon_),
	invert_helix_(src.invert_helix()),
	delta_t_(src.delta_t()),
	allow_dihedrals_(src.allow_dihedrals()),
	allow_bondangles_(src.allow_bondangles()),
	allow_bondlengths_(src.allow_bondlengths())
{
}

BundleParameters::~BundleParameters() {}


/// @brief make a copy of this residue( allocate actual memory for it )
///
core::conformation::parametric::ParametersOP
BundleParameters::clone() const
{
	return ParametersOP( new BundleParameters( *this ) );
}

/// @brief Get a summary of this ParametersSet object, for output to remark lines of a PDB file.
///
void BundleParameters::get_pdb_remark(std::stringstream &remark) const {
	remark.setf( std::ios::fixed, std::ios::floatfield );
	remark.precision(8);
	remark << " MAJOR HELIX PARAMETERS:" << std::endl;
	remark << "   Radius (r0,Angstroms): " << r0() << std::endl;
	remark << "   Twist (omega0,radians/residue): " << omega0() << std::endl;
	remark << "   Rotational offset (delta_omega0,radians): " << delta_omega0() << std::endl;
	remark << "   Axial offset (z0_offset,Angstroms): " << z0_offset() << std::endl;
	remark << "   Lateral squash factor (epsilon,dimensionless): " << epsilon() << std::endl;
	remark << "   Invert helix: " << (invert_helix() ? "true" : "false") << std::endl;
	remark << "   Dihedrals set by generator: " << (allow_dihedrals() ? "true" : "false") << std::endl;
	remark << "   Bond angles set by generator: " << (allow_bondangles() ? "true" : "false") << std::endl;
	remark << "   Bond lengths set by generator: " << (allow_bondlengths() ? "true" : "false") << std::endl;
	remark << " MINOR HELIX PARAMETERS (that are typically sampled):" << std::endl;
	remark << "   Roll about axis (delta_omega1,radians): " << delta_omega1_all() << std::endl;
	remark << "   Axial offset (z1_offset,vert. Angstroms): " << z1_offset() << std::endl;
	remark << "   Registry shift (delta_t,residues): " << delta_t() << std::endl;

	remark << " OTHER MINOR HELIX PARAMETERS (fixed):" << std::endl;
	remark << "   Residues/repeat: " << residues_per_repeat() << std::endl;
	for ( core::Size i=1, imax=residues_per_repeat(); i<=imax; ++i ) {
		remark << "   Atoms/residue" << i << ": " << atoms_per_residue(i) << std::endl;
	}
	remark << "   Repeat unit offset: " << repeating_unit_offset() << std::endl;
	remark << "   Twist (omega1,radians/residue): " << omega1() << std::endl;
	remark << "   Rise (z1,Angstroms/residue): " << z1() << std::endl;
	for ( core::Size i=1, imax=r1_.size(); i<=imax; ++i ) {
		remark << "   MAINCHAIN ATOM #" << i << ":" << std::endl;
		remark << "   Radius (r1,Angstroms): " << r1(i) << std::endl;
		remark << "   Axial offset (delta_z1,Angstroms): " << delta_z1(i) << std::endl;
	}
}

} // namespace parameters
} // namespace helical_bundle
} // namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::helical_bundle::parameters::BundleParameters::save( Archive & arc ) const {
	arc( cereal::base_class< core::conformation::parametric::Parameters >( this ) );
	arc( CEREAL_NVP( r0_ ) ); // core::Real
	arc( CEREAL_NVP( omega0_ ) ); // core::Real
	arc( CEREAL_NVP( delta_omega0_ ) ); // core::Real
	arc( CEREAL_NVP( residues_per_repeat_ ) ); // core::Size
	arc( CEREAL_NVP( repeating_unit_offset_ ) ); // core::Size
	arc( CEREAL_NVP( atoms_per_residue_ ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( r1_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( omega1_ ) ); // core::Real
	arc( CEREAL_NVP( delta_omega1_all_ ) ); // core::Real
	arc( CEREAL_NVP( z1_ ) ); // core::Real
	arc( CEREAL_NVP( delta_omega1_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( delta_z1_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( z1_offset_ ) ); // core::Real
	arc( CEREAL_NVP( z0_offset_ ) ); // core::Real
	arc( CEREAL_NVP( epsilon_ ) ); // core::Real
	arc( CEREAL_NVP( invert_helix_ ) ); // _Bool
	arc( CEREAL_NVP( delta_t_ ) ); // core::Real
	arc( CEREAL_NVP( allow_dihedrals_ ) ); // _Bool
	arc( CEREAL_NVP( allow_bondangles_ ) ); // _Bool
	arc( CEREAL_NVP( allow_bondlengths_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::helical_bundle::parameters::BundleParameters::load( Archive & arc ) {
	arc( cereal::base_class< core::conformation::parametric::Parameters >( this ) );
	arc( r0_ ); // core::Real
	arc( omega0_ ); // core::Real
	arc( delta_omega0_ ); // core::Real
	arc( residues_per_repeat_ ); // core::Size
	arc( repeating_unit_offset_ ); // core::Size
	arc( atoms_per_residue_ ); // utility::vector1<core::Size>
	arc( r1_ ); // utility::vector1<core::Real>
	arc( omega1_ ); // core::Real
	arc( delta_omega1_all_ ); // core::Real
	arc( z1_ ); // core::Real
	arc( delta_omega1_ ); // utility::vector1<core::Real>
	arc( delta_z1_ ); // utility::vector1<core::Real>
	arc( z1_offset_ ); // core::Real
	arc( z0_offset_ ); // core::Real
	arc( epsilon_ ); // core::Real
	arc( invert_helix_ ); // _Bool
	arc( delta_t_ ); // core::Real
	arc( allow_dihedrals_ ); // _Bool
	arc( allow_bondangles_ ); // _Bool
	arc( allow_bondlengths_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::helical_bundle::parameters::BundleParameters );
CEREAL_REGISTER_TYPE( protocols::helical_bundle::parameters::BundleParameters )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_helical_bundle_parameters_BundleParameters )
#endif // SERIALIZATION
