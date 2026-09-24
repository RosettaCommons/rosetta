// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/parameters/BarrelParametersSet.cc
/// @brief  Container for all strands in a parametric beta-barrel.
/// @author Andy Watkins

#include <protocols/beta_barrel/parameters/BarrelParametersSet.hh>
#include <protocols/beta_barrel/parameters/BarrelParameters.hh>

#include <core/conformation/parametric/ParametersSet.fwd.hh>

#include <basic/Tracer.hh>

#ifdef    SERIALIZATION
#include <utility/serialization/serialization.hh>
#include <cereal/types/polymorphic.hpp>
#endif

namespace protocols {
namespace beta_barrel {
namespace parameters {

static basic::Tracer TR( "protocols.beta_barrel.parameters.BarrelParametersSet" );

BarrelParametersSet::BarrelParametersSet() :
	n_strands_(0),
	shear_number_(0),
	antiparallel_(false),
	barrel_radius_(0.0)
{}

BarrelParametersSet::BarrelParametersSet( BarrelParametersSet const & /*src*/ ) = default;

BarrelParametersSet::~BarrelParametersSet() = default;

core::conformation::parametric::ParametersSetOP
BarrelParametersSet::clone() const
{
	return utility::pointer::make_shared< BarrelParametersSet >( *this );
}

void BarrelParametersSet::get_pdb_remark( std::stringstream & remark ) const {
	remark << "----------------------------------------" << std::endl;
	remark << "BETA-BARREL PARAMETER INFORMATION" << std::endl;
	remark << "Number of strands: " << n_strands() << std::endl;
	remark << "Shear number: " << shear_number() << std::endl;
	remark << "Topology: " << ( antiparallel() ? "antiparallel" : "all-parallel" ) << std::endl;
	remark << "Barrel radius: " << barrel_radius() << " Angstroms" << std::endl;
	if ( n_strands() > 0 ) {
		for ( core::Size i = 1, imax = n_parameters(); i <= imax; ++i ) {
			BarrelParametersCOP cur_strand(
				utility::pointer::static_pointer_cast< BarrelParameters const >( parameters(i) ) );
			remark << "---STRAND " << i << " PARAMETERS:---" << std::endl;
			cur_strand->get_pdb_remark( remark );
		}
	}
	remark << "----------------------------------------" << std::endl;
}

} // namespace parameters
} // namespace beta_barrel
} // namespace protocols


#ifdef    SERIALIZATION

template< class Archive >
void
protocols::beta_barrel::parameters::BarrelParametersSet::save( Archive & arc ) const {
	arc( cereal::base_class< core::conformation::parametric::ParametersSet >( this ) );
	arc( CEREAL_NVP( n_strands_ ) );
	arc( CEREAL_NVP( shear_number_ ) );
	arc( CEREAL_NVP( antiparallel_ ) );
	arc( CEREAL_NVP( barrel_radius_ ) );
}

template< class Archive >
void
protocols::beta_barrel::parameters::BarrelParametersSet::load( Archive & arc ) {
	arc( cereal::base_class< core::conformation::parametric::ParametersSet >( this ) );
	arc( n_strands_ );
	arc( shear_number_ );
	arc( antiparallel_ );
	arc( barrel_radius_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::beta_barrel::parameters::BarrelParametersSet );
CEREAL_REGISTER_TYPE( protocols::beta_barrel::parameters::BarrelParametersSet )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_beta_barrel_parameters_BarrelParametersSet )
#endif
