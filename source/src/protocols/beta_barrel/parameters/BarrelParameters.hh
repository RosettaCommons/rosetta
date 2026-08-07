// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/parameters/BarrelParameters.hh
/// @brief  Parameters for a single strand in a parametric beta-barrel.
/// @author Andy Watkins

#ifndef INCLUDED_protocols_beta_barrel_parameters_BarrelParameters_hh
#define INCLUDED_protocols_beta_barrel_parameters_BarrelParameters_hh

#include <protocols/beta_barrel/parameters/BarrelParameters.fwd.hh>
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/ParametersSet.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.hh>

#ifdef    SERIALIZATION
#include <cereal/types/polymorphic.fwd.hpp>
#endif

namespace protocols {
namespace beta_barrel {
namespace parameters {

class BarrelParameters : public core::conformation::parametric::Parameters
{
public:

	typedef core::conformation::parametric::Parameters Parameters;
	typedef core::conformation::parametric::ParametersOP ParametersOP;

public:

	BarrelParameters();
	BarrelParameters( BarrelParameters const & src );
	~BarrelParameters() override;

	ParametersOP clone() const override;

	virtual void get_pdb_remark( std::stringstream & remark ) const;

private:

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif

}; // class BarrelParameters

} // namespace parameters
} // namespace beta_barrel
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_beta_barrel_parameters_BarrelParameters )
#endif

#endif
