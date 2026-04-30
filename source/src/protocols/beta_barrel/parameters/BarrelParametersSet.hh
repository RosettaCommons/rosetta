// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/parameters/BarrelParametersSet.hh
/// @brief  Container for all strands in a parametric beta-barrel.
/// @author Andy Watkins

#ifndef INCLUDED_protocols_beta_barrel_parameters_BarrelParametersSet_hh
#define INCLUDED_protocols_beta_barrel_parameters_BarrelParametersSet_hh

#include <protocols/beta_barrel/parameters/BarrelParametersSet.fwd.hh>
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/ParametersSet.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.hh>
#include <core/types.hh>

#ifdef    SERIALIZATION
#include <cereal/types/polymorphic.fwd.hpp>
#endif

namespace protocols {
namespace beta_barrel {
namespace parameters {

class BarrelParametersSet : public core::conformation::parametric::ParametersSet
{
public:

	typedef core::conformation::parametric::Parameters Parameters;
	typedef core::conformation::parametric::ParametersOP ParametersOP;
	typedef core::conformation::parametric::ParametersSet ParametersSet;
	typedef core::conformation::parametric::ParametersSetOP ParametersSetOP;

public:

	BarrelParametersSet();
	BarrelParametersSet( BarrelParametersSet const & src );
	~BarrelParametersSet() override;

	ParametersSetOP clone() const override;

public: // Getters

	core::Size n_strands() const { return n_strands_; }
	core::Size shear_number() const { return shear_number_; }
	bool antiparallel() const { return antiparallel_; }
	core::Real barrel_radius() const { return barrel_radius_; }

public: // Setters

	void set_n_strands( core::Size const val ) { n_strands_ = val; }
	void set_shear_number( core::Size const val ) { shear_number_ = val; }
	void set_antiparallel( bool const val ) { antiparallel_ = val; }
	void set_barrel_radius( core::Real const val ) { barrel_radius_ = val; }

	void get_pdb_remark( std::stringstream & remark ) const override;

private:

	core::Size n_strands_;
	core::Size shear_number_;
	bool antiparallel_;
	core::Real barrel_radius_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif

}; // class BarrelParametersSet

} // namespace parameters
} // namespace beta_barrel
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_beta_barrel_parameters_BarrelParametersSet )
#endif

#endif
