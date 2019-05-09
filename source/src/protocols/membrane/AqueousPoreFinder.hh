// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane/AqueousPoreFinder.hh
/// @brief Compute the center, major axis and minor axis of an ellipsoidal aqueous pore
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_AqueousPoreFinder_HH
#define INCLUDED_protocols_membrane_AqueousPoreFinder_HH

// Unit headers
#include <protocols/membrane/AqueousPoreFinder.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/conformation/membrane/AqueousPoreParameters.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <numeric/MathMatrix.hh>
#include <numeric/linear_algebra/EllipseParameters.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <numeric/cubic_polynomial.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace membrane {

///@brief Compute the center, major axis and minor axis of an ellipsoidal aqueous pore
class AqueousPoreFinder : public protocols::moves::Mover {

	typedef utility::vector1< numeric::CubicPolynomial > piecewise_poly;

public:

	/// @brief Default constructor
	AqueousPoreFinder();

	/// @brief Copy constructor (not needed unless you need deep copies)
	AqueousPoreFinder( AqueousPoreFinder const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~AqueousPoreFinder() override;

public:

	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	/// @brief Return a rotation matrix that results in no rotation
	numeric::MathMatrix< core::Real > get_dummy_rotation_matrix() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private: // methods

	core::Real
	find_max_z( utility::vector1< numeric::xyzVector< core::Real > > coords ) const;

	core::Real
	find_min_z( utility::vector1< numeric::xyzVector< core::Real > > coords ) const;

	core::Real
	get_rotation_angle( numeric::linear_algebra::EllipseParametersOP ellipse ) const;

	utility::vector1< numeric::CubicPolynomial >
	generate_piecewise_cubic_poly_func(
		utility::vector1< core::Real > ellipse_locations,
		utility::vector1< core::Real > ellipse_parameters,
		core::Real const min_value,
		core::Real const max_value
	) const;

	numeric::CubicPolynomial
	generate_cubic_polynomial_func(   core::Real const start_coord,
		core::Real const start_value,
		core::Real const start_deriv,
		core::Real const end_coord,
		core::Real const end_value,
		core::Real const end_deriv
	) const;

	utility::vector1< utility::vector1< numeric::xyzVector< core::Real > > >
	distribute_coords_into_bins(
		utility::vector1< numeric::xyzVector< core::Real > > aqueous_coords,
		core::Real const binsize
	) const;

	core::conformation::membrane::AqueousPoreParametersOP
	construct_aqueous_pore(
		core::Real const z_axis_buffer,
		core::Real const ellipse_radius_buffer,
		utility::vector1< core::Real > ellipse_locations,
		utility::vector1< numeric::linear_algebra::EllipseParametersOP > ellipse_parameters
	) const;

private: // data

	core::Real tolerance_;

};

std::ostream &
operator<<( std::ostream & os, AqueousPoreFinder const & mover );

} //protocols
} //membrane

#endif //protocols_membrane_AqueousPoreFinder_HH
