// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane/AqueousPoreFinder.cc
/// @brief Compute the center, major radius and minor radius of an ellipsoidal aqueous pore
/// @author Rebecca Alford (rfalford12@gmail.com)

// Unit headers
#include <protocols/membrane/AqueousPoreFinder.hh>
#include <protocols/membrane/AqueousPoreFinderCreator.hh>

// Core headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/ImplicitLipidInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/AqueousPoreParameters.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <core/chemical/AA.hh>

// Numeric Headers
#include <numeric/MathMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/linear_algebra/minimum_bounding_ellipse.hh>
#include <numeric/linear_algebra/EllipseParameters.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>
#include <numeric/interpolation/spline/SimpleInterpolator.hh>
#include <numeric/cubic_polynomial.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ Headers
#include <ostream>
#include <fstream>

static basic::Tracer TR( "protocols.membrane.AqueousPoreFinder" );
typedef utility::vector1< numeric::CubicPolynomial > piecewise_poly;

namespace protocols {
namespace membrane {

/// @brief Default constructor
AqueousPoreFinder::AqueousPoreFinder():
	protocols::moves::Mover( AqueousPoreFinder::mover_name() ),
	tolerance_( 0.1 )
{}

/// @brief Copy constructor
AqueousPoreFinder::AqueousPoreFinder( AqueousPoreFinder const & src ):
	protocols::moves::Mover( src ),
	tolerance_( src.tolerance_ )
{}

/// @brief Destructor (important for properly forward-declaring smart-pointer members)
AqueousPoreFinder::~AqueousPoreFinder(){}

core::Real
AqueousPoreFinder::find_max_z( utility::vector1< numeric::xyzVector< core::Real > > coords ) const {
	core::Real max_z = 0;
	for ( core::Size ii = 1; ii <= coords.size(); ++ii ) {
		if ( max_z < coords[ii].z() ) {
			max_z = coords[ii].z();
		}
	}
	return max_z;
}

core::Real
AqueousPoreFinder::find_min_z( utility::vector1< numeric::xyzVector< core::Real > > coords ) const {
	core::Real min_z = 0;
	for ( core::Size ii = 1; ii <= coords.size(); ++ii ) {
		if ( min_z > coords[ii].z() ) {
			min_z = coords[ii].z();
		}
	}
	return min_z;
}

core::Real
AqueousPoreFinder::get_rotation_angle( numeric::linear_algebra::EllipseParametersOP ellipse ) const {
	core::Real a11 = ellipse->rotation()(0,0);
	return std::acos( a11 );
}

numeric::CubicPolynomial
AqueousPoreFinder::generate_cubic_polynomial_func(
	core::Real const start_coord,
	core::Real const start_value,
	core::Real const start_deriv,
	core::Real const end_coord,
	core::Real const end_value,
	core::Real const end_deriv
) const {

	using namespace numeric;
	using namespace numeric::interpolation::spline;

	SplineGenerator gen_spline( start_coord, start_value, start_deriv, end_coord, end_value, end_deriv );
	InterpolatorOP interp( gen_spline.get_interpolator() );
	SimpleInterpolatorOP sinterp = utility::pointer::dynamic_pointer_cast< SimpleInterpolator > ( interp );
	if ( ! sinterp ) {
		utility_exit_with_message( "Created non-simple-interpolator in initialization of pore parameters" );
	}

	numeric::SplineParameters sp;
	sp.ylo = sinterp->y()[1];
	sp.yhi = sinterp->y()[2];
	sp.y2lo = sinterp->ddy()[1];
	sp.y2hi = sinterp->ddy()[2];

	return cubic_polynomial_from_spline( start_coord, end_coord, sp );

}

utility::vector1< numeric::CubicPolynomial >
AqueousPoreFinder::generate_piecewise_cubic_poly_func(
	utility::vector1< core::Real > ellipse_locations,
	utility::vector1< core::Real > ellipse_parameters,
	core::Real const min_value,
	core::Real const max_value
) const {

	using namespace numeric;

	// Initialize a list of polynomials for the piecewise function
	utility::vector1< numeric::CubicPolynomial > list_of_polynomials;

	// Create the left boundary cubic polynomial func
	CubicPolynomial min_value_poly( generate_cubic_polynomial_func(
		ellipse_locations.front(), min_value, 0,
		ellipse_locations[2], ellipse_parameters.front(), 0 ));
	list_of_polynomials.push_back( min_value_poly );

	// Calculate splines that connect all of the boundary ellipse locations
	for ( core::Size ii = 2; ii < ellipse_locations.size()-1; ++ii ) {

		CubicPolynomial middle_values_poly( generate_cubic_polynomial_func(
			ellipse_locations[ii], ellipse_parameters[ii-1], 0,
			ellipse_locations[ii+1], ellipse_parameters[ii], 0 ) );
		list_of_polynomials.push_back( middle_values_poly );

	}

	// Create the right boundary cubic polynomial func
	CubicPolynomial max_value_poly( generate_cubic_polynomial_func(
		ellipse_locations[ ellipse_locations.size() - 1 ], ellipse_parameters.back(), 0,
		ellipse_locations.back(), max_value, 0 ));
	list_of_polynomials.push_back( max_value_poly );

	return list_of_polynomials;

}


utility::vector1< utility::vector1< numeric::xyzVector< core::Real > > >
AqueousPoreFinder::distribute_coords_into_bins(
	utility::vector1< numeric::xyzVector< core::Real > > aqueous_coords,
	core::Real const binsize
) const {

	using namespace numeric;

	// Calculate the minimum and maximum z-coordinate for aqueous facing atoms
	core::Real min_z = std::floor( find_min_z( aqueous_coords ) );
	core::Real max_z = std::ceil( find_max_z( aqueous_coords ) );

	// Calculate the bin size for sorting
	core::Size numbins = std::ceil( (max_z-min_z)/binsize );

	// Create a vector1 of vector1's of coordinates for data storage and return
	utility::vector1< utility::vector1< xyzVector< core::Real > > > sorted_coordinates;

	// Iterate through the aqueous coordinate files
	core::Real current_location( min_z );
	core::Real next_location( min_z + binsize );
	for ( core::Size ii = 1; ii <= numbins; ++ii ) {
		utility::vector1< xyzVector< core::Real > > coords_for_current_bin;
		for ( core::Size jj = 1; jj <= aqueous_coords.size(); ++jj ) {
			core::Real zcoord( aqueous_coords[jj].z() );
			if ( zcoord >= current_location && zcoord < next_location ) {
				coords_for_current_bin.push_back( aqueous_coords[jj] );
			}
		}
		sorted_coordinates.push_back( coords_for_current_bin );
		current_location = next_location;
		next_location = next_location + binsize;
	}

	return sorted_coordinates;
}

/// TODO: need a sub-function to construct this long polynomal
core::conformation::membrane::AqueousPoreParametersOP
AqueousPoreFinder::construct_aqueous_pore(
	core::Real const z_axis_buffer,
	core::Real const ellipse_radius_buffer,
	utility::vector1< core::Real > ellipse_locations,
	utility::vector1< numeric::linear_algebra::EllipseParametersOP > ellipse_parameters
) const {

	using namespace numeric;
	using namespace numeric::interpolation::spline;
	using namespace core::conformation::membrane;

	debug_assert( ellipse_locations.size() > 0 );
	debug_assert( ellipse_parameters.size() > 0 );

	// Generate edge coordinate parameters
	core::Real min_zcoord( ellipse_locations.front() - z_axis_buffer );
	core::Real max_zcoord( ellipse_locations.back() + z_axis_buffer );

	// Create an extended list of ellipse_locations
	utility::vector1< core::Real > extended_ellipse_locations;
	extended_ellipse_locations.push_back( min_zcoord );
	for ( core::Size ii = 1; ii <= ellipse_locations.size(); ++ii ) {
		extended_ellipse_locations.push_back( ellipse_locations[ii] );
	}
	extended_ellipse_locations.push_back( max_zcoord );

	// Extrapolate parameters from Ellipse parameters
	utility::vector1< core::Real > pore_center_x_values;
	utility::vector1< core::Real > pore_center_y_values;
	utility::vector1< core::Real > major_radius_values;
	utility::vector1< core::Real > minor_radius_values;
	utility::vector1< core::Real > rotation_angle_values;
	for ( core::Size ii = 1; ii <= ellipse_parameters.size(); ++ii ) {

		pore_center_x_values.push_back( ellipse_parameters[ii]->center_h() );
		pore_center_y_values.push_back( ellipse_parameters[ii]->center_k() );
		major_radius_values.push_back( ellipse_parameters[ii]->major_radius() );
		minor_radius_values.push_back( ellipse_parameters[ii]->minor_radius() );
		rotation_angle_values.push_back( get_rotation_angle( ellipse_parameters[ii] ) );

	}

	// Generate functions describing motion in the pore center
	core::Real min_center_x( ellipse_parameters.front()->center_h() );
	core::Real max_center_x( ellipse_parameters.back()->center_h() );
	piecewise_poly center_x_poly( generate_piecewise_cubic_poly_func(
		extended_ellipse_locations,
		pore_center_x_values,
		min_center_x,
		max_center_x
		));

	core::Real min_center_y( ellipse_parameters.front()->center_k() );
	core::Real max_center_y( ellipse_parameters.back()->center_k() );
	piecewise_poly center_y_poly( generate_piecewise_cubic_poly_func(
		extended_ellipse_locations,
		pore_center_y_values,
		min_center_y,
		max_center_y
		));

	// Generate functions describing the major and minor radius
	core::Real min_z_major_radius( ellipse_parameters.front()->major_radius() + ellipse_radius_buffer );
	core::Real max_z_major_radius( ellipse_parameters.back()->major_radius() + ellipse_radius_buffer );
	piecewise_poly major_radius_poly( generate_piecewise_cubic_poly_func(
		extended_ellipse_locations,
		major_radius_values,
		min_z_major_radius,
		max_z_major_radius
		));

	core::Real min_z_minor_radius( ellipse_parameters.front()->minor_radius() + ellipse_radius_buffer );
	core::Real max_z_minor_radius( ellipse_parameters.back()->minor_radius() + ellipse_radius_buffer );
	piecewise_poly minor_radius_poly( generate_piecewise_cubic_poly_func(
		extended_ellipse_locations,
		minor_radius_values,
		min_z_minor_radius,
		max_z_minor_radius
		));

	// Generate edge rotation angle parameters
	core::Real min_z_rotation_angle( get_rotation_angle( ellipse_parameters.front() ) );
	core::Real max_z_rotation_angle( get_rotation_angle( ellipse_parameters.back() ) );
	piecewise_poly rotation_angle_poly( generate_piecewise_cubic_poly_func(
		extended_ellipse_locations,
		rotation_angle_values,
		min_z_rotation_angle,
		max_z_rotation_angle
		));

	// Store all of the data in AqueousPoreParameters
	AqueousPoreParametersOP aqueous_pore( new AqueousPoreParameters(
		min_center_x, max_center_x, min_center_y, max_center_y, min_z_major_radius, max_z_major_radius,
		min_z_minor_radius, max_z_minor_radius, min_z_rotation_angle, max_z_rotation_angle,
		extended_ellipse_locations, center_x_poly, center_y_poly, major_radius_poly, minor_radius_poly, rotation_angle_poly ));

	return aqueous_pore;
}

/// @brief Apply the mover
void
AqueousPoreFinder::apply( core::pose::Pose & pose ) {

	using namespace utility;
	using namespace numeric;
	using namespace numeric::linear_algebra;
	using namespace numeric::interpolation::spline;
	using namespace core::conformation::membrane;

	// Incclude a case where there are less than 3 transmembrane spans
	if ( pose.conformation().membrane_info()->spanning_topology()->nspans() >= 3 ) {

		TR << "Initializing a custom aqueous pore boundary" << std::endl;

		// Detect aqueous-facing residues in the membrane - collect the xyz coordinates
		utility::vector1< xyzVector< core::Real > > aqueous_coords;
		for ( core::Size ii = 1; ii <= nres_protein( pose ); ++ii ) {

			// Add CA coordinate
			core::Size CA( 2 );
			core::Real ca_acc( pose.conformation().membrane_info()->implicit_lipids()->per_atom_lipid_accessibility(ii,CA) );
			bool ca_in_memb( pose.conformation().membrane_info()->in_membrane( pose.conformation(), ii ) );
			if ( ca_acc == 0.0 && ca_in_memb ) {
				aqueous_coords.push_back( pose.residue( ii ).atom( "CA" ).xyz() );
			}

			// Add CB coordinate
			if ( pose.residue(ii).aa() != core::chemical::AA::aa_gly ) {
				core::Size CB( pose.residue(ii).first_sidechain_atom() );
				core::Real cb_acc( pose.conformation().membrane_info()->implicit_lipids()->per_atom_lipid_accessibility(ii,CB) );
				bool cb_in_memb( pose.conformation().membrane_info()->in_membrane( pose.conformation(), ii ) );
				if ( cb_acc == 0.0 && cb_in_memb ) {
					aqueous_coords.push_back( pose.residue( ii ).atom( "CB" ).xyz() );
				}
			}
		}

		// Sort coordinates into bins
		core::Real binsize(5);
		utility::vector1< utility::vector1< numeric::xyzVector< core::Real > > > sorted_coordinates = distribute_coords_into_bins( aqueous_coords, binsize );

		// Calculate an ellipse for the coordinates in each bin
		utility::vector1< EllipseParametersOP > ellipse_parameters;
		for ( core::Size ii = 1; ii <= sorted_coordinates.size(); ++ii ) {

			// Initialize a blank ellipse
			EllipseParametersOP current_ellipse = EllipseParametersOP( new EllipseParameters() );

			// If there is more than 10 points, build a minimum bounding ellipse
			if ( sorted_coordinates[ii].size() >= 10 ) {

				EllipseParametersOP ellipse = minimum_bounding_ellipse( sorted_coordinates[ii], tolerance_ );

				// Widen the ellipse for alpha helical proteins
				if ( pose.conformation().membrane_info()->implicit_lipids()->is_helical() ) {
					ellipse->add_buffer( 3, 3 );
				}

				ellipse_parameters.push_back( ellipse );
				current_ellipse = ellipse;

			} else {

				ellipse_parameters.push_back( current_ellipse );

			}

		}

		// Calculate the minimum and maximum z-coordinate for aqueous facing atoms
		core::Real min_z = std::floor( find_min_z( aqueous_coords ) );
		core::Real max_z = std::ceil( find_max_z( aqueous_coords ) );

		// Calculate the bin size for sorting
		core::Size numbins = std::ceil( (max_z-min_z)/binsize );

		// Make a list of ellipse locations
		utility::vector1< core::Real > ellipse_locations;
		core::Real current_location = min_z;
		for ( core::Size ii = 1; ii <= numbins; ++ii ) {
			ellipse_locations.push_back( current_location );
			current_location = current_location + binsize;
		}

		// Once you have your array of ellipses, write the parameters to a file
		TR << "Writing skeleton ellipse layers to an output file" << std::endl;
		utility::vector1< std::string > temp( utility::string_split( pose.pdb_info()->name(), '/') );
		std::string tempstr = temp[ temp.size() ].substr(0, temp[ temp.size() ].size()-4 );
		std::string filename( tempstr + "_" + pose.conformation().membrane_info()->implicit_lipids()->lipid_composition_name() + "_ellipse.dat" );
		utility::io::ozstream output( filename );

		// Write header
		output << "ellipse_no z_loc h k a b r00 r01 r10 r11" << std::endl;

		for ( core::Size ii = 1; ii <= ellipse_parameters.size(); ++ii ) {
			output << ii << " " << ellipse_locations[ii] << " " << ellipse_parameters[ii]->center_h() << " " << ellipse_parameters[ii]->center_k() << " " << ellipse_parameters[ii]->major_radius() << " " << ellipse_parameters[ii]->minor_radius();
			output << " " << ellipse_parameters[ii]->rotation()(0,0);
			output << " " << ellipse_parameters[ii]->rotation()(0,1);
			output << " " << ellipse_parameters[ii]->rotation()(1,0);
			output << " " << ellipse_parameters[ii]->rotation()(1,1) << std::endl;
		}

		// Construct piecewise cubic polynomials to approximate the changing elliptical pore shape
		core::Real z_axis_buffer(binsize/2);
		core::Real ellipse_radius_buffer(5);

		// Make an object w/ AqueousPore information
		AqueousPoreParametersOP aqueous_pore(
			construct_aqueous_pore( z_axis_buffer, ellipse_radius_buffer, ellipse_locations, ellipse_parameters ) );

		// Here is the step where AqueousPoreParameters will be set within ImplicitLipidInfo
		TR << "Setting aqueous pore parameters in implicit lipid info" << std::endl;
		pose.conformation().membrane_info()->implicit_lipids()->set_aqueous_pore_parameters( aqueous_pore  );
		pose.conformation().membrane_info()->implicit_lipids()->has_pore( true );

	} else {

		TR << "Pose has too few TM spans for a pore. Setting an empty pore parameter set in implicit lipid info" << std::endl;

		// Approximate the minimum bounding ellipse as an extremely small center point
		utility::vector1< core::Real > empty_boundaries;
		utility::vector1< CubicPolynomial > empty_poly;
		AqueousPoreParametersOP aqueous_pore_zero( new AqueousPoreParameters(
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			empty_boundaries, empty_poly, empty_poly, empty_poly, empty_poly, empty_poly ));

		// Set parameters
		pose.conformation().membrane_info()->implicit_lipids()->set_aqueous_pore_parameters( aqueous_pore_zero  );
		pose.conformation().membrane_info()->implicit_lipids()->has_pore( false );
	}
}


/// @brief: Return a rotation matrix that results in no rotation
numeric::MathMatrix< core::Real >
AqueousPoreFinder::get_dummy_rotation_matrix() const {

	using namespace numeric;
	MathMatrix< core::Real > rot_matrix( 2, 2 );
	rot_matrix(0,0) = 1;
	rot_matrix(0,1) = 0;
	rot_matrix(1,0) = 0;
	rot_matrix(1,1) = 1;
	return rot_matrix;

}

/// @brief Show the contents of the Mover
void
AqueousPoreFinder::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
AqueousPoreFinder::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	if ( tag->hasOption( "tolerance" ) ) {
		tolerance_ = tag->getOption< core::Real >( "tolerance" );
	}
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
AqueousPoreFinder::fresh_instance() const
{
	return protocols::moves::MoverOP( new AqueousPoreFinder );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
AqueousPoreFinder::clone() const
{
	return protocols::moves::MoverOP( new AqueousPoreFinder( *this ) );
}

std::string AqueousPoreFinder::get_name() const {
	return mover_name();
}

std::string AqueousPoreFinder::mover_name() {
	return "AqueousPoreFinder";
}

void AqueousPoreFinder::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "tolerance", xsct_real, "Tolerance for minimum bounding ellipse estimation" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Mover to calculate the parameters of an elliptical membrane pore", attlist );
}

// Creators
protocols::moves::MoverOP
AqueousPoreFinderCreator::create_mover() const
{
	return protocols::moves::MoverOP( new AqueousPoreFinder );
}

std::string
AqueousPoreFinderCreator::keyname() const
{
	return AqueousPoreFinder::mover_name();
}

void AqueousPoreFinderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AqueousPoreFinder::provide_xml_schema( xsd );
}



std::ostream &
operator<<( std::ostream & os, AqueousPoreFinder const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //membrane
