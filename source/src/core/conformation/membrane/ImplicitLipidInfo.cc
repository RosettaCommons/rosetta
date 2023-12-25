// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/membrane/ImplicitLipidInfo.hh
/// @brief An Implicit Lipid Membrane Model
/// @version menv-franklin2018
/// @cite Alford et al. (2019) "A thermodynamically-inspired, biologically-driven energy function for membrane protein modeling and design"
///
/// @detail This class defines physical and chemical properties of the membrane environment. It
/// is mainly used by the membrane energy function.
///  1. Parameters of the hydration function (e.g. thickness, rate, of transition, pore size)
///  2. Lipid composition details
///  3. Hydration function smoothing parameters
///  4. Structure-based lipid accessibiity information
///
/// @author Rebecca Alford (rfalford12@gmail.com)

#include <core/conformation/membrane/ImplicitLipidInfo.hh>

// Pakcage Headers
#include <core/conformation/membrane/AqueousPoreParameters.hh>

#include <core/types.hh>


#include <numeric/MathMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/cubic_polynomial.hh>

// Utility Headers
#include <basic/database/open.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

#include <cmath>

static basic::Tracer TR( "core.conformation.membrane.ImplicitLipidInfo" );
typedef utility::vector1< numeric::CubicPolynomial > piecewise_poly;

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace conformation {
namespace membrane {

// privatized default constructor
ImplicitLipidInfo::ImplicitLipidInfo() :
	utility::VirtualBase(),
	water_thickness_( 0.0 ),
	change_in_water_density_( 0.0 ),
	transformed_water_thickness_( 0.0 ),
	chain_type_( "" ),
	headgroup_type_( "" ),
	lipid_composition_name_( "" ),
	lipid_composition_name_long_( "" ),
	degrees_of_saturation_( 0 ),
	temperature_( 0.0 ),
	is_helical_( true ),
	pore_params_()
{}

ImplicitLipidInfo::ImplicitLipidInfo(
	std::string lipid_composition_name,
	core::Real temperature
) : utility::VirtualBase(),
	water_thickness_( 0.0 ),
	change_in_water_density_( 0.0 ),
	transformed_water_thickness_( 0.0 ),
	chain_type_( "" ),
	headgroup_type_( "" ),
	lipid_composition_name_( lipid_composition_name ),
	lipid_composition_name_long_( "" ),
	degrees_of_saturation_( 0 ),
	temperature_( temperature ),
	is_helical_( true ),
	pore_params_()
{

	// Initialize lipid composition specific parameters
	initialize_implicit_lipid_parameters();
	initialize_implicit_lipid_electricfield_parameters();
}

ImplicitLipidInfo::~ImplicitLipidInfo() {}

ImplicitLipidInfoOP
ImplicitLipidInfo::clone() const {
	return ImplicitLipidInfoOP( new ImplicitLipidInfo( *this ) );
}

void
ImplicitLipidInfo::show( std::ostream & output ) const {

	// Show information specific to Lipid Membrane Info
	output << "Information about the Lipid Composition: " << std::endl;
	output << "Water thickness: " << water_thickness_ << std::endl;
	output << "Change in water density: " << change_in_water_density_ << std::endl;
	output << "Transformed water thickness: " << transformed_water_thickness_ << std::endl;
	output << "Lipid chain type: " << chain_type_ << std::endl;
	output << "Lipid headgroup type: " << headgroup_type_ << std::endl;
	output << "Lipid composition name (short-form): " << lipid_composition_name_ << std::endl;
	output << "Lipid composition name (long-form): " << lipid_composition_name_long_ << std::endl;
	output << "Degrees of saturation: " << degrees_of_saturation_ << std::endl;
	output << "Temperature (celcius): " << temperature_ << std::endl;
	output << "Is this an alpha helical protein?: ";

	if ( is_helical_ ) {
		output << "yes" << std::endl;
	} else {
		output << "no" << std::endl;
	}

}

// Getters for information about the lipid composition

/// @brief Number of carbons and degrees of saturation in the chains
std::string
ImplicitLipidInfo::chain_type() const {
	return chain_type_;
}

/// @brief Chemical name of the headgroup
std::string
ImplicitLipidInfo::headgroup_type() const {
	return headgroup_type_;
}

/// @brief Abbreviated name of the lipid composition
std::string
ImplicitLipidInfo::lipid_composition_name() const {
	return lipid_composition_name_;
}

/// @brief Full name of the lipid composiiton
std::string
ImplicitLipidInfo::lipid_composition_name_long() const {
	return lipid_composition_name_long_;
}

/// @brief Number of degrees of saturation in the overall lipid
core::Real
ImplicitLipidInfo::degrees_of_saturation() const {
	return degrees_of_saturation_;
}

/// @brief Temperature at which the Db parameter was measured/calculated (in celcius)
core::Real
ImplicitLipidInfo::temperature() const {
	return temperature_;
}


/// @brief Is the protein alpha helical or beta barrel
bool
ImplicitLipidInfo::is_helical() const {
	return is_helical_;
}

void
ImplicitLipidInfo::is_helical( bool const is_helical ) {
	is_helical_ = is_helical;
}


/// @brief Make empty pore parameters
void
ImplicitLipidInfo::make_no_pore_parameters() {

	// Make a zero pore AqueousPoreParams object
	using namespace numeric;
	utility::vector1< core::Real > empty_boundaries;
	utility::vector1< CubicPolynomial > empty_poly;
	AqueousPoreParametersOP aqueous_pore_zero( new AqueousPoreParameters(
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		empty_boundaries, empty_poly, empty_poly, empty_poly, empty_poly, empty_poly ));
	pore_params_ = aqueous_pore_zero;
}

/// @brief Membrane Aqueous pore center - h parameter
core::Real
ImplicitLipidInfo::pore_center_x( core::Real const zcoord ) const {
	return pore_params_->pore_center_x( zcoord );
}

/// @brief Membrane aqueous pore center - k parameter
core::Real
ImplicitLipidInfo::pore_center_y( core::Real const zcoord ) const {
	return pore_params_->pore_center_y( zcoord );
}

/// @brief Membrane aqueous pore - major radius
core::Real
ImplicitLipidInfo::pore_major_radius( core::Real const zcoord ) const {
	return pore_params_->pore_major_radius( zcoord );
}

/// @brief Membrane aqueous pore - minor radius
core::Real
ImplicitLipidInfo::pore_minor_radius( core::Real const zcoord ) const {
	return pore_params_->pore_minor_radius( zcoord );
}

/// @brief Membrane aqueous pore - rotation matrix
numeric::MathMatrix< core::Real >
ImplicitLipidInfo::pore_rotation( core::Real const zcoord ) const {
	return pore_params_->pore_rotation( zcoord );
}


/// @brief Set membrane aqueous pore parameters
void
ImplicitLipidInfo::set_aqueous_pore_parameters(
	AqueousPoreParametersOP aqueous_pore
) {
	pore_params_ = aqueous_pore;
}


// Chemical information about this membrane
/// @brief Water thickness of the membrane
/// @detail Thickness of the membrane, defined by the pair of z coordinates
/// with a water density of 50% (SAXS Db value)
core::Real
ImplicitLipidInfo::water_thickness() const {
	return water_thickness_;
}

/// @brief Change in water density from membrane core to water bulk water
/// @detail Steepness defined by the number of waters lost per increase in z
/// from the membrane center (s = b)
core::Real
ImplicitLipidInfo::water_steepness() const {
	return change_in_water_density_;
}

/// @brief Pseudo-thickness parameter
/// @details A parameter in the expotential membrane definition
/// t = -(1/b) * ln(1/A)
core::Real
ImplicitLipidInfo::water_pseudo_thickness() const {
	return transformed_water_thickness_;
}

/// @brief Overall hydration given the atomic depth and cavity structure
core::Real
ImplicitLipidInfo::f_elec_field( core::Real z ) const {
	core::Real phi_z( 0.0 );

	if ( z>=-1.0*std::abs(center_root_) && z<=std::abs(center_root_) ) {
		///phi_z = center_a_ * std::exp(-1.0*(z-center_b_)*(z-center_b_)/(center_c_*center_c_));
		phi_z = ( center_a_ * z * z * z * z ) + ( center_b_ * z * z * z ) + ( center_c_ * z * z ) + ( center_d_ * z ) + center_e_ ;

	} else {
		phi_z = pp_a_ + (( pp_a_-pp_b_ )/( 1.0 + std::exp((std::abs(z)-pp_c_)/pp_d_) ));
	}
	return phi_z;
}

/// @brief Gradient
core::Real
ImplicitLipidInfo::f_elec_field_gradient( core::Real z ) const {

	core::Real grad_phi_z( 0.0 );
	//bringing the derivative as close as possible. 0.70 came with trial and error/
	if ( z>=-1.0*std::abs(center_root_+0.70) && z<=std::abs(center_root_+0.70) ) {
		/// grad_phi_z = -2.0 * (z-center_b_) * (center_a_/(center_c_*center_c_)) * std::exp(-1.0*(z-center_b_)*(z-center_b_)/(center_c_*center_c_));
		grad_phi_z = ( 4.0 * center_a_ * z * z * z ) + ( 3.0 * center_b_ * z * z ) + ( 2.0 * center_c_ * z ) + center_d_ ;

	} else {
		grad_phi_z = -1.0 * (z/std::abs(z)) * (( pp_a_-pp_b_ )/pp_d_)*(std::exp((std::abs(z)-pp_c_)/pp_d_)/( (1.0 + std::exp((std::abs(z)-pp_c_)/pp_d_) )*(1.0 + std::exp((std::abs(z)-pp_c_)/pp_d_) )));
	}

	return grad_phi_z;
}


// Private helper functions for initializing polynomials and parameters

/// @brief Helper function to initialize the lipid composiiton data
void
ImplicitLipidInfo::initialize_implicit_lipid_parameters() {

	using namespace basic;
	std::string dbfile( "membrane/implicit_lipid_parameters.txt" );
	utility::io::izstream infile;

	database::open( infile, dbfile );
	if ( !infile.good() ) {
		utility_exit_with_message( "Unable to open database file containing implicit lipid parameters" );
	}

	std::string lipid_composition_name(""), lipid_composition_name_long(""), headgroup_type(""), chain_type("");
	std::string water_thickness(""), transformed_water_thickness(""), change_in_water_density(""), temperature(""), degrees_of_saturation("");

	bool lipid_composition_found(false);
	bool temperature_found(false);

	std::string line;
	getline( infile, line ); // throw out header line
	while ( getline( infile, line ) ) {
		utility::trim(line, "\t\n" ); // remove leading and trailing spaces
		if ( line.empty() > 0 ) continue;
		if ( line.find("#",0) == 0 ) continue; // skip comment lines
		if ( line.find("LIPID_NAME",0) == 0 ) continue; // slip header line

		std::istringstream l( line );
		l >> lipid_composition_name >> lipid_composition_name_long >> headgroup_type >> chain_type >> water_thickness >> transformed_water_thickness >> change_in_water_density >> temperature >> degrees_of_saturation;
		if ( l.fail() ) {
			utility_exit_with_message( "bad line: " + line );
		}

		if ( lipid_composition_name == lipid_composition_name_ ) {
			lipid_composition_found = true;
			core::Real temperature2 = utility::from_string( temperature, core::Real(0.0) );
			if ( temperature2 == temperature_ ) {
				temperature_found = true;
				lipid_composition_name_long_ = lipid_composition_name_long;
				headgroup_type_ = headgroup_type;
				chain_type_ = chain_type;
				water_thickness_ = utility::from_string( water_thickness, core::Real(0.0) );
				transformed_water_thickness_ = utility::from_string( transformed_water_thickness, core::Real(0.0) );
				change_in_water_density_ = utility::from_string( change_in_water_density, core::Real(0.0) );
				temperature_ = utility::from_string( temperature, core::Real(0.0) );
				degrees_of_saturation_ = utility::from_string( degrees_of_saturation, core::Real(0.0) );
			}
		}
	}

	// Check that a lipid composition was found
	if ( !lipid_composition_found ) {
		std::string msg = "Lipid composition parameters for " + lipid_composition_name_ + " were not found!";
		utility_exit_with_message( msg );
	} else if ( !temperature_found ) {
		std::string msg = "Parameters for lipid composition " + lipid_composition_name_ + " at temp " + std::to_string( temperature_ ) + " were not found";
		utility_exit_with_message( msg);
	}
}

// Private helper functions for initializing polynomials and parameters

/// @brief Helper function to initialize the lipid composiiton data
void
ImplicitLipidInfo::initialize_implicit_lipid_electricfield_parameters() {

	std::string dbfile( "membrane/lipid_electric_field_params.txt" );
	using namespace basic;
	using namespace core;

	bool lipid_composition_found(false);
	utility::io::izstream infile;
	TR << "Reading electric field fitting parameters from the database" << std::endl;
	database::open( infile, dbfile );
	if ( !infile.good() ) {
		utility_exit_with_message( "Unable to open database file containing electric field fitting parameters" );
	}

	std::string line;
	getline( infile, line );
	while ( getline( infile, line ) ) {
		utility::trim(line, "\t\n" );
		if ( line.empty() > 0 ) continue;
		if ( line.find("#",0) == 0 ) continue;

		std::istringstream l( line );

		std::string lipid_type(""), center_a(""), center_b("");
		std::string center_c(""), center_d(""), center_e(""), center_root("");
		std::string pp_a(""), pp_b(""), pp_c(""), pp_d("");
		l >> lipid_type >> center_a >> center_b >> center_c >> center_d >> center_e >> center_root >> pp_a >> pp_b >> pp_c >> pp_d ;
		//TR << "Reading the line is done" <<std::endl;
		if ( l.fail() ) {
			utility_exit_with_message( "bad input line: " + line );
		}

		if ( lipid_type == lipid_composition_name_ ) {
			lipid_composition_found = true;
			//should there be more meaningful names? will have to think more.
			if ( center_a == " " || center_b == " " || center_c == " " || center_d == " " || center_e == " " || center_root == " " ) {
				std::string msg = "center_variables for " + lipid_composition_name_ + " is missing!";
				utility_exit_with_message( msg );
			} else {


				center_a_ = utility::from_string( center_a, core::Real(0.0) );
				center_b_ = utility::from_string( center_b, core::Real(0.0) );
				center_c_ = utility::from_string( center_c, core::Real(0.0) );
				center_d_ = utility::from_string( center_d, core::Real(0.0) );
				center_e_ = utility::from_string( center_e, core::Real(0.0) );
				center_root_ = utility::from_string( center_root, core::Real(0.0) );

			}

			//should there be more meaningful names? will have to think more.
			if ( pp_a == " " || pp_b == " " || pp_c == " " || pp_d == " " ) {
				std::string msg = "pp_variables for " + lipid_composition_name_ + " is missing!";
				utility_exit_with_message( msg );
			} else {

				TR << "lipid parameters are read" <<std::endl;
				//TR << "pp_a:"<< pp_a <<"pp_b:"<< pp_b <<"pp_c:"<< pp_c <<std::endl;

				//pp_a_ = std::stof( pp_a );
				//pp_b_ = std::stof( pp_b );
				//pp_c_ = std::stof( pp_c );
				//pp_d_ = std::stof( pp_d );

				pp_a_ = utility::from_string( pp_a, core::Real(0.0) );
				pp_b_ = utility::from_string( pp_b, core::Real(0.0) );
				pp_c_ = utility::from_string( pp_c, core::Real(0.0) );
				pp_d_ = utility::from_string( pp_d, core::Real(0.0) );

				// TR << "pp_a:"<< pp_a_ <<"pp_b:"<< pp_b_ <<"pp_c:"<< pp_c_ <<std::endl;


			}
		}

	}
	// Check that a lipid composition was found
	if ( !lipid_composition_found ) {
		std::string msg = "Lipid composition parameters for " + lipid_composition_name_ + " were not found!";
		utility_exit_with_message( msg );
	}
}


} //core
} //conformation
} //membrane


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::membrane::ImplicitLipidInfo::save( Archive & arc ) const {
	arc( CEREAL_NVP( water_thickness_ ) ); // core::Real
	arc( CEREAL_NVP( change_in_water_density_ ) ); // core::Real
	arc( CEREAL_NVP( transformed_water_thickness_ ) ); // core::Real
	arc( CEREAL_NVP( chain_type_ ) ); // std::string
	arc( CEREAL_NVP( headgroup_type_ ) ); // std::string
	arc( CEREAL_NVP( lipid_composition_name_ ) ); // std::string
	arc( CEREAL_NVP( lipid_composition_name_long_ ) ); // std::string
	arc( CEREAL_NVP( degrees_of_saturation_ ) ); // core::Real
	arc( CEREAL_NVP( temperature_ ) ); // core::Real
	arc( CEREAL_NVP( is_helical_ ) ); // _Bool
	arc( CEREAL_NVP( pore_params_ ) ); // AqueousPoreParametersOP

	arc( CEREAL_NVP( center_a_ ) ); //core::Real
	arc( CEREAL_NVP( center_b_ ) ); //core::Real
	arc( CEREAL_NVP( center_c_ ) ); //core::Real
	arc( CEREAL_NVP( center_d_ ) ); //core::Real
	arc( CEREAL_NVP( center_e_ ) ); //core::Real
	arc( CEREAL_NVP( center_root_ ) ); //core::Real
	arc( CEREAL_NVP( pp_a_ ) ); //core::Real
	arc( CEREAL_NVP( pp_b_ ) ); //core::Real
	arc( CEREAL_NVP( pp_c_ ) ); //core::Real
	arc( CEREAL_NVP( pp_d_ ) ); //core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::membrane::ImplicitLipidInfo::load( Archive & arc ) {
	arc( water_thickness_ ); // core::Real
	arc( change_in_water_density_ ); // core::Real
	arc( transformed_water_thickness_ ); // core::Real
	arc( chain_type_ ); // std::string
	arc( headgroup_type_ ); // std::string
	arc( lipid_composition_name_ ); // std::string
	arc( lipid_composition_name_long_ ); // std::string
	arc( degrees_of_saturation_ ); // core::Real
	arc( temperature_ ); // core::Real
	arc( is_helical_ ); // _Bool
	arc( pore_params_ ); // AqueousPoreParametersOP

	arc( center_a_ ); //core::Real
	arc( center_b_ ); //core::Real
	arc( center_c_ ); //core::Real
	arc( center_d_ ); //core::Real
	arc( center_e_ ); //core::Real
	arc( center_root_ ); //core::Real
	arc( pp_a_ ); //core::Real
	arc( pp_b_ ); //core::Real
	arc( pp_c_ ); //core::Real
	arc( pp_d_ ); //core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::membrane::ImplicitLipidInfo );
CEREAL_REGISTER_TYPE( core::conformation::membrane::ImplicitLipidInfo )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_membrane_ImplicitLipidInfo )
#endif // SERIALIZATION
