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
	pore_params_(),
	pore_transition_steepness_( 0.0 )
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
	pore_params_(),
	pore_transition_steepness_( 10.0 )
{

	// Initialize lipid composition specific parameters
	initialize_implicit_lipid_parameters();
}

ImplicitLipidInfo::ImplicitLipidInfo( ImplicitLipidInfo const & src ) :
	utility::VirtualBase( src ),
	water_thickness_( src.water_thickness_ ),
	change_in_water_density_( src.change_in_water_density_ ),
	transformed_water_thickness_( src.transformed_water_thickness_ ),
	chain_type_( src.chain_type_ ),
	headgroup_type_( src.headgroup_type_ ),
	lipid_composition_name_( src.lipid_composition_name_ ),
	lipid_composition_name_long_( src.lipid_composition_name_long_ ),
	degrees_of_saturation_( src.degrees_of_saturation_ ),
	temperature_( src.temperature_ ),
	is_helical_( src.is_helical_ ),
	pore_transition_steepness_( src.pore_transition_steepness_ )
{}

ImplicitLipidInfo &
ImplicitLipidInfo::operator=( ImplicitLipidInfo const & src ) {

	// Abort self-assignment
	if ( this == &src ) {
		return *this;
	}

	// make a deep copy of everything
	this->water_thickness_ = src.water_thickness_;
	this->change_in_water_density_ = src.change_in_water_density_;
	this->transformed_water_thickness_ = src.transformed_water_thickness_;
	this->chain_type_ = src.chain_type_;
	this->headgroup_type_ = src.headgroup_type_;
	this->lipid_composition_name_ = src.lipid_composition_name_;
	this->lipid_composition_name_long_ = src.lipid_composition_name_long_;
	this->degrees_of_saturation_ = src.degrees_of_saturation_;
	this->is_helical_ = src.is_helical_;
	this->temperature_ = src.temperature_;
	this->pore_transition_steepness_ = src.pore_transition_steepness_;

	return *this;
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
	arc( CEREAL_NVP( pore_transition_steepness_ ) ); // core::Real
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
	arc( pore_transition_steepness_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::membrane::ImplicitLipidInfo );
CEREAL_REGISTER_TYPE( core::conformation::membrane::ImplicitLipidInfo )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_membrane_ImplicitLipidInfo )
#endif // SERIALIZATION
