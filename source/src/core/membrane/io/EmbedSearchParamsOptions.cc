// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedSearchParamsOptions.cc
///
/// @brief      Embedding Search Parameters Options class  - Contains options for membrane search and score
/// @details    Membrane proteins in rosetta use the membrane scoring function and an MCM embedidng search
///             which can be tuned and adjusted using the following options.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedSearchParamsOptions_cc
#define INCLUDED_core_membrane_io_EmbedSearchParamsOptions_cc

// Unit Headers
#include <core/membrane/io/EmbedSearchParamsOptions.hh>
#include <core/membrane/io/EmbedSearchParamsOptionsCreator.hh>

// Project Headers
#include <core/conformation/membrane/Exceptions.hh>

// Utility Headers
#include <utility/tag/Tag.hh>

// C++ Headers
#include <cstdlib>

using namespace core::conformation::membrane;

namespace core {
namespace membrane {
namespace io {

/// @brief Constructor
EmbedSearchParamsOptions::EmbedSearchParamsOptions() : basic::resource_manager::ResourceOptions(),
		normal_search_(false),
		normal_delta_angle_(0),
		normal_start_angle_(0),
		normal_max_angle_(0),
		center_search_(false),
		center_max_delta_(0),
		center_mag_(0),
		normal_mag_(0),
		normal_cycles_(0),
		penalties_(false),
		no_interpolate_mpair_(false)
{}

/// @brief Destructor
EmbedSearchParamsOptions::~EmbedSearchParamsOptions() {}

/// @brief Parse .xml file for embedding opts
void
EmbedSearchParamsOptions::parse_my_tag(
							  utility::tag::TagCOP tag
							  )
{

	// Read in search parameters for the normal
	set_normal_search( tag->getOption< bool >( "normal_search", true ));
	set_normal_start_angle( tag->getOption< core::Real >( "normal_start_angle", true ));
	set_normal_max_angle( tag->getOption< core::Real >( "normal_max_angle", true ));
	set_normal_delta_angle( tag->getOption< core::Real >( "normal_delta_angle", true ));

	// Read in search parameters for the center
	set_center_search( tag->getOption< bool >( "center_search", true ));
	set_center_max_delta( tag->getOption< core::Real >( "center_max_delta", true ));

	// Read in normal/center mag
	set_normal_mag( tag->getOption< core::Size >( "normal_mag", true ));
	set_center_mag( tag->getOption< core::Size >( "center_mag", true ));

	// Read in normal cycles
	set_normal_cycles( tag->getOption< core::Size >( "normal_cycles", true ));

	// Read in new params
	set_penalties( tag->getOption< bool >( "penalties", true) );
	set_no_interpolate_mpair( tag->getOption< bool >( "no_interpolate_Mpair", true));

}

/// @brief Return Options Class Type - Embedding Search Parameters
std::string
EmbedSearchParamsOptions::type() const
{
	return "EmbedSearchParamsOptions";
}


// Setters and Getters

/// @brief Getter/Setter Pair for normal search
/// @details Param Type: Boolean
bool EmbedSearchParamsOptions::normal_search() const { return normal_search_; }
void EmbedSearchParamsOptions::set_normal_search( bool setting ) { normal_search_ = setting; }

/// @brief Getter/Setter Pair for Normal Start Angle
/// @details Param Type: abs(int)
/// Preconditions: 0 <= setting && setting <= 360
core::Real EmbedSearchParamsOptions::normal_start_angle() const { return normal_start_angle_; }
void EmbedSearchParamsOptions::set_normal_start_angle( core::Real setting )
{
	if ( setting < 0 || setting > 360 ) {
		throw EXCN_Illegal_Arguments("Invlaid argument for normal start angle. Parameter must be >= 0");
	}
	normal_start_angle_ = setting;
}

/// @brief Getter/Setter Pair for Normal Max Angle
/// @details Param Type: abs(int)
/// Precodnition  0 <= setting && setting <= 360
core::Real EmbedSearchParamsOptions::normal_max_angle() const { return normal_max_angle_; }
void EmbedSearchParamsOptions::set_normal_max_angle( core::Real setting ) {

	if ( setting < 0 || setting > 360 ) {
		throw EXCN_Illegal_Arguments("Invlaid argument for normal max angle. Parameter must be >= 0");
	}
	normal_max_angle_ = setting;
}

/// @brief Getter/Setter Pair for Normal delta Angle
/// @details Param Type: abs(int)
/// Precodnition  0 <= setting && setting <= 360
core::Real EmbedSearchParamsOptions::normal_delta_angle() const { return normal_delta_angle_; }
void EmbedSearchParamsOptions::set_normal_delta_angle( core::Real setting ) {

	if ( setting < 0 || setting > 360 ) {
		throw EXCN_Illegal_Arguments("Invlaid argument for normal delta angle. Parameter must be >= 0");
	}

	normal_delta_angle_ = setting;
}

/// @brief Getter/Setter Pair for center search
/// @details Param Type: Boolean
bool EmbedSearchParamsOptions::center_search() const { return center_search_; }
void EmbedSearchParamsOptions::set_center_search( bool setting ) { center_search_ = setting; }

/// @brief Getter/Setter Pair for Center max delta
/// @details Param Type: abs(int)
/// Precodnition  0 <= setting
core::Real EmbedSearchParamsOptions::center_max_delta() const { return center_max_delta_; }
void EmbedSearchParamsOptions::set_center_max_delta( core::Real setting ) {

	if ( setting < 0 ) {
		throw EXCN_Illegal_Arguments("Invlaid argument for center max delta. Parameter must be >= 0");
	}

	center_max_delta_ = setting;
}

/// @brief Getter/Setter Pair for Center Mag
/// @details Param Type: abs(int)
/// Precodnition setting >=0
core::Size EmbedSearchParamsOptions::center_mag() const { return center_mag_; }
void EmbedSearchParamsOptions::set_center_mag( core::Size setting ) {
	center_mag_ = setting;
}

/// @brief Getter/Setter Pair for Normal Mag
/// @details Param Type: abs(int)
/// Precodnition setting >=0
core::Size EmbedSearchParamsOptions::normal_mag() const { return normal_mag_; }
void EmbedSearchParamsOptions::set_normal_mag( core::Size setting ) {
	normal_mag_ = setting;
}

/// @brief Getter/Setter Pair for Normal Cycles
/// @details Param Type: abs(int)
/// Precodnition setting >=0
core::Size EmbedSearchParamsOptions::normal_cycles() const { return normal_cycles_; }
void EmbedSearchParamsOptions::set_normal_cycles( core::Size setting ) { normal_cycles_ = setting; }

/// @brief Getter/Setter pair for penalties
/// @details Param tyupe: bool
bool EmbedSearchParamsOptions::penalties() const { return penalties_; }
void EmbedSearchParamsOptions::set_penalties( bool setting ) { penalties_ = setting; }

/// @brief Getter/Setter pair for interpolating mpair term
/// @details Param Type: bool
bool EmbedSearchParamsOptions::no_interpolate_mpair() const { return no_interpolate_mpair_; }
void EmbedSearchParamsOptions::set_no_interpolate_mpair( bool setting ) { no_interpolate_mpair_ = setting; }

/// @brief Creator Class - Return options class type
std::string
EmbedSearchParamsOptionsCreator::options_type() const { return "EmbedSearchParamsOptions"; }

/// @brief Creator class - return options class
basic::resource_manager::ResourceOptionsOP
EmbedSearchParamsOptionsCreator::create_options() const { return new EmbedSearchParamsOptions; }

} // io
} // memrbane
} // core

#endif // INCLUDED_core_membrane_io_EmbedSearchParamsOptions_cc


