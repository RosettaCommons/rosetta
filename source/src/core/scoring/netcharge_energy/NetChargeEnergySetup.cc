// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/netcharge_energy/NetChargeEnergySetup.cc
/// @brief A helper for the EnergyMethod that penalizes deviation from a desired net charge.
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

// Unit headers
#include <core/scoring/netcharge_energy/NetChargeEnergySetup.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>

// Options system
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// File I/O
#include <basic/database/open.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>

// Other Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <cmath>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace netcharge_energy {

static basic::Tracer TR("core.scoring.netcharge_energy.NetChargeEnergySetup");


/**************************************************
NetChargeEnergySetup class:
**************************************************/

/// @brief Default constructor for NetChargeEnergySetup.
///
NetChargeEnergySetup::NetChargeEnergySetup() :
	utility::pointer::ReferenceCount(),
	desired_charge_(0),
	charge_penalties_range_( std::pair< signed long int, signed long int>(0, 0) ),
	penalties_(),
	tailfunctions_( std::pair < TailFunction, TailFunction >( tf_quadratic, tf_quadratic ) )
{}


/// @brief Default destructor for NetChargeEnergySetup.
///
NetChargeEnergySetup::~NetChargeEnergySetup() = default;

/// @brief Clone: create a copy of this object, and return an owning pointer
/// to the copy.
NetChargeEnergySetupOP NetChargeEnergySetup::clone() const {
	return NetChargeEnergySetupOP( new NetChargeEnergySetup(*this) );
}

/// @brief Reset all data in this data storage object.
///
void NetChargeEnergySetup::reset() {
	desired_charge_ = 0;
	charge_penalties_range_ = std::pair< signed long int, signed long int>(0,0);
	penalties_.clear();
	tailfunctions_ = std::pair< TailFunction, TailFunction >( tf_quadratic, tf_quadratic );
}

/// @brief Initialize from a .charge file.
///
void NetChargeEnergySetup::initialize_from_file( std::string const &filename ) {
	using namespace utility::io;

	//First, reset all data:
	reset();

	std::string filename_formatted = filename;
	if ( utility::file::file_extension(filename_formatted)!="charge" ) filename_formatted+= ".charge";

	izstream infile;
	infile.open( filename_formatted );
	if ( !infile.good() ) {
		filename_formatted = "scoring/score_functions/netcharge/" + utility::file::file_basename(filename_formatted) + ".charge";
		basic::database::open( infile, filename_formatted );
		runtime_assert_string_msg( infile.good(), "Error in core::scoring::netcharge_energy::NetChargeEnergySetup::initialize_from_file():  Unable to open .charge file for read!" );
	}

	if ( TR.visible() ) TR << "Reading netcharge scoring term setup data from " << filename_formatted << "." << std::endl;
	std::string curline(""); //Buffer for current line.
	utility::vector1< std::string > lines; //Storing all lines

	//Read the file:
	while ( getline(infile, curline) ) {
		if ( curline.size() < 1 ) continue; //Ignore blank lines.
		//Find and process comments:
		std::string::size_type pound = curline.find('#', 0);
		if ( pound == std::string::npos ) {
			lines.push_back( ObjexxFCL::strip( curline, " \t\n" ) );
		} else {
			curline = curline.substr(0,pound);
			lines.push_back( ObjexxFCL::strip( curline, " \t\n" ) );
		}
	}
	infile.close();

	if ( TR.Debug.visible() ) TR.Debug << "Read complete.  Parsing charge penalty definitions." << std::endl;
	parse_a_penalty_definition( lines );

	check_data();

	return;
}

/// @brief Initialize from a string in the format of a .charge file.
/// @details Allows external code to initialize object without having it read
/// directly from disk.
void
NetChargeEnergySetup::initialize_from_file_contents(
	std::string const &filecontents
) {
	reset();

	std::istringstream filestream(filecontents);

	std::string curline(""); //Buffer for current line.
	utility::vector1< std::string > lines; //Storing all lines

	//Split the string into lines:
	while ( getline(filestream, curline) ) {
		if ( curline.size() < 1 ) continue; //Ignore blank lines.
		//Find and process comments:
		std::string::size_type pound = curline.find('#', 0);
		if ( pound == std::string::npos ) {
			lines.push_back( ObjexxFCL::strip( curline, " \t\n" ) );
		} else {
			curline = curline.substr(0,pound);
			lines.push_back( ObjexxFCL::strip( curline, " \t\n" ) );
		}
	}

	if ( TR.Debug.visible() ) TR.Debug << "Initial processing of file contents string complete.  Parsing penalty definitions." << std::endl;
	parse_a_penalty_definition( lines );

	check_data();
}

/// @brief Get tail function name from enum.
///
std::string NetChargeEnergySetup::get_tailfunction_name( TailFunction const tf ) const
{
	switch(tf) {
	case tf_linear :
		return "LINEAR";
	case tf_quadratic :
		return "QUADRATIC";
	case tf_constant :
		return "CONSTANT";
	default :
		break;
	}

	return "UNKNOWN";
}

/// @brief Get tail function enum from name.
/// @details This is slow; it calls get_tailfunction_name repeatedly.  Intended only for use during setup.
/// Returns tf_unknown if tail function couldn't be parsed.
TailFunction NetChargeEnergySetup::get_tailfunction_from_name( std::string const &name ) const {
	for ( core::Size i=1; i<tf_end_of_list; ++i ) {
		if ( get_tailfunction_name( static_cast<TailFunction>(i) ) == name ) return static_cast<TailFunction>(i);
	}
	return tf_unknown;
}

/// @brief Get a summary of the data stored in this object
///
std::string NetChargeEnergySetup::report() const {
	std::stringstream output("");

	output << "desired_charge_:\t" << desired_charge_ << std::endl;

	output << "charge_penatly_range_:\t" << charge_penalties_range_.first << ", " << charge_penalties_range_.second << std::endl;

	output << "penalties_:\t";
	for ( core::Size i=1, imax=penalties_.size(); i<=imax; ++i ) {
		output << penalties_[i];
		if ( i<imax ) output << "\t";
	}
	output << std::endl;

	output << "tailfunctions_:\t" << get_tailfunction_name( tailfunctions_.first ) << ", " << get_tailfunction_name( tailfunctions_.second ) << std::endl;

	return output.str();
}


/// @brief Parse out penalty definition from a single block of lines from file.
/// @details The lines vector should be the lines BETWEEN the PENALTY_DEFINITION and END_PENALTY_DEFINITION lines.
void NetChargeEnergySetup::parse_a_penalty_definition( utility::vector1 < std::string > const & lines ) {

	static const std::string errmsg("Error in core::scoring::netcharge_energy::NetChargeEnergySetup::parse_a_penalty_definition():  ");

	// Bools for whether we've found all the lines we're looking for:
	bool netchargefound(false);
	bool rangefound(false);
	bool penaltiesfound(false);

	signed long int desired_charge(0), charge_penalties_range_min(0), charge_penalties_range_max(0);
	utility::vector1< core::Real > penalties;
	TailFunction beforefxn( tf_unknown ); //TailFunction is an enum defined in NetChargeEnergySetup.hh
	TailFunction afterfxn( tf_unknown );

	core::Size const nlines( lines.size() ); //Number of lines we'll be going through.

	for ( core::Size i=1; i<=nlines; ++i ) { //Loop through all lines
		std::istringstream curline(lines[i]);
		std::string oneword("");

		curline >> oneword;
		if ( curline.fail() ) continue; //A blank line.

		if ( !oneword.compare( "DESIRED_CHARGE" ) ) {
			runtime_assert_string_msg( !netchargefound, errmsg + "A \"DESIRED_CHARGE\" line can only be present only once in a netcharge energy setup file." );
			curline >> desired_charge;
			runtime_assert_string_msg( !curline.fail(), errmsg + "Could not parse DESIRED_CHARGE line." );
			netchargefound=true;
		} else if ( !oneword.compare( "PENALTIES_CHARGE_RANGE" ) ) {
			runtime_assert_string_msg( !rangefound, errmsg + "A \"PENALTIES_CHARGE_RANGE\" line can only be present only once in a netcharge energy setup file." );
			curline >> charge_penalties_range_min >> charge_penalties_range_max;
			runtime_assert_string_msg( !curline.fail(), errmsg + "Could not parse PENALTIES_CHARGE_RANGE line." );
			rangefound=true;
		} else if ( !oneword.compare( "PENALTIES" ) ) {
			runtime_assert_string_msg( !penaltiesfound, errmsg + "A \"PENALTIES\" line can only be present only once in a netcharge energy setup file." );
			bool at_least_one_found(false);
			core::Real penval;
			do {
				curline >> penval;
				if ( curline.fail() ) {
					at_least_one_found = false;
					break;
				}
				at_least_one_found = true;
				penalties.push_back(penval);
			} while( !curline.eof() );
			runtime_assert_string_msg( at_least_one_found, errmsg + "Could not parse PENALTIES line." );
			penaltiesfound=true;
		} else if ( !oneword.compare( "BEFORE_FUNCTION" ) ) {
			runtime_assert_string_msg( beforefxn == tf_unknown, errmsg + "A \"BEFORE_FUNCTION\" line can only be present only once in a netcharge energy setup file." );
			std::string beforefxn_str;
			curline >> beforefxn_str;
			runtime_assert_string_msg( !curline.fail(), errmsg + "Could not parse BEFORE_FUNCTION line." );
			beforefxn = get_tailfunction_from_name( beforefxn_str );
			runtime_assert_string_msg( beforefxn != tf_unknown, errmsg + "Could not parse BEFORE_FUNCTION line." );
		} else if ( !oneword.compare( "AFTER_FUNCTION" ) ) {
			runtime_assert_string_msg( afterfxn == tf_unknown, errmsg + "An \"AFTER_FUNCTION\" line can only be present only once in a netcharge energy setup file." );
			std::string afterfxn_str;
			curline >> afterfxn_str;
			runtime_assert_string_msg( !curline.fail(), errmsg + "Could not parse AFTER_FUNCTION line." );
			afterfxn = get_tailfunction_from_name( afterfxn_str );
			runtime_assert_string_msg( afterfxn != tf_unknown, errmsg + "Could not parse AFTER_FUNCTION line." );
		}
	} //End loop through all lines

	runtime_assert_string_msg( netchargefound, errmsg + "No DESIRED_CHARGE line was found." );
	runtime_assert_string_msg( rangefound, errmsg + "No PENALTIES_CHARGE_RANGE line was found." );
	runtime_assert_string_msg( penaltiesfound, errmsg + "No PENALTIES line was found." );

	// If a BEFORE_FUNCTION line was not provided, default to quadratic:
	if ( beforefxn == tf_unknown ) { beforefxn = tf_quadratic;}
	// If an AFTER_FUNCTION line was not provided, default to quadratic:
	if ( afterfxn == tf_unknown ) { afterfxn = tf_quadratic;}

	//Store everything that we've parsed:
	desired_charge_ = desired_charge;
	charge_penalties_range_ = std::pair< signed long int, signed long int >( charge_penalties_range_min, charge_penalties_range_max );
	penalties_ = penalties;
	tailfunctions_ = std::pair< TailFunction, TailFunction >( beforefxn, afterfxn );

	check_data();

	return;
}

/// @brief Do some final checks to ensure that data were loaded properly.
///
void NetChargeEnergySetup::check_data() const {
	if ( TR.Debug.visible() ) TR.Debug << "Checking loaded data." << std::endl;
	static const std::string errmsg("Error in core::scoring::netcharge_energy::NetChargeEnergySetup::check_data():  ");

	runtime_assert_string_msg( penalties_.size() >=3, errmsg + "At least three penalty values must be specified.  Too few were found." );
	runtime_assert_string_msg( charge_penalties_range_.first < desired_charge_, errmsg + "The start of the charge penalties range must be less than the desired charge." );
	runtime_assert_string_msg( charge_penalties_range_.second > desired_charge_, errmsg + "The end of the charge penalties range must be greater than the desired charge." );
	runtime_assert_string_msg( charge_penalties_range_.second - charge_penalties_range_.first + 1 == static_cast< signed long int>( penalties_.size() ), errmsg + "The number of penalty values provided must equal the number of charge values in the range." );

	if ( TR.Debug.visible() ) {
		TR.Debug << "Data checks passed!" << std::endl;
		TR.Debug << report() << std::endl;
	}
	return;
}

} // netcharge_energy
} // scoring
} // core

#ifdef    SERIALIZATION
template< class Archive >
void
core::scoring::netcharge_energy::NetChargeEnergySetup::save( Archive & arc ) const {
	arc( CEREAL_NVP( desired_charge_ ) );
	arc( CEREAL_NVP( charge_penalties_range_ ) );
	arc( CEREAL_NVP( penalties_ ) );
	arc( CEREAL_NVP( tailfunctions_ ) );
}

template< class Archive >
void
core::scoring::netcharge_energy::NetChargeEnergySetup::load( Archive & arc ) {
	arc( desired_charge_ );
	arc( charge_penalties_range_ );
	arc( penalties_ );
	arc( tailfunctions_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::netcharge_energy::NetChargeEnergySetup );
CEREAL_REGISTER_TYPE( core::scoring::netcharge_energy::NetChargeEnergySetup )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_netcharge_energy_NetChargeEnergySetup )
#endif // SERIALIZATION
