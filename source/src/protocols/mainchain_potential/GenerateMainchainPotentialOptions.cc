// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/mainchain_potential/GenerateMainchainPotentialOptions.cc
/// @brief Options container for the generator for mainchain potentials.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project includes:
#include <protocols/mainchain_potential/GenerateMainchainPotentialOptions.hh>

// Basic includes:
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/make_mainchain_potential.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// Utility includes:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.mainchain_potential.GenerateMainchainPotentialOptions" );


/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// STATIC FUNCTIONS //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

/// @brief Static function used to register releavnt options.
void
protocols::mainchain_potential::GenerateMainchainPotentialOptions::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( basic::options::OptionKeys::make_mainchain_potential::do_minimization );
	option.add_relevant( basic::options::OptionKeys::run::min_type );
	option.add_relevant( basic::options::OptionKeys::run::min_tolerance );
	option.add_relevant( basic::options::OptionKeys::make_mainchain_potential::mainchain_potential_points_per_dimension );
	option.add_relevant( basic::options::OptionKeys::make_mainchain_potential::mainchain_torsions_covered );
	option.add_relevant( basic::options::OptionKeys::make_mainchain_potential::make_pre_proline_potential );
	option.add_relevant( basic::options::OptionKeys::make_mainchain_potential::output_filename );
	option.add_relevant( basic::options::OptionKeys::make_mainchain_potential::residue_name );
	option.add_relevant( basic::options::OptionKeys::make_mainchain_potential::symmetrize_output );
	option.add_relevant( basic::options::OptionKeys::make_mainchain_potential::write_potentials_for_individual_scoreterms );
	option.add_relevant( basic::options::OptionKeys::score::weights );
	option.add_relevant( basic::options::OptionKeys::packing::ex1::ex1 );
	option.add_relevant( basic::options::OptionKeys::packing::ex1::level );
	option.add_relevant( basic::options::OptionKeys::packing::ex2::ex2 );
	option.add_relevant( basic::options::OptionKeys::packing::ex2::level );
	option.add_relevant( basic::options::OptionKeys::packing::ex3::ex3 );
	option.add_relevant( basic::options::OptionKeys::packing::ex3::level );
	option.add_relevant( basic::options::OptionKeys::packing::ex4::ex4 );
	option.add_relevant( basic::options::OptionKeys::packing::ex4::level );
	option.add_relevant( basic::options::OptionKeys::packing::extrachi_cutoff);
}




namespace protocols {
namespace mainchain_potential {

/// @brief Default constructor.
/// @details "True" indicates that we should initialize from the global options system.  "False" by default,
/// meaning that all options are set to default values.
GenerateMainchainPotentialOptions::GenerateMainchainPotentialOptions( bool const initialize_from_globals /* = false */ ):
	utility::pointer::ReferenceCount(),
	residue_name_( "ALA" ),
	scorefxn_filename_(""), //Empty string means use default.
	make_pre_proline_potential_(false),
	dimensions_(utility::vector1< core::Size >({36,36})),
	mainchain_torsions_covered_(),
	kbt_(0.63),
	do_minimization_(true),
	minimization_type_("dfpmin"),
	minimization_threshold_(1.0e-7),
	output_filename_("generated_mainchain_potential.txt"),
	write_potentials_for_individual_scoreterms_(false),
	symmetrize_output_(false)
{
	if ( initialize_from_globals ) {
		do_initialization_from_globals();
	}
}

/// @brief Destructor.
GenerateMainchainPotentialOptions::~GenerateMainchainPotentialOptions(){}

/// @brief Clone function: create a copy of this object, and return an owning pointer to the copy.
GenerateMainchainPotentialOptionsOP
GenerateMainchainPotentialOptions::clone() const {
	return GenerateMainchainPotentialOptionsOP( utility::pointer::make_shared< GenerateMainchainPotentialOptions >( *this ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// PUBLIC FUNCTIONS //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

/// @brief Set the name of the residue type for which we'll be generating a mainchain potential.
void
GenerateMainchainPotentialOptions::set_residue_name(
	std::string const &name_in
) {
	runtime_assert_string_msg( !name_in.empty(), "Error in GenerateMainchainPotentialOptions::set_residue_name(): The name cannot be empty!" );
	residue_name_ = name_in;
}

/// @brief Set whether we're generating a pre-proline potential.
/// @details If true, then we generate a potential for a position before a proline, sarcosine, or other N-substituted building block (e.g. a peptoid).  If false,
/// then the next residue's nitrogen is unsubstituted.
/// @note Determines the patch applied to the C-terminus.
void
GenerateMainchainPotentialOptions::set_make_pre_proline_potential(
	bool const setting
) {
	make_pre_proline_potential_ = setting;
}

/// @brief Set the vector of number of points in each dimension for the mainchain potential that we're generating.
/// @details This takes a vector of ints, not Sizes, because that's what comes from the global options system.
void
GenerateMainchainPotentialOptions::set_dimensions(
	utility::vector1< int > const & dimensions_in
) {
	static const std::string errmsg( "Error when parsing the make_mainchain_potential application's \"-mainchain_potential_points_per_dimension\" flag:" );
	core::Size const dimensionality( dimensions_in.size() );
	runtime_assert_string_msg( dimensionality != 0, errmsg + "The flag was used, but no values were specified.  One value must be provided for each dimension of the mainchain potential to be generated.  The mainchain potential cannot be zero-dimensional!" );
	dimensions_.resize( dimensionality );
	for ( core::Size i(1); i<=dimensionality; ++i ) {
		runtime_assert_string_msg( dimensions_in[i] > 0, errmsg + "All entries provided with this flag must be positive.  The number of points for a given dimension cannot be zero (or negative)." );
		dimensions_[i] = static_cast<core::Size>(dimensions_in[i]); //Cast to lose sign.
	}
}

/// @brief Set which mainchain torsions are covered.
/// @details If this is an empty vector, it indicates that all mainchain torsions are covered.  This takes a vector of
/// ints, not Sizes, because that's what comes from the global options sytem.
/// @note Pass this an empty vector to set the default (all mainchain torsions for the residue type.)
void
GenerateMainchainPotentialOptions::set_mainchain_torsions_covered(
	utility::vector1< int > const & mainchain_torsions_in
) {
	static const std::string errmsg( "Error in GenerateMainchainPotentialOptions::set_mainchain_torsions_covered(): ");
	core::Size const ntors( mainchain_torsions_in.size() );
	mainchain_torsions_covered_.clear();
	if ( ntors == 0 ) return;
	mainchain_torsions_covered_.reserve(ntors);
	for ( core::Size i(1); i<=ntors; ++i ) {
		runtime_assert_string_msg( mainchain_torsions_in[i] > 0, errmsg + "A mainchain torsion index was specified that was not positive." );
		runtime_assert_string_msg( !mainchain_torsions_covered_.has( static_cast<core::Size>(mainchain_torsions_in[i]) ), errmsg + "A duplicated torsion index was found in the list of torsion indices." );
		mainchain_torsions_covered_.push_back( static_cast< core::Size >( mainchain_torsions_in[i] ) );
	}
}

/// @brief Set name of the scorefunction that we'll use.  An empty string indicates that the default should be used.
void
GenerateMainchainPotentialOptions::set_scorefxn_filename(
	std::string const & filename_in
) {
	scorefxn_filename_ = filename_in;
}

/// @brief Set the Boltzmann temperature.
void
GenerateMainchainPotentialOptions::set_kbt(
	core::Real const & kbt_in
) {
	runtime_assert_string_msg( kbt_in >= 0.0, "Error in GenerateMainchainPotentialOptions::set_kbt(): The Boltzmann temperature cannot be negative." );
	kbt_ = kbt_in;
}

/// @brief Set whether we're minimizing each rotamer.
void
GenerateMainchainPotentialOptions::set_do_minimization( bool const setting ) {
	do_minimization_ = setting;
}

/// @brief Set the minimization type.
void
GenerateMainchainPotentialOptions::set_minimization_type( std::string const & type_in ) {
	minimization_type_ = type_in;
}

/// @brief Set the minimization threshold.
void
GenerateMainchainPotentialOptions::set_minimization_threshold( core::Real const & threshold_in ) {
	runtime_assert_string_msg( threshold_in > 0, "Error in GenerateMainchainPotentialOptions::set_minimization_threshold(): The minimization threshold must be greater than zero." );
	minimization_threshold_ = threshold_in;
}

/// @brief Set the mainchain potential filename to write.
void
GenerateMainchainPotentialOptions::set_output_filename(
	std::string const & filename_in
) {
	runtime_assert_string_msg( !filename_in.empty(), "Error in GenerateMainchainPotentialOptions::set_output_filename(): The output filename cannot be empty!" );
	TR << "Set output filename to \"" << filename_in << "\"." << std::endl;
	output_filename_ = filename_in;
}

/// @brief Set whether we are writing mainchain potentials for individual scoreterms.
void
GenerateMainchainPotentialOptions::set_write_potentials_for_individual_scoreterms(
	bool const setting
) {
	if ( setting ) {
		TR << "Potentials will be written for individual scoreterms.  Note that these will be unnormalized, and score terms will be treated as though they have a weight of 1.0." << std::endl;
	}
	write_potentials_for_individual_scoreterms_ = setting;
}

/// @brief Set whether the output should be made symmetric.
void
GenerateMainchainPotentialOptions::set_symmetrize_output(
	bool const setting
) {
	if ( setting ) {
		TR << "The generated mainchain potential will be made symmetric." << std::endl;
	}
	symmetrize_output_ = setting;
}

//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// PRIVATE FUNCTIONS /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

/// @brief Initialize this object from the global options system.
void
GenerateMainchainPotentialOptions::do_initialization_from_globals() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static const std::string errmsg( "Error in make_mainchain_potential application input options: " );

	runtime_assert_string_msg( !option[ make_mainchain_potential::residue_name ]().empty(), errmsg + "The residue type string cannot be empty.  Please specify a residue type with the \"-make_mainchain_potential:residue_name\" flag." );
	residue_name_ = option[ make_mainchain_potential::residue_name ]();

	make_pre_proline_potential_ = option[make_mainchain_potential::make_pre_proline_potential].value();

	runtime_assert_string_msg( option[ make_mainchain_potential::mainchain_potential_points_per_dimension ].user(), errmsg + "The \"-mainchain_potential_points_per_dimension\" flag was not provided.  Please indicate the number of points for each dimension of the mainchain potential to be generated.  For example, if we were making a 2D mainchain potential with 72 points in the first dimension and 36 in the second, the input would be \"-mainchain_potential_points_per_dimension 72 36\"." );
	set_dimensions( option[make_mainchain_potential::mainchain_potential_points_per_dimension].value() );

	if ( option[make_mainchain_potential::mainchain_torsions_covered].user() ) {
		runtime_assert_string_msg( !option[make_mainchain_potential::mainchain_torsions_covered]().empty(), errmsg + "The \"-mainchain_torsions_covered\" option was specified, but no torsion indices were given." );
		set_mainchain_torsions_covered( option[make_mainchain_potential::mainchain_torsions_covered]() );
	} else {
		mainchain_torsions_covered_.clear();
	}

	set_output_filename( option[make_mainchain_potential::output_filename]() );
	if ( !option[ make_mainchain_potential::output_filename ].user() ) {
		TR.Warning << TR.Red << "Warning!  No output filename was specified.  Using the default (\"" << output_filename() << "\") instead.  It is strongly recommended that you provide a custom output filename!" << TR.Reset << std::endl;
	}

	set_do_minimization( option[make_mainchain_potential::do_minimization ].value() );
	set_minimization_threshold( option[run::min_tolerance ]() );
	set_minimization_type( option[run::min_type]() );

	if ( option[score::weights].user() ) {
		set_scorefxn_filename( option[score::weights]() );
	} else {
		set_scorefxn_filename(""); //Empty string means use default.
	}

	set_write_potentials_for_individual_scoreterms( option[make_mainchain_potential::write_potentials_for_individual_scoreterms]() );

	set_symmetrize_output( option[make_mainchain_potential::symmetrize_output]() );
}


} //protocols
} //mainchain_potential
