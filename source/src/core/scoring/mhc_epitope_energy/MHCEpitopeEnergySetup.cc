// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/mhc_epitope_energy/MHCEpitopeEnergySetup.cc
/// @brief A helper for the MHCEpitopeEnergy
/// @details Follows analogous file for Vikram K. Mulligan's NetChargeEnergy
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

// Unit headers
#include <core/scoring/mhc_epitope_energy/MHCEpitopeEnergySetup.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictor.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorExternal.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorMatrix.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorPreLoaded.hh>

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
#include <utility/pointer/memory.hh>
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
namespace mhc_epitope_energy {

static basic::Tracer TR("core.scoring.mhc_epitope_energy.MHCEpitopeEnergySetup");


/**************************************************
MHCEpitopeEnergySetup class:
**************************************************/

/// @brief Default constructor for MHCEpitopeEnergySetup.
///
MHCEpitopeEnergySetup::MHCEpitopeEnergySetup() :
	utility::pointer::ReferenceCount()
{}


/// @brief Default destructor for MHCEpitopeEnergySetup.
///
MHCEpitopeEnergySetup::~MHCEpitopeEnergySetup() = default;

/// @brief Clone: create a copy of this object, and return an owning pointer
/// to the copy.
MHCEpitopeEnergySetupOP MHCEpitopeEnergySetup::clone() const {
	return utility::pointer::make_shared< MHCEpitopeEnergySetup >(*this);
}

/// @brief Reset all data in this data storage object.
///
void MHCEpitopeEnergySetup::reset() {
}

/// @brief Initialize from a .mhc file.
///
void MHCEpitopeEnergySetup::initialize_from_file( std::string const &filename ) {
	using namespace utility::io;

	//First, reset all data:
	reset();

	std::string filename_formatted = filename;
	if ( utility::file::file_extension(filename_formatted)!="mhc" ) filename_formatted+= ".mhc";

	izstream infile;
	infile.open( filename_formatted );
	if ( !infile.good() ) {
		TR << "Cannot find .mhc file " << filename_formatted << ".  Checking the Rosetta database." << std::endl;
		filename_formatted = "scoring/score_functions/mhc_epitope/" + utility::file::file_basename(filename_formatted) + ".mhc";
		basic::database::open( infile, filename_formatted );
		runtime_assert_string_msg( infile.good(), "Error in core::scoring::mhc_epitope_energy::MHCEpitopeEnergySetup::initialize_from_file():  Unable to open .mhc file for read!" );
	}

	if ( TR.visible() ) TR << "Reading mhc_epitope scoring term setup data from " << filename_formatted << "." << std::endl;
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

	if ( TR.Debug.visible() ) TR.Debug << "Read complete.  Parsing epitope prediction specs." << std::endl;
	parse_specs( lines );

	return;
}

void
MHCEpitopeEnergySetup::initialize_from_file_contents(
	std::string const &filecontents
) {
	reset();

	std::istringstream filestream(filecontents);
	if ( TR.visible() ) TR << "Reading the following mhc_epitope scoring term setup data:" << std::endl << filecontents << std::endl;

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

	if ( TR.Debug.visible() ) TR.Debug << "Initial processing of file contents string complete.  Parsing epitope prediction specs." << std::endl;
	parse_specs( lines );
}

/// @brief The MHC epitope score for the peptide, as returned by the predictor
core::Real MHCEpitopeEnergySetup::raw_score(std::string peptide) const {
	if ( is_default() ) { //If we haven't initalized (i.e. this is the "default"), print a warning and return 0.
		TR.Warning << "Attempting to generate an mhc_epitope score without a setup file.  Did you forget to include a setup .mhc file?" << std::endl;
		TR.Warning << "mhc_epitope will return a score of 0." << std::endl;
		return 0;
	}

	return predictor_->score(peptide);
}

/// @brief Get a summary of the data stored in this object
///
std::string MHCEpitopeEnergySetup::report() const {
	std::stringstream output("");

	// Predictor
	output << predictor_->report();

	// Score transform
	output << "; xform ";
	if ( !relative_ ) output << "raw ";
	else {
		if ( relative_additive_ ) output << "additive ";
		else output << "multiplicative ";
	}
	if ( apply_offset_ ) output << score_offset_;

	output << std::endl;

	return output.str();
}

/// @details
void MHCEpitopeEnergySetup::parse_specs( utility::vector1 < std::string > const & lines ) {

	static const std::string errmsg("Error in core::scoring::mhc_epitope_energy::MHCEpitopeEnergySetup::parse_specs():  ");

	core::Size const nlines( lines.size() ); //Number of lines we'll be going through.

	TR.Debug << "parsing" << std::endl;

	// TODO: very fragile -- do more error checking and friendly error reporting

	// The first line of the .mhc file contains the method keyword, the method type, and for the filename of the reference data.

	// Get method from first line
	std::istringstream firstline(lines[1]);
	std::string keyword, method, filename;
	firstline >> keyword;
	firstline >> method;
	TR.Debug << "method: '" << method  << "'" << std::endl;

	// Hold on to predictor as its subclass, to set options without typecasting
	MHCEpitopePredictorExternalOP pred_external;
	MHCEpitopePredictorMatrixOP pred_matrix;
	MHCEpitopePredictorPreLoadedOP pred_pre;
	if ( method == "matrix" ) {
		// For matrix methods, read the filename containing the matrix, and use it to create a new MHCEpitopePredictorMatrix object.
		firstline >> filename;
		pred_matrix = MHCEpitopePredictorMatrixOP(utility::pointer::make_shared<MHCEpitopePredictorMatrix>(filename));
		predictor_ = pred_matrix;
	} else if ( method == "external" ) {
		// We don't know how the database will work in a multi-threading environment, so we are disabling it in case.
#ifdef MULTI_THREADED // Once this issue is resolved, the #ifndef blocks in the unit test suite should be removed (MHCEpitopePredictorExternal.cxxtest.hh, MHCEpitopeEnergy.cxxtest.hh, MHCEpitopeEnergySetup.cxxtest.hh)
		utility_exit_with_message( "The MHCEpitopePredictorExternal is not supported in a multi-threading environment at this time.  Contact cbk@cs.dartmouth.edu and brahm.yachnin@rutgers.edu if you would like to help test it with multi-threading." );
#endif
		// For external database methods, read the filename containing the database, and use it to create a new MHCEpitopePredictorExternal object.
		firstline >> filename;
		pred_external = MHCEpitopePredictorExternalOP(utility::pointer::make_shared<MHCEpitopePredictorExternal>(filename));
		predictor_ = pred_external;
	} else if ( method == "preloaded" ) {
		pred_pre = MHCEpitopePredictorPreLoadedOP(utility::pointer::make_shared<MHCEpitopePredictorPreLoaded>());
		predictor_ = pred_pre;
		std::string filetype;
		firstline >> filetype >> filename;
		if ( filetype == "db" ) {
			pred_pre->load_database(filename);
		} else if ( filetype == "csv" ) {
			pred_pre->load_csv(filename);
		} else {
			utility_exit_with_message("ERROR: Unknown preloaded file type " + filetype);
		}
	} else {
		utility_exit_with_message("ERROR: Unknown epitope prediction method " + method);
	}

	for ( core::Size i=2; i<=nlines; ++i ) { //Loop through all remaining lines, containing options
		TR.Debug << "option: '" << lines[i] << "'" << std::endl;
		std::istringstream line(lines[i]);
		std::string option;
		line >> option;
		if ( option == "alleles" ) {
			// Use only the alleles listed after the alleles keyword.
			// Use * to use all alleles in the matrix file (default behaviour).
			// TODO -- allow restriction
		} else if ( option == "threshold" && method=="matrix" ) {
			// If the threshold keyword is used, set the threshold in the predictor.
			core::Real thresh;
			line >> thresh;
			pred_matrix->set_thresh(thresh);
		} else if ( option == "unseen" && (method=="external" || method=="preloaded") ) {
			// How to handle unseen epitopes.
			// The line gives the type of handler and its parameters
			std::string handler;
			line >> handler;
			core::Size unseen(0);
			if ( handler == "ignore" ) {
				// Useful with experimental data, where if you haven't seen it, it's a good thing
				unseen = 0; // already was, but now it's doubly so
			} else if ( handler == "penalize" || handler == "score" ) {
				line >> unseen;
			} else {
				TR.Error << "Unknown unseen handler " << handler << "; ignoring: " << lines[i] << std::endl;
			}
			if ( method=="external" ) pred_external->set_unseen_score(unseen);
			else if ( method=="preloaded" ) pred_pre->set_unseen_score(unseen);
			// those are the only two possibilities now, but put the else if for future expansion
		} else if ( option == "xform" ) {
			// A transformation to the score returned by the predictor
			std::string xform;
			line >> xform;
			if ( xform == "raw" ) {
				// Use the score as given, subject to a scoring offset
				relative_ = false;
				line >> score_offset_;
				apply_offset_ = true;
			} else if ( xform == "relative+" ) {
				// Relative to the native score, plus a scoring offset
				relative_ = true;
				relative_additive_ = true;
				line >> score_offset_;
				apply_offset_ = true;
			} else if ( xform == "relative*" ) {
				// Relative to the native score, with multiplicative offset (i.e., native * score_offset_)
				relative_ = true;
				relative_additive_ = false;
				line >> score_offset_;
				apply_offset_ = true;
			} else {
				TR.Error << "Unknown xform; ignoring:" << lines[i] << std::endl;
			}
		} else {
			TR.Error << "Unknown option; ignoring:" << lines[i] << std::endl;
		}
	} //End loop through all lines

	TR << "Parsing complete.  MHC setup settings: " << std::endl;
	TR << report();
	TR.flush();
}

core::Real MHCEpitopeEnergySetup::xform(core::Real raw, core::Real native) const
{
	if ( !relative_ ) {
		// score as is; check offset
		if ( apply_offset_ && raw < score_offset_ ) return 0;
		return raw - score_offset_;
	}

	// score relative to native
	// If we aren't setting a offset, just return raw - native
	if ( !apply_offset_ ) {
		return raw - native;
	}

	core::Real xformed;
	// check additive or multiplicative offset
	if ( relative_additive_ ) {
		// If we're in additive mode, calculate xformed appropriately.
		xformed = raw - native + score_offset_;
	} else { // if not additive, then multiplicative
		// Calculate xformed in multiplicative mode
		xformed = raw - native * score_offset_;
	}

	// If xformed is negative, return 0.  Otherwise, return xformed.
	if ( xformed < 0 ) {
		return 0;
	} else {
		return xformed;
	}
}

bool MHCEpitopeEnergySetup::operator == ( MHCEpitopeEnergySetup const & other ) const
{
	// Check member variables
	if ( relative_ != other.get_relative() ) return false;
	if ( apply_offset_ != other.get_apply_offset() ) return false;
	if ( relative_additive_ != other.get_relative_additive() ) return false;
	if ( score_offset_ != other.get_score_offset() ) return false;
	if ( predictor_ == NULL && other.predictor_ != NULL ) return false;
	if ( predictor_ != NULL && !(*predictor_ == *other.predictor_) ) return false;
	return true;
}

} // mhc_epitope_energy
} // scoring
} // core

#ifdef    SERIALIZATION
template< class Archive >
void
core::scoring::mhc_epitope_energy::MHCEpitopeEnergySetup::save( Archive & arc ) const {
	arc( CEREAL_NVP( predictor_ ) );
	arc( CEREAL_NVP( relative_ ) );
	arc( CEREAL_NVP( apply_offset_ ) );
	arc( CEREAL_NVP( relative_additive_ ) );
	arc( CEREAL_NVP( score_offset_ ) );
}

template< class Archive >
void
core::scoring::mhc_epitope_energy::MHCEpitopeEnergySetup::load( Archive & arc ) {
	arc( predictor_ );
	arc( relative_ );
	arc( apply_offset_ );
	arc( relative_additive_ );
	arc( score_offset_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::mhc_epitope_energy::MHCEpitopeEnergySetup );
CEREAL_REGISTER_TYPE( core::scoring::mhc_epitope_energy::MHCEpitopeEnergySetup )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_mhc_epitope_energy_MHCEpitopeEnergySetup )
#endif // SERIALIZATION
