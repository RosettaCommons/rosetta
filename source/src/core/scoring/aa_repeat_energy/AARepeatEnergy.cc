// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/aa_repeat_energy/AARepeatEnergy.cc
/// @brief An EnergyMethod that penalizes stretches of a repeating amino acid (e.g. poly-Q sequences).
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

// Unit headers
#include <core/scoring/aa_repeat_energy/AARepeatEnergy.hh>
#include <core/scoring/aa_repeat_energy/AARepeatEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/EnergyMethod.hh>
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

namespace core {
namespace scoring {
namespace aa_repeat_energy {

using namespace core::scoring::methods;
static THREAD_LOCAL basic::Tracer TR("core.scoring.methods.AARepeatEnergy");

/// @brief This must return a fresh instance of the AARepeatEnergy class, never an instance already in use.
///
methods::EnergyMethodOP
AARepeatEnergyCreator::create_energy_method( methods::EnergyMethodOptions const & ) const
{
	return methods::EnergyMethodOP( new AARepeatEnergy );
}

/// @brief Defines the score types that this energy method calculates.
///
ScoreTypes
AARepeatEnergyCreator::score_types_for_method() const
{
	ScoreTypes sts;
	sts.push_back( aa_repeat );
	return sts;
}

/// @brief Default constructor.
///
AARepeatEnergy::AARepeatEnergy() :
	parent1( methods::EnergyMethodCreatorOP( new AARepeatEnergyCreator ) ),
	parent2(),
	penalties_()
{
	using namespace basic::options;
	load_penalty_table_from_file( option[OptionKeys::score::aa_repeat_energy_penalty_file]() );
}

/// @brief Copy constructor.
///
AARepeatEnergy::AARepeatEnergy( AARepeatEnergy const &src ) :
	parent1( methods::EnergyMethodCreatorOP( new AARepeatEnergyCreator ) ),
	parent2( src ),
	penalties_(src.penalties_)
{
}

/// @brief Default destructor.
///
AARepeatEnergy::~AARepeatEnergy() {}

/// @brief Clone: create a copy of this object, and return an owning pointer
/// to the copy.
EnergyMethodOP AARepeatEnergy::clone() const {
	return EnergyMethodOP( new AARepeatEnergy(*this) );
}

/// @brief AARepeatEnergy is context-independent and thus indicates that no context graphs need to be maintained by
/// class Energies.
void AARepeatEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const
{
	//Do nothing.
	return;
}

/// @brief AARepeatEnergy is version 1.0 right now.
///
core::Size AARepeatEnergy::version() const
{
	return 1; // Initial versioning
}

/// @brief Actually calculate the total energy
/// @details Called by the scoring machinery.
void AARepeatEnergy::finalize_total_energy( core::pose::Pose & pose, ScoreFunction const &, EnergyMap & totals ) const
{
	//Number of residues:
	core::Size const nres( pose.n_residue() );

	//Vector of residue owning pointers:
	utility::vector1< core::conformation::ResidueCOP > resvector;
	resvector.reserve(nres);

	//Populate the vector with const owning pointers to the residues:
	for ( core::Size ir=1; ir<=nres; ++ir ) {
		resvector.push_back( pose.residue(ir).get_self_ptr() );
	}

	totals[ aa_repeat ] += calculate_energy( resvector ); //Using the vector of residue owning pointers, calculate the repeat energy (unweighted) and set the aa_repeat_energy to this value.

	return;
}

/// @brief Calculate the total energy given a vector of const owning pointers to residues.
/// @details Called by finalize_total_energy().
core::Real AARepeatEnergy::calculate_energy( utility::vector1< core::conformation::ResidueCOP > const &resvect ) const
{
	core::Size const nres( resvect.size() );
	if ( nres==0 ) return 0.0;
	if ( nres==1 ) return penalties(1);

	//A counter (for number of residues in a row):
	core::Size counter(1);
	//An accumulator for the score:
	core::Real accumulator(0.0);

	for ( core::Size ir=2; ir<=nres; ++ir ) { //Loop through all residues, starting with the second.
		if ( //If this residue lacks a lower connection, or it's not connected to the previous residue normally, add the penalty to the accumulator and reset the counter.
				!resvect[ir]->has_lower_connect() || //Residue lacks a lower connection slot
				resvect[ir]->connection_incomplete( resvect[ir]->type().lower_connect_id() ) || //Incomplete lower connection.
				resvect[ir]->connected_residue_at_resconn( resvect[ir]->type().lower_connect_id() )!=ir-1 || //Not connected to the previous normally
				resvect[ir-1]->connected_residue_at_resconn( resvect[ir-1]->type().upper_connect_id() )!=ir //Previous not connected to this normally
				) {
			accumulator += penalties(counter); //Score the last stretch.
			counter=1; //Reset the counter.
			continue;
		}
		//If this residue is the same as the previous, increment the counter:
		if ( resvect[ir]->type().name3() == resvect[ir-1]->type().name3() ) ++counter;
		else { //If it's not, score the last stretch and reset the counter:
			accumulator += penalties(counter);
			counter=1;
		}
	}

	//The last stretch will not have been counted at this point.  We need to count it.
	accumulator += penalties(counter);

	return accumulator;
}

/// @brief Return a penalty for N residues in a row, from a lookup table.
/// @details The last entry in the lookup table is the penalty to return for more residues
/// than there are entries in the lookup table.  Returns 0 if nres is zero.
inline core::Real AARepeatEnergy::penalties( core::Size const nres ) const {
	if ( nres==0 || penalties_.empty() ) return 0.0;
	if ( nres > penalties_.size() ) {
		return penalties_[penalties_.size()];
	}
	return penalties_[nres];
}

/// @brief Read a penalty data file and load the penalties into the energy method.  Called once by the constructor.
/// @details Comment lines are ignored in the penalty data file.  The file should have one line that's a whitespace-separated
/// row of numbers.  The numbers represent the penalty for having a stretch of 1, 2, 3, 4, etc. of the same residue as a repeating
/// stretch.  The function looks in the working directory and in database/scoring/score_functions/aa_repeat_energy/ for the penalty
/// file.
void AARepeatEnergy::load_penalty_table_from_file( std::string const &filename ) {
	using namespace utility::io;

	std::string filename_formatted = filename;
	if ( utility::file::file_extension(filename_formatted)!="rpt_pen" ) filename_formatted+= ".rpt_pen";

	izstream infile;
	infile.open( filename_formatted );
	if ( !infile.good() ) {
		filename_formatted = "scoring/score_functions/aa_repeat_energy/" + utility::file::file_basename(filename_formatted) + ".rpt_pen";
		basic::database::open( infile, filename_formatted );
		runtime_assert_string_msg( infile.good(), "Error in AARepeatEnergy::load_penalty_table_from_file():  Unable to open .rpt_pen file for read!" );
	}

	if ( TR.Debug.visible() ) TR.Debug << "Reading repeat penalty energies from " << filename_formatted << std::endl;

	bool non_comment_line_present(false); //Have we found an actual data line in the file?
	std::string curline(""); //Buffer for current line.
	utility::vector1< std::string > lines; //Storing all lines

	//Read the file:
	while ( getline(infile, curline) ) {
		if ( curline.size() < 1 ) continue; //Ignore blank lines.
		//Find and process comments:
		std::string::size_type pound = curline.find('#', 0);
		if ( pound == std::string::npos ) {
			lines.push_back( curline );
		} else {
			lines.push_back(curline.substr(0, pound));
		}
	}
	infile.close();

	//Parse the lines from the file:
	for ( core::Size i=1, imax=lines.size(); i<=imax; ++i ) {
		if ( non_comment_line_present ) {
			runtime_assert_string_msg( !parse_line( lines[i], penalties_ ), "Error in AARepeatEnergy::load_penalty_table_from_file():  More than one line containing data was found in the .rpt_pen file." );
		} else { //If we have not yet found the data
			if ( parse_line( lines[i], penalties_) ) non_comment_line_present=true;
		}
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Done reading penalties from " << filename_formatted << std::endl;
		TR.Debug.flush();
	}
	return;
}

/// @brief Parse a series of floats from a string and store them in the penalties vector.  Return true if successful and false otherwise.
/// @details Only modifies the penalties vector if successful.
bool AARepeatEnergy::parse_line( std::string const &line, utility::vector1< core::Real > &penalties ) const {
	std::istringstream l(line);
	core::Real curval;
	utility::vector1 < core::Real > outvect;
	bool good_parse(false);

	while ( !l.eof() ) {
		l >> curval;
		if ( !l.fail() ) {
			good_parse=true;
			outvect.push_back(curval);
		} else {
			good_parse=false;
			break;
		}
	}

	if ( good_parse ) penalties=outvect;
	return good_parse;
}

} // aa_repeat
} // scoring
} // core
