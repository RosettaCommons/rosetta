// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/bin_transitions/BinTransitionCalculator.cc
/// @brief  Class member functions for BinTransitionCalculator class.
/// @details This class loads data associated with transitions from one mainchain torsion bin to another (e.g. ABEGO bins, OO-ABBA bins, etc.)
/// and provides methods for generating bin sequences randomly, perturbing sequences, etc.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit Headers
#include <core/scoring/bin_transitions/BinTransitionCalculator.hh>
#include <core/scoring/bin_transitions/BinTransitionData.hh>
#include <core/scoring/bin_transitions/BinTransitionData.fwd.hh>

// File I/O
#include <basic/database/open.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>

//Random
#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>

// Other Headers
#include <basic/Tracer.hh>

namespace core {
namespace scoring {
namespace bin_transitions {

static THREAD_LOCAL basic::Tracer TR( "core.scoring.bin_transitions.BinTransitionCalculator" );

/// @brief Default constructor for BinTransitionCalculator
///
BinTransitionCalculator::BinTransitionCalculator(): //TODO -- initialize variables here:
	bin_params_loaded_(false),
	bin_params_file_(""),
	bin_transition_data_()
{}

/// @brief Copy constructor for BinTransitionCalculator
///
BinTransitionCalculator::BinTransitionCalculator( BinTransitionCalculator const &src ): //TODO -- copy variables here:
	utility::pointer::ReferenceCount(),
	bin_params_loaded_(src.bin_params_loaded_),
	bin_params_file_(src.bin_params_file_),
	bin_transition_data_()
{
	bin_transition_data_.clear();
	for ( core::Size i=1, imax=src.bin_transition_data_.size(); i<=imax; ++i ) bin_transition_data_.push_back( src.bin_transition_data_[i]->clone() ); //Clone all of the BinTransitionData objects in src.
}

/// @brief Default destructor for BinTransitionCalculator
///
BinTransitionCalculator::~BinTransitionCalculator() {}

/// @brief Clone operation for BinTransitionCalculator.
/// @details Returns an owning pointer to a copy of this object.
BinTransitionCalculatorOP BinTransitionCalculator::clone() const
{ return BinTransitionCalculatorOP( new BinTransitionCalculator( *this ) ); }

/// @brief Set the bin params file, and load all the bin params and transition probabilities.
///
void BinTransitionCalculator::load_bin_params( std::string const &filename ) {
	using namespace utility::io;

	//Check that this object has not already been initialized.
	runtime_assert_string_msg( !bin_params_loaded(), "In core::scoring::bin_transitions::BinTransitionCalculator::load_bin_params(): A bin params file was already loaded by this object!  BinTransitionCalculator objects are intended to be single-use (as far as bin params are concerned)." );
	bin_params_loaded_=true;

	//Format the filename:
	bin_params_file_ = filename;
	if ( utility::file::file_extension(bin_params_file_)!="bin_params" ) bin_params_file_ += ".bin_params";
	izstream infile;
	infile.open( bin_params_file_ );
	if ( !infile.good() ) {
		bin_params_file_ = "protocol_data/generalizedKIC/bin_params/" + utility::file::file_basename(bin_params_file_) + ".bin_params";
		basic::database::open( infile, bin_params_file_ );
		runtime_assert_string_msg( infile.good(),
			"In core::scoring::bin_transitions::BinTransitionCalculator::load_bin_params(): Unable to open .bin_params file for read!  The file " + filename + " does not exist or could not be read!" );
	}
	if ( TR.visible() ) TR << "Opened file " << bin_params_file_ << " for read." << std::endl;

	std::string curline;
	utility::vector1< std::string > lines;

	while ( getline(infile, curline) ) {

		if ( TR.Debug.visible() ) TR.Debug << curline << std::endl;

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

	runtime_assert_string_msg( lines.size() > 0, "In core::scoring::bin_transitions::BinTransitionCalculator::load_bin_params(): loaded no lines from the specified bin_params file!" );

	//Dump out the lines to parse:
	if ( TR.Debug.visible() ) { //CONVERT TO DEBUG MODE LATER!
		TR.Debug << std::endl << "Lines to parse:" << std::endl;
		for ( core::Size i=1, imax=lines.size(); i<=imax; ++i ) TR.Debug << lines[i] << std::endl;
		TR.Debug << std::endl;
	}

	//Parse the lines in chunks, where each chunk is between a pair of lines that read "BEGIN" and "END":
	core::Size firstline(0), lastline(0);
	while ( has_another_matrix(lines, firstline, lastline) ) {
		BinTransitionDataOP curdata( add_bin_transition_data() );
		parse_mainchain_torsions( lines, firstline, lastline, curdata ); //Get the number of mainchain torsions in the ith and i+1st residues for the probability matrix defined in this block.
		initialize_matrix( lines, firstline, lastline, curdata ); //Get the dimensions of the probability matrix, and initialize the storage space in memory (in the BinTransitionData object).
		set_up_bins( lines, firstline, lastline, curdata); //Get the angle ranges for each bin (ith and i+1st residue), as well as the bin names.
		populate_matrix( lines, firstline, lastline, curdata ); //Actually read in the transition probabilities and populate the transition matrix.
		store_properties( lines, firstline, lastline, curdata ); //Store the properties associated with this transition matrix, for the ith and i+1st residues.
		store_residue_identities( lines, firstline, lastline, curdata); //Store the allowed residue identities for the ith and i+1st residues.
		curdata->finalize(); //Do the final post-load calculations.
	}

	if ( TR.Debug.visible() ) TR.Debug << summarize_stored_data( true );

	if ( TR.visible() ) {
		TR << "Finished read of " << bin_params_file_ << "." << std::endl;
		TR.flush();
	}
	return;
} //load_bin_params()

/// @brief Given a bin name and a residue, find a BinTransitionsData object describing that residue,
/// and the index of the bin within that object.
/// @details data_index and bin_index are outputs.  Both are set to 0 if the search fails.  Everything
/// else is const input.  If bin_name is set to "", then the bin index that is returned is the index
/// corresponding to the residue's mainchain torsion vector.
void BinTransitionCalculator::find_data_and_bin(
	std::string const &bin_name,
	core::conformation::Residue const &res,
	core::Size &data_index,
	core::Size &bin_index,
	bool const use_iplus1
) const {
	//Initialize data_index and bin_index
	data_index=0;
	bin_index=0;

	if ( n_bin_transition_data()==0 ) return; //If there are no BinTransitionData objects, we've already failed.

	bool const use_tors_vect( bin_name=="" ); //If the bin_name is set to "", use the residue's mainchain torsion vector to fish out the bin.

	//Loop through all BinTransitionData objects:
	for ( core::Size i=1, imax=n_bin_transition_data(); i<=imax; ++i ) {
		//Initialize data_index and bin_index
		data_index=0;
		bin_index=0;

		if ( !use_iplus1 ) {
			if ( !bin_transition_data(i)->criteria_match_i( res ) ) continue;
			data_index=i;
			bin_index = ( use_tors_vect ?  bin_transition_data(i)->which_bin_i( res.mainchain_torsions() ) : bin_transition_data(i)->binname_index_from_string_i(bin_name) );
			if ( bin_index==0 ) continue;
		} else {
			if ( !bin_transition_data(i)->criteria_match_iplus1( res ) ) continue;
			data_index=i;
			bin_index = ( use_tors_vect ?  bin_transition_data(i)->which_bin_iplus1( res.mainchain_torsions() ) : bin_transition_data(i)->binname_index_from_string_iplus1(bin_name) );
			if ( bin_index==0 ) continue;
		}
		//If we get to here, then we've found a BinTransitionData object and the bin index,
		//and data_index and bin_index have been set.  So we're done.
		if ( data_index!=0 && bin_index!=0 /*Redundant check*/ ) return;
	}

	//If we get to here, then we have failed and should return 0 for data_index and bin_index.
	data_index=0; bin_index=0;
	return;
} //find_data_and_bin

/// @brief Given residues at positions i and iplus1, find a bin transition probability data object that
/// describes the pair.
/// @details Inputs are res_i and res_iplus1.  Outputs are data_i and data_iplus1 (BinTransitionData
/// indices).  Function returns true if data are found successfully, false if no BinTransitionData object
/// could be found describing the residues in question.
bool BinTransitionCalculator::find_data(
	core::conformation::Residue const &res_i,
	core::conformation::Residue const &res_iplus1,
	core::Size &data_index
) const {
	if ( n_bin_transition_data()==0 ) return 1; //Fail immediately if there are no defined BinTranstionData objects.

	for ( core::Size i=1, imax=n_bin_transition_data(); i<=imax; ++i ) {
		if ( !bin_transition_data(i)->criteria_match_i( res_i ) ) continue;
		if ( !bin_transition_data(i)->criteria_match_iplus1( res_iplus1 ) ) continue;
		//At this point, we've found a match.
		data_index=i;
		return true;
	}

	//At this point, we've failed to find a match.
	data_index=0;
	return false;
}

/// @brief Given two residues (rsd_i and rsd_iplus1 at positions i and i+1, respectively), give me the
/// probability of seeing rsd_iplus1 in its bin given that rsd_i is in its bin.
/// @details The probability value is set to a number from 0 to 1.  Inputs are rsd_i and rsd_iplus1.
/// Function returns true if data are found successfully, false if no BinTransitionData object
/// could be found describing the residues in question.
bool BinTransitionCalculator::p_iplus1_given_i(
	core::conformation::Residue const &res_i,
	core::conformation::Residue const &res_iplus1,
	core::Real &probability
) const {
	if ( !bin_params_loaded() ) {
		utility_exit_with_message(
			"In core::scoring::bin_transitions::BinTransitionCalculator::p_iplus1_given_i(): Bin transition parameters have not yet been loaded!"
		);
	}

	core::Size data_index(0);
	core::Size bin_i(0);
	core::Size bin_iplus1(0);

	if ( !find_data( res_i, res_iplus1, data_index ) || data_index==0 ) {
		if ( TR.Warning.visible() ) {
			TR.Warning << "No bin transition probability data could be found for residues " << res_i.seqpos() << " and " << res_iplus1.seqpos() << "." << std::endl;
			TR.Warning.flush();
		}
		return false;
	}

	bin_i = bin_transition_data(data_index)->which_bin_i( res_i.mainchain_torsions() );
	bin_iplus1 = bin_transition_data(data_index)->which_bin_iplus1( res_iplus1.mainchain_torsions() );

	probability = bin_transition_data(data_index)->probability_matrix( bin_i, bin_iplus1 ) / bin_transition_data(data_index)->binsums_i( bin_i );

	return true;
}

/// @brief Given two residues (rsd_i and rsd_iplus1 at positions i and i+1, respectively), give me the
/// probability of seeing rsd_i in its bin given that rsd_iplus1 is in its bin.
/// @details The probability value is set to a number from 0 to 1.  Inputs are rsd_i and rsd_iplus1.
/// Function returns true if data are found successfully, false if no BinTransitionData object
/// could be found describing the residues in question.
bool BinTransitionCalculator::p_i_given_iplus1(
	core::conformation::Residue const &res_i,
	core::conformation::Residue const &res_iplus1,
	core::Real &probability
) const {
	if ( !bin_params_loaded() ) {
		utility_exit_with_message(
			"In core::scoring::bin_transitions::BinTransitionCalculator::p_iplus1_given_i(): Bin transition parameters have not yet been loaded!"
		);
	}

	core::Size data_index(0);
	core::Size bin_i(0);
	core::Size bin_iplus1(0);

	if ( !find_data( res_i, res_iplus1, data_index ) || data_index==0 ) {
		if ( TR.Warning.visible() ) {
			TR.Warning << "No bin transition probability data could be found for residues " << res_i.seqpos() << " and " << res_iplus1.seqpos() << "." << std::endl;
			TR.Warning.flush();
		}
		return false;
	}

	bin_i = bin_transition_data(data_index)->which_bin_i( res_i.mainchain_torsions() );
	bin_iplus1 = bin_transition_data(data_index)->which_bin_iplus1( res_iplus1.mainchain_torsions() );

	probability = bin_transition_data(data_index)->probability_matrix( bin_i, bin_iplus1 ) / bin_transition_data(data_index)->binsums_iplus1( bin_iplus1 );

	return true;
}

/// @brief Is the given residue in the given bin?
/// @details For the bin definitions, this uses the first BinTransitionsData object that it finds where residue i matches the properties of the given
/// residue.  Checks the i+1 definitions if nothing is found for the i transitions.  Fails if no bin definitions for the given residue type are found.
bool BinTransitionCalculator::is_in_bin (
	core::conformation::Residue const &res,
	std::string const &bin_name
) const {
	if ( bin_name=="" ) utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionCalculator::is_in_bin(): Received an empty string for the bin name." );

	core::Size data_index(0);
	core::Size bin_index(0);
	bool use_iplus1(false);

	find_data_and_bin( bin_name, res, data_index, bin_index, use_iplus1 );
	if ( data_index==0 ) {
		use_iplus1=true;
		find_data_and_bin( bin_name, res, data_index, bin_index, use_iplus1 );
	}
	if ( data_index==0 ) utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionCalculator::is_in_bin(): Could not find suitable bin definitions for the given residue type." );

	if ( !use_iplus1 ) {
		return bin_transition_data(data_index)->in_bin_i( bin_index, res );
	} else {
		return bin_transition_data(data_index)->in_bin_iplus1( bin_index, res );
	}

	return true;
}

/// @brief Initialize a string of residues to a bunch of random bins, based on bin transition probabilities; then draw random mainchain torsion angles from those bins.
/// @details Takes a const conformation and a const list of residue indices as input; the conformation is just for checking residues types, numbers of mainchain torsions, etc.
/// The residue indices must be in order, defining a contiguous chain (running backwards or forwards).  Output is the mainchain_torsions vector of vectors (reset and
/// replaced by this operation).  The distribution WITHIN the bin depends on the BinTransitionData object and what was specified in the bin_params file.  Default is uniform
/// within each bin, though Ramachandran-biased distributions are also permitted for alpha-amino acids.
void BinTransitionCalculator::random_mainchain_torsions_from_bins(
	core::conformation::Conformation const &conformation,
	utility::vector1 <core::Size> const &res_indices,
	utility::vector1 < utility::vector1 < core::Real > > &mainchain_torsions
) const {
	//Get the number of residues that we'll be setting:
	core::Size const nres(res_indices.size());
	runtime_assert_string_msg(nres > 0, "In core::scoring::bin_transitions::BinTransitionCalculator::random_mainchain_torsions_from_bins(): We need at least one residue to operate on!");

	//Initialize the mainchain_torsions vector:
	mainchain_torsions.clear();
	for ( core::Size ir=1; ir<=nres; ++ir ) {
		utility::vector1 < core::Real > newvect;
		newvect.resize( conformation.residue(res_indices[ir]).mainchain_torsions().size(), 0.0);
		mainchain_torsions.push_back(newvect);
	}

	bool iplus1(false); //Are we drawing from the bins for the i+1st residue, or the ith?
	core::Size data_index(0); //Which BinTransitionData object has criteria matching this residue?

	//Step 1: draw a random bin for the first residue, based on the bin probability distribution.
	core::Size bin_index_1( random_bin(conformation.residue(res_indices[1]), iplus1, data_index )); //iplus1 and data_index are set by this function

	//Step 2: draw random mainchain torsions from that bin.
	random_mainchain_torsions_from_bin( bin_index_1, iplus1, data_index, mainchain_torsions[1] ); //iplus1 and data_index are used by this function

	if ( nres > 1 ) {
		for ( core::Size ir=2; ir<=nres; ++ir ) {
			//Step 3: draw a random bin for the second residue, based on the probability of transitioning from the first.
			core::Size bin_index( random_bin_based_on_previous( conformation.residue(res_indices[ir]), conformation.residue(res_indices[ir-1]), mainchain_torsions[ir-1], data_index ) ); //data_index is reset by this function (recycling the variable).

			//Step 4: draw random mainchain torsions for the second residue's bin.
			random_mainchain_torsions_from_bin( bin_index, true, data_index, mainchain_torsions[ir] );

		} //Step 5: repeat 3 and 4 for all subsequent residues.
	}

	return;
} //random_mainchain_torsions_from_bins

/// @brief Draw random mainchain torsion values for a set of residues, given a bin from which the values should be drawn.
/// @details Takes a bin name, a const conformation, and a const list of residue indices as input; the conformation is just for checking residues types, numbers of mainchain
/// torsions, etc.  Output is the mainchain_torsions vector of vectors (reset and replaced by this operation).  The distribution WITHIN the bin depends on the BinTransitionData
/// object and what was specified in the bin_params file.  Default is uniform within each bin, though Ramachandran-biased distributions are also permitted for alpha-amino acids.
/// Note that this function uses bins for residue i, and only checks i+1 if no suitable data are found for i.
void BinTransitionCalculator::random_mainchain_torsions_from_bin(
	std::string const &bin_name,
	core::conformation::Conformation const &conformation,
	utility::vector1 <core::Size> const &res_indices,
	utility::vector1 < utility::vector1 < core::Real > > &mainchain_torsions
) const {
	//Initial checks:
	if ( bin_name=="" ) {
		utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionCalculator::random_mainchain_torsions_from_bin(): No bin name was specified!" );
	}

	//Get the number of residues that we'll be setting:
	core::Size const nres(res_indices.size());
	if ( nres < 1 ) utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionCalculator::random_mainchain_torsions_from_bin(): We need at least one residue to operate on!");

	//Storage slots for the index of the BinTransitionData object that we'll use, and for the bin index:
	core::Size data_index(0);
	core::Size bin_index(0);

	//Generate the mainchain_torsions vector:
	mainchain_torsions.clear();
	for ( core::Size ir=1; ir<=nres; ++ir ) {
		utility::vector1 < core::Real > newvect;
		newvect.resize( conformation.residue(res_indices[ir]).mainchain_torsions().size(), 0.0);

		//Find the data index and bin index:
		data_index=0;
		bin_index=0;
		bool use_iplus1=false;
		find_data_and_bin( bin_name, conformation.residue(res_indices[ir]), data_index, bin_index, use_iplus1 );
		if ( data_index==0 ) {
			use_iplus1=true;
			find_data_and_bin( bin_name, conformation.residue(res_indices[ir]), data_index, bin_index, use_iplus1 ); //Use the i+1st bin if nothing is found for bin i.
		}
		if ( data_index==0 ) utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionCalculator::random_mainchain_torsions_from_bin(): No BinTransitionData object could be found for the specified bin and residue type!" );

		//Draw random mainchain torsions from that bin:
		random_mainchain_torsions_from_bin( bin_index, use_iplus1, data_index, newvect  );

		//Add the new mainchain torsions to the mainchain torsion vector:
		mainchain_torsions.push_back(newvect);
	}

	//And we're done!
	return;
} //random_mainchain_torsions_from_bin

/// @brief Randomly pick mainchain torsions for a residue based on the torsion bins of its i+1 and i-1 neighbours.
/// @details Takes a const conformatoin, a const residue index, and a boolean valueas input.  The conformation is for checking residue
/// types, numbers of mainchain torsions, etc.  The boolean determines whether this residue should be allowed to stay in its own bin
/// or be required to switch to another bin.  Output is the mainchain_torsions vector of Reals, with one entry for each mainchain torsion
/// of the residue.  Note that random mainchain torsions are drawn using the sub-bins for the case where the current residue is the ith
/// residue.
void BinTransitionCalculator::random_mainchain_torsions_using_adjacent_bins(
	core::conformation::Conformation const &conformation,
	core::Size const res_index,
	bool const must_switch_bins,
	utility::vector1 < core::Real > &mainchain_torsions
) const {
	core::Size data_index_i(0); //Which BinTransitionData object has criteria matching this residue for the ith position?
	core::Size data_index_iplus1(0); //Which BinTransitionData object has criteria matching this residue for the i+1st position?

	//Step 1: Draw a random bin for this residue, based on the bins of the previous and next residues.
	core::Size bin_index( random_bin_based_on_prev_and_next( conformation, res_index, must_switch_bins, data_index_i, data_index_iplus1 ) );

	//Step 2: Draw random mainchain torsions from that bin.
	random_mainchain_torsions_from_bin( bin_index, false, data_index_i, mainchain_torsions );

	return;
} //random_mainchain_torsions_using_adjacent_bins

/// @brief Randomly pick a bin, given an input residue.
/// @details The residue's properties and name are used to find the first BinTransitionData object matching the properties and name.
/// A random bin is then picked based on the frequency with which that bin is observed for that class of residues.  If the bin ends
/// up being from the i+1st residue, iplus1 is set to "true"; otherwise, it remains "false".  The data_index variable is set by this
/// function to the index of the BinTransitionData object that provides the transition probability data.
core::Size BinTransitionCalculator::random_bin(
	core::conformation::Residue const &res,
	bool &iplus1,
	core::Size &data_index
) const {
	core::Size const ndata( n_bin_transition_data() );
	if ( !bin_params_loaded() || ndata < 1 ) utility_exit_with_message( "In BinTransitionCalculator::random_bin(): no transition probability matrices have been defined!\n" );

	for ( core::Size i=1; i<=ndata; ++i ) { //Loop through all BinTransitionData objects (i.e. all transition probability matrices) and use the first with criteria matching this residue.
		if ( bin_transition_data(i)->criteria_match_i( res ) ) {
			iplus1=false;
			data_index=i;
			return bin_transition_data(i)->random_bin_i();
		} else if ( bin_transition_data(i)->criteria_match_iplus1( res ) ) {
			iplus1=true;
			data_index=i;
			return bin_transition_data(i)->random_bin_iplus1();
		}
	}

	utility_exit_with_message( "In BinTransitionCalculator::random_bin(): No transition probability data were found for a " + res.name3() + " residue!\n" );

	return 0;
} //random_bin

/// @brief Randomly pick mainchain torsion values from within a bin
/// @details The residue's properties and name are used to find the first BinTransitionData object matching the properties and name.
/// Mainchain torsion values are then picked randomly from within the bin.  The distribution WITHIN the bin depends on the BinTransitionData
/// object and what was specified in the bin_params file.  Default is uniform within each bin, though Ramachandran-biased distributions are
/// also permitted for alpha-amino acids.  If iplus1 is true, we draw from the bins for the i+1st residue; otherwise, we draw from the bins
/// for the ith residue.  The data_index value tells this function which BinTransitionData object to use.
void BinTransitionCalculator::random_mainchain_torsions_from_bin(
	//core::conformation::Residue const &res,
	core::Size const bin_index,
	bool const iplus1,
	core::Size const data_index,
	utility::vector1 < core::Real > &mainchain_torsions //output
) const {
	if ( data_index < 1 || data_index > n_bin_transition_data() ) utility_exit_with_message( "In BinTransitionCalculator::random_mainchain_torsions_from_bin(): the index of the BinTransitionData object is out of range.\n" );
	if ( bin_index < 1 || (!iplus1 && bin_index > bin_transition_data(data_index)->n_bins_i()) || (iplus1 && bin_index > bin_transition_data(data_index)->n_bins_iplus1()) ) {
		utility_exit_with_message( "In BinTransitionCalculator::random_mainchain_torsions_from_bin(): the index of the bin is out of range.\n" );
	}

	BTSB_SUBBIN_TYPE bintype = (!iplus1 ? bin_transition_data(data_index)->subbin_type_i() : bin_transition_data(data_index)->subbin_type_iplus1());

	switch( bintype ) {
	case BTSB_L_AA:
	case BTSB_D_AA:
	case BTSB_GLY :
		biased_mainchain_torsions( /*res,*/ bin_index, iplus1, data_index, mainchain_torsions );
		break;
	default :
		uniform_mainchain_torsions( /*res,*/ bin_index, iplus1, data_index, mainchain_torsions);
		break;
	}

	return;
} // random_mainchain_torsions_from_bin

/// @brief Given the current residue and the previous residue, as well as the mainchain torsion values for the previous residue,
/// pick a bin for the current residue based on bin transition probabilities.
/// @details This function takes thisres (the current residue) and prevres (the previous residue) as inputs, using them for their
/// properties and identities to find a BinTransitionData object that matches the criteria at the i and i+1st residues.  It then
/// sets data_index to the index of that BinTransitionData object before returning the index of a bin chosen for the current (i+1st)
/// residue at random, based on the relative probabilities of bins given the previous (ith) residue.
core::Size BinTransitionCalculator::random_bin_based_on_previous (
	core::conformation::Residue const &thisres,
	core::conformation::Residue const &prevres,
	utility::vector1 <core::Real> const &prev_torsions,
	core::Size &data_index
) const {
	//First we must find a BinTransitionData object describing this pair of residues:
	core::Size const ndata( n_bin_transition_data() );
	data_index=0; //This will be set to the index if the BinTransitionData object.  We set it to zero here to use to test for whether the object was found.
	if ( !bin_params_loaded() || ndata < 1 ) {
		utility_exit_with_message(
			"In BinTransitionCalculator::random_bin_based_on_previous(): No transition probability matrices have been defined!\n" );
	}

	//Check that the residues are normally bonded to one another:
	bool normal_bond_found(are_normally_bonded( thisres, prevres ));
	if ( !normal_bond_found ) {
		utility_exit_with_message(
			"In BinTransitionCalculator::random_bin_based_on_previous(): The two residues passed to this function are not properly connected." );
	}

	for ( core::Size i=1; i<=ndata; ++i ) { //Loop through all BinTransitionData objects (i.e. all transition probability matrices) and use the first with criteria matching this residue.
		if ( bin_transition_data(i)->criteria_match_i( prevres ) && bin_transition_data(i)->criteria_match_iplus1(thisres) ) { //If we've found a BinTransitionData object describing this pair of residues
			data_index=1;
			break;
		}
	}

	//Check that a BinTransitionData object was found, and create an owning pointer to it.
	if ( data_index==0 ) {
		utility_exit_with_message(
			"In BinTransitionCalculator::random_bin_based_on_previous(): No transition probability matrix was found describing a " + prevres.name3() + "->" + thisres.name3() + " transition.\n" );
	}
	BinTransitionDataCOP curdata( bin_transition_data(data_index) );

	//Next we must figure out what bin the previous residue was in, according to the current BinTransitionData object:
	core::Size prevbin( curdata->which_bin_i( prev_torsions ) );

	return curdata->random_bin_iplus1_given_i( prevbin );
} //random_bin_based_on_previous

/// @brief Given the current residue and its conformation, pick a bin for the current residue based on the bins of the previous and
/// next residues, and the bin transition probabilities.
/// @details Bin transition probabilities are assumed to be independent.  That is, P(ABC) = P(B | A & C) = P (B | A ) * P(B | C).
/// This function also finds and returns the indices of the BinTransitionData objects that describe the AB transition (data_index_iplus1)
/// and the BC transition (data_index_i).  So conformation and res_index are const inputs, and data_index_i and data_index_iplus1 are outputs.
/// If must_switch_bins is true, then the current bin for the current residue is never chosen; otherwise, it is in the pool to be picked.
core::Size BinTransitionCalculator::random_bin_based_on_prev_and_next(
	core::conformation::Conformation const &conformation,
	core::Size const res_index,
	bool const must_switch_bins,
	core::Size &data_index_i,
	core::Size &data_index_iplus1
) const {
	//Check that res_index is reasonable:
	if ( res_index==0 || res_index>conformation.size() ) {
		utility_exit_with_message(
			"In BinTransitionCalculator::random_bin_based_on_prev_and_next(): The residue index passed to this function is not within the current conformation's range of residues.");
	}

	//The current residue:
	core::conformation::Residue const &thisres( conformation.residue(res_index) );

	//Check that the residue has upper and lower connections.
	if ( !thisres.has_lower_connect() || !thisres.has_upper_connect() ) {
		utility_exit_with_message(
			"In BinTransitionCalculator::random_bin_based_on_prev_and_next(): The residue passed to this function lacks lower or upper connections.");
	}

	//Get the residues that this one is bound to:
	core::Size const prev_index( thisres.connected_residue_at_resconn( thisres.type().lower_connect_id() ) );
	core::Size const next_index( thisres.connected_residue_at_resconn( thisres.type().upper_connect_id() ) );
	if ( prev_index==0 || next_index==0 ) {
		utility_exit_with_message(
			"In BinTransitionCalculator::random_bin_based_on_prev_and_next(): The residue passed to this function lacks residues that it's connected to at its lower or upper connections.");
	}
	core::conformation::Residue const &prevres( conformation.residue(prev_index) );
	core::conformation::Residue const &nextres( conformation.residue(next_index) );

	//Check that the residues are normally bonded to one another:
	bool normal_bond_found( are_normally_bonded( thisres, prevres ) && are_normally_bonded( nextres, thisres ) );
	if ( !normal_bond_found ) {
		utility_exit_with_message(
			"In BinTransitionCalculator::random_bin_based_on_prev_and_next(): The residue passed to this function is not properly connected to the previous or next residue." );
	}

	//Get the number of BinTransitionData objects:
	core::Size const ndata( n_bin_transition_data() );

	//Get the BinTransitionData indices describing the prevres->thisres and thisres->nextres transitions.
	data_index_i=0;
	data_index_iplus1=0;
	for ( core::Size i=1; i<=ndata; ++i ) {
		if ( data_index_i==0 && bin_transition_data(i)->criteria_match_i( thisres ) && bin_transition_data(i)->criteria_match_iplus1(nextres) ) { //We've found a this->next transition
			data_index_i=i;
		}
		if ( data_index_iplus1==0 && bin_transition_data(i)->criteria_match_i( prevres ) && bin_transition_data(i)->criteria_match_iplus1(thisres) ) { //We've found a prev->this transition
			data_index_iplus1=i;
		}
		if ( data_index_i != 0 && data_index_iplus1 != 0 ) break;
	}

	if ( data_index_i == 0 ) {
		utility_exit_with_message(
			"In BinTransitionCalculator::random_bin_based_on_prev_and_next(): Could not find a suitable transition probability matrix (BinTransitionData object) describing the transition from the current residue to the next residue.");
	}
	if ( data_index_iplus1 == 0 ) {
		utility_exit_with_message(
			"In BinTransitionCalculator::random_bin_based_on_prev_and_next(): Could not find a suitable transition probability matrix (BinTransitionData object) describing the transition from the previous residue to the current residue.");
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Number of bins for this residue based on this->next BinTransitionData object: " << bin_transition_data(data_index_i)->n_bins_i() << std::endl;
		TR.Debug << "Number of bins for this residue based on prev->this BinTransitionData object: " << bin_transition_data(data_index_iplus1)->n_bins_iplus1() << std::endl;
	}
	if ( bin_transition_data(data_index_i)->n_bins_i() != bin_transition_data(data_index_iplus1)->n_bins_iplus1() ) {
		utility_exit_with_message( "In BinTransitionCalculator::random_bin_based_on_prev_and_next(): The number of bins for this residue when considering the prev->this transition does not match the number of bins when considering the this->next transition.  Are bin definitions consistent in the bin_params file?" );
	} else {
		if ( TR.Debug.visible() ) TR.Debug << "Bin counts match.  Assuming bin identities are identical, too." << std::endl;
	}

	//Figure out which bin the residue is currently in.
	core::Size const curbin_index_i( bin_transition_data(data_index_i)->which_bin_i( thisres.mainchain_torsions() ) );
	core::Size const curbin_index_iplus1( bin_transition_data(data_index_iplus1)->which_bin_iplus1( thisres.mainchain_torsions() ) );
	if ( curbin_index_i != curbin_index_iplus1 ) {
		utility_exit_with_message(
			"In BinTransitionCalculator::random_bin_based_on_prev_and_next(): The bin indices for this residue when it is at the i or i+1st position don't match!  Are bin definitions consistent in the bin_params file?");
	}

	//Figure out which bin the previous and next residues are currently in.
	core::Size const prevbin_index( bin_transition_data(data_index_iplus1)->which_bin_i( prevres.mainchain_torsions() ) );
	core::Size const nextbin_index( bin_transition_data(data_index_i)->which_bin_iplus1( nextres.mainchain_torsions() ) );

	//Construct the probability vector.  Each entry corresponds to a bin at this position.
	utility::vector1 < core::Real > probvect;
	probvect.resize( bin_transition_data(data_index_i)->n_bins_i() , 0.0 );
	core::Real probvect_sum(0.0); //Accumulator for normalization
	for ( core::Size i=1, imax=probvect.size(); i<=imax; ++i ) { //Loop through and populate the probability vector:
		if ( must_switch_bins && i==curbin_index_i ) probvect[i]=0; //Set the probability of the current bin to zero if must_switch_bins is true.
		else { //Calculate the probability of the current bin as the product of probabilities:
			probvect[i] = bin_transition_data(data_index_iplus1)->probability_matrix(prevbin_index, i) / bin_transition_data(data_index_iplus1)->binsums_i(prevbin_index)
				* bin_transition_data(data_index_i)->probability_matrix(i, nextbin_index) / bin_transition_data(data_index_i)->binsums_iplus1(nextbin_index);
			probvect_sum += probvect[i];
		}
	}

	//Normalize the probability vector and construct the cdf:
	utility::vector1 < core::Real > probvect_cdf;
	probvect_cdf.resize( probvect.size(), 0.0 );
	for ( core::Size i=1, imax=probvect.size(); i<=imax; ++i ) {
		probvect[i] = probvect[i]/probvect_sum;
		if ( i<imax ) probvect_cdf[i+1]=probvect[i]+probvect_cdf[i]; //Note that Rosetta uses weird cdfs, where the P(i) = cdf(i+1)-cdf(i).
	}

	//Now use the cdf to pick a bin:
	core::Size const random_bin (numeric::random::pick_random_index_from_cdf( probvect_cdf, numeric::random::rg() ) );

	if ( TR.visible() ) TR.flush();
	if ( TR.Debug.visible() ) TR.Debug.flush();
	return random_bin;
} //random_bin_based_on_prev_and_next

/// @brief Given a bin, generate phi and psi values from within the bin, biased by the sub-bin cumulative probability distribution function.
/// @details bin_index, iplus1, and data_index are inputs, mainchain_torsions vector is cleared and set by this function (i.e. it's output).
void BinTransitionCalculator::biased_mainchain_torsions(
	//core::conformation::Residue const &res,
	core::Size const bin_index,
	bool const iplus1,
	core::Size const data_index,
	utility::vector1 <core::Real> &mainchain_torsions
) const {
	core::Size const ndata( n_bin_transition_data() );
	if ( !bin_params_loaded() || ndata < 1 ) utility_exit_with_message( "In BinTransitionCalculator::biased_mainchain_torsions(): no transition probability matrices have been defined!\n" );
	if ( data_index < 1 || data_index > ndata ) utility_exit_with_message( "In BinTransitionCalculator::biased_mainchain_torsions(): the index of the BinTransitionData object is out of range.\n" );
	if ( bin_index < 1 || (!iplus1 && bin_index > bin_transition_data(data_index)->n_bins_i()) || (iplus1 && bin_index > bin_transition_data(data_index)->n_bins_iplus1()) ) {
		utility_exit_with_message( "In BinTransitionCalculator::biased_mainchain_torsions(): the index of the bin is out of range.\n" );
	}

	BinTransitionDataCOP thisres_data( bin_transition_data_[data_index] ); //Pointer to the transition probability data for this residue.
	core::Size const ntors( (iplus1 ? thisres_data->n_mainchain_torsions_iplus1() : thisres_data->n_mainchain_torsions_i()) ); //Number of mainchain torsions.
	assert(ntors>0); //Should already be guaranteed.

	//Initialize mainchain torsions vector:
	mainchain_torsions.clear();
	mainchain_torsions.resize(ntors, 0.0);

	//Randomly pick a sub-bin based on the probability distribution:
	core::Size subbin (numeric::random::pick_random_index_from_cdf( ( iplus1 ? thisres_data->subbin_cdf_iplus1(bin_index) : thisres_data->subbin_cdf_i(bin_index) ), numeric::random::rg() ));
	for ( core::Size i=1; i<=ntors; ++i ) {
		//Get the ranges of this sub-bin for this mainchain torsion
		core::Real start( (iplus1 ? thisres_data->subbin_start_iplus1(bin_index, subbin, i) : thisres_data->subbin_start_i(bin_index, subbin, i)) );
		core::Real end( (iplus1 ? thisres_data->subbin_end_iplus1(bin_index, subbin, i) : thisres_data->subbin_end_i(bin_index, subbin, i)) );
		if ( start > end ) end +=360.0; //If this is a wrap-around bin, we'll pick numbers that go past 180, then re-adjust so that the result is between -180 and 180.
		mainchain_torsions[i] = start + ( end-start )*numeric::random::rg().uniform(); //Pick a random angle within this sub-bin
		if ( mainchain_torsions[i] > 180.0 ) mainchain_torsions[i] -= 360.0; //Bring it back to the -180, 180 range.
		assert(mainchain_torsions[i] <= 180.0 && mainchain_torsions[i] >= -180.0); //Should be true.
	}

	return;
} //biased_mainchain_torsions

/// @brief Given a bin, generate phi and psi values randomly from within the bin, with uniform distribution.
/// @details bin_index, iplus1, and data_index are inputs, mainchain_torsions vector is cleared and set by this function (i.e. it's output).
void BinTransitionCalculator::uniform_mainchain_torsions(
	//core::conformation::Residue const &res,
	core::Size const bin_index,
	bool const iplus1,
	core::Size const data_index,
	utility::vector1 <core::Real> &mainchain_torsions
) const
{
	core::Size const ndata( n_bin_transition_data() );
	if ( !bin_params_loaded() || ndata < 1 ) utility_exit_with_message( "In BinTransitionCalculator::uniform_mainchain_torsions(): no transition probability matrices have been defined!\n" );
	if ( data_index < 1 || data_index > ndata ) utility_exit_with_message( "In BinTransitionCalculator::uniform_mainchain_torsions(): the index of the BinTransitionData object is out of range.\n" );
	if ( bin_index < 1 || (!iplus1 && bin_index > bin_transition_data(data_index)->n_bins_i()) || (iplus1 && bin_index > bin_transition_data(data_index)->n_bins_iplus1()) ) {
		utility_exit_with_message( "In BinTransitionCalculator::uniform_mainchain_torsions(): the index of the bin is out of range.\n" );
	}

	BinTransitionDataCOP thisres_data( bin_transition_data_[data_index] ); //Pointer to the transition probability data for this residue.
	core::Size const ntors( (iplus1 ? thisres_data->n_mainchain_torsions_iplus1() : thisres_data->n_mainchain_torsions_i()) ); //Number of mainchain torsions.
	assert(ntors>0); //Should already be guaranteed.

	//Initialize mainchain torsions vector:
	mainchain_torsions.clear();
	mainchain_torsions.resize(ntors, 0.0);

	for ( core::Size i=1; i<=ntors; ++i ) {
		core::Real range_start = (iplus1 ? thisres_data->bin_boundaries_iplus1(bin_index,i).first : thisres_data->bin_boundaries_i(bin_index,i).first );
		core::Real range_end = (iplus1 ? thisres_data->bin_boundaries_iplus1(bin_index,i).second : thisres_data->bin_boundaries_i(bin_index,i).second );
		if ( range_start > range_end ) range_end+=360.0; //For wrap-around bins, we'll pick something in the trans-180 range, then wrap it back.
		assert( range_start < range_end ); //Should be true, now.
		mainchain_torsions[i]=numeric::random::rg().uniform()*( range_end-range_start )+range_start; //Pick a random number in the range
		if ( mainchain_torsions[i] > 180.0 ) mainchain_torsions[i]-=360.0; //Bring it back into the -180, 180 range.
		assert(mainchain_torsions[i] <= 180.0 && mainchain_torsions[i] >= -180.0); //Should be true.
	}

	return;
}

/// @brief Get the first and last lines of the next set of lines flanked by "BEGIN" and "END" lines.
/// @details This function returns "true" if another BEGIN/END pair is found, false otherwise.  It starts searching at
/// lastline+1, and continues to the end of lines.  If it finds a BEGIN/END pair, it sets firstline to the position of
/// the BEGIN line and lastline to the position of the END line.  It throws an error if a BEGIN line does not have a
/// corresponding END line.
bool BinTransitionCalculator::has_another_matrix( utility::vector1 <std::string> const &lines, core::Size &firstline, core::Size &lastline ) const {
	bool has_another(false);

	core::Size const n_lines( lines.size() ); //Number of lines

	runtime_assert_string_msg(lastline <= n_lines,
		"In core::scoring::bin_transitions::BinTransitionCalculator::has_another_matrix(): Somehow, the starting line is past the last line.  This should be impossible.  Consult a mortitian or a developer." );
	if ( lastline == n_lines ) return false; //No additional matrices if we've already reached the end of the lines.


	bool begin_found(false), end_found(false);
	for ( core::Size i=lastline+1; i<=n_lines; ++i ) {
		if ( !begin_found ) {
			if ( lines[i].substr(0,5)=="BEGIN" ) {
				begin_found=true;
				firstline=i;
				if ( TR.Debug.visible() ) TR.Debug << "BEGIN found at line " << firstline << "." << std::endl; //DELETE ME
				continue;
			}
		} else {
			if ( lines[i].substr(0,3)=="END" ) {
				end_found=true;
				has_another=true;
				lastline=i;
				if ( TR.Debug.visible() ) TR.Debug << "END found at line " << lastline << "." << std::endl; //DELETE ME
				break;
			}
		}
	}

	if ( begin_found ) {
		runtime_assert_string_msg(end_found,
			"In core::scoring::bin_transitions::BinTransitionCalculator::has_another_matrix():  A \"BEGIN\" statement was found without a corresponding \"END\" statement in the bin_params file.");
	}

	return has_another;
} //has_another_matrix()

/// @brief Parse the number of mainchain torsions for the ith and i+1st residues defined in the current block of lines
/// from the bin_params file.
void BinTransitionCalculator::parse_mainchain_torsions(
	utility::vector1 < std::string > const &lines,
	core::Size const firstline,
	core::Size const lastline,
	BinTransitionDataOP curdata
) {

	bool mainchain_torsions_i_found(false);
	bool mainchain_torsions_iplus1_found(false);

	for ( core::Size i=firstline+1, imax=lastline-1; i<=imax; ++i ) { //Loop through the lines between BEGIN and END lines.
		std::istringstream l(lines[i]);
		std::string tag;

		l >> tag;
		if ( tag == "MAINCHAIN_TORSIONS_I" ) {
			runtime_assert_string_msg( !mainchain_torsions_i_found,
				"In core::scoring::bin_transitions::BinTransitionCalculator::parse_mainchain_torsions(): Duplicate \"MAINCHAIN_TORSIONS_I\" entries found in a BEGIN/END block in the bin_params file." );
			mainchain_torsions_i_found=true;
			core::Size n_mainchain_torsions_i(0);
			l >> n_mainchain_torsions_i;
			curdata->set_n_mainchain_torsions_i( n_mainchain_torsions_i );
			continue;
		} else if ( tag == "MAINCHAIN_TORSIONS_IPLUS1" ) {
			runtime_assert_string_msg( !mainchain_torsions_iplus1_found,
				"In core::scoring::bin_transitions::BinTransitionCalculator::parse_mainchain_torsions(): Duplicate \"MAINCHAIN_TORSIONS_IPLUS1\" entries found in a BEGIN/END block in the bin_params file." );
			mainchain_torsions_iplus1_found=true;
			core::Size n_mainchain_torsions_iplus1(0);
			l >> n_mainchain_torsions_iplus1;
			curdata->set_n_mainchain_torsions_iplus1( n_mainchain_torsions_iplus1 );
			continue;
		} else continue; //Go on to the next line if this isn't a MAINCHAIN_TORSIONS line.
	}

	runtime_assert_string_msg(mainchain_torsions_i_found,
		"In core::scoring::bin_transitions::BinTransitionCalculator::parse_mainchain_torsions(): No \"MAINCHAIN_TORSIONS_I\" entry found in a BEGIN/END block in the bin_params file.");
	runtime_assert_string_msg(mainchain_torsions_iplus1_found,
		"In core::scoring::bin_transitions::BinTransitionCalculator::parse_mainchain_torsions(): No \"MAINCHAIN_TORSIONS_IPLUS1\" entry found in a BEGIN/END block in the bin_params file.");

	return;
} //parse_mainchain_torsions()

/// @brief Parse the number of torsion bins for the ith and i+1st residues as defined in the current block of lines,
/// and initialize the probability matrix to a grid of zeroes of appropiate dimensions.
void BinTransitionCalculator::initialize_matrix( utility::vector1 < std::string > const &lines, core::Size const firstline, core::Size const lastline, BinTransitionDataOP curdata)
{
	core::Size bins_i(0), bins_iplus1(0);

	for ( core::Size i=firstline+1, imax=lastline-1; i<=imax; ++i ) { //Loop through the lines between BEGIN and END lines.
		std::istringstream l(lines[i]);
		std::string tag;

		l >> tag;
		if ( tag == "BIN_COUNT_I" ) {
			runtime_assert_string_msg(bins_i==0,
				"In core::scoring::bin_transitions::BinTransitionCalculator::initialize_matrix(): Dublicate \"BIN_COUNT_I\" lines found in a BEGIN/END block in the bin_params file.");
			l >> bins_i;
			continue;
		} else if ( tag == "BIN_COUNT_IPLUS1" ) {
			runtime_assert_string_msg(bins_iplus1==0,
				"In core::scoring::bin_transitions::BinTransitionCalculator::initialize_matrix(): Dublicate \"BIN_COUNT_IPLUS1\" lines found in a BEGIN/END block in the bin_params file.");
			l >> bins_iplus1;
			continue;
		} else continue;
	}

	runtime_assert_string_msg(bins_i!=0,
		"In core::scoring::bin_transitions::BinTransitionCalculator::initialize_matrix(): No \"BIN_COUNT_I\" line found in a BEGIN/END block in the bin_params file.");
	runtime_assert_string_msg(bins_iplus1!=0,
		"In core::scoring::bin_transitions::BinTransitionCalculator::initialize_matrix(): No \"BIN_COUNT_IPLUS1\" line found in a BEGIN/END block in the bin_params file.");

	//Set the number of bins and initialize the matrix in the data object:
	curdata->set_n_bins( bins_i, bins_iplus1 );

	return;
}//initialize_matrix()

/// @brief Parse the names and angle ranges for each bin for the ith and i+1st residue.
///
void BinTransitionCalculator::set_up_bins( utility::vector1 < std::string > const &lines, core::Size const firstline, core::Size const lastline, BinTransitionDataOP curdata )
{
	core::Size bins_i_found(0), bins_iplus1_found(0); //How many bin names have we found so far?
	core::Size const n_bins_i( curdata->n_bins_i() ); core::Size const n_bins_iplus1( curdata->n_bins_iplus1() ); //How many bins should there be?
	core::Size const n_torsions_i( curdata->n_mainchain_torsions_i() ); core::Size const n_torsions_iplus1( curdata->n_mainchain_torsions_iplus1() ); //How many mainchain torsions should there be?

	for ( core::Size i=firstline+1, imax=lastline-1; i<=imax; ++i ) { //Loop through lines between BEGIN and END lines
		std::istringstream l(lines[i]);
		std::string tag;
		l >> tag;

		if ( tag=="I_BIN" ) {
			if ( TR.Debug.visible() ) TR.Debug << "Found bin: \"" << lines[i] << "\"" << std::endl; //DELETE ME

			++bins_i_found; //We've found a bin definition for residue i.
			runtime_assert_string_msg(!l.eof(), "In core::scoring::bin_transitions::BinTransitionCalculator::set_up_bins(): could not get bin name!");
			l >> tag; //Read the bin name.
			curdata->set_binname_i(bins_i_found, tag);

			std::string errmsg("In core::scoring::bin_transitions::BinTransitionCalculator::set_up_bins(): could not get all bin torsion ranges!"); //The string for error messages -- create it once.

			for ( core::Size curtors=1; curtors<=n_torsions_i; ++curtors ) {
				core::Real rangestart(0.0), rangeend(0.0);
				runtime_assert_string_msg(!l.eof(), errmsg);
				l >> rangestart;
				runtime_assert_string_msg(!l.eof(), errmsg);
				l >> rangeend;
				curdata->set_binrange_i( bins_i_found, curtors, rangestart, rangeend );
			}
		} else if ( tag=="IPLUS1_BIN" ) {
			if ( TR.Debug.visible() ) TR.Debug << "Found bin: \"" << lines[i] << "\"" << std::endl; //DELETE ME

			++bins_iplus1_found; //We've found a bin definition for residue iplus1.
			runtime_assert_string_msg(!l.eof(), "In core::scoring::bin_transitions::BinTransitionCalculator::set_up_bins(): could not get bin name!");
			l >> tag; //Read the bin name.
			curdata->set_binname_iplus1(bins_iplus1_found, tag);

			std::string errmsg("In core::scoring::bin_transitions::BinTransitionCalculator::set_up_bins(): could not get all bin torsion ranges!"); //The string for error messages -- create it once.

			for ( core::Size curtors=1; curtors<=n_torsions_iplus1; ++curtors ) {
				core::Real rangestart(0.0), rangeend(0.0);
				runtime_assert_string_msg(!l.eof(), errmsg);
				l >> rangestart;
				runtime_assert_string_msg(!l.eof(), errmsg);
				l >> rangeend;
				curdata->set_binrange_iplus1( bins_iplus1_found, curtors, rangestart, rangeend );
			}
		} else if ( tag=="SUB_BINS_I" ) {
			l >> tag; //Read the sub-bin type.
			curdata->set_subbin_type_i( tag );
		} else if ( tag=="SUB_BINS_IPLUS1" ) {
			l >> tag; //Read the sub-bin type.
			curdata->set_subbin_type_iplus1( tag );
		} else continue;
	}

	//At this point, all bins for residue i should be defined:
	if ( TR.Debug.visible() ) { //DELETE ME
		TR.Debug << "bins_i_found=" << bins_i_found << "\tn_bins_i=" << n_bins_i << std::endl; //DELETE ME
		TR.Debug << "bins_iplus1_found=" << bins_iplus1_found << "\tn_bins_iplus1=" << n_bins_iplus1 << std::endl; //DELETE ME
		TR.Debug.flush(); //DELETE ME
	}
	runtime_assert_string_msg( bins_i_found == n_bins_i,
		"In core::scoring::bin_transitions::BinTransitionCalculator::set_up_bins(): not all bins for residue i were defined in the bin_params file!" );
	//Bins for residue i+1 have either all been defined, or they're undefined and are to be copied from residue i.
	if ( bins_iplus1_found==0 ) {
		bool copy_bins(false);
		for ( core::Size i=firstline+1, imax=lastline-1; i<=imax; ++i ) {
			std::istringstream l(lines[i]);
			std::string tag;
			l >> tag;
			if ( tag=="IPLUS1_BINS_COPY_I" ) { //The bin_params file specifies that the bins for residue i+1 are identical to those for residue i.
				copy_bins=true;
				curdata->copy_i_bins_to_iplus1(); //Copy the bins for residue i to residue i+1
				break;
			}
		}
		runtime_assert_string_msg( copy_bins,
			"In core::scoring::bin_transitions::BinTransitionCalculator::set_up_bins(): no bins were defined for the i+1 residue.  Did you forget to include \"IPLUS1_BINS_COPY_I\" in the bin_params file?" );
	} else {
		runtime_assert_string_msg( bins_iplus1_found == n_bins_iplus1,
			"In core::scoring::bin_transitions::BinTransitionCalculator::set_up_bins(): not all bins for residue i+1 were defined in the bin_params file!" );
	}

	return;
} //set_up_bins()

/// @brief Parse the transition probabilities in the probability matrix, and store the values.
///
void BinTransitionCalculator::populate_matrix( utility::vector1 < std::string > const &lines, core::Size const firstline, core::Size const lastline, BinTransitionDataOP curdata )
{
	runtime_assert_string_msg( curdata->matrix_initialized(),
		"In core::scoring::bin_transitions::BinTransitionCalculator::populate_matrix(): The matrix has not yet been initialized!" );

	core::Size const n_bins_i( curdata->n_bins_i() );
	core::Size const n_bins_iplus1( curdata->n_bins_iplus1() );
	core::Size matrix_lines_found(0);

	for ( core::Size iline=firstline+1,iline_max=lastline-1; iline<=iline_max; ++iline ) { //Loop from line after BEGIN to line before END.
		std::istringstream l(lines[iline]);
		std::string tag;
		l >> tag;

		if ( tag=="MATRIX" ) {
			if ( TR.Debug.visible() ) TR.Debug << "Found a matrix line:\t\"" << lines[iline] << "\"" << std::endl; //DELETE ME
			matrix_lines_found++;
			for ( core::Size j=1; j<=n_bins_iplus1; ++j ) { //Loop through the horizontal bins corresponding to the i+1 position
				runtime_assert_string_msg(!l.eof(), "In core::scoring::bin_transitions::BinTransitionCalculator::populate_matrix(): too few entries in a \"MATRIX\" line in the bin_params file.");
				core::Real matrixval(0.0);
				l >> matrixval;
				curdata->set_matrix_entry( matrix_lines_found, j, matrixval ); //Store the value
			}
		} else continue;
	}

	runtime_assert_string_msg( matrix_lines_found==n_bins_i,
		"In core::scoring::bin_transitions::BinTransitionCalculator::populate_matrix(): Too few \"MATRIX\" lines found in the bin_params file." );

	//TEMPORARY -- write out the probability matrix:
	if ( TR.Debug.visible() ) { //DELETE ME
		TR.Debug << "Set matrix to the following:" << std::endl;
		for ( core::Size i=1; i<=n_bins_i; ++i ) {
			for ( core::Size j=1; j<=n_bins_iplus1; ++j ) {
				TR.Debug << curdata->probability_matrix(i,j) << "\t";
			}
			TR.Debug << std::endl;
		}
		TR.Debug.flush();
	}

	return;
} //populate_matrix

/// @brief Parse the required and prohibited properties for the ith and i+1st residues.
///
void BinTransitionCalculator::store_properties( utility::vector1 < std::string > const &lines, core::Size const firstline, core::Size const lastline, BinTransitionDataOP curdata ) {
	for ( core::Size i=firstline+1,imax=lastline-1; i<=imax; ++i ) { //Loop through all lines between BEGIN and END lines.
		std::istringstream l(lines[i]);
		std::string tag;
		l >> tag;
		if ( tag=="PROPERTIES_I" ) {
			while ( !l.eof() ) {
				l >> tag;
				curdata->add_property_i(tag);
			}
		} else if ( tag == "PROPERTIES_IPLUS1" ) {
			while ( !l.eof() ) {
				l >> tag;
				curdata->add_property_iplus1(tag);
			}
		} else if ( tag == "NOT_PROPERTIES_I" ) {
			while ( !l.eof() ) {
				l >> tag;
				curdata->prohibit_property_i(tag);
			}
		} else if ( tag == "NOT_PROPERTIES_IPLUS1" ) {
			while ( !l.eof() ) {
				l >> tag;
				curdata->prohibit_property_iplus1(tag);
			}
		} else continue;
	}

	curdata->check_property_overlap(); //Ensure that no properties are both required and prohibited.

	return;
} //store_properties

/// @brief Parse the required and prohibited residue types (three-letter codes) for the ith and i+1st residues.
///
void BinTransitionCalculator::store_residue_identities( utility::vector1 < std::string > const &lines, core::Size const firstline, core::Size const lastline, BinTransitionDataOP curdata )
{
	for ( core::Size i=firstline+1,imax=lastline-1; i<=imax; ++i ) { //Loop through all lines between BEGIN and END lines.
		std::istringstream l(lines[i]);
		std::string tag;

		std::string const errmsg("In core::scoring::bin_transitions::BinTransitionCalculator::store_residue_identities(): Residue identities must be provided as space-separated three-letter codes."); //Define this error message string ONCE.

		l >> tag;
		if ( tag=="RES_I" ) {
			while ( !l.eof() ) {
				l >> tag;
				runtime_assert_string_msg( tag.length()<=3, errmsg );
				curdata->add_res_identity_i(tag);
			}
		} else if ( tag == "RES_IPLUS1" ) {
			while ( !l.eof() ) {
				l >> tag;
				runtime_assert_string_msg( tag.length()<=3, errmsg );
				curdata->add_res_identity_iplus1(tag);
			}
		} else if ( tag == "NOT_RES_I" ) {
			while ( !l.eof() ) {
				l >> tag;
				runtime_assert_string_msg( tag.length()<=3, errmsg );
				curdata->prohibit_res_identity_i(tag);
			}
		} else if ( tag == "NOT_RES_IPLUS1" ) {
			while ( !l.eof() ) {
				l >> tag;
				runtime_assert_string_msg( tag.length()<=3, errmsg );
				curdata->prohibit_res_identity_iplus1(tag);
			}
		} else continue;
	}

	curdata->check_residentity_overlap(); //Ensure that no residue identities are both required and prohibited.

	return;
} //store_residue_identities

/// @brief Is thisres connected to prevres at the lower connection of thisres and the upper connection of prevres?
/// @details Returns true or false based on whether this is the case.
bool BinTransitionCalculator::are_normally_bonded(
	core::conformation::Residue const &thisres,
	core::conformation::Residue const &prevres
) const {
	if ( thisres.has_lower_connect() && prevres.has_upper_connect() ) {
		core::Size const lower_id_this( thisres.type().lower_connect_id() );
		core::Size const upper_id_that( prevres.type().upper_connect_id() );
		if ( thisres.residue_connection_partner(lower_id_this) == prevres.seqpos() && prevres.residue_connection_partner(upper_id_that) == thisres.seqpos() ) return true;
	}
	return false;
} //are_normally_bonded

/// @brief Prints a report summarizing the data stored in this object, including
/// all of the BinTransitionData objects.
/// @details If verbose is true, the full sub-bin information is printed, too.
std::string BinTransitionCalculator::summarize_stored_data( bool const verbose ) const
{
	std::ostringstream outstream;

	outstream << std::endl << "********************" << std::endl;
	outstream << "Summary of data stored in the BinTransitionCalculator" << std::endl << std::endl;
	outstream << "bin_params_file: " << bin_params_file_ << std::endl;

	for ( core::Size i=1, imax=bin_transition_data_.size(); i<=imax; ++i ) { //Loop through all BinTransitionData objects.
		BinTransitionDataOP curdata( bin_transition_data_[i] );
		outstream << "----------" << std::endl << "BinTransitionData " << i << std::endl;
		outstream << curdata->summarize_data(verbose) << "----------" << std::endl;
	}

	outstream << "End of report." << std::endl << "********************" << std::endl;

	return outstream.str();

} //summarize_stored_data


} //namespace bin_transitions
}//namespace scoring
}//namespace core
