// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/bin_transitions/BinTransitionCalculator.hh
/// @brief  Headers for BinTransitionCalculator class.
/// @details This class loads data associated with transitions from one mainchain torsion bin to another (e.g. ABEGO bins, OO-ABBA bins, etc.)
/// and provides methods for generating bin sequences randomly, perturbing sequences, etc.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_scoring_bin_transitions_BinTransitionCalculator_hh
#define INCLUDED_core_scoring_bin_transitions_BinTransitionCalculator_hh

//BinTransitionCalculator owning pointers header:
#include <core/scoring/bin_transitions/BinTransitionCalculator.fwd.hh>
#include <core/scoring/bin_transitions/BinTransitionData.fwd.hh>
#include <core/scoring/bin_transitions/BinTransitionData.hh>

//Other headers:
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

///////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace bin_transitions {

class BinTransitionCalculator : public utility::pointer::ReferenceCount {
public: //Constructors and destructors:
	/// @brief Default constructor for BinTransitionCalculator
	///
	BinTransitionCalculator();

	/// @brief Copy constructor for BinTransitionCalculator
	///
	BinTransitionCalculator( BinTransitionCalculator const &src );

	/// @brief Default destructor for BinTransitionCalculator
	///
	~BinTransitionCalculator();

	/// @brief Clone operation for BinTransitionCalculator.
	/// @details Returns an owning pointer to a copy of this object.
	BinTransitionCalculatorOP clone() const;

public: //Public functions -- information:

	/// @brief Prints a report summarizing the data stored in this object, including
	/// all of the BinTransitionData objects.
	/// @details If verbose is true, the full sub-bin information is printed, too.
	std::string summarize_stored_data( bool const verbose ) const;

public: //Public functions -- calculators and initializers:

	/// @brief Set the bin params file, and load all the bin params and transition probabilities.
	///
	void load_bin_params( std::string const &filename );

	/// @brief Given a bin name and a residue, find a BinTransitionsData object describing that residue,
	/// and the index of the bin within that object.
	/// @details data_index and bin_index are outputs.  Both are set to 0 if the search fails.  Everything
	/// else is const input.  If bin_name is set to "", then the bin index that is returned is the index
	/// corresponding to the residue's mainchain torsion vector.
	void find_data_and_bin(
		std::string const &bin_name,
		core::conformation::Residue const &res,
		core::Size &data_index,
		core::Size &bin_index,
		bool const use_iplus1
	) const;

	/// @brief Given residues at positions i and iplus1, find a bin transition probability data object that
	/// describes the pair.
	/// @details Inputs are res_i and res_iplus1.  Outputs are data_i and data_iplus1 (BinTransitionData
	/// indices).  Function returns true if data are found successfully, false if no BinTransitionData object
	/// could be found describing the residues in question.
	bool find_data(
		core::conformation::Residue const &res_i,
		core::conformation::Residue const &res_iplus1,
		core::Size &data_index
	) const;

	/// @brief Given two residues (rsd_i and rsd_iplus1 at positions i and i+1, respectively), give me the
	/// probability of seeing rsd_iplus1 in its bin given that rsd_i is in its bin.
	/// @details The probability value is set to a number from 0 to 1.  Inputs are rsd_i and rsd_iplus1.
	/// Function returns true if data are found successfully, false if no BinTransitionData object
	/// could be found describing the residues in question.
	bool p_iplus1_given_i(
		core::conformation::Residue const &rsd_i,
		core::conformation::Residue const &rsd_iplus1,
		core::Real &probability
	) const;

	/// @brief Given two residues (rsd_i and rsd_iplus1 at positions i and i+1, respectively), give me the
	/// probability of seeing rsd_i in its bin given that rsd_iplus1 is in its bin.
	/// @details The probability value is set to a number from 0 to 1.  Inputs are rsd_i and rsd_iplus1.
	/// Function returns true if data are found successfully, false if no BinTransitionData object
	/// could be found describing the residues in question.
	bool p_i_given_iplus1(
		core::conformation::Residue const &rsd_i,
		core::conformation::Residue const &rsd_iplus1,
		core::Real &probability
	) const;

	/// @brief Is the given residue in the given bin?
	/// @details For the bin definitions, this uses the first BinTransitionsData object that it finds where residue i matches the properties of the given
	/// residue.  Checks the i+1 definitions if nothing is found for the i transitions.  Fails if no bin definitions for the given residue type are found.
	bool is_in_bin (
		core::conformation::Residue const &res,
		std::string const &bin_name
	) const;

	/// @brief Is the given residue in the bin given by a data object index and a bin index?
	///
	bool is_in_bin (
		core::conformation::Residue const &res,
		core::Size const data_index,
		core::Size const bin_index,
		bool const use_iplus1
	) const;

	/// @brief Initialize a string of residues to a bunch of random bins, based on bin transition probabilities; then draw random mainchain torsion angles from those bins.
	/// @details Takes a const conformation and a const list of residue indices as input; the conformation is just for checking residues types, numbers of mainchain torsions, etc.
	/// The residue indices must be in order, defining a contiguous chain (running backwards or forwards).  Output is the mainchain_torsions vector of vectors (reset and
	/// replaced by this operation).  The distribution WITHIN the bin depends on the BinTransitionData object and what was specified in the bin_params file.  Default is uniform
	/// within each bin, though Ramachandran-biased distributions are also permitted for alpha-amino acids.
	void random_mainchain_torsions_from_bins(
		core::conformation::Conformation const &conformation,
		utility::vector1 <core::Size> const &res_indices,
		utility::vector1 < utility::vector1 < core::Real > > &mainchain_torsions
	) const;

	/// @brief Draw random mainchain torsion values for a set of residues, given a bin from which the values should be drawn.
	/// @details Takes a bin name, a const conformation, and a const list of residue indices as input; the conformation is just for checking residues types, numbers of mainchain
	/// torsions, etc.  Output is the mainchain_torsions vector of vectors (reset and replaced by this operation).  The distribution WITHIN the bin depends on the BinTransitionData
	/// object and what was specified in the bin_params file.  Default is uniform within each bin, though Ramachandran-biased distributions are also permitted for alpha-amino acids.
	/// Note that this function uses bins for residue i, and only checks i+1 if no suitable data are found for i.
	void random_mainchain_torsions_from_bin(
		std::string const &bin_name,
		core::conformation::Conformation const &conformation,
		utility::vector1 <core::Size> const &res_indices,
		utility::vector1 < utility::vector1 < core::Real > > &mainchain_torsions
	) const;

	/// @brief Randomly pick mainchain torsions for a residue based on the torsion bins of its i+1 and i-1 neighbours.
	/// @details Takes a const conformatoin, a const residue index, and a boolean valueas input.  The conformation is for checking residue
	/// types, numbers of mainchain torsions, etc.  The boolean determines whether this residue should be allowed to stay in its own bin
	/// or be required to switch to another bin.  Output is the mainchain_torsions vector of Reals, with one entry for each mainchain torsion
	/// of the residue.
	void random_mainchain_torsions_using_adjacent_bins(
		core::conformation::Conformation const &conformation,
		core::Size const res_index,
		bool const must_switch_bins,
		utility::vector1 < core::Real > &mainchain_torsions
	) const;

	/// @brief Randomly pick a bin, given an input residue.
	/// @details The residue's properties and name are used to find the first BinTransitionData object matching the properties and name.
	/// A random bin is then picked based on the frequency with which that bin is observed for that class of residues.  If the bin ends
	/// up being from the i+1st residue, iplus1 is set to "true"; otherwise, it remains "false".  The data_index variable is set by this
	/// function to the index of the BinTransitionData object that provides the transition probability data.
	core::Size random_bin(
		core::conformation::Residue const &res,
		bool &iplus1,
		core::Size &data_index
	) const;

	/// @brief Randomly pick mainchain torsion values from within a bin
	/// @details The residue's properties and name are used to find the first BinTransitionData object matching the properties and name.
	/// Mainchain torsion values are then picked randomly from within the bin.  The distribution WITHIN the bin depends on the BinTransitionData
	/// object and what was specified in the bin_params file.  Default is uniform within each bin, though Ramachandran-biased distributions are
	/// also permitted for alpha-amino acids.  If iplus1 is true, we draw from the bins for the i+1st residue; otherwise, we draw from the bins
	/// for the ith residue.  The data_index value tells this function which BinTransitionData object to use.
	void random_mainchain_torsions_from_bin(
		core::Size const bin_index,
		bool const iplus1,
		core::Size const data_index,
		utility::vector1 < core::Real > &mainchain_torsions //output
	) const;

	/// @brief Given the current residue and the previous residue, as well as the mainchain torsion values for the previous residue,
	/// pick a bin for the current residue based on bin transition probabilities.
	/// @details This function takes thisres (the current residue) and prevres (the previous residue) as inputs, using them for their
	/// properties and identities to find a BinTransitionData object that matches the criteria at the i and i+1st residues.  It then
	/// sets data_index to the index of that BinTransitionData object before returning the index of a bin chosen for the current (i+1st)
	/// residue at random, based on the relative probabilities of bins given the previous (ith) residue.
	core::Size random_bin_based_on_previous (
		core::conformation::Residue const &thisres,
		core::conformation::Residue const &prevres,
		utility::vector1 <core::Real> const &prev_torsions,
		core::Size &data_index
	) const;

	/// @brief Given the current residue and its conformation, pick a bin for the current residue based on the bins of the previous and
	/// next residues, and the bin transition probabilities.
	/// @details Bin transition probabilities are assumed to be independent.  That is, P(ABC) = P(B | A & C) = P (B | A ) * P(B | C).
	/// This function also finds and returns the indices of the BinTransitionData objects that describe the AB transition (data_index_iplus1)
	/// and the BC transition (data_index_i).  So conformation and res_index are const inputs, and data_index_i and data_index_iplus1 are outputs.
	/// If must_switch_bins is true, then the current bin for the current residue is never chosen; otherwise, it is in the pool to be picked.
	core::Size random_bin_based_on_prev_and_next(
		core::conformation::Conformation const &conformation,
		core::Size const res_index,
		bool const must_switch_bins,
		core::Size &data_index_i,
		core::Size &data_index_iplus1
	) const;

public: //Public functions -- getters:

	/// @brief Return whether the bin params file has been loaded.
	///
	bool bin_params_loaded() const { return bin_params_loaded_; }

	/// @brief Number of BinTransitionData objects stored (i.e. number of transition probability matrices stored).
	///
	core::Size n_bin_transition_data() const { return bin_transition_data_.size(); }

	/// @brief Is a particular bin defined for at least one residue type?
	///
	bool bin_definition_exists( std::string const &name ) const;

public: //Public functions -- setters:

private: //Private functions:

	/// @brief Given a bin, generate phi and psi values from within the bin, biased by the sub-bin cumulative probability distribution function.
	/// @details bin_index, iplus1, and data_index are inputs, mainchain_torsions vector is cleared and set by this function (i.e. it's output).
	void biased_mainchain_torsions(
		//core::conformation::Residue const &res,
		core::Size const bin_index,
		bool const iplus1,
		core::Size const data_index,
		utility::vector1 <core::Real> &mainchain_torsions
	) const;


	/// @brief Given a bin, generate phi and psi values randomly from within the bin, with uniform distribution.
	/// @details bin_index, iplus1, and data_index are inputs, mainchain_torsions vector is cleared and set by this function (i.e. it's output).
	void uniform_mainchain_torsions(
		//core::conformation::Residue const &res,
		core::Size const bin_index,
		bool const iplus1,
		core::Size const data_index,
		utility::vector1 <core::Real> &mainchain_torsions
	) const;

	/// @brief Get the first and last lines of the next set of lines flanked by "BEGIN" and "END" lines.
	/// @details This function returns "true" if another BEGIN/END pair is found, false otherwise.  It starts searching at
	/// lastline+1, and continues to the end of lines.  If it finds a BEGIN/END pair, it sets firstline to the position of
	/// the BEGIN line and lastline to the position of the END line.  It throws an error if a BEGIN line does not have a
	/// corresponding END line.
	bool has_another_matrix( utility::vector1 <std::string> const &lines, core::Size &firstline, core::Size &lastline ) const;

	/// @brief Add a BinTransitionData object to the list of BinTransitionData objects, and return an owning pointer to the
	/// newly-added object.
	/// @details Each BinTransitionData object stores the BinTransition matrix for a particular i / i+1 pair.  For example,
	/// there might be one matrix for i=L-amino acid, i+1=L-amino acid, and another for i=D-amino acid, i+1=D-amino acid.
	BinTransitionDataOP add_bin_transition_data() {
		bin_transition_data_.push_back( BinTransitionDataOP( new BinTransitionData ) );
		return bin_transition_data_[bin_transition_data_.size()];
	}

	/// @brief Parse the number of mainchain torsions for the ith and i+1st residues defined in the current block of lines
	/// from the bin_params file.
	void parse_mainchain_torsions( utility::vector1 < std::string > const &lines, core::Size const firstline, core::Size const lastline, BinTransitionDataOP curdata );

	/// @brief Parse the number of torsion bins for the ith and i+1st residues as defined in the current block of lines,
	/// and initialize the probability matrix to a grid of zeroes of appropiate dimensions.
	void initialize_matrix( utility::vector1 < std::string > const &lines, core::Size const firstline, core::Size const lastline, BinTransitionDataOP curdata );

	/// @brief Parse the names and angle ranges for each bin for the ith and i+1st residue.
	///
	void set_up_bins( utility::vector1 < std::string > const &lines, core::Size const firstline, core::Size const lastline, BinTransitionDataOP curdata );

	/// @brief Parse the transition probabilities in the probability matrix, and store the values.
	///
	void populate_matrix( utility::vector1 < std::string > const &lines, core::Size const firstline, core::Size const lastline, BinTransitionDataOP curdata );

	/// @brief Parse the required and prohibited properties for the ith and i+1st residues.
	///
	void store_properties( utility::vector1 < std::string > const &lines, core::Size const firstline, core::Size const lastline, BinTransitionDataOP curdata );

	/// @brief Parse the required and prohibited residue types (three-letter codes) for the ith and i+1st residues.
	///
	void store_residue_identities( utility::vector1 < std::string > const &lines, core::Size const firstline, core::Size const lastline, BinTransitionDataOP curdata );

	/// @brief Access one of the bin_transition_data objects (const-access).
	///
	BinTransitionDataCOP bin_transition_data( core::Size const index ) const {
		if ( index < 1 || index > bin_transition_data_.size() ) utility_exit_with_message( "In BinTransitionCalculator::bin_transition_data(): index out of range!" );
		return( bin_transition_data_[index] );
	}

	/// @brief Is thisres connected to prevres at the lower connection of thisres and the upper connection of prevres?
	/// @details Returns true or false based on whether this is the case.
	bool are_normally_bonded(
		core::conformation::Residue const &thisres,
		core::conformation::Residue const &prevres
	) const;

private: //Private variables:

	/// @brief Has a bin params file been loaded?
	/// @details Default false; set to true by the load_bin_params() function.
	bool bin_params_loaded_;

	/// @brief The filename for the bin params file that defines the bins and their transition probabilities.
	///
	std::string bin_params_file_;

	/// @brief Vector1 of BinTransitionData objects.
	/// @details Each BinTransitionData object stores the BinTransition matrix for a particular i / i+1 pair.  For example,
	/// there might be one matrix for i=L-amino acid, i+1=L-amino acid, and another for i=D-amino acid, i+1=D-amino acid.
	BinTransitionDataOPs bin_transition_data_;

}; //BinTransitionCalculator class

} //namespace bin_transitions
}//namespace scoring
}//namespace core

#endif //INCLUDED_core_scoring_bin_transitions_BinTransitionCalculator_hh
