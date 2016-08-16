// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/antibody/grafting/scs_blast.hh
/// @brief Structural Component Selector (SCS) implementation with NCBI-BLAST+
/// @author Sergey Lyskov


/// @details (Brian/Sergey correct me if I'm wrong) - JRJ
/// SCS_Result is a general struct containing a vector of PDBs. These PDBs could come from anywhere.
///
/// SCS_BlastResult is a special case where these PDBs come from a BLAST+ alignment,
/// The SCS_BlastResult struct inherits from SCS_Result, so it also contains further information
/// about the aligned template, not the query -- i.e. the sequence similar to the one you are looking for.
///
/// SCS_ResultVector contains a vector of SCS_Result structs such as SCS_BlastResult.
/// It does not contain a list of SCS_Results!
///
/// SCS_Results contains alignment results for each antibody region.
/// Each alignment result is itself an SCS_ResultVector.
///
/// SCS_ResultSet is a slice of the n-th alignment of each region.
///
/// SCS_Base is base class containing a vector of filters and a sorter.
/// It has adders/setters for the filters/sorter.
/// More importantly, it has two functions: select / raw_select which take a query sequence.
/// The raw_select function should be overwritten for the specific alignment method. (see SCS_BlastPlus)
/// The select function calls raw_select then filters and sorts (if these have been applied).
///
/// SCS_BlastPlus is where the actual BLASTing occurs!
/// It overwrites raw_select to use BLAST+.

#ifndef INCLUDED_protocols_antibody_grafting_scs_blast_hh
#define INCLUDED_protocols_antibody_grafting_scs_blast_hh

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__


#include <protocols/antibody/grafting/scs_blast.fwd.hh>

#include <protocols/antibody/grafting/antibody_sequence.hh>
#include <protocols/antibody/grafting/chothia_numberer.hh>
#include <protocols/antibody/grafting/scs_functor.fwd.hh>
#include <protocols/antibody/grafting/exception.hh>

#include <core/types.hh>

#include <basic/report.hh>

#include <map>
#include <string>

namespace protocols {
namespace antibody {
namespace grafting {

/// @brief Base struct for any SCS single result. Hold information of selected template
struct SCS_Result
{
	virtual ~SCS_Result() {}  // This datatype is pure struct but we need to add at least one virtual method so we can use dynamic_cast later

	std::string pdb;

	// true if this result originated from results-padding instead of 'normal' tempalte selection
	bool padded = false;
};

struct SCS_Antibody_Database_Result : public SCS_Result
{
	core::Real resolution;
	std::string bio_type, light_type, struct_source;

	/// sequences of selected template (do not confuse with query sequences! 'sequence' will hold value corresponding to region results column)
	std::string /*sequence, */h1, h2, h3, frh, l1, l2, l3, frl;
};

struct SCS_BlastResult : public SCS_Antibody_Database_Result
{
	int alignment_length;
	core::Real bit_score, identity;
};



struct SCS_ResultVector : public utility::vector0< SCS_ResultOP >
{
	//std::string name, sequence;  // query info: name of the region and it sequence
};


struct SCS_Results
{
	//AntibodySequence const antibody_sequence; // our initial query for reference

	SCS_ResultVector h1, h2, h3, l1, l2, l3, frh, frl, orientation;

	SCS_Results& operator=(SCS_Results&&) = default;


	/// @brief Create result set from 'row' element of each SCS_ResultVector vector
	///        if strict is true the throw if no resuls found otherwise use empty OP
	///
	/// @throw std::out_of_range if for some of the row is not present and strict==true
	SCS_ResultSet get_result_set(uint row, bool strict=true);
};


// 'Horisontal' cut though SCS results. We storing OP here to allow fields to be optional and to accomodate for various results type (and possible mix of them)
struct SCS_ResultSet
{
	SCS_ResultOP h1, h2, h3, l1, l2, l3, frh, frl, orientation;
};


class SCS_Base : public basic::Reporter {
	std::vector<SCS_FunctorCOP> filters_;
	SCS_FunctorCOP sorter_;

public:
	typedef std::string string;

	using Reporter::Reporter;

	virtual ~SCS_Base() {}

	/// Results post-processing options: use filter to filter out resutls and sorting for chnaging results priority
	//void clear_filters();
	void add_filter(SCS_FunctorCOP filter);
	void set_sorter(SCS_FunctorCOP sorter);

	/// @brief Select CDR's template without filtering or sorting. In general you probably need to call select(...) instead
	/// @throw _AE_scs_failed_ on failure
	virtual SCS_ResultsOP raw_select(AntibodySequence const &) = 0;

	/// @brief Select CDR's template, filter it and sort. Try to provide at least 'n' templates if possible
	/// @throw _AE_scs_failed_ on failure
	virtual SCS_ResultsOP select(uint n, AntibodySequence const &);

	/// @brief Pad results vectors for each region (if possible) by adding arbitraty but compatible templates so at least n templates for each region is avalible
	virtual void pad_results(uint n, AntibodySequence const &, SCS_Results &) = 0;

protected:
	/// @brief output summary to 'Report'
	void report(SCS_ResultsOP r, uint n);
};


class SCS_LoopOverSCs : public SCS_Base
{
protected:
	typedef utility::vector0< std::map<string, string> > Antibody_SCS_Database;

private:
	std::string database_path_; // path to antibody grafting database root, should point to tools/antibody

	Antibody_SCS_Database antibody_scs_database_;

public:
	using SCS_Base::SCS_Base;


	/// @brief Select CDR's template
	/// @throw _AE_scs_failed_ on failure
	SCS_ResultsOP raw_select(AntibodySequence const &) override;  //, AntibodyNumbering const &);


	/// @brief set working dir/output-prefix for intermediate files based on command-line options
	void init_from_options();

	/// @brief set path to antibody grafting database root, should point to tools/antibody
	void set_database_path(string const &database) { database_path_ = database; }

	// /// @brief return database path
	// string database_path(void) const { return database_; }

protected:
	struct Result { std::string name, sequence; SCS_ResultVector &results; };

	virtual void select_template(Result & j,
								 std::string const & db_to_query,
								 std::map< std::string, std::map< std::string, std::string> > const & ab_db ) const = 0;

	/// Return Antibody_SCS_Database object by ether reading it from database or from cache if it was read before
	Antibody_SCS_Database &antibody_scs_database();

};


class SCS_BlastPlus : public SCS_LoopOverSCs
{
	typedef std::string string;

public:
	using SCS_LoopOverSCs::SCS_LoopOverSCs;

	/// @brief Select CDR's template
	/// @throw _AE_scs_failed_ on failure
	// SCS_ResultsOP raw_select(AntibodySequence const &) override;  //, AntibodyNumbering const &);

	/// Pad results vectors for each region (if possible) by adding arbitraty but compatible templates so at least n templates for each region is avalible
	void pad_results(uint n, AntibodySequence const &, SCS_Results &) override;

	/// @brief set working dir/output-prefix for intermediate files based on command-line options
	void init_from_options();

	/// @brief set custom working dir/output-prefix for intermediate files
	void set_output_prefix(std::string const &prefix) { prefix_ = prefix; }

	/// @brief set path to NCBI-Blast+ executable
	void set_blastp_executable(std::string const &blastp) { blastp_ = blastp; }

private:
	std::string prefix_; // output prefix for intermediate files

	std::string blastp_ = "blastp"; // path to NCBI-Blast+ executable

	void select_template( Result & j,
	                      std::string const & db_to_query, std::map< std::string,
						  std::map< std::string, std::string> > const & ab_db ) const override;
};


struct FRH_FRL
{
	std::string frh1, frh2, frh3, frh4,
				frl1, frl2, frl3, frl4;

	// struct Residue {
	// 	char aa;
	// 	std::string number;
	// };
	// utility::vector0< Residue > heavy, light;
};

/// @brief Calculate antibody framework regions (frh1...frl4) and Trim it by removing less-preserved elements
void trim_framework(AntibodySequence const &A, AntibodyFramework &heavy_fr, AntibodyFramework &light_fr);


/// @brief Calculate frh+frl pair by trimming fr* and removing 'less-preserved' regions
FRH_FRL calculate_frh_frl(AntibodySequence const &);

void populate_results_from_db( SCS_Antibody_Database_ResultOP const & result,
                               std::map< std::string, std::map< std::string, std::string > > const & db );


} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__


#endif // INCLUDED_protocols_antibody_grafting_scs_blast_hh
