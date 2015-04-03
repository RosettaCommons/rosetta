// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/pmut_scan/Mutant.hh
/// @brief A helper class for the point mut scan protocol.
/// @author Ron Jacak

#ifndef INCLUDED_protocols_pmut_scan_Mutant_HH
#define INCLUDED_protocols_pmut_scan_Mutant_HH

// Project Headers
#include <protocols/pmut_scan/Mutant.fwd.hh>

#include <core/types.hh>


// Utility headers

// ObjexxFCL header

// C++

#include <utility/vector1.hh>


namespace protocols {
namespace pmut_scan {

class MutationData {

public:

	MutationData( char wt_residue, char mut_residue, core::Size pose_resnum, core::Size pdb_resnum, char icode, char chain );
	~MutationData();

	std::string mutation_string() const;
	std::string mutation_string_PDB_numbering() const;

	char mut_residue() const;
	core::Size pose_resnum() const;

	char pdb_chain() const { return chain_; }

	// print function for mutation data class; only used by unit tests
	void print_mutation_data( MutationData & md );

	// function which tests two mutation_data objects for equality; only used by unit tests
	bool operator==( const MutationData & md_other ) const;

	friend std::ostream & operator<< ( std::ostream & os, const MutationData & md );

private:

	char wt_residue_;
	char mut_residue_;
	core::Size pose_resnum_;
	core::Size pdb_resnum_; // position in PDB numbering
	char icode_; // insertion code
	char chain_; // insertion code


friend class Mutant;
friend class PointMutScanDriver;

}; // class MutationData


class Mutant {

public:

	Mutant();
	~Mutant();

	core::Size n_mutations() const;
	void add_mutation( MutationData md );

	utility::vector1< MutationData >::const_iterator mutations_begin() const;
	utility::vector1< MutationData >::const_iterator mutations_end() const;

	// function which tests two mutation_data objects for equality; only used by unit tests
	bool operator==( const Mutant & m_other ) const;

	MutationData pop_mutation();

	friend std::ostream & operator<< ( std::ostream & os, const Mutant & m );

private:

	utility::vector1< MutationData > mutations_;

friend class PointMutScanDriver;

}; // class Mutant


} // namespace pmut_scan
} // namespace protocols

#endif //INCLUDED_protocols_pmut_scan_Mutant_HH
