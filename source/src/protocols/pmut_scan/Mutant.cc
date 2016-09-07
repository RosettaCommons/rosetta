// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pmut_scan/Mutant.cc
/// @brief A helper class for the point mut scan protocol which holds data about mutants.
/// @author Ron Jacak


// Project Headers
#include <protocols/pmut_scan/Mutant.hh>

#include <core/types.hh>


// Utility headers

// ObjexxFCL header

// C++
#include <iostream>

#include <utility/vector1.hh>
#include <sstream>


namespace protocols {
namespace pmut_scan {

///
/// @brief
/// Constructor for MutationData objects. A mutation holds the wt and mutant amino acid type, pdb resnum, pose resnum, chain and icode.
///
MutationData::MutationData( char wt_residue, char mut_residue, core::Size pose_resnum, core::Size pdb_resnum, char icode, char chain ) :
	wt_residue_( wt_residue ),
	mut_residue_( mut_residue ),
	pose_resnum_( pose_resnum ),
	pdb_resnum_( pdb_resnum ),
	icode_( icode ),
	chain_( chain )
{}

///
/// @brief
/// Destructor for MutationData objects. No dynamically allocated memory held in MutationData objects so nothing to do here.
///
MutationData::~MutationData() = default;

///
/// @brief
/// Returns a string representation of this mutation.
///
std::string MutationData::mutation_string() const {
	std::stringstream out;
	out << chain_ << "-" << wt_residue_ << pose_resnum_ << mut_residue_;
	return out.str();
}

///
/// @brief
/// Returns a string representation of this mutation using PDB not pose numbering.
///
std::string MutationData::mutation_string_PDB_numbering() const {
	std::stringstream out;
	out << chain_ << "-" << wt_residue_ << pdb_resnum_;
	if ( icode_ != ' ' ) {
		out << icode_;
	}
	out << mut_residue_;
	return out.str();
}

///
/// @brief
/// Accessor for the mut_residue member variable.  Needed by the function make_mutant_structure() in the PointMutScanDriver class.
///
char MutationData::mut_residue() const {
	return mut_residue_;
}

///
/// @brief
/// Accessor for the pose_resnum member variable.  Needed by the function make_mutant_structure() in the PointMutScanDriver class.
///
core::Size MutationData::pose_resnum() const {
	return pose_resnum_;
}

///
/// @brief
/// print function for MutationData class; only used by unit tests
///
void MutationData::print_mutation_data( MutationData & md ) {
	std::cout << md.wt_residue_ << md.pose_resnum_ << md.mut_residue_ << " (pdb chain/res: " << md.chain_ << "/" << md.pdb_resnum_ << ", icode: '" << md.icode_ << "')";
}

///
/// @brief
/// function which tests two MutationData objects for equality; only used by unit tests
///
bool MutationData::operator==( const MutationData & md_other ) const {
	if ( wt_residue_ == md_other.wt_residue_ && mut_residue_ == md_other.mut_residue_ &&
			pose_resnum_ == md_other.pose_resnum_ && pdb_resnum_ == md_other.pdb_resnum_ &&
			icode_ == md_other.icode_ && chain_ == md_other.chain_ ) {
		return true;
	}
	return false;
}


/// @brief
/// Mutant class constructor
Mutant::Mutant() {}

/// @brief
/// Mutant class destructor
Mutant::~Mutant() = default;

///
/// @brief
/// Returns the number of mutations in this mutant.
///
core::Size Mutant::n_mutations() const {
	return mutations_.size();
}

///
/// @brief
/// Adds the passed in mutation to the class member list.
///
void Mutant::add_mutation( MutationData md ) {
	mutations_.push_back( md );
}

///
/// @brief
/// Returns a const iterator to beginning of the mutations vector
///
utility::vector1< MutationData >::const_iterator Mutant::mutations_begin() const {
	return mutations_.begin();
}

///
/// @brief
/// Returns a const iterator to end of the mutations vector
///
utility::vector1< MutationData >::const_iterator Mutant::mutations_end() const {
	return mutations_.end();
}

///
/// @brief
/// Function which tests two Mutant objects for equality; only used by unit tests
///
bool Mutant::operator==( const Mutant & m_other ) const {

	if ( this->n_mutations() != m_other.n_mutations() ) { return false; }
	utility::vector1< MutationData >::const_iterator this_iter, m_other_iter;

	for ( this_iter = mutations_begin(), m_other_iter = m_other.mutations_begin(); ((this_iter != mutations_end()) && (m_other_iter != m_other.mutations_end())); ++this_iter, ++m_other_iter ) {

		if ( !( (*this_iter).operator==(*m_other_iter) ) ) {
			return false;
		}
	}

	return true;
}

///
/// @brief
/// Sets the passed in reference to the first element of the mutations_ vector, and removes that element from the vector.
///
MutationData Mutant::pop_mutation() {

	MutationData md = mutations_.front(); // uses MutationData assignment operator!  does this need to be defined?
	mutations_.erase( mutations_.begin() );

	return md;
}

std::ostream & operator<< ( std::ostream & os, const MutationData & md ) {
	os << md.chain_ << "/" << md.wt_residue_ << md.pose_resnum_ << md.mut_residue_ << " (pdb chain/res/icode: " << md.chain_ << "/" << md.pdb_resnum_ << md.icode_ << ")";
	return os;
}


std::ostream & operator<< ( std::ostream & os, const Mutant & m ) {

	utility::vector1< MutationData >::const_iterator m_iter;
	for ( m_iter = m.mutations_begin(); m_iter != m.mutations_end(); ++m_iter ) {
		// don't print " + " token on first time through
		if ( m_iter != m.mutations_begin() ) {
			os << " + ";
		}
		os << *(m_iter);
	}
	return os;
}


} // namespace pmut_scan
} // namespace protocols

