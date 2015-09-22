// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/DFIRE_Potential.cc
/// @author James Thompson

#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/scoring/methods/dfire/DFIRE_Potential.hh>

#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <boost/algorithm/string/trim.hpp>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {
namespace dfire {

// @brief Auto-generated virtual destructor
DFIRE_Potential::~DFIRE_Potential() {}

static THREAD_LOCAL basic::Tracer TR( "core.scoring.methods.DFIRE_Potential" );

static Size const INVALID(999);

DFIRE_Potential::DFIRE_Potential()
: potential_is_loaded_(false)
{}

std::string const get_joint_id(
	std::string const & res1,
	std::string const & atom1,
	std::string const & res2,
	std::string const & atom2
) {
	std::string atom1_trimmed(atom1);
	std::string atom2_trimmed(atom2);
	boost::algorithm::trim(atom1_trimmed);
	boost::algorithm::trim(atom2_trimmed);
	using std::string;
	if ( res1 < res2 ) {
		string const joint_id1( res1 + " " + atom1_trimmed + " " + res2 + " " + atom2_trimmed);
		return joint_id1;
	}
	string const joint_id2( res2 + " " + atom2_trimmed + " " + res1 + " " + atom1_trimmed );
	return joint_id2;
}

void
DFIRE_Potential::read_potential(std::string const & fn) {
	using core::Size;
	using core::Real;
	using std::string;
	using utility::vector1;

	//std::cout << "reading potential!" << std::endl;

	utility::io::izstream input;
	basic::database::open( input, fn );
	if ( !input.good() ) {
		std::string const msg( "Error opening file: " + fn );
		utility_exit_with_message( msg );
	}
	string line;
	while ( getline(input,line) ) {
		if ( line.substr(0,1) != "#" ) {
			std::istringstream ss(line);
			string res1(""), atom1(""), res2(""), atom2("");
			ss >> res1 >> atom1 >> res2 >> atom2;

			//Size const res_idx1(res_index(res1)), res_idx2(res_index(res2));
			//Size const atom_idx1(atom_index(atom1)), atom_idx2(atom_index(atom2));

			//std::cout << "line = " << line << std::endl;
			//if ( res_idx1 != INVALID && res_idx2 != INVALID && atom_idx1 != INVALID && atom_idx2 != INVALID ) {
			vector1< Real > pair_potential; // for (res1,atom1,res2,atom2)
			Real value(-999);
			ss >> value;
			while ( !ss.fail() ) {
				pair_potential.push_back(value);
				ss >> value;
				//std::cout << "value = " << value << std::endl;
			}

			potential_.push_back( pair_potential );
			string const joint_id(get_joint_id(res1,atom1,res2,atom2));
			atom_res_idx_[joint_id] = potential_.size();
			//std::cout << "joint_id = (" << joint_id << ")" << std::endl;
			//std::cout << "potential size = " << potential_.size() << std::endl;
			//std::cout << "pair_potential size = " << pair_potential.size() << std::endl;
			//std::cout << joint_id;
			//for ( Size ii = 1; ii <= pair_potential.size(); ++ii ) {
			// std::cout << " " << pair_potential[ii];
			//}
			//std::cout << std::endl;
			//}
		}
	} // getline

	potential_is_loaded_ = true;
} // read_potential


core::Size DFIRE_Potential::res_index(
	std::string const & input
) const {
	std::string res_name(input);
	boost::algorithm::trim(res_name);
	if      ( res_name == "ASP" ) return  1;
	else if ( res_name == "PRO" ) return  2;
	else if ( res_name == "LYS" ) return  3;
	else if ( res_name == "ILE" ) return  4;
	else if ( res_name == "TRP" ) return  5;
	else if ( res_name == "CYS" ) return  6;
	else if ( res_name == "GLY" ) return  7;
	else if ( res_name == "PHE" ) return  8;
	else if ( res_name == "GLN" ) return  9;
	else if ( res_name == "SER" ) return 10;
	else if ( res_name == "ASN" ) return 11;
	else if ( res_name == "LEU" ) return 12;
	else if ( res_name == "VAL" ) return 13;
	else if ( res_name == "TYR" ) return 14;
	else if ( res_name == "GLU" ) return 15;
	else if ( res_name == "ARG" ) return 16;
	else if ( res_name == "THR" ) return 17;
	else if ( res_name == "ALA" ) return 18;
	else if ( res_name == "MET" ) return 19;
	else if ( res_name == "HIS" ) return 20;

	return INVALID;
}

core::Size DFIRE_Potential::atom_index(
	std::string const & atom_name
) const {
	if      ( atom_name ==  "N" ) return 1;
	else if ( atom_name == "CA" ) return 2;
	else if ( atom_name ==  "C" ) return 3;
	else if ( atom_name ==  "O" ) return 4;
	else if ( atom_name == "CB" ) return 5;

	return INVALID;
}

core::Real
DFIRE_Potential::eval_dfire_pair_energy(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2
) const {

	if ( rsd1.seqpos() == rsd2.seqpos() ) return 0.0;

	using core::Size;
	using std::string;
	Size const natoms1( rsd1.natoms() );
	Size const natoms2( rsd2.natoms() );

	Real score = 0;
	//std::cout << "calling score with " << natoms1 << "," << natoms2 << std::endl;
	for ( Size ii = 1; ii <= natoms1; ++ii ) {
		for ( Size jj = 1; jj <= natoms2; ++jj ) {
			Real const dist = rsd1.xyz(ii).distance( rsd2.xyz(jj) );
			Size const dist_bin_idx = (Size) (2*dist + 0.5);
			//std::cout << "bin(" << dist << ") = " << dist_bin_idx << std::endl;
			if ( dist_bin_idx < 30 ) {
				string const atom_id1( rsd1.type().atom_name(ii) );
				string const atom_id2( rsd2.type().atom_name(jj) );
				string const res_id1 ( rsd1.type().name3() );
				string const res_id2 ( rsd2.type().name3() );
				string const joint_bin_id( get_joint_id( res_id1, atom_id1, res_id2, atom_id2 ) );
				boost::unordered_map< string, Size >::const_iterator it = atom_res_idx_.find(joint_bin_id);
				utility::vector1< Real > const & scores( potential_[it->second] );;
				if ( it != atom_res_idx_.end() && dist_bin_idx <= scores.size() ) {
					score += scores[dist_bin_idx];
					//for ( Size kk = 1; kk <= scores.size(); ++kk ) {
					// std::cout << " " << scores[kk];
					//}
					//std::cout << std::endl;
					//std::cout << "score(" << joint_bin_id << "," << dist << "," << dist_bin_idx <<  ") = " << score << std::endl;
				}
			}
		}
	}
	return score/100;
}

bool DFIRE_Potential::is_loaded() const {
	return potential_is_loaded_;
}

DFIRE_Potential & get_DFIRE_potential() {
	static DFIRE_Potential potential;
	if ( !potential.is_loaded() ) {
		potential.read_potential("scoring/dfire/dfire_pair.lib.txt");
	}
	return potential;
}

} // database
} // methods
} // scoring
} // core
