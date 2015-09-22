// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMBondAngleResidueTypeParam.cc
/// @brief  Class to store bond angle parameters for a particular ResidueType
/// @author Colin A. Smith (colin.smith@ucsf.edu)

// Unit headers
#include <core/scoring/mm/MMBondAngleResidueTypeParam.hh>

// Project headers
#include <core/chemical/Atom.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.hh>
#include <core/scoring/mm/MMBondAngleLibrary.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>

// Utility header
#include <numeric/conversions.hh>
#include <utility/string_util.hh>

#include <core/chemical/MMAtomType.hh>
#include <utility/vector1.hh>


// C++ headers

static THREAD_LOCAL basic::Tracer TR( "core.mm.MMBondAngleResidueTypeParam" );


namespace core {
namespace scoring {
namespace mm {

MMBondAngleResidueTypeParam::MMBondAngleResidueTypeParam()
{ }

bool
score_atom_centrally(
	core::chemical::ResidueType const & restype,
	utility::vector1<std::string> const & central_atoms_to_score,
	Size atomno
)
{
	if ( central_atoms_to_score.size() == 0 ) return true;

	for ( Size ii = 1; ii <= central_atoms_to_score.size(); ++ii ) {
		if ( utility::same_ignoring_spaces( restype.atom_name( atomno ), central_atoms_to_score[ ii ] ) ) {
			//std::cout << "Counting atom " << atomno << " " << restype.name() << " " << restype.atom_name( atomno ) << std::endl;
			return true;
		}
	}
	return false;
}

void
MMBondAngleResidueTypeParam::init(
	core::chemical::ResidueType const & residue_type,
	MMBondAngleLibrary const & mm_bondangle_library,
	bool use_residue_type_theta0,
	utility::vector1<std::string> const & central_atoms_to_score
)
{
	// add data for intraresidue bond angles
	bondangle_atom_sets_.clear();
	Ktheta_.clear();
	theta0_.clear();

	for ( core::Size i = 1; i <= residue_type.num_bondangles(); ++i ) {

		three_atom_set const & atom_set(residue_type.bondangle(i));
		if ( !score_atom_centrally(residue_type, central_atoms_to_score, atom_set.key2()) ) continue;

		core::Real residue_type_theta0(numeric::angle_radians(residue_type.atom(atom_set.key1()).ideal_xyz(), residue_type.atom(atom_set.key2()).ideal_xyz(), residue_type.atom(atom_set.key3()).ideal_xyz()));

		std::string const & type1(residue_type.mm_atom_type(atom_set.key1()).name());
		std::string const & type2(residue_type.mm_atom_type(atom_set.key2()).name());
		std::string const & type3(residue_type.mm_atom_type(atom_set.key3()).name());

		//TR << residue_type.atom_name(atom_set.key1()) << "-"
		//  << residue_type.atom_name(atom_set.key2()) << "-"
		//  << residue_type.atom_name(atom_set.key3())
		//   << " (" << type1 << "-" << type2 << "-" << type3 << ")";

		core::scoring::mm::mm_bondangle_library_citer_pair mm_pair(mm_bondangle_library.lookup(type1, type2, type3));

		// make sure there is at least one set of parameters defined
		debug_assert(mm_pair.first != mm_pair.second);

		core::Real mm_Ktheta((mm_pair.first->second).key1());
		core::Real mm_theta0((mm_pair.first->second).key2());

		// make sure there was only one set of parameters defined
		debug_assert(++mm_pair.first == mm_pair.second);

		//TR << " mm_Ktheta: " << mm_Ktheta << " mm_theta0: " << numeric::conversions::degrees(mm_theta0)
		//   << " rt_theta0: " << numeric::conversions::degrees(residue_type_theta0);

		if ( mm_Ktheta ) {
			bondangle_atom_sets_.push_back(atom_set);
			Ktheta_.push_back(mm_Ktheta);
			theta0_.push_back(use_residue_type_theta0 ? residue_type_theta0 : mm_theta0);
		} else {
			//TR << " Ignoring";
		}

		//TR << std::endl;
	}

	// update bondangles_for_atom_ and bondangle_index_
	bondangles_for_atom_.clear();
	bondangles_for_atom_.resize( residue_type.natoms() );
	bondangle_index_.clear();

	for ( core::Size i = 1; i <= bondangle_atom_sets_.size(); ++i ) {

		three_atom_set const & bondangle_atom_set( bondangle_atom_sets_[i] );
		bondangles_for_atom_[bondangle_atom_set.key1()].push_back(i);
		bondangles_for_atom_[bondangle_atom_set.key2()].push_back(i);
		bondangles_for_atom_[bondangle_atom_set.key3()].push_back(i);

		// forward mapping
		bondangle_index_[bondangle_atom_set] = i;
		// reverse mapping
		bondangle_index_[three_atom_set(bondangle_atom_set.key3(), bondangle_atom_set.key2(), bondangle_atom_set.key1())] = i;
	}

	// add data for interresidue bond angles
	connection_atom_sets_.clear();
	connection_theta0_.clear();
	connection_use_theta0_.clear();
	connection_index_.clear();
	connection_atom_sets_.resize(residue_type.n_residue_connections());
	connection_theta0_.resize(residue_type.n_residue_connections());
	connection_use_theta0_.resize(residue_type.n_residue_connections());
	connection_index_.resize(residue_type.n_residue_connections());

	for ( core::Size i = 1; i <= residue_type.n_residue_connections(); ++i ) {

		core::chemical::ResidueConnection const & residue_connection(residue_type.residue_connection(i));
		core::Size const connection_atomno(residue_connection.atomno());
		if ( !score_atom_centrally(residue_type, central_atoms_to_score, connection_atomno) ) {
			continue;
		}
		core::Vector external_xyz(residue_connection.icoor().build(residue_type));
		core::chemical::AtomIndices const & bonded_neighbors(residue_type.bonded_neighbor(connection_atomno));

		for ( core::Size j = 1; j <= bonded_neighbors.size(); ++j ) {

			core::Real residue_type_theta0(numeric::angle_radians(residue_type.atom(bonded_neighbors[j]).ideal_xyz(), residue_type.atom(connection_atomno).ideal_xyz(), external_xyz));

			//std::string const & type1(residue_type.mm_atom_type(bonded_neighbors[j]).name());
			//std::string const & type2(residue_type.mm_atom_type(connection_atomno).name());

			//TR << residue_type.atom_name(bonded_neighbors[j]) << "-"
			//  << residue_type.atom_name(connection_atomno) << "-?"
			//   << " (" << type1 << "-" << type2 << "-?)"
			//   << " rt_theta0: " << numeric::conversions::degrees(residue_type_theta0) << std::endl;

			two_atom_set const connection_atom_set(bonded_neighbors[j], connection_atomno);
			connection_atom_sets_[i].push_back(connection_atom_set);
			connection_theta0_[i].push_back(residue_type_theta0);
			connection_use_theta0_[i].push_back(use_residue_type_theta0);

			connection_index_[i][connection_atom_set] = j;
		}
	}
}

/// @brief stream << MMBondAngleResidueTypeParam
std::ostream &
operator <<(
	std::ostream & os,
	MMBondAngleResidueTypeParam const & residue_type_param
)
{
	os << "Intraresidue Bond Angles:" << std::endl;
	for ( core::Size i = 1; i <= residue_type_param.bondangle_atom_sets_.size(); ++i ) {
		os << residue_type_param.bondangle_atom_sets_[i].key1() << "-"
			<< residue_type_param.bondangle_atom_sets_[i].key2() << "-"
			<< residue_type_param.bondangle_atom_sets_[i].key3()
			<< " Ktheta: " << residue_type_param.Ktheta_[i]
			<< " theta0: " << numeric::conversions::degrees(residue_type_param.theta0_[i]) << std::endl;
	}

	for ( core::Size i = 1; i <= residue_type_param.connection_atom_sets_.size(); ++i ) {
		os << "Connection " << i << " Bond Angles:" << std::endl;
		for ( core::Size j = 1; j <= residue_type_param.connection_atom_sets_[i].size(); ++j ) {
			os << residue_type_param.connection_atom_sets_[i][j].key1() << "-"
				<< residue_type_param.connection_atom_sets_[i][j].key2() << "-?"
				<< " rt_theta0: " << numeric::conversions::degrees(residue_type_param.connection_theta0_[i][j]) << " "
				<< (residue_type_param.connection_use_theta0_[i][j] ? "use_residue_type_theta0" : "use_mm_theta0") << std::endl;
		}
	}

	return os;
}

} // namespace mm
} // namespace scoring
} // namespace core
