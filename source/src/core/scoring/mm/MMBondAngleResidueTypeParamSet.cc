// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMBondAngleLibary.cc
/// @brief  Class to store bond angle parameters for a set of ResidueTypes
/// @author Colin A. Smith (colin.smith@ucsf.edu)

// Unit headers
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.hh>

// Project headers
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/methods/MMBondAngleEnergy.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers

// C++ headers
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <utility/assert.hh>

#include <core/scoring/mm/MMBondAngleResidueTypeParam.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace mm {

static thread_local basic::Tracer TR( "core.mm.MMBondAngleResidueTypeParamSet" );

MMBondAngleResidueTypeParamSet::MMBondAngleResidueTypeParamSet() :
	mm_bondangle_library_( &scoring::ScoringManager::get_instance()->get_MMBondAngleLibrary() ),
	use_residue_type_theta0_( false )
{ }

/// @brief copy ctor
MMBondAngleResidueTypeParamSet::MMBondAngleResidueTypeParamSet(
	MMBondAngleResidueTypeParamSet const & src
) :
	utility::pointer::ReferenceCount( src )
{
	*this = src;
}

MMBondAngleResidueTypeParamSet::~MMBondAngleResidueTypeParamSet() {}

/// @brief lookup a param object for a given ResidueType
MMBondAngleResidueTypeParam const &
MMBondAngleResidueTypeParamSet::get(
	core::chemical::ResidueType const & residue_type
)
{
	std::map<std::string, MMBondAngleResidueTypeParam>::iterator iter(reside_type_param_map_.find(residue_type.name()));

	if ( iter == reside_type_param_map_.end() ) {

		MMBondAngleResidueTypeParam new_param;
		new_param.init(residue_type, *mm_bondangle_library_, use_residue_type_theta0_, central_atoms_to_score_);
		reside_type_param_map_[residue_type.name()] = new_param;

		return reside_type_param_map_[residue_type.name()];
	}

	return iter->second;
}

/// @brief lookup a param object for a given ResidueType, does not auto-create
MMBondAngleResidueTypeParam const *
MMBondAngleResidueTypeParamSet::get(
	core::chemical::ResidueType const & residue_type
) const
{
	std::map<std::string, MMBondAngleResidueTypeParam>::const_iterator iter(reside_type_param_map_.find(residue_type.name()));

	if ( iter == reside_type_param_map_.end() ) {

		return NULL;
	}

	return &(iter->second);
}

/// @brief get the indices of the connections that connects one of my atoms to an atom on another residue
bool
connection_indices(
	core::conformation::Residue const & residue,
	core::conformation::Residue const & other_residue,
	core::Size const my_atomno,
	core::Size const other_atomno,
	core::Size & my_connection,
	core::Size & other_connection
)
{
	for ( core::Size i = 1; i <= residue.n_residue_connections(); ++i ) {
		if ( residue.residue_connection(i).atomno() == (signed)my_atomno ) {
			if ( residue.connect_map(i).resid() == other_residue.seqpos() ) {
				my_connection = i;
				other_connection = residue.connect_map(i).connid();
				if ( other_residue.residue_connection(other_connection).atomno() == (signed)other_atomno ) {
					return true;
				}
			}
		}
	}

	// didn't find a connection between the two atoms
	my_connection = 0;
	other_connection = 0;
	return false;
}

/// @brief lookup Ktheta and theta0 for any bond angle in a conformation
void
MMBondAngleResidueTypeParamSet::lookup(
	core::conformation::Conformation const & conformation,
	core::id::AtomID const & atomid1,
	core::id::AtomID const & atomid2,
	core::id::AtomID const & atomid3,
	core::Real & Ktheta,
	core::Real & theta0
)
{
	MMBondAngleResidueTypeParam const & residue_type_param(get(conformation.residue_type(atomid2.rsd())));

	lookup(conformation, atomid1, atomid2, atomid3, residue_type_param, Ktheta, theta0);
}


/// @brief lookup Ktheta and theta0 for any bond angle in a conformation
void
MMBondAngleResidueTypeParamSet::lookup(
	core::conformation::Conformation const & conformation,
	core::id::AtomID const & atomid1,
	core::id::AtomID const & atomid2,
	core::id::AtomID const & atomid3,
	core::Real & Ktheta,
	core::Real & theta0
) const
{
	MMBondAngleResidueTypeParam const * residue_type_param(get(conformation.residue_type(atomid2.rsd())));

	if ( residue_type_param ) {
		lookup(conformation, atomid1, atomid2, atomid3, *residue_type_param, Ktheta, theta0);
	} else {
		// this is slow, but shouln't be called after scoring caches all MMBondAngleResidueTypeParam objects
		MMBondAngleResidueTypeParam new_param;
		new_param.init(conformation.residue_type(atomid2.rsd()), *mm_bondangle_library_, use_residue_type_theta0_, central_atoms_to_score_);
		lookup(conformation, atomid1, atomid2, atomid3, new_param, Ktheta, theta0);
	}
}

/// @brief lookup Ktheta and theta0 for any bond angle in a conformation
void
MMBondAngleResidueTypeParamSet::lookup(
	core::conformation::Conformation const & conformation,
	core::id::AtomID const & atomid1,
	core::id::AtomID const & atomid2,
	core::id::AtomID const & atomid3,
	MMBondAngleResidueTypeParam const & residue_type_param,
	core::Real & Ktheta,
	core::Real & theta0
) const
{
	Ktheta = 0;
	theta0 = 0;

	if ( atomid1.rsd() == atomid2.rsd() && atomid1.rsd() == atomid3.rsd() ) {

		// simple case, all atoms are in the same residue
		three_atom_set atom_set(atomid1.atomno(), atomid2.atomno(), atomid3.atomno());
		core::Size bondangle_index(residue_type_param.bondangle_index(atom_set));

		if ( bondangle_index != 0 ) {
			// use residue_type_param data
			Ktheta = residue_type_param.Ktheta(bondangle_index);
			theta0 = residue_type_param.theta0(bondangle_index);
		}

		return;
	}

	int const mmtype1(conformation.residue_type(atomid1.rsd()).atom(atomid1.atomno()).mm_atom_type_index());
	int const mmtype2(conformation.residue_type(atomid2.rsd()).atom(atomid2.atomno()).mm_atom_type_index());
	int const mmtype3(conformation.residue_type(atomid3.rsd()).atom(atomid3.atomno()).mm_atom_type_index());

	mm_bondangle_library_citer_pair params(mm_bondangle_library_->lookup(mmtype1, mmtype2, mmtype3));

	Ktheta = (params.first->second).key1();
	theta0 = (params.first->second).key2();

	// the rest is all to determine if theta0 should be derived from the ResidueType geometry

	core::Size primary_rsd(atomid2.rsd());
	core::Size secondary_rsd(0);
	core::Size atomno1(0);
	core::Size atomno2(atomid2.atomno());
	core::Size atomno3(0);

	if ( atomid1.rsd() == atomid2.rsd() ) {
		// more complicated case, an interresidue bond angle with atoms 1 & 2 in the same residue
		secondary_rsd = atomid3.rsd();
		atomno1 = atomid1.atomno();
		atomno3 = atomid3.atomno();
	} else if ( atomid2.rsd() == atomid3.rsd() ) {
		// more complicated case, an interresidue bond angle with atoms 2 & 3 in the same residue
		secondary_rsd = atomid1.rsd();
		atomno1 = atomid3.atomno();
		atomno3 = atomid1.atomno();
	} else {
		// unsupported case, an interresidue bond angle with atom 2 in its own residue distinct from 1 & 3
		utility_exit_with_message("Lookup of three residue bond angle parameters not supported");
	}

	core::Size primary_connection;
	core::Size secondary_connection;
	// lookup the residue connection ids
	if ( connection_indices(conformation.residue(primary_rsd), conformation.residue(secondary_rsd), atomno2, atomno3,
			primary_connection, secondary_connection) ) {
		two_atom_set const primary_atom_set(atomno1, atomno2);
		core::Size connection_index(residue_type_param.connection_index(primary_connection, primary_atom_set));
		if ( !connection_index ) {
			Ktheta = 0;
			return;
		}
		if ( residue_type_param.connection_use_theta0(primary_connection, connection_index) ) {
			theta0 = residue_type_param.connection_theta0(primary_connection, connection_index);
		}
	}
}

MMBondAngleResidueTypeParamSetCOP
mm_bond_angle_residue_type_param_set(
	core::scoring::ScoreFunction const & scorefxn
)
{
	for ( core::scoring::ScoreFunction::TWO_B_MethodIterator iter(scorefxn.ci_2b_intrares_begin());
			iter != scorefxn.ci_2b_intrares_end(); ++iter ) {

		core::scoring::methods::MMBondAngleEnergyCOP bond_angle_energy_method =
			utility::pointer::dynamic_pointer_cast<core::scoring::methods::MMBondAngleEnergy const>( *iter );

		if ( bond_angle_energy_method ) {
			return bond_angle_energy_method->residue_type_param_set();
		}
	}

	return NULL;
}


} // namespace mm
} // namespace scoring
} // namespace core
