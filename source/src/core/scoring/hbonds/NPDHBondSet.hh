// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbonds/NPDHBondSet.hh
/// @brief  Container for the non-pairwise-decomposable hydrogen-bond set
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_hbonds_NPDHBondSet_hh
#define INCLUDED_core_scoring_hbonds_NPDHBondSet_hh


// Unit Headers
#include <core/scoring/hbonds/NPDHBondSet.fwd.hh>

// Package headers
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
#include <core/scoring/EnergyMap.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <basic/datacache/CacheableData.hh>


// Utility headers
//#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

// C++
#include <map> // what is the right header for std::pair ?

#include <core/id/AtomID.fwd.hh>

//#include <set>




#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace hbonds {


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// NPDHBondSet
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief A class that holds Hbond objects and helps setup Hbonds for scoring
///
/// @details For general hydrogen bond information, either use the default or option constructor,
/// then use the fill methods in hbonds.hh OR use the convenience constructors to detect all Hbonds.
/// Use the copy constructors to fill NPDHBondSets with the Hydrogen bonds you are interested in.
///
///
class NPDHBondSet : public HBondSet {

	typedef id::AtomID AtomID;
	typedef HBondSet parent;

public:
	NPDHBondSet();
	~NPDHBondSet();

	NPDHBondSet( Size const nres );
	NPDHBondSet( HBondOptions const & options );
	NPDHBondSet( HBondOptions const & options, Size const nres);


	////////////////////////////////////////////////////////////////////////////
	// Convenience Constructors - Use fill methods in hbond.hh for more options
	//

	/// @brief convenience constructor: Find all the hbonds in the pose. BB only default.
	NPDHBondSet(
		pose::Pose & pose,
		EnergyMap const & weights
	);

	/// @brief convenience constructor: Find all the hbonds in the pose. BB only default.
	NPDHBondSet(
		HBondOptions const & options,
		pose::Pose & pose,
		EnergyMap const & wights
	);

	////////////////////////////////////////////////////////////////////////////
	// Copy Constructors
	// These will have to be thought about more thuroughly; does it mae sense in an NPD context
	// to exclude hydrogen bonds to parts of the pose outside of the desired section?

	NPDHBondSet( NPDHBondSet const & src );

	NPDHBondSet( NPDHBondSet const & src, utility::vector1< core::Size > const & exclude_list );

	NPDHBondSet( NPDHBondSet const & src, utility::vector1< bool > const & residue_mask );

	NPDHBondSet( NPDHBondSet const & src, Size seqpos );


public:

	void
	setup_for_residue_pair_energies(
		pose::Pose const & pose,
		bool const calculate_derivative = false
	);

	////////////////////////////////////////////////////////////////////////////

	/// @brief Clone this object
	basic::datacache::CacheableDataOP
	clone() const;

	friend
	std::ostream &
	operator<< ( std::ostream & out, const NPDHBondSet & hbond_set );

	/// @brief equality operator
	friend
	bool
	operator==(NPDHBondSet const & a, NPDHBondSet const & b);

	////////////////////////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////////////////////////
	// Show Methods
	//

	/// @brief Print just the information stored in each individual
	/// hbond.
	void show(std::ostream & out) const;

	// PyRosetta friendly version
	void show() const { show(std::cout); };


	/// @brief Print nicely formated summary of the hbonds and their geometry in the pose.
	void show(
		pose::Pose const & pose,
		bool const print_header,
		std::ostream & out) const;

	// PyRosetta friendly version
	void show(
		pose::Pose const & pose,
		bool const print_header=true) const { show(pose, print_header, std::cout); }


	/// @brief Print nicely formated summary of all the hbonds to a
	/// specific residue
	void show(
		pose::Pose const & pose,
		Size const residue,
		bool const print_header,
		std::ostream & out) const;

	// PyRosetta friendly version
	void show(
		pose::Pose const & pose,
		Size const residue,
		bool const print_header=true) const { show(pose, residue, print_header, std::cout); }

	/// @brief Return the non-pairwise-decomposable weight for a particular atom for a particular
	/// hbond (by index for that atom -- stored in HBond::acc_index() or HBond::don_index())
	Real hbond_weight( id::AtomID const & id, Size hbond_index ) const;

	/// @brief Return the derivative of the non-pairwise-decomposable weight for a particular atom for a particular
	/// hbond (by index for that atom ) as a function of the change in energy of another (possibly the same) hbond (by index)
	/// given the AtomID of the acceptor or the donor.
	Real d_hbond_weight_dE( id::AtomID const & id, Size hbond_fixed_index, Size hbond_changing_index ) const;


	/// @brief Return the derivative for the system induced by a change in energy of a particular hbond
	/// as caused by movement by its atoms
	Real dEtot_dEhb( Size hbond_index ) const;


private:

	void derive_per_hbond_donor_and_acceptor_weights( pose::Pose const & pose );

	void
	get_weights_for_one_partner_hbonder( id::AtomID jj_id );

	void
	get_weights_for_two_partner_hbonder( id::AtomID jj_id );

	void
	sum_dEtot_dEhb_for_atom( id::AtomID const & jj_id );

private:

	/// @brief The weights (in particular, for the npd_hbond_* terms) that will be used to
	/// determine the donor and acceptor weights
	EnergyMap sfxn_weights_;

	/// @brief One (unweighted) energy per hbond per atom
	id::AtomID_Map< utility::vector1< Real > > atom_hbond_energies_;

	/// @brief One score-function weighted energy per hbond per atom
	id::AtomID_Map< utility::vector1< Real > > atom_hbond_sfxn_wtd_energies_;

	/// @brief One weight per hbond per atom
	id::AtomID_Map< utility::vector1< Real > > atom_hbond_weights_;

	///// @brief hbond energy sums for two-partner hbonders
	//id::AtomID_Map< utility::vector1< Real > > atom_hbond_esum_;

	/// @brief d wt / d E so that this needs to be calculated only once
	/// per hydrogen bond
	id::AtomID_Map< utility::vector1< utility::vector1< Real > > > dwt_dE_for_hbond_;

	utility::vector1< Real > dEtot_dEhb_;
	//utility::vector1< Real > don_dEtot_dEhb_;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

/// @brief Given a vector of energies for hydrogen bonds for a particular atom, compute the weights
/// assigned to them. The residue need not have been seen by the NPDHBondSet prior to this function
/// evaluation.
void
weights_for_hbonds(
	conformation::Residue const & res,
	Size atom,
	utility::vector1< Real > const & energies,
	utility::vector1< Real > & weights
);

///// @brief Given a vector of energies for hydrogen bonds for a particular atom, compute the weights
///// assigned to them, and the change in the weights wrt a change in the hbonds' energies.
///// The residue need not have been seen by the NPDHBondSet prior to this function
///// evaluation.
//void
//weights_for_hbonds(
//	conformation::Residue const & res,
//	Size atom,
//	utility::vector1< Real > const & energies,
//	utility::vector1< Real > & weights,
//	utility::vector1< utility::vector1< Real > > & d_wt_dE
//);

void
get_weights_for_one_partner_hbonder(
 	utility::vector1< Real > const & energies,
	utility::vector1< Real > & weights,
	utility::vector1< utility::vector1< Real > > & dwt_dE
);

/// @brief Helper function used by the weights_for_hbonds function above.
Real
get_weights_for_one_partner_hbonder(
 	utility::vector1< Real > const & energies,
	utility::vector1< Real > & weights
);

void
get_weights_for_two_partner_hbonder(
 	utility::vector1< Real > const & energies,
	utility::vector1< Real > & weights,
	utility::vector1< utility::vector1< Real > > & dwt_dE
);

/// @brief Helper function used by the weights_for_hbonds function above.
void
get_weights_for_two_partner_hbonder(
 	utility::vector1< Real > const & energies,
	utility::vector1< Real > & weights
);


} // namespace hbonds
} // namespace scoring
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_hbonds_NPDHBondSet )
#endif // SERIALIZATION


#endif // INCLUDED_core_scoring_ScoreFunction_HH
