// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author


#ifndef INCLUDED_core_scoring_NeighborList_hh
#define INCLUDED_core_scoring_NeighborList_hh

// Unit Headers
#include <core/scoring/NeighborList.fwd.hh>

// Package headers

// AUTO-REMOVED #include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/scoring/etable/EtableEnergy.fwd.hh>


// Project Headers
#include <core/id/AtomID.hh>
// AUTO-REMOVED #include <core/kinematics/DomainMap.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <numeric/xyzVector.hh>


// Utility Headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/kinematics/DomainMap.fwd.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.hh>



namespace core {
namespace scoring {

	/// an atom-atom neighborlist object

	/**
		 The neighborlist is used during minimization to speed atom-atom energy
		 calculations. It stores a list of potentially interacting neighbor atoms
		 for each atom in the system.

		 The logic for using the nblist is tricky.

		 Tentative scheme:
				turn on nblist scoring at start of minimization

				// at this point, want pose to be fully scored
				// so perhaps a call to scorefxn(pose) ??
				// Real const start_score( scorefxn( pose ) );

				pose.energies().setup_use_nblist( true );

				Real const start_func( func( vars ) ); // nblist setup inside this call

				now require that all energy evaluations have an identical set of moving
				dofs (guaranteed if all energy calculations are inside function
				evaluations). This is checked inside each scorecaln using the
				nblist.

				when using the nblist, the rsd-rsd neighbor information is not
				updated. This will probably be a good thing in that it will smooth
				the energy landscape during minimization...

				in a nblist score calculation, we do two things: recover cached
				energies for non-pair-moved positions, and get atom-atom energies
				for the pairs that are on the nblist. We don't cache 2d energies
				for moving positions, since we are not looping over rsd nbr links
				for that score calculation so the caching would be pretty time-
				consuming I think.

				The nblist has the count_pair weights stored, so no calls to
				count_pair !!

				turn off nblist scoring at the end of minimization. Since we have not
				been updating rsd-pair energies for moving pairs, and the rsd-rsd
				nblist is potentially out of data, we reset the neighborgraph at this
				point to ensure a complete score calculation next time.

	**/

// move to separate file
class AtomNeighbor
{
public:
	AtomNeighbor() : rsd_( 0 ), atomno_( 0 ), path_dist_( 0 ), weight_( 0.0 ), weight_func_( 0.0 ) {}

	AtomNeighbor(
		int const rsd_in,
		int const atomno_in,
		Size const path_dist_in,
		Real const weight_in,
		Real const weight_func_in = 1
	):
		rsd_(rsd_in),
		atomno_(atomno_in),
		path_dist_(path_dist_in),
		weight_(weight_in),
		weight_func_(weight_func_in)
	{}


	///
	int
	rsd() const
	{
		return rsd_;
	}

	///
	int
	atomno() const
	{
		return atomno_;
	}

	Size
	path_dist() const
	{
		return path_dist_;
	}

	///
	Real
	weight() const
	{
		return weight_;
	}

	///fpd
	Real
	weight_func() const
	{
		return weight_func_;
	}


	Real & temp1() const { return temp1_; }
	Real & temp2() const { return temp2_; }
	Real & temp3() const { return temp3_; }
	Real & temp4() const { return temp4_; }


private:
	int rsd_;
	int atomno_;
	Size path_dist_;
	Real weight_;
	Real weight_func_;  ///fpd

	mutable Real temp1_;
	mutable Real temp2_;
	mutable Real temp3_;
	mutable Real temp4_;

};

typedef utility::vector1< AtomNeighbor > AtomNeighbors;

///////////////////////////////////////////////////////////////////////////////
class NeighborList : public utility::pointer::ReferenceCount
{
public:
	NeighborList(
		kinematics::DomainMap const & domain_map,
		Real const XX_cutoff,
		Real const XH_cutoff,
		Real const HH_cutoff
	);

	virtual
	~NeighborList();

	// Clone method used in copy ctors of classes that contain NeighborListOP's
	// like, for example, Energies.
	NeighborListOP clone() const { return NeighborListOP( new NeighborList( *this ) ); }

	///
	AtomNeighbors const &
	atom_neighbors(
		int const pos,
		int const atomno
	) const
	{
		return nblist_[ pos ][ atomno ];
	}

	///
	AtomNeighbors const &
	atom_neighbors(
		id::AtomID const & id
	) const
	{
		return nblist_[ id.rsd() ][ id.atomno() ];
	}

	AtomNeighbors const &
	upper_atom_neighbors(
		int const pos,
		int const atomno
	) const
	{
		return upper_nblist_[ pos ][ atomno ];
	}

	///
	AtomNeighbors const &
	upper_atom_neighbors(
		id::AtomID const & id
	) const
	{
		return upper_nblist_[ id.rsd() ][ id.atomno() ];
	}

	AtomNeighbors const &
	intrares_upper_atom_neighbors(
		int const pos,
		int const atomno
	) const
	{
		return intrares_upper_nblist_[ pos ][ atomno ];
	}

	///
	AtomNeighbors const &
	intrares_upper_atom_neighbors(
		id::AtomID const & id
	) const
	{
		return intrares_upper_nblist_[ id.rsd() ][ id.atomno() ];
	}



	/// @brief Initialize the nblist so that it reflects the current coordinates in the pose.
	template < class T_Etable >
	void
	setup(
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		T_Etable const & etable_method
	) const;

	///
	void
	check_domain_map(
		kinematics::DomainMap const & domain_map_in
	) const;


	///
	void
	clear()
	{
		nblist_.clear();
	}

	///
	kinematics::DomainMap const &
	domain_map() const
	{
		return domain_map_;
	}

	/// @brief If auto_update_, ensure that no atom in the pose has not moved too much
	/// since the last time the neighborlist was updated.  The neighborlist
	/// tracks the starting coords for all atoms, and then updates
	template < class T_Etable >
	void
	prepare_for_scoring(
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		T_Etable const & etable_method
	) const;

	void
	set_auto_update( Distance move_tolerance );

	void
	disable_auto_update();

private:

	void
	update_from_wide_nblist( pose::Pose const & pose ) const;

	inline
	DistanceSquared
	atom_pair_cutoff( bool atom1_is_hydrogen, bool atom2_is_hydrogen ) const
	{
		return ( ( atom1_is_hydrogen && atom2_is_hydrogen ) ?
			HH_cutoff_ : ( ( atom1_is_hydrogen || atom2_is_hydrogen ) ?
			XH_cutoff_ : XX_cutoff_ ) );
	}

protected:

	void
	declare_atoms_neighbors( id::AtomID at1, id::AtomID at2, Size path_dist, Real weight, Real weight_func = 1.0 ) const;

	void
	declare_atom_neighbor_1sided( id::AtomID at1, id::AtomID at2, Size path_dist, Real weight, Real weight_func = 1.0 ) const;

private:

	bool auto_update_;

	//square of how far can any atom move before the nblist needs to be updated
	DistanceSquared move_tolerance_sqr_;

	// How far out should the wide nblist reach beyond XX_cutoff_
	Distance const wide_nblist_extension_;

	// square of the coordinate movement before the wide_nblist's data becomes stale
	// == ( wide_nblist_extension - sqrt( move_tolerance_) ) ^ 2
	DistanceSquared wide_move_tolerance_sqr_;

	mutable utility::vector1< utility::vector1< AtomNeighbors > > nblist_;
	mutable utility::vector1< utility::vector1< AtomNeighbors > > upper_nblist_;
	mutable utility::vector1< utility::vector1< AtomNeighbors > > intrares_upper_nblist_;
	mutable utility::vector1< utility::vector1< AtomNeighbors > > wide_nblist_;
	mutable utility::vector1< utility::vector1< Vector > > reference_coords_;
	mutable utility::vector1< utility::vector1< Vector > > wide_reference_coords_;

	kinematics::DomainMap const domain_map_;
	DistanceSquared const XX_cutoff_;
	DistanceSquared const XH_cutoff_;
	DistanceSquared const HH_cutoff_;

	Distance const sqrt_XX_cutoff_;
	Distance const sqrt_XH_cutoff_;
	Distance const sqrt_HH_cutoff_;


	/// Separation square distance for atom pairs in kept in the wide neighbor list
	DistanceSquared const XX_cutoff_wide_;
	DistanceSquared const XH_cutoff_wide_;
	DistanceSquared const HH_cutoff_wide_;

	/// Variables for updating the nblist from the wide nblist
	/// don't use v1< bool >: too slow
	mutable utility::vector1< utility::vector1< Size > > atom_needs_update_from_wide_;
	//mutable utility::vector1< utility::vector1< Size > > atom_has_been_updated_from_wide_;
	mutable utility::vector1< id::AtomID > atoms_to_update_;

	mutable Size n_prepare_for_scorings_;
	mutable Size n_update_from_wide_;
	mutable Size n_full_updates_;
};

typedef utility::pointer::shared_ptr< NeighborList > NeighborListOP;

} // namespace scoring
} // namespace core


#endif
