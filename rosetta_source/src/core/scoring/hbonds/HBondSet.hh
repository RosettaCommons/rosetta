// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hbonds/HBondSet.hh
/// @brief  Hydrogen bond set class declaration
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_hbonds_HBondSet_hh
#define INCLUDED_core_scoring_hbonds_HBondSet_hh


// Unit Headers
#include <core/scoring/hbonds/HBondSet.fwd.hh>


// Package headers
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
// AUTO-REMOVED #include <core/id/AtomID.hh>
#include <basic/datacache/CacheableData.hh>
// AUTO-REMOVED #include <core/pose/datacache/CacheableDataType.hh>


// Utility headers
//#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

// C++
#include <map> // what is the right header for std::pair ?

#include <core/id/AtomID.fwd.hh>

//#include <set>



namespace core {
namespace scoring {
namespace hbonds {

class HBond : public utility::pointer::ReferenceCount {

public:
	HBond(
		Size const dhatm,
		bool const dhatm_is_protein_backbone,
		bool const dres_is_protein,
		bool const dres_is_dna,
		bool const dhatm_is_backbone,
		Size const dres,
		Size const aatm,
		bool const aatm_is_protein_backbone,
		bool const ares_is_protein,
		bool const ares_is_dna,
		bool const aatm_is_backbone,
		Size const ares,
		HBEvalType const hbe_type,
		Real const energy_in, // unweighted
		Real const weight_in,
		HBondDerivs const & derivs_in
	);

	///
	Size
	don_res() const;

	///
	Size
	don_hatm() const;

	/// needed for silly allow logic
	bool
	don_hatm_is_protein_backbone() const;

	bool
	don_res_is_protein() const;

	bool
	don_res_is_dna() const;

	/// needed for silly allow logic
	bool
	don_hatm_is_backbone() const;

	///
	Size
	acc_res() const;

	///
	Size
	acc_atm() const;

	/// needed for silly allow logic
	bool
	acc_atm_is_protein_backbone() const;

	bool
	acc_res_is_protein() const;

	bool
	acc_res_is_dna() const;

	/// needed for silly allow logic
	bool
	acc_atm_is_backbone() const;

	/// NOTE: this is unweighted energy, see weight() for the weight
	Real
	energy() const;

	///
	Real
	weight() const;

	///
	HBondDerivs const &
	derivs() const;

	///
	HBEvalType const &
	eval_type() const;

	///
	bool
	atom_is_donorH( id::AtomID const & atom ) const;

	///
	bool
	atom_is_acceptor( id::AtomID const & atom ) const;

	///@breif a bare bones description of the data contained in the hbond object
	void
	show( std::ostream & out ) const;

	///@brief a prettier, more interpretable description of an hbond,
	///including pdb identified residues and the geometric dimensions of
	///the hydrogen bond.
	void show(
		pose::Pose const & pose,
		bool const print_header,
		std::ostream & out) const;

	// PyRosetta friendly version
	void show(pose::Pose const & pose, bool const print_header=true) const  { show(pose, print_header, std::cout); }

	friend
	std::ostream &
	operator<< ( std::ostream & out, const HBond & hbond );

	///
	friend
	bool
	operator==(HBond const & a, HBond const & b);

	///////////////////////////////////////////////////////////////////////////////
	static
	bool
	hbond_energy_comparer( HBondCOP a, HBondCOP b );

	///////
	// data
private:

	Size don_hatm_;
	bool don_hatm_is_protein_backbone_;
	bool don_res_is_protein_;
	bool don_res_is_dna_;
	bool don_hatm_is_backbone_;
	Size don_res_;
	Size acc_atm_;
	bool acc_atm_is_protein_backbone_;
	bool acc_res_is_protein_;
	bool acc_res_is_dna_;
	bool acc_atm_is_backbone_;
	Size acc_res_;
	HBEvalType eval_type_;

	Real energy_;
	Real weight_;
	HBondDerivs derivs_;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// HBondSet
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

class HBondSet : public basic::datacache::CacheableData {

public:
	HBondSet();
	~HBondSet();

	HBondSet( Size const nres );
	HBondSet( HBondOptions const & options );
	HBondSet( HBondOptions const & options, Size const nres);

	///@brief convenience constructor: Find all the hbonds in the pose.
	HBondSet(
		pose::Pose & pose);

	///@brief convenience constructor: Find all the hbonds in the pose.
	HBondSet(
		HBondOptions const & options,
		pose::Pose & pose);


	HBondSet( HBondSet const & src );

	HBondSet( HBondSet const & src, utility::vector1< core::Size > exclude_list );

	HBondSet( HBondSet const & src, utility::vector1< bool > residue_mask );

	HBondSet( HBondSet const & src, Size seqpos );

	// typedefs
	typedef id::AtomID AtomID;

public:

	/// @brief Clone this object
	basic::datacache::CacheableDataOP
	clone() const;

	/// \brief  Number of hbonds
	Size
	nhbonds() const;

	// for accessing the nbrs, allows hacky non-tenA_nbr count
	int
	nbrs( Size const seqpos ) const;

	void
	set_nbrs( Size const seqpos, Size value );

	/// \brief  Access hbond
	HBond const &
	hbond( Size const number ) const;

	/// \brief  Add a new hbond to the list
	/// updates the "hbchk" array as necessary
	void
	append_hbond(
		Size const dhatm,
		conformation::Residue const & don_rsd,
		Size const aatm,
		conformation::Residue const & acc_rsd,
		HBEvalType const & hbe_type,
		Real const energy,
		Real const weight,
		HBondDerivs const & deriv
	);

	/// \brief  Is this hbond allowed under the bb-bb exclusion scheme?
	bool
	allow_hbond( Size const index ) const;

	bool
	allow_hbond( HBond const & hbond ) const;

	/// @brief is the backbone bone acceptor group in a bb/bb hydrogen bond?
	bool
	acc_bbg_in_bb_bb_hbond( Size const residue ) const;

	/// @brief is the backbone bone donor group in a bb/bb hydrogen bond?
	bool
	don_bbg_in_bb_bb_hbond( Size const residue ) const;

	/// @Manually set the state of backbone-backbone acceptor. Used for symmetry.
	void
	set_backbone_backbone_acceptor( Size const residue, bool state )
	{
		backbone_backbone_acceptor_[ residue ] = state;
	}

	/// @Manually set the state of backbone-backbone donor. Used for symmetry.
	void
	set_backbone_backbone_donor( Size const residue, bool state )
	{
		backbone_backbone_donor_[ residue ] = state;
	}

	/// @brief  Setup the mapping from atoms to hbonds
	void
	setup_atom_map() const;

	/// \brief  Get a vector of all the hbonds involving this atom
	utility::vector1< HBondCOP > const &
	atom_hbonds( AtomID const & atom ) const;

	/// \brief  Delete all the data
	void
	clear();

	/// \brief Resize bb info arrays
	void
	resize_bb_donor_acceptor_arrays( Size const new_dimension );

	void
	copy_bb_donor_acceptor_arrays( HBondSet const & src );

	///
	void
	setup_for_residue_pair_energies(
		pose::Pose const & pose,
		bool const calculate_derivative = false,
		bool const backbone_only = true
	);

	void
	sort_by_weighted_energy();

	/// @brief Read access to the stored hbond options
	HBondOptions const &
	hbond_options() const;

	/// @breif set the hbond options for this hbond set; clears all hbonds already stored
	void set_hbond_options( HBondOptions const & options );

	friend
	std::ostream &
	operator<< ( std::ostream & out, const HBondSet & hbond_set );

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


	/// @brief equality operator
	friend
	bool
	operator==(HBondSet const & a, HBondSet const & b);

private:
	typedef std::map< AtomID, utility::vector1< HBondCOP > > HBondAtomMap;

private:
	////////
	// data ---- IF YOU ADD DATA ALSO ADD IT TO THE COPY C-TOR

	HBondOptionsCOP options_;

	utility::vector1< HBondOP > hbonds_;
	utility::vector1< bool > backbone_backbone_donor_;
	utility::vector1< bool > backbone_backbone_acceptor_;
	utility::vector1< int > nbrs_;
	mutable HBondAtomMap atom_map_;
	mutable bool atom_map_init_;
	static utility::vector1< HBondCOP > empty_list_of_hbonds_;
};

} // namespace hbonds
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
