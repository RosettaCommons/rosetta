// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/fldsgn/topology/StrandPairing.hh
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

#ifndef INCLUDED_protocols_fldsgn_topology_StrandPairing_hh
#define INCLUDED_protocols_fldsgn_topology_StrandPairing_hh

// Unit headers
#include <protocols/fldsgn/topology/StrandPairing.fwd.hh>

// Project headers
#include <core/types.hh>
#include <protocols/fldsgn/topology/DimerPairing.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>
#include <utility/vector1.hh>
#include <map>
#include <set>
#include <string>


namespace protocols {
namespace fldsgn {
namespace topology {

class StrandPairing : public utility::pointer::ReferenceCount {
public:

	typedef std::string String;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef protocols::fldsgn::topology::SS_Info2_OP SS_Info2_OP;
	typedef protocols::fldsgn::topology::SS_Info2_COP SS_Info2_COP;

public:// construct/destruct


	/// @brief default constructor
	StrandPairing();

	/// @brief value constructor
	StrandPairing(
		Size const s1,
		Size const s2,
		Size const b1,
		Size const b2,
		Size const p,
		Real const rs,
		char const o
	);

	/// @brief value constructor
	StrandPairing(
		Size const s1,
		Size const s2,
		Real const rs,
		char const o
	);

	/// @brief copy constructor
	StrandPairing( String const & spair);

	/// @brief default destructor
	virtual ~StrandPairing();

	/// @brief clone this object
	StrandPairingOP clone();

	/// @brief return strand pairing
	friend
	std::ostream & operator<<(std::ostream & out, const StrandPairing &sp);


public: //accessors


	/// @brief the strand number of the 1st strand in strand pairing
	inline Size s1() const { return s1_; }

	/// @brief the strand number of the 2nd strand in strand pairing
	inline Size s2() const { return s2_; }

	/// @brief the residue number of the beginning of 1st strand
	inline Size begin1() const { return begin1_; }

	/// @brief the residue number of the end of 1st strand
	inline Size end1() const { return end1_; }

	/// @brief the residue number of the beginning of 2nd strand
	inline Size begin2() const { return begin2_; }

	/// @brief the residue number of the end of 2nd strand
	inline Size end2() const { return end2_; }

	/// @brief the number of register shift between the strands
	inline Real rgstr_shift() const { return rgstr_shift_; }

	/// @brief the pleating at the begining of strand_pairing
	inline utility::vector1< Size > pleats1() const { return pleats1_; }

	/// @brief the pleating at the end of strand_pairing
	inline utility::vector1< Size > pleats2() const { return pleats2_; }

	/// @brief orientation, parallel or anti-parallel, of strand pairing
	inline char orient() const { return orient_; }

	/// @brief whether the strand pairing have bulge or not
	inline bool has_bulge() const { return ( !bulges1_.empty() ) || ( !bulges2_.empty() ); }

	/// @brief whethe the given residue number is a bulge
	bool is_bulge( Size const resid ) const { return ( bulges1_.find( resid ) != bulges1_.end() ) || ( bulges2_.find( resid ) != bulges2_.end() ); }

	/// @brief StrandPairing is descripbed as s1()-s2().orient().rgstr_shift()
	/// For example, 2-3.A.1 means 2nd and 3rd strands make anti-parallel strand_pairing with register shift 1
	inline String name() const { return name_; }

	/// @brief return residue pairing
	bool has_paired_residue( Size const res ) const;

	/// @brief residue pair
	Size residue_pair( Size const res ) const;


public:


	/// @brief elongate strand pairings
	bool elongate( Size const r1, Size const r2, Size const p1, Size const p2 );

	/// @brief
	bool add_pair( Size const r1, Size const r2, char const orient, Real const rgstr );

	/// @brief the length of 1st strand
	Size size1() const;

	/// @brief the length of 2nd strand
	Size size2() const;

	/// @brief is parallel
	bool is_parallel() const;

	/// @brief whether input residue is included in this StrandPairinge or not
	bool is_member( Size const res );

	/// @brief reset begin1_, begin2_, and end1_, end2_ based on ssinfo
	/// @detailed abego is used for obtaining bulge information if it is non-empty
	void redefine_begin_end( SS_Info2_COP const ss_info, utility::vector1< String > const & abego );


private: // initialize


	void initialize();


private:  // data


	/// @brief Strand number of first strand in the strand pair
	Size s1_;

	/// @brief Strand number of second strand in the strand pair
	Size s2_;

	/// @brief end resides of first and second strands
	Size begin1_, end1_, begin2_, end2_;

	/// @brief pleats of end residues
	utility::vector1< Size > pleats1_, pleats2_;

	/// @brief register shift between two strands
	Real rgstr_shift_;

	/// @brief two strands make a sheet by parallel, "P", anti parallel, "A", and
	/// if not defined, "N"
	char orient_;

	/// @brief strand_pairing as in the style: s1_-s2_.orient_.rgstr_shift_
	String name_;

	/// @brief residue pair
	std::map< Size, Size > residue_pair_;

	/// @brief list of bulges in strand 1
	std::set< Size > bulges1_;

	/// @brief list of bulges in strand 2
	std::set< Size > bulges2_;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////
class StrandPairingSet : public utility::pointer::ReferenceCount {
public: // typedef


	typedef std::string String;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef utility::vector1< Size > VecSize;
	typedef protocols::fldsgn::topology::SS_Info2 SS_Info2;
	typedef protocols::fldsgn::topology::SS_Info2_COP SS_Info2_COP;
	typedef protocols::fldsgn::topology::DimerPairing DimerPairing;
	typedef protocols::fldsgn::topology::DimerPairings DimerPairings;


public:// construct/destruct


	/// @brief default constructor
	StrandPairingSet();

	/// @brief value constructor
	StrandPairingSet( StrandPairings const & strand_pairings );

	/// @brief value constructor
	StrandPairingSet( SS_Info2 const & ssinfo, DimerPairings const & dimer_pairs );

	/// @brief value constructor
	StrandPairingSet( String const & spairstring, SS_Info2_COP const ssinfo = NULL );

	/// @brief value constructor
	StrandPairingSet( String const & spairstring, SS_Info2_COP const ssinfo, utility::vector1< String > const & abego );

	/// @brief default destructor
	virtual ~StrandPairingSet();

	/// @brief clone this object
	StrandPairingSetOP clone() const;

	/// @brief return strand pairing
	friend std::ostream & operator<<( std::ostream & out, const StrandPairingSet &s );


public: // mutators


	/// @brief add StrandPairingOP to StandPairingSet
	void push_back( StrandPairingOP const sop );

	/// @brief add StrandPairingOP to StandPairingSet
	void push_back_and_finalize( StrandPairingOP const sop );

	/// @brief clear data of this StrandPairingSet
	void clear();


public: // accessors


	/// @brief return begin of iterator of strand_pairings_
	StrandPairings::const_iterator begin() { return strand_pairings_.begin(); }

	/// @brief return end of iterator of strand_pairings_
	StrandPairings::const_iterator end() { return strand_pairings_.end(); }

	/// @brief
	Size size() const;


public: //  accessors

	/// @brief
	bool finalized() const { return finalized_; }

	/// @brief return all strand pairings
	Size num_strands() const { return num_strands_; };

	/// @brief return all strand pairings
	StrandPairings const & strand_pairings() const;

	/// @brief return one of the strand_pairings given a number
	StrandPairingOP strand_pairing( Size const s ) const;

	/// @brief return the pointer of the strand pairing, given the two strand numbers of s1, and s2
	/// if s1 and s2 does not make pairing, return 0
	StrandPairingOP strand_pairing( Size const s1, Size const s2 ) const;

	/// @brief return strand number of neighbor strands of a input strand
	VecSize const & neighbor_strands( Size const s ) const;


public:


	/// @brief the name of StrandPairingSet is expressed by the combination of strand pairings
	/// For example, 2kl8 of ferredoxin-like fold is described as 1-3.A.0;2-3.A.0;1-4.A.0
	String name() const;

	/// @brief the name of StrandPairingSet without register shift
	/// For example, 2kl8 of ferredoxin-like fold is described as 1-3.A;2-3.A;1-4.A
	String name_wo_rgstr() const;


public:


	/// @brief make the number of strand pairing as two
	void make_strand_neighbor_two();

	/// @brief remove a set of strand pairings from datay
	void drop_strand_pairs( StrandPairings const & drop_spairs );

	/// @brief finalize this and create_map_strand_pairings
	void finalize();

protected:
	void initialize_by_sspair_string( String const & spairstring, SS_Info2_COP const ssinfo );
	void initialize_by_sspair_string_and_abego(
		String const & spairstring,
		SS_Info2_COP const ssinfo,
		utility::vector1< String > const & abego );

private:
	void initialize_by_dimer_pairs( SS_Info2 const & ssinfo, DimerPairings const & dimer_pairs );


private:// data


	/// @brief vector1 including owning pointers of StrandPairing
	StrandPairings strand_pairings_;

	/// @brief the name of StrandPairingSet is expressed by the combination of strand pairings
	/// For example, 2kl8 of ferredoxin-like fold is given as 1-3.A.0;2-3.A.0;1-4.A.0
	String spairset_name_;

	/// @brief the total number of strands included in StrandPairingSet
	Size num_strands_;

	/// @brief whether the map_strand_pairings_ is initialized or not
	bool finalized_;

	/// @brief 2D table of the pointers of strand pairings, which are sorted by the strand number
	utility::vector1< utility::vector1< StrandPairingOP > > map_strand_pairings_;

	// @brief neighbor strasnd
	mutable std::map< Size, VecSize > neighbor_strands_;

	// @brief
	StrandPairingOP empty_;


}; // StrandPairingSet

} // namespace topology
} // namespace fldsgn
} // namespace protocols

#endif
