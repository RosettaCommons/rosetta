// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// @brief
// @author olange: ported from original bblum-rosetta++ version $


#ifndef INCLUDED_core_scoring_dssp_StrandPairing_HH
#define INCLUDED_core_scoring_dssp_StrandPairing_HH

// Unit Headers
#include <core/scoring/dssp/StrandPairing.fwd.hh>
#include <core/scoring/dssp/PairingsList.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>


// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2.fwd.hh>

// C++ headers
#include <string>

#include <utility/exit.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace dssp {

//////////////////////////////////////////////////////////////////
// StrandPairing
// This class encapsulates a pairing between two beta strands.
// It is designed to be initialized with a single residue-residue
// pairing and then extended, residue by residue, along the
// strands.  It MUST be extended sequentially along the length
// of the strands.  It can accomodate beta bulges of length at
// most 1 on one strand and 4 on the other (inspired by dssp
// guidelines).
//////////////////////////////////////////////////////////////////
class StrandPairing {
	typedef utility::vector1< core::Size > SizeList;
public:
	StrandPairing();
	StrandPairing(core::Size res1, core::Size res2, bool antiparallel, core::Size pleating);
	~StrandPairing();

	core::Size operator==(const StrandPairing &rhs) const;
	// core::Size operator!=(const StrandPairing &rhs) const;
	core::Size operator<(const StrandPairing &rhs) const;

	core::Size size() const {
		return end1_ - begin1_ + 1;
	}

	core::Size begin1() const {
		return begin1_;
	}

	core::Size end1() const {
		return end1_;
	}

	bool range_check() const;

	core::Size contact_order() const;

	core::Size get_register() const;
	void get_all_register_and_bulges( SizeList& regs, SizeList& bulges ) const;
	std::size_t hash_value() const;
	bool contains(core::Size res) const;
	bool is_bulge(core::Size res) const;
	bool is_ladder() const;
	core::Size get_pair(core::Size res)const ;
	bool check_pleat() const;
	core::Size get_pleating(core::Size res)const ;
	bool extend(core::Size res, core::Size res2, bool antiparallel, core::Size pleating);
	void extend_to(core::Size res);
	bool antiparallel() const;
	bool has_pairing( core::scoring::dssp::Pairing const& ) const;
	bool has_common_pairing( StrandPairing const& other ) const;
	bool paired( core::Size res1, core::Size res2, bool anti ) const;
	void get_beta_pairs( core::scoring::dssp::PairingList& ) const;
	bool merge( StrandPairing const& other, bool domerge = false);
	bool mergeable( StrandPairing const& other ) const;
	void show_internals( std::ostream& out ) const;
	bool valid_ends() const;
	friend std::ostream & operator<<(std::ostream & out, StrandPairing const & sp);
	friend std::istream & operator>>( std::istream &in, StrandPairing &sp );
	static core::Size BIG_BULGE_LIMIT;
	static core::Size SMALL_BULGE_LIMIT;

private:
	//@brief first pairing of antipar. strand is begin1_ --> end2_
	//       last  pairing of antipar. strand is end1_ --> begin2_
	// hence strand1 goes from begin1-->end1 and strand2 goes from begin2-->end2
	core::Size begin1_, end1_, begin2_, end2_;

	// vector listing which residues are paired to the residues
	// in strand 1.  0 indicates beta bulge (unpaired).
	// one entry for pos begin1_...end1_
	std::vector<core::Size> pairing1;
	std::vector<core::Size> pleating1;
	// similar to pairing1 but for strand 2.
	std::vector<core::Size> pairing2;
	//bool pleat_weird;
	bool antipar;
};

std::ostream & operator<<(std::ostream & out, StrandPairing const& sp);
//////////////////////////////////////////////////////////////////
// StrandPairingSet
// This class maintains a set of strand pairings and provides
// access functions to determine whether a particular residue
// participates in any of them, and in what capacity.
//
//////////////////////////////////////////////////////////////////
class StrandPairingSet : public utility::pointer::ReferenceCount {
	typedef utility::vector1< StrandPairing > StrandPairings;
	typedef StrandPairings::iterator iterator;
public:
	typedef StrandPairings::const_iterator const_iterator;

	StrandPairingSet() {};
	StrandPairingSet( ObjexxFCL::FArray2_float const& hbonds, float threshold, core::pose::Pose const& );
	StrandPairingSet( core::pose::Pose const&, core::Real threshold = -0.5 );
	StrandPairingSet( core::scoring::dssp::PairingList const& );
	virtual ~StrandPairingSet();

	//void add_decoy( core::Size dec );
	bool check_pleat() const;
	char dssp_state( core::Size res ) const;
	char featurizer_state( core::Size res ) const;
	bool paired( core::Size res1, core::Size res2, bool antiparallel ) const;
	bool has_pairing( core::scoring::dssp::Pairing const& ) const;
	bool has_pairing( StrandPairing const& ) const;
	void get_beta_pairs( core::scoring::dssp::PairingList& ) const;
	bool merge( const StrandPairingSet &other, bool domerge = false );

	friend std::ostream & operator<<(std::ostream & out, const StrandPairingSet &sp);
	friend std::istream & operator>>(std::istream & is, StrandPairingSet &sp);

	const_iterator begin() const { return pairings_.begin(); }
	const_iterator end() const { return pairings_.end(); }
	Size size() const { return pairings_.size(); }

	StrandPairing const& strand_pairing( Size i ) const {
		runtime_assert( i <= pairings_.size() );
		return pairings_[ i ];
	}

	void push_back( StrandPairing const& sp ) {
		pairings_.push_back( sp );
	}
private:
	void add_pairing( core::Size res1, core::Size res2, bool antiparallel, core::Size pleating );
	void add_pairing( core::scoring::dssp::Pairing const& pairing );
	void selfmerge();
	void compute( ObjexxFCL::FArray2_float const& hbonds, float threshold, core::pose::Pose const& );

	StrandPairings pairings_;
};

std::ostream & operator<<(std::ostream & out, const StrandPairingSet &sp);

}
}
}

#endif
