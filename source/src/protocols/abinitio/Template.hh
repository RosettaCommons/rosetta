// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Oliver Lange

#ifndef INCLUDED_protocols_abinitio_Template_hh
#define INCLUDED_protocols_abinitio_Template_hh

// Unit Headers
#include <protocols/abinitio/Template.fwd.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/FrameList.fwd.hh>
#include <core/fragment/SingleResidueFragData.fwd.hh>

#include <core/sequence/DerivedSequenceMapping.hh>

#ifdef  WIN32
#include <core/scoring/constraints/NamedAtomPairConstraint.hh>
#endif


// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>


//// C++ headers
#include <string>

#include <core/scoring/constraints/AtomPairConstraint.fwd.hh>
#include <core/scoring/constraints/NamedAtomPairConstraint.fwd.hh>
#include <core/scoring/dssp/PairingsList.fwd.hh>
#include <core/scoring/dssp/StrandPairing.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace abinitio {

class Template : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~Template();
	typedef utility::vector1< core::scoring::constraints::AtomPairConstraintOP > AtomPairConstraintList;
	typedef utility::vector1< core::scoring::constraints::Obsolet_NamedAtomPairConstraintOP > NamedAtomPairConstraintList;

	static void register_options();

	Template( std::string const& name, core::pose::PoseCOP, core::sequence::DerivedSequenceMapping const& mapping);
	Template( std::string const& name, core::pose::PoseCOP, std::string const& map_file, int offset, core::Real score );

	//@brief pick fragments of <length> from Template according to mapping, return nr of Frames
	/// good for continuous fragments: no torsions from insertions are mapped
	//  Size pick_frags( core::fragment::FragSet&, core::fragment::FragDataCOP frag_type ) const;


	//@brief steals frames as templated in target_frames and accumulates them in FragSet
	// steal frags does not take care of insertion/deletions ( fine for e.g., jump-frames )
	Size steal_frags( core::fragment::FrameList const& target_frames, core::fragment::FragSet& accumulator, Size ncopies = 1 ) const;

	//@brief maps pairings from target to template ( or the reverse )
	void map_pairings2template( core::scoring::dssp::PairingList const& in, core::scoring::dssp::PairingList& out ) const;
	void map_pairings2target( core::scoring::dssp::PairingList const& in, core::scoring::dssp::PairingList& out ) const;

	//@brief generate new list of frames aligned to the target --> only alignable frames in target_frames
	void map2target( core::fragment::FrameList const& template_frames, core::fragment::FrameList& target_frames ) const;

	//@brief in-place mapping ---> breaks if not alignable
	void map2target( core::fragment::FrameList& ) const;

	//@brief generate new list of frame aligned to the template --> only alignable frames in template_frames
	void map2template( core::fragment::FrameList const& target_frames, core::fragment::FrameList& template_frames ) const;

	//@brief in-place mapping ----> breaks if not alignable
	void map2template( core::fragment::FrameList& ) const;

	Size pick_large_frags( core::fragment::FragSet& frag_set, core::fragment::SingleResidueFragDataOP frag_type, core::Size ncopies = 1 ) const;

	//@brief map constraints:
	// in the moment only support for AtomPairConstraints... otherwise would have to add clone() method to Constraint
	//  and write-accessors for AtomID to all derived classes
	void map2target( NamedAtomPairConstraintList const&, NamedAtomPairConstraintList& ) const;

	void map2template( NamedAtomPairConstraintList const&, NamedAtomPairConstraintList& ) const;

	//  core::scoring::constraints::NamedAtomPairConstraintOP
	//  map2template( core::scoring::constraints::NamedAtomPairConstraint const& ) const;

	void cull_violators( NamedAtomPairConstraintList const& target_list, NamedAtomPairConstraintList& culled_list ) const;

	void read_constraints( std::string const& cst_file );

	//@brief constraints for this tempalte are present
	bool has_constraints() const {
		return cstfile_.size();
		//		return cstset_.size();
	}

	core::scoring::dssp::StrandPairingSet const& strand_pairings() const {
		return *strand_pairings_;
	}

	//@brief constraints -- template numbering
	NamedAtomPairConstraintList const& constraints() const {
		if ( cstset_.size() == 0 ) _read_constraints( cstfile_ );
		return cstset_;
	}

  //void map2template( FrameList const& target_frames, FrameList const&  template_frame ) const;

  std::string const& name() const {
    return name_;
  }

	core::Real external_score() const {
		return score_;
	}

	std::string query_sequence() const {
		return mapping_.seq1();
	}

	core::Real topology_score() const {
		return topol_score_;
	}

	void topology_score( core::Real setting ) {
		topol_score_ = setting;
	}

  bool is_good() const { return good_; }

private:
	bool map_pairing( core::scoring::dssp::Pairing const&, core::scoring::dssp::Pairing&, core::sequence::DerivedSequenceMapping const& map ) const;
	void _read_constraints( std::string const& cst_file ) const;

	core::sequence::DerivedSequenceMapping mapping_; //target2template
	core::sequence::DerivedSequenceMapping reverse_mapping_; //template2target
	//  std::string aligned_seq_;
	//  std::string target_seq_;
	core::pose::PoseCOP pose_;
	std::string name_;
	core::Real score_; //the svr_alignment score
	core::Real topol_score_; //the pairing-stat score
	//constraints template numbering
	mutable NamedAtomPairConstraintList cstset_; //because of lazy read

	core::scoring::dssp::StrandPairingSetOP strand_pairings_;

	std::string cstfile_;

	bool good_;
};

} //abinitio
} //protocols

#endif
