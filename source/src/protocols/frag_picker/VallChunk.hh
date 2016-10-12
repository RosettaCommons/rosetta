// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/VallChunk.hh
/// @brief  a contiguous chunk of residues taken from a vall.
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


#ifndef INCLUDED_protocols_frag_picker_VallChunk_hh
#define INCLUDED_protocols_frag_picker_VallChunk_hh

// unit headers
#include <protocols/frag_picker/VallChunk.fwd.hh>

// package headers
#include <protocols/frag_picker/VallProvider.fwd.hh>
#include <protocols/frag_picker/VallResidue.hh>

// mini
#include <core/sequence/SequenceProfile.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

namespace protocols {
namespace frag_picker {

/// @brief  represents a chunk of residues extracted from a vall.
/// @details VallChunk contains a vector of VallResidue objects and provides a basic ways to access them
class VallChunk: public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< VallChunk >
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~VallChunk() override;

	VallChunk(VallProviderAP provider);

	inline VallChunkCOP get_self_ptr() const { return shared_from_this(); }
	inline VallChunkOP  get_self_ptr() { return shared_from_this(); }
	inline VallChunkCAP get_self_weak_ptr() const { return VallChunkCAP( shared_from_this() ); }
	inline VallChunkAP  get_self_weak_ptr() { return VallChunkAP( shared_from_this() ); }

	/// @brief  returns a PDB id (a string of four letters, e.g. "4mba")
	inline std::string get_pdb_id() const {
		return at(1)->id().substr(0, 4);
	}

	/// @brief  returns protein chain ID
	inline char get_chain_id() const {
		return at(1)->id()[4];
	}

	/// @brief  returns integer key of this chunk, which is the key of this chunk's first residue
	inline core::Size key() const {
		return at(1)->key();
	}

	/// @brief  returns the size of this chunk i.e. the number of residues stored in there
	inline core::Size size() const {
		return residues_.size();
	}

	inline core::Size vall_key() const {
		return vall_key_;
	}

	void vall_key(core::Size key) {
		vall_key_ = key;
	}

	/// @brief  returns i-th residue form this chunk. The first residue has index 1
	inline VallResidueOP at(core::Size index) const {
		runtime_assert( index <= residues_.size() );
		runtime_assert( index >= 1 );
		return residues_.at(index);
	}

	/// @brief  appends a residue to this chunk
	inline void push_back(VallResidueOP what) {
		residues_.push_back(what);
	}

	/// @brief  returns amino acid sequence of this chunk
	std::string& get_sequence();

	/// @brief  returns amino acid profile of this chunk
	/// @details the profile object is created when this function is called for the first time
	/// and then cached within a VallProvider object.
	/// Every time this method is called for a new chunk, VallProvider caches new data
	core::sequence::SequenceProfileOP get_profile();

	/// @brief  returns a pose created for this chunk
	/// @details the pose object is created when this function is called for the first time
	/// and then cached within a VallProvider object
	/// Every time this method is called for a new chunk, VallProvider caches new data
	core::pose::PoseOP get_pose();

	/// @brief returns the Vall file used for this chunk
	std::string get_vall_filename();
	/// @brief returns the Vall Provider
	VallProviderOP get_vall_provider();

	/// @brief  returns a string that is unique for each chunk in vall
	std::string & chunk_key() { if ( !has_key_ ) { create_key(); } return chunk_key_; }

private:
	utility::vector1<VallResidueOP> residues_;
	std::string sequence_;
	VallProviderAP my_provider_;
	std::string chunk_key_;
	core::Size vall_key_;
	bool has_key_;

	void create_key();
};

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_VallChunk_HH */
