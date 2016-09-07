// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/VallProvider.hh
/// @brief  reads a vall library and serves the data in chunks
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_VallProvider_hh
#define INCLUDED_protocols_frag_picker_VallProvider_hh

// type headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>

// unit headers
#include <protocols/frag_picker/VallProvider.fwd.hh>

// package headers
#include <protocols/frag_picker/VallChunk.hh>

// mini headers
#include <core/sequence/SequenceProfile.hh>
#include <core/pose/Pose.hh>


namespace protocols {
namespace frag_picker {

/// @brief  a vector of vall chunks
class VallProvider: public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< VallProvider >
{
public:

	/// @brief
	VallProvider() {
		cached_profile_id_ = "";
		cached_pose_id_ = "";
		largest_chunk_size_ = 0;
	}
	~VallProvider() override;

	inline VallProviderCOP get_self_ptr() const { return shared_from_this(); }
	inline VallProviderOP  get_self_ptr() { return shared_from_this(); }
	inline VallProviderCAP get_self_weak_ptr() const { return VallProviderCAP( shared_from_this() ); }
	inline VallProviderAP  get_self_weak_ptr() { return VallProviderAP( shared_from_this() ); }

	/// @brief Vall reader
	/// THe defaults should ensure that the file is fully read if startline and endline ar not specified. endline = 0 means read to the end.
	Size vallChunksFromLibrary(std::string const & filename, Size startline = 1, Size endline = 0 );

	Size vallChunksFromLibraries( utility::vector1< std::string > const & fns );

	/// @brief Runs through the Vall and stores number of lines
	Size vallNumLines(std::string const & filename);

	/// @brief says how many chunks do we have
	inline Size size() {
		return chunks_.size();
	}

	/// @brief returns a certain chunk (starts from 1)
	inline VallChunkOP at(Size index) {
		return chunks_.at(index);
	}

	/// @brief adds a new chunk to this provider
	inline void push_back(VallChunkOP what) {
		return chunks_.push_back(what);
	}

	/// @brief says what is the length of the largest chunk known to this
	/// provider
	inline Size get_largest_chunk_size() {
		return largest_chunk_size_;
	}

	inline Size get_vall_count() {
		return vall_keys_.size();
	}

	inline std::string get_vall_by_key(Size index) {
		return vall_keys_.at(index);
	}

	inline Size get_vall_start_line_by_key(Size index) {
		return vall_start_line_.at(index);
	}

	inline Size get_vall_end_line_by_key(Size index) {
		return vall_end_line_.at(index);
	}

	inline Size get_vall_last_residue_key_by_key(Size index) {
		return vall_last_residue_key_.at(index);
	}

	/// @brief tries to find a chunk defined by PDB id, chain id and a residue
	/// sequence id @details If this VallProvider does not contain a desired
	/// chunk, 0 is returned.
	VallChunkOP find_chunk(std::string, char, Size);

	/// @brief cache a sequence profile for a given chunk
	core::sequence::SequenceProfileOP cache_profile(VallChunkOP source_chunk);

	/// @brief cache a pose for a given chunk
	core::pose::PoseOP cache_pose(VallChunkOP source_chunk);

private:
	utility::vector1<VallChunkOP> chunks_;
	utility::vector1<std::string> vall_keys_;
	utility::vector1<Size> vall_start_line_;
	utility::vector1<Size> vall_end_line_;
	utility::vector1<Size> vall_last_residue_key_;
	Size largest_chunk_size_;
	std::string cached_profile_id_;
	core::sequence::SequenceProfileOP cached_profile_;
	core::pose::PoseOP cached_pose_;
	std::string cached_pose_id_;
	std::string poly_A_seq_;
};

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_VallProvider_HH */
