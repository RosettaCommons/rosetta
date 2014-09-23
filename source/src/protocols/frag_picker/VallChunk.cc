// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/VallChunk.cc
/// @brief  a contiguous chunk of residues taken from a vall.
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


// unit headers
#include <protocols/frag_picker/VallChunk.hh>

// package headers
#include <protocols/frag_picker/VallResidue.hh>
#include <protocols/frag_picker/VallProvider.hh>

// mini
#include <core/sequence/SequenceProfile.hh>
#include <core/pose/Pose.hh>

// C++ stuff
#include <string>
#include <sstream>


namespace protocols {
namespace frag_picker {

/// @details Auto-generated virtual destructor
VallChunk::~VallChunk() {}

VallChunk::VallChunk(VallProviderAP provider) {
	sequence_ = "";
	my_provider_ = provider;
	has_key_ = false;
}


void VallChunk::create_key() {

	std::stringstream out;
	out << get_pdb_id() << get_chain_id() << ':' << at(1)->resi();
	chunk_key_ =  out.str();
	has_key_ = true;
}

std::string& VallChunk::get_sequence() {
	if (sequence_.length() == 0) {
		for (Size i = 1; i <= residues_.size(); i++) {
			char next = residues_.at(i)->aa();
			sequence_ += next;
		}
	}
	return sequence_;
}

core::sequence::SequenceProfileOP VallChunk::get_profile() {

	VallProviderOP my_provider( my_provider_ );
	return my_provider->cache_profile( get_self_ptr() );
}

core::pose::PoseOP VallChunk::get_pose() {

	VallProviderOP my_provider( my_provider_ );
	return my_provider->cache_pose( get_self_ptr() );
}

} // frag_picker
} // protocols
