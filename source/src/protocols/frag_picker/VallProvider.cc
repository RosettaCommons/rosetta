// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/VallProvider.cc
/// @brief  reads a vall library and serves the data in chunks
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


// type headers
#include <core/types.hh>

#include <protocols/frag_picker/VallProvider.hh>

// package headers
#include <protocols/frag_picker/VallResidue.hh>
#include <protocols/frag_picker/VallChunk.hh>

// project headers
#include <basic/Tracer.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ChemicalManager.hh>


// utility headers
#include <utility/exit.hh>
#include <utility/io/izstream.hh>

// C++ headers
#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif


#include <sstream>
#include <string>

//Auto Headers
#include <core/pose/annotated_sequence.hh>


// Tracer instance for this file
// Named after the original location of this code
static thread_local basic::Tracer TR( "protocols.frag_picker.VallProvider" );

namespace protocols {
namespace frag_picker {

using namespace core;

VallProvider::~VallProvider() {}

VallChunkOP VallProvider::find_chunk(std::string pdb_id, char chain_id, Size residue_id) {

	for (Size i = 1; i <= chunks_.size(); ++i) {
		VallChunkOP c = chunks_[i];
		if (c->get_pdb_id().compare(pdb_id) != 0)
			continue;
		if (c->get_chain_id() != chain_id)
			continue;
		for (Size j = 1; j < c->size(); ++j) {
			if (c->at(j)->resi() == residue_id)
				return c;
		}
	}

	TR.Warning << "Can't find chunk that contains the residue number "
			<< residue_id << " in a protein " << pdb_id << " in a chain "
			<< chain_id << std::endl;

	return 0;
}

core::pose::PoseOP VallProvider::cache_pose(VallChunkOP source_chunk) {

	std::string key = source_chunk->chunk_key();
	TR.Debug << "looking for a pose for " << key << std::endl;
	if (key.compare(cached_pose_id_) == 0)
		return cached_pose_;

	TR.Debug << "caching a pose for >" << key << "< the previous key was: >"
			<< cached_pose_id_ << "<  (" << key.compare(cached_pose_id_)
			<< ")... ";

	for (Size i = 1; i <= source_chunk->size(); ++i) {
		VallResidueOP r = source_chunk->at(i);
		cached_pose_->set_phi(i, r->phi());
		cached_pose_->set_psi(i, r->psi());
		cached_pose_->set_omega(i, r->omega());
	}
	TR.Debug << " has " << source_chunk->size() << " residues." << std::endl;

	cached_pose_id_.clear();
	cached_pose_id_.assign(key);

//	pose->dump_pdb("dump-"+chunk->get_pdb_id()+".pdb");

	return cached_pose_;
}

core::sequence::SequenceProfileOP VallProvider::cache_profile(VallChunkOP source_chunk) {

	std::string key = source_chunk->chunk_key();
	TR.Debug << "looking for a profile for " << key << std::endl;
	if (key.compare(cached_profile_id_) == 0)
		return cached_profile_;
	TR.Debug << "caching a sequence profile for >" << key
			<< "< the previous key was: >" << cached_profile_id_ << "< ("
			<< key.compare(cached_pose_id_) << ")... ";

	utility::vector1<utility::vector1<core::Real> > prof;
	for (Size i = 1; i <= source_chunk->size(); ++i) {
		prof.push_back(source_chunk->at(i)->profile());
	}
	cached_profile_id_.clear();
	cached_profile_id_.assign(key);

	if (cached_profile_ == 0)
		cached_profile_ = new core::sequence::SequenceProfile();
	cached_profile_->profile(prof);
	cached_profile_->sequence(source_chunk->get_sequence());

	TR.Debug << " has " << cached_profile_->length() << " columns."
			<< std::endl;

	return cached_profile_;
}

Size VallProvider::vallChunksFromLibraries(utility::vector1< std::string > const & fns) {

	using std::string;
	using utility::vector1;

	typedef vector1< string >::const_iterator iter;
	for ( iter it = fns.begin(), end = fns.end(); it != end; ++it ) {
		vallChunksFromLibrary( *it );
	}

	return 0;
}
Size VallProvider::vallNumLines(std::string const & filename) {
	Size num_lines;
		utility::io::izstream stream(filename);
	if (!stream)
		utility_exit_with_message( "can't open file: " + filename );
	num_lines = 0;
	std::string line;
	getline(stream, line);
	while(line[0] == '#') {
		getline(stream, line);
	}
	while(getline(stream,line)) {
		num_lines++;
	}
	return num_lines;
}

Size VallProvider::vallChunksFromLibrary(std::string const & filename, core::Size startline, core::Size endline) {

	TR.Info << "vallChunksFromLibrary" << std::endl;

	utility::io::izstream stream(filename);
	if (!stream)
		utility_exit_with_message( "can't open file: " + filename );

	// statistics
	Size n_lines = 0;

	TR.Info << "Reading Vall library from " << filename << " ... startline: " << startline << "  endline: " << endline << std::endl;

	time_t time_start = time(NULL);

	Size last_key = 0;
	if(chunks_.size() > 0) {
	    VallChunkOP last_chunk = chunks_[ chunks_.size() ];
	    last_key = last_chunk->at( last_chunk->size() )->key();
	}
	std::string prior_id = "";
	Size prior_resi = 0;
	VallChunkOP current_section = new VallChunk(get_self_weak_ptr());
	std::string line;
	getline(stream, line);
	while(line[0] == '#') getline(stream, line);

	VallResidueOP firstRes = new VallResidue();
	firstRes->key(1);

	//Decides if vall is in the old nnmake format
	//if not it tries to use the nnmake + chemical shift's format
	if (line.length() > 300)
		firstRes->fill_from_string_version1(line);
	else if (line.length() > 240)
		firstRes->fill_from_string_cs(line);
	else if (line.length() > 110)
		firstRes->fill_from_string(line);
	else
		firstRes->fill_from_string_residue_depth_version1(line);

	current_section->push_back(firstRes);
	prior_id = firstRes->id();
	prior_resi = firstRes->resi();
	// parse Vall from file
	n_lines = 1;

	while (getline(stream, line)) {
		++n_lines;
		if( n_lines < startline ) continue;

		// If endline is 0, just read to the end of the file
		if( endline != 0 && n_lines > endline ) break;

		VallResidueOP current_residue = new VallResidue();

		if (line.length() > 300)
			current_residue->fill_from_string_version1(line);
		else if (line.length() > 240)
			current_residue->fill_from_string_cs(line);
		else if (line.length() > 110)
			current_residue->fill_from_string(line);
		else
			current_residue->fill_from_string_residue_depth_version1(line);

		current_residue->key( n_lines + last_key );
		// check for start of new continuous stretch
		if ( (current_residue->resi() != prior_resi + 1)
				|| (current_residue->id() != prior_id)) {
			push_back(current_section);
			Size t = current_section->size();
			if (t > largest_chunk_size_)
			    largest_chunk_size_ = t;
			TR.Debug << "Created a new chunk for : " << current_section->get_pdb_id()
					<< " having " << current_section->size() << " residues "
					<< " at index " << size() << ". The largest chunk's size is: "
					<<largest_chunk_size_<<std::endl;
			current_section = new VallChunk(get_self_weak_ptr());
			prior_id = current_residue->id();
		}
		prior_resi = current_residue->resi();
		current_section->push_back(current_residue);
//		if (n_lines % 100000 == 0) {
//			TR.Info << "   " << n_lines << std::endl;
//			TR.flush();
//		}

	} // line loop

	push_back(current_section);
	TR.Debug << "Created a new chunk for : " << current_section->get_pdb_id()
			<< " having " << current_section->size() << " residues "
			<< " at index " << size() << std::endl;
	Size t = current_section->size();
	if (t > largest_chunk_size_) largest_chunk_size_ = t;

	time_t time_end = time(NULL);

	TR.Info << "... done.  Read " << n_lines << " lines.  Time elapsed: "
			<< (time_end - time_start) << " seconds." << std::endl;

	TR.Info << "Total chunks: " << size() << std::endl;
	TR.Info << "Largest chunk: " << largest_chunk_size_ << " aa" << std::endl;

	// create cached pose
	for (Size i = 1; i <= largest_chunk_size_; i++) poly_A_seq_ += "A";
	cached_pose_ = new core::pose::Pose();
	core::pose::make_pose_from_sequence(*cached_pose_, poly_A_seq_,
			*(chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard")));
	TR.flush();
	return n_lines;
}

} // frag_picker
} // protocols
