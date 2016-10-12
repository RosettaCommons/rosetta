// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
static THREAD_LOCAL basic::Tracer TR( "protocols.frag_picker.VallProvider" );

namespace protocols {
namespace frag_picker {

using namespace core;

VallProvider::~VallProvider() = default;

VallChunkOP VallProvider::find_chunk(std::string pdb_id, char chain_id, core::Size residue_id) {

	for ( core::Size i = 1; i <= chunks_.size(); ++i ) {
		VallChunkOP c = chunks_[i];
		if ( c->get_pdb_id().compare(pdb_id) != 0 ) {
			continue;
		}
		if ( c->get_chain_id() != chain_id ) {
			continue;
		}
		for ( core::Size j = 1; j < c->size(); ++j ) {
			if ( c->at(j)->resi() == residue_id ) {
				return c;
			}
		}
	}

	TR.Warning << "Can't find chunk that contains the residue number "
		<< residue_id << " in a protein " << pdb_id << " in a chain "
		<< chain_id << std::endl;

	return nullptr;
}

core::pose::PoseOP VallProvider::cache_pose(VallChunkOP source_chunk) {

	std::string key = source_chunk->chunk_key();
	TR.Debug << "looking for a pose for " << key << std::endl;
	if ( key.compare(cached_pose_id_) == 0 ) {
		return cached_pose_;
	}

	TR.Debug << "caching a pose for >" << key << "< the previous key was: >"
		<< cached_pose_id_ << "<  (" << key.compare(cached_pose_id_)
		<< ")... ";

	for ( core::Size i = 1; i <= source_chunk->size(); ++i ) {
		VallResidueOP r = source_chunk->at(i);
		cached_pose_->set_phi(i, r->phi());
		cached_pose_->set_psi(i, r->psi());
		cached_pose_->set_omega(i, r->omega());
	}
	TR.Debug << " has " << source_chunk->size() << " residues." << std::endl;

	cached_pose_id_.clear();
	cached_pose_id_.assign(key);

	// pose->dump_pdb("dump-"+chunk->get_pdb_id()+".pdb");

	return cached_pose_;
}

core::sequence::SequenceProfileOP VallProvider::cache_profile(VallChunkOP source_chunk) {

	std::string key = source_chunk->chunk_key();
	TR.Debug << "looking for a profile for " << key << std::endl;
	if ( key.compare(cached_profile_id_) == 0 ) {
		return cached_profile_;
	}
	TR.Debug << "caching a sequence profile for >" << key
		<< "< the previous key was: >" << cached_profile_id_ << "< ("
		<< key.compare(cached_pose_id_) << ")... ";

	utility::vector1<utility::vector1<core::Real> > prof;
	for ( core::Size i = 1; i <= source_chunk->size(); ++i ) {
		prof.push_back(source_chunk->at(i)->profile());
	}
	cached_profile_id_.clear();
	cached_profile_id_.assign(key);

	if ( cached_profile_ == nullptr ) {
		cached_profile_ = core::sequence::SequenceProfileOP( new core::sequence::SequenceProfile() );
	}
	cached_profile_->profile(prof);
	cached_profile_->sequence(source_chunk->get_sequence());

	TR.Debug << " has " << cached_profile_->length() << " columns."
		<< std::endl;

	return cached_profile_;
}

core::Size VallProvider::vallChunksFromLibraries(utility::vector1< std::string > const & fns) {

	using std::string;
	using utility::vector1;

	for ( auto const & fn : fns ) {
		vallChunksFromLibrary( fn );
	}

	return 0;
}
core::Size VallProvider::vallNumLines(std::string const & filename) {
	core::Size num_lines;
	utility::io::izstream stream(filename);
	if ( !stream ) {
		utility_exit_with_message( "can't open file: " + filename );
	}
	num_lines = 0;
	std::string line;
	getline(stream, line);
	while ( line[0] == '#' ) {
		getline(stream, line);
	}
	while ( getline(stream,line) ) {
		num_lines++;
	}
	return num_lines;
}

core::Size VallProvider::vallChunksFromLibrary(std::string const & filename, core::Size startline, core::Size endline) {

	TR.Info << "vallChunksFromLibrary" << std::endl;

	utility::io::izstream stream(filename);
	if ( !stream ) {
		utility_exit_with_message( "can't open file: " + filename );
	}

	// save Vall filename
	vall_keys_.push_back(filename);

	// statistics
	core::Size n_lines = 0;

	TR.Info << "Reading Vall library from " << filename << " ... startline: " << startline << "  endline: " << endline << std::endl;

	time_t time_start = time(nullptr);

	// get last residue key
	core::Size last_key = 0;
	if ( chunks_.size() > 0 ) {
		VallChunkOP last_chunk = chunks_[ chunks_.size() ];
		last_key = last_chunk->at( last_chunk->size() )->key();
	}
	vall_last_residue_key_.push_back(last_key);

	std::string prior_id = "";
	core::Size prior_resi = 0;
	VallChunkOP current_section( new VallChunk(get_self_weak_ptr()) );
	current_section->vall_key(vall_keys_.size());


	std::string line;
	while ( getline(stream, line) ) {
		if ( line[0] == '#' ) continue;
		++n_lines;
		if ( n_lines < startline ) continue;

		// If endline is 0, just read to the end of the file
		if ( endline != 0 && n_lines > endline ) break;

		VallResidueOP firstRes( new VallResidue() );

		// Chunk residue key is the last vall residue key + the vall line number for the residue
		firstRes->key(n_lines + last_key);

		//Decides if vall is in the old nnmake format
		//if not it tries to use the nnmake + chemical shift's format
		if ( line.length() > 300 ) {
			firstRes->fill_from_string_version1(line);
		} else if ( line.length() > 240 ) {
			firstRes->fill_from_string_cs(line);
		} else if ( line.length() > 110 ) {
			firstRes->fill_from_string(line);
		} else {
			firstRes->fill_from_string_residue_depth_version1(line);
		}

		current_section->push_back(firstRes);
		prior_id = firstRes->id();
		prior_resi = firstRes->resi();
		vall_start_line_.push_back(n_lines);
		break;
	}

	// parse Vall from file
	core::Size end_line = 0;
	while ( getline(stream, line) ) {
		if ( line[0] == '#' ) continue;
		++n_lines;
		if ( n_lines < startline ) continue;

		// If endline is 0, just read to the end of the file
		if ( endline != 0 && n_lines > endline ) break;

		VallResidueOP current_residue( new VallResidue() );

		if ( line.length() > 300 ) {
			current_residue->fill_from_string_version1(line);
		} else if ( line.length() > 240 ) {
			current_residue->fill_from_string_cs(line);
		} else if ( line.length() > 110 ) {
			current_residue->fill_from_string(line);
		} else {
			current_residue->fill_from_string_residue_depth_version1(line);
		}

		// Chunk residue key is the last vall residue key + the vall line number for the residue
		current_residue->key( n_lines + last_key );

		// check for start of new continuous stretch
		if ( (current_residue->resi() != prior_resi + 1)
				|| (current_residue->id() != prior_id) ) {
			push_back(current_section);
			core::Size t = current_section->size();
			if ( t > largest_chunk_size_ ) {
				largest_chunk_size_ = t;
			}
			TR.Debug << "Created a new chunk for : " << current_section->get_pdb_id()
				<< " having " << current_section->size() << " residues "
				<< " at index " << size() << ". The largest chunk's size is: "
				<<largest_chunk_size_<<std::endl;
			current_section = VallChunkOP( new VallChunk(get_self_weak_ptr()) );
			current_section->vall_key(vall_keys_.size());
			prior_id = current_residue->id();
		}
		prior_resi = current_residue->resi();
		current_section->push_back(current_residue);
		end_line = n_lines;
		//  if (n_lines % 100000 == 0) {
		//   TR.Info << "   " << n_lines << std::endl;
		//   TR.flush();
		//  }

	} // line loop


	vall_end_line_.push_back(end_line);

	push_back(current_section);
	TR.Debug << "Created a new chunk for : " << current_section->get_pdb_id()
		<< " having " << current_section->size() << " residues "
		<< " at index " << size() << std::endl;
	core::Size t = current_section->size();
	if ( t > largest_chunk_size_ ) largest_chunk_size_ = t;

	time_t time_end = time(nullptr);

	if ( vall_keys_.size() > 0 && vall_keys_.size() == vall_end_line_.size() && vall_keys_.size() == vall_start_line_.size() && vall_keys_.size() == vall_last_residue_key_.size() ) {
		TR.Debug << "Vall key: " << vall_keys_.size() << std::endl;
	} else {
		utility_exit_with_message( "There was an error reading the Vall: " + filename );
	}

	TR.Info << "... done.  Read " << n_lines << " lines.  Time elapsed: "
		<< (time_end - time_start) << " seconds." << std::endl;

	TR.Info << "Total chunks: " << size() << std::endl;
	TR.Info << "Largest chunk: " << largest_chunk_size_ << " aa" << std::endl;

	// create cached pose
	for ( core::Size i = 1; i <= largest_chunk_size_; i++ ) poly_A_seq_ += "A";
	cached_pose_ = core::pose::PoseOP( new core::pose::Pose() );
	core::pose::make_pose_from_sequence(*cached_pose_, poly_A_seq_,
		*(chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard")));
	TR.flush();
	return n_lines;
}

} // frag_picker
} // protocols
