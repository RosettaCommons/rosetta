// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragSet.cc
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @author James Thompson (tex@u.washington.edu)
/// @author Christopher Miles (cmiles@uw.edu)

// Unit Headers
#include <core/fragment/FragmentIO.hh>
#ifdef WIN32
#include <core/fragment/FragID.hh>
#endif

// Factory Headers
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/SecstructSRFD.hh>
#include <core/fragment/JumpSRFD.hh>
#include <core/fragment/JumpingFrame.hh>

// Packet Headers
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/OrderedFragSet.hh>

// Project Headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>
#include <ostream>

// Auto Headers
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FrameIterator.hh>

#include <utility/vector1.hh>


namespace core {
namespace fragment {

static thread_local basic::Tracer tr( "core.io.fragments" );

FragFactory FragmentIO::frag_factory_;
FragmentIO::FragFileCache FragmentIO::frag_cache_;

void
FragFactory::add_frag_type( std::string const&  type_name, SingleResidueFragDataOP new_frag ) {
	frag_types_[ type_name ] = new_frag;
}

void
FragFactory::add_frame_type( std::string const& type_name, FrameOP new_frame ) {
	frame_types_[ type_name ] = new_frame;
}

FragFactory::FragFactory(void) {
	// initialization of fragtions which this factory knows how to instantiate
	FragFactory::add_frag_type( BBTorsionSRFD::_static_type_name() , new BBTorsionSRFD );
	FragFactory::add_frag_type( DownJumpSRFD::_static_type_name(), new DownJumpSRFD );
	FragFactory::add_frag_type( UpJumpSRFD::_static_type_name(), new UpJumpSRFD );
	FragFactory::add_frag_type( SecstructSRFD::_static_type_name(), new SecstructSRFD );

	FragFactory::add_frame_type( Frame::_static_type_name(), new Frame );
	FragFactory::add_frame_type( JumpingFrame::_static_type_name(), new JumpingFrame );
}

FrameOP FragFactory::frame( std::string const& type ) const {
	FrameTypes::const_iterator const iter ( frame_types_.find( type ) );
	if ( iter != frame_types_.end() ) {
		return iter->second->clone();
	}
	return NULL;
}

SingleResidueFragDataOP
FragFactory::frag_type( std::string const& frag_name ) const {
	SRFD_Types::const_iterator const iter ( frag_types_.find( frag_name ) );
	if ( iter != frag_types_.end() ) {
		return iter->second->clone();
	}
	return NULL;
}

FragFactory & FragmentIO::get_frag_factory(void) {
	return frag_factory_;
}

void FragmentIO::read_next_frames( std::istream& data, std::string& next_line, FrameList &next_frames ) {
	//read lines defining the FRAMES for a given block of Fragments
	do {
		if ( !next_line.size() ) continue; //skip empty lines
		std::string tag;
		std::istringstream line_stream( next_line );
		line_stream >> tag;
		FrameOP new_frame = frag_factory_.frame( tag );
		if ( new_frame ) {
			new_frame->read( line_stream );
			tr.Trace << "read frame " << *new_frame << std::endl;
			next_frames.push_back( new_frame );
		} else break; // not a Frame tag: line means something else .. finished reading frames for this block
		if ( line_stream.fail() ) data.setstate( std::ios_base::failbit );
	} while ( getline( data, next_line ) );
}

void FragmentIO::read_frag_data( std::istream& data, std::string& next_line, FrameList &new_frames ) {
	if ( !new_frames.size() ) return;
	Size const length( new_frames[ 1 ]->length() );
	tr.Trace << "frags in frame at beginning: " << new_frames[ 1 ]->nr_frags() << std::endl;

	FragDataOP current_fragment = NULL;//new FragData;
	do {
		if ( !next_line.size() ) continue; //skip empty lines
		if ( next_line == "ENDFRAME" ) { //can probably vanish soon
			getline( data, next_line );
			break;
		}
		if ( next_line == "FRAME" ) { // shouldn't work, since next_line will be FRAME     12 20
			break;
		}
		std::string tag;
		Size pos;
		Size pdb_pos;
		std::string pdb_id;

		std::istringstream line_stream( next_line + " " );
		line_stream >> pos >> pdb_pos >> pdb_id >> tag;
		if ( line_stream.fail() ) break; //e.g., we have a line that starts with a tag

		//	bIgnoreFragments if we have already top_ -- but don't forget to eat all lines until next "FRAME"
		if ( top_ && (*new_frames.begin())->nr_frags() >= top_ * ncopies_ ) continue;

		if ( !current_fragment ) {
			if ( pdb_pos && bAnnotate_ ) current_fragment = new AnnotatedFragData( pdb_id, pdb_pos );
			else current_fragment = new FragData;
		}

		SingleResidueFragDataOP new_srfd = frag_factory_.frag_type( tag );
		if ( !new_srfd ) {
			tr.Error << "ERROR: tag " << tag << " not known by frag_factory" << std::endl;
			break;
		}
		line_stream >> *new_srfd;
		if ( line_stream.fail() ) {
			tr.Error << "ERROR: reading of " << tag << " failed in line: " << next_line << std::endl;
			data.setstate( std::ios_base::failbit );
			break;
		}
		current_fragment->add_residue( new_srfd );
		if ( current_fragment->size() == length ) {
			current_fragment->set_valid();
			tr.Trace << "read FragData" << std::endl << *current_fragment << std::endl;
			for ( FrameList::iterator frame = new_frames.begin(), eframe = new_frames.end();
						frame != eframe; ++frame ) {

				//jump out if top_ fragments have been read
				if ( !(*frame)->add_fragment( current_fragment ) ) {
					tr.Error << "failed to add " << *current_fragment << " to frame " << (**frame) << std::endl;
					data.setstate( std::ios_base::failbit );
				} else {
					for ( Size i = 2; i<=ncopies_; i++ ) {
						(*frame)->add_fragment( current_fragment );
					}
				}
			} //frames
			current_fragment = NULL;
		}
	} while ( getline( data, next_line ) );
	if ( data.fail() && data.eof() && !data.bad() ) { //stupid getline sets fail bit instead of eof if last line is read.
		tr.Debug << "End of fragment file during fragment reading... set eofbit" << std::endl;
		data.clear( std::ios_base::eofbit );
	}
	if ( current_fragment ) {
		tr.Error << "ERROR: read incomplete fragment at line" << next_line << std::endl;
		data.setstate( std::ios_base::failbit );
	}
}

void FragmentIO::read_data( std::istream& data, FrameList& all_frames ) {
	tr.Info << " read fragments options: top = " << top_ << " copies: " << ncopies_
					<< " annotate: " << ( bAnnotate_ ? "yes" : "no" ) << std::endl;

	std::string next_line;
	while ( data.good() ) {
		FrameList new_frames;

		read_next_frames( data, next_line, new_frames );
		if ( new_frames.size() == 0 ) tr.Debug << "found no new FRAME entry ... finished. " << std::endl;
		if ( data.eof() ) tr.Debug << " reached EOF during reading of FRAME entries " << std::endl;
		if ( !data || new_frames.size() == 0 ) break;

		read_frag_data( data, next_line, new_frames );
		if ( data.fail() ) break;

		for ( FrameList::iterator frame = new_frames.begin(), eframe = new_frames.end();
					frame != eframe; ++frame ) {
			if ( (*frame)->nr_frags() ) {
				all_frames.push_back( *frame );
			}
		}
	}
	if ( data.fail() || data.bad() ) {
		tr.Error << "reading failed at line " << next_line << std::endl;
	}
}

FragSetOP FragmentIO::read_data( std::string const& filename ) {
	//detect format:
	std::string line;

	if ( frag_cache_.find( filename ) == frag_cache_.end() ) {
		//okay we are about to read a new fragment file...
		// at this moment we loose trust in fragments that are only referenced by the frag_cache..
		clean_frag_cache();

		// format check...
		utility::io::izstream data( filename );
		if ( !data.good() )
			utility_exit_with_message("ERROR: FragmentIO: could not open file " + filename );

		tr.Info << "reading fragments from file: " << filename << " ... " << std::endl;
		getline( data, line );
		std::istringstream str( line );
		std::string tag;
		str >> tag;
		if ( tag == "position:" ) {
			tr.Info << "rosetta++ fileformat detected! Calling legacy reader... "
							<< std::endl;

			ConstantLengthFragSetOP frags = new ConstantLengthFragSet;
			frags->read_fragment_file( filename, top_, ncopies_, bAnnotate_ );
			frag_cache_[ filename ] = frags;
			return frags;
		}

		FrameList frames;
		read_data( filename, frames );

		//find out if ConstantLengthFragSet is sufficient...
		//  for now always use OrderedFragSet
		FragSetOP frags = new OrderedFragSet();
		frags->add( frames );
		frag_cache_[ filename ] = frags;
	}
	return frag_cache_[ filename ];
}

void FragmentIO::clean_frag_cache() {
	/// FIX ME: this will not work with std::shared_ptr as the reference
	/// count is not stored with the object, but with the owning_ptr.
	/// Commenting out for now.
	/*
	for ( FragFileCache::iterator it = frag_cache_.begin(); it != frag_cache_.end();  ) {
		if ( it->second->ref_count() == 1 ) {
			tr.Info << "GARBAGE COLLECTION: remove " << it->first << " from frag_cache " << std::endl;
			frag_cache_.erase( it++ ); //doesn't return iterator ... help with post-fix iterator
		} else ++it;
	}
	*/
}

void FragmentIO::read_data( std::string const& filename, FrameList& frames ) {
	utility::io::izstream data( filename );
	if ( !data.good() ) utility_exit_with_message("ERROR: FragmentIO: could not open file " + filename );
	tr.Info << "reading fragment frames from file: " << filename << " ... " << std::endl;

	read_data( data, frames);
	if ( data.fail() ) {
		utility_exit_with_message( "failed to read fragments from file " + filename );
	}
}

void FragmentIO::write_data( std::string const& file, FragSet const& frags )  {
	utility::io::ozstream out( file );
	for ( ConstFrameIterator it=frags.begin(), eit=frags.end(); it!=eit; ++it ) {
		(*it)->show( out );
	}
}

} //fragment
} //core
