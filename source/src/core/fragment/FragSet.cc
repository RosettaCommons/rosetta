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
/// @date   Wed Oct 20 12:08:31 2007
/// @author Roland A. Pache, PhD
///

// Unit Headers
#include <core/fragment/FragSet.hh>

// Package Headers
// AUTO-REMOVED #include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FragData.hh>


// Project Headers
#include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/types.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
#include <utility/pointer/owning_ptr.hh>

#include <basic/Tracer.hh>

#include <ostream>
// AUTO-REMOVED #include <set>

#include <core/fragment/FragID.hh>
#include <core/fragment/FrameIterator.hh>
#include <utility/vector1.hh>


namespace core {
namespace fragment {

using namespace kinematics;

static thread_local basic::Tracer tr( "core.fragment" );

///@brief return a list of frames that all sample the specified region, assume all motions are allowed
Size
FragSet::region_all(
	core::Size start,
	core::Size end,
	core::Size min_overlap,
	core::Size min_length,
	FrameList &frames
) const {
	kinematics::MoveMap move_map;
	move_map.set_bb( true );
	move_map.set_chi( true );
	move_map.set_jump( true );
	return region( move_map, start, end, min_overlap, min_length, frames );
}

///@brief returns the number and list of all fragment alignment frames that somehow overlap with the given region
///(also allows those frames that start before the region and reach into it)
Size FragSet::overlapping_with_region(
    kinematics::MoveMap const& mm,
    core::Size start,
    core::Size end,
    core::Size min_overlap,
    core::Size min_length,
    FrameList &frames
) const {
    //call the region method by default (method can be overloaded to change behavior in child classes of FragSet)
    return region( mm, start, end, min_overlap, min_length, frames );
}
    

	// put all fragments in FragID_list into this FragSet.
	// this function has the following effect:
	//      fragments that belong to the same frame are copied into a new frame
	//      the frame gets added. If all fragments of a frame are in the list, the frame is just added as is
	//
void FragSet::insert_fragID_list( FragID_List& list ) {
	utility::vector1< bool > handled(list.size(), false );
	Size pos( 1 );
	for ( FragID_List::iterator it=list.begin(), eit=list.end(); it!=eit; ++it,++pos ) {
		if ( !handled[ pos ] ) {
			FrameOP new_frame = it->frame().clone();
			Size spos( pos );
			for ( FragID_List::iterator sit=it; sit!=eit; ++sit, spos++ ) {
				if ( sit->frame_ptr() == it->frame_ptr() ) {
					handled[ spos ] = true;
					new_frame->add_fragment( &( sit->fragment() ) );
					new_frame->clone_cache_data( it->frame(), it->id(), new_frame->nr_frags() /* last added fragment */ );
				};
			}

			tr.Debug << pos << ": add frame " << new_frame << " " << new_frame->nr_frags() << std::endl;
			add( new_frame );
		} // handled
	}
}

void FragSet::add( FragID const& frag_id ) {
	Frame const& aFrame( frag_id.frame() );
	//std::cerr << "FragSet::add_frame " << std::endl;
	runtime_assert( aFrame.nr_frags() ); // do not allow insertion of empty frames --> makes frag_id iterator sooo much easier
	Size start ( aFrame.start() );
	Size end ( aFrame.end() );
	Size length( aFrame.length() );

	if ( min_pos() > start ) {
		set_min_pos( start );
	};

	if ( max_pos() < end ) {
		set_max_pos( end );
	};

	//	tr.Trace << "frag length " << length << " ( " << max_frag_length() <<  " ) " << std::endl;
	if ( length > max_frag_length() ) {
		tr.Trace << "set max frag length " << length << std::endl;
		set_max_frag_length( length );
	}

	// now add frame:
	FrameList present_frames;
	Size nr_present = frames( start, present_frames );
	if ( nr_present ) {
		for ( FrameList::iterator it = present_frames.begin(),
						eit = present_frames.end(); it!=eit; ++it ) {
			if ( (*it)->is_mergeable( aFrame ) ) {
				Size const new_id( (*it)->add_fragment( const_cast< FragData*>( frag_id.fragment_ptr().get() ) ) );
				(*it)->clone_cache_data( aFrame, frag_id.id(), new_id );
				return; //finished early
			}
		}
	}
	//didn't found mergable frames at this sequence position
	// make a new empty frame for this Fragment
	FrameOP new_frame = aFrame.clone();
	Size const new_id( new_frame->add_fragment( const_cast< FragData*>( frag_id.fragment_ptr().get() ) ) );
	new_frame->clone_cache_data( aFrame, frag_id.id(), new_id );
	add_( new_frame );
}

void
FragSet::generate_insert_map( MoveMap const& mm, InsertMap &insert_map, InsertSize &insert_size ) const {
	tr.Debug << "generate insert map from Movemap:\n";
	for ( Size i = 1; i<= max_pos(); i++) {
		if ( mm.get_bb( i ) ) tr.Debug << "*";
		else tr.Debug << "x";
	}
	tr.Debug << std::endl;
	typedef std::map< Size, Size> InsertSet;
	InsertSet insert_set;

	for ( ConstFrameIterator it=begin(), eit=end(); it!=eit; ++it ) {
		Size size ( it->is_valid() ? it->is_applicable( mm ) : 0 );
		if ( size ) {
			if ( insert_set[ it->start() ] < size ) insert_set[ it->start() ] = size;
		}
	}

	insert_map.clear();
	insert_size.clear();
	if( insert_set.size() != 0){
		insert_size.resize( insert_set.rbegin()->first ); //the largest residue in insert_map
	}
	// now copy it into a simple vector of numbers
	//for ( std::set< std::pair< Size, Size > >::const_iterator it=insert_set.begin(), eit=insert_set.end();
	for ( InsertSet::const_iterator it=insert_set.begin(), eit=insert_set.end();
				it!=eit; ++it ) {
		insert_map.push_back( it->first );
		insert_size[ it->first ] = it->second;
	}
}

Size FragSet::size() const {
	Size tot = 0;
	for ( ConstFrameIterator it=begin(), eit=end(); it!=eit; ++it ) {
		tot+=(*it)->nr_frags();
	}
	return tot; //frames_.size();
}

Size FragSet::nr_frames() const {
	Size tot = 0;
	for ( ConstFrameIterator it=begin(), eit=end(); it!=eit; ++it ) {
		tot+=1;
	}
	return tot; //frames_.size();
}

void
FragSet::add( FrameCOP aFrame ) {
	//std::cerr << "FragSet::add_frame " << std::endl;
	runtime_assert( aFrame->nr_frags() ); // do not allow insertion of empty frames --> makes frag_id iterator sooo much easier
	Size start ( aFrame->start() );
	Size end ( aFrame->end() );
	Size length( aFrame->length() );

	if ( min_pos() > start ) {
		set_min_pos( start );
	};

	if ( max_pos() < end ) {
		set_max_pos( end );
	};



	//	tr.Trace << "frag length " << length << " ( " << max_frag_length() <<  " ) " << std::endl;
	if ( length > max_frag_length() ) {
		tr.Trace << "set max frag length " << length << std::endl;
		set_max_frag_length( length );
	}

	// now add frame:
	FrameList present_frames;
	Size nr_present = frames( start, present_frames );
	if ( !nr_present ) {
		add_( aFrame->clone_with_frags() );
	} else {
		for ( FrameList::iterator it = present_frames.begin(),
						eit = present_frames.end(); it!=eit; ++it ) {
			if ( (*it)->is_mergeable( *aFrame ) ) {
				(*it)->merge( *aFrame );
				return; //finished early
			}
		}
		//didn't found mergable frames at this sequence position
		add_( aFrame->clone_with_frags() );
	}
}

void
FragSet::add( FrameList const& frames ) {
	for ( FrameList::const_iterator it=frames.begin(), eit=frames.end(); it!=eit; ++it ) {
		add( *it );
	}
}

void
FragSet::add( FragSet const& cframes ) {
	for ( ConstFrameIterator it=cframes.begin(), eit=cframes.end(); it!=eit; ++it ) {
		add( *it );
	}
}

std::ostream& operator<< (std::ostream& out, FragSet const& cfrags ) {
	for ( ConstFrameIterator it = cfrags.begin(), eit = cfrags.end(); it!=eit; ++it ) {
		out << *(*it);
	}
	return out;
}

void FragSet::shift_by( int offset ) {
	if ( offset != 0 ) {
    min_pos_ += offset;
    max_pos_ += offset;
    for ( FrameIterator it=nonconst_begin(), eit=nonconst_end(); it!=eit; ++it ) {
      //this should read ( *it != NULL ) but for some-reason gcc throws a warning that NULL needs to get converted to an unsigned int
      //somehow the returned FrameOP cannot be compared with NULL. calling .get() on the owning_pointer interface retrieves the naked pointer
      if ( (*it).get() != NULL ) {
        it->shift_by( offset );
      }
    }
  } else {
    tr.Debug << "Attempt to shift offset by " << offset << " ignored because it's 0." << std::endl;
  }
}

void FragSet::global_offset( int offset ){
	if ( offset != global_offset_ ) {
    shift_by( offset - global_offset_ );
    global_offset_ = offset;
    tr.Debug << "Shifted FragSet relative to fragment file position by " << offset << std::endl;
  } else {
    tr.Debug << "FragSets have not been shifted as they are already in place (offset: " << offset << ")." << std::endl;
  }
}


FragSetOP FragSet::clone_shifted( int offset ) const {
	FragSetOP newFragSet = empty_clone();

	for ( ConstFrameIterator it=begin(), eit=end(); it!=eit; ++it ) {
		FrameOP newFrame = (*it)->clone_with_frags();
		newFragSet->add( newFrame );
	}

  newFragSet->global_offset( offset );

  assert( ( newFragSet->max_pos() - newFragSet->min_pos() ) == ( this->max_pos() - this->min_pos() ) );

	return newFragSet;
}


} //fragment
} //core
