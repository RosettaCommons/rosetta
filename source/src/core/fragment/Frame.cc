// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/Frame.cc
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
///

// Unit Headers
#include <core/fragment/Frame.hh>


// Package Headers

// Project Headers
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/id/SequenceMapping.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <basic/Tracer.hh>

#include <core/fragment/FragData.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

namespace core {
namespace fragment {

/// @details Auto-generated virtual destructor
Frame::~Frame() {}

using namespace kinematics;


static basic::Tracer tr("core.fragment");
//pose::PoseOP Frame::my_static_pose_for_testing_ = NULL;
void
make_pose_from_sequence_(
	std::string sequence,
	chemical::ResidueTypeSet const& residue_set,
	pose::Pose& pose
) {
	using namespace chemical;
	// clear all of the old data in the pose
	pose.clear();

	// setup the pose by appending the appropriate residues residues
	for ( Size seqpos = 1; seqpos <= sequence.length(); ++seqpos ) {
		char aa = sequence[seqpos-1]; // string indexing is zero-based!
		//		std::cerr<<  aa << " aminoacid requested" << std::endl;
		AA my_aa = aa_from_oneletter_code( aa );
		ResidueTypeCOPs const & rsd_type_list( residue_set.aa_map( my_aa ) );
		if ( rsd_type_list.size() == 0 ) {
			std::cerr << " cannot find aminoacid " << my_aa << " from sequence " << sequence << " in Frame::fragment_from_pose() " << std::endl;
		}
		runtime_assert( rsd_type_list.size() );
		Size best_index = 1;
		ResidueType const & rsd_type( *(rsd_type_list[ best_index ]) );
		conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( rsd_type ) );
		if ( seqpos == 1 ) {
			pose.append_residue_by_jump( *new_rsd, 1 );
		} else {
			pose.append_residue_by_bond( *new_rsd, true );
		}
	} // for seqpos

	//	if ( pose.total_residue() > 1)  pose.conformation().insert_chain_ending( pose.total_residue() - 1 );		// probably not necessary

} // make_pose_match_sequence_


Frame::Frame() :
	start_( 0 ),
	end_( 0 ),
	nr_res_( 0 )
{}

Frame::Frame( core::Size begin, core::Size end, core::Size nr_res ) :
	start_( begin ),
	end_( end ),
	nr_res_( nr_res )
{}

Frame::Frame( core::Size start, core::Size length ) :
	start_( start),
	end_ ( start + length - 1 ),
	nr_res_( length )
{}

Frame::Frame( core::Size start ) :
	start_( start),
	end_ ( 0 ), //needs to be determined later
	nr_res_( 0 )
{}

Frame::Frame( core::Size start, FragDataCOP const& frag1 ) :
	start_( start ),
	end_( start + frag1->size() - 1 ),
	nr_res_( frag1->size() )
{
	add_fragment( frag1 );
}

Frame::Frame( core::Size start, core::Size length, SingleResidueFragDataOP srfd )
	: start_( start ),
		end_( start + length - 1),
		nr_res_( length )
{
	add_fragment( new FragData( srfd, length ) );
}

/// @brief clone method, new frame with same alignment position, fragments are not copied!
FrameOP Frame::clone() const {
	return new Frame( start(), end(), length() );
}

/// @brief clone method, new frame with same alignment position, fragments are not copied!
FrameOP Frame::clone_with_frags() const {
	FrameOP newFrame = clone();// new Frame( start(), end(), length() );
	*newFrame = *this; //usually that is enough
	return newFrame;
}

/// @brief clone method, new frame with same alignment position, one fragments is copied as template ( valid() == false )
FrameOP Frame::clone_with_template() {
	FrameOP newFrame = clone();// new Frame( start(), end(), length() );
	if ( nr_frags() ) {
		FragDataOP new_frag_data=frag_list_[ 1 ]->clone();
		new_frag_data->set_valid( false );
		newFrame->frag_list_.push_back( new_frag_data );
	}
	return newFrame;
}

///@brief type() is specifying the output name of the Frame in FragmentIO ("FRAME", "JUMPFRAME", etc)
std::string Frame::type() const {
	return _static_type_name();
}

std::string Frame::_static_type_name() {
	return "FRAME";
}


/// @brief accesors for underlying FragData
FragData const & Frame::fragment( core::Size frag_num ) const {
	return *frag_list_[ frag_num ];
}

/// @brief accessor for underlying FragData
//FragData const& Frame::fragment( core::Size frag_num ) {
//	return *frag_list_[ frag_num ];
//}

/// @brief accessor for underlying FragData as owning ptr
FragDataCOP Frame::fragment_ptr( core::Size frag_num ) const {
	return frag_list_[ frag_num ];
}

/// @brief accessor for underlying FragData as owning ptr
//FragDataOP Frame::fragment_ptr( core::Size frag_num ) {
	//return frag_list_[ frag_num ];
//}

/// @brief a frame is considered valid if at least one fragment is contained and this fragment is also valid
/// (not an empty template fragment)
bool Frame::is_valid() const {
 // why check both 1 and 2?
	return ( nr_frags() >= 2 && fragment( 2 ).is_valid() )
		||   ( nr_frags() >= 1 && fragment( 1 ).is_valid() );
}

/// @brief translate intra-frame position into sequence position. (trivial for base-class)
core::Size Frame::seqpos( core::Size intra_pos ) const { // BaseClass --> continuous frames
	return start()+intra_pos-1;
}

/// @brief  a unique id for every fragment in the list.
/// his is silly, but would enable later on to avoid cache_clearence on deletion of FragData entries
/// in this case, we would require that the ID of a certain fragment never changes, even if the position in FragList changes
core::Size Frame::frag_id( core::Size frag_num ) const
{ return frag_num; }


/// @brief true if frame is continuous (always true for base class)
bool Frame::is_continuous() const
{ return true; } // base class can only handle continuous frames

/// @brief number of fragments attached to this frame
core::Size Frame::nr_frags() const {
	return frag_list_.size();
}


/// @brief whether this fragment contains a certain position
bool Frame::contains_seqpos( core::Size seqpos ) const {
	return (seqpos >= start_ ) && (seqpos <= end_ ); }

/// @brief first seqpos of this frame
core::Size Frame::start() const {
	return start_;
}


/// @brief last sequence position affected by this frame
core::Size Frame::end() const {
	return end_;
}

///	///@brief set stop position
//	core::Size stop( core::Size setting );

/// @brief last sequence position affected by this frame
core::Size Frame::stop() const {
	return end_;
}


bool Frame::moves_residue( core::Size pos ) const {
	return pos <= stop() && pos >= start();
}
/// @brief number of residues affected by this frame
core::Size Frame::nr_res_affected( kinematics::MoveMap const& mm ) const {
	return is_applicable( mm );
}

/// @brief number of residues in this frame ( for continuous frames it is the same as end()-start() + 1 )
core::Size Frame::length() const {
	return nr_res_;
}


// can we make this private and accessible via
// a friend BaseCacheUnit statement --> are derived classes of a friend still a friend
/// @brief return handle to cached data stored under "tag"
/// shouldn't be called directly
BaseCacheUnit& Frame::cache( std::string tag, BaseCacheUnitOP const& new_cache )  const {
	CacheMap::const_iterator iter( cache_.find( tag ) );
	if ( iter == cache_.end() ) {
		return *( cache_[ tag ] = new_cache->clone() );
	}
	return *( cache_[ tag ] );
}

/// @brief copies all entries in the Frame-Cache for fragment "sid" of Frame "source" to fragment "nid" of "this" frame
void Frame::clone_cache_data( Frame const& source, core::Size sid, core::Size nid ) {
	for ( CacheMap::iterator it=source.cache_.begin(), eit=source.cache_.end(); it!=eit; ++it ) {
		cache( it->first, it->second  ).remap_value( *it->second, sid, nid );
	}
}


//@brief prints frame to stream -- multiline object


bool Frame::is_mergeable( Frame const& other ) const {
	if ( !( other.length() == length() && other.start() == start() && other.stop() == stop() )) return false;
	if ( other.nr_frags()==0 || nr_frags()==0 ) return true;
	if ( other.fragment( 1 ).is_compatible( fragment( 1 ) )) return true;
	return false;
}


/// @brief is a FragData object compatible with the already stored ones ?
/// @detail you can only add instances of FragData to the same Frame that are compatible, i.e., that contain the same
/// class of FragData, e.g., based on BBTorsionSRFD,
/// if you want to have different fragment for other dof's at the same sequence position create a new Frame.
/// Users of the fragment-core are aware that multiple Frames for the same sequence position may exist.
bool Frame::is_compatible( FragDataCOP new_frag ) const {
	if ( nr_frags() ) {
		return new_frag->is_compatible( fragment( 1 ) );
	} else {
		if ( length() == 0 ) {
			nr_res_ = new_frag->size();
			end_ = start_ + nr_res_ - 1;
			return true;
		} else return (new_frag->size() == length() );
	}
	/*			if ( frag_list_.size() ) {
			kinematics::MoveMap new_move_map;
			new_frag->generate_move_map( new_move_map, *this );
			return true;
			//			return move_map_ == new_move_map;
		} else {
			new_frag->generate_move_map( move_map_, *this );
			return true;
		}
	*/
}

void Frame::init_length( core::Size start, core::Size end, core::Size length ) {
	start_ = start;
	end_ = end;
	nr_res_ = length;
}







////////////////////////////////////////// I M P L E M E N T A T I O N S ///////////////////////////////////////////////

core::Size Frame::add_fragment( FragDataCOP new_frag ) {
	assert( new_frag );
	bool success ( is_compatible( new_frag ) );
	if ( success ) frag_list_.push_back( new_frag );
	return frag_list_.size();
}


//
// void Frame::delete_fragment( core::Size frag_num )
// {
// 	std::cerr << "SPEED WARNING: fragment deleted -- cache invalidated"
// 						<< "for now we  use frag_num as frag_id, so all cached data is invalidated" << std::endl;
// 	cache_.clear();

// 	//and now delete fragment...
// 	std::cerr << "delete stubbed out" << frag_num << std::endl;
// }

core::Size Frame::is_applicable( kinematics::MoveMap const& mm ) const {
	//	if ( nr_frags()==0 ) return true;
	if ( is_continuous() ) {
		return fragment( 1 ).is_applicable( mm, start(), end() );
	} else {
		return fragment( 1 ).is_applicable( mm, *this );
	}
}

core::Size Frame::apply( kinematics::MoveMap const& mm, core::Size frag_num, pose::Pose & pose ) const {
	if ( is_continuous() ) {
		return fragment( frag_num ).apply( mm, pose, start(), end() );
	} else {
		return fragment( frag_num ).apply( mm, pose, *this );
	}
}

core::Size Frame::apply( core::Size frag_num, pose::Pose & pose ) const {
	if ( is_continuous() ) {
		return fragment( frag_num ).apply( pose, start(), end() );
	} else {
		return fragment( frag_num ).apply( pose, *this );
	}
}



core::Size Frame::apply_ss( kinematics::MoveMap const& mm, core::Size frag_num, std::string& ss ) const {
	return fragment( frag_num ).apply_ss( mm, ss, *this );
}


void Frame::shift_to( core::Size setting ) {
	using core::Size;
	// set end first
	if ( start_ < setting ) {
		end_ += ( setting - start_ );
	} else { // start >= setting
		end_ -= ( start_ - setting );
	}
	// set start second
	start_ = setting;
	assert( nr_res_ == ( end_ - start_ + 1 ) ); //OL: changed assert to use == instead of "=",  much better now
}

void Frame::shift_by( int offset ) {
	if ( start_ + offset < 1 ) {
		std::ostringstream msg;
		msg << "offset " << offset << " would shift Frame " << *this << " to negative or zero starting position " << std::endl;
		throw utility::excn::EXCN_RangeError( msg.str() );
	}
	start_ += offset;
	end_ += offset;
}



// void show( std::ostream& out ) {
//
// }




void Frame::fragment_as_pose(
	Size frag_num,
	pose::Pose & pose,
	chemical::ResidueTypeSetCAP restype_set ) const
{
	//	if (!my_static_pose_for_testing_ ) {
	//		std::string const pdbfile ( "protocols/abinitio/2GB3.pdb" );
	//my_static_pose_for_testing_=new pose::Pose;
	pose.clear();
	make_pose_from_sequence_( frag_list_[ frag_num ]->sequence(),
		*restype_set,
		pose );
		//		core::import_pose::pose_from_pdb( *my_static_pose_for_testing_, pdbfile );
	fragment( frag_num ).apply( pose, 1, length() );
}

bool Frame::add_fragment( FragDataCOPs new_frag ) {
	for ( FragDataCOPs::const_iterator it=new_frag.begin(), eit=new_frag.end();
				it!=eit; ++it ) {
		bool success ( add_fragment( *it ) );
		if ( !success ) return false;
	}
	return true;
}

void Frame::show_classic( std::ostream &out ) const {
	using namespace ObjexxFCL::format;
	Size position = start();
	out << "position: " << RJ( 10, position) << " neighbors: " << RJ( 10, nr_frags() ) << std::endl << std::endl;
	for ( Size nr = 1; nr <= nr_frags(); nr ++ ) {
		fragment( nr ).show( out );
		out << std::endl;
	}
}

void Frame::show( std::ostream &out ) const {
// 	using namespace ObjexxFCL::format;
// 	out << "FRAME " << " " << RJ( 10, start() ) << RJ( 10, end() ) << std::endl;
	show_header( out );
	show_fragments( out );
}

void Frame::show_header( std::ostream &out ) const {
	using namespace ObjexxFCL::format;
	out << "FRAME " << " " << RJ( 3, start() ) << " " << RJ( 3, end() ) << std::endl;
}

void Frame::read( std::istream &in ) {
	//because of one FragData - multiple frames.. reading of fragments is outsourced to FragmentIO
	//	std::string tag;
	in >> start_ >> end_;
	nr_res_ = end_ - start_ + 1;
}

void Frame::show_fragments( std::ostream& out ) const {
	for ( Size nr = 1; nr <= nr_frags(); nr ++ ) {
		runtime_assert( fragment_ptr( nr ) );
		//		std::cerr << "FRAME::show_fragments " << nr << " seqpos:" << start() << std::endl;
		fragment( nr ).show( out, *this );
		out << std::endl << std::endl;
		//		std::cerr << "FRAME::show_fragments done" << std::endl;
	}
}

bool Frame::merge( Frame const& other ) {
	if ( !is_mergeable( other) ) return false;
	//append fragdata
	//copy cached data
	Size insert_pos = frag_list_.size()+1;
	Size other_frag_num = 1;
	for ( FragDataCOPs::const_iterator it = other.frag_list_.begin(),
					eit = other.frag_list_.end(); it!=eit; ++it ) {
		frag_list_.push_back( *it );
		clone_cache_data( other, other_frag_num++, insert_pos++ );
	}
	return true;
}

///@brief change frames residue numbers according to map
bool Frame::align( core::id::SequenceMapping const& map) {
	Size s( map[ start() ] );
  Size e( map[ end() ] );
	bool success = ( s>0 && e>0 );
	if ( is_continuous() && success ) {
		success = ( (e - s + 1 ) == length() );
		if ( !success ) tr.Trace << "failed to align frame " << start() << " " << end() << " different size after alignment" << std::endl;
		//check that all residues within frame can be mapped -- or is that an unwanted restriction ?
		for ( Size pos = start(); success && pos<=end(); pos ++ ) { // this should only happen if we hae insertion and deletion cancelling each other  ?
			success = map[ pos ] > 0;
			if ( !success ) tr.Trace << "failed to align frame " << start() << " " << end() << " due to pos " << pos << std::endl;
		}
	}
	if ( success ) {
		start_ = s;
		end_ = e;
	}
	return success;
}

FrameOP Frame::generate_sub_frame(  Size length, Size start /* = 1*/ ) const {
	if ( !is_continuous() ) return NULL;
	FrameOP sub_frame( clone() ); //creates empty frame with same type
	sub_frame->start_ += start - 1;
	sub_frame->nr_res_ = length;
	sub_frame->end_ = sub_frame->start_ + length - 1;
	for ( FragDataCOPs::const_iterator it = frag_list_.begin(), eit = frag_list_.end();
				it!=eit; ++it ) {
		sub_frame->add_fragment( (*it)->generate_sub_fragment( start, start + length - 1 ) );
	}
	runtime_assert( nr_frags() == sub_frame->nr_frags() );
	return sub_frame;
}

bool Frame::steal( pose::Pose const& pose) {
	runtime_assert( nr_frags() );
	FragDataOP new_frag = frag_list_.front()->clone();
	bool success ( new_frag->steal( pose, *this ) );
	if ( success ) {
		// check if first fragment was just an unitialized template...
		if ( frag_list_.front()->is_valid() ) {
			frag_list_.push_back( new_frag );
		} else {  // remove unitialized template and replace it with the new fragment
			tr.Trace << "remove unitialized template FragData and replace it with stolen one" << std::endl;
			runtime_assert( nr_frags() == 1 );
			frag_list_.clear();
			frag_list_.push_back( new_frag );
		}
	}
	return success;
}

void Frame::clear() {
	cache_.clear();
	if ( nr_frags() ) {
		FragDataOP frag_template = frag_list_.front()->clone();
		frag_list_.clear();
		frag_template->set_valid( false ); //make it a template
		add_fragment( frag_template );
	}
}

}
}
