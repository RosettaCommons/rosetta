// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/Frame.hh
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @author James Thompson
/// @date   Wed Oct 20 12:08:31 2007
///
#ifndef INCLUDED_core_fragment_Frame_HH
#define INCLUDED_core_fragment_Frame_HH

// Unit Headers
#include <core/fragment/Frame.fwd.hh>

// Package Headers
#include <core/fragment/BaseCacheUnit.hh>

// Project Headers
#include <core/fragment/SingleResidueFragData.fwd.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
//#include <core/io/pdb/pose_io.hh> //needed for hack in fragment_as_pose
#include <core/id/SequenceMapping.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#ifdef BOINC
#include <utility/fix_boinc_read.hh>

#endif

#ifdef WIN32
	#include <string>
	#include <core/fragment/FragData.hh>
#endif

// C++ STL Headers
#include <map>

#include <core/types.hh>
#include <core/fragment/FragData.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace fragment {

/* How to make a top_N_frag thingy, i.e., get a Frame that does only use a certain subset of the
	 available FragData
	 clearly that will be a copy of a Frame with less elements.
		tex: not necessarily. we could just set weights of certain fragments to zero and do weight-based
		sampling.

	 Frame contains only FragDataCOPs so it can refer to the same elements but less of them ...

	 ERROR Handling: have to think about that. So far I often return a bool
	 Raising Exceptions would be better. But I think we are still free of exceptions...
*/
//WARNING: when fragments inside of a frame are deleted the FragID of that frame will be fucked up


// if deleting of fragments in a frame is really needed a lot one could
// introduce a second class with a slightly different implementation
// DeletableFrame that maps frag_id --> frag_nr
// but for most purposes it will probably suffice to make a new frame that contains all wanted fragments


/// @brief Frame couples a list of FragData-instances to a certain alignment frame, i.e., position in sequence space
/// @detail A frame may be continous, i.e., its fragment data will be applied to all residues between start() and end()
/// or ( in a specialized version inheriting this interface) may contain a loose list of sequence positions.
///
/// Instances of FragData (e.g., the Torsion-Angles of a specific Fragment) can be
/// added and deleted via add_- or delete_fragment() methods.
///
/// a fragment is inserted into the structure via the apply() method
/// and the steal() method is its inverse: an instance of FragData is created from the structure of the pose
/// at the dof's specified in the other FragData's of the Frame. The new FragData is added as last fragment to the frame.
///
/// accessors to underlying FragData-instances are available as fragment() or fragment_ptr().
///
/// Frame supports the FragCache --> see FragCache for documentation.
///
/// MoveMaps: It should be possible to specify which dofs can be moved. In this case fragment-insertions should
/// a) only be allowed if all dofs affected are moveable or b) only change those dofs that are moveable.
/// this information could be handled by an extended version of the kinematics::movemap.
/// The movemap would then have to functions: 1) which dofs are changed by a certain FragData/Frame
///                                           2) which dofs are moveable
///                             comparison of movemaps, i.e., M1 <= M2 could tell us if dofs tagged in M1 are also tagged in M2:
///                       i.e., the M_fragdata<= M_pose would tell us that the fragment is applicable to the Pose.


class Frame : public utility::pointer::ReferenceCount {
	typedef	std::map<std::string, BaseCacheUnitOP > CacheMap;
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~Frame();

	Frame();

	Frame( core::Size begin, core::Size end, core::Size nr_res );

	Frame( core::Size start, core::Size length );

	Frame( core::Size start );

	Frame( core::Size start, FragDataCOP const& frag1 );

	Frame( core::Size start, core::Size length, SingleResidueFragDataOP srfd );

	/// @brief clone method, new frame with same alignment position, fragments are not copied!
	virtual FrameOP clone() const;

	/// @brief clone method, new frame with same alignment position, fragments are not copied!
	virtual FrameOP clone_with_frags() const;

	/// @brief clone method, new frame with same alignment position, one fragments is copied as template ( valid() == false )
	virtual FrameOP clone_with_template();

	/// @brief type() is specifying the output name of the Frame in FragmentIO ("FRAME", "JUMPFRAME", etc)
	virtual std::string type() const;

	static std::string _static_type_name();


	/// @brief add a fragment .. return new frag_nr
	core::Size add_fragment( FragDataCOP new_frag );

	/// @brief add all fragments in list
	bool add_fragment( FragDataCOPs new_frags );

	/// @brief delete a fragment: Attention: all data in the FragCache is invalidated ( and deleted )
	/// it would be complicated to change this behaviour. Thus, it is desirable to avoid using delete_fragment() altogether.
	/// Faster: Make a new Frame and add all fragments you are interested in.
	//void delete_fragment( core::Size frag_num );

	/// @brief accesors for underlying FragData
	FragData const & fragment( core::Size frag_num ) const;

	/// @brief accessor for underlying FragData
	//FragData & fragment( core::Size frag_num );

	/// @brief accessor for underlying FragData as owning ptr
	FragDataCOP fragment_ptr( core::Size frag_num ) const;

	/// @brief accessor for underlying FragData as owning ptr
	//FragDataOP fragment_ptr( core::Size frag_num );

	/// @brief a frame is considered valid if at least one fragment is contained and this fragment is also valid
	/// (not an empty template fragment)
	bool is_valid() const;


	/// @brief insert fragment frag_num into pose
	core::Size apply( kinematics::MoveMap const&, core::Size frag_num, pose::Pose & pose ) const;

	/// @brief insert fragment frag_num into pose --- ignore movemap
	core::Size apply( core::Size frag_num, pose::Pose & pose ) const;

	/// @brief change ss-string according to sec-struct info in fragments
	core::Size apply_ss( kinematics::MoveMap const&, core::Size frag_num, std::string& ss ) const;

	/// @brief obtain fragment from pose at frame position
	// ATTENTION: functionality of steal is different for FRAMES than for FragData.
	// since many frames can own the same FragData, we make a new clone of FragData for the stolen frame
	bool steal( pose::Pose const& pose);

	/// @brief is the Frame applicable to the pose with the given movemap?
	core::Size is_applicable( kinematics::MoveMap const& /*trial_move_map*/ ) const;


	/// @brief translate intra-frame position into sequence position. (trivial for base-class)
	virtual core::Size seqpos( core::Size intra_pos ) const; // BaseClass --> continuous frames

	/// @brief  a unique id for every fragment in the list.
	/// his is silly, but would enable later on to avoid cache_clearence on deletion of FragData entries
	/// in this case, we would require that the ID of a certain fragment never changes, even if the position in FragList changes
	core::Size frag_id( core::Size frag_num ) const;

	/// @brief returns a (small) pose with fragment ( continous from seqpos 1 ... nr_res_affected() )
	void fragment_as_pose(
		core::Size frag_num,
		pose::Pose & pose,
		chemical::ResidueTypeSetCAP restype_set
	) const;

	////////// properties of the Frame///////////////////

	/// @brief true if frame is continuous (always true for base class)
	virtual bool is_continuous() const; // base class can only handle continuous frames

	/// @brief number of fragments attached to this frame
	core::Size nr_frags() const;

	/// @brief remove all valid fragments, a template fragment will be left over
	void clear();

	/// @brief whether this fragment contains a certain position
	bool contains_seqpos( core::Size seqpos ) const;

	/// @brief first seqpos of this frame
	core::Size start() const;

	/// @brief set start position
	//	inline
	//	void start( core::Size setting );

	/// @brief shift to new start position ( change end accordingly )
	virtual void shift_to( core::Size setting );

	/// @brief shift frame by offset relative to current start position ( change end accordingly )
	virtual void shift_by( int offset );

	/// @brief last sequence position affected by this frame
	core::Size end() const;

	///	/// @brief set stop position
	//	core::Size stop( core::Size setting );

	/// @brief last sequence position affected by this frame
	core::Size stop() const;

	virtual bool moves_residue( core::Size pos ) const;

/// @brief number of residues affected by this frame
	core::Size nr_res_affected( kinematics::MoveMap const& mm ) const;

	/// @brief number of residues in this frame ( for continuous frames it is the same as end()-start() + 1 )
	core::Size length() const;


	// can we make this private and accessible via
	// a friend BaseCacheUnit statement --> are derived classes of a friend still a friend
	/// @brief return handle to cached data stored under "tag"
	/// shouldn't be called directly
	BaseCacheUnit& cache( std::string tag, BaseCacheUnitOP const& new_cache )  const;

	/// @brief copies all entries in the Frame-Cache for fragment "sid" of Frame "source" to fragment "nid" of "this" frame
	void clone_cache_data( Frame const& source, core::Size sid, core::Size nid );

	virtual
	void show_classic( std::ostream& out ) const;

	//@brief prints frame to stream -- multiline object
	virtual
	void show( std::ostream& out ) const;

	virtual
	void show_header( std::ostream& out ) const;

	virtual
	void read( std::istream& in );

	bool is_mergeable( Frame const& other ) const;

	bool merge( Frame const& other );

	/// @brief change frames residue numbers accoriding to map
	virtual
	bool align( core::id::SequenceMapping const& map );

	/// @brief generate_sub_frame of length from start ( internal numbers )
	FrameOP generate_sub_frame( Size length, Size start = 1 ) const;

	/// @brief NOT IMPLEMENTED YET: generate_sub_frame according to mapping ( residue numbers ) returns NULL if mapping invalid
	// Commenting out to make Python bindings compile
	//FrameOP generate_sub_frame( core::id::SequenceMapping const& map );

protected:

	//@brief called by show() to iterate over fragments
	virtual
	void show_fragments( std::ostream& out ) const;

	/// @brief is a FragData object compatible with the already stored ones ?
	/// @detail you can only add instances of FragData to the same Frame that are compatible, i.e., that contain the same
	/// class of FragData, e.g., based on BBTorsionSRFD,
	/// if you want to have different fragment for other dof's at the same sequence position create a new Frame.
	/// Users of the fragment-core are aware that multiple Frames for the same sequence position may exist.
	virtual bool is_compatible( FragDataCOP new_frag ) const;

	void init_length( core::Size start, core::Size end, core::Size length );


private:
	// first seqpos of frame
	core::Size start_;

	// last seqpos of frame
	mutable core::Size end_;

	// nr of residues affected ( start_ - end_ + 1, for continuous frames )
	mutable core::Size nr_res_;

	// not used right now..
	//kinematics::MoveMap move_map_;

	// cache data for FragCache and FragStore functionality
	mutable CacheMap cache_;

	// a list of fragments for this frame
	FragDataCOPs frag_list_;

	//static pose::PoseOP my_static_pose_for_testing_; //replace that with something more sensible ...
};

inline std::ostream& operator<< ( std::ostream& out, Frame const& frame ) {
	frame.show( out );
	return out;
}


} //fragment
} //core

#endif
