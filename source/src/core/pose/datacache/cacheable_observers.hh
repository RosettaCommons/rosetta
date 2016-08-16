// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/datacache/cacheable_observers.hh
/// @brief  file that has class definitions for a bunch of general CacheableObserver implementations
///
/// @author Florian Richter (floric@u.washington.edu)

#ifndef INCLUDED_core_pose_datacache_cacheable_observers_hh
#define INCLUDED_core_pose_datacache_cacheable_observers_hh

// unit headers
#include <core/pose/datacache/cacheable_observers.fwd.hh>
#include <core/pose/datacache/CacheableObserver.hh>

// project headers
#include <core/conformation/signals/LengthEvent.hh>
//#include <core/pose/Pose.fwd.hh>


// utility headers
#include <utility/signals/Link.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace datacache {


/// @brief a cacheable observer that keeps track of what length events occured
class LengthEventCollector : public core::pose::datacache::CacheableObserver {


public: // typedefs


	typedef utility::signals::Link Link;


public: // construct/destruct


	/// @brief default constructor
	LengthEventCollector();


	/// @brief copy constructor
	/// @warning Subject being observed (represented by Link/pointer) is not copied!
	LengthEventCollector( LengthEventCollector const & rval );


	/// @brief default destructor
	/// @remarks detaches during destruction
	//virtual
	~LengthEventCollector();


public: // assignment


	/// @brief copy assignment
	/// @warning Subject being observed (represented by Link/pointer) is not copied!
	LengthEventCollector & operator =( LengthEventCollector const & rval );


public: // virtual constructors


	/// @brief clone this object
	/// @warning Subject (represented by Link/pointer) is not copied!
	pose::datacache::CacheableObserverOP clone();


	/// @brief create a new instance of this object
	pose::datacache::CacheableObserverOP create();


public: // interface


	/// @brief is this observer attached to a Pose/Conformation?
	bool is_attached() const {
		return length_event_link_.valid(); }

	void
	clear_events(){
		length_events_.clear(); }

	utility::vector1< core::conformation::signals::LengthEvent > const &
	events() const;


protected: // virtual observer interface


	/// @brief attach to Pose/Conformation
	virtual
	void attach_impl( pose::Pose & pose );


	/// @brief detach from Pose/Conformation
	virtual
	void detach_impl();

	void
	on_length_change( conformation::signals::LengthEvent const & event );

	//data
private:

	void
	copy_length_events( utility::vector1< core::conformation::signals::LengthEvent > const & events );

	utility::vector1< core::conformation::signals::LengthEvent > length_events_;

	Link length_event_link_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


/// @brief observer that tracks the fate of a one or more segments (i.e. pose
/// residues) of interest.
/// note: the convention should be that a segment.second marks the end of the segment
/// but is not part of it, i.e. the last position of a segment is segment.second - 1
/// reason: some peculiar stuff regarding the meaning of length events
class SpecialSegmentsObserver : public core::pose::datacache::CacheableObserver {

public: // typedefs
	typedef core::pose::datacache::CacheableObserver Parent;
	typedef utility::signals::Link Link;
	typedef core::Size Size;
	typedef std::pair< Size, Size > Segment;

public: // construct/destruct, assignment, virtual constructors

	SpecialSegmentsObserver();

	SpecialSegmentsObserver( SpecialSegmentsObserver const & rval );

	~SpecialSegmentsObserver();

	SpecialSegmentsObserver & operator =( SpecialSegmentsObserver const & rval );

	/// @brief clone this object
	/// @warning Subject (represented by Link/pointer) is not copied!
	pose::datacache::CacheableObserverOP clone();

	/// @brief create a new instance of this object
	pose::datacache::CacheableObserverOP create();

public: //interface

	/// @brief is this observer attached to a Pose/Conformation?
	bool is_attached() const {
		return length_event_link_.valid(); }

	utility::vector1< Segment > const &
	segments() const {
		return segments_;
	}

	void
	clear();

	void
	add_segment(
		Size begin,
		Size end
	);

	/// @brief utility function that sets all elements found in the
	/// SpecialSegmentsObserver in the pose to value
	static
	void
	set_farray_from_sso(
		ObjexxFCL::FArray1D_bool & array,
		core::pose::Pose const & pose,
		bool const value
	);

protected: // virtual observer interface

	void
	on_length_change( conformation::signals::LengthEvent const & event );

	/// @brief attach to Pose/Conformation
	virtual
	void attach_impl( pose::Pose & pose );


	/// @brief detach from Pose/Conformation
	virtual
	void detach_impl();

private: //data

	utility::vector1< Segment > segments_;

	Link length_event_link_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //SpecialSegmentsObserver

} // namespace datacache
} // namespace pose
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pose_datacache_cacheable_observers )
#endif // SERIALIZATION


#endif /* INCLUDED_core_pose_datacache_CacheableObserver_HH */
