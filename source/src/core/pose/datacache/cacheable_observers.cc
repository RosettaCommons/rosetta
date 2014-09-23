// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/datacache/cacheable_observers.cc
/// @brief  a bunch of cacheable observers
///
/// @author Florian Richter ( floric@u.washington.edu)

// Unit headers
#include <core/pose/datacache/cacheable_observers.hh>

// Package headers
// AUTO-REMOVED #include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>

// Project headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/pose/Pose.hh>

#include <ObjexxFCL/FArray1D.hh>

//c++ headers
#include <set>

#include <utility/vector1.hh>
#include <boost/foreach.hpp>

//Auto Headers
namespace core {
namespace pose {
namespace datacache {


/// @brief default constructor
LengthEventCollector::LengthEventCollector() :
	CacheableObserver()
{
	length_events_.clear();
}


/// @brief copy constructor
LengthEventCollector::LengthEventCollector( LengthEventCollector const & rval ) :
	CacheableObserver( rval )
{
	copy_length_events( rval.length_events_ );
}

LengthEventCollector::~LengthEventCollector() {
	detach_from();
}


/// @brief copy assignment
LengthEventCollector &
LengthEventCollector::operator =( LengthEventCollector const & rval ) {
	if ( this != &rval ) {
		CacheableObserver::operator =( rval );

		copy_length_events( rval.length_events_ );
	}

	return *this;
}

utility::vector1< core::conformation::signals::LengthEvent > const &
LengthEventCollector::events() const{
	return length_events_;
}


pose::datacache::CacheableObserverOP
LengthEventCollector::create()
{
	return pose::datacache::CacheableObserverOP( new LengthEventCollector() );
}

pose::datacache::CacheableObserverOP
LengthEventCollector::clone()
{
	return pose::datacache::CacheableObserverOP( new LengthEventCollector( *this ) );
}


void
LengthEventCollector::copy_length_events(
	utility::vector1< core::conformation::signals::LengthEvent > const & events
){

	using namespace core::conformation::signals;
	length_events_.clear();

	BOOST_FOREACH(LengthEvent event, events){
		length_events_.push_back( LengthEvent( event ) );
	}

}

/// @details all this class does is keep track of length events
void
LengthEventCollector::on_length_change( conformation::signals::LengthEvent const & event ){

	length_events_.push_back( conformation::signals::LengthEvent( event ) );
}

void
LengthEventCollector::attach_impl( pose::Pose & pose ){

	length_event_link_ = pose.conformation().attach_length_obs( &LengthEventCollector::on_length_change, this );

}

void
LengthEventCollector::detach_impl(){

	length_event_link_.invalidate();

}


SpecialSegmentsObserver::SpecialSegmentsObserver()
: Parent() {}

SpecialSegmentsObserver::SpecialSegmentsObserver( SpecialSegmentsObserver const & rval )
	: Parent( rval ), segments_(rval.segments_ ) {}

SpecialSegmentsObserver::~SpecialSegmentsObserver() {
	detach_from(); }

/// @brief copy assignment
SpecialSegmentsObserver &
SpecialSegmentsObserver::operator =( SpecialSegmentsObserver const & rval ) {
	if ( this != &rval ) {
		CacheableObserver::operator =( rval );

		segments_ = rval.segments_;
	}

	return *this;
}

pose::datacache::CacheableObserverOP
SpecialSegmentsObserver::create()
{
	return pose::datacache::CacheableObserverOP( new SpecialSegmentsObserver() );
}

pose::datacache::CacheableObserverOP
SpecialSegmentsObserver::clone()
{
	return pose::datacache::CacheableObserverOP( new SpecialSegmentsObserver( *this ) );
}

void
SpecialSegmentsObserver::clear()
{
	segments_.clear();
}

void
SpecialSegmentsObserver::add_segment(
	Size begin,
	Size end
){
	//dummy check
	if( begin > end ) utility_exit_with_message("Negative length segment added to SpecialSegmentsObserver.");

	segments_.push_back( std::pair< Size, Size >(begin, end ) );
}


void
SpecialSegmentsObserver::set_farray_from_sso(
	ObjexxFCL::FArray1D_bool & array,
	core::pose::Pose const & pose,
	bool const value
){

	if( pose.observer_cache().has( core::pose::datacache::CacheableObserverType::SPECIAL_SEGMENTS_OBSERVER) ){
		utility::vector1< std::pair< core::Size, core::Size > > const & segments = utility::pointer::static_pointer_cast< core::pose::datacache::SpecialSegmentsObserver const >(pose.observer_cache().get_const_ptr( core::pose::datacache::CacheableObserverType::SPECIAL_SEGMENTS_OBSERVER ) )->segments();


		for( core::Size i = 1; i <= segments.size(); ++i ){
			for( core::Size j = segments[i].first; j < segments[i].second; ++j ) array[ j - 1 ] = value;
		}
	}
}

/// @details figure out where the loops are after the length has changed
/// note: in case the segment got deleted, it will also be removed from
/// the observer
void
SpecialSegmentsObserver::on_length_change( conformation::signals::LengthEvent const & event ){

	std::set< Size > deleted_segments; //in case stuff gets deleted, we have to rearrange the segment vector
	//std::cerr << "SSO getting hit with event at pos " << event.position << " of change " << event.length_change << "   ";
	for( Size i(1); i<= segments_.size(); ++i ){
		//3 possibilites:
		//1. event downstream of segment, segment remains unchanged
		//2. event upstream of segment, segment gets shifted or (partially) deleted
		//3. event in segment, segment gets lengthened/shortened
		//std::cerr << "(" << segments_[i].first << "," << segments_[i].second << ")" << "moves to ";

		//1.
		if( event.position > segments_[i].second ){} //event happened downstream of the segment, i.e. segment unchanged

		//2.
		else if( event.position < segments_[i].first ){

			if( event.length_change < 0 ) { //deletion?

				if( int( event.position - event.length_change ) > int(segments_[i].second) ){ //complete segment deleted?
					deleted_segments.insert(i);
				}
				else if(  int( event.position - event.length_change ) > int(segments_[i].first) ){ //partial segment deleted?
					segments_[i].first = event.position - event.length_change;
					segments_[i].second += event.length_change;
				}
				else{ //deletion before segment
					segments_[i].first += event.length_change;
					segments_[i].second += event.length_change;
				}
			} //if length_change < 0 //deletion
			else{ //insertion before segment
				segments_[i].first += event.length_change;
				segments_[i].second += event.length_change;
			}
		} //if event.position < segments_[i].first

		//3. event in segment
		else{
			if( event.length_change < 0 ){ //deletion?
				if( int( event.position - event.length_change ) > int(segments_[i].second) ){ //deletion past end of segment?
					segments_[i].second = event.position;
				}
				else{ //deletion within segment
					segments_[i].second += event.length_change;
				}
			} //if deletion
			else{ //insertion in segment
				segments_[i].second += event.length_change;
			}
		}
		//std::cerr << "(" << segments_[i].first << "," << segments_[i].second << ")      ";
	}//loop over segments
	//std::cerr << std::endl;

	//in case complete segments have been deleted, we need to restructure the vector
	if( deleted_segments.size() > 0 ){
		utility::vector1< Segment > new_segments;
		for( Size i = 1; i <= segments_.size(); ++i){
			if( deleted_segments.find( i ) != deleted_segments.end() ) new_segments.push_back( segments_[i] );
		}
		segments_ = new_segments;
	}
}


void
SpecialSegmentsObserver::attach_impl( pose::Pose & pose ){

	length_event_link_ = pose.conformation().attach_length_obs( &SpecialSegmentsObserver::on_length_change, this );

}

void
SpecialSegmentsObserver::detach_impl(){

	length_event_link_.invalidate();

}



} // namespace datacache
} // namespace pose
} // namespace core
