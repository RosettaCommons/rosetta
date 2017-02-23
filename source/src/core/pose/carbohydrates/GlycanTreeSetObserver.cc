// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/carbohydrates/GlycanTreeSetObserver.cc
/// @brief The CacheablePoseObserver version of GlycanTreeSet that will react to pose length changes..
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <core/pose/carbohydrates/GlycanTreeSetObserver.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>


#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "core.pose.carbohydrates.GlycanTreeSetObserver" );


namespace core {
namespace pose {
namespace carbohydrates {

	
	using namespace core::conformation::carbohydrates;
	using namespace core::pose::datacache;
	using namespace core::conformation::signals;
	using namespace core::pose::datacache;


///@brief
///
///  Get the GlycanTreeSet from the pose.
///  Returns NULLPTR if not setup or the data is not valid!
///
///@details
///
/// This is so we can get the glycan_tree_set using a const pose!
///
GlycanTreeSetCOP
get_glycan_tree_set(Pose const & pose){


	//const access: if cacheable observer hasn't been set, return NULL pointer
	if ( !pose.observer_cache().has( core::pose::datacache::GLYCAN_TREE_OBSERVER  ) ) {
		return nullptr;
	}
	else {
		
		//Not sure if this is slow or not and how it can be sped up.
		CacheableObserverCOP obs = pose.observer_cache().get_const_ptr( core::pose::datacache::CacheableObserverType::GLYCAN_TREE_OBSERVER );
		GlycanTreeSetObserverCOP observer = utility::pointer::static_pointer_cast< GlycanTreeSetObserver const >( obs );
		return observer->get_glycan_tree_set();

	}

}


GlycanTreeSetObserver::GlycanTreeSetObserver():
	CacheableObserver()
{

}

GlycanTreeSetObserver::GlycanTreeSetObserver( conformation::Conformation const & conf ):
	CacheableObserver()
{
	glycan_tree_set_ = GlycanTreeSetOP( new GlycanTreeSet( conf ) );
}

GlycanTreeSetObserver::GlycanTreeSetObserver( core::pose::Pose & pose ):
	CacheableObserver()
{
	glycan_tree_set_ = GlycanTreeSetOP( new GlycanTreeSet( pose.conformation() ) );
	attach_impl( pose );
}

GlycanTreeSetObserver::~GlycanTreeSetObserver(){
	detach_from();
}

GlycanTreeSetObserver::GlycanTreeSetObserver( GlycanTreeSetObserver const & observer):
	CacheableObserver( observer ),
	glycan_tree_set_( new GlycanTreeSet( *observer.glycan_tree_set_))

{

}

GlycanTreeSetObserverOP
GlycanTreeSetObserver::clone() const {
	return GlycanTreeSetObserverOP( new GlycanTreeSetObserver( *this ));
}

core::pose::datacache::CacheableObserverOP
GlycanTreeSetObserver::clone()
{
	return core::pose::datacache::CacheableObserverOP( new GlycanTreeSetObserver( *this ) );
}

core::pose::datacache::CacheableObserverOP
GlycanTreeSetObserver::create()
{
	return core::pose::datacache::CacheableObserverOP( new GlycanTreeSetObserver() );
}

GlycanTreeSetCOP
GlycanTreeSetObserver::get_glycan_tree_set() const {
	return glycan_tree_set_;
}

void
GlycanTreeSetObserver::on_length_change( core::conformation::signals::LengthEvent const & event ){
	glycan_tree_set_->on_length_change( event );
}




bool
GlycanTreeSetObserver::is_attached() const {
	return length_event_link_.valid();
}

void
GlycanTreeSetObserver::attach_impl( core::pose::Pose & pose ){
	
	length_event_link_ = pose.conformation().attach_length_obs( &GlycanTreeSetObserver::on_length_change, this );

}

void
GlycanTreeSetObserver::detach_impl(){
	length_event_link_.invalidate();
}

} //core
} //pose
} //carbohydrates






