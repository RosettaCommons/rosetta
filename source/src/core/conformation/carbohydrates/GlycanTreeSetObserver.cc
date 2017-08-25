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

#include <core/conformation/carbohydrates/GlycanTreeSetObserver.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/conformation/signals/LengthEvent.hh>


#include <basic/Tracer.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION



static THREAD_LOCAL basic::Tracer TR( "core.conformation.carbohydrates.GlycanTreeSetObserver" );


namespace core {
namespace conformation {
namespace carbohydrates {


using namespace core::conformation::carbohydrates;
using namespace core::conformation::signals;





GlycanTreeSetObserver::GlycanTreeSetObserver():
	utility::pointer::ReferenceCount()
{

}

GlycanTreeSetObserver::GlycanTreeSetObserver( conformation::Conformation const & conf ):
	utility::pointer::ReferenceCount()
{
	glycan_tree_set_ = GlycanTreeSetOP( new GlycanTreeSet( conf ) );
}

/*
GlycanTreeSetObserver::GlycanTreeSetObserver( core::pose::Pose & pose ):
	CacheableObserver()
{
	glycan_tree_set_ = GlycanTreeSetOP( new GlycanTreeSet( pose.conformation() ) );
	attach_impl( pose );
}
*/


GlycanTreeSetObserver::~GlycanTreeSetObserver(){
	detach_impl();
}

GlycanTreeSetObserver::GlycanTreeSetObserver( GlycanTreeSetObserver const & observer):
	utility::pointer::ReferenceCount(),
	glycan_tree_set_( new GlycanTreeSet( *observer.glycan_tree_set_))

{
	//Detach on clone
	detach_impl();
}

GlycanTreeSetObserverOP
GlycanTreeSetObserver::clone() const {
	return GlycanTreeSetObserverOP( new GlycanTreeSetObserver( *this ));
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
GlycanTreeSetObserver::attach_to( Conformation & conf ){
	detach_impl();
	attach_impl( conf );
}
	
	
void
GlycanTreeSetObserver::attach_impl( Conformation & conf ){
	length_event_link_ = conf.attach_length_obs( &GlycanTreeSetObserver::on_length_change, this );

}

void
GlycanTreeSetObserver::detach_impl(){
	length_event_link_.invalidate();
}

} //core
} //conformation
} //carbohydrates

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::carbohydrates::GlycanTreeSetObserver::save( Archive & arc ) const {
	// EXEMPT length_event_link_
	arc( CEREAL_NVP( glycan_tree_set_ ) ); // core::Real

}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::carbohydrates::GlycanTreeSetObserver::load( Archive & arc ) {
	// EXEMPT length_event_link_
	arc( glycan_tree_set_ ); // core::Real

}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::carbohydrates::GlycanTreeSetObserver );
CEREAL_REGISTER_TYPE( core::conformation::carbohydrates::GlycanTreeSetObserver )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_carbohydrates_GlycanTreeSetObserver )
#endif // SERIALIZATION




