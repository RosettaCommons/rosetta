// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/denovo_design/components/StructureDataObserver.cc
/// @brief  Observes a pose and updates StructureData accordingly
///
/// @author Tom Linsky (tlinsky at uw dot edu)

// Unit headers
#include <protocols/denovo_design/components/StructureDataObserver.hh>

// Core headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// Boost/ObjexxFCL headers
#include <ObjexxFCL/FArray1D.hh>
#include <boost/foreach.hpp>

// C++ headers
#include <set>

//Auto Headers
#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.components.StructureDataObserver" );

namespace protocols {
namespace denovo_design {
namespace components {

StructureDataObserver::StructureDataObserver():
	core::pose::datacache::CacheableObserver(),
	sd_(),
	length_event_link_()
{
}

/// @brief default constructor
StructureDataObserver::StructureDataObserver( StructureDataOP sd_op ):
	core::pose::datacache::CacheableObserver(),
	sd_( sd_op ),
	length_event_link_()
{
}

StructureDataObserver::~StructureDataObserver()
{
	detach_from();
}

core::pose::datacache::CacheableObserverOP
StructureDataObserver::create()
{
	return core::pose::datacache::CacheableObserverOP( new StructureDataObserver() );
}

core::pose::datacache::CacheableObserverOP
StructureDataObserver::clone()
{
	return core::pose::datacache::CacheableObserverOP( new StructureDataObserver( *this ) );
}


bool
StructureDataObserver::is_attached() const
{
	//TR.Debug << "StructureDataObserver::is_attached(): " << length_event_link_.valid() << std::endl;
	return length_event_link_.valid();
}

void
StructureDataObserver::on_length_change( core::conformation::signals::LengthEvent const & event )
{
	if ( event.tag == core::conformation::signals::LengthEvent::RESIDUE_APPEND ) {
		on_residue_append( event );
	} else if ( event.tag == core::conformation::signals::LengthEvent::RESIDUE_PREPEND ) {
		on_residue_prepend( event );
	} else if ( event.tag == core::conformation::signals::LengthEvent::RESIDUE_DELETE ) {
		on_residue_delete( event );
	} else if ( event.tag == core::conformation::signals::LengthEvent::EMPTY ) {
		TR.Debug << "Empty signal" << std::endl;
	} else if ( event.tag == core::conformation::signals::LengthEvent::INVALIDATE ) {
		TR.Debug << "Invalidate signal" << std::endl;
	} else {
		TR.Debug << "Unknown signal" << std::endl;
	}
}

void
StructureDataObserver::on_residue_append( core::conformation::signals::LengthEvent const & ASSERT_ONLY(event) )
{
	debug_assert( event.tag == core::conformation::signals::LengthEvent::RESIDUE_APPEND );
}

void
StructureDataObserver::on_residue_prepend( core::conformation::signals::LengthEvent const & ASSERT_ONLY(event) )
{
	debug_assert( event.tag == core::conformation::signals::LengthEvent::RESIDUE_PREPEND );
}

void
StructureDataObserver::on_residue_delete( core::conformation::signals::LengthEvent const & event )
{
	debug_assert( event.tag == core::conformation::signals::LengthEvent::RESIDUE_DELETE );
	TR << "On_residue_delete(): " << event.position << " " << event.length_change
		<< " " << event.conformation_size << std::endl;
	debug_assert( event.length_change < 0 );
	core::Size const nres = -event.length_change;
	for ( core::Size resid=1; resid<=nres; ++resid ) {
		sd_->delete_residue( event.position );
	}
}

/*
/// @brief residue position where the event happened
Size position;

/// @brief The length of the conformation after the length event
/// @details This is needed for things like regenerating mappings
/// after multiple applications of events which may invalidate the confromation.
Size conformation_size;

/// @brief overall length change of the conformation
int length_change;

/// @brief direct access to residue
/// @remarks Almost always want to use this to access the residue instead of
///  the conformation.  Calling Conformation::residue() can cause an internal
///  update/re-sync inside Pose, which may have consequences if you're depending
///  upon multiple residue operations to be setup (such as bond angle/length
///  changes) prior to an internal update.
Residue const * residue;
*/

void
StructureDataObserver::attach_impl( core::pose::Pose & pose )
{
	length_event_link_ = pose.conformation().attach_length_obs( &StructureDataObserver::on_length_change, this );
	/* SET
	core::pose::datacache::LengthEventCollectorOP lencollect( new core::pose::datacache::LengthEventCollector() );
	pose.observer_cache().set( core::pose::datacache::CacheableObserverType::LENGTH_EVENT_COLLECTOR, lencollect );
	*/
	/* GET
	pose::datacache::CacheableObserverCOP len_obs = pose.observer_cache().get_const_ptr( pose::datacache::CacheableObserverType::LENGTH_EVENT_COLLECTOR );
	pose::datacache::LengthEventCollectorCOP lencollect( utility::pointer::static_pointer_cast< pose::datacache::LengthEventCollector const >( len_obs ) );
	*/
}

void
StructureDataObserver::detach_impl()
{
	length_event_link_.invalidate();
}

/*
StructureData const &
StructureDataObserver::sd() const
{
debug_assert( sd_ );
return *sd_;
}

StructureDataCOP
StructureDataObserver::sd_ptr() const
{
return sd_;
}
*/

} // namespace components
} // namespace denovo_design
} // namespace protocols

#ifdef    SERIALIZATION

/// @brief Serialization method for StructureDataObserver
/// @details This does not serialize either the length_event_ nor length_event_link_
/// data members as these hold raw pointers and would be invalid following deserialization.
template< class Archive >
void
protocols::denovo_design::components::StructureDataObserver::save( Archive & arc ) const
{
	arc( cereal::base_class< core::pose::datacache::CacheableObserver >( this ) );
	arc( CEREAL_NVP( sd_ ) );
	// EXEMPT length_event_link_
}

/// @brief Deserialization method
/// @brief This does not deserialize either the length_event_ vector nor the length_event_link_
/// data member because both of these classes hold raw pointers and would be invalid following
/// deserialization.
template< class Archive >
void
protocols::denovo_design::components::StructureDataObserver::load( Archive & arc )
{
	arc( cereal::base_class< core::pose::datacache::CacheableObserver >( this ) );
	StructureDataOP local_sd_ptr;
	arc( local_sd_ptr ); // StructureDataCOP
	sd_ = local_sd_ptr;
	// EXEMPT length_event_link_
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::denovo_design::components::StructureDataObserver );
CEREAL_REGISTER_TYPE( protocols::denovo_design::components::StructureDataObserver )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_denovo_design_components_StructureDataObserver )
#endif // SERIALIZATION
