// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/network/ui_mover.cc
/// @brief: UIMover, apply to send in-progress results to UI
///
/// @author Sergey Lyskov

#include <protocols/network/ui_mover.hh>

#include <protocols/network/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/ObserverCache.hh>

#include <utility/json_utilities.hh>

// #ifdef    SERIALIZATION
// #include <utility/serialization/serialization.hh>
// #include <cereal/types/polymorphic.hpp>
// #endif // SERIALIZATION


namespace protocols {
namespace network {

using std::string;

#if defined(ZEROMQ)

UIMover::UIMover() : hal(zmq_context(), ZMQ_DEALER)
{
	hal.connect(_hal_address_);
}

UIMover::UIMover(UIMover const &) : UIMover()
{
}

UIMover::~UIMover() = default;

UIMover & UIMover::operator= (UIMover const &)
{
	return *this;
}

void UIMover::apply(Pose & pose)
{
	apply( static_cast<Pose const &>(pose) );
}


void UIMover::apply(Pose const & pose)
{
	auto pose_binary = protocols::network::pose_to_bytes(pose);

	nlohmann::json result;

	result["pose"] = pose_binary;

	string binary_result;
	nlohmann::json::basic_json::to_msgpack(result, binary_result);

	send_message(hal, _m_progress_, ZMQ_SNDMORE);
	send_message(hal, binary_result);
}


#else // defined(ZEROMQ)

UIMover::UIMover() = default;

UIMover::UIMover(UIMover const &) = default;

UIMover::~UIMover() = default;

void UIMover::apply(Pose &) {}
void UIMover::apply(Pose const &) {}

UIMover & UIMover::operator= (UIMover const &) = default;

#endif // defined(ZEROMQ)




UIObserver::UIObserver(): CacheableObserver()
{
}

UIObserver::UIObserver(UIObserver const & rval) :
	CacheableObserver( rval ),
	type_( rval.type_ ),
	ui_(rval.ui_)
	// Do NOT copy the *_event_link_s
{
}

UIObserver::~UIObserver() {
	detach_from();
}

UIObserver & UIObserver::operator= (UIObserver const &other) {
	if ( this != &other ) {
		core::pose::datacache::CacheableObserver::operator=( other );
		type_ = other.type_;
		ui_ = other.ui_;
		// Do NOT copy the *_event_link_s
	}
	return *this;
}

core::pose::datacache::CacheableObserverOP
UIObserver::clone() {
	return core::pose::datacache::CacheableObserverOP( new UIObserver( *this ) );
}

core::pose::datacache::CacheableObserverOP
UIObserver::create() {
	return core::pose::datacache::CacheableObserverOP( new UIObserver );
}

void
UIObserver::set_type( ObserverType setting ) {
	type_ = setting;
	// We don't have a pose, so wait until we get attached
}

void
UIObserver::add_type( ObserverType setting ) {
	type_ = type_ | setting;
	// We don't have a pose, so wait until we get attached
}

void UIObserver::attach(core::pose::Pose &p)
{
	attach_to(p);
}

void UIObserver::detach(core::pose::Pose & /*p*/)
{
	detach_from();
}

bool
UIObserver::is_attached() const {
	return general_event_link_.valid() || energy_event_link_.valid() || conformation_event_link_.valid();
}

void
UIObserver::attach_impl(core::pose::Pose & pose) {
	general_event_link_.invalidate();
	energy_event_link_.invalidate();
	conformation_event_link_.invalidate();

	if ( has_observer_of_type(type_, ObserverType::general) ) general_event_link_ = pose.attach_general_obs( &UIObserver::generalEvent, this );

	if ( has_observer_of_type(type_, ObserverType::energy) ) energy_event_link_ = pose.attach_energy_obs( &UIObserver::energyEvent, this );

	if ( has_observer_of_type(type_, ObserverType::conformation) ) conformation_event_link_ = pose.attach_conformation_obs( &UIObserver::conformationEvent, this );
}

void UIObserver::detach_impl() {
	general_event_link_.invalidate();
	energy_event_link_.invalidate();
	conformation_event_link_.invalidate();
}


void UIObserver::generalEvent( core::pose::signals::GeneralEvent const & ev) {
	ui_.apply( *ev.pose );
}

void UIObserver::energyEvent( core::pose::signals::EnergyEvent const & ev) {
	ui_.apply( *ev.pose );
}

void UIObserver::conformationEvent( core::pose::signals::ConformationEvent const & ev) {
	ui_.apply( *ev.pose);
}


UIObserverOP get_ui_observer(core::pose::Pose & pose)
{
	using namespace core::pose::datacache;

	if ( !pose.observer_cache().has( CacheableObserverType::UI_OBSERVER) ) {
		UIObserverOP obs( new UIObserver );
		pose.observer_cache().set( CacheableObserverType::UI_OBSERVER, obs, /*autoattach*/ false );
	}
	CacheableObserverOP obs = pose.observer_cache().get_ptr( core::pose::datacache::CacheableObserverType::UI_OBSERVER);
	return utility::pointer::dynamic_pointer_cast< UIObserver >( obs );
}

UIObserverOP AddUIObserver(core::pose::Pose &p) //, bool keep_history, core::Real update_interval)
{
	UIObserverOP o( get_ui_observer(p) );
	//o->ui().keep_history(keep_history);
	//o->ui().update_interval(update_interval);
	o->add_type( UIObserver::ObserverType::general );
	o->attach(p);
	return o;
}

UIObserverOP AddUIObserver_to_energies(core::pose::Pose &p) //, bool keep_history, core::Real update_interval)
{
	UIObserverOP o( get_ui_observer(p) );
	//o->pymol().keep_history(keep_history);
	//o->pymol().update_interval(update_interval);
	o->add_type( UIObserver::ObserverType::energy );
	o->attach(p);
	return o;
}

UIObserverOP AddUIObserver_to_conformation(core::pose::Pose &p) //, bool keep_history, core::Real update_interval)
{
	UIObserverOP o( get_ui_observer(p) );
	//o->pymol().keep_history(keep_history);
	//o->pymol().update_interval(update_interval);
	o->add_type( UIObserver::ObserverType::conformation );
	o->attach(p);
	return o;
}


} // namespace network
} // namespace protocols



// #ifdef    SERIALIZATION

// /// @brief Automatically generated serialization method
// template< class Archive >
// void
// protocols::network::UIObserver::save( Archive & ) const
// {
//  // no-op
// }

// /// @Brief Automatically generated deserialization method
// template< class Archive >
// void
// protocols::network::UIObserver::load( Archive & )
// {
//  // no-op, and set observation type to `none` to avoid UI self recursion
//  type_ = ObserverType::none;
// }

// SAVE_AND_LOAD_SERIALIZABLE( protocols::network::UIObserver );
// CEREAL_REGISTER_TYPE( protocols::network::UIObserver )

// CEREAL_REGISTER_DYNAMIC_INIT( protocols_network_UIObserver )
// #endif // SERIALIZATION
