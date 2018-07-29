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

#include <utility/thread/shared_thread_local_data.impl.hh>

#include <utility/json_utilities.hh>

// #ifdef    SERIALIZATION
// #include <utility/serialization/serialization.hh>
// #include <cereal/types/polymorphic.hpp>
// #endif // SERIALIZATION


namespace protocols {
namespace network {

using std::string;


#if defined(ZEROMQ)

auto const _hal_socket_high_water_mark_  = 4;


UIMover::HalSocket::HalSocket() : socket(zmq_context(), ZMQ_DEALER)
{
	//std::cout << "UIMover::HalSocket()" << std::endl;

	socket.setsockopt(ZMQ_LINGER, 100);  // 100 milliseconds

	socket.connect(_hal_address_);

	socket.setsockopt(ZMQ_SNDHWM, _hal_socket_high_water_mark_);
	socket.setsockopt(ZMQ_RCVHWM, _hal_socket_high_water_mark_);
}



UIMover::UIMover() = default;

UIMover::UIMover(UIMover const &) = default;

UIMover::~UIMover() = default;


// UIMover::UIMover(UIMover const &other) : hal_(other.hal_)
// {
//  std::cout << "UIMover::UIMover(UIMover const &)" << std::endl;
// }

// UIMover::~UIMover()
// {
//  std::cout << "UIMover::~UIMover()" << std::endl;
// }


UIMover & UIMover::operator= (UIMover const &)
{
	//std::cout << "IMover & UIMover::operator= (UIMover const &)" << std::endl;
	return *this;
}

void UIMover::apply(Pose & pose)
{
	apply( static_cast<Pose const &>(pose) );
}


void UIMover::apply(Pose const & pose)
{
	sleep_if_paused();

	zmq::socket_t & hal = hal_.get().socket;

	auto pose_binary = protocols::network::pose_to_bytes(pose);

	nlohmann::json result;

	result[_f_pose_] = pose_binary;

	string binary_result;
	nlohmann::json::basic_json::to_msgpack(result, binary_result);

	send_message(hal, _m_progress_, ZMQ_SNDMORE);
	send_message(hal, binary_result);

	sleep_if_paused();
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
	//std::cout << "UIObserver::UIObserver()" << std::endl;
}

UIObserver::UIObserver(UIObserver const & rval) :
	CacheableObserver( rval ),
	type_( rval.type_ ),
	ui_(rval.ui_)
	// Do NOT copy the *_event_link_s
{
	//std::cout << "UIObserver::UIObserver(UIObserver const &)" << std::endl;
}

UIObserver::~UIObserver()
{
	//std::cout << "UIObserver::~UIObserver()" << std::endl;
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
	//std::cout << "UIObserverOP AddUIObserver(core::pose::Pose &)" << std::endl;

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
