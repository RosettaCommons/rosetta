// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/network/ui_mover.hh
/// @brief: UIMover, apply to send in-progress results to UI
///
/// @author Sergey Lyskov

#pragma once

#include <protocols/network/ui_mover.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/pose/signals/GeneralEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/datacache/CacheableObserver.hh>
#include <protocols/moves/Mover.hh>

#include <utility/signals/Link.hh>

#ifdef ZEROMQ
#include <libzmq/include/zmq.hpp>
#endif // ZEROMQ

// #ifdef    SERIALIZATION
// #include <cereal/types/polymorphic.fwd.hpp>
// #endif // SERIALIZATION


namespace protocols {
namespace network {


/// Mover to send in-progress data to UI client
/// Note: you will need to build Rosetta with extra=zmq for this Mover to do any real work
class UIMover : public protocols::moves::Mover
{
public:
	/// @brief ctor
	explicit UIMover();

	/// @brief cctor
	explicit UIMover(UIMover const &);

	/// @brief dtor
	~UIMover() override;

	void apply(Pose &) override;
	virtual void apply(Pose const &);

	std::string get_name() const override { return "UIMover"; };

	UIMover & operator= (UIMover const &other);

#ifdef ZEROMQ

private:
	mutable zmq::socket_t hal;

#endif // ZEROMQ

};




/// @brief Special Observer which apply UIObserver if Pose is changed
class UIObserver : public core::pose::datacache::CacheableObserver
{
public:
	// This is set up to allow multiple settings with bit twiddling
	enum class ObserverType {
		general      = 1,
		energy       = 2,
		conformation = 4,

		none         = 0,
	};

	UIObserver();
	UIObserver(UIObserver const &);
	~UIObserver() override;

	UIObserver & operator= (UIObserver const &);


	core::pose::datacache::CacheableObserverOP clone() override;


	core::pose::datacache::CacheableObserverOP create() override;

	void set_type(ObserverType setting);

	void add_type(ObserverType setting);

	ObserverType get_type() const { return type_; };

	bool is_attached() const override;

	virtual void generalEvent( core::pose::signals::GeneralEvent const & ev);
	virtual void energyEvent( core::pose::signals::EnergyEvent const & ev);
	virtual void conformationEvent( core::pose::signals::ConformationEvent const & ev);

	UIMover & ui_mover() { return ui_; };

	/// Attach observer to the pose object
	void attach(core::pose::Pose &p);

	/// Detach observer from the pose object
	void detach(core::pose::Pose &p);

protected:

	void
	attach_impl(core::pose::Pose &pose) override;

	void
	detach_impl() override;

	void
	update_links();

private:

	ObserverType type_ = ObserverType::none;
	UIMover ui_;

	utility::signals::Link general_event_link_;
	utility::signals::Link energy_event_link_;
	utility::signals::Link conformation_event_link_;


	// #ifdef    SERIALIZATION
	// public:
	//  template< class Archive > void save( Archive & arc ) const;
	//  template< class Archive > void load( Archive & arc );
	// #endif // SERIALIZATION

};


inline UIObserver::ObserverType operator|(UIObserver::ObserverType lhs, UIObserver::ObserverType rhs)
{
	using I = std::underlying_type<UIObserver::ObserverType>::type;
	return static_cast<UIObserver::ObserverType>(static_cast<I>(lhs) | static_cast<I>(rhs));
}

inline UIObserver::ObserverType operator&(UIObserver::ObserverType lhs, UIObserver::ObserverType rhs)
{
	using I = std::underlying_type<UIObserver::ObserverType>::type;
	return static_cast<UIObserver::ObserverType>(static_cast<I>(lhs) & static_cast<I>(rhs));
}


inline bool has_observer_of_type(UIObserver::ObserverType flag, UIObserver::ObserverType observer)
{
	return (flag & observer) == observer;
}



/// @brief (Internal) helper function to create a UIObserver and add it to the given pose
/// NOTE: You NEED to adjust the observer type and call attach() on the return - by default a new UIObserver isn't attached/observing.
UIObserverOP get_ui_observer(core::pose::Pose & pose);

/// @brief Helper function that create UIObserver Object and add it to the give Pose.
///        This is the most likely the only function that you need to call...
UIObserverOP AddUIObserver(core::pose::Pose &p);//, bool keep_history=false, core::Real update_interval=0);

/// @brief Helper function that create UIObserver Object and add it to the give Pose energies object so pymol only updates on energy changes.
UIObserverOP AddUIObserver_to_energies(core::pose::Pose & p);//, bool keep_history=false, core::Real update_interval=0);

/// @brief Helper function that create UIObserver Object and add it to the give Pose conformation object so pymol only updates on conformation changes.
UIObserverOP AddUIObserver_to_conformation(core::pose::Pose & p);//, bool keep_history = false, core::Real update_interval = 0);


} // namespace network
} // namespace protocols

// #ifdef    SERIALIZATION
// CEREAL_FORCE_DYNAMIC_INIT( protocols_network_UIObserver )
// #endif // SERIALIZATION
