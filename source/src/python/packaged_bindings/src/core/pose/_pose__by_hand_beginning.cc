// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
///
/// @author Sergey Lyskov


#include "boost/python.hpp"

#include <core/pose/Pose.hh>
#include <core/pose/signals/GeneralEvent.hh>

namespace bp = boost::python;

class PosePyObserver : public utility::pointer::ReferenceCount
{
public:
	PosePyObserver() {};
	virtual ~PosePyObserver() {};


	void add_observer(core::pose::Pose &p) {
		p.attach_general_obs(&PosePyObserver::generalEvent, this);
	}

	void remove_observer(core::pose::Pose &p) {
		p.detach_general_obs(&PosePyObserver::generalEvent, this);
	}

	virtual void generalEvent( core::pose::signals::GeneralEvent const & ) { };

private:
	core::pose::PoseOP pose_; // we want to keep OP to linked pose, so it cannot got deleted by Python
};

struct Wrapper_PosePyObserver : public PosePyObserver, bp::wrapper<PosePyObserver>
{
    void generalEvent(core::pose::signals::GeneralEvent const &ev)
    {
    	bp::override f = this->get_override("generalEvent");

    	if( f ) {
        f(ev);
      } else {
        this->PosePyObserver::generalEvent(ev);
      }
    }

    void default_generalEvent( core::pose::signals::GeneralEvent const & ev) {
    	this->PosePyObserver::generalEvent(ev);
    }
};

void __pose_by_hand_beginning__()
{
	boost::python::class_<Wrapper_PosePyObserver, utility::pointer::shared_ptr<Wrapper_PosePyObserver>, boost::noncopyable>( "PosePyObserver" )
 		.def("add_observer", &PosePyObserver::add_observer)
		.def("remove_observer", &PosePyObserver::remove_observer)
		.def("generalEvent", &PosePyObserver::generalEvent, &Wrapper_PosePyObserver::default_generalEvent)
    ;
}
