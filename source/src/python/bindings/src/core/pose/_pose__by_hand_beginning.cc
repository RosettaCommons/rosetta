// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
///
/// @author Sergey Lyskov
///

#include "boost/python.hpp"

//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/scoring/NeighborList.hh>
//#include <core/scoring/SecondaryStructureWeights.hh>
//#include <core/scoring/Energies.hh>


#include <core/pose/Pose.hh>
#include <core/pose/signals/GeneralEvent.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/io/pdb/pose_io.hh>
#include "core/import_pose/import_pose.hh"


//#include <core/io/pdb/pose_io.hh>
//#include <core/scoring/ScoreFunctionInfo.hh>

namespace bp = boost::python;


class PosePyObserver : public utility::pointer::ReferenceCount
{
public:
	PosePyObserver() {};
	virtual ~PosePyObserver() {
		std::cout << "~PosePyObserver..." << std::endl;
	};


	void add_observer(core::pose::Pose &p) {
		p.attach_general_obs(&PosePyObserver::generalEvent, this);
		//pose_ = p;
	}

	void remove_observer(core::pose::Pose &p) {
		p.detach_general_obs(&PosePyObserver::generalEvent, this);
		//pose_ = 0;
	}

	virtual void generalEvent( core::pose::signals::GeneralEvent const & ) {
		std::cout << "PosePyObserver::generalEvent... C++ version" << std::endl;
	};

private:
	core::pose::PoseOP pose_; // we want to keep OP to linked pose, so it cannot got deleted by Python
};
/*
struct Wrapper_PosePyObserver : public PosePyObserver
{
    Wrapper_PosePyObserver(PyObject* self_) : self(self_) {}

    void generalEvent(core::pose::signals::GeneralEvent const &ev)
    {
        return boost::python::call_method<void>(self, "generalEvent", ev);
    }

    PyObject* self;
};*/

struct Wrapper_PosePyObserver : public PosePyObserver, bp::wrapper<PosePyObserver>
{
    void generalEvent(core::pose::signals::GeneralEvent const &ev)
    {
    	bp::override f = this->get_override("generalEvent");
    	if( f ) f(ev);
    	else this->PosePyObserver::generalEvent(ev);
    }

    void default_generalEvent( core::pose::signals::GeneralEvent const & ev) {
    	this->PosePyObserver::generalEvent(ev);
    }
};


class PosePyObserver2 : public utility::pointer::ReferenceCount
{
public:
	PosePyObserver2() {};
	virtual ~PosePyObserver2() {
		std::cout << "~PosePyObserver..." << std::endl;
	};


	void link(core::pose::Pose &p) {
		p.attach_general_obs(&PosePyObserver2::generalEvent, this);
		//pose_ = p;
	}

	void unlink(core::pose::PoseOP &p) {
		p->detach_general_obs(&PosePyObserver2::generalEvent, this);
		//pose_ = 0;
	}

	virtual void generalEvent( core::pose::signals::GeneralEvent const & ) {
		std::cout << "PosePyObserver::generalEvent... C++ version" << std::endl;
	};

private:
	core::pose::PoseOP pose_; // we want to keep OP to linked pose, so it cannot got deleted by Python
};


// dummy function to test all events link
//void PosePyObserverTesterFunction(PosePyObserver &observer)
//void PosePyObserverTesterFunction(core::pose::Pose &__pose, core::scoring::ScoreFunction &__scorefxn)
void PosePyObserverTesterFunction()
{
	//core::pose::signals::GeneralEvent ge;
	//observer.generalEvent( ge );

	core::pose::Pose pose;
	core::import_pose::pose_from_pdb(pose, "test_in.pdb");


	std::cout << "Creating PosePyObserver...\n";
	PosePyObserver2 PO;

	std::cout << "Linkin... PosePyObserver...\n";
	PO.link(pose);

	core::scoring::ScoreFunction scorefxn;

	std::cout << "Scoring pose agagin...\n";
	scorefxn(pose);

	//PO.unlink(pose);

	std::cout << "Yeah... we still alive...\n";
}



void __pose_by_hand_beginning__()
{   /*
	bp::class_< PosePyObserver >("PosePyObserverBase")
		.def("link", &PosePyObserver::link)
		.def("unlink", &PosePyObserver::unlink)

		.def("generalEvent", &PosePyObserver::generalEvent)
    ; */

    //boost::python::class_<PosePyObserver, Wrapper_PosePyObserver, boost::noncopyable>( "PosePyObserver" )
    //boost::python::class_<Wrapper_PosePyObserver, boost::noncopyable>( "PosePyObserver" )
	boost::python::class_<Wrapper_PosePyObserver, utility::pointer::owning_ptr<Wrapper_PosePyObserver>, boost::noncopyable>( "PosePyObserver" )
 		.def("add_observer", &PosePyObserver::add_observer)
		.def("remove_observer", &PosePyObserver::remove_observer)

		.def("generalEvent", &PosePyObserver::generalEvent, &Wrapper_PosePyObserver::default_generalEvent)
    ;

    bp::def("QQQ_PosePyObserverTesterFunction", PosePyObserverTesterFunction);

}
