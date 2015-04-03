// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#include "boost/python.hpp"


#include <protocols/moves/Mover.hh>

//#include <protocols/moves/TrialMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
//#include <protocols/moves/MonteCarlo.hh>
#include <core/pose/Pose.fwd.hh>

#include <iostream>

//#include <boost/python/suite/indexing/map_indexing_suite.hpp>

namespace bp = boost::python;

/*
class PyMover : public protocols::moves::Mover
{
public:
	virtual std::string get_name() const { return "PyMover"; };

	virtual void apply( Pose & ) {
	    	//std::cout << "PyMover:apply..." << std::endl;
	};
};

struct Wrapper_PyMover : public PyMover, bp::wrapper<PyMover>
{

	void apply( Pose & pose )
    {
    	//std::cout << "Wrapper_Mover:apply..." << std::endl;
    	bp::override f = this->get_override("apply");
    	if( f ) f(pose);
    	else this->PyMover::apply(pose);
    }

    void default_apply( Pose & pose ) {
	    //std::cout << "Wrapper PyMover default_apply..." << std::endl;
    	this->PyMover::apply(pose);
    }
};
*/

#include <protocols/simple_moves/MinMover.hh>

protocols::moves::MoverOP QQQ_CreateMinMover()
{
	return new protocols::simple_moves::MinMover();
}

void QQQ_SubclassTester(protocols::moves::MoverOP m)
{
	std::cout << "QQQ_SubclassTester..." << m->get_name() << std::endl;
}

void QQQ_SubclassTester2(protocols::moves::MoverOP m, core::pose::PoseOP pose)
{
	std::cout << "QQQ_SubclassTester...2" << m->get_name() << std::endl;
	m->apply(*pose);
}

void __moves_by_hand_beginning__()
{
	// Debug
	bp::def("QQQ_CreateMinMover", QQQ_CreateMinMover);
	bp::def("QQQ_SubclassTester", QQQ_SubclassTester);
	bp::def("QQQ_SubclassTester2", QQQ_SubclassTester2);

    // bp::implicitly_convertible< utility::pointer::owning_ptr< ::protocols::moves::Mover >
    //                           , utility::pointer::owning_ptr< ::protocols::moves::Mover const > >();

	/*bp::implicitly_convertible< utility::pointer::owning_ptr<  PyMover >
                              , utility::pointer::owning_ptr< ::protocols::moves::Mover > >();


	bp::implicitly_convertible< utility::pointer::owning_ptr<  Wrapper_PyMover >
                              , utility::pointer::owning_ptr< ::protocols::moves::Mover > >();

	bp::implicitly_convertible< utility::pointer::owning_ptr< Wrapper_PyMover >
                              , utility::pointer::owning_ptr< Wrapper_PyMover const > >();


    //boost::python::class_<Wrapper_Mover, bp::bases< ::protocols::moves::Mover >, utility::pointer::owning_ptr<PyMover>, boost::noncopyable>( "PyMover" )
    boost::python::class_<Wrapper_PyMover, utility::pointer::owning_ptr<Wrapper_PyMover>, boost::noncopyable>( "PyMover" )
		.def("apply", &PyMover::apply, &Wrapper_PyMover::default_apply)
    ;
	*/

  /*  bp::def("toMover", toMover);
    bp::def("toMC", toMC);

    bp::implicitly_convertible< utility::pointer::owning_ptr< ::protocols::moves::TrialMover >
                              , utility::pointer::owning_ptr< ::protocols::moves::Mover > >();

    bp::implicitly_convertible< utility::pointer::owning_ptr< ::protocols::simple_moves::SmallMover >
                              , utility::pointer::owning_ptr< ::protocols::moves::Mover > >();

    bp::implicitly_convertible< utility::pointer::owning_ptr< ::protocols::simple_moves::ShearMover >
                              , utility::pointer::owning_ptr< ::protocols::moves::Mover > >();
*/

	/*
	bp::class_< ::protocols::moves::PyMolMover::ColorMap >("ColorMap")
		.def( bp::map_indexing_suite< ::protocols::moves::PyMolMover::ColorMap >())
		; */
}
