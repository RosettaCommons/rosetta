#include <boost/python.hpp>

#include <core/scoring/constraints/AtomPairConstraint.hh>

namespace bp = boost::python;



core::scoring::constraints::ConstraintCOP rtest_COP( core::scoring::constraints::ConstraintCOP cst )
{
	return cst;
}


void __constraints_by_hand_beginning__()
{
    bp::class_<core::scoring::constraints::Constraint,
               utility::pointer::owning_ptr< core::scoring::constraints::Constraint >, boost::noncopyable >
      ("Constraint", bp::no_init)
      /*
        .def("clone", bp::pure_virtual(&core::scoring::constraints::Constraint::clone))
		.def("natoms", bp::pure_virtual(&core::scoring::constraints::Constraint::natoms))
		.def("atom", bp::pure_virtual(&core::scoring::constraints::Constraint::atom))
		.def("score", bp::pure_virtual(&core::scoring::constraints::Constraint::score))
		.def("fill_f1_f2", bp::pure_virtual(&core::scoring::constraints::Constraint::fill_f1_f2))
*/
      ;

    bp::class_< utility::pointer::owning_ptr< core::scoring::constraints::Constraint const > >("ConstraintCOP");

    bp::implicitly_convertible< utility::pointer::owning_ptr< ::core::scoring::constraints::Constraint >
                              , utility::pointer::owning_ptr< ::core::scoring::constraints::Constraint const > >();

    //bp::implicitly_convertible< utility::pointer::owning_ptr< ::core::scoring::constraints::AtomPairConstraint >
    //                          , utility::pointer::owning_ptr< ::core::scoring::constraints::Constraint > >();



    //bp::implicitly_convertible< ::core::scoring::constraints::AtomPairConstraint
    //                          , ::core::scoring::constraints::Constraint > ();


    /*bp::implicitly_convertible< utility::pointer::owning_ptr< ::core::scoring::constraints::AtomPairConstraint const >
                              , utility::pointer::owning_ptr< ::core::scoring::constraints::Constraint const > >();


*/
    //bp::converter::implicit< ::core::scoring::constraints::AtomPairConstraint
    //                          , ::core::scoring::constraints::Constraint >();

    //bp::implicitly_convertible< ::core::scoring::constraints::Constraint const
      //                        , utility::pointer::owning_ptr< ::core::scoring::constraints::Constraint const > >();
}
