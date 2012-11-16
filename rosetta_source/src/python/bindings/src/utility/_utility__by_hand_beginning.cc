// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include "boost/python.hpp"
// #include "boost/python/suite/indexing/indexing_suite.hpp"
// #include "boost/python/suite/indexing/vector_indexing_suite.hpp"
// #include "boost/python/suite/indexing/map_indexing_suite.hpp"

#include "utility/vector1.hh"
#include "utility/pointer/access_ptr.hh"

#include "core/chemical/AtomType.hh"
#include "core/chemical/AtomTypeSet.hh"
#include "core/chemical/MMAtomType.hh"
#include "core/chemical/MMAtomTypeSet.hh"
#include "core/chemical/ResidueType.hh"
#include "core/chemical/ResidueTypeSet.hh"

#include "core/conformation/Atom.hh"

#include "core/coarse/Translator.hh"
#include "core/coarse/CoarseEtable.hh"

// for AP wrapping
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>

#include <core/pack/task/PackerTask.hh>

#include <core/scoring/ScoreType.hh>

#include <utility/stream_util.hh>


#include "utility/exit.hh"

#include <ostream>
#include <set>


namespace bp = boost::python;
using namespace std;
using namespace utility;


#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>

// Some subclassing testing functions
void Q_Test_CI1B(core::scoring::methods::ContextIndependentOneBodyEnergyOP )
{
    std::cout << "Q_Test_CI1B!" << std::endl;
}

// Some subclassing testing functions
void Q_Test_EnergyMethodCreator(core::scoring::methods::EnergyMethodCreatorOP cr)
{
    std::cout << "Q_Test_EnergyMethodCreator..." << std::endl;

    core::scoring::methods::EnergyMethodOptions options;
    cr->create_energy_method( options );
    std::cout << "Q_Test_EnergyMethodCreator... Done!" << std::endl;
}




template< class T >
T * getCAP( pointer::access_ptr<T> rs ) {
  T & rs_ref( *rs );
  T * rs_ptr = &rs_ref;
  return rs_ptr;
}

/*
template <class T>
std::ostream& operator <<(std::ostream &os, utility::vector1<T> const & v)
{
    os << "[";
    for(unsigned int i=1; i<=v.size(); i++) {
        os << v[i] << ", ";
    }
    os << "]";
    return os;
} */

template <class T>
std::string vector1_repr(utility::vector1<T> const & v)
{
    std::ostringstream os;

    os << "[";
    for(unsigned int i=1; i<=v.size(); i++) {
        os << v[i] << ", ";
    }
    os << "]";
    return os.str();
}

template< class TT > inline void vector1_set( vector1<TT> & v, size_t const & i, TT const & val ) { v[i] = val; }
template< class TT > inline std::size_t vector1_len( vector1<TT> & v ) { return v.size(); }

template< class TT > inline std::string vector1_str( vector1<TT> & v ) { std::ostringstream s; s<<v; return s.str(); }

template< class TT > inline typename vector1<TT>::iterator vector1_begin( vector1<TT> & v ) { return v.begin(); }
template< class TT > inline typename vector1<TT>::iterator vector1_end  ( vector1<TT> & v ) { return v.end(); }

template< class TT > inline void vector1_reserve( vector1<TT> & v, std::size_t n) { v.reserve(n); }
template< class TT > inline void vector1_resize( vector1<TT> & v, std::size_t n) { v.resize(n); }

template< class Htype, class CP, class CP_const>
void wrap_vector1(char * name) {
  typedef vector1<Htype> Ttype;
  typedef vectorL<1,Htype,allocator<Htype> > Btype;
  typedef vector<Htype> Vtype;
  bp::class_<Ttype>(name)
    .def( bp::init< size_t >() )
    .def( bp::init< vector1<Htype> const & >() )
    // .def( bp::init< size_t, TT >() )
    .def("__getitem__"
        , (Htype const & (Ttype::*)(size_t const) const)( &Ttype::at )
        , CP_const()    )
    .def("__getitem__"
        , (Htype & (Ttype::*)(size_t const))( &Ttype::at )
        , CP()        )
    .def("__setitem__"
        , &vector1_set<Htype> )
    .def("append"
        , (Btype & (Btype::*)(Htype const &))( &Btype::add_back )
        , bp::return_value_policy< bp::reference_existing_object >()        )
    .def("__len__", & vector1_len<Htype> )
    .def("__iter__", bp::range(&vector1_begin<Htype>,&vector1_end<Htype>) )

    //.def("__str__", & vector1_str<Htype> )
    .def("__str__", & vector1_repr<Htype> )
    //.def( bp::self_ns::str( bp::self ) )

    .def("reserve", &vector1_reserve<Htype> )
    .def("resize", &vector1_resize<Htype> )

  ;
}

template< class Htype, class CP, class CP_const>
void wrap_vector1_part(char * name) {
  typedef vector1<Htype> Ttype;
  typedef vectorL<1,Htype,allocator<Htype> > Btype;
  typedef vector<Htype> Vtype;
  bp::class_<Ttype>(name)
    .def("__getitem__"
        , (Htype const & (Ttype::*)(size_t const) const)( &Ttype::at )
        , CP_const()    )
    .def("__getitem__"
        , (Htype & (Ttype::*)(size_t const))( &Ttype::at )
        , CP()        )
    .def("__setitem__"
        , &vector1_set<Htype> )
    .def("append"
        , (Btype & (Btype::*)(Htype const &))( &Btype::add_back )
        , bp::return_value_policy< bp::reference_existing_object >()        )
    .def("__len__", & vector1_len<Htype> )
    .def("__iter__", bp::range(&vector1_begin<Htype>,&vector1_end<Htype>) )
  ;
}

// std::map --------------------------------------------------------------------------------------------------------------------
template< class T >
void wrap_std_map(char * name)
{
	bp::class_< T >(name)
		.def( bp::map_indexing_suite< T >())
		;
}


// std::set --------------------------------------------------------------------------------------------------------------------
template< class T > void add_to_set(std::set<T> & s, const T & v) { s.insert(v); }
template< class T > void erase_from_set(std::set<T> & s, const T & v) { s.erase(v); }

template< class TT > inline typename std::set<TT>::iterator set_begin( std::set<TT> & v ) { return v.begin(); }
template< class TT > inline typename std::set<TT>::iterator set_end  ( std::set<TT> & v ) { return v.end(); }

template <class T>
std::string set_repr(std::set<T> const & s)
{
	typedef std::set<T> Stype;
	typedef typename std::set<T>::iterator Stype_iterator;

    std::ostringstream os;
    os << "<set>[";
    for(Stype_iterator p=s.begin(); p!=s.end(); ++p) os << *p << ", ";
    os << "]";
    return os.str();
}


template< class Htype, class CP, class CP_const>
void wrap_std_set(char * name)
{
	typedef std::set<Htype> Ttype;
	bp::class_<Ttype>(name)
	.def( bp::init< >() )
	.def( bp::init< std::set<Htype> const & >() )

	.def("__contains__", &std::set<Htype>::count )
	.def("add", &add_to_set<Htype> )
	.def("erase",  &erase_from_set<Htype> )
	.def("__len__",  &std::set<Htype>::size )
	.def("__iter__", bp::range(&set_begin<Htype>, &set_end<Htype> ) )
	.def("__str__", &set_repr<Htype> )
	;
}

template< class Type >
void wrap_owning_pointer(char * name)
{
    bp::class_<Type>(name)
        //.def("get", &Type::get)
        .def("reset_to_null", &Type::reset_to_null)
    ;
}



#ifndef _MSC_VER
	template< class T >  T * wrap_access_pointer_get_function( pointer::access_ptr<T> rs ) {  return rs.get(); }
#else
	template< class T >  T * wrap_access_pointer_get_function( pointer::access_ptr<T> const & rs ) {  return rs.get(); }
#endif

/*
template< class T >
void wrap_access_pointer(std::string class_name)
{
    boost::python::implicitly_convertible< utility::pointer::access_ptr< T >
                                         , utility::pointer::access_ptr< T const > >();

    bp::class_< utility::pointer::access_ptr< T > >( std::string(class_name+"AP").c_str() )
        .def("get", (  T * (*)( utility::pointer::access_ptr<T> )  )( & wrap_access_pointer_get_function<T> )
             , bp::return_value_policy< bp::reference_existing_object >() );

    bp::class_< utility::pointer::access_ptr< T const > >( std::string(class_name+"CAP").c_str() )
        .def("get", (  T const * (*)( utility::pointer::access_ptr<T const > )  )( & wrap_access_pointer_get_function<T const> )
             , bp::return_value_policy< bp::reference_existing_object >() );
}
*/

// .def("__iter__", bp::range( &core::pose::Pose::res_begin, &core::pose::Pose::res_end));

inline bool vector1_bool_get ( vector1<bool> & v, int i ) { if(v[i]) return true; else return false; }
inline void vector1_bool_push( vector1<bool> & v, bool h ) { return v.push_back(h); }

void pyexit_callback(void)
{
    throw "RosettaException";
}

void set_pyexit_callback(void)
{
    set_main_exit_callback(pyexit_callback);
}


// Python Derived class demo -----------------------------------------------------------------------
// An abstract base class...
class DemoBase //: public boost::noncopyable
{
public:

    DemoBase() {};
    ~DemoBase() {};

    virtual std::string testMethod() { return "testMethod() of a C++ Base"; };
};

// A derived class...
class DemoDerived : public DemoBase
{
public:

    DemoDerived() {}
    ~DemoDerived() {}

    std::string testMethod()
    {
        return "testMethod() of a C++ Derived object called!";
    }
};

std::string DemoTesterFunction(DemoBase* p) { return p->testMethod(); }

// Boost.Python wrapper class for Base
struct BaseWrap : public DemoBase
{
    BaseWrap(PyObject* self_) : self(self_) {}

    std::string testMethod()
    {
        return boost::python::call_method<std::string>(self, "testMethod");
    }

    PyObject* self;
};
/* Python part of the demo
class PythonDerived( Base ):
    def testMethod( self ):
        return "testMethod() of a PythonDerived object called!"

p = PythonDerived()
DemoTesterFunction(p)
*/

class OOO
{
public:
};

void __utility_by_hand_beginning__()
{
    // Testing functions
    bp::def("Q_Test_CI1B", Q_Test_CI1B);
    bp::def("Q_Test_EnergyMethodCreator", Q_Test_EnergyMethodCreator);

    bp::def("set_pyexit_callback", set_pyexit_callback);
    //bp::def("set_main_exit_callback", set_main_exit_callback);
    bp::def("pyexit_callback", pyexit_callback);

    // Wraping for derived Demo
    bp::class_< DemoBase >("DemoBase")
    .def("testMethod", &DemoBase::testMethod)
    ;

    bp::class_< DemoDerived, bp::bases<DemoBase> >("DemoDerived")
    .def("testMethod", &DemoDerived::testMethod)
    ;

    boost::python::class_<DemoBase, BaseWrap, boost::noncopyable>( "Base" );

    bp::def("DemoTesterFunction", DemoTesterFunction);

    //wrap_owning_pointer<core::pack::task::PackerTaskOP>("PackerTaskOP");


    // OStringStream like wrappers ---------------------------------------------------------------------
    //bp::implicitly_convertible< utility::OStringStream & , std::ostringstream & >();
    //bp::implicitly_convertible< utility::OStringStream2 * , std::ostream * >();
    //bp::implicitly_convertible< utility::OStringStream2 *,  std::ostringstream * >();

    bp::class_< std::ostream, boost::noncopyable >("OStream", bp::no_init);

    typedef void ( std::ostringstream::*str1_function_type )( std::string const & );
    typedef string ( std::ostringstream::*str2_function_type )( ) const;

    bp::class_< std::ostringstream, bp::bases<std::ostream>, boost::noncopyable >("OStringStream")
        .def("str", str1_function_type( &::std::ostringstream::str ) )
        .def("str", str2_function_type( &::std::ostringstream::str ) )

      //.def("str", (void ( ::std::ostringstream::str::* )( std::string const & ) )( &::std::ostringstream::str )

    ;

    //, (void ( ::core::scoring::constraints::AtomPairConstraint::* )( ::core::id::AtomID const &,::core::conformation::Conformation const &,::core::Vector &,::core::Vector &,::core::scoring::EnergyMap const & ) const)( &::core::scoring::constraints::AtomPairConstraint::fill_f1_f2 )


    //     , (::core::scoring::constraints::ConstraintOP ( ::core::scoring::constraints::AtomPairConstraint::* )(  ) const)( &::core::scoring::constraints::AtomPairConstraint::clone )


    using namespace pointer;
    typedef bp::return_value_policy< bp::reference_existing_object > CP_REF;
    typedef bp::return_value_policy< bp::copy_const_reference >      CP_CCR;
    typedef bp::return_value_policy< bp::copy_non_const_reference >  CP_CNCR;



    // bp::class_< vector1<vector1<size_t> > >("utility___vec1_vec1_size")
    //   .def(bp::vector_indexing_suite< vector1<vector1<size_t> > >() );

    bp::class_< access_ptr< core::chemical::AtomTypeSet const   > >("core___chemical___AtomTypeSetCAP");
    bp::class_< access_ptr< core::chemical::ResidueType const   > >("core___chemical___ResidueTypeCAP");
    bp::class_< access_ptr< core::chemical::ResidueTypeSet const> >("core___chemical___ResidueTypeSetCAP");
    bp::class_< access_ptr< core::chemical::MMAtomTypeSet const > >("core___chemical___MMAtomTypeSetCAP");
    bp::class_< access_ptr< core::coarse::Translator const      > >("core___coarse___TranslatorCAP");
    bp::class_< access_ptr< core::coarse::CoarseEtable const    > >("core___coarse___CoarseEtableCAP");

    using namespace core::chemical;
    using namespace core::coarse;

    // old code - only for compatibility with previous verisons - deprecated, will be removed in the future...
    bp::def("utility___getCAP"
         , (  ResidueTypeSet const * (*)( access_ptr<ResidueTypeSet const> )  )( & getCAP<ResidueTypeSet const> )
         , bp::return_value_policy< bp::reference_existing_object >() );
    bp::def("utility___getCAP"
         , (  AtomTypeSet const * (*)( access_ptr<AtomTypeSet const> )  )( & getCAP<AtomTypeSet const> )
         , bp::return_value_policy< bp::reference_existing_object >() );
    bp::def("utility___getCAP"
         , (  ResidueType const * (*)( access_ptr<ResidueType const> )  )( & getCAP<ResidueType const> )
         , bp::return_value_policy< bp::reference_existing_object >() );
    bp::def("utility___getCAP"
         , (  MMAtomTypeSet const * (*)( access_ptr<MMAtomTypeSet const> )  )( & getCAP<MMAtomTypeSet const> )
         , bp::return_value_policy< bp::reference_existing_object >() );
    bp::def("utility___getCAP"
         , (  Translator const * (*)( access_ptr<Translator const> )  )( & getCAP<Translator const> )
         , bp::return_value_policy< bp::reference_existing_object >() );
    bp::def("utility___getCAP"
         , (  CoarseEtable const * (*)( access_ptr<CoarseEtable const> )  )( & getCAP<CoarseEtable const> )
         , bp::return_value_policy< bp::reference_existing_object >() );
    /*
    //wrap_access_pointer<  >("AP");
    wrap_access_pointer< core::chemical::AtomTypeSet >("core_chemical_AtomTypeSet");
    wrap_access_pointer< core::chemical::ResidueType >("core_chemical_ResidueType");
    wrap_access_pointer< core::chemical::ResidueTypeSet >("core_chemical_ResidueTypeSet");
    wrap_access_pointer< core::chemical::MMAtomTypeSet >("core_chemical_MMAtomTypeSet");
    wrap_access_pointer< core::coarse::Translator >("core_coarse_Translator");
    wrap_access_pointer< core::coarse::CoarseEtable >("core_coarse_CoarseEtable");

    // Adding AP for subclassin in Python
    wrap_access_pointer< core::pose::Pose >("core_pose_Pose");  // for PyMover

    wrap_access_pointer< core::conformation::Residue >("core_conformation_Residue");  // Energy methods
    wrap_access_pointer< core::scoring::EnergyMap >("core_scoring_EnergyMap");
    wrap_access_pointer< core::scoring::ScoreFunction >("core_scoring_ScoreFunction");
    wrap_access_pointer< core::id::DOF_ID >("core_id_DOF_ID_");
    wrap_access_pointer< core::id::TorsionID >("core_id_TorsionID_");
    */
    wrap_access_pointer< utility::vector1< bool > >("utility_vector1_bool_");


    boost::python::class_<OOO> my_obj("OOO");
    //boost::python::object my_obj;

    boost::python::scope within(my_obj);
    bp::def("Q_Test_CI1B", Q_Test_CI1B);
    bp::def("Q_Test_EnergyMethodCreator", Q_Test_EnergyMethodCreator);

}
