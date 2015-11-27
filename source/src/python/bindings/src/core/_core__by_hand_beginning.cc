// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include "utility/vector1.hh"
#include "utility/pointer/access_ptr.hh"

#include "core/chemical/AtomType.hh"
#include "core/chemical/AtomTypeSet.hh"
#include "core/chemical/MMAtomType.hh"
#include "core/chemical/MMAtomTypeSet.hh"
#include "core/chemical/ResidueType.hh"
#include "core/chemical/ResidueTypeSet.hh"

#include "core/conformation/Atom.hh"

//#include "core/coarse/Translator.hh"
//#include "core/coarse/CoarseEtable.hh"

// for AP wrapping
#include <core/kinematics/Edge.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>

#include <core/conformation/Residue.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>

#include <core/pack/task/PackerTask.hh>

#include <core/scoring/ScoreType.hh>

#include <utility/stream_util.hh>
#include "utility/exit.hh"

#include <platform/types.hh>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>


// #include <iostream>
// #include <ostream>
#include <istream>
#include <sstream>
#include <set>
#include <map>

#ifdef DEBUG
#include <cstdio>
#endif


#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>


// Includes for dummy bindings to simplify import orders
// #include <core/chemical/ChemicalManager.hh>
// #include <core/scoring/ScoringManager.hh>
// #include <core/scoring/constraints/ConstraintIO.hh>
// #include <core/scoring/constraints/ConstraintFactory.hh>
// #include <core/io/silent/SilentStructFactory.hh>
// #include <core/pack/task/operation/ResLvlTaskOperationFactory.hh>


namespace bp = boost::python;

#ifdef PYROSETTA_NUMPY

// Numpy API has strange import behavior
// Must define PY_ARRAY_UNIQUE_SYMBOL for all cpp files within a module which include numpy/arrayobject.h
// Must then call import_array to load the numpy module
// import_array1 returns given value if the module fails to load, may be used to raise importerrors, or failfast as in this case
#define PY_ARRAY_UNIQUE_SYMBOL PYROSETTA_ARRAY_API
#include <numpy/arrayobject.h>

namespace
{
  static struct pyrosetta_numpy_importer
  {
    static bool do_import_array()
    {
      // import_array1 is a macro defined in numpy/arrayobject.h which returns with its argument if import fails.
#ifdef DEBUG
      std::cerr << "Calling numpy import_array." << std::endl;
#endif
      import_array1(false);
      // Will not return if import fails
      return true;
    }

    pyrosetta_numpy_importer()
    {
      if (!do_import_array())
      {
        throw std::runtime_error("numpy failed to initialize.");
      }
    }
  } _array_importer;
}

// Arrayobject typedefs the lowercase npy_* to c types
// and defines the NPY_TYPES enum w/ corrosponding uppercase numpy
// type ids.
inline NPY_TYPES get_typenum(bool) { return NPY_BOOL; }
// inline NPY_TYPES get_typenum(npy_bool) { return NPY_BOOL; }
inline NPY_TYPES get_typenum(npy_byte) { return NPY_BYTE; }
inline NPY_TYPES get_typenum(npy_ubyte) { return NPY_UBYTE; }
inline NPY_TYPES get_typenum(npy_short) { return NPY_SHORT; }
inline NPY_TYPES get_typenum(npy_ushort) { return NPY_USHORT; }
inline NPY_TYPES get_typenum(npy_int) { return NPY_INT; }
inline NPY_TYPES get_typenum(npy_uint) { return NPY_UINT; }
inline NPY_TYPES get_typenum(npy_long) { return NPY_LONG; }
inline NPY_TYPES get_typenum(npy_ulong) { return NPY_ULONG; }
inline NPY_TYPES get_typenum(npy_longlong) { return NPY_LONGLONG; }
inline NPY_TYPES get_typenum(npy_ulonglong) { return NPY_ULONGLONG; }
inline NPY_TYPES get_typenum(npy_float) { return NPY_FLOAT; }
inline NPY_TYPES get_typenum(npy_double) { return NPY_DOUBLE; }
inline NPY_TYPES get_typenum(npy_cfloat) { return NPY_CFLOAT; }
inline NPY_TYPES get_typenum(npy_cdouble) { return NPY_CDOUBLE; }
inline NPY_TYPES get_typenum(std::complex<float>) { return NPY_CFLOAT; }
inline NPY_TYPES get_typenum(std::complex<double>) { return NPY_CDOUBLE; }
#if HAVE_LONG_DOUBLE && (NPY_SIZEOF_LONGDOUBLE > NPY_SIZEOF_DOUBLE)
inline NPY_TYPES get_typenum(npy_longdouble) { return NPY_LONGDOUBLE; }
inline NPY_TYPES get_typenum(npy_clongdouble) { return NPY_CLONGDOUBLE; }
inline NPY_TYPES get_typenum(std::complex<long double>) { return NPY_CLONGDOUBLE; }
#endif


// array scalars ------------------------------------------------------------
template <class T>
struct array_scalar_converter
{
  static const PyTypeObject * get_array_scalar_typeobj()
  {
    return (PyTypeObject *) PyArray_TypeObjectFromType(get_typenum(T()));
  }

  static PyArray_Descr * get_descr()
  {
    return PyArray_DescrFromType(get_typenum(T()));
  }

  static void * check_array_scalar(PyObject *obj)
  {
    if (obj->ob_type == get_array_scalar_typeobj())
    {
      // Check if the object type is the array scalar type corrosponding to the
      // input C type.
      return obj;
    }
    else if(PyArray_CheckScalar(obj) && !PyArray_IsZeroDim(obj))
    {
      // Convert to array scalar to check casting to handle
      // signed/unsigned conversion.
      if(check_zero_dim_array(PyArray_FromScalar(obj, NULL)) != 0)
      {
        return obj;
      }
      else
      {
        return 0;
      }
    }
    else
    {
      return 0;
    }
  }

  static void convert_array_scalar(
    PyObject* obj,
    bp::converter::rvalue_from_python_stage1_data* data)
  {
#ifdef DEBUG
      std::cerr << "array_scalar_converter converting. typenum: " << get_typenum(T()) << "\n";
#endif

    void* storage = ((bp::converter::rvalue_from_python_storage<T>*)data)->storage.bytes;

    // no constructor needed, only dealing with POD types
    PyArray_CastScalarToCtype(obj, storage, get_descr());

    // record successful construction
    data->convertible = storage;
  }

  static void * check_zero_dim_array(PyObject *obj)
  {
    // Check if the object is an PyArray_Type of
    // dimensionality 0
    // See:
    // http://docs.scipy.org/doc/numpy/reference/arrays.scalars.html
    // and
    // http://docs.scipy.org/doc/numpy/reference/c-api.array.html#general-check-of-python-type
    //
    // Convert to an array scalar to perform C-type conversion check.
    //
    // Check if object is an array and a scalar
    // if so, the object must be a zero-length array.
    if(PyArray_IsZeroDim(obj) &&
       PyArray_CanCastArrayTo((PyArrayObject*) obj, get_descr(), NPY_SAME_KIND_CASTING))
    {
      return obj;
    }
    else
    {
      return 0;
    }
  }

  static void convert_zero_dim_array(
    PyObject* obj,
    bp::converter::rvalue_from_python_stage1_data* data)
  {
    convert_array_scalar(PyArray_ToScalar(PyArray_DATA(obj), obj), data);
  }
};

#endif


//template< class T >
// T * getCAP( utility::pointer::access_ptr<T> rs ) {
//   T & rs_ref( *rs );
//   T * rs_ptr = &rs_ref;
//   return rs_ptr;
// }

// std::pair ---------------------------------------------------------------------------------------------------


//vector1 wrapper requires ostream '<<' operator

template< class T1, class T2 >
std::ostream& operator<<(std::ostream& strm, const std::pair< T1, T2>& kvPair)
{
  strm << "(" << kvPair.first << ", " << kvPair.second << ")";
  return strm;
}

template <class T1, class T2>
std::string pair_repr(std::pair<T1, T2> const & v)
{
    std::ostringstream os;

    os << v;

    return os.str();
}

template< class T1, class T2 >
void wrap_std_pair(std::string name)
{
    bp::class_< std::pair< T1, T2 > >(name.c_str())
			.def( bp::init< T1 const &, T2 const & >(( bp::arg("__a"),
bp::arg("__b") )))
			.def_readwrite( "first", &std::pair< T1, T2 >::first)
			.def_readwrite( "second", &std::pair< T1, T2 >::second)
      .def("__str__", &pair_repr<T1, T2> );
}


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

template< class TT > inline void vector1_set( utility::vector1<TT> & v, platform::Size const & i, TT const & val ) { v[i] = val; }
template< class TT > inline platform::Size vector1_len( utility::vector1<TT> & v ) { return v.size(); }

template< class TT > inline std::string vector1_str( utility::vector1<TT> & v ) { std::ostringstream s; s<<v; return s.str(); }

template< class TT > inline typename utility::vector1<TT>::iterator vector1_begin( utility::vector1<TT> & v ) { return v.begin(); }
template< class TT > inline typename utility::vector1<TT>::iterator vector1_end  ( utility::vector1<TT> & v ) { return v.end(); }

template< class TT > inline void vector1_reserve( utility::vector1<TT> & v, platform::Size n) { v.reserve(n); }
template< class TT > inline void vector1_resize( utility::vector1<TT> & v, platform::Size n) { v.resize(n); }


template< class Htype, class CP, class CP_const>
void wrap_vector1(std::string name) {
  typedef utility::vector1<Htype> Ttype;
  typedef utility::vectorL<1,Htype, std::allocator<Htype> > Btype;
  typedef std::vector<Htype> Vtype;
  bp::class_<Ttype>(name.c_str())
    .def( bp::init< platform::Size >() )
    .def( bp::init< utility::vector1<Htype> const & >() )
    // .def( bp::init< platform::Size, TT >() )
    .def("__getitem__"
        , (Htype const & (Ttype::*)(platform::Size const) const)( &Ttype::at )
        , CP_const()    )
    .def("__getitem__"
        , (Htype & (Ttype::*)(platform::Size const))( &Ttype::at )
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
void wrap_vector1_part(const char * name) {
  typedef utility::vector1<Htype> Ttype;
  typedef utility::vectorL<1,Htype, std::allocator<Htype> > Btype;
  typedef std::vector<Htype> Vtype;
  bp::class_<Ttype>(name)
    .def("__getitem__"
        , (Htype const & (Ttype::*)(platform::Size const) const)( &Ttype::at )
        , CP_const()    )
    .def("__getitem__"
        , (Htype & (Ttype::*)(platform::Size const))( &Ttype::at )
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

template< class Key, class Val >
struct map_wrapper
{
  typedef std::map<Key,Val> Map;

  static boost::python::list keys(Map const& self)
  {
    boost::python::list t;

    for(typename Map::const_iterator it = self.begin(); it!=self.end(); ++it)
    {
      t.append(it->first);
    }

    return t;
  }

  static boost::python::list values(Map const& self)
  {
    boost::python::list t;

    for(typename Map::const_iterator it = self.begin(); it!=self.end(); ++it)
    {
      t.append(it->second);
    }

    return t;
  }

  static boost::python::list items(Map const& self)
  {
    boost::python::list t;

    for(typename Map::const_iterator it = self.begin(); it!=self.end(); ++it)
    {
        t.append( boost::python::make_tuple(it->first, it->second) );
    }

    return t;
  }
};

template< class Key, class Val >
void wrap_std_map(std::string name)
{
	bp::class_< std::map< Key, Val > >(name.c_str())
		.def( bp::map_indexing_suite< std::map< Key, Val > >())
    .def( "keys", &map_wrapper<Key, Val>::keys )
    .def( "values", &map_wrapper<Key, Val>::values )
    .def( "items", &map_wrapper<Key, Val>::items )
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
void wrap_std_set(std::string name)
{
	typedef std::set<Htype> Ttype;
	bp::class_<Ttype>(name.c_str())
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


template <class T>
void expose_number_type(std::string name)
{
    typedef bp::return_value_policy< bp::reference_existing_object > CP_REF;
    typedef bp::return_value_policy< bp::copy_const_reference >      CP_CCR;
    typedef bp::return_value_policy< bp::copy_non_const_reference >  CP_CNCR;

#ifdef PYROSETTA_NUMPY
    // conversion of array scalars
    //
#ifdef DEBUG
    std::cerr << "Registering array scalar wrapper. name: " << name << " typenum: " << get_typenum(T()) << std::endl;
#endif

    bp::converter::registry::push_back(
        array_scalar_converter<T>::check_array_scalar,
        array_scalar_converter<T>::convert_array_scalar,
        bp::type_id<T>());

    bp::converter::registry::push_back(
        array_scalar_converter<T>::check_zero_dim_array,
        array_scalar_converter<T>::convert_zero_dim_array,
        bp::type_id<T>());
#endif

    wrap_vector1<numeric::xyzVector<T>,    CP_CNCR, CP_CCR>("vector1_xyzVector_" + name);
}

template <class T1, class T2>
void expose_pair_types(std::string name_1, std::string name_2)
{
  typedef bp::return_value_policy< bp::reference_existing_object > CP_REF;
  typedef bp::return_value_policy< bp::copy_const_reference >      CP_CCR;
  typedef bp::return_value_policy< bp::copy_non_const_reference >  CP_CNCR;

  wrap_std_pair<T1, T2>("pair_" + name_1 + "_" + name_2);
  wrap_vector1<std::pair<T1, T2>, CP_CNCR, CP_CCR>("vector1_pair_" + name_1 + "_" + name_2);
  wrap_std_map< T1, T2 >("map_" + name_1 + "_" + name_2);
}


// Some subclassing testing functions ------------------------------------------
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


// Python Char/String arguments overload demo ------------------------------------------------------
void _test_char_string_args_(int c)
{
    std::cout << "_test_char_string_args_::char_version: c = " << c << std::endl;
}

void _test_char_string_args_(std::string s)
{
    std::cout << "_test_char_string_args_::string_version: s = " << s << std::endl;
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


void __core_by_hand_beginning__()
{
    // Testing functions
    bp::def("Q_Test_CI1B", Q_Test_CI1B);
    bp::def("Q_Test_EnergyMethodCreator", Q_Test_EnergyMethodCreator);

    using namespace utility::pointer;
    typedef bp::return_value_policy< bp::reference_existing_object > CP_REF;
    typedef bp::return_value_policy< bp::copy_const_reference >      CP_CCR;
    typedef bp::return_value_policy< bp::copy_non_const_reference >  CP_CNCR;

    // bp::class_< vector1<vector1<platform::Size> > >("utility___vec1_vec1_size")
    //   .def(bp::vector_indexing_suite< vector1<vector1<platform::Size> > >() );

    // bp::class_< access_ptr< core::chemical::AtomTypeSet const   > >("core___chemical___AtomTypeSetCAP");
    // bp::class_< access_ptr< core::chemical::ResidueType const   > >("core___chemical___ResidueTypeCAP");
    // bp::class_< access_ptr< core::chemical::ResidueTypeSet const> >("core___chemical___ResidueTypeSetCAP");
    // bp::class_< access_ptr< core::chemical::MMAtomTypeSet const > >("core___chemical___MMAtomTypeSetCAP");
    // bp::class_< access_ptr< core::coarse::Translator const      > >("core___coarse___TranslatorCAP");
    // bp::class_< access_ptr< core::coarse::CoarseEtable const    > >("core___coarse___CoarseEtableCAP");

    using namespace core::chemical;
    //using namespace core::coarse;

    // // old code - only for compatibility with previous verisons - deprecated, will be removed in the future...
    // bp::def("utility___getCAP"
    //      , (  ResidueTypeSet const * (*)( access_ptr<ResidueTypeSet const> )  )( & getCAP<ResidueTypeSet const> )
    //      , bp::return_value_policy< bp::reference_existing_object >() );
    // bp::def("utility___getCAP"
    //      , (  AtomTypeSet const * (*)( access_ptr<AtomTypeSet const> )  )( & getCAP<AtomTypeSet const> )
    //      , bp::return_value_policy< bp::reference_existing_object >() );
    // bp::def("utility___getCAP"
    //      , (  ResidueType const * (*)( access_ptr<ResidueType const> )  )( & getCAP<ResidueType const> )
    //      , bp::return_value_policy< bp::reference_existing_object >() );
    // bp::def("utility___getCAP"
    //      , (  MMAtomTypeSet const * (*)( access_ptr<MMAtomTypeSet const> )  )( & getCAP<MMAtomTypeSet const> )
    //      , bp::return_value_policy< bp::reference_existing_object >() );
    // bp::def("utility___getCAP"
    //      , (  Translator const * (*)( access_ptr<Translator const> )  )( & getCAP<Translator const> )
    //      , bp::return_value_policy< bp::reference_existing_object >() );
    // bp::def("utility___getCAP"
    //      , (  CoarseEtable const * (*)( access_ptr<CoarseEtable const> )  )( & getCAP<CoarseEtable const> )
    //      , bp::return_value_policy< bp::reference_existing_object >() );

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


	// Some wrapping test funtions/demos -----------------------------------------
    typedef void ( * _test_char_string_args_int_)(int);
    typedef void ( * _test_char_string_args_string_)(std::string);
    bp::def("Q_test_char_string_args_", _test_char_string_args_int_( &_test_char_string_args_) );
    bp::def("Q_test_char_string_args_", _test_char_string_args_string_( &_test_char_string_args_) );


    // Wraping for derived Demo
    bp::class_< DemoBase >("DemoBase")
    .def("testMethod", &DemoBase::testMethod)
    ;

    bp::class_< DemoDerived, bp::bases<DemoBase> >("DemoDerived")
    .def("testMethod", &DemoDerived::testMethod)
    ;

    boost::python::class_<DemoBase, BaseWrap, boost::noncopyable>( "Base" );

    bp::def("DemoTesterFunction", DemoTesterFunction);

    boost::python::class_<OOO> my_obj("OOO");
    //boost::python::object my_obj;

    boost::python::scope within(my_obj);
    bp::def("Q_Test_CI1B", Q_Test_CI1B);
    bp::def("Q_Test_EnergyMethodCreator", Q_Test_EnergyMethodCreator);


	// Dummy bindings to simplify import orders
	// boost::python::class_< utility::SingletonBase<core::chemical::ChemicalManager>, boost::noncopyable >( "__utility_SingletonBase_core_chemical_ChemicalManager__");
	// boost::python::class_< utility::SingletonBase<core::scoring::ScoringManager>, boost::noncopyable >( "__utility_SingletonBase_core_scoring_ScoringManager__");
	// boost::python::class_< utility::SingletonBase<core::scoring::constraints::ConstraintIO>, boost::noncopyable >( "__utility_SingletonBase_core_scoring_constraints_ConstraintIO__");
	// boost::python::class_< utility::SingletonBase<core::scoring::constraints::ConstraintFactory>, boost::noncopyable >( "__utility_SingletonBase_core_scoring_constraints_ConstraintFactory__");
	// boost::python::class_< utility::SingletonBase<core::io::silent::SilentStructFactory>, boost::noncopyable >( "__utility_SingletonBase_core_io_silent_SilentStructFactory__");
	// boost::python::class_< utility::SingletonBase<core::pack::task::operation::ResLvlTaskOperationFactory>, boost::noncopyable >( "__utility_SingletonBase_core_pack_task_operation_ResLvlTaskOperationFactory__");
}
