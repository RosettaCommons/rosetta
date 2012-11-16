
void __utility_by_hand_ending__()
{
    using namespace pointer;
    typedef bp::return_value_policy< bp::reference_existing_object > CP_REF;
    typedef bp::return_value_policy< bp::copy_const_reference >      CP_CCR;
    typedef bp::return_value_policy< bp::copy_non_const_reference >  CP_CNCR;

    wrap_vector1< vector1< std::size_t >,CP_REF,CP_REF>("vec1_vec1_Size");

    wrap_vector1<string,       CP_CNCR, CP_CCR>("vector1_string");
    wrap_vector1<int,          CP_CNCR, CP_CCR>("vector1_int");
    wrap_vector1<long,         CP_CNCR, CP_CCR>("vector1_long");
    wrap_vector1<unsigned long,CP_CNCR, CP_CCR>("vector1_ulong");
    wrap_vector1<unsigned int, CP_CNCR, CP_CCR>("vector1_uint");
    wrap_vector1<char,         CP_CNCR, CP_CCR>("vector1_char");
    wrap_vector1<float,        CP_CNCR, CP_CCR>("vector1_float");
    wrap_vector1<double,       CP_CNCR, CP_CCR>("vector1_double");
    wrap_vector1<unsigned,     CP_CNCR, CP_CCR>("vector1_unsigned");

    wrap_vector1<core::Size,   CP_CNCR, CP_CCR>("vector1_Size");
    wrap_vector1<core::SSize,  CP_CNCR, CP_CCR>("vector1_SSize");

    wrap_vector1<numeric::xyzVector_int,    CP_CNCR, CP_CCR>("vector1_xyzVector_int");
    wrap_vector1<numeric::xyzVector_uint,   CP_CNCR, CP_CCR>("vector1_xyzVector_uint");
    wrap_vector1<numeric::xyzVector_long,   CP_CNCR, CP_CCR>("vector1_xyzVector_long");
    wrap_vector1<numeric::xyzVector_ulong,  CP_CNCR, CP_CCR>("vector1_xyzVector_ulong");
    wrap_vector1<numeric::xyzVector_double, CP_CNCR, CP_CCR>("vector1_xyzVector_float");
    wrap_vector1<numeric::xyzVector_double, CP_CNCR, CP_CCR>("vector1_xyzVector_double");
    wrap_vector1<numeric::xyzVector_char,   CP_CNCR, CP_CCR>("vector1_xyzVector_char");
    wrap_vector1<numeric::xyzVector_uchar,  CP_CNCR, CP_CCR>("vector1_xyzVector_uchar");

    wrap_vector1<core::scoring::ScoreType,  CP_CNCR, CP_CCR>("vector1_ScoreType");

    wrap_vector1<core::id::AtomID, CP_REF,CP_REF>("vector1_AtomID");

    wrap_vector1<numeric::xyzVector_bool, CP_CNCR, CP_CCR>("vector1_xyzVector_bool");

    wrap_vector1<core::pose::PoseOP,  CP_CNCR, CP_CCR>("vector1_PoseOP");


    //if( true ) { // This code never executed, but we need it so compilere create operator ostream << vector1[core::conformation::Atom]
    //	cout << vector1< core::conformation::Atom > () ;
    //}
    //wrap_vector1_part<core::conformation::Atom, CP_REF, CP_REF>("vector1_core_conformation_Atom");
    wrap_vector1_part<core::conformation::Atom, CP_CNCR, CP_CCR>("vector1_core_conformation_Atom");

    // wrap_vector1<bool,         CP_CNCR,CP_CCR>("utility___vector1_bool");
    bp::class_< vector1<bool> >("vector1_bool")
    .def("__len__",&vector1_len<bool> )
    .def("append", &vector1_bool_push )
    .def("__getitem__", &vector1_bool_get )
    .def("__iter__",bp::range(&vector1_begin<bool>,&vector1_end<bool>))
    ;

    wrap_std_map< std::map<int, int> >("map_int_int");

    //std::set<int> s_int;  set_repr(s_int);
    wrap_std_set< int,           CP_CNCR, CP_CCR >("set_int");
    wrap_std_set< unsigned int,  CP_CNCR, CP_CCR >("set_uint");
    wrap_std_set< unsigned long, CP_CNCR, CP_CCR >("set_ulong");
    wrap_std_set< double,        CP_CNCR, CP_CCR >("set_double");
    wrap_std_set< char,          CP_CNCR, CP_CCR >("set_char");
    wrap_std_set< std::string,   CP_CNCR, CP_CCR >("set_string");
    wrap_std_set< core::Size,    CP_CNCR, CP_CCR >("set_Size");
    wrap_std_set< core::SSize,   CP_CNCR, CP_CCR >("set_SSize");
}
