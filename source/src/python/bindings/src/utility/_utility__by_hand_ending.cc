void __utility_by_hand_ending__()
{
    using namespace utility::pointer;
    typedef bp::return_value_policy< bp::reference_existing_object > CP_REF;
    typedef bp::return_value_policy< bp::copy_const_reference >      CP_CCR;
    typedef bp::return_value_policy< bp::copy_non_const_reference >  CP_CNCR;


    expose_basic_type<std::string>("string");

    expose_basic_type<int>("int");
    // conflict with SSize expose_basic_type<long>("long");
    // conflict with Size expose_basic_type<unsigned long>("ulong");
    expose_basic_type<unsigned int>("uint");
    expose_basic_type<char>("char");
    expose_basic_type<float>("float");
    expose_basic_type<double>("double");
    // conflict with Real expose_basic_type<double>("double");
    //expose_basic_type<unsigned>("unsigned");
    expose_basic_type<core::Size>("Size");
    expose_basic_type<core::SSize>("SSize");
    expose_basic_type<core::Real>("Real");

    //Need overload of bool_get so no:
    // expose_basic_type<bool>("bool");
    bp::class_< utility::vector1<bool> >("vector1_bool")
    .def("__len__",&vector1_len<bool> )
    .def("append", &vector1_bool_push )
    .def("__getitem__", &vector1_bool_get )
    .def("__iter__",bp::range(&vector1_begin<bool>,&vector1_end<bool>))
    ;
    wrap_vector1< utility::vector1<bool>, CP_REF, CP_REF >("vec1_vec1_bool");
}
