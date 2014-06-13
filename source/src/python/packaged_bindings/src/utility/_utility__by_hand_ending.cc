
void __utility_by_hand_ending__()
{
    using namespace pointer;
    typedef bp::return_value_policy< bp::reference_existing_object > CP_REF;
    typedef bp::return_value_policy< bp::copy_const_reference >      CP_CCR;
    typedef bp::return_value_policy< bp::copy_non_const_reference >  CP_CNCR;

    expose_basic_type<string>("string");

    expose_basic_type<int>("int");
    // conflict with SSize expose_basic_type<long>("long");
    // conflict with Size expose_basic_type<unsigned long>("ulong");
    expose_basic_type<unsigned int>("uint");
    expose_basic_type<char>("char");
    expose_basic_type<float>("float");
    // conflict with Real expose_basic_type<double>("double");
    //expose_basic_type<unsigned>("unsigned");
    expose_basic_type<core::Size>("Size");
    expose_basic_type<core::SSize>("SSize");
    expose_basic_type<core::Real>("Real");

    //Need overload of bool_get so no:
    // expose_basic_type<bool>("bool");
    bp::class_< vector1<bool> >("vector1_bool")
    .def("__len__",&vector1_len<bool> )
    .def("append", &vector1_bool_push )
    .def("__getitem__", &vector1_bool_get )
    .def("__iter__",bp::range(&vector1_begin<bool>,&vector1_end<bool>))
    ;
    wrap_vector1< vector1<bool>, CP_REF, CP_REF >("vec1_vec1_bool");


    expose_number_type<int>("int");
    // conflict with SSize expose_number_type<long>("long");
    // conflict with Size expose_number_type<unsigned long>("ulong");
    expose_number_type<unsigned int>("uint");
    expose_number_type<char>("char");
    expose_number_type<float>("float");
    // conflict with Real expose_number_type<double>("double");
    //expose_number_type<unsigned>("unsigned");
    expose_number_type<bool>("bool");

    expose_number_type<core::Size>("Size");
    expose_number_type<core::SSize>("SSize");
    expose_number_type<core::Real>("Real");

    expose_pair_types<std::string, std::string>("string", "string");

    expose_pair_types<int, std::string>("int", "string");
    expose_pair_types<std::string, int>("string", "int");

    expose_pair_types<core::Size, std::string>("Size", "string");
    expose_pair_types<std::string, core::Size>("string", "Size");

    expose_pair_types<core::SSize, std::string>("SSize", "string");
    expose_pair_types<std::string, core::SSize>("string", "SSize");

    expose_pair_types<core::Real, std::string>("Real", "string");
    expose_pair_types<std::string, core::Real>("string", "Real");

    wrap_vector1<core::scoring::ScoreType,  CP_CNCR, CP_CCR>("vector1_ScoreType");
    wrap_vector1<core::id::AtomID, CP_REF,CP_REF>("vector1_AtomID");

    wrap_vector1<core::pose::PoseOP,  CP_CNCR, CP_CCR>("vector1_PoseOP");
    wrap_vector1<core::pose::Pose,  CP_REF, CP_REF>("vector1_Pose");

    //if( true ) { // This code never executed, but we need it so compiler create operator ostream << vector1[core::conformation::Atom]
    //	cout << vector1< core::conformation::Atom > () ;
    //}
    //wrap_vector1_part<core::conformation::Atom, CP_REF, CP_REF>("vector1_core_conformation_Atom");
    wrap_vector1_part<core::conformation::Atom, CP_CNCR, CP_CCR>("vector1_core_conformation_Atom");
}
