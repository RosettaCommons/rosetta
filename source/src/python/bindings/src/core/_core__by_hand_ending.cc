
void __core_by_hand_ending__()
{
    using namespace utility::pointer;
    typedef bp::return_value_policy< bp::reference_existing_object > CP_REF;
    typedef bp::return_value_policy< bp::copy_const_reference >      CP_CCR;
    typedef bp::return_value_policy< bp::copy_non_const_reference >  CP_CNCR;

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

	expose_pair_types<int, int>("int", "int");
	expose_pair_types<core::Size, core::Size>("Size", "Size");
	expose_pair_types<core::SSize, core::SSize>("SSize", "SSize");
    expose_pair_types<core::Real, core::Real>("Real", "Real");

    expose_pair_types<int, std::string>("int", "string");
    expose_pair_types<std::string, int>("string", "int");

    expose_pair_types<core::Size, std::string>("Size", "string");
    expose_pair_types<std::string, core::Size>("string", "Size");

    expose_pair_types<core::SSize, std::string>("SSize", "string");
    expose_pair_types<std::string, core::SSize>("string", "SSize");

    expose_pair_types<core::Real, std::string>("Real", "string");
    expose_pair_types<std::string, core::Real>("string", "Real");


    expose_pair_types<core::Size, core::Real>("Size", "Real");


	// Some types for core/conformation/symmetry/SymmetryInfo.hh templates
    expose_pair_types<core::Size, std::string>("Size", "string");
	expose_pair_types<core::Size, core::conformation::symmetry::SymDof>("Size", "SymDof");
	//expose_pair_types<core::Size, core::conformation::symmetry::SymmetryInfo::WtedClones>("Size", "SymmetryInfo_WtedClones");
	expose_pair_types<char, std::pair<core::Size, core::Size> >("char", "pair_Size_Size");


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
