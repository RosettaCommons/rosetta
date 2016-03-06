namespace bp = boost::python;

// core/io/Remarks -----------------------------------------------------------------------------------------------------------------
core::io::RemarkInfo& remark_vector_at(core::io::Remarks & v, platform::Size i ) { return v.at(i); }

void remark_vector_set(core::io::Remarks & v, platform::Size const & i, core::io::RemarkInfo const & val ) { v[i] = val; }

platform::Size remark_vector_len(core::io::Remarks & v ) { return v.size(); }

void remark_vector_append(core::io::Remarks & v, core::io::RemarkInfo const & val) { v.push_back(val); }

std::string remark_vector_repr(core::io::Remarks & v ) { std::ostringstream s; s<<v; return s.str(); }

core::io::Remarks::iterator remark_vector_begin(core::io::Remarks & v ) { return v.begin(); }
core::io::Remarks::iterator remark_vector_end  (core::io::Remarks & v ) { return v.end(); }

void remark_vector_reserve(core::io::Remarks & v, platform::Size n) { v.reserve(n); }
void remark_vector_resize(core::io::Remarks & v, platform::Size n) { v.resize(n); }

void wrap_remarks()
{
    typedef bp::return_value_policy< bp::reference_existing_object > CP_REF;

	bp::class_<core::io::Remarks>("Remarks")
		.def( bp::init<>() )

		.def("__getitem__",  &remark_vector_at, CP_REF() )
		.def("__setitem__",  &remark_vector_set)

		.def("__len__", &remark_vector_len)

		.def("append", &remark_vector_append)

		.def("__str__", &remark_vector_repr)

		.def("__iter__", bp::range(&remark_vector_begin,&remark_vector_end) )


		.def("reserve", &remark_vector_reserve)
		.def("resize", &remark_vector_resize)
  ;
}


void __io_by_hand_ending__()
{
	wrap_remarks();
}
