// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta_binders.hpp
/// @brief  Support for custom binders for some Rosetta template classes
/// @author Sergey Lyskov


#ifndef _INCLUDED_rosetta_binders_hpp_
#define _INCLUDED_rosetta_binders_hpp_

#include <pybind11/operators.h>
#include <pybind11/stl_binders.h>

#include <utility/vectorL.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <type_traits>
#include <sstream>


// Adding Rosetta-specific cases to has_equal_operator_s template specialization
namespace pybind11 { namespace detail {

template <typename A>
struct has_equal_operator_s< utility::vector0<A> >
{
	static const bool value = has_equal_operator_s<A>::value;
};
template <typename A>
struct has_equal_operator_s< utility::vector1<A> >
{
	static const bool value = has_equal_operator_s<A>::value;
};

} }


namespace rosetta_binders {


template <typename T, class Allocator>
void vector_binder(pybind11::module &m, char const *name, char const * /*allocator name*/) {
	pybind11::vector_binder<T, Allocator, std::shared_ptr< std::vector<T, Allocator> > >(m, name);
}


template<typename Vector, platform::SSize L, typename T, typename Allocator>
class utility_vector_binder
{
	using SizeType = typename Vector::size_type;

	using Class_ = pybind11::class_<Vector, std::shared_ptr< Vector > >;


	// template<typename U = T, typename std::enable_if< std::is_constructible<U>{} >::type * = nullptr>
	// void maybe_constructible(Class_ &cl) {
	// 	cl.def(pybind11::init<>());
	// }
	// template<typename U = T, typename std::enable_if< !std::is_constructible<U>{} >::type * = nullptr>
	// void maybe_constructible(Class_ &cl) {}

	template<typename U = T, typename std::enable_if< std::is_default_constructible<U>::value >::type * = nullptr>
	void maybe_default_constructible(Class_ &cl) {
		cl.def(pybind11::init<SizeType>());
		cl.def("resize", (void (Vector::*)(SizeType count)) &Vector::resize, "changes the number of elements stored");

		/// Slicing protocol
		cl.def("__getitem__", [](Vector const &v, pybind11::slice slice) -> Vector * {
				pybind11::ssize_t start, stop, step, slicelength;
				if(!slice.compute(v.size(), &start, &stop, &step, &slicelength))
					throw pybind11::error_already_set();
				Vector *seq = new Vector(slicelength);
				for (int i=0; i<slicelength; ++i) {
					(*seq)[i] = v[start]; start += step;
				}
				return seq;
			});
	}
	template<typename U = T, typename std::enable_if< !std::is_default_constructible<U>::value >::type * = nullptr>
	void maybe_default_constructible(Class_ &) {}


	template<typename U = T, typename std::enable_if< std::is_copy_constructible<U>::value >::type * = nullptr>
	void maybe_copy_constructible(Class_ &cl) {
		cl.def(pybind11::init< Vector const &>());
	}
	template<typename U = T, typename std::enable_if< !std::is_copy_constructible<U>::value >::type * = nullptr>
	void maybe_copy_constructible(Class_ &) {}


	template<typename U = T, typename std::enable_if< pybind11::detail::has_equal_operator_s<U>::value >::type * = nullptr>
	void maybe_has_equal_operator(Class_ &cl) {
	    cl.def(pybind11::self == pybind11::self);
	    cl.def(pybind11::self != pybind11::self);

		cl.def("count", [](Vector const &v, T const & value) { return std::count(v.begin(), v.end(), value); }, "counts the elements that are equal to value");

		cl.def("remove", [](Vector &v, T const &t) {
				auto p = std::find(v.begin(), v.end(), t);
				if(p != v.end()) v.erase(p);
				else throw pybind11::value_error();
			}, "Remove the first item from the list whose value is x. It is an error if there is no such item.");

		cl.def("__contains__", [](Vector const &v, T const &t) { return std::find(v.begin(), v.end(), t) != v.end(); }, "return true if item in the container");
	}
	template<typename U = T, typename std::enable_if< !pybind11::detail::has_equal_operator_s<U>::value >::type * = nullptr>
	void maybe_has_equal_operator(Class_ &) {}


	template<typename U = T, typename std::enable_if< pybind11::detail::has_insertion_operator_s<U>::value >::type * = nullptr>
	void maybe_has_insertion_operator(std::string const &name, Class_ &cl) {
		cl.def("__repr__", [name](Vector &v) {
				std::ostringstream s;
				s << name << '[';
				for(SizeType i=v.l(); i<=v.u(); ++i) {
					s << v[i];
					if(i != v.u()) s << ", ";
				}
				s << ']';
 				return s.str();
			});

	}
	template<typename U = T, typename std::enable_if< !pybind11::detail::has_insertion_operator_s<U>::value >::type * = nullptr>
	void maybe_has_insertion_operator(std::string const &, Class_ &) {}


public:
	utility_vector_binder(pybind11::module &m, std::string const &lower_index, std::string const &name, std::string const & /*allocator name*/) {
		Class_ cl(m, std::string("vector"+ lower_index +'_' + name).c_str() );

		cl.def(pybind11::init<>());

		//maybe_constructible(cl);
		maybe_default_constructible(cl);
		maybe_copy_constructible(cl);

		// Element access
		cl.def("front", [](Vector &v) {
				if(v.size()) return v.front();
				else throw pybind11::index_error();
			}, "access the first element");
		cl.def("back", [](Vector &v) {
				if(v.size()) return v.back();
				else throw pybind11::index_error();
			}, "access the last element ");
		// Not needed, the operator[] is already providing bounds checking cl.def("at", (T& (Vector::*)(SizeType i)) &Vector::at, "access specified element with bounds checking");

		// Capacity, C++ style
		cl.def("max_size",      &Vector::max_size,      "returns the maximum possible number of elements");
		cl.def("reserve",       &Vector::reserve,       "reserves storage");
		cl.def("capacity",      &Vector::capacity,      "returns the number of elements that can be held in currently allocated storage");
		cl.def("shrink_to_fit", &Vector::shrink_to_fit, "reduces memory usage by freeing unused memory");

		// Modifiers, C++ style
		cl.def("clear", &Vector::clear, "clears the contents");

		// Modifiers, Python style
		cl.def("append", (void (Vector::*)(const T&)) &Vector::push_back, "adds an element to the end");
		cl.def("insert", [](Vector &v, SizeType i, const T&t) {v.insert(v.begin()+i, t);}, "insert an item at a given position");
		cl.def("extend", [](Vector &v, Vector &src) { v.reserve( v.size() + src.size() ); v.insert(v.end(), src.begin(), src.end()); }, "extend the list by appending all the items in the given vector");
		cl.def("pop", [](Vector &v) {
				if(v.size()) {
					T t = v.back();
					v.pop_back();
					return t;
				}
				else throw pybind11::index_error();
			}, "remove and return last item");

		cl.def("pop", [](Vector &v, SizeType i) {
				if(i >= v.size()) throw pybind11::index_error();
				T t = v[i];
				v.erase(v.begin() + i);
				return t;
			}, "remove and return item at index");

		cl.def("erase", [](Vector &v, SizeType i) {
				if(i >= v.size()) throw pybind11::index_error();
				v.erase(v.begin() + i);
			}, "erases element at index");


		// Python friendly bindings
		#ifdef PYTHON_ABI_VERSION // Python 3+
			cl.def("__bool__",    [](Vector &v) -> bool { return v.size() != 0; }); // checks whether the container has any elements in it
		#else
			cl.def("__nonzero__", [](Vector &v) -> bool { return v.size() != 0; }); // checks whether the container has any elements in it
		#endif

		cl.def("__getitem__", [](Vector const &v, SizeType i) {
				if(i >= v.size()) throw pybind11::index_error();
				return v[i];
			});

		cl.def("__setitem__", [](Vector &v, SizeType i, T const & t) {
				if(i >= v.size()) throw pybind11::index_error();
				v[i] = t;
			});

		cl.def("__len__", &Vector::size);

		cl.def("__iter__", [](Vector &v) { return pybind11::make_iterator(v.begin(), v.end()); },
			   pybind11::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */);

		/// Slicing protocol
		cl.def("__setitem__", [](Vector &v, pybind11::slice slice,  Vector const &value) {
				pybind11::ssize_t start, stop, step, slicelength;
				if(!slice.compute(v.size(), &start, &stop, &step, &slicelength))
					throw pybind11::error_already_set();
				if((size_t) slicelength != value.size())
					throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
				for(int i=0; i<slicelength; ++i) {
					v[start] = value[i]; start += step;
				}
			});

		cl.def("l", (SizeType (Vector::*)() const) &Vector::l, "lower index");
		cl.def("u", (SizeType (Vector::*)() const) &Vector::u, "upper index");

		// fixme later: add separate binding for to/from std::vector  cl.def("swap",   &Vector::swap, "swaps the contents");


		// Comparisons
		maybe_has_equal_operator(cl);

		// Printing
		maybe_has_insertion_operator(name, cl);
	}

	utility_vector_binder(pybind11::module &m, std::string const &name, std::string const & allocator) : utility_vector_binder(m, std::to_string(L), name, allocator) {}
};


template<platform::SSize L, typename T, typename A> using vectorL_binder = utility_vector_binder<utility::vectorL<L, T, A>, L, T, A>;
template< typename T, typename A > using vector0_binder = utility_vector_binder<utility::vector0<T, A>, 0, T, A>;
template< typename T, typename A > using vector1_binder = utility_vector_binder<utility::vector1<T, A>, 1, T, A>;


//using Vector = utility::vectorL<L, T, Allocator>;
//template< typename T, typename Allocator > using vector0_binder = vectorL_binder<0, T, Allocator>;
//template< typename T, typename Allocator > using vector1_binder = vectorL_binder<1, T, Allocator>;


// template< typename T, typename Allocator >
// class vector0_binder : public vectorL_binder<0, T, Allocator>
// {
// public:
// 	using vectorL_binder<0, T, Allocator>::vectorL_binder;
// 	vector0_binder(pybind11::module &m, std::string const &name, std::string const & allocator) : vector0_binder(m, "0", name, allocator) {}
// };


// template< typename T, typename Allocator >
// class vector1_binder : public vectorL_binder<1, T, Allocator>
// {
// public:
// 	using vectorL_binder<1, T, Allocator>::vectorL_binder;
// 	vector1_binder(pybind11::module &m, std::string const &name, std::string const & allocator) : vector1_binder(m, "1", name, allocator) {}
// };



		// cl.def("vector", (const class std::vector<class core::scoring::rna::data::RNA_Reactivity, class std::allocator<class core::scoring::rna::data::RNA_Reactivity> > & (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)() const) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::vector, "doc");
		// cl.def("vector", (class std::vector<class core::scoring::rna::data::RNA_Reactivity, class std::allocator<class core::scoring::rna::data::RNA_Reactivity> > & (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)()) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::vector, "doc");
		// cl.def("append", (class utility::vectorL<1, class core::scoring::rna::data::RNA_Reactivity, class std::allocator<class core::scoring::rna::data::RNA_Reactivity> > & (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)(const class utility::vectorL<1, class core::scoring::rna::data::RNA_Reactivity, class std::allocator<class core::scoring::rna::data::RNA_Reactivity> > &)) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::append, "doc", pybind11::arg("v"));
		// cl.def("add_back", (class utility::vectorL<1, class core::scoring::rna::data::RNA_Reactivity, class std::allocator<class core::scoring::rna::data::RNA_Reactivity> > & (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)(const class core::scoring::rna::data::RNA_Reactivity &)) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::add_back, "doc", pybind11::arg("t"));
		// cl.def("remove_back", (class utility::vectorL<1, class core::scoring::rna::data::RNA_Reactivity, class std::allocator<class core::scoring::rna::data::RNA_Reactivity> > & (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)()) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::remove_back, "doc");
		// cl.def("shrink", (void (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)()) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::shrink, "doc");
		// cl.def("has", (bool (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)(const unsigned long) const) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::has, "doc", pybind11::arg("i"));
		// cl.def("at", (const class core::scoring::rna::data::RNA_Reactivity & (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)(const unsigned long) const) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::at, "doc", pybind11::arg("i"));
		// cl.def("at", (class core::scoring::rna::data::RNA_Reactivity & (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)(const unsigned long)) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::at, "doc", pybind11::arg("i"));
		// cl.def("swap", (void (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)(class utility::vectorL<1, class core::scoring::rna::data::RNA_Reactivity, class std::allocator<class core::scoring::rna::data::RNA_Reactivity> > &)) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::swap, "doc", pybind11::arg("v"));



} // namespace rosetta_binders

#endif // _INCLUDED_rosetta_binders_hpp_
