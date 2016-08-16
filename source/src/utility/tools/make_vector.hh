// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file   utility/tools/make_vector.hh
/// @brief  Common function to build vector, vector0, vector1, map.
/// @author Sergey Lyskov

#ifndef INCLUDED_utility_tools_make_vector_hh
#define INCLUDED_utility_tools_make_vector_hh

#include <vector>


namespace utility {
namespace tools {

template<typename T>
std::vector<T> make_vector(const T & i0)
{
	std::vector<T> v;
	v.push_back(i0);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7, const T & i8)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7); v.push_back(i8);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7, const T & i8, const T & i9)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7); v.push_back(i8); v.push_back(i9);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7, const T & i8, const T & i9, const T & i10)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7); v.push_back(i8); v.push_back(i9); v.push_back(i10);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7, const T & i8, const T & i9, const T & i10, const T & i11)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7); v.push_back(i8); v.push_back(i9); v.push_back(i10); v.push_back(i11);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7, const T & i8, const T & i9, const T & i10, const T & i11, const T & i12)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7); v.push_back(i8); v.push_back(i9); v.push_back(i10); v.push_back(i11); v.push_back(i12);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7, const T & i8, const T & i9, const T & i10, const T & i11, const T & i12, const T & i13)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7); v.push_back(i8); v.push_back(i9); v.push_back(i10); v.push_back(i11); v.push_back(i12); v.push_back(i13);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7, const T & i8, const T & i9, const T & i10, const T & i11, const T & i12, const T & i13, const T & i14)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7); v.push_back(i8); v.push_back(i9); v.push_back(i10); v.push_back(i11); v.push_back(i12); v.push_back(i13); v.push_back(i14);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7, const T & i8, const T & i9, const T & i10, const T & i11, const T & i12, const T & i13, const T & i14, const T & i15)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7); v.push_back(i8); v.push_back(i9); v.push_back(i10); v.push_back(i11); v.push_back(i12); v.push_back(i13); v.push_back(i14); v.push_back(i15);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7, const T & i8, const T & i9, const T & i10, const T & i11, const T & i12, const T & i13, const T & i14, const T & i15, const T & i16)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7); v.push_back(i8); v.push_back(i9); v.push_back(i10); v.push_back(i11); v.push_back(i12); v.push_back(i13); v.push_back(i14); v.push_back(i15); v.push_back(i16);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7, const T & i8, const T & i9, const T & i10, const T & i11, const T & i12, const T & i13, const T & i14, const T & i15, const T & i16, const T & i17)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7); v.push_back(i8); v.push_back(i9); v.push_back(i10); v.push_back(i11); v.push_back(i12); v.push_back(i13); v.push_back(i14); v.push_back(i15); v.push_back(i16); v.push_back(i17);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7, const T & i8, const T & i9, const T & i10, const T & i11, const T & i12, const T & i13, const T & i14, const T & i15, const T & i16, const T & i17, const T & i18)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7); v.push_back(i8); v.push_back(i9); v.push_back(i10); v.push_back(i11); v.push_back(i12); v.push_back(i13); v.push_back(i14); v.push_back(i15); v.push_back(i16); v.push_back(i17); v.push_back(i18);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7, const T & i8, const T & i9, const T & i10, const T & i11, const T & i12, const T & i13, const T & i14, const T & i15, const T & i16, const T & i17, const T & i18, const T & i19)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7); v.push_back(i8); v.push_back(i9); v.push_back(i10); v.push_back(i11); v.push_back(i12); v.push_back(i13); v.push_back(i14); v.push_back(i15); v.push_back(i16); v.push_back(i17); v.push_back(i18); v.push_back(i19);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7, const T & i8, const T & i9, const T & i10, const T & i11, const T & i12, const T & i13, const T & i14, const T & i15, const T & i16, const T & i17, const T & i18, const T & i19, const T & i20)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7); v.push_back(i8); v.push_back(i9); v.push_back(i10); v.push_back(i11); v.push_back(i12); v.push_back(i13); v.push_back(i14); v.push_back(i15); v.push_back(i16); v.push_back(i17); v.push_back(i18); v.push_back(i19); v.push_back(i20);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7, const T & i8, const T & i9, const T & i10, const T & i11, const T & i12, const T & i13, const T & i14, const T & i15, const T & i16, const T & i17, const T & i18, const T & i19, const T & i20, const T & i21)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7); v.push_back(i8); v.push_back(i9); v.push_back(i10); v.push_back(i11); v.push_back(i12); v.push_back(i13); v.push_back(i14); v.push_back(i15); v.push_back(i16); v.push_back(i17); v.push_back(i18); v.push_back(i19); v.push_back(i20); v.push_back(i21);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7, const T & i8, const T & i9, const T & i10, const T & i11, const T & i12, const T & i13, const T & i14, const T & i15, const T & i16, const T & i17, const T & i18, const T & i19, const T & i20, const T & i21, const T & i22)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7); v.push_back(i8); v.push_back(i9); v.push_back(i10); v.push_back(i11); v.push_back(i12); v.push_back(i13); v.push_back(i14); v.push_back(i15); v.push_back(i16); v.push_back(i17); v.push_back(i18); v.push_back(i19); v.push_back(i20); v.push_back(i21); v.push_back(i22);
	return(v);
}

template<typename T>
std::vector<T> make_vector(const T & i0, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6, const T & i7, const T & i8, const T & i9, const T & i10, const T & i11, const T & i12, const T & i13, const T & i14, const T & i15, const T & i16, const T & i17, const T & i18, const T & i19, const T & i20, const T & i21, const T & i22, const T & i23)
{
	std::vector<T> v;
	v.push_back(i0); v.push_back(i1); v.push_back(i2); v.push_back(i3); v.push_back(i4); v.push_back(i5); v.push_back(i6); v.push_back(i7); v.push_back(i8); v.push_back(i9); v.push_back(i10); v.push_back(i11); v.push_back(i12); v.push_back(i13); v.push_back(i14); v.push_back(i15); v.push_back(i16); v.push_back(i17); v.push_back(i18); v.push_back(i19); v.push_back(i20); v.push_back(i21); v.push_back(i22); v.push_back(i23);
	return(v);
}


} // namespace utility
} // namespace tools

#endif // INCLUDED_utility_tools_make_vector_hh

