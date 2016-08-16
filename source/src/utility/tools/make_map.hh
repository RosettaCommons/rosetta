// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file   utility/tools/make_map.hh
/// @brief  Common function to build vector, vector0, vector1, map.
/// @author Sergey Lyskov

#ifndef INCLUDED_utility_tools_make_map_hh
#define INCLUDED_utility_tools_make_map_hh

#include <map>


namespace utility {
namespace tools {

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0)
{
	std::map<T1, T2> m;
	m[f0]=s0;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7, const T1 & f8, const T2 & s8)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7; m[f8]=s8;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7, const T1 & f8, const T2 & s8, const T1 & f9, const T2 & s9)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7; m[f8]=s8; m[f9]=s9;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7, const T1 & f8, const T2 & s8, const T1 & f9, const T2 & s9, const T1 & f10, const T2 & s10)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7; m[f8]=s8; m[f9]=s9; m[f10]=s10;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7, const T1 & f8, const T2 & s8, const T1 & f9, const T2 & s9, const T1 & f10, const T2 & s10, const T1 & f11, const T2 & s11)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7; m[f8]=s8; m[f9]=s9; m[f10]=s10; m[f11]=s11;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7, const T1 & f8, const T2 & s8, const T1 & f9, const T2 & s9, const T1 & f10, const T2 & s10, const T1 & f11, const T2 & s11, const T1 & f12, const T2 & s12)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7; m[f8]=s8; m[f9]=s9; m[f10]=s10; m[f11]=s11; m[f12]=s12;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7, const T1 & f8, const T2 & s8, const T1 & f9, const T2 & s9, const T1 & f10, const T2 & s10, const T1 & f11, const T2 & s11, const T1 & f12, const T2 & s12, const T1 & f13, const T2 & s13)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7; m[f8]=s8; m[f9]=s9; m[f10]=s10; m[f11]=s11; m[f12]=s12; m[f13]=s13;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7, const T1 & f8, const T2 & s8, const T1 & f9, const T2 & s9, const T1 & f10, const T2 & s10, const T1 & f11, const T2 & s11, const T1 & f12, const T2 & s12, const T1 & f13, const T2 & s13, const T1 & f14, const T2 & s14)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7; m[f8]=s8; m[f9]=s9; m[f10]=s10; m[f11]=s11; m[f12]=s12; m[f13]=s13; m[f14]=s14;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7, const T1 & f8, const T2 & s8, const T1 & f9, const T2 & s9, const T1 & f10, const T2 & s10, const T1 & f11, const T2 & s11, const T1 & f12, const T2 & s12, const T1 & f13, const T2 & s13, const T1 & f14, const T2 & s14, const T1 & f15, const T2 & s15)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7; m[f8]=s8; m[f9]=s9; m[f10]=s10; m[f11]=s11; m[f12]=s12; m[f13]=s13; m[f14]=s14; m[f15]=s15;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7, const T1 & f8, const T2 & s8, const T1 & f9, const T2 & s9, const T1 & f10, const T2 & s10, const T1 & f11, const T2 & s11, const T1 & f12, const T2 & s12, const T1 & f13, const T2 & s13, const T1 & f14, const T2 & s14, const T1 & f15, const T2 & s15, const T1 & f16, const T2 & s16)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7; m[f8]=s8; m[f9]=s9; m[f10]=s10; m[f11]=s11; m[f12]=s12; m[f13]=s13; m[f14]=s14; m[f15]=s15; m[f16]=s16;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7, const T1 & f8, const T2 & s8, const T1 & f9, const T2 & s9, const T1 & f10, const T2 & s10, const T1 & f11, const T2 & s11, const T1 & f12, const T2 & s12, const T1 & f13, const T2 & s13, const T1 & f14, const T2 & s14, const T1 & f15, const T2 & s15, const T1 & f16, const T2 & s16, const T1 & f17, const T2 & s17)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7; m[f8]=s8; m[f9]=s9; m[f10]=s10; m[f11]=s11; m[f12]=s12; m[f13]=s13; m[f14]=s14; m[f15]=s15; m[f16]=s16; m[f17]=s17;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7, const T1 & f8, const T2 & s8, const T1 & f9, const T2 & s9, const T1 & f10, const T2 & s10, const T1 & f11, const T2 & s11, const T1 & f12, const T2 & s12, const T1 & f13, const T2 & s13, const T1 & f14, const T2 & s14, const T1 & f15, const T2 & s15, const T1 & f16, const T2 & s16, const T1 & f17, const T2 & s17, const T1 & f18, const T2 & s18)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7; m[f8]=s8; m[f9]=s9; m[f10]=s10; m[f11]=s11; m[f12]=s12; m[f13]=s13; m[f14]=s14; m[f15]=s15; m[f16]=s16; m[f17]=s17; m[f18]=s18;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7, const T1 & f8, const T2 & s8, const T1 & f9, const T2 & s9, const T1 & f10, const T2 & s10, const T1 & f11, const T2 & s11, const T1 & f12, const T2 & s12, const T1 & f13, const T2 & s13, const T1 & f14, const T2 & s14, const T1 & f15, const T2 & s15, const T1 & f16, const T2 & s16, const T1 & f17, const T2 & s17, const T1 & f18, const T2 & s18, const T1 & f19, const T2 & s19)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7; m[f8]=s8; m[f9]=s9; m[f10]=s10; m[f11]=s11; m[f12]=s12; m[f13]=s13; m[f14]=s14; m[f15]=s15; m[f16]=s16; m[f17]=s17; m[f18]=s18; m[f19]=s19;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7, const T1 & f8, const T2 & s8, const T1 & f9, const T2 & s9, const T1 & f10, const T2 & s10, const T1 & f11, const T2 & s11, const T1 & f12, const T2 & s12, const T1 & f13, const T2 & s13, const T1 & f14, const T2 & s14, const T1 & f15, const T2 & s15, const T1 & f16, const T2 & s16, const T1 & f17, const T2 & s17, const T1 & f18, const T2 & s18, const T1 & f19, const T2 & s19, const T1 & f20, const T2 & s20)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7; m[f8]=s8; m[f9]=s9; m[f10]=s10; m[f11]=s11; m[f12]=s12; m[f13]=s13; m[f14]=s14; m[f15]=s15; m[f16]=s16; m[f17]=s17; m[f18]=s18; m[f19]=s19; m[f20]=s20;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7, const T1 & f8, const T2 & s8, const T1 & f9, const T2 & s9, const T1 & f10, const T2 & s10, const T1 & f11, const T2 & s11, const T1 & f12, const T2 & s12, const T1 & f13, const T2 & s13, const T1 & f14, const T2 & s14, const T1 & f15, const T2 & s15, const T1 & f16, const T2 & s16, const T1 & f17, const T2 & s17, const T1 & f18, const T2 & s18, const T1 & f19, const T2 & s19, const T1 & f20, const T2 & s20, const T1 & f21, const T2 & s21)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7; m[f8]=s8; m[f9]=s9; m[f10]=s10; m[f11]=s11; m[f12]=s12; m[f13]=s13; m[f14]=s14; m[f15]=s15; m[f16]=s16; m[f17]=s17; m[f18]=s18; m[f19]=s19; m[f20]=s20; m[f21]=s21;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7, const T1 & f8, const T2 & s8, const T1 & f9, const T2 & s9, const T1 & f10, const T2 & s10, const T1 & f11, const T2 & s11, const T1 & f12, const T2 & s12, const T1 & f13, const T2 & s13, const T1 & f14, const T2 & s14, const T1 & f15, const T2 & s15, const T1 & f16, const T2 & s16, const T1 & f17, const T2 & s17, const T1 & f18, const T2 & s18, const T1 & f19, const T2 & s19, const T1 & f20, const T2 & s20, const T1 & f21, const T2 & s21, const T1 & f22, const T2 & s22)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7; m[f8]=s8; m[f9]=s9; m[f10]=s10; m[f11]=s11; m[f12]=s12; m[f13]=s13; m[f14]=s14; m[f15]=s15; m[f16]=s16; m[f17]=s17; m[f18]=s18; m[f19]=s19; m[f20]=s20; m[f21]=s21; m[f22]=s22;
	return(m);
}

template<typename T1, typename T2>
std::map<T1, T2> make_map(const T1 &f0, const T2 & s0, const T1 & f1, const T2 & s1, const T1 & f2, const T2 & s2, const T1 & f3, const T2 & s3, const T1 & f4, const T2 & s4, const T1 & f5, const T2 & s5, const T1 & f6, const T2 & s6, const T1 & f7, const T2 & s7, const T1 & f8, const T2 & s8, const T1 & f9, const T2 & s9, const T1 & f10, const T2 & s10, const T1 & f11, const T2 & s11, const T1 & f12, const T2 & s12, const T1 & f13, const T2 & s13, const T1 & f14, const T2 & s14, const T1 & f15, const T2 & s15, const T1 & f16, const T2 & s16, const T1 & f17, const T2 & s17, const T1 & f18, const T2 & s18, const T1 & f19, const T2 & s19, const T1 & f20, const T2 & s20, const T1 & f21, const T2 & s21, const T1 & f22, const T2 & s22, const T1 & f23, const T2 & s23)
{
	std::map<T1, T2> m;
	m[f0]=s0; m[f1]=s1; m[f2]=s2; m[f3]=s3; m[f4]=s4; m[f5]=s5; m[f6]=s6; m[f7]=s7; m[f8]=s8; m[f9]=s9; m[f10]=s10; m[f11]=s11; m[f12]=s12; m[f13]=s13; m[f14]=s14; m[f15]=s15; m[f16]=s16; m[f17]=s17; m[f18]=s18; m[f19]=s19; m[f20]=s20; m[f21]=s21; m[f22]=s22; m[f23]=s23;
	return(m);
}


} // namespace utility
} // namespace tools

#endif // INCLUDED_utility_tools_make_map_hh

