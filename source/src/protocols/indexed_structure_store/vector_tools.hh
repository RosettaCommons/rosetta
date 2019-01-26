// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/vector_tools.hh
/// @brief Support templates for working with vectors
/// @author Alex Ford (fordas@uw.edu)
#pragma once

#include <vector>
#include <algorithm>
#include <ndarray.h>
#include <boost/range.hpp>

template <typename To, typename Collection, typename unop>
std::vector<To> v_map(Collection src, unop op) {
	std::vector<To> result;
	std::transform(
		src.begin(),src.end(),std::back_inserter(result), op);
	return result;
}

template<typename T, typename unop>
void v_erase_if(std::vector<T> vec, unop op) {
	vec.erase(
		std::remove_if(vec.begin(), vec.end(), op),
		vec.end());
}

template<typename T, typename Collection, typename unop>
std::vector<T> v_copy_if(Collection src, unop op) {
	std::vector<T> result;
	std::copy_if(
		src.begin(),src.end(),std::back_inserter(result), op);
	return result;
}

template<typename T>
ndarray::Array<T, 1, 1>
v_to_a(std::vector<T> vec) {
	if ( vec.empty() ) {
		return ndarray::external(
			(T*)( nullptr ),
			ndarray::makeVector(vec.size()), ndarray::makeVector(0)
		);
	}
	return ndarray::external(
		&vec.front(),
		ndarray::makeVector(vec.size()), ndarray::makeVector(1)
	);
}

template<typename T, typename Range>
std::vector<T>
r_to_v(Range r){
	return std::vector<T>( boost::begin(r), boost::end(r) );
}

template<typename T>
std::vector<T>
a_to_v(ndarray::Array<T, 1> arr){
	return r_to_v<T>(arr);
}
