// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// Copyright (c) 2016 Sergey Lyskov <sergey.lyskov@jhu.edu>
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.

/// @file   binder/test/T42.stl.names.hpp
/// @brief  Binder self-test file. Tests for bindings of STL classes.
/// @author Sergey Lyskov

#ifndef _INCLUDED_T42_stl_names_hpp_
#define _INCLUDED_T42_stl_names_hpp_

#include <set>
#include <map>
#include <list>
#include <deque>
#include <vector>
#include <forward_list>
#include <unordered_set>
#include <unordered_map>

void foo(std::set<float> ) {}
void foo(std::list<double> ) {}
void foo(std::forward_list<float> ) {}
void foo(std::map<float,int> ) {}
void foo(std::vector<double> ) {}
void foo(std::deque<int> ) {}

void foo(std::multimap<float, int> ) {}
void foo(std::multiset<float> ) {}

void foo(std::unordered_set<float> ) {}
void foo(std::unordered_map<float,int> ) {}

void foo(std::unordered_multiset<float> ) {}
void foo(std::unordered_multimap<float, int> ) {}

#endif // _INCLUDED_T42_stl_names_hpp_
