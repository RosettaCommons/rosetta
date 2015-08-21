// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Svm_rosetta.hh
/// @brief
/// @details
///
/// @author TJ Brunette


#ifndef INCLUDED_utility_libsvm_Svm_rosetta_hh
#define INCLUDED_utility_libsvm_Svm_rosetta_hh
#include <utility/libsvm/Svm.hh>
#include <utility/libsvm/Svm_rosetta.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
// type headers
#include <platform/types.hh>
#include <string>

namespace utility {
namespace libsvm {

using utility::vector1;
using std::string;
//using platform::Size;
//using platform::Real;

class Svm_node_rosetta: public utility::pointer::ReferenceCount {
public:
	Svm_node_rosetta(platform::Size index, platform::Real value);
	~Svm_node_rosetta();
	void set_index(platform::Size index){
		index_ = index;
	};
	void set_value(platform::Real value){
		value_ = value;
	};
	platform::Size index(){
		return(index_);
	};
	platform::Real value(){
		return(value_);
	};
private:
	platform::Size index_;
	platform::Real value_;
};


class Svm_rosetta: public utility::pointer::ReferenceCount {
public:
	Svm_rosetta(string model_filename);
	~Svm_rosetta();
	platform::Size get_nr_class();
	vector1< platform::Real > predict_probability(vector1< Svm_node_rosettaOP > & features);
	platform::Real predict( const vector1< Svm_node_rosettaOP > & features);
private:
	svm_model *svm_model_;
};
}//libsvm
}//utility

#endif
