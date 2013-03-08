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
/// @detailed
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
#include <core/types.hh>
#include <string>

namespace utility {
namespace libsvm {

using utility::vector1;
using core::Real;
using core::Size;
using std::string;

class Svm_node_rosetta: public utility::pointer::ReferenceCount {
 public:
	Svm_node_rosetta(Size index, Real value);
	~Svm_node_rosetta();
	void set_index(Size index){
		index_ = index;
	};
	void set_value(Real value){
		value_ = value;
	};
	Size index(){
		return(index_);
	};
	Real value(){
		return(value_);
	};
 private:
	Size index_;
	Real value_;
};


class Svm_rosetta: public utility::pointer::ReferenceCount {
 public:
  Svm_rosetta(string model_filename);
  ~Svm_rosetta();
  Size get_nr_class();
	vector1< Real > predict_probability(vector1< Svm_node_rosettaOP > features);
 private:
	svm_model *svm_model_;
};
}//libsvm
}//utility

#endif
