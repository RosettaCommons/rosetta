// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/tensorflow_manager/TFDataTypeDetector.hh
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
/// @author Jack Maguire, jackmaguire1444@gmail.com


#ifndef INCLUDED_basic_tensorflow_manager_TFDataTypeDetector_hh
#define INCLUDED_basic_tensorflow_manager_TFDataTypeDetector_hh

#include <utility/exit.hh>

// External headers
#include <tensorflow/c/c_api.h>

namespace basic {
namespace tensorflow_manager {

template< typename T >
struct TFDataTypeDetector {
	//This is only called if there is no special instantiation below
	TF_DataType value = TF_FLOAT;

	static void validate(){
		utility_exit_with_message( "TFDataTypeDetector has no instantiation for this T" );
	}
};

template<>
struct TFDataTypeDetector< float > {
	TF_DataType value = TF_FLOAT;

	static void validate(){}
};

template<>
struct TFDataTypeDetector< double > {
	TF_DataType value = TF_DOUBLE;

	static void validate(){}
};

template<>
struct TFDataTypeDetector< bool > {
	TF_DataType value = TF_BOOL;

	static void validate(){}
};

template<>
struct TFDataTypeDetector< int > {
	TF_DataType value = TF_INT32;

	static void validate(){ runtime_assert( sizeof(int) == 4 ); }
};

template<>
struct TFDataTypeDetector< unsigned int > {
	TF_DataType value = TF_UINT32;

	static void validate(){ runtime_assert( sizeof(unsigned int) == 4 ); }
};

template<>
struct TFDataTypeDetector< long int > {
	TF_DataType value = TF_INT64;

	static void validate(){ runtime_assert( sizeof(long int) == 8 ); }
};

template<>
struct TFDataTypeDetector< unsigned long int > {
	TF_DataType value = TF_UINT64;

	static void validate(){ runtime_assert( sizeof(unsigned long int) == 8 ); }
};

} //tensorflow_manager
} //basic

#endif //INCLUDED_basic_tensorflow_manager_TFDataTypeDetector_hh
