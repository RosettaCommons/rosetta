// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
#pragma once

#include <ndarray.h>
#include <protocols/indexed_structure_store/Datatypes.hh>

namespace protocols { namespace indexed_structure_store {

inline ndarray::ArrayRef<float, 4, 1>
orient_array(ndarray::Array<ResidueEntry, 2> res){
	typedef ndarray::detail::ArrayAccess< ndarray::ArrayRef<float,4,1> > Access;

	ndarray::Vector< ndarray::Size, 4> shape = ndarray::makeVector(
		(int)res.getSize<0>(), (int)res.getSize<1>(), 4, 3);
	ndarray::Vector< ndarray::Offset, 4> strides = ndarray::makeVector(
		(int)(res.getStride<0>() * sizeof(ResidueEntry) / sizeof(float)),
		(int)(res.getStride<1>() * sizeof(ResidueEntry) / sizeof(float)),
		3, 1);

	return Access::construct(
		&res.getData()->orient.N[0],
		Access::Core::create(shape, strides, res.getManager())
	);
}

inline ndarray::ArrayRef<float, 3, 1>
orient_array(ndarray::Array<ResidueEntry, 1> res){
	typedef ndarray::detail::ArrayAccess< ndarray::ArrayRef<float,3,1> > Access;

	ndarray::Vector< ndarray::Size, 3> shape = ndarray::makeVector((int)res.size(), 4, 3);
	ndarray::Vector< ndarray::Offset, 3> strides = ndarray::makeVector(
		(int)(res.getStride<0>() * sizeof(ResidueEntry) / sizeof(float)), 3, 1);

	return Access::construct(
		&res[0].orient.N[0],
		Access::Core::create(shape, strides, res.getManager())
	);
}

inline ndarray::ArrayRef<float, 2, 1>
orient_array(ResidueEntry & src) {
	return ndarray::external(
		&src.orient.N[0],
		ndarray::makeVector(4, 3),
		ndarray::makeVector(3, 1)
	);
}
} }
