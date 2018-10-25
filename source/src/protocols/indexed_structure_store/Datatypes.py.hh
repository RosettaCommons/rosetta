// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/indexed_structure_store/Datatypes.py.hh
/// @brief Pybind11 converter support for indexed_structure_store datatypes.
/// @author Alex Ford <fordas@uw.edu>
//
#pragma once

#include "stdint.h"
#include <iostream>

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

#include "ndarray.h"
#include "ndarray/pybind11.h"

#include <protocols/indexed_structure_store/Datatypes.hh>
#include <protocols/indexed_structure_store/nonpod_dtype_support.h>

DECLARE_NONPOD_DTYPE_SUPPORT(protocols::indexed_structure_store::ResidueBackboneEntry);
DECLARE_NONPOD_DTYPE_SUPPORT(protocols::indexed_structure_store::ResidueSidechainEntry);
DECLARE_NONPOD_DTYPE_SUPPORT(protocols::indexed_structure_store::ResidueOrientEntry);
DECLARE_NONPOD_DTYPE_SUPPORT(protocols::indexed_structure_store::ResidueEntry);
