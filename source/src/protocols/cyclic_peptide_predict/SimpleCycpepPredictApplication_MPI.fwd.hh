// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication_MPI.fwd.hh
/// @brief Wrapper for SimpleCycpepPredictApplication that allows the app to use hierarchical MPI/pthreads based job distribution.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifdef USEMPI

#ifndef INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_MPI_fwd_hh
#define INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_MPI_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace cyclic_peptide_predict {

class SimpleCycpepPredictApplication_MPI;

typedef utility::pointer::shared_ptr< SimpleCycpepPredictApplication_MPI > SimpleCycpepPredictApplication_MPIOP;
typedef utility::pointer::shared_ptr< SimpleCycpepPredictApplication_MPI const > SimpleCycpepPredictApplication_MPICOP;

} //protocols
} //cyclic_peptide_predict

#endif //INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_MPI_fwd_hh

#endif //USEMPI
