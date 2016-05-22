// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication_MPI_JobResultsSummary.fwd.hh
/// @brief  Defines owning pointers for a small helper class used by the SimpleCycpepPredict_MPI class.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_MPI_JobResultsSummary_fwd_hh
#define INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_MPI_JobResultsSummary_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace cyclic_peptide_predict {

class SimpleCycpepPredictApplication_MPI_JobResultsSummary; // fwd declaration
typedef utility::pointer::shared_ptr< SimpleCycpepPredictApplication_MPI_JobResultsSummary > SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP;
typedef utility::pointer::shared_ptr< SimpleCycpepPredictApplication_MPI_JobResultsSummary const > SimpleCycpepPredictApplication_MPI_JobResultsSummaryCOP;
typedef utility::vector1<SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP> SimpleCycpepPredictApplication_MPI_JobResultsSummaryOPs;
typedef utility::vector1<SimpleCycpepPredictApplication_MPI_JobResultsSummaryCOP> SimpleCycpepPredictApplication_MPI_JobResultsSummaryCOPs;

} // namespace cyclic_peptide_predict
} // namespace protocols

#endif // INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_MPI_JobResultsSummary_fwd_hh

