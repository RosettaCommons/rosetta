// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/JobOutputIndex.hh
/// @brief  Definition of the %JobOutputIndex class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_JobOutputIndex_HH
#define INCLUDED_protocols_jd3_JobOutputIndex_HH

// Unit headers
#include <protocols/jd3/JobOutputIndex.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {

/// @brief The %JobOutputIndex holds the four relevant counts for describing
/// the "name" for an output (e.g. a Pose) that is going to be written out to disk
struct JobOutputIndex : utility::pointer::ReferenceCount
{
	// There might be several jobs that ran with the same "job tag" (e.g. the same starting structure), but
	// perhaps their outputs are not all going to be written out to disk. The prime_output_index
	// could be the nstruct id if you had one InnerLarvalJob and nstruct LarvalJobs that
	// pointed at it, but if you should have several InnerLarvalJobs for a single starting
	// structure, you may wish to number your output Poses differently (because nstruct id is insufficient).
	// The prime_output_index lets you assign an index for a Pose (or a group of Poses) that came out of
	// a single job.
	core::Size primary_output_index = 1;

	// For zero-padding at the beginning of the primary_output_index. You do not need to know exactly
	// how many primary outputs for a job there are so long as you can give an upper bound;
	// A default of 1000 (equivalent to 9999 in terms of the number of leading zeros) is commonly
	// used in Rosetta.
	core::Size n_primary_outputs = 1000;

	// If a job produces more than one output (e.g. several poses), you can group them together and label
	// them with the same "output index for job" id, and then differentiate them by
	// a secondary index.  E.g. If job 50 starting from 1ubq produces 3 output PDBs, they these
	// can be written out as 1ubq_0050_0001.pdb, 1ubq_0050_0002.pdb, etc.
	core::Size secondary_output_index = 1;

	// For zero padding; an upper bound on the number of local output indices for a job. If there
	// is only a single ouput Pose, e.g. for job 48 of 1ubq and local_output_index & n_local_outputs
	// are both set to 1, then the Pose will be written out to 1ubq_0048.pdb instead of 1ubq_0048_0001.pdb (assuming
	// that you are writing Poses as PDBs).
	// It is entirely reasonable for primary index 53 to have a different number of secondary outputs from
	// primary index 54; 53 could have 10, 54 could have 100.
	core::Size n_secondary_outputs = 1;
#ifdef    SERIALIZATION
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_JobOutputIndex )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_JobOutputIndex_HH
