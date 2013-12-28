// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_BasePairClassifier.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_RNA_BasePairClassifier_hh
#define INCLUDED_protocols_rna_RNA_BasePairClassifier_hh

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <core/scoring/rna/RNA_BaseDoubletClasses.hh>

#include <utility/vector1.hh>


//// C++ headers


namespace protocols {
namespace farna {

void
classify_base_pairs( core::pose::Pose const & pose, utility::vector1< core::scoring::rna::Base_pair> & base_pair_list, utility::vector1< bool > & is_bulged  );

core::Size
get_number_base_stacks(
	 core::pose::Pose const & pose_input
);

} //farna
} //protocols

#endif
