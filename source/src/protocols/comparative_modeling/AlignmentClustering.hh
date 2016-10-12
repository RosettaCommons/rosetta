// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
///
/// @brief  Divide input alns into clusters based on gdtmm comparison of partial models. Uses the top e-value or hh-value sample as the cluster center.
/// @author  TJ Bcunette

#ifndef INCLUDED_protocols_comparative_modeling_AlignmentClustering_hh
#define INCLUDED_protocols_comparative_modeling_AlignmentClustering_hh

#include <protocols/comparative_modeling/AlignmentClustering.fwd.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <map>

#include <core/sequence/SequenceAlignment.fwd.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.fwd.hh>


namespace protocols {
namespace comparative_modeling {

class AlignmentCluster : public utility::pointer::ReferenceCount {
private:
	utility::vector1<core::sequence::SequenceAlignment> alns;
public:
	AlignmentCluster(core::sequence::SequenceAlignment & aln_in);
	~AlignmentCluster() override;
	void add_aln(core::sequence::SequenceAlignment & aln_in);
	core::sequence::SequenceAlignment get_aln(core::Size index);
	core::Real size();
	core::sequence::SequenceAlignment get_clusterCenter();
	void output(std::ostream & alignment_out);
	void merge(AlignmentClusterOP cluster_in);
	core::Real overlap(AlignmentClusterOP cluster_in);
};
class AlignmentClustering : public utility::pointer::ReferenceCount {
public:
	AlignmentClustering();
	~AlignmentClustering() override;
private:
	utility::vector1<AlignmentClusterOP> cluster(utility::vector1< utility::vector1< core::Real > > & gdtmms, utility::vector1<core::sequence::SequenceAlignment> & rankedAlignments, core::Real threshold_gdt);
	std::map< std::string, core::pose::Pose > poses_from_cmd_line(utility::vector1< std::string > const & fn_list);
	utility::vector1< core::sequence::SequenceAlignment > generateRankedAlignments(std::map <std::string,core::sequence::SequenceAlignment> & alns, core::Real THRESHOLD_FOR_E_VAL);
};
} // comparative_modeling
} // protocols
#endif

