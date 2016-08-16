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

using utility::vector1;
using core::Size;
using core::Real;
using std::string;
using std::map;
using core::pose::Pose;
// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core::sequence;

class AlignmentCluster : public utility::pointer::ReferenceCount {
private:
	vector1<SequenceAlignment> alns;
public:
	AlignmentCluster(SequenceAlignment & aln_in);
	virtual ~AlignmentCluster();
	void add_aln(SequenceAlignment & aln_in);
	SequenceAlignment get_aln(Size index);
	Real size();
	SequenceAlignment get_clusterCenter();
	void output(std::ostream & alignment_out);
	void merge(AlignmentClusterOP cluster_in);
	Real overlap(AlignmentClusterOP cluster_in);
};
class AlignmentClustering : public utility::pointer::ReferenceCount {
public:
	AlignmentClustering();
	virtual ~AlignmentClustering();
private:
	vector1<AlignmentClusterOP> cluster(vector1< vector1< Real > > & gdtmms, vector1<SequenceAlignment> & rankedAlignments, Real threshold_gdt);
	map< string, Pose > poses_from_cmd_line(utility::vector1< std::string > const & fn_list);
	vector1< SequenceAlignment > generateRankedAlignments(map <string,SequenceAlignment> & alns, Real THRESHOLD_FOR_E_VAL);
};
} // comparative_modeling
} // protocols
#endif

