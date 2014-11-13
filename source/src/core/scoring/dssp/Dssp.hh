// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
// @brief
// @author olange: ported from original bblum-rosetta++ version $


#ifndef INCLUDED_core_scoring_dssp_Dssp_hh
#define INCLUDED_core_scoring_dssp_Dssp_hh

// Unit Headers
//#include <core/scoring/dssp/Dssp.fwd.hh>

// Package Headers
// AUTO-REMOVED #include <core/scoring/dssp/StrandPairing.hh>


// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
// AUTO-REMOVED #include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/FArray1A.fwd.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

//// C++ headers
//#include <cstdlib>
#include <string>
#include <cmath>

#include <core/scoring/dssp/StrandPairing.fwd.hh>
#include <utility/vector1.hh>

//#include <vector>

namespace core {
namespace scoring {
namespace dssp {

///@details Non-protein residues get a DSSP value of ' '
class Dssp {
public:

	Dssp( core::pose::Pose const& );
	~Dssp();

	void dssp_reduced_IG_as_L_if_adjcent_H( ObjexxFCL::FArray1_char &secstruct );
	void dssp_reduced_IG_as_L( ObjexxFCL::FArray1_char &secstruct );

	void dssp_reduced( ObjexxFCL::FArray1_char &secstruct );
	void dssp_reduced();
	void dssp_featurizer( ObjexxFCL::FArray1_char &secstruct );
	void dssp( ObjexxFCL::FArray1_char &dssp_secstruct );
	bool paired( core::Size res1, core::Size res2, bool antiparallel );

	StrandPairingSet const& strand_pairing_set() {
		return *pair_set_;
	}

	void insert_ss_into_pose(             core::pose::Pose & pose, bool recompute=true );
	void insert_dssp_ss_into_pose(        core::pose::Pose & pose, bool recompute=true );
	void insert_edge_ss_into_pose(        core::pose::Pose & pose, bool recompute=true );
	void insert_ss_into_pose_no_IG_helix( core::pose::Pose & pose, bool recompute=true );

	char get_dssp_secstruct( core::Size resid );

	std::string get_dssp_secstruct();
	std::string get_dssp_reduced_IG_as_L_secstruct();

	float bb_pair_score( Size res1, Size res2 );

	Size num_pairings(Size resi ) const;
	bool in_paired_strands(Size res1, Size res2 ) const;

private:
	void compute( core::pose::Pose const& );


	ObjexxFCL::FArray1D_char dssp_secstruct_;
	StrandPairingSetOP pair_set_;

	ObjexxFCL::FArray2D_float hbond_bb_pair_score_;
};

extern void fill_hbond_bb_pair_score_dssp( core::pose::Pose const&, ObjexxFCL::FArray2D_float &hbond_bb_pair_score_ );

} //dssp
} //scoring
} //core

#endif
