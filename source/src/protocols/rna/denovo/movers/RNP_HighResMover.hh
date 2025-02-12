// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RNP_HighResMover.hh
/// @brief
/// @details
///
/// @author Kalli Kappel


#ifndef INCLUDED_protocols_rna_RNP_HighResMover_hh
#define INCLUDED_protocols_rna_RNP_HighResMover_hh

#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>

//// C++ headers
#include <string>

#include <protocols/rna/denovo/movers/RNA_FragmentMover.fwd.hh>
#include <protocols/rna/denovo/movers/RNA_Minimizer.fwd.hh>
#include <core/import_pose/options/RNA_FragmentMonteCarloOptions.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>
#include <protocols/rna/movers/RNA_LoopCloser.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh> // AUTO IWYU For ScoreFunctionCOP


namespace protocols {
namespace rna {
namespace denovo {
namespace movers {


typedef utility::pointer::shared_ptr< RNA_FragmentMover > RNA_FragmentMoverOP;
typedef utility::pointer::shared_ptr< RNA_Minimizer > RNA_MinimizerOP;

/// @brief The RNA de novo structure modeling protocol
class RNP_HighResMover: public protocols::moves::Mover {
public:
	/// @brief Construct the protocol object given
	/// the RNA fragment library to use.
	RNP_HighResMover( RNA_FragmentMoverOP rna_fragment_mover,
		protocols::rna::movers::RNA_LoopCloserOP rna_loop_closer,
		core::import_pose::options::RNA_FragmentMonteCarloOptionsCOP options );

	~RNP_HighResMover() override;

	void initialize( core::pose::Pose const & pose );

	void apply( core::pose::Pose & pose ) override;

	std::string get_name() const override;

private:
	//private methods
	void protein_sidechain_packing( core::pose::Pose & pose,
		core::scoring::ScoreFunctionCOP const & packer_scorefxn ) const;

	void
	rnp_docking( core::pose::Pose & pose ) const;

	//data
	RNA_FragmentMoverOP rna_fragment_mover_;
	protocols::rna::movers::RNA_LoopCloserOP rna_loop_closer_;
	core::import_pose::options::RNA_FragmentMonteCarloOptionsCOP options_;

	core::Size const fragment_cycles_;
	core::Size const frag_size_;
	core::Real protein_sc_pack_cutoff_;
	bool is_init_;
	core::kinematics::FoldTree rnp_docking_ft_;
	protocols::rigid::RigidBodyPerturbMoverOP  rnp_docking_mover_;
	bool docking_;

}; // class RNP_HighResMover


} //movers
} //denovo
} //rna
} //protocols

#endif
