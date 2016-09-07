// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/seeded_abinitio/GrowPeptides.cc
/// @author Eva-Maria Strauch (evas01@u.washington.edu)

#ifndef INCLUDED_protocols_seeded_abinitio_GrowPeptides_hh
#define INCLUDED_protocols_seeded_abinitio_GrowPeptides_hh

#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/string_util.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <set>

namespace protocols {
namespace seeded_abinitio {

class GrowPeptides : public protocols::moves::Mover
{
public:
	GrowPeptides();

	~GrowPeptides() override;
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(     utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;

	// Undefined, commenting out to fix PyRosetta build   void add_chainbreakterm( bool ac );

	// Undefined, commenting out to fix PyRosetta build   bool add_chainbreakterm();

	bool ddg();

private: ///functions

	void append_residues_nterminally (
		core::Size seq_register,
		core::Size res_pos,
		core::Size stop,
		std::string & nat_seq ,
		core::pose::Pose & target_seeds
	);

	void append_residues_cterminally (
		core::Size seq_register,
		core::Size res_pos,
		core::Size stop,
		std::string & nat_seq ,
		core::pose::Pose & target_seeds
	);

	void grow_from_seeds(
		core::pose::Pose curr_pose, //the pose that gets changed
		std::string sequence,
		protocols::loops::Loops & seeds,
		utility::vector1< Size > cutpoints//,
	);


	void grow_from_verteces(
		core::pose::Pose & curr_pose,
		std::string sequence,
		protocols::loops::Loops & seeds,
		std::set< core::Size > vertex_set
	);

	void setup_cached_observers( core::pose::Pose & pose );

private: /// data

	core::Size extend_nterm;
	core::Size extend_cterm;
	std::string nsequence_;
	std::string csequence_;
	std::string seq_;
	bool all_ala_N;
	bool all_ala_C;
	//bool add_chainbreakterm_;

	bool ddg_;
	/// template pose to derrive a sequence from
	core::pose::PoseOP template_pdb_;
	bool template_presence;
	bool anchor_specified_;
	utility::vector1< core::Size > anchors_;
	core::scoring::ScoreFunctionOP scorefxn_;
	bool output_centroid;
	/// vector of sequence pieces
	utility::vector1< std::string > sequence_chunks_;
	protocols::loops::Loops all_seeds_;
	bool fetch_foldtree;
	///foldtree after growing out peptide pieces
	core::kinematics::FoldTreeOP seed_foldtree_;
	/// pointer for the current pose
	core::pose::PoseOP curr_pose_;
	std::set< core::Size > verteces_;
	//bool add_chainbreakterm_;

};
}//end seeded_abinitio
}//end protocols

#endif
