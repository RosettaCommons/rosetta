// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_moves/StructProfileMover.hh
/// @brief Quickly generates a structure profile

#ifndef INCLUDED_protocols_simple_moves_StructProfileMover_hh
#define INCLUDED_protocols_simple_moves_StructProfileMover_hh

#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/StructProfileMover.fwd.hh>

#include <core/indexed_structure_store/ABEGOHashedFragmentStore.hh>
#include <core/pose/Pose.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>

// C++ Headers
#include <string>
#include <map>
// Utility Headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {
// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;
using utility::vector1;

typedef  core::Real  Probability;

class StructProfileMover : public protocols::moves::Mover {
public:
	StructProfileMover();
	StructProfileMover(Real rmsThreshold,Size consider_topN_frags, Real burialWt, bool only_loops , Real allowed_deviation, Real allowed_deviation_loops, bool eliminate_background, bool outputProfile, bool add_csts_to_pose, bool ignore_terminal_res);
	Size ss_type_convert(char ss_type);
	void read_P_AA_SS_cen6();
	vector1<std::string> get_closest_sequence_at_res(core::pose::Pose const pose, Size res,vector1<Real> cenList);
	vector1<vector1<std::string> > get_closest_sequences(core::pose::Pose const pose,vector1<Real> cenList);
	vector1<vector1<Size> >generate_counts(vector1<vector1<std::string> > top_frag_sequences,core::pose::Pose const pose);
	vector1<vector1<Real> >generate_profile_score(vector1<vector1<Size> > res_per_pos,core::pose::Pose const pose);
	vector1<vector1<Real> >generate_profile_score_wo_background(vector1<vector1<Size> > res_per_pos, vector1<Real> cenList, core::pose::Pose const pose);
	void save_MSAcst_file(vector1<vector1<Real> > profile_score,core::pose::Pose const pose);
	void add_MSAcst_to_pose(vector1<vector1<Real> > profile_score,core::pose::Pose & pose);
	Real get_cen_deviation(std::vector<Real> cenListFrag,vector1<Real> cenListModel);
	vector1< Real> calc_cenlist(core::pose::Pose const pose);
	void apply( Pose & pose ) override;
	std::string get_name() const override;
	moves::MoverOP clone() const override { return moves::MoverOP( new StructProfileMover( *this ) ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
private:
	Real rmsThreshold_;
	std::string aa_order_;
	Size consider_topN_frags_;
	Real burialWt_;
	bool outputProfile_;
	bool add_csts_to_pose_;
	Size cenType_;
	core::indexed_structure_store::ABEGOHashedFragmentStore * ABEGOHashedFragmentStore_;
	typedef utility::vector1< utility::vector1< utility::vector1< Probability > > > Probability_AA_n_n;
	Probability_AA_n_n P_AA_SS_burial_;
	Real allowed_deviation_;
	Real allowed_deviation_loops_;
	bool only_loops_;
	bool eliminate_background_;
	bool ignore_terminal_res_;

};


} // simple_moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_StructProfileMover_hh
