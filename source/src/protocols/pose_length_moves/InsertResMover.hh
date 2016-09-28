// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/pose_length_moves/InsertResMover.hh
/// @brief inserts ideal residues into pose. Useful for extending helices
///
/// @author TJ Brunette tjbrunette@gmail.com

#ifndef INCLUDED_protocols_pose_length_moves_InsertResMover_hh
#define INCLUDED_protocols_pose_length_moves_InsertResMover_hh


#include <protocols/moves/Mover.hh>

#include <protocols/pose_length_moves/InsertResMover.fwd.hh>

#include <core/scoring/ScoreFunction.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>

// C++ Headers
#include <string>
#include <map>
// Utility Headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <numeric/xyzVector.hh>


namespace protocols {
namespace pose_length_moves {
// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;
using utility::vector1;



class InsertResMover : public protocols::moves::Mover {
public:
	InsertResMover();
	numeric::xyzVector<core::Real> center_of_mass(core::pose::Pose const & pose);
	void extendRegion(core::pose::PoseOP poseOP, Size chain_id, Size length);
	void apply( Pose & pose ) override;
	std::string get_name() const override;
	moves::MoverOP clone() const override { return moves::MoverOP( new InsertResMover( *this ) ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	core::pose::PoseOP get_additional_output() override;
private:
	vector1<core::pose::PoseOP> posesToOutput_;
	vector1<bool> posesOutputed_;
	std::string chain_;
	Size residue_;
	bool grow_toward_Nterm_;
	bool ideal_;
	Size lowAddRes_;
	Size highAddRes_;
	std::string resType_;
	bool useInputAngles_;
	Real phi_;
	Real psi_;
	Real omega_;
	core::scoring::ScoreFunctionOP scorefxn_;
	Size steal_angles_from_res_;
};


} // pose_length_moves
} // protocols

#endif
