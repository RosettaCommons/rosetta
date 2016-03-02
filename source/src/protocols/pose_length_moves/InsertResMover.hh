// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/pose_length_moves/InsertResMover.hh
/// @brief

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
using namespace core;
using namespace std;
using utility::vector1;



class InsertResMover : public protocols::moves::Mover {
public:
	InsertResMover();
	numeric::xyzVector<core::Real> center_of_mass(core::pose::Pose const & pose);
	void extendRegion(core::pose::PoseOP poseOP, Size chain_id, Size length);
	virtual void apply( Pose & pose );
	virtual std::string get_name() const;
	moves::MoverOP clone() const { return moves::MoverOP( new InsertResMover( *this ) ); }
	virtual void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	core::pose::PoseOP get_additional_output();
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
};


} // pose_length_moves
} // protocols

#endif
