// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/pose_length_moves/AnalyzeLoopModeling.hh
/// @details connects chains using a very fast RMSD lookback. only works for chains <5 residues. Designed to make loops look within .4 RMSD to naturally occuring loops
/// @author TJ Brunette tjbrunette@gmail.com


#ifndef INCLUDED_protocols_pose_length_moves_AnalyzeLoopModeling_hh
#define INCLUDED_protocols_pose_length_moves_AnalyzeLoopModeling_hh


#include <protocols/moves/Mover.hh>

#include <protocols/pose_length_moves/AnalyzeLoopModeling.fwd.hh>

#include <core/pose/Pose.fwd.hh>

//#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/indexed_structure_store/SSHashedFragmentStore.hh>
// C++ Headers
#include <string>
#include <map>
// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <ctime>

namespace protocols {
namespace pose_length_moves {


class AnalyzeLoopModeling : public protocols::moves::Mover {
public:
	AnalyzeLoopModeling();
	core::Real rmsd_between_coordinates(std::vector< numeric::xyzVector<core::Real> > fragCoordinates,std::vector< numeric::xyzVector<core::Real> > coordinates);
	core::Real get_loop_rmsd(core::pose::Pose native_pose,core::pose::Pose designed_pose, core::Size loopStart, core::Size loopEnd);
	core::Size get_valid_resid(int resid,core::pose::Pose const pose);
	core::Real generate_lookback_rmsd(core::pose::Pose pose, core::Size position);
	void apply( Pose & pose ) override;
	moves::MoverOP clone() const override { return moves::MoverOP( new AnalyzeLoopModeling( *this ) ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();


	//static
	//utility::tag::XMLSchemaRestriction
	//define_size_cspair();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	protocols::loops::Loops get_loops(core::pose::Pose const & pose);
	core::Size loopLengthRangeLow_;
	core::Size loopLengthRangeHigh_;

	core::indexed_structure_store::SSHashedFragmentStore * SSHashedFragmentStore_;
};


} // pose_length_moves
} // protocols

#endif
