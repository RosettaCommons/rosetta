// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Yifan Song
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_comparative_modeling_hybridize_CartesianHybridize_hh
#define INCLUDED_protocols_comparative_modeling_hybridize_CartesianHybridize_hh

#include <protocols/comparative_modeling/hybridize/InsertChunkMover.hh>
#include <protocols/comparative_modeling/hybridize/CartesianHybridize.fwd.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/util/kinematics_util.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>

#include <core/scoring/ScoreFunction.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <protocols/moves/Mover.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <boost/unordered/unordered_map.hpp>

namespace protocols {
namespace comparative_modeling {
namespace hybridize {

using namespace core;
using namespace protocols::moves;
using namespace protocols::loops;
	
class CartesianHybridize: public protocols::moves::Mover {
public:
	CartesianHybridize() : ncycles_(400) {};

	CartesianHybridize(
		utility::vector1 < core::pose::PoseOP > const & templates_in,
		utility::vector1 < core::Real > const & template_wts_in,
		utility::vector1 < protocols::loops::Loops > const & template_chunks_in, 
		utility::vector1 < protocols::loops::Loops > const & template_contigs_in,
		core::fragment::FragSetOP fragments9_in );

	// run the protocol
	void apply(core::pose::Pose & pose);

	// set the centroid scorefunction
	void set_scorefunction(core::scoring::ScoreFunctionOP scorefxn_in);

	//
	std::string	get_name() const { return "CartesianHybridize"; }

protected:
	// apply a homologue fragment
	void apply_frag( core::pose::Pose &pose, core::pose::Pose &templ, protocols::loops::Loop &frag, bool superpose=true);

	// apply a sequence fragment
	void apply_frame( core::pose::Pose & pose, core::fragment::Frame &frame );

private:
	// parameters
	core::Real increase_cycles_;
	core::Size ncycles_;

	// fragments
	utility::vector1 < core::pose::PoseOP > templates_;
	utility::vector1 < core::Real > template_wts_;
	utility::vector1 < protocols::loops::Loops > template_contigs_;
	core::fragment::FragSetOP fragments9_;
	boost::unordered_map<core::Size, core::fragment::Frame> library_;

	// scorefunctions
	core::scoring::ScoreFunctionOP lowres_scorefxn_, min_scorefxn_, bonds_scorefxn_, nocst_scorefxn_;
}; //class CartesianHybridize
	
} // hybridize 
} // comparative_modeling 
} // protocols

#endif
