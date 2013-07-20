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


#ifndef INCLUDED_protocols_hybridization_CartesianSampler_hh
#define INCLUDED_protocols_hybridization_CartesianSampler_hh

#include <protocols/hybridization/InsertChunkMover.hh>
#include <protocols/hybridization/CartesianSampler.fwd.hh>

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
#include <numeric/random/WeightedSampler.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <boost/unordered/unordered_map.hpp>

namespace protocols {
//namespace comparative_modeling {
namespace hybridization {

using namespace core;
using namespace protocols::moves;
using namespace protocols::loops;

class CartesianSampler: public protocols::moves::Mover {
public:
	CartesianSampler();
	CartesianSampler( utility::vector1<core::fragment::FragSetOP> fragments_in );

	// initialize options to defaults
	void init();

	// run the protocol
	void apply(core::pose::Pose & pose);

	// set the centroid scorefunction
	void set_scorefunction(core::scoring::ScoreFunctionOP scorefxn_in) { scorefxn_=scorefxn_in; }

	// set options
	void set_ncycles(core::Size ncycles_in) { ncycles_=ncycles_in; }
	void set_overlap(core::Size overlap_in) { overlap_=overlap_in; }

	std::string	get_name() const { return "CartesianSampler"; }

	void parse_my_tag(
		utility::tag::TagPtr const tag,
		moves::DataMap & data,
		filters::Filters_map const & ,
		moves::Movers_map const & ,
		core::pose::Pose const & pose );

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

protected:
	// apply a sequence fragment
	void apply_frame( core::pose::Pose & pose, core::fragment::Frame &frame );

	//
	void compute_fragment_bias( core::pose::Pose & pose );

	//
	void update_fragment_library_pointers( );

	//
	void apply_constraints( core::pose::Pose & pose );

private:
	// parameters
	core::Size ncycles_, overlap_;

	// fragments
	utility::vector1<core::fragment::FragSetOP> fragments_;
	utility::vector1<boost::unordered_map<core::Size, core::fragment::Frame> > library_;

	// fragment bias
	std::string fragment_bias_strategy_;
	utility::vector1<numeric::random::WeightedSampler> frag_bias_;
	std::set<core::Size> user_pos_;

	// reference model
	core::pose::Pose ref_model_;
	bool input_as_ref_;

	// scorefunctions
	core::scoring::ScoreFunctionOP scorefxn_;
}; //class CartesianSampler

} // hybridize
//} // comparative_modeling
} // protocols

#endif
