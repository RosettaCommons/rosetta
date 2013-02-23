// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/loophash_loopclosure/LoopHashLoopClosureMover.hh
/// @brief Close loops using loophash library
/// @author Sachko Honda (honda@apl.washington.edu)

#ifndef INCLUDED_devel_loophash_loopclosure_LoopHashLoopClosureMover_HH
#define INCLUDED_devel_loophash_LoopHashLoopClosureMover_HH

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/ddG.fwd.hh>
#include <protocols/forge/remodel/RemodelMover.fwd.hh>

#include <devel/loophash_loopclosure/LoopHashLoopClosureMover.fwd.hh>

#include <utility/tag/Tag.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/vector1.hh>

#include <string>
#include <map>
#include <vector>

namespace devel {
namespace loophash_loopclosure {

class LoopHashLoopClosureMover : public protocols::moves::Mover {
public:

  LoopHashLoopClosureMover();
  virtual ~LoopHashLoopClosureMover();
  virtual void apply( core::pose::Pose & pose );
  virtual std::string get_name() const;
  virtual void parse_my_tag( utility::tag::TagPtr const,
			     protocols::moves::DataMap &,
			     protocols::filters::Filters_map const &,
			     protocols::moves::Movers_map const &,
			     core::pose::Pose const &);
  virtual protocols::moves::MoverOP fresh_instance() const;
  protocols::moves::MoverOP clone() const;

private:
	protocols::forge::remodel::RemodelMover_OP remodel_;
	const std::string make_blueprint(const core::pose::Pose& pose, const std::string & loop_insert_instruction) const;
};
} // loophash_loopclosure
} // devel
#endif
