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
#include <protocols/loops/Loops.fwd.hh>
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
class MyLoop{
public:
	unsigned int r1_, r2_;
  	char c1_, c2_;
	unsigned int minn_, maxn_; // potential length range
	MyLoop(){}
  MyLoop(unsigned int r1, char c1, unsigned int minn, unsigned int maxn, unsigned int r2, char c2) 
		: r1_(r1), r2_(r2), c1_(c1), c2_(c2), minn_(minn), maxn_(maxn) {}
};

class LoopHashLoopClosureMover : public protocols::moves::Mover {
public:

  LoopHashLoopClosureMover();
  virtual ~LoopHashLoopClosureMover();
  virtual void apply( core::pose::Pose & pose );
  virtual std::string get_name() const;
  virtual void parse_my_tag( utility::tag::TagCOP const,
			     basic::datacache::DataMap &,
			     protocols::filters::Filters_map const &,
			     protocols::moves::Movers_map const &,
			     core::pose::Pose const &);
  virtual protocols::moves::MoverOP fresh_instance() const;
  protocols::moves::MoverOP clone() const;

private:
	protocols::forge::remodel::RemodelMover_OP remodel_;
	void make_blueprint( const core::pose::Pose& pose, 
											 const std::string & loop_insert_instruction,
											 const std::string & bpname ) const;
	void make_blueprint( const core::pose::Pose& pose, 
											 const std::vector<MyLoop> & loops,
											 const std::string & bpname ) const;
	const std::vector<std::string> tokenize( const std::string& in_str, 
       					const std::string& delimiters ) const;
	const std::vector<MyLoop> make_loops(const std::string & loop_insert_instruction) const;

};
} // loophash_loopclosure
} // devel
#endif
