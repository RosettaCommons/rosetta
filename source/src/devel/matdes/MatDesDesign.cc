// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file 
/// @brief 
/// @author Neil King ( neilking@uw.edu )
/// @author Javier Castellanos ( javiercv@uw.edu )

// Unit headers
#include <devel/matdes/MatDesDesign.hh>
#include <devel/matdes/MatDesDesignMoverCreator.hh>

// Package headers

// project headers
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/jd2/parser/BluePrint.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

namespace devel {
namespace matdes {

using namespace core;
using namespace utility;

// -------------  Mover Creator -------------
std::string
MatDesDesignCreator::keyname() const
{
	return MatDesDesignCreator::mover_name();
}

MoverOP
MatDesDesignCreator::create_mover() const {
	return new ConstrainedDesign;
}

std::string
MatDesDesignCreator::mover_name()
{
	return "MatDesDesign";
}
// -------------  Mover Creator -------------


MatDesDesign::MatDesDesign() {
} //MatDesDesign()


MatDesDesign::MatDesDesign(const MatDesDesign& rval) {
} // MatDesDesign(const MatDesDesign& rval) 

void 
MatDesDesign::apply(Pose& pose) {
	vector1<Size> design_pos = pick_design_position(pose, contact_dist_, bblock_dist_, probe_radius_);
	if(native_bias_ > 0) 
		add_native_bias_constraints(pose, native_bias, design_pos);
	
  SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);

	// do n_cycles_ cycles of desing
	for(Size iteration = 1; iteration <= n_cycles_; ++iteration) {
    PackerTaskOP task(TaskFactory::create_packer_task(pose));
    if(resfile_ != "") {
    	core::pack::task::parse_resfile(pose,*task, resfile_);
      // Get from the packer task the mutalyze positions
      Size nres_monomer = sym_info->num_independent_residues();
      for(Size i=1; i<=nres_monomer; ++i) {
      	if(! task->residue_task(i).has_behavior("AUTO")) {
      		design_pos.insert(i);
      		task->nonconst_residue_task(i).initialize_from_command_line();
      	}
      }
    }

		// Repack and minimize using score12
		// design
  	protocols::moves::MoverOP packer = new protocols::simple_moves::symmetry::SymPackRotamersMover(scorefxn_design_, task);
  	packer->apply(pose);
		minimize(pose, scorefxn_design_, design_pos, min_bb_, min_sc_, min_rb_);
	}
	repack(pose, scorefxn_minimize_, task);
	minimize(pose, scorefxn_minimize_, design_pos, min_bb_, min_sc_, min_rb_);

	scorefxn_minimize_->score(pose);
	
} // apply


MoverOP MatDesDesign::clone() const {
} // clone


MoverOP
MatDesDesign::fresh_instance() const {
} // fresh_instance

void 
MatDesDesign::minimize(Pose & pose, ScoreFunctionOP scorefxn, std::set<Size> design_pos, bool min_bb, bool min_sc, bool min_rb) {
	// Initialize a MoveMap
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_jump(false);
	movemap->set_bb(false);
	movemap->set_chi(false);
	
	// Set allowable move types at interface positions
	// Currently, only sc moves allowed
	for (std::set<Size>::iterator i = design_pos.begin(); i != design_pos.end(); i++) {
		movemap->set_bb (*i, min_bb);
		movemap->set_chi(*i, min_sc);
	}
	// Make MoveMap symmetric, apply it to minimize the pose
	core::pose::symmetry::make_symmetric_movemap( pose, *movemap );
	protocols::simple_moves::symmetry::SymMinMover m1( movemap, scorefxn, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );
	m1.apply(pose);

	if (min_rb) {
		movemap->set_jump(true);
		core::pose::symmetry::make_symmetric_movemap( pose, *movemap );
 	  protocols::simple_moves::symmetry::SymMinMover m2( movemap, scorefxn, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );
		m2.apply(pose);
	}
}

void
MatDesDesign::repack(Pose & pose, ScoreFunctionOP sf, std::set<Size> design_pos) {

  using namespace core;
  using namespace pack;
  using namespace task;
  using namespace conformation;
  using namespace conformation::symmetry;
  using namespace scoring;
  using namespace chemical;

  // Get the symmetry info and make the packer task
  SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
  PackerTaskOP task( TaskFactory::create_packer_task( pose ));

  // Set which residues can be repacked
  for (Size i=1; i<=pose.n_residue(); i++) {
    if (!sym_info->bb_is_independent(i)) {
      task->nonconst_residue_task(i).prevent_repacking();
    } else if (!design_pos.count(i)) {
      task->nonconst_residue_task(i).prevent_repacking();
    } else {
			vector1<bool> allowed_aas(20, false);
      allowed_aas[pose.residue(i).aa()] = true;
      task->nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_aas);
			task->nonconst_residue_task(i).initialize_from_command_line();
    }
  }

  // Actually repack.
  protocols::moves::MoverOP packer = new protocols::simple_moves::symmetry::SymPackRotamersMover(sf, task);
  packer->apply(pose);

}


void 
MatDesDesign::parse_my_tag( TagCOP const tag,
										 basic::datacache::DataMap & data,
										 Filters_map const &,
										 Movers_map const &,
										 Pose const & ) {
} // parse_my_tag


} // devel
} // matdes
