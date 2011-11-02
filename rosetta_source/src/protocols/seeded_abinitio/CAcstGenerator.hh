// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
///
/// @file protocols/seeded_abinitio/CAcstGenerator.hh
/// @brief Bruno Correia's CA cst generator
/// @author Eva-Maria Strauch (evas01@u.washington.edu)

#ifndef INCLUDED_protocols_seeded_abinitio_CAcstGenerator_hh
#define INCLUDED_protocols_seeded_abinitio_CAcstGenerator_hh

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <utility/string_util.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace seeded_abinitio {

class CAcstGenerator : public protocols::moves::Mover {
 public:
  typedef core::pose::Pose Pose;

  CAcstGenerator();

  // undefined, commenting out to fix PyRosetta build  bool is_seed  ( protocols::loops::Loops & loops, core::Size & residue );

  void apply( core::pose::Pose & pose );

  virtual std::string get_name() const;

  void parse_my_tag( utility::tag::TagPtr const tag,
                     protocols::moves::DataMap &,
                     protocols::filters::Filters_map const &,
                     protocols::moves::Movers_map const &,
                     core::pose::Pose const & );

  protocols::moves::MoverOP clone() const { return( protocols::moves::MoverOP( new CAcstGenerator( *this ) ) ); }
  protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new CAcstGenerator ); }

  virtual ~CAcstGenerator();

 private:
  core::scoring::constraints::ConstraintSetOP ca_cst_;

  ///determines whether constraints for the areas should be which will be "replaced" by the seeds
  bool add_cst_seed_;

  ///container for the cutpoints, since there shouldnt be constraints between two cutpoints
  utility::vector1< core::Size > cut_points_;

  ///stddeviation for the harmonic CA constraints
  core::Real stddev_;

  ///container that has the seed information
  protocols::loops::Loops all_seeds_;

  ///area for which no constraints should be derrived
  protocols::loops::Loops clear_seeds_;

  ///user specified which chain to gather the constraints from
  core::Size from_chain_;

  ///user specified to which chain of the input chain is applied to
  core::Size to_chain_;

  ///user specified a template
  bool template_presence_;

  ///template pdb
  core::pose::PoseOP template_pdb_;

  ///the chain/pose that the user actually wants to read the constraints from
  core::pose::PoseOP curr_pose_;

  ///replace constraints or add onto them
  bool replace_;

  //core::pose::PoseOP target_chain;
  /// sequence separation after which the pair constraints are added
  core::Size seq_separation_;
};

}
}

#endif
