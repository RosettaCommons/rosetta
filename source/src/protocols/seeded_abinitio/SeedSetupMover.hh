// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file protocols/seeded_abinitio/SeedSetupMover
/// @author Eva-Maria Strauch (evas01@u.washington.edu)

#ifndef INCLUDED_protocols_seeded_abinitio_SeedSetupMover_hh
#define INCLUDED_protocols_seeded_abinitio_SeedSetupMover_hh

#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/Loops.hh>



namespace protocols {
namespace seeded_abinitio {

class SeedSetupMover : public protocols::moves::Mover {
 public:
  SeedSetupMover();
  virtual ~SeedSetupMover();
  void apply( core::pose::Pose & pose );

  virtual std::string get_name() const;
  virtual protocols::moves::MoverOP clone() const;
  virtual protocols::moves::MoverOP fresh_instance() const;

  void parse_my_tag( utility::tag::TagCOP tag,
                     basic::datacache::DataMap &,
                     protocols::filters::Filters_map const &,
                     protocols::moves::Movers_map const &,
                     core::pose::Pose const & );

  void clear_task_factory();
  void task_factory( core::pack::task::TaskFactoryOP tf );

  core::pack::task::TaskFactoryOP & task_factory();

 private:
  void set_packerTasks_target_and_seeds (
      core::pose::Pose & pose ,
      protocols::loops::Loops & seeds,
      utility::vector1< core::Size > & designable_residues,
      utility::vector1< core::Size > & norepack_res,
        core::pack::task::TaskFactoryOP & tf );

  /// this taskfactory contains all the information about which residues should be repacked and designed
  /// it assumes the fold pose is the last chain. Seeds are parse at runtime and incorporated into this taskfactory
  core::pack::task::TaskFactoryOP task_factory_;
  core::pack::task::PackerTaskOP task_;

  /// movemap
  core::kinematics::MoveMapOP movemap_;

  /// switches for folding and desing behavior:
  bool repack_target_;
  bool repack_foldpose_;
  bool allow_all_aas_;
  bool design_target_;
  bool design_foldpose_;

  //seeds
  utility::vector1< std::pair < std::string,std::string > > seed_vector_;

  ///designable residues
  std::string design_res_;

  ///dont repack residues
  std::string norepack_res_;

  //seed info
  protocols::loops::Loops all_seeds_;

  /// in case someone wanted to instantly design, but really this is only for debugging purposes in there
  core::scoring::ScoreFunctionOP scorefxn_repack_ ;
  core::scoring::ScoreFunctionOP scorefxn_minimize_;
  bool design_;

  /// let chi angles of chain 1 move
  bool chi_chain1_;
  /// let chi angles of chain 2 move
  bool chi_chain2_;

	///interface chi angles be movable
	bool interface_chi1_;
	bool interface_chi2_;
  core::Real interface_distance_cutoff_;
  /// set helper taskfactory
  bool packtask_;
};

}//end seeded_abinitio
}//end protocols

#endif
