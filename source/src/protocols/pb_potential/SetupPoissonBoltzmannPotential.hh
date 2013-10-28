// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/SetupAPBSMover.hh
/// @brief Setup for Poisson-Boltzmann energy term use in scorefunction
/// This mover assumes you have a path to the APBS (Adaptive Poisson-Boltzmann Solver) program
/// in the system PATH.  The executable name must be "apbs", case-sensitive.
/// @author Sachko Honda (honda@apl.washington.edu)

#ifndef INCLUDED_protocols_pb_potential_SetupPoissonBoltzmannPotential_hh
#define INCLUDED_protocols_pb_potential_SetupPoissonBoltzmannPotential_hh

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
#include <protocols/pb_potential/SetupPoissonBoltzmannPotential.fwd.hh>

#include <utility/tag/Tag.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/vector1.hh>

#include <string>
#include <map>
#include <vector>

namespace protocols {
namespace pb_potential {

class SetupPoissonBoltzmannPotential : public moves::Mover {
public:

  static const std::string APBS_CONFIG_EXT;
  static const std::string APBS_PQR_EXT;
  static const std::string APBS_DX_EXT;
  static const std::string DEFAULT_STATE;
	static const std::string DEFAULT_APBS_PATH;

  SetupPoissonBoltzmannPotential();
  virtual ~SetupPoissonBoltzmannPotential();
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
	protocols::simple_moves::ddGOP ddg_;  // used as pre-scoring & caching state data
};
}
}
#endif
