// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/*
 * DumpConstraints.cc
 *
 *  Created on: Apr 21, 2009
 *      Author: dgront
 */

#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <core/scoring/func/Func.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/ConstraintIO.hh>

#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


static thread_local basic::Tracer tr( "DisulfideTest" );

void register_options() {

  OPT( in::file::native );
}

int main(int argc, char * argv[]) {
  try {
  using namespace core;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace core::scoring::constraints;

  register_options();
  devel::init(argc, argv);

  pose::Pose init_pose;
  core::import_pose::pose_from_pdb(init_pose, option[in::file::native]());

  std::cerr<<init_pose.sequence()<<"\n";

  ConstraintSetOP cstset;
  if (option[constraints::cst_file].user()) {
    cstset = ConstraintIO::get_instance()->read_constraints(core::scoring::constraints::get_cst_file_option(), new ConstraintSet,
        init_pose);
  }
  // read in constraints
  //  std::string cstfile = option[basic::options::OptionKeys::constraints::cst_file]()[1];
  //  cstset = ConstraintIO::get_instance()->read_constraints(cstfile, new ConstraintSet, init_pose);

  utility::vector1<ConstraintCOP> csts = cstset->get_all_constraints();

  for (utility::vector1<ConstraintCOP>::iterator it = csts.begin(), end = csts.end(); it != end; ++it) {
    (*it)->show_def(std::cout, init_pose);
  }
  } catch ( utility::excn::EXCN_Base const & e ) {
                            std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                }
  return 0;

}
