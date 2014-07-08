// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief A very simple case of SAXS score that reads a PDB and prints the score value
/// @author Dominik Gront

#include <devel/init.hh>
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/rms_util.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/options/keys/casp.OptionKeys.gen.hh>

#include <core/scoring/saxs/SAXSEnergyFA.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>


#include <utility/file/FileName.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>
#include <string>

#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <utility/excn/Exceptions.hh>

basic::Tracer TR("call_saxs_score");

using namespace core;

void register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  OPT(casp::repack);
  OPT(casp::sc_min);
  OPT(in::file::native);
  OPT(in::file::residue_type_set);
  OPT(score::saxs::q_min);
  OPT(score::saxs::q_max);
  OPT(score::saxs::q_step);
  OPT(score::saxs::custom_ff);
  OPT(score::saxs::ref_cen_spectrum);
  OPT(score::saxs::ref_fa_spectrum);
//  OPT(relax::fast);
}

core::Real saxs_energy(core::pose::Pose & a_pose) {

    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
    if ( option[ casp::repack ].user() ) {
	core::pack::task::PackerTaskOP task(
                            core::pack::task::TaskFactory::create_packer_task( a_pose ));
	task->initialize_from_command_line();
        task->restrict_to_repacking();
        core::pack::pack_rotamers(a_pose, (*scorefxn), task);
    }
    if ( option[ casp::sc_min ] ) {
//	std::string const min_type( option[ run::min_type ]() );
	std::string min_type("dfpmin");
        core::kinematics::MoveMap final_mm;
        final_mm.set_chi( true );
        final_mm.set_bb( false );
        core::optimization::AtomTreeMinimizer().run( a_pose, final_mm, *scorefxn,
    		core::optimization::MinimizerOptions( min_type, 0.001, true ) );
    }
    core::scoring::saxs::SAXSEnergyFA en;
    return en.total_energy(a_pose);
}


////////////////////////////////////////////////////////
int main( int argc, char * argv [] ) {
    try {
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    register_options();
    devel::init( argc, argv );

    core::pose::Pose native_pose;
    bool if_native = false;
//------------- Read the native pose  ----------
    if ( option[ in::file::native ].user() )
        core::import_pose::pose_from_pdb( native_pose, option[ in::file::native ]().name() );

//------------- Read the pose for scoring ----------
//    core::chemical::ResidueTypeSetCAP rsd_set_fa = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
    core::pose::Pose fa_pose;
    utility::vector1<utility::file::FileName> s = option[in::file::s]();
    for(Size i=1;i<=s.size();i++) {
	core::import_pose::pose_from_pdb( fa_pose, s[i].name());
	std::cout<< saxs_energy(fa_pose);
	if( if_native )
	    std::cout << ' ' << core::scoring::CA_rmsd( native_pose, fa_pose ) << std::endl;
	else
	    std::cout << std::endl;
    }
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }

    return 0;
}

