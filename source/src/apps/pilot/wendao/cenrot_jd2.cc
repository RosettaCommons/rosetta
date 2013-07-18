// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file fa_to_cenrot_score.cc
/// @brief test centroid rot model
/// @author Yuan Liu


// Core Headers
#include <core/chemical/AA.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>

#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/TenANeighborGraph.hh>

// 
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/dunbrack/cenrot/SingleResidueCenrotLibrary.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
//#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/Minimizer.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/Mover.fwd.hh>

//sampling
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/BBGaussianMover.hh>

#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include <map>

/////////////////////////////////////////////////////////////////
using namespace core;
using namespace core::import_pose;
using namespace core::io::pdb;
using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace core::pose;
using namespace core::chemical;
using namespace core::conformation;
using namespace core::scoring;
using namespace core::kinematics;
using namespace protocols::moves;

using namespace core::pack::task;
using namespace core::pack::dunbrack;
using namespace core::pack::dunbrack::cenrot;

//////////////////////////////////////////////////////////////////
static numeric::random::RandomGenerator RG(62331911);
basic::Tracer TR("pilot.wendao.cenrot");
std::map<std::string, Real> masslst;

utility::vector1<core::Size> nrecovery(20);
//utility::vector1<core::Real> err_buried(20);
utility::vector1<core::Size> n_total(20);

///////////////////////////////////////////////////////////////////
OPT_KEY(Boolean, fit_best_rotamer);
OPT_KEY(Boolean, switch_to_centroid);
OPT_KEY(Boolean, input_cenrot_pdb);
OPT_KEY(Boolean, output_cenrot_intcoord);
OPT_KEY(Boolean, output_bestrot_err);
OPT_KEY(Boolean, output_cenrot_pdb);
OPT_KEY(Boolean, repack_cenrot);
OPT_KEY(Boolean, relax_cenrot);
OPT_KEY(Boolean, opt_after_relax);
OPT_KEY(Integer, repack_buried_cutoff);
OPT_KEY(Real, repack_bfactor_cutoff);
OPT_KEY(Integer, repack_ncycle);
OPT_KEY(String, output_cenrot_score);
OPT_KEY(String, output_cenrot_dir);
OPT_KEY(String, output_cenrot_prefix);

OPT_KEY(Real, relax_temp);
OPT_KEY(Integer, relax_step_per_cycle);
OPT_KEY(Integer, relax_cycle_number);

int main( int argc, char * argv [] ) {
	NEW_OPT(fit_best_rotamer, "fit the exact centroid to the closest in lib", false);
	NEW_OPT(switch_to_centroid, "switch the fa pdb to the old centroid", false);
	NEW_OPT(output_cenrot_intcoord, "output the internal coordinate, for building lib", false);
	NEW_OPT(output_bestrot_err, "output the distance between the native rot and best fitted one", false);
	NEW_OPT(input_cenrot_pdb, "input centroid pdbs for scoring", false);
	NEW_OPT(output_cenrot_pdb, "output centroid pdbs for building database", false);
	NEW_OPT(output_cenrot_score, "score the centroid pdbs", "cenrot_score.out");
	NEW_OPT(output_cenrot_dir, "dir for output centroid pdbs", ".");
	NEW_OPT(output_cenrot_prefix, "prefix for pdbs", "idealized_");
	NEW_OPT(repack_cenrot, "repack the centroid rotamer model", false);
	NEW_OPT(repack_bfactor_cutoff, "count repack side-chain with Bfactor lower than default 100", 100.0);
	NEW_OPT(repack_buried_cutoff, "count repack side-chain with buried cutoff", 0);
	NEW_OPT(repack_ncycle, "how many times to repack", 1);

	NEW_OPT(relax_cenrot, "relax the centroid rotamer model", false);
	NEW_OPT(opt_after_relax, "opt after relax", false);
	NEW_OPT(relax_temp, "temp", 1.0);
	NEW_OPT(relax_step_per_cycle, "step", 100);
	NEW_OPT(relax_cycle_number, "cycle", 1);

	devel::init(argc, argv);

	//element mass
	masslst["C"]=12.0107;
	masslst["O"]=15.9994;
	masslst["S"]=32.066;
	masslst["N"]=14.00674;
	masslst["H"]=1.00794;

	return 0;
}

