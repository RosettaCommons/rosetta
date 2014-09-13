// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

///
/// @brief Performs constrained minimization on a pose
///
/// @author Andrea Bazzoli (bazzoli@ku.edu)
///
/// @params[in] -cst_coord_sdv X, where X is the standard deviation (Angstrom?)
/// 	of coordinate	displacements (?) (default=0.5)
/// @params[in] -cst_anchor Y, where Y is the pose index of the residue chosen
/// 	to be the constraint anchor (default=1)

#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>

using basic::options::option;
using std::string;
using core::Size;
using core::Real;

OPT_KEY( Real, cst_coord_sdv )
OPT_KEY( Integer, cst_anchor )
using basic::options::OptionKeys::cst_coord_sdv;
using basic::options::OptionKeys::cst_anchor;

int main( int argc, char * argv [] )
{
	NEW_OPT( cst_coord_sdv, "Standard deviation of allowed coordinate displacement (?)", 0.5);
	NEW_OPT( cst_anchor, "Constraint anchor residue index", 1);

	devel::init(argc, argv);

	// load pose from pdb file
	core::pose::Pose ps;
	std::string const input_pdb_name( basic::options::start_file() );
	core::import_pose::pose_from_pdb( ps, input_pdb_name );

	// prepare prefix for output file names
	string::size_type pfxsta = input_pdb_name.find_last_of('/');
	pfxsta = (pfxsta==string::npos) ? 0 : pfxsta+1;

	string::size_type ifnlen = input_pdb_name.length();
	string end4 = input_pdb_name.substr(ifnlen-4, 4);
	string::size_type pfxend = (end4 == ".pdb") ? ifnlen-5 : ifnlen-1;

	string outpfx = input_pdb_name.substr(pfxsta, pfxend-pfxsta+1);

	// create scoring function
	core::scoring::ScoreFunctionOP scorefxn(core::scoring::getScoreFunction());

	// impose coordinate constraint
	Size const NRES = ps.total_residue();
	Real const COORD_SDV = option[cst_coord_sdv];
	Size const ANCHOR = option[cst_anchor];
	for(Size i=1; i<=NRES; ++i) {
		core::conformation::Residue const& rsd(ps.residue(i));
		for(Size ii=1; ii<=rsd.nheavyatoms(); ++ii)
			ps.add_constraint(new core::scoring::constraints::CoordinateConstraint(
				core::id::AtomID(ii,i), core::id::AtomID(1,ANCHOR), rsd.xyz(ii),
					new core::scoring::constraints::HarmonicFunc(0.0, COORD_SDV)));
		}

	scorefxn->set_weight(core::scoring::coordinate_constraint, 1.0);

	// create minimizer
	core::optimization::AtomTreeMinimizer minimizer;
	core::optimization::MinimizerOptions min_options( "dfpmin", 0.00001, true);
	min_options.nblist_auto_update(true);

	core::kinematics::MoveMap mm;
	mm.set_chi( true );
	mm.set_bb( true );
	mm.set_jump( true );

	// minimize pose
	minimizer.run( ps, mm, *scorefxn, min_options);

	// output score
	scorefxn->set_weight(core::scoring::coordinate_constraint, 0);
	Real nrg = (*scorefxn)(ps);
	std::cout << "Energy after minimization: " << nrg << std::endl;

	// output scored pose
	ps.dump_scored_pdb(outpfx + "_" + "0001.pdb", *scorefxn);
}
