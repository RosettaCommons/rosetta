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

/// @brief Computes the packing score of a ligand bound to a protein.
/// @param[in] -s <PDBFIL>, where <PDBFIL> is the path to the PDB file
/// 	describing the protein+ligand complex.
/// @param[in] -extra_res_fa <PARFIL>, where <PARFIL> is the path to the params
/// 	file containing the Rosetta description of the ligand.
/// @param[in] -neigh. This flag is optional; if specified, the program computes
/// 	the overall packing score of the ligand and its neighbor residues.
/// @param[in] -nevals <E>, where <E> is the number of evaluations of packing
/// 	over which to compute the final, average score.
/// @details It is assumed that the ligand is the last residue of the pose.
/// @author Andrea Bazzoli (bazzoli@ku.edu)

#include <core/scoring/packstat/compute_sasa.hh>
#include <protocols/neighbor/Neighborhood.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/packstat.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>
#include <string>


using namespace basic::options::OptionKeys;
using basic::options::option;
using utility::vector1;
using core::Size;
using core::Real;

OPT_KEY( Boolean, neigh)
OPT_KEY( Integer, nevals)

static basic::Tracer TR( "apps.pilot.lig_packsco.main" );

using protocols::neighbor::Neighborhood;
using basic::options::OptionKeys::packstat::oversample;
using core::scoring::packstat::compute_residue_packing_scores;

////////////////////////////////////////////////////////////////////////////////
//                                    MAIN                                    //
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char * argv [] )
{
	NEW_OPT( neigh, "Include the ligand's neighbors in the computation of the packing score", false );
	NEW_OPT( nevals, "Number of evaluations of ligand packing", 10);

	devel::init(argc, argv);

	// load pose from pdb file
	core::pose::Pose ps;
	std::string const input_pdb_name( basic::options::start_file() );
	core::import_pose::pose_from_pdb( ps, input_pdb_name );

	// score pose
	core::scoring::ScoreFunctionOP scorefxn(
		core::scoring::getScoreFunction());

	(*scorefxn)(ps);

	// prepare set of residues to compute packing score for
	Size const N = ps.total_residue();
	vector1<Size> ligvec(1, N);
	vector1<Size> packres;
	if(option[neigh]) {
		Neighborhood neigh(ligvec, ps, protocols::neighbor::in_ngbat_sphere);
		packres = neigh.get();
	}
	packres.push_back(N);
	const Size NPCK = packres.size();

	std::cout << "Computing packing score over the following residue set: " << std::endl;
	for(Size i=1; i<=NPCK; ++i)
		std::cout << ps.pdb_info()->chain(packres[i]) << "-" << ps.pdb_info()->number(packres[i]) << std::endl;
	std::cout << std::endl;

	// get average packing score over desired number of evaluations
	int const E = option[nevals];
	std::cout << "Evaluating packing score " << E << " times" << std::endl;
	Real toti=0;
	for(int i=0; i<E; ++i) {

		vector1<Real> ressco = compute_residue_packing_scores(ps, option[oversample]);

		Real totj=0;
		for(Size j=1; j<=NPCK; ++j)
			totj += ressco[packres[j]];
		totj /= NPCK;
		std::cout << "eval " << i << ": " << totj << std::endl;

		toti += totj;
	}
	std::cout << std::endl;

	toti /= E;
	std::cout << "### Average packing score: " << toti << std::endl;
}
