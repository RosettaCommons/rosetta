// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief protocols for folding into density
/// @details
/// @author Frank DiMaio

#include <protocols/electron_density/util.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <core/scoring/dssp/Dssp.hh>


#include <core/id/AtomID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Residue.hh>

// Symmetry
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/pose/symmetry/util.hh>

#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/rigid/RB_geometry.hh>


#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/options/option.hh>

// option key includes
#include <basic/options/keys/edensity.OptionKeys.gen.hh>


#include <basic/Tracer.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>

using basic::T;
using basic::Error;
using basic::Warning;


namespace protocols {
namespace electron_density {

static THREAD_LOCAL basic::Tracer TR( "protocols.electron_density.util" );

using namespace protocols;
using namespace core;


protocols::loops::Loops findLoopFromDensity( core::pose::Pose & pose, core::Real frac, int max_helix_melt, int max_strand_melt ) {
	int nres = pose.size();
	while ( !pose.residue(nres).is_polymer() ) nres--;

	// get dssp parse
	core::scoring::dssp::Dssp secstruct( pose );
	ObjexxFCL::FArray1D< char > dssp_pose( nres );
	secstruct.dssp_reduced (dssp_pose);
	utility::vector1< core::Real > perResCC( nres ), smoothPerResCC( nres );

	// align to map
	protocols::electron_density::SetupForDensityScoringMoverOP dockindens( new protocols::electron_density::SetupForDensityScoringMover );
	dockindens->apply( pose );

	// per-res score
	core::scoring::electron_density::getDensityMap().set_nres( nres );

	//////////////////////
	// now do the matching!
	for ( int r=1; r<=nres; ++r ) {
		// check for missing density ... how? look for atoms > 40A apart
		bool isMissingDens = false;
		for ( int j=1; j<(int)pose.residue(r).natoms(); ++j ) {
			for ( int i=j+1; j<=(int)pose.residue(r).natoms(); ++j ) {
				if ( (pose.residue(r).atom(i).xyz() - pose.residue(r).atom(j).xyz()).length() > 40 ) {
					isMissingDens = true;
					break;
				}
			}
		}

		if ( isMissingDens ) {
			perResCC[r] = 0.0;
		} else {
			perResCC[r] =
				std::max(
				core::scoring::electron_density::getDensityMap().matchRes( r , pose.residue(r), pose, nullptr , false),
				0.001);
		}
	}

	// sort by filtered smoothed CC
	utility::vector1< bool > loopMarker( nres, false );

	// smooth CCs
	smoothPerResCC[1] = 0.67*perResCC[1] + 0.33*perResCC[2];
	for ( int r=2; r<=nres-1; ++r ) {
		smoothPerResCC[r] = 0.25*perResCC[r+1] + 0.5*perResCC[r] + 0.25*perResCC[r];
	}
	smoothPerResCC[nres] = 0.67*perResCC[nres] + 0.33*perResCC[nres-1];

	// missing dens
	for ( int r=1; r<=nres; ++r ) {
		if ( perResCC[r] == 0 ) smoothPerResCC[r] = 0.0;
		if ( r != 1    && perResCC[r-1] == 0 ) smoothPerResCC[r] = 0.0;
		if ( r != nres && perResCC[r+1] == 0 ) smoothPerResCC[r] = 0.0;
	}

	// don't eat into sec struct elements too much
	// filter by setting their CCs artificially high
	for ( int r=1; r<=nres; ++r ) {
		bool in_strand = (max_strand_melt >= 0), in_helix = (max_helix_melt >= 0);
		//bool in_strand = true, in_helix = true;
		for ( int i=std::max(1,r-max_strand_melt), i_end=std::min(nres,r+max_strand_melt); i<= i_end; ++i ) {
			in_strand &= (dssp_pose(i) == 'E');
		}
		for ( int i=std::max(1,r-max_helix_melt), i_end=std::min(nres,r+max_helix_melt); i<= i_end; ++i ) {
			in_helix &= (dssp_pose(i) == 'H');
		}
		if ( in_strand || in_helix ) {
			if ( smoothPerResCC[ r ] > 0.0 ) {  // missing dens
				smoothPerResCC[ r ] = 2;
			}
		}
	}

	utility::vector1< core::Real > sortPerResCC = smoothPerResCC;
	std::sort( sortPerResCC.begin(), sortPerResCC.end() );
	core::Real CCcutoff = sortPerResCC[ (int)std::floor(frac*nres + 0.5) ];
	for ( int r=1; r<=nres; ++r ) {
		loopMarker[r] = (smoothPerResCC[r] < CCcutoff);
	}

	// remove singletons
	for ( int r=2; r<=nres-1; ++r ) {
		if ( loopMarker[r+1] && loopMarker[r-1] ) {
			loopMarker[r] = true;
		}
		if ( !loopMarker[r+1] && !loopMarker[r-1] ) {
			if ( perResCC[r] > 0 ) {
				loopMarker[r] = false;
			} else {
				loopMarker[r+1] = loopMarker[r-1] = true;  // r is missing density; force rebuild
			}
		}
	}

	// fix termini
	// if _ANY_ of the four terminal residues should be rebuilt, build them all
	utility::vector1< int > cuts = pose.fold_tree().cutpoints();
	for ( int i=1; i<=(int)cuts.size(); ++i ) {
		int j = cuts[i];
		if ( j <= nres-4 && (loopMarker[j+4] || loopMarker[j+3] || loopMarker[j+2] || loopMarker[j+1] ) ) {
			loopMarker[j+4] = loopMarker[j+3] = loopMarker[j+2] = loopMarker[j+1] = true;
		}
		if ( j >= 3 && (loopMarker[j] || loopMarker[j-3] || loopMarker[j-2] || loopMarker[j-1]) ) {
			loopMarker[j] = loopMarker[j-3] = loopMarker[j-2] = loopMarker[j-1] = true;
		}
	}

	// finally ... write the loopfile
	protocols::loops::Loops retval;
	core::Size loop_start=0, loop_end=0;
	bool inloop=false;
	for ( int r=1; r<=nres; ++r ) {
		if ( loopMarker[r] && !inloop ) {
			loop_start = r;
			loop_end = r;
			inloop = true;
		} else if ( loopMarker[r] && inloop ) {
			loop_end = r;
		} else if ( !loopMarker[r] && inloop ) {
			//out << "LOOP  " << loop_start << " " << loop_end << "  0 0" << std::endl;
			retval.push_back( protocols::loops::Loop( loop_start, loop_end ) );
			inloop = false;
		}

		// force loop to end at cutpoint
		if ( pose.fold_tree().is_cutpoint( r ) && inloop ) {
			// above cases handle singletons
			retval.push_back( protocols::loops::Loop( loop_start, loop_end ) );
			inloop = false;
		}
	}

	if ( inloop ) {
		retval.push_back( protocols::loops::Loop( loop_start, loop_end ) );
	}

	return retval;
}


core::Real dockPoseIntoMap( core::pose::Pose & pose, std::string align_in /* ="" */ ) {
	using namespace basic::options;

	std::string align = align_in;
	if ( align.length() == 0 ) {
		align = option[ OptionKeys::edensity::realign ]();
	}

	if ( align == "no" ) {
		return 0.0; // do nothing
	}

	// minimization
	core::scoring::ScoreFunctionOP scorefxn_dens( new core::scoring::ScoreFunction() );
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		scorefxn_dens = core::scoring::symmetry::symmetrize_scorefunction( *scorefxn_dens );
	}
	core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn_dens );

	core::scoring::ScoreFunctionOP scorefxn_input = core::scoring::get_score_function();

	core::Real dens_wind = scorefxn_input->get_weight( core::scoring::elec_dens_window );
	core::Real dens_allca = scorefxn_input->get_weight( core::scoring::elec_dens_whole_structure_ca );
	core::Real dens_allatom = scorefxn_input->get_weight( core::scoring::elec_dens_whole_structure_allatom );
	core::Real dens_fast = scorefxn_input->get_weight( core::scoring::elec_dens_fast );

	if ( dens_wind>0 )    scorefxn_dens->set_weight( core::scoring::elec_dens_window, dens_wind );
	if ( dens_allca>0 )   scorefxn_dens->set_weight( core::scoring::elec_dens_whole_structure_ca, dens_allca );
	if ( dens_allatom>0 ) scorefxn_dens->set_weight( core::scoring::elec_dens_whole_structure_allatom, dens_allatom );
	if ( dens_fast>0 )    scorefxn_dens->set_weight( core::scoring::elec_dens_fast, dens_fast );

	// make sure at least 1 density term is on
	core::Real dens_score_sum =
		scorefxn_dens->get_weight( core::scoring::elec_dens_window ) +
		scorefxn_dens->get_weight( core::scoring::elec_dens_whole_structure_ca ) +
		scorefxn_dens->get_weight( core::scoring::elec_dens_whole_structure_allatom ) +
		scorefxn_dens->get_weight( core::scoring::elec_dens_fast );

	if ( align != "no" && dens_score_sum==0 ) {
		scorefxn_dens->set_weight( core::scoring::elec_dens_fast, 10.0 );
	}

	core::Real dens_score = 0.0;

	if ( align == "min" ) {
		bool isSymm = core::pose::symmetry::is_symmetric(pose);

		int root = pose.fold_tree().root();
		utility::vector1< core::kinematics::Edge > root_edges = pose.fold_tree().get_outgoing_edges (root);

		core::kinematics::MoveMapOP rbmm( new core::kinematics::MoveMap );
		rbmm->set_bb( false ); rbmm->set_chi( false );
		TR << "RB minimizing pose into density alongs jump(s)";
		for ( core::Size i=1; i<=root_edges.size(); ++i ) {
			TR << "  " << root_edges[i].label();
			rbmm->set_jump ( root_edges[i].label() , true );
		}
		TR << std::endl;

		if ( isSymm ) {
			core::scoring::ScoreFunctionOP symmscorefxn_dens = core::scoring::symmetry::symmetrize_scorefunction( *scorefxn_dens );
			core::pose::symmetry::make_symmetric_movemap( pose, *rbmm );
			moves::MoverOP min_mover( new simple_moves::symmetry::SymMinMover( rbmm, symmscorefxn_dens,  "lbfgs_armijo_nonmonotone", 1e-5, true ) );
			min_mover->apply( pose );
			symmscorefxn_dens->show( TR, pose ); TR<<std::endl;
			dens_score = (*symmscorefxn_dens)( pose );
		} else {
			moves::MoverOP min_mover( new protocols::simple_moves::MinMover( rbmm, scorefxn_dens, "lbfgs_armijo_nonmonotone", 1e-5, true ) );
			min_mover->apply( pose );
			dens_score = (*scorefxn_dens)( pose );
		}
	}
	return dens_score;
}



}
}

