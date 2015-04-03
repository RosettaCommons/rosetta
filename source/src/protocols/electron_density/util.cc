// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

static thread_local basic::Tracer TR( "protocols.electron_density.util" );

using namespace protocols;
using namespace core;

// find stratch of residues worst agreeing with patterson map
// for each possible segment, remove and rescore
protocols::loops::Loops findLoopFromPatterson( core::pose::Pose & pose, core::Size N, core::Size nloops, bool allow_termini ) {
	int nres = pose.total_residue();
	while (!pose.residue(nres).is_protein()) nres--;

	runtime_assert( nres > (int)N );

	utility::vector1< core::Real > scores(nres, -1);
	int start_res = allow_termini? 1 : 6;
	int stop_res = allow_termini? nres : nres-5;

	// SLOW
	for (int i=start_res+1; i<= stop_res-(int)N; i++) {
		// check for cutpoints in this segment
		bool contains_cut = false;
		for (int j=i; j<i+(int)N; ++j)
			contains_cut |= pose.fold_tree().is_cutpoint( j );
		if (contains_cut) continue;

		core::scoring::electron_density::getDensityMap().clearMask();

		// mask i->i+N-1
		for (int j=i; j<i+(int)N; ++j)
			core::scoring::electron_density::getDensityMap().maskResidues( j );

		// score
		scores[i] = core::scoring::electron_density::getDensityMap().matchPoseToPatterson( pose, false );
		//std::cerr << "res " << i << " to " << i+N-1 << " --- " << scores[i] << std::endl;
	}

	protocols::loops::Loops retval;
	for (int i=0; i<(int)nloops; ++i) {
		core::Real worst_match = -1;
		int worst_match_idx = -1;
		for (int j=0; j<nres; ++j) {
			if (scores[j] > worst_match) {
				worst_match = scores[j];
				worst_match_idx = j;
			}
		}

		retval.push_back( protocols::loops::Loop( worst_match_idx, worst_match_idx+N-1 ) );
		for (int j=worst_match_idx; j<worst_match_idx+(int)N; ++j)
			scores[j] = -1;
	}

	return ( retval );
}


protocols::loops::Loops findLoopFromDensity( core::pose::Pose & pose, core::Real frac, int max_helix_melt, int max_strand_melt ) {
	int nres = pose.total_residue();
	while (!pose.residue(nres).is_polymer()) nres--;

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
	for (int r=1; r<=nres; ++r) {
		// check for missing density ... how? look for atoms > 40A apart
		bool isMissingDens = false;
		for (int j=1; j<(int)pose.residue(r).natoms(); ++j) {
			for (int i=j+1; j<=(int)pose.residue(r).natoms(); ++j) {
				if ( (pose.residue(r).atom(i).xyz() - pose.residue(r).atom(j).xyz()).length() > 40 ) {
					isMissingDens = true;
					break;
				}
			}
		}

		if (isMissingDens) {
			perResCC[r] = 0.0;
		} else {
			perResCC[r] =
				std::max(
					core::scoring::electron_density::getDensityMap().matchRes( r , pose.residue(r), pose, NULL , false),
					0.001);
		}
	}

	// sort by filtered smoothed CC
	utility::vector1< bool > loopMarker( nres, false );

	// smooth CCs
	smoothPerResCC[1] = 0.67*perResCC[1] + 0.33*perResCC[2];
	for (int r=2; r<=nres-1; ++r) {
		smoothPerResCC[r] = 0.25*perResCC[r+1] + 0.5*perResCC[r] + 0.25*perResCC[r];
	}
	smoothPerResCC[nres] = 0.67*perResCC[nres] + 0.33*perResCC[nres-1];

	// missing dens
	for (int r=1; r<=nres; ++r) {
		if (perResCC[r] == 0) smoothPerResCC[r] = 0.0;
		if (r != 1    && perResCC[r-1] == 0) smoothPerResCC[r] = 0.0;
		if (r != nres && perResCC[r+1] == 0) smoothPerResCC[r] = 0.0;
	}

	// don't eat into sec struct elements too much
	// filter by setting their CCs artificially high
	for (int r=1; r<=nres; ++r) {
		bool in_strand = (max_strand_melt >= 0), in_helix = (max_helix_melt >= 0);
		//bool in_strand = true, in_helix = true;
		for (int i=std::max(1,r-max_strand_melt), i_end=std::min(nres,r+max_strand_melt); i<= i_end; ++i) {
			in_strand &= (dssp_pose(i) == 'E');
		}
		for (int i=std::max(1,r-max_helix_melt), i_end=std::min(nres,r+max_helix_melt); i<= i_end; ++i) {
			in_helix &= (dssp_pose(i) == 'H');
		}
		if ( in_strand || in_helix ) {
			if ( smoothPerResCC[ r ] > 0.0 )  // missing dens
				smoothPerResCC[ r ] = 2;
		}
	}

	utility::vector1< core::Real > sortPerResCC = smoothPerResCC;
	std::sort( sortPerResCC.begin(), sortPerResCC.end() );
	core::Real CCcutoff = sortPerResCC[ (int)std::floor(frac*nres + 0.5) ];
	for (int r=1; r<=nres; ++r) {
		loopMarker[r] = (smoothPerResCC[r] < CCcutoff);
	}

	// remove singletons
	for (int r=2; r<=nres-1; ++r) {
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
	utility::vector1< int >	cuts = pose.fold_tree().cutpoints();
	for (int i=1; i<=(int)cuts.size(); ++i) {
		int j = cuts[i];
		if (j <= nres-4 && (loopMarker[j+4] || loopMarker[j+3] || loopMarker[j+2] || loopMarker[j+1] ) ) {
			loopMarker[j+4] = loopMarker[j+3] = loopMarker[j+2] = loopMarker[j+1] = true;
		}
		if (j >= 3 && (loopMarker[j] || loopMarker[j-3] || loopMarker[j-2] || loopMarker[j-1]) ) {
			loopMarker[j] = loopMarker[j-3] = loopMarker[j-2] = loopMarker[j-1] = true;
		}
	}

	// finally ... write the loopfile
	protocols::loops::Loops retval;
	core::Size loop_start=0, loop_end=0;
	bool inloop=false;
	for (int r=1; r<=nres; ++r) {
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
		if (pose.fold_tree().is_cutpoint( r ) && inloop) {
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


// dock pose into a density map
//    respect recentering flags
//    should we ensure we are set up for density scoring?????
core::Real dockPoseIntoMap( core::pose::Pose & pose, std::string align_in /* ="" */ ) {
	using namespace basic::options;

	std::string align = align_in;
	if (align.length() == 0)
		align = option[ OptionKeys::edensity::realign ]();

	// minimization
	core::scoring::ScoreFunctionOP scorefxn_dens( new core::scoring::ScoreFunction() );
	if ( core::pose::symmetry::is_symmetric(pose) )
		scorefxn_dens = core::scoring::symmetry::symmetrize_scorefunction( *scorefxn_dens );
	core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn_dens );

	core::scoring::ScoreFunctionOP scorefxn_input = core::scoring::get_score_function();
	core::Real dens_patt = scorefxn_input->get_weight( core::scoring::patterson_cc );
	core::Real dens_wind = scorefxn_input->get_weight( core::scoring::elec_dens_window );
	core::Real dens_allca = scorefxn_input->get_weight( core::scoring::elec_dens_whole_structure_ca );
	core::Real dens_allatom = scorefxn_input->get_weight( core::scoring::elec_dens_whole_structure_allatom );
	core::Real dens_fast = scorefxn_input->get_weight( core::scoring::elec_dens_fast );

	if (dens_patt>0)    scorefxn_dens->set_weight( core::scoring::patterson_cc, dens_patt );
	if (dens_wind>0)    scorefxn_dens->set_weight( core::scoring::elec_dens_window, dens_wind );
	if (dens_allca>0)   scorefxn_dens->set_weight( core::scoring::elec_dens_whole_structure_ca, dens_allca );
	if (dens_allatom>0) scorefxn_dens->set_weight( core::scoring::elec_dens_whole_structure_allatom, dens_allatom );
	if (dens_fast>0)    scorefxn_dens->set_weight( core::scoring::elec_dens_fast, dens_fast );

	// make sure at least 1 density term is on
	core::Real dens_score_sum =
		scorefxn_dens->get_weight( core::scoring::patterson_cc ) +
		scorefxn_dens->get_weight( core::scoring::elec_dens_window ) +
		scorefxn_dens->get_weight( core::scoring::elec_dens_whole_structure_ca ) +
		scorefxn_dens->get_weight( core::scoring::elec_dens_whole_structure_allatom );
		scorefxn_dens->get_weight( core::scoring::elec_dens_fast );

	if (align != "no" && align != "random" && dens_score_sum==0 ) {
		scorefxn_dens->set_weight( core::scoring::elec_dens_fast, 1.0 );
	}

	core::Real dens_score = 0.0;

	if (align == "no") {
		//dens_score = (*scorefxn_dens)( pose );
		return dens_score; // do nothing
	}
	if (align == "random") {
		// get jump index of root jump
		int root = pose.fold_tree().root();
		utility::vector1< core::kinematics::Edge > root_edges = pose.fold_tree().get_outgoing_edges (root);
		numeric::xyzMatrix< Real > rot = protocols::geometry::random_reorientation_matrix();
		for (int i=1; i<=(int)root_edges.size(); ++i) {
			core::kinematics::Jump flexible_jump = pose.jump( i );
			flexible_jump.set_rotation( rot * flexible_jump.get_rotation() );
			pose.set_jump( i, flexible_jump );
		}
		dens_score = (*scorefxn_dens)( pose );
		return dens_score;
	}

	/////
	// initial alignment
	if ( align.substr(0,8) == "membrane") {
		// align centers of mass
		fastTransAlignPose( pose );

		// rotation
		fast2DRotAlignPose( pose, option[ OptionKeys::edensity::membrane_axis ]());
	}

	/////
	// minimization
	// special case for symmetric poses
	if ( align.length() >= 3 && align.substr( align.length()-3 ) == "min" ) {
		bool isSymm = core::pose::symmetry::is_symmetric(pose);

		// get jump index of root jump
		int root = pose.fold_tree().root();
		utility::vector1< core::kinematics::Edge > root_edges = pose.fold_tree().get_outgoing_edges (root);

		core::kinematics::MoveMapOP rbmm( new core::kinematics::MoveMap );
		rbmm->set_bb( false ); rbmm->set_chi( false );
		// TODO? a flag which toggles minimization of:
		//  a) all rigid-body DOFs
		//  b) all symmetric DOFs
		// HOWEVER, if this is done then fa_rep / vdw should be turned on
 		TR << "RBminimizing pose into density alongs jump(s)";
 		for (core::Size i=1; i<=root_edges.size(); ++i) {
 			TR << "  " << root_edges[i].label();
 			rbmm->set_jump ( root_edges[i].label() , true );
 		}
 		TR << std::endl;

		if (isSymm) {
			core::scoring::ScoreFunctionOP symmscorefxn_dens = core::scoring::symmetry::symmetrize_scorefunction( *scorefxn_dens );

			core::pose::symmetry::make_symmetric_movemap( pose, *rbmm );
			moves::MoverOP min_mover( new simple_moves::symmetry::SymMinMover( rbmm, symmscorefxn_dens,  "dfpmin_armijo_nonmonotone", 1e-5, true ) );

			bool densInMinimizer = core::scoring::electron_density::getDensityMap().getUseDensityInMinimizer();
			core::scoring::electron_density::getDensityMap().setUseDensityInMinimizer( true );
			min_mover->apply( pose );
			core::scoring::electron_density::getDensityMap().setUseDensityInMinimizer( densInMinimizer );

			symmscorefxn_dens->show( TR, pose ); TR<<std::endl;
			dens_score = (*symmscorefxn_dens)( pose );
		} else {
			moves::MoverOP min_mover( new protocols::simple_moves::MinMover( rbmm, scorefxn_dens, "dfpmin_armijo_nonmonotone", 1e-5, true ) );

			bool densInMinimizer = core::scoring::electron_density::getDensityMap().getUseDensityInMinimizer();
			core::scoring::electron_density::getDensityMap().setUseDensityInMinimizer( true );
			min_mover->apply( pose );
			core::scoring::electron_density::getDensityMap().setUseDensityInMinimizer( densInMinimizer );
			dens_score = (*scorefxn_dens)( pose );
		}
	}
	return dens_score;
}


// align pose CoM to density CoM
core::Real fastTransAlignPose(core::pose::Pose & pose) {
	// align CoM of map and fragment
	int  nres( pose.total_residue() ), nAtms = 0;
  numeric::xyzVector<core::Real> massSum(0,0,0), poseCoM, mapCoM, mapOri;
	for ( int i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		if (rsd.aa() == core::chemical::aa_vrt) continue;

		for ( Size j=1; j<= rsd.nheavyatoms(); ++j ) {
			conformation::Atom const & atom( rsd.atom(j) );
			massSum += atom.xyz();
			nAtms++;
		}
	}
	poseCoM = massSum / (core::Real)nAtms;

	mapCoM = core::scoring::electron_density::getDensityMap().getCoM();
	//mapOri = core::scoring::electron_density::getDensityMap().getOrigin();

	// grid coords -> cart coords
	numeric::xyzVector< core::Real > mapCoMxyz;
	core::scoring::electron_density::getDensityMap().idx2cart( mapCoM , mapCoMxyz );

  numeric::xyzVector<core::Real> translation = mapCoMxyz-poseCoM;
	//TR << "Aligning inital pose to density ... " << std::endl;
	//TR << "   pdb center-of-mass = " << poseCoM << std::endl;
	//TR << "   map center-of-mass = " << mapCoM << std::endl;
	//TR << "   translation (structure -> map) = " << translation  << std::endl;

	// apply translation to each residue in the pose
	for ( int i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		if (rsd.type().aa() == core::chemical::aa_vrt)
			continue;
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			numeric::xyzVector<core::Real> atom_ij = pose.xyz( id::AtomID(j,i) );
			pose.set_xyz( id::AtomID(j,i) ,  (atom_ij+translation) );
		}
	}
	return 0.0;
}

core::Real fast2DRotAlignPose( core::pose::Pose & pose , std::string axis ) {
	if (axis != "X" && axis != "Y" && axis != "Z") {
		TR.Error << "Unrecognized membrane-normal axis '" << axis << "'" << std::endl;
		TR.Error << "Not aligning!" << std::endl;
		return 0.0;
	}

	int  nres( pose.total_residue() ), nAtms = 0;
  numeric::xyzVector<core::Real> massSum(0,0,0), poseCoM, mapCoM, mapOri;
	for ( int i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		if (rsd.aa() == core::chemical::aa_vrt) continue;
		for ( Size j=1; j<= rsd.nheavyatoms(); ++j ) {
			conformation::Atom const & atom( rsd.atom(j) );
			massSum += atom.xyz();
			nAtms++;
		}
	}
	poseCoM = massSum / (core::Real)nAtms;

	numeric::xyzMatrix< core::Real > R
		= core::scoring::electron_density::getDensityMap().rotAlign2DPose( pose, axis );

	// apply translation to each residue in the pose
	for ( int i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		if (rsd.type().aa() == core::chemical::aa_vrt) continue;
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			numeric::xyzVector<core::Real> atom_ij = pose.xyz( id::AtomID(j,i) );
			pose.set_xyz( id::AtomID(j,i) ,  R*(atom_ij-poseCoM)+poseCoM );
		}
	}
	return 0.0;
}


}
}

