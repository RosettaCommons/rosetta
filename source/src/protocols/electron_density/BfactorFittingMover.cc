// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Set up morphing with electron density map
/// @author Yifan Song

#include <protocols/electron_density/BfactorFittingMover.hh>
#include <protocols/electron_density/BfactorFittingMoverCreator.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/electron_density/xray_scattering.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/scoring/cryst/util.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/kinematics/MoveMap.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/cryst.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>


#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>

static basic::Tracer TR( "protocols.electron_density.BfactorFittingMover" );

namespace protocols {
namespace electron_density {

using core::scoring::electron_density::poseCoords;
using core::scoring::electron_density::poseCoord;

/// helper function
void
symmetrizeBfactors( core::pose::Pose & pose ) {

  if ( !core::pose::symmetry::is_symmetric(pose) )
		return;

	core::conformation::symmetry::SymmetryInfoOP symm_info;
  core::conformation::symmetry::SymmetricConformation & SymmConf (
		dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
	symm_info = SymmConf.Symmetry_Info();

	for (core::Size resid=1; resid<=pose.total_residue(); ++resid) {
		if ( symm_info && !symm_info->bb_is_independent( resid ) ) continue;
		core::conformation::Residue const & rsd( pose.residue(resid) );
		if ( rsd.aa() == core::chemical::aa_vrt ) continue;

		for ( core::conformation::symmetry::SymmetryInfo::Clones::const_iterator pos=symm_info->bb_clones( resid ).begin(),
		      epos=symm_info->bb_clones( resid ).end(); pos != epos; ++pos ) {
			core::Size natoms = rsd.natoms();
			for (core::Size atmid = 1; atmid <= natoms; ++atmid) {
				core::Real B = pose.pdb_info()->temperature( resid, atmid );
				pose.pdb_info()->temperature( *pos, atmid , B );
			}
		}
	}
}


///
BfactorMultifunc::BfactorMultifunc(
		core::pose::Pose & pose_in,
		core::Real wt_adp,
		core::Real wt_dens,
		core::Real rmax,
		core::Real radius_exp,
		core::Real scorescale,
		bool exact, bool verbose, bool deriv_check ) :
	pose_(pose_in), wt_adp_(wt_adp), wt_dens_(wt_dens), rmax_(rmax), radius_exp_(radius_exp), scorescale_( scorescale ),
	exact_( exact ), verbose_(verbose), deriv_check_( deriv_check )
{
	B_EPS=0.0001;

	// map AtomIDs <-> indices
	core::pose::initialize_atomid_map( atom_indices_, pose_in );
	core::conformation::symmetry::SymmetryInfoOP symm_info;
  if ( core::pose::symmetry::is_symmetric(pose_in) ) {
    core::conformation::symmetry::SymmetricConformation & SymmConf (
      dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose_in.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
  }

	for (Size resid=1; resid<=pose_in.total_residue(); ++resid) {
		if ( symm_info && !symm_info->bb_is_independent( resid ) ) continue;
		core::conformation::Residue const & rsd( pose_in.residue(resid) );
		if ( rsd.aa() == core::chemical::aa_vrt ) continue;

		core::Size natoms = rsd.nheavyatoms();
		for (core::Size atmid = 1; atmid <= natoms; ++atmid) {
			if (rsd.atom_type( atmid ).is_virtual()) continue;
			core::id::AtomID a_i( atmid, resid );
			moving_atoms_.push_back(a_i);
			atom_indices_[ a_i ] = moving_atoms_.size();
		}
	}
}


///
///  pose with B factors --> multivec
void
BfactorMultifunc::poseBfacts2multivec( core::pose::Pose & pose, core::optimization::Multivec &y ) const {
	runtime_assert( pose.pdb_info() );  // must exist

	core::Size natoms = moving_atoms_.size();
	y.resize(natoms);

	for (core::Size i = 1; i <= natoms; ++i) {
		core::id::AtomID a_i = moving_atoms_[i];
		core::Real B = pose.pdb_info()->temperature( a_i.rsd(), a_i.atomno() );
		y[i] = log( B );
	}
}



///
///   multivec ---> pose B factors
void
BfactorMultifunc::multivec2poseBfacts( core::optimization::Multivec const &y, core::pose::Pose & pose ) const {
	runtime_assert( pose.pdb_info() );  // must exist

	core::Size natoms = moving_atoms_.size();
	runtime_assert( y.size() == natoms );

	for (core::Size i = 1; i <= natoms; ++i) {
		core::id::AtomID a_i = moving_atoms_[i];
		// fpd! upper limit must match what is in xray_scattering.hh
		//  THIS SHOULD BE SELECTABLE
		pose.pdb_info()->temperature( a_i.rsd(), a_i.atomno(), std::min(std::max( exp( y[i] ) ,0.0),600.0) );
	}
}

///
///  B factor optimization multifunc
core::Real
BfactorMultifunc::operator ()( core::optimization::Multivec const & vars ) const {
	// [[1]] density score
	core::pose::Pose pose_copy = pose_;  // inefficient?
	multivec2poseBfacts( vars, pose_copy );

	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose_copy) ) {
		core::conformation::symmetry::SymmetricConformation const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose_copy.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}

	core::Real dens_score = 0;
	core::Size natom=0;
	if (!exact_) {
		for ( Size i = 1; i <= pose_copy.total_residue(); ++i ) {
			if ( symm_info && !symm_info->bb_is_independent( i ) ) continue;
			core::conformation::Residue const & rsd ( pose_copy.residue(i) );
			if ( rsd.aa() == core::chemical::aa_vrt ) continue;
			dens_score -= core::scoring::electron_density::getDensityMap().matchResFast( i, rsd, pose_copy, NULL );
		}
	} else {
		core::scoring::electron_density::poseCoords litePose;
		for (uint i = 1; i <= pose_copy.total_residue(); ++i) {
			if (symm_info && !symm_info->bb_is_independent( i ) ) continue;
			core::conformation::Residue const & rsd_i ( pose_copy.residue(i) );
			if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;

			core::Size natoms = rsd_i.nheavyatoms();
			for (uint j = 1; j <= natoms; ++j) {
				if ( rsd_i.atom_type(j).is_virtual() ) continue;
				core::conformation::Atom const &atom_j( rsd_i.atom(j) );
				core::chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );

				core::scoring::electron_density::poseCoord coord_j;
				coord_j.x_ = rsd_i.xyz( j );
				coord_j.B_ = pose_copy.pdb_info()->temperature( i, j );
				coord_j.elt_ = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
				litePose.push_back( coord_j );
			}
		}
		core::scoring::electron_density::getDensityMap().calcRhoC( litePose );
		dens_score -= (moving_atoms_.size()/10.0) * core::scoring::electron_density::getDensityMap().getRSCC(litePose);
		natom = litePose.size();
	}


	// [[2]]  constraint score
	core::Real cst_score=0;
	core::Size nedge=0;
	if (wt_adp_ != 0) {
		core::scoring::EnergyGraph const & energy_graph( pose_copy.energies().energy_graph() );
		Size const nres( energy_graph.num_nodes() );

		for ( Size i = 1; i <= pose_copy.total_residue(); ++i ) {
			if ( symm_info && !symm_info->bb_is_independent( i ) ) continue;
			core::conformation::Residue const & rsd1 ( pose_copy.residue(i) );
			if ( rsd1.aa() == core::chemical::aa_vrt ) continue;

			core::Size natoms = rsd1.nheavyatoms();
			for (uint k = 1; k <= natoms; ++k) {
				core::Real B_k = pose_copy.pdb_info()->temperature( i, k );
				if (atom_indices_[core::id::AtomID(k,i)] == 0) continue;

				for (uint l = k + 1; l <= natoms; ++l) {
					core::Real B_l = pose_copy.pdb_info()->temperature( i, l );
					if (atom_indices_[core::id::AtomID(l,i)] == 0) continue;

					core::Real dist_kl = (rsd1.atom( k ).xyz() - rsd1.atom( l ).xyz()).length();
					if (radius_exp_ != 1.0) dist_kl = std::pow( dist_kl, radius_exp_ );

					if (dist_kl < rmax_) {
						cst_score += (B_k-B_l)*(B_k-B_l) / ( dist_kl*(B_k+B_l+B_EPS) );
						nedge++;
					}
				}
			}

			for ( core::graph::Graph::EdgeListConstIter
					iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
					irue = energy_graph.get_node(i)->const_upper_edge_list_end();
					iru != irue; ++iru ) {
				core::scoring::EnergyEdge const * edge( static_cast< core::scoring::EnergyEdge const *> (*iru) );
				core::Size const j( edge->get_second_node_ind() );
				core::conformation::Residue const & rsd2 ( pose_copy.residue(j) );
				if ( rsd2.aa() == core::chemical::aa_vrt ) continue;
				core::Size natoms2 = rsd2.nheavyatoms();

				// inter-atom
				for (uint k = 1; k <= natoms; ++k) {
					core::Real B_k = pose_copy.pdb_info()->temperature( i, k );
					if (atom_indices_[core::id::AtomID(k,i)] == 0) continue;
					for (uint l = 1; l <= natoms2; ++l) {
						core::Real B_l = pose_copy.pdb_info()->temperature( j, l );
						if (atom_indices_[core::id::AtomID(l,j)] == 0) continue;

						core::Real dist_kl = (rsd1.atom( k ).xyz() - rsd2.atom( l ).xyz()).length();
						if (radius_exp_ != 1.0) dist_kl = std::pow( dist_kl, radius_exp_ );

						if (dist_kl < rmax_) {
							cst_score += (B_k-B_l)*(B_k-B_l) / ( dist_kl*(B_k+B_l+B_EPS) );
							nedge++;
						}
					}
				}
			}
		}
	}

	if (verbose_) {
		//core::optimization::Multivec varsCopy = vars;
		//std::cerr << "[score] evaluated " << natom << " atoms and " << nedge << " edges" << std::endl;
		std::cerr << "[score] score = " << scorescale_* (wt_dens_*dens_score + wt_adp_*cst_score) << " [ " << dens_score << " , " << cst_score << " ] " << std::endl;
	}

	return ( scorescale_* (wt_dens_*dens_score + wt_adp_*cst_score) );
}

void
BfactorMultifunc::dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const {
	core::pose::Pose pose_copy = pose_;

	multivec2poseBfacts( vars, pose_copy );
	dE_dvars.resize( vars.size(), 0.0 );

	if (!exact_) {
		// [[1]] density deriv
		for (core::Size i = 1; i <= vars.size(); ++i) {
			core::id::AtomID a_i = moving_atoms_[i];
			core::Real dCCdb = core::scoring::electron_density::getDensityMap().dCCdB_fastRes(
				a_i.atomno(), a_i.rsd(), pose_copy.residue( a_i.rsd() ), pose_copy );
			dE_dvars[ i ] = -scorescale_*wt_dens_*dCCdb * exp(vars[i]);
		}
	} else {
		core::scoring::electron_density::getDensityMap().dCCdBs( pose_copy, dE_dvars );
		for (core::Size i = 1; i <= vars.size(); ++i)
			dE_dvars[i] *= -(moving_atoms_.size()/10.0)*scorescale_*wt_dens_ * exp(vars[i]);
	}

	// [[2]]  constraint deriv
	core::Size nedge=0;
	if (wt_adp_ != 0) {
		core::conformation::symmetry::SymmetryInfoOP symm_info;
		if ( core::pose::symmetry::is_symmetric(pose_copy) ) {
			core::conformation::symmetry::SymmetricConformation & SymmConf (
				dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose_copy.conformation()) );
			symm_info = SymmConf.Symmetry_Info();
		}

		core::scoring::EnergyGraph const & energy_graph( pose_copy.energies().energy_graph() );
		core::Size const nres( energy_graph.num_nodes() );

		for ( core::Size i = 1; i <= pose_copy.total_residue(); ++i ) {
			if ( symm_info && !symm_info->bb_is_independent( i ) ) continue;
			core::conformation::Residue const & rsd1 ( pose_copy.residue(i) );
			if ( rsd1.aa() == core::chemical::aa_vrt ) continue;

			// intra-res
			core::Size natoms = rsd1.nheavyatoms();
			for (uint k = 1; k <= natoms; ++k) {
				core::Real B_k = pose_copy.pdb_info()->temperature( i, k );
				if (atom_indices_[core::id::AtomID(k,i)] == 0) continue;
				for (uint l = k + 1; l <= natoms; ++l) {
					core::Real B_l = pose_copy.pdb_info()->temperature( i, l );
					if (atom_indices_[core::id::AtomID(l,i)] == 0) continue;

					core::Real dist_kl = (rsd1.atom( k ).xyz() - rsd1.atom( l ).xyz()).length();
					if (radius_exp_ != 1.0) dist_kl = std::pow( dist_kl, radius_exp_ );

					if (dist_kl < rmax_) {
						core::Size atomK = atom_indices_[ core::id::AtomID( k,i ) ];
						core::Size atomL = atom_indices_[ core::id::AtomID( l,i ) ];

						core::Real cst_kl = (B_k-B_l)*(B_k-B_l) / ( dist_kl*(B_k+B_l+B_EPS)*(B_k+B_l+B_EPS) );
						core::Real dcst_dbk = -cst_kl + 2*(B_k-B_l) / ( dist_kl*(B_k+B_l+B_EPS) );
						core::Real dcst_dbl = -cst_kl - 2*(B_k-B_l) / ( dist_kl*(B_k+B_l+B_EPS) );

						dE_dvars[ atomK ] += scorescale_*wt_adp_*dcst_dbk * exp(vars[atomK]);
						dE_dvars[ atomL ] += scorescale_*wt_adp_*dcst_dbl * exp(vars[atomL]);
						nedge++;
					}
				}
			}

			for ( core::graph::Graph::EdgeListConstIter
					iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
					irue = energy_graph.get_node(i)->const_upper_edge_list_end();
					iru != irue; ++iru ) {
				core::scoring::EnergyEdge const * edge( static_cast< core::scoring::EnergyEdge const *> (*iru) );
				core::Size const j( edge->get_second_node_ind() );
				core::conformation::Residue const & rsd2 ( pose_copy.residue(j) );
				if ( rsd2.aa() == core::chemical::aa_vrt ) continue;
				core::Size natoms2 = rsd2.nheavyatoms();

				// inter-res
				for (uint k = 1; k <= natoms; ++k) {
					core::Real B_k = pose_copy.pdb_info()->temperature( i, k );
					if (atom_indices_[core::id::AtomID(k,i)] == 0) continue;

					for (uint l = 1; l <= natoms2; ++l) {
						core::Real B_l = pose_copy.pdb_info()->temperature( j, l );
						if (atom_indices_[core::id::AtomID(l,j)] == 0) continue;

						core::Real dist_kl = (rsd1.atom( k ).xyz() - rsd2.atom( l ).xyz()).length();
						if (radius_exp_ != 1.0) dist_kl = std::pow( dist_kl, radius_exp_ );

						if (dist_kl < rmax_) {
							core::Size atomK = atom_indices_[ core::id::AtomID( k,i ) ];
							core::Size atomL = atom_indices_[ core::id::AtomID( l,j ) ];

							core::Real cst_kl = (B_k-B_l)*(B_k-B_l) / ( dist_kl*(B_k+B_l+B_EPS)*(B_k+B_l+B_EPS) );
							core::Real dcst_dbk = -cst_kl + 2*(B_k-B_l) / ( dist_kl*(B_k+B_l+B_EPS) );
							core::Real dcst_dbl = -cst_kl - 2*(B_k-B_l) / ( dist_kl*(B_k+B_l+B_EPS) );

							dE_dvars[ atomK ] += scorescale_*wt_adp_*dcst_dbk * exp(vars[atomK]);
							dE_dvars[ atomL ] += scorescale_*wt_adp_*dcst_dbl * exp(vars[atomL]);
							nedge++;
						}
					}
				}
			}
		}
	}

	if (verbose_) {
		core::Real score = (*this)(vars);
		//std::cerr << "[deriv] evaluated " << vars.size() << " atoms and " << nedge << " edges" << std::endl;
		std::cerr << "[deriv] score = " << score << std::endl;
	}
}

void
BfactorMultifunc::dump( core::optimization::Multivec const & x1, core::optimization::Multivec const & x2 ) const {
	// debug
	if (deriv_check_) {
 		 core::optimization::Multivec varsCopy = x1;
 		 core::optimization::Multivec dE_dvars;
 		 this->dfunc(x1,dE_dvars);

 		 core::Real score = (*this)(x1);
 		 std::cerr << "[deriv] score = " << score << std::endl;
 		 for (uint i = 1; i <= x1.size(); ++i) {
 		 	varsCopy[i]+=0.0001;
 		 	core::Real scorep = (*this)(varsCopy);
 		 	varsCopy[i]-=0.0002;
 		 	core::Real scoren = (*this)(varsCopy);
 		 	varsCopy[i]+=0.0001;
 		 	std::cerr << "[deriv] x:" << i << "  B: " << x1[i] << "  A: " << dE_dvars[i] << "  N: "
 				<< (scorep-scoren)/0.0002 << "    r: " << 0.0002*dE_dvars[i]/(scorep-score) << std::endl;
 		}
	}
}



BfactorFittingMover::BfactorFittingMover() : moves::Mover() {
	init();
}

void BfactorFittingMover::init() {
	wt_adp_ = 0.001;
	wt_dens_ = 1.0;
	rmax_ = 5.0;
	radius_exp_ = 1;
	max_iter_ = 200;
	scorescale_ = 5000.0;
	//minimizer_ = "dfpmin_armijo";
	minimizer_ = "lbfgs_armijo_nonmonotone_atol";
	init_ = true;
	exact_ = false;
	opt_to_fsc_ = false;
	verbose_ = false;
	deriv_check_ = false;

	weights_to_scan_.clear();
	weights_to_scan_.push_back(10);
	weights_to_scan_.push_back(20);
	weights_to_scan_.push_back(40);
	weights_to_scan_.push_back(60);
	weights_to_scan_.push_back(80);
	weights_to_scan_.push_back(100);
	weights_to_scan_.push_back(120);
	weights_to_scan_.push_back(140);
	weights_to_scan_.push_back(160);
	weights_to_scan_.push_back(180);
	weights_to_scan_.push_back(200);
	weights_to_scan_.push_back(300);
}

void BfactorFittingMover::apply(core::pose::Pose & pose) {
	if (!basic::options::option[ basic::options::OptionKeys::cryst::crystal_refine ]()) {
		utility_exit_with_message( "Bfactor optimization _must_ be run with the -cryst::crystal_refine option" );
	}

	core::optimization::MinimizerOptions options( minimizer_, 1e-4, true, false, false );
	options.max_iter(max_iter_);

	// set up optimizer
	BfactorMultifunc f_b( pose, wt_adp_, wt_dens_, rmax_, radius_exp_, scorescale_, exact_, verbose_, deriv_check_);
	core::optimization::Minimizer minimizer( f_b, options );

	// make sure interaction graphs are set up
	core::scoring::ScoreFunctionOP sf;
	if ( core::pose::symmetry::is_symmetric(pose) )
		sf = new core::scoring::symmetry::SymmetricScoreFunction();
	else
		sf = new core::scoring::ScoreFunction();
	sf->set_weight( core::scoring::fa_rep, 1.0 );
	(*sf)(pose);

	// ... and run
	core::optimization::Multivec y;
	f_b.poseBfacts2multivec( pose, y );

	// if b=0 set b to some small value
	for (uint i = 1; i <= y.size(); ++i) if (y[i] < 1e-6) y[i] = 5;

	core::scoring::electron_density::getDensityMap().rescale_fastscoring_temp_bins( pose, init_ );

	// optionally grid search for best starting value
	core::Size nvalues = weights_to_scan_.size();
	core::Real bestVal, bestScore;

	if (opt_to_fsc_) {
		// pose->poseCoords
		poseCoords litePose;
		for (uint i = 1; i <= pose.total_residue(); ++i) {
			core::conformation::Residue const & rsd_i ( pose.residue(i) );
			if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;
			core::Size natoms = rsd_i.nheavyatoms();
			for (uint j = 1; j <= natoms; ++j) {
				core::conformation::Atom const &atom_j( rsd_i.atom(j) );
				core::chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );
				poseCoord coord_j;
				coord_j.x_ = rsd_i.xyz( j );
				coord_j.elt_ = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
				litePose.push_back( coord_j );
			}
		}

		for (core::Size j=1; j<=nvalues; ++j) {
			for (uint i = 1; i <= litePose.size(); ++i) litePose[i].B_ = weights_to_scan_[j];

			utility::vector1< core::Real > modelmapFSC(nresbins_,1.0), modelmapError(nresbins_,1.0);
			core::scoring::electron_density::getDensityMap().getFSC(
				litePose, nresbins_, 1.0/res_low_, 1.0/res_high_, modelmapFSC );
			core::Real thisScore = 0;
			for (Size i=1; i<=modelmapFSC.size(); ++i) thisScore+=modelmapFSC[i];
			thisScore /= modelmapFSC.size();

			TR << "B=" << weights_to_scan_[j] << ": FSC=" << thisScore << std::endl;

			if (thisScore>bestScore || j==1) {
				bestScore = thisScore;
				bestVal = weights_to_scan_[j];
			}
		}

		for (uint i = 1; i <= y.size(); ++i) y[i] = bestVal;
	} else if (init_) {
		for (core::Size j=1; j<=nvalues; ++j) {
			for (uint i = 1; i <= y.size(); ++i)
				y[i] = log( weights_to_scan_[j] );
			core::Real thisScore = f_b(y);
			TR << "B=" << weights_to_scan_[j] << ": score=" << thisScore << std::endl;
			if (thisScore<bestScore || j==1) {
				bestScore = thisScore;
				bestVal = weights_to_scan_[j];
			}
		}
		for (uint i = 1; i <= y.size(); ++i) y[i] = log(bestVal);
	}

	if (!opt_to_fsc_) {
		minimizer.run( y );
	}

	f_b.multivec2poseBfacts( y, pose );
	core::scoring::cryst::fix_bfactorsH( pose ); // fix hydrogen B's to riding atoms
	symmetrizeBfactors( pose );
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
BfactorFittingMover::parse_my_tag(
								   TagCOP const tag,
								   basic::datacache::DataMap & datamap,
								   Filters_map const & filters,
								   moves::Movers_map const & movers,
								   Pose const & pose
							   )
{
	if ( tag->hasOption("wt_adp") ) {
		wt_adp_ = tag->getOption<core::Real>("wt_adp");
	}
	if ( tag->hasOption("wt_dens") ) {
		wt_dens_ = tag->getOption<core::Real>("wt_dens");
	}
	if ( tag->hasOption("rmax") ) {
		rmax_ = tag->getOption<core::Real>("rmax");
	}
	if ( tag->hasOption("scale") ) {
		scorescale_ = tag->getOption<core::Real>("scale");
	}
	if ( tag->hasOption("radius_exp") ) {
		radius_exp_ = tag->getOption<core::Size>("radius_exp");
	}
	if ( tag->hasOption("max_iter") ) {
		max_iter_ = tag->getOption<core::Size>("max_iter");
	}
	if ( tag->hasOption("minimizer") ) {
		minimizer_ = tag->getOption<std::string>("minimizer");
	}
	if ( tag->hasOption("init") ) {
		init_ = tag->getOption<bool>("init");
	}
	if ( tag->hasOption("exact") ) {
		exact_ = tag->getOption<bool>("exact");
	}
	if ( tag->hasOption("verbose") ) {
		verbose_ = tag->getOption<bool>("verbose");
	}

	if ( tag->hasOption("deriv_check") ) {
		deriv_check_ = tag->getOption<bool>("deriv_check");
	}

	//
	if ( tag->hasOption("opt_to_fsc") ) {
		opt_to_fsc_ = tag->getOption<bool>("opt_to_fsc");
	}
	if ( tag->hasOption("res_low") ) {
		res_low_ = tag->getOption<core::Real>("res_low");
	}
	if ( tag->hasOption("res_high") ) {
		res_high_ = tag->getOption<core::Real>("res_high");
	}
	if ( tag->hasOption("nresbins") ) {
		nresbins_ = tag->getOption<core::Size>("nresbins");
	}

	if ( tag->hasOption("weights_to_scan") ) {
		std::string weight_string = tag->getOption<std::string>("weights_to_scan");
		utility::vector1< std::string > tokens ( utility::string_split( weight_string, ',' ) );
		weights_to_scan_.clear();
		for (uint j = 1; j <= tokens.size(); ++j)
			weights_to_scan_.push_back( atof(tokens[j].c_str()) );
	}
}

protocols::moves::MoverOP
BfactorFittingMoverCreator::create_mover() const {
	return new BfactorFittingMover;
}

std::string
BfactorFittingMoverCreator::keyname() const
{
	return BfactorFittingMoverCreator::mover_name();
}

std::string
BfactorFittingMoverCreator::mover_name()
{
	return "BfactorFitting";
}

}
}
