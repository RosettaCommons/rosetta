// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.electron_density.BfactorFittingMover" );

namespace protocols {
namespace electron_density {

using core::scoring::electron_density::poseCoords;
using core::scoring::electron_density::poseCoord;

/// helper function
void
symmetrizeBfactors( core::pose::Pose & pose ) {

	if ( !core::pose::symmetry::is_symmetric(pose) ) {
		return;
	}

	core::conformation::symmetry::SymmetryInfoOP symm_info;
	auto & SymmConf (
		dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
	symm_info = SymmConf.Symmetry_Info();

	for ( core::Size resid=1; resid<=pose.size(); ++resid ) {
		if ( symm_info && !symm_info->bb_is_independent( resid ) ) continue;
		core::conformation::Residue const & rsd( pose.residue(resid) );
		if ( rsd.aa() == core::chemical::aa_vrt ) continue;

		for ( auto pos=symm_info->bb_clones( resid ).begin(),
				epos=symm_info->bb_clones( resid ).end(); pos != epos; ++pos ) {
			core::Size natoms = rsd.natoms();
			for ( core::Size atmid = 1; atmid <= natoms; ++atmid ) {
				core::Real B = pose.pdb_info()->temperature( resid, atmid );
				pose.pdb_info()->temperature( *pos, atmid , B );
			}
		}
	}
}


BfactorMultifunc::BfactorMultifunc(
	core::pose::Pose & pose_in,
	core::Real wt_adp,
	core::Real wt_dens,
	core::Real rmax,
	core::Real radius_exp,
	core::Real scorescale,
	bool exact, bool verbose, bool deriv_check ) :
	pose_(pose_in),
	scorescale_( scorescale ),
	wt_adp_(wt_adp),
	wt_dens_(wt_dens),
	rmax_(rmax),
	radius_exp_(radius_exp),
	exact_( exact ),
	verbose_(verbose),
	deriv_check_( deriv_check )
{
	B_EPS=0.001;

	// map AtomIDs <-> indices
	core::pose::initialize_atomid_map( atom_indices_, pose_in );
	core::conformation::symmetry::SymmetryInfoOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose_in) ) {
		auto & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose_in.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}

	for ( Size resid=1; resid<=pose_in.size(); ++resid ) {
		if ( symm_info && !symm_info->bb_is_independent( resid ) ) continue;
		core::conformation::Residue const & rsd( pose_in.residue(resid) );
		if ( rsd.aa() == core::chemical::aa_vrt ) continue;

		core::Size natoms = rsd.nheavyatoms();
		for ( core::Size atmid = 1; atmid <= natoms; ++atmid ) {
			if ( rsd.atom_type( atmid ).is_virtual() ) continue;
			core::id::AtomID a_i( atmid, resid );
			moving_atoms_.push_back(a_i);
			atom_indices_[ a_i ] = moving_atoms_.size();
		}
	}
}


///  pose with B factors --> multivec
void
BfactorMultifunc::poseBfacts2multivec( core::pose::Pose & pose, core::optimization::Multivec &y ) const {
	runtime_assert( pose.pdb_info() != nullptr );  // must exist

	core::Size natoms = moving_atoms_.size();
	y.resize(natoms);

	for ( core::Size i = 1; i <= natoms; ++i ) {
		core::id::AtomID a_i = moving_atoms_[i];
		core::Real B = pose.pdb_info()->temperature( a_i.rsd(), a_i.atomno() );
		y[i] = log( B );
	}
}


///   multivec ---> pose B factors
void
BfactorMultifunc::multivec2poseBfacts( core::optimization::Multivec const &y, core::pose::Pose & pose ) const {
	runtime_assert( pose.pdb_info() != nullptr );  // must exist

	core::Size natoms = moving_atoms_.size();
	runtime_assert( y.size() == natoms );

	for ( core::Size i = 1; i <= natoms; ++i ) {
		core::id::AtomID a_i = moving_atoms_[i];
		// fpd! upper limit must match what is in xray_scattering.hh
		//  THIS SHOULD BE SELECTABLE
		pose.pdb_info()->temperature( a_i.rsd(), a_i.atomno(), std::min(std::max( exp( y[i] ) ,0.0),600.0) );
	}
}


///  B factor optimization multifunc
core::Real
BfactorMultifunc::operator ()( core::optimization::Multivec const & vars ) const {
	// [[1]] density score
	core::pose::Pose pose_copy = pose_;  // inefficient?
	multivec2poseBfacts( vars, pose_copy );

	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose_copy) ) {
		auto const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose_copy.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}

	core::Real dens_score = 0;
	if ( !exact_ ) {
		for ( Size i = 1; i <= pose_copy.size(); ++i ) {
			if ( symm_info && !symm_info->bb_is_independent( i ) ) continue;
			core::conformation::Residue const & rsd ( pose_copy.residue(i) );
			if ( rsd.aa() == core::chemical::aa_vrt ) continue;
			dens_score -= core::scoring::electron_density::getDensityMap().matchResFast( i, rsd, pose_copy, nullptr );
		}
	} else {
		core::scoring::electron_density::poseCoords litePose;
		for ( core::Size i = 1; i <= pose_copy.size(); ++i ) {
			if ( symm_info && !symm_info->bb_is_independent( i ) ) continue;
			core::conformation::Residue const & rsd_i ( pose_copy.residue(i) );
			if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;

			core::Size natoms = rsd_i.nheavyatoms();
			for ( core::Size j = 1; j <= natoms; ++j ) {
				if ( rsd_i.atom_type(j).is_virtual() ) continue;
				core::chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );

				core::scoring::electron_density::poseCoord coord_j;
				coord_j.x_ = rsd_i.xyz( j );
				coord_j.B_ = pose_copy.pdb_info()->temperature( i, j );
				coord_j.elt_ = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
				litePose.push_back( coord_j );
			}
		}

		core::scoring::electron_density::getDensityMap().calcRhoC( litePose, 0.0, rhoC_, rhoMask_, 600.0 );
		dens_score -= (moving_atoms_.size()/10.0) * core::scoring::electron_density::getDensityMap().getRSCC(rhoC_, rhoMask_);
	}

	// [[2]]  constraint score
	core::Real cst_score=0;
	core::Size nedge=0;
	if ( wt_adp_ != 0 ) {
		core::scoring::EnergyGraph const & energy_graph( pose_copy.energies().energy_graph() );

		for ( Size i = 1; i <= pose_copy.size(); ++i ) {
			if ( symm_info && !symm_info->bb_is_independent( i ) ) continue;
			core::conformation::Residue const & rsd1 ( pose_copy.residue(i) );
			if ( rsd1.aa() == core::chemical::aa_vrt ) continue;

			core::Size natoms = rsd1.nheavyatoms();
			for ( core::Size k = 1; k <= natoms; ++k ) {
				core::Real B_k = pose_copy.pdb_info()->temperature( i, k );
				if ( atom_indices_[core::id::AtomID(k,i)] == 0 ) continue;

				for ( core::Size l = k + 1; l <= natoms; ++l ) {
					core::Real B_l = pose_copy.pdb_info()->temperature( i, l );
					if ( atom_indices_[core::id::AtomID(l,i)] == 0 ) continue;

					core::Real dist_kl = (rsd1.atom( k ).xyz() - rsd1.atom( l ).xyz()).length();
					if ( radius_exp_ != 1.0 ) dist_kl = std::pow( dist_kl, radius_exp_ );

					if ( dist_kl < rmax_ ) {
						cst_score += (B_k-B_l)*(B_k-B_l) / ( dist_kl*(B_k+B_l+B_EPS) );
						nedge++;
					}
				}
			}

			for ( utility::graph::Graph::EdgeListConstIter
					iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
					irue = energy_graph.get_node(i)->const_upper_edge_list_end();
					iru != irue; ++iru ) {
				auto const * edge( static_cast< core::scoring::EnergyEdge const *> (*iru) );
				core::Size const j( edge->get_other_ind(i) );
				core::Size jasu = j;

				if ( symm_info && !symm_info->bb_is_independent( j ) ) {
					jasu = symm_info->bb_follows( j );
				}

				core::conformation::Residue const & rsd2 ( pose_copy.residue(j) );
				if ( rsd2.aa() == core::chemical::aa_vrt ) continue;
				core::Size natoms2 = rsd2.nheavyatoms();

				// inter-atom
				for ( core::Size k = 1; k <= natoms; ++k ) {
					core::Real B_k = pose_copy.pdb_info()->temperature( i, k );
					if ( atom_indices_[core::id::AtomID(k,i)] == 0 ) continue;
					for ( core::Size l = 1; l <= natoms2; ++l ) {
						core::Real B_l = pose_copy.pdb_info()->temperature( jasu, l );
						if ( atom_indices_[core::id::AtomID(l,jasu)] == 0 ) continue;

						core::Real dist_kl = (rsd1.atom( k ).xyz() - rsd2.atom( l ).xyz()).length();
						if ( radius_exp_ != 1.0 ) dist_kl = std::pow( dist_kl, radius_exp_ );

						if ( dist_kl < rmax_ ) {
							cst_score += (B_k-B_l)*(B_k-B_l) / ( dist_kl*(B_k+B_l+B_EPS) );
							nedge++;
						}
					}
				}
			}
		}
	}

	if ( verbose_ ) {
		//core::optimization::Multivec varsCopy = vars;
		//std::cerr << "[score] evaluated " << natom << " atoms and " << nedge << " edges" << std::endl;
		std::cerr << "[score] score = " << scorescale_* (wt_dens_*dens_score + wt_adp_*cst_score)
			<< " [ " << dens_score*10.0/moving_atoms_.size() << " , " << cst_score << " ] " << std::endl;
	}

	return ( scorescale_* (wt_dens_*dens_score + wt_adp_*cst_score) );
}

void
BfactorMultifunc::dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const {
	core::pose::Pose pose_copy = pose_;

	multivec2poseBfacts( vars, pose_copy );
	dE_dvars.clear();
	dE_dvars.resize( vars.size(), 0.0 );

	if ( !exact_ ) {
		// [[1]] density deriv
		for ( core::Size i = 1; i <= vars.size(); ++i ) {
			core::id::AtomID a_i = moving_atoms_[i];
			core::Real dCCdb = core::scoring::electron_density::getDensityMap().dCCdB_fastRes(
				a_i.atomno(), a_i.rsd(), pose_copy.residue( a_i.rsd() ), pose_copy );
			dE_dvars[ i ] = -scorescale_*wt_dens_*dCCdb * exp(vars[i]);
		}
	} else {
		core::scoring::electron_density::getDensityMap().dCCdBs( pose_copy, dE_dvars, rhoMask_ );
		for ( core::Size i = 1; i <= vars.size(); ++i ) {
			dE_dvars[i] *= -(moving_atoms_.size()/10.0)*scorescale_*wt_dens_ * exp(vars[i]);
		}
	}

	// [[2]]  constraint deriv
	core::Size nedge=0;
	if ( wt_adp_ != 0 ) {
		core::conformation::symmetry::SymmetryInfoOP symm_info;
		if ( core::pose::symmetry::is_symmetric(pose_copy) ) {
			auto & SymmConf (
				dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose_copy.conformation()) );
			symm_info = SymmConf.Symmetry_Info();
		}

		core::scoring::EnergyGraph const & energy_graph( pose_copy.energies().energy_graph() );

		for ( core::Size i = 1; i <= pose_copy.size(); ++i ) {
			if ( symm_info && !symm_info->bb_is_independent( i ) ) continue;
			core::conformation::Residue const & rsd1 ( pose_copy.residue(i) );
			if ( rsd1.aa() == core::chemical::aa_vrt ) continue;

			// intra-res
			core::Size natoms = rsd1.nheavyatoms();
			for ( core::Size k = 1; k <= natoms; ++k ) {
				core::Real B_k = pose_copy.pdb_info()->temperature( i, k );
				if ( atom_indices_[core::id::AtomID(k,i)] == 0 ) continue;
				for ( core::Size l = k + 1; l <= natoms; ++l ) {
					core::Real B_l = pose_copy.pdb_info()->temperature( i, l );
					if ( atom_indices_[core::id::AtomID(l,i)] == 0 ) continue;

					core::Real dist_kl = (rsd1.atom( k ).xyz() - rsd1.atom( l ).xyz()).length();
					if ( radius_exp_ != 1.0 ) dist_kl = std::pow( dist_kl, radius_exp_ );

					if ( dist_kl < rmax_ ) {
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

			for ( utility::graph::Graph::EdgeListConstIter
					iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
					irue = energy_graph.get_node(i)->const_upper_edge_list_end();
					iru != irue; ++iru ) {
				auto const * edge( static_cast< core::scoring::EnergyEdge const *> (*iru) );
				core::Size const j( edge->get_other_ind(i) );
				core::Size jasu = j;

				if ( symm_info && !symm_info->bb_is_independent( j ) ) {
					jasu = symm_info->bb_follows( j );
				}

				core::conformation::Residue const & rsd2 ( pose_copy.residue(j) );
				if ( rsd2.aa() == core::chemical::aa_vrt ) continue;
				core::Size natoms2 = rsd2.nheavyatoms();

				// inter-res
				for ( core::Size k = 1; k <= natoms; ++k ) {
					core::Real B_k = pose_copy.pdb_info()->temperature( i, k );
					if ( atom_indices_[core::id::AtomID(k,i)] == 0 ) continue;

					for ( core::Size l = 1; l <= natoms2; ++l ) {
						core::Real B_l = pose_copy.pdb_info()->temperature( jasu, l );
						if ( atom_indices_[core::id::AtomID(l,jasu)] == 0 ) continue;

						core::Real dist_kl = (rsd1.atom( k ).xyz() - rsd2.atom( l ).xyz()).length();
						if ( radius_exp_ != 1.0 ) dist_kl = std::pow( dist_kl, radius_exp_ );

						if ( dist_kl < rmax_ ) {
							core::Size atomK = atom_indices_[ core::id::AtomID( k,i ) ];
							core::Size atomL = atom_indices_[ core::id::AtomID( l,jasu ) ];

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

	if ( verbose_ ) {
		core::Real score = (*this)(vars);
		std::cerr << "[deriv] score = " << score << std::endl;
	}
}

void
BfactorMultifunc::dump( core::optimization::Multivec const & x1, core::optimization::Multivec const & ) const {
	// debug
	if ( deriv_check_ ) {
		core::optimization::Multivec varsCopy = x1;

		core::optimization::Multivec dE_dvars;
		this->dfunc(x1,dE_dvars);

		core::Real score = (*this)(x1);
		std::cerr << "[deriv] score = " << score << std::endl;
		for ( core::Size i = 1; i <= x1.size(); ++i ) {
			varsCopy[i]+=B_EPS;
			core::Real scorep = (*this)(varsCopy);
			varsCopy[i]-=2*B_EPS;
			core::Real scoren = (*this)(varsCopy);
			varsCopy[i]+=B_EPS;
			std::cerr << "[deriv] x:" << i << "  B: " << x1[i] << "  A: " << dE_dvars[i] << "  N: "
				<< (scorep-scoren)/2*B_EPS << "    r: " << 2*B_EPS*dE_dvars[i]/(scorep-scoren) << std::endl;
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
	//minimizer_ = "lbfgs_armijo";
	minimizer_ = "lbfgs_armijo_nonmonotone_atol";
	init_ = true;
	exact_ = false;
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
	if ( !basic::options::option[ basic::options::OptionKeys::cryst::crystal_refine ]() ) {
		utility_exit_with_message( "Bfactor optimization _must_ be run with the -cryst::crystal_refine option" );
	}

	// make sure pose has pdbinfo
	if ( ! pose.pdb_info() ) {
		core::pose::PDBInfoOP newinfo( new core::pose::PDBInfo(pose) );
		pose.pdb_info( newinfo );
	}

	core::optimization::MinimizerOptions options( minimizer_, 1e-4, true, deriv_check_, deriv_check_ );
	options.max_iter(max_iter_);

	// set up optimizer
	BfactorMultifunc f_b( pose, wt_adp_, wt_dens_, rmax_, radius_exp_, scorescale_, exact_, verbose_, deriv_check_);
	core::optimization::Minimizer minimizer( f_b, options );

	// make sure interaction graphs are set up
	core::scoring::ScoreFunctionOP sf;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		sf = core::scoring::ScoreFunctionOP( new core::scoring::symmetry::SymmetricScoreFunction() );
	} else {
		sf = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() );
	}
	sf->set_weight( core::scoring::fa_rep, 1.0 );
	(*sf)(pose);

	// ... and run
	core::optimization::Multivec y;
	f_b.poseBfacts2multivec( pose, y );

	// if b=0 set b to some small value
	for ( core::Size i = 1; i <= y.size(); ++i ) if ( y[i] < 1e-6 ) y[i] = 5;

	core::scoring::electron_density::getDensityMap().rescale_fastscoring_temp_bins( pose, init_ );

	// optionally grid search for best starting value
	core::Size nvalues = weights_to_scan_.size();
	core::Real bestVal = 0, bestScore = 0;

	if ( init_ ) {
		for ( core::Size j=1; j<=nvalues; ++j ) {
			for ( core::Size i = 1; i <= y.size(); ++i ) {
				y[i] = log( weights_to_scan_[j] );
			}
			core::Real thisScore = f_b(y);
			TR << "B=" << weights_to_scan_[j] << ": score=" << thisScore << std::endl;
			if ( thisScore<bestScore || j==1 ) {
				bestScore = thisScore;
				bestVal = weights_to_scan_[j];
			}
		}
		for ( core::Size i = 1; i <= y.size(); ++i ) y[i] = log(bestVal);
	}

	// MAIN LOOP
	minimizer.run( y );

	f_b.multivec2poseBfacts( y, pose );
	core::scoring::cryst::fix_bfactorsH( pose ); // fix hydrogen B's to riding atoms
	symmetrizeBfactors( pose );
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
BfactorFittingMover::parse_my_tag(TagCOP const tag, basic::datacache::DataMap &, Filters_map const &,
	moves::Movers_map const &, Pose const &)
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
		for ( core::Size j = 1; j <= tokens.size(); ++j ) {
			weights_to_scan_.push_back( atof(tokens[j].c_str()) );
		}
	}
}

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP BfactorFittingMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new BfactorFittingMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP BfactorFittingMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return BfactorFittingMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP BfactorFittingMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "BfactorFitting";
// XRW TEMP }

std::string BfactorFittingMover::get_name() const {
	return mover_name();
}

std::string BfactorFittingMover::mover_name() {
	return "BfactorFitting";
}

void BfactorFittingMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute("wt_adp", xsct_real, "Weight on enforcing nearby atoms to take the same B factors. "
		"0.0005 is generally well-behaved, default=0.001")
		+ XMLSchemaAttribute("wt_dens", xsct_real, "Weight on density contraints, default=1.0.")
		+ XMLSchemaAttribute("rmax", xsct_real, "XRW TO DO, default=5.0")
		+ XMLSchemaAttribute("scale", xsct_real, "Scaling factor. Is multiplied with the weighted scores "
		"after wt_adp and wt_dens are applied, default=5000")
		+ XMLSchemaAttribute("radius_exp", xsct_non_negative_integer, "XRW TO DO, default=1")
		+ XMLSchemaAttribute("max_iter", xsct_non_negative_integer, "Maximum nr of round before exiting, default=200.")
		+ XMLSchemaAttribute("minimizer", xs_string, "Minimizer used for energy minimization. default=lbfgs_armijo_nonmonotone_atol")
		+ XMLSchemaAttribute("init", xsct_rosetta_bool,
		"init=1 means to do a quick scan of overall B factors before beginning refinement; "
		"if there is more than one call to this mover in a single trajectory, "
		"then only the first needs to have init=1")
		+ XMLSchemaAttribute("exact", xsct_rosetta_bool,
		"Exact B-factor refinement - recommended. default=false "
		"Approximate version can behave poorly and is memory intensive.")
		+ XMLSchemaAttribute("verbose", xsct_rosetta_bool, "Increase level of output information. default=false")
		+ XMLSchemaAttribute("deriv_check", xsct_rosetta_bool, "XRW TO DO, default=false")
		+ XMLSchemaAttribute("res_low", xsct_real, "Low resultion cut off (Angstroms)")
		+ XMLSchemaAttribute("res_high", xsct_real, "High resolution cut off (Angstrom)")
		+ XMLSchemaAttribute("nresbins", xsct_non_negative_integer, "This option is not used.")
		+ XMLSchemaAttribute("weights_to_scan", xsct_real_cslist, "List of ??? weights to scan. Finds the optimum one to be used. Default=\"10,20,...,300\"");
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Fit atomic B factors to maximize model-map correlation.",
		attlist );
}

std::string BfactorFittingMoverCreator::keyname() const {
	return BfactorFittingMover::mover_name();
}

protocols::moves::MoverOP
BfactorFittingMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new BfactorFittingMover );
}

void BfactorFittingMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BfactorFittingMover::provide_xml_schema( xsd );
}


}
}
