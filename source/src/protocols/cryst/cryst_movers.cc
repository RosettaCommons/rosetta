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

// keep first
#include <core/scoring/cryst/PhenixInterface.hh>

#include <protocols/cryst/cryst_movers.hh>
#include <protocols/cryst/cryst_movers_creator.hh>

//#include <protocols/comparative_modeling/coord_util.hh>
#include <protocols/jd2/util.hh>

#include <core/util/cryst_util.hh>
#include <core/scoring/cryst/util.hh>
#include <core/scoring/cryst/XtalMLEnergy.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/coulomb/Coulomb.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/elec/FA_ElecEnergy.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/lkball/LK_BallEnergy.hh>

#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/Remarks.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/select/jump_selector/JumpIndexSelector.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>

/////////////
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/types.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/AtomTreeMultifunc.hh>

#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/symmetry/SymMinimizerMap.hh>
#include <core/optimization/symmetry/SymAtomTreeMultifunc.hh>

#include <core/optimization/CartesianMinimizerMap.hh>
#include <core/optimization/CartesianMultifunc.hh>
#include <core/optimization/CartesianMinimizer.hh>

///////////

#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>
#include <numeric/model_quality/rms.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <fstream>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using basic::Error;
using basic::Warning;

namespace protocols {
namespace cryst {

static basic::Tracer TR( "protocols.cryst.cryst_movers" );

using namespace protocols;
using namespace core;
using namespace kinematics;
using namespace scoring;


//////////////////
//////////////////
/// creators


// XRW TEMP std::string
// XRW TEMP ReportGradientsMoverCreator::keyname() const {
// XRW TEMP  return ReportGradientsMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP ReportGradientsMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new ReportGradientsMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ReportGradientsMover::mover_name() {
// XRW TEMP  return "ReportGradients";
// XRW TEMP }

///

// XRW TEMP std::string
// XRW TEMP SetCrystWeightMoverCreator::keyname() const {
// XRW TEMP  return SetCrystWeightMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SetCrystWeightMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SetCrystWeightMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SetCrystWeightMover::mover_name() {
// XRW TEMP  return "SetCrystWeight";
// XRW TEMP }


// XRW TEMP std::string
// XRW TEMP RecomputeDensityMapMoverCreator::keyname() const {
// XRW TEMP  return RecomputeDensityMapMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP RecomputeDensityMapMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new RecomputeDensityMapMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP RecomputeDensityMapMover::mover_name() {
// XRW TEMP  return "RecomputeDensityMap";
// XRW TEMP }


// XRW TEMP std::string
// XRW TEMP LoadDensityMapMoverCreator::keyname() const {
// XRW TEMP  return LoadDensityMapMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP LoadDensityMapMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new LoadDensityMapMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP LoadDensityMapMover::mover_name() {
// XRW TEMP  return "LoadDensityMap";
// XRW TEMP }


// XRW TEMP std::string
// XRW TEMP FitBfactorsMoverCreator::keyname() const {
// XRW TEMP  return FitBfactorsMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP FitBfactorsMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new FitBfactorsMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP FitBfactorsMover::mover_name() {
// XRW TEMP  return "FitBfactors";
// XRW TEMP }


// XRW TEMP std::string
// XRW TEMP UpdateSolventMoverCreator::keyname() const {
// XRW TEMP  return UpdateSolventMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP UpdateSolventMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new UpdateSolventMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP UpdateSolventMover::mover_name() {
// XRW TEMP  return "UpdateSolvent";
// XRW TEMP }


// XRW TEMP std::string
// XRW TEMP TagPoseWithRefinementStatsMoverCreator::keyname() const {
// XRW TEMP  return TagPoseWithRefinementStatsMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP TagPoseWithRefinementStatsMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new TagPoseWithRefinementStatsMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP TagPoseWithRefinementStatsMover::mover_name() {
// XRW TEMP  return "TagPoseWithRefinementStats";
// XRW TEMP }


// XRW TEMP std::string
// XRW TEMP SetRefinementOptionsMoverCreator::keyname() const {
// XRW TEMP  return SetRefinementOptionsMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SetRefinementOptionsMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SetRefinementOptionsMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SetRefinementOptionsMover::mover_name() {
// XRW TEMP  return "SetRefinementOptions";
// XRW TEMP }

//////////////////
//////////////////
/// movers

void ReportGradientsMover::apply( core::pose::Pose & pose ) {
	compute (pose);
}

core::Real ReportGradientsMover::compute(core::pose::Pose & pose ) {
	using namespace core::scoring;
	using namespace core::optimization;
	using namespace core::optimization::symmetry;

	core::scoring::ScoreFunctionOP reference_scorefxn( score_function_->clone() );
	core::scoring::ScoreFunctionOP working_scorefxn( new core::scoring::ScoreFunction() );

	working_scorefxn->set_weight( core::scoring::xtal_ml, 0.0 );
	reference_scorefxn->set_weight( core::scoring::xtal_ml, 0.0 );

	core::kinematics::MoveMap move_map;
	move_map.set_bb  ( true );
	move_map.set_chi ( true );
	move_map.set_jump( true );

	// additional setup for symmetric poses
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		auto const & symm_conf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation() ) );
		core::conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

		// symmetrize scorefunct & movemap
		working_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *working_scorefxn );
		core::pose::symmetry::make_symmetric_movemap( pose, move_map );
	}

	// compute gradients using both scorefunctions
	CartesianMinimizerMap min_map;
	min_map.setup( pose, move_map );

	Multivec vars( min_map.ndofs() );
	min_map.copy_dofs_from_pose( pose, vars );

	Multivec dEros_dvars;
	(*reference_scorefxn)(pose);  // score pose first
	reference_scorefxn->setup_for_minimizing( pose, min_map );
	CartesianMultifunc f_ros( pose, min_map, *reference_scorefxn, false, false );
	f_ros.dfunc( vars, dEros_dvars );

	// lkball setup
	core::scoring::lkball::LK_BallEnergy lkb( reference_scorefxn->energy_method_options() );
	lkb.setup_for_scoring(pose, *reference_scorefxn);

	utility::vector1< Multivec > dEros_i_dvars;
	utility::vector1< core::scoring::ScoreType > term_i;

	// per score-term gradients
	// used for:
	//  a) verbose output
	//  b) normalization for non-etable derivs

	for ( int ii=1; ii<=(int)core::scoring::n_score_types; ++ii ) {
		auto st_ii = (core::scoring::ScoreType)ii;

		if ( reference_scorefxn->get_weight( (core::scoring::ScoreType)ii ) == 0.0 ) continue;

		term_i.push_back( (core::scoring::ScoreType)ii );
		dEros_i_dvars.push_back( Multivec() );

		working_scorefxn->set_weight( (core::scoring::ScoreType)ii , reference_scorefxn->get_weight( (core::scoring::ScoreType)ii )  );

		// lk ball hack
		if ( st_ii == lk_ball_wtd || st_ii == lk_ball || st_ii == lk_ball_iso || st_ii == lk_ball_bridge || st_ii == lk_ball_bridge_uncpl ) {
			working_scorefxn->set_weight( fa_atr, 1e-9 );
		}


		(*working_scorefxn)(pose);
		working_scorefxn->setup_for_minimizing( pose, min_map );
		CartesianMultifunc f_ros_i( pose, min_map, *working_scorefxn, false, false );
		f_ros_i.dfunc( vars, dEros_i_dvars[ dEros_i_dvars.size() ] );

		working_scorefxn->set_weight( (core::scoring::ScoreType)ii , 0.0  );
		if ( st_ii == lk_ball_wtd || st_ii == lk_ball || st_ii == lk_ball_iso || st_ii == lk_ball_bridge || st_ii == lk_ball_bridge_uncpl ) {
			working_scorefxn->set_weight( fa_atr, 0 );
		}
	}

	// compute the total gradient sum
	core::Size natoms = dEros_dvars.size()/3;  // integer division, numerator is always divisible by 3
	core::Real gradsum( 0.0 );

	for ( core::Size ii_atm=0; ii_atm<natoms; ++ii_atm ) {
		id::AtomID id = min_map.get_atom( ii_atm+1 );
		numeric::xyzVector< core::Real > grad_i(dEros_dvars[3*ii_atm+1], dEros_dvars[3*ii_atm+2], dEros_dvars[3*ii_atm+3]);

		if ( grad_i.length()<1e-6 ) continue;

		// etable terms
		core::Real norm = normalization(pose,  min_map.get_atom( ii_atm+1 ), reference_scorefxn );

		// other terms
		for ( int ii=1; ii<=(int)term_i.size(); ++ii ) {
			if ( term_i[ii] == fa_rep || term_i[ii] == fa_atr || term_i[ii] == fa_sol ) continue;
			if ( term_i[ii] == fa_elec ) continue;
			if ( term_i[ii] == fa_intra_rep || term_i[ii] == fa_intra_atr || term_i[ii] == fa_intra_sol ) continue;
			if ( term_i[ii] == fa_intra_rep_xover4 || term_i[ii] == fa_intra_atr_xover4 || term_i[ii] == fa_intra_sol_xover4 ) continue;

			numeric::xyzVector< core::Real > grad_ij (dEros_i_dvars[ii][3*ii_atm+1], dEros_i_dvars[ii][3*ii_atm+2], dEros_i_dvars[ii][3*ii_atm+3]);
			norm += grad_ij.length();
		}

		if ( norm<1e-6 ) {
			TR << "ERROR! at atom " << id << " got grad = " << grad_i.length() << " with norm " << norm << std::endl;
		}

		gradsum += grad_i.length() / norm;

	}
	gradsum /= (core::Real) natoms;
	TR << "[grad] average gradient " << gradsum << std::endl;

	// (optionally) report per-atom/per-scoreterm gradients
	if ( verbose_ ) {
		// header
		TR << "[grad] restype resid chnid atomname grad_unnorm norm grad_norm ";
		for ( int ii=1; ii<=(int)term_i.size(); ++ii ) {
			TR << term_i[ii] << " ";
		}
		TR << std::endl;

		for ( core::Size ii_atm=0; ii_atm<natoms; ++ii_atm ) {
			id::AtomID id = min_map.get_atom( ii_atm+1 );
			core::conformation::Residue const & rsd_i = pose.residue( id.rsd() );

			// ASSUMES PDBINFO IS PRESENT!!!!
			core::Size pdb_res = pose.pdb_info()->number( id.rsd() );
			char pdb_chain = pose.pdb_info()->chain( id.rsd() );

			if ( id.atomno() <= rsd_i.natoms() ) {
				numeric::xyzVector< core::Real > grad_i(dEros_dvars[3*ii_atm+1], dEros_dvars[3*ii_atm+2], dEros_dvars[3*ii_atm+3]);
				TR << "[grad] " << rsd_i.name3() << " " << pdb_res << " " << pdb_chain<< " " << rsd_i.atom_name( id.atomno() );

				// etable terms
				core::Real norm = normalization(pose,  min_map.get_atom( ii_atm+1 ), reference_scorefxn );

				// other terms
				for ( int jj=1; jj<=(int)term_i.size(); ++jj ) {
					if ( term_i[jj] == fa_rep || term_i[jj] == fa_atr || term_i[jj] == fa_sol ) continue;
					if ( term_i[jj] == fa_elec ) continue;
					if ( term_i[jj] == fa_intra_rep || term_i[jj] == fa_intra_atr || term_i[jj] == fa_intra_sol ) continue;
					if ( term_i[jj] == fa_intra_rep_xover4 || term_i[jj] == fa_intra_atr_xover4 || term_i[jj] == fa_intra_sol_xover4 ) continue;

					numeric::xyzVector< core::Real > grad_ij(dEros_i_dvars[jj][3*ii_atm+1], dEros_i_dvars[jj][3*ii_atm+2], dEros_i_dvars[jj][3*ii_atm+3]);
					norm += grad_ij.length();
				}

				TR << " " << grad_i.length() << " " << norm << " " << grad_i.length()/norm;
				for ( int jj=1; jj<=(int)term_i.size(); ++jj ) {
					numeric::xyzVector< core::Real > grad_ij(dEros_i_dvars[jj][3*ii_atm+1], dEros_i_dvars[jj][3*ii_atm+2], dEros_i_dvars[jj][3*ii_atm+3]);
					TR << " " << grad_ij.length();
				}
				TR << std::endl;
			}
		}
	}

	if ( outfile_.length() > 0 ) {
		std::ofstream OUTF(outfile_.c_str());
		OUTF << gradsum << std::endl;
	}
	return gradsum;
}

core::Real ReportGradientsMover::normalization(core::pose::Pose & pose, core::id::AtomID atmid, core::scoring::ScoreFunctionOP sfxn ) {
	using namespace core;
	using namespace core::scoring;
	using namespace core::scoring::etable;
	using namespace core::scoring::etable::coulomb;
	using namespace core::scoring::etable::count_pair;
	using namespace core::scoring::methods;

	core::Real symmscale=1.0;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		auto const & symm_conf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation() ) );
		core::conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
		symmscale = symm_info->score_multiply_factor();
	}

	(*sfxn)(pose); // make sure energies obj initialized

	core::Real norm=0.0;

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	methods::EnergyMethodOptions e_opts = sfxn->energy_method_options();
	Etable etable = *( ScoringManager::get_instance()->etable( e_opts ).lock() ); // copy

	core::Real fa_atr_wt = sfxn->get_weight( core::scoring::fa_atr );
	core::Real fa_rep_wt = sfxn->get_weight( core::scoring::fa_rep );
	core::Real fa_sol_wt = sfxn->get_weight( core::scoring::fa_sol );
	core::Real fa_elec_wt = sfxn->get_weight( core::scoring::fa_elec );

	core::Real fa_intra_atr_wt = sfxn->get_weight( core::scoring::fa_intra_atr );
	core::Real fa_intra_rep_wt = sfxn->get_weight( core::scoring::fa_intra_rep );
	core::Real fa_intra_sol_wt = sfxn->get_weight( core::scoring::fa_intra_sol );

	core::Real fa_intra_atr_x4_wt = sfxn->get_weight( core::scoring::fa_intra_atr_xover4 );
	core::Real fa_intra_rep_x4_wt = sfxn->get_weight( core::scoring::fa_intra_rep_xover4 );
	core::Real fa_intra_sol_x4_wt = sfxn->get_weight( core::scoring::fa_intra_sol_xover4 );

	Size ires = atmid.rsd();
	Size iatm = atmid.atomno();
	core::conformation::Residue const & rsd1( pose.residue(ires) );

	// needed for electrostatic calcs
	core::scoring::elec::FA_ElecEnergy faelec( e_opts );
	Coulomb coulomb( sfxn->energy_method_options() );
	coulomb.initialize();

	// get inter-res etable energies
	if ( fa_atr_wt>0 || fa_rep_wt>0 || fa_sol_wt>0 || fa_elec_wt>0 ) {
		for ( utility::graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(ires)->const_edge_list_begin(),
				irue = energy_graph.get_node(ires)->const_edge_list_end();
				iru != irue; ++iru ) {
			auto const * edge( static_cast< EnergyEdge const *> (*iru) );
			Size jres=edge->get_other_ind(ires);
			core::conformation::Residue const &rsd2( pose.residue(jres) );

			if ( ires == jres ) continue;

			// count pair - assume 4 for now (actually 3 for ligands/RNA or if rama is off)
			CountPairFunctionOP cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

			for ( Size jatm=1; jatm<= rsd2.natoms(); ++jatm ) {
				core::Real weight=1.0;
				core::Size path_dist;
				core::Real d_faatr=0, d_farep=0, d_fasol=0, d_faelec=0, invD; //, dummy;

				etable.analytic_etable_derivatives(rsd1.atom(iatm), rsd2.atom(jatm),  d_faatr, d_farep, d_fasol, invD );

				Real const dis = (rsd1.xyz(iatm)-rsd2.xyz(jatm)).length();
				Real const dis2 = dis*dis;

				d_faelec = coulomb.eval_dfa_elecE_dr_over_r( dis2, rsd1.atomic_charge(iatm), rsd2.atomic_charge(jatm) );

				// fa_rep/atr/sol
				if ( cpfxn->count( iatm, jatm, weight, path_dist ) ) {
					norm += symmscale * weight *(
						fa_atr_wt*std::abs(d_faatr) +
						fa_rep_wt*std::abs(d_farep) +
						fa_sol_wt*std::abs(d_fasol));
				}

				// fa_elec
				core::Size iirep = faelec.get_countpair_representative_atom(rsd1.type(), iatm);
				core::Size jjrep = faelec.get_countpair_representative_atom(rsd2.type(), jatm);

				if ( cpfxn->count( iirep, jjrep, weight, path_dist ) ) {
					norm += symmscale * weight *(dis*fa_elec_wt*std::abs(d_faelec));
				}

			}
		}
	}

	// add fa_intra_* contribution
	if ( fa_intra_atr_wt>0 || fa_intra_rep_wt>0 || fa_intra_sol_wt>0 ) {
		CountPairFunctionOP cpfxn =
			CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_3 );
		for ( Size jatm=1; jatm<= rsd1.natoms(); ++jatm ) {
			core::Real weight=1.0;
			core::Size path_dist;
			core::Real d_faatr, d_farep, d_fasol, invD;
			if ( cpfxn->count( iatm, jatm, weight, path_dist ) ) {
				etable.analytic_etable_derivatives(rsd1.atom(iatm), rsd1.atom(jatm),  d_faatr, d_farep, d_fasol, invD );
				//Real const dis = (rsd1.xyz(iatm)-rsd1.xyz(jatm)).length();
				norm += symmscale * weight* (
					fa_intra_atr_wt*std::abs(d_faatr) +
					fa_intra_rep_wt*std::abs(d_farep) +
					fa_intra_sol_wt*std::abs(d_fasol) );
			}
		}
	}

	// add fa_intra_*_xover4 contribution
	if ( fa_intra_atr_x4_wt>0 || fa_intra_rep_x4_wt>0 || fa_intra_sol_x4_wt>0 ) {
		CountPairFunctionOP cpfxn =
			CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_4 );
		for ( Size jatm=1; jatm<= rsd1.natoms(); ++jatm ) {
			core::Real weight=1.0;
			core::Size path_dist;
			core::Real d_faatr, d_farep, d_fasol, invD;
			if ( cpfxn->count( iatm, jatm, weight, path_dist ) ) {
				etable.analytic_etable_derivatives(rsd1.atom(iatm), rsd1.atom(jatm),  d_faatr, d_farep, d_fasol, invD );
				//Real const dis = (rsd1.xyz(iatm)-rsd1.xyz(jatm)).length();
				norm += symmscale * weight*(
					fa_intra_atr_x4_wt*std::abs(d_faatr) +
					fa_intra_rep_x4_wt*std::abs(d_farep) +
					fa_intra_sol_x4_wt*std::abs(d_fasol) );
			}
		}
	}

	return norm;
}

void ReportGradientsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &data,
	filters::Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/ ) {
	score_function_ = protocols::rosetta_scripts::parse_score_function( tag, data );
	verbose_ = tag->getOption<bool>("verbose", false);
	outfile_ = tag->getOption<std::string>("outfile", "");
}

std::string ReportGradientsMover::get_name() const {
	return mover_name();
}

std::string ReportGradientsMover::mover_name() {
	return "ReportGradients";
}

void ReportGradientsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "verbose", xsct_rosetta_bool, "report per-atom/per-scoreterm gradients", "false" );
	attlist + XMLSchemaAttribute( "outfile", xs_string, "output file for gradients" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "PHENIX crystallographic refinement interface code.  Reports to file per score-term gradients, used for: a) verbose output; b) normalization for non-etable derivs", attlist );
}

std::string ReportGradientsMoverCreator::keyname() const {
	return ReportGradientsMover::mover_name();
}

protocols::moves::MoverOP
ReportGradientsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ReportGradientsMover );
}

void ReportGradientsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ReportGradientsMover::provide_xml_schema( xsd );
}


////

void SetCrystWeightMover::apply( core::pose::Pose & pose ) {
	if ( autoset_wt_ ) {
		// if autoset_wt_ guess at cryst wt
		//    use score_function_ref_ to do the scaling
		core::Real auto_weight = 0;
		if ( cartesian_ ) {
			if ( mmf_ ) {
				core::kinematics::MoveMapOP mm( mmf_->create_movemap_from_pose( pose ) );
				auto_weight = core::util::getMLweight_cart( *score_function_ref_, pose, *mm );
			} else {
				auto_weight = core::util::getMLweight_cart( *score_function_ref_, pose );
			}
		} else {
			if ( mmf_ ) {
				core::kinematics::MoveMapOP mm( mmf_->create_movemap_from_pose( pose ) );
				auto_weight = core::util::getMLweight( *score_function_ref_, pose, *mm );
			} else {
				auto_weight = core::util::getMLweight( *score_function_ref_, pose );
			}
		}

		// stay above min weight
		auto_weight = std::max( auto_weight, weight_min_ );

		score_function_->set_weight( core::scoring::xtal_ml, auto_weight * weight_scale_ );
		TR << "set xtal_weight to " << score_function_->get_weight( core::scoring::xtal_ml ) << std::endl;
	} else {
		TR << "override automatic weight with wt=" << weight_ * weight_scale_ << std::endl;
		score_function_->set_weight( core::scoring::xtal_ml, weight_ * weight_scale_ );
	}
}

void SetCrystWeightMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &data,
	filters::Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	core::pose::Pose const & ) {

	score_function_ = protocols::rosetta_scripts::parse_score_function( tag, data );

	if ( tag->hasOption("scorefxn_ref") ) {
		std::string const scorefxn_key( tag->getOption<std::string>("scorefxn_ref") );
		if ( ! data.has( "scorefxns", scorefxn_key ) ) {
			utility_exit_with_message("ScoreFunction " + scorefxn_key + " not found in basic::datacache::DataMap.");
		}
		score_function_ref_ = data.get_ptr< ScoreFunction >( "scorefxns", scorefxn_key );
	} else {
		score_function_ref_ = score_function_->clone();
	}

	mmf_ = protocols::rosetta_scripts::parse_movemap_factory_legacy( tag, data );

	// also can specify defaults
	if ( tag->hasOption("jump") ) {
		if ( ! mmf_ ) mmf_ = core::select::movemap::MoveMapFactoryOP( new core::select::movemap::MoveMapFactory );
		utility::vector1<std::string> jumps = utility::string_split( tag->getOption<std::string>( "jump" ), ',' );
		// string 'ALL' makes all jumps movable
		if ( jumps.size() == 1 && (jumps[1] == "ALL" || jumps[1] == "All" || jumps[1] == "all") ) {
			mmf_->all_jumps( true );
		} else {
			for ( std::vector<std::string>::const_iterator it = jumps.begin(); it != jumps.end(); ++it ) {
				Size const value = std::atoi( it->c_str() ); // convert to C string, then convert to integer, then set a Size (phew!)
				core::select::jump_selector::JumpIndexSelectorOP jumpselect( new core::select::jump_selector::JumpIndexSelector( value ) );
				mmf_->add_jump_action( core::select::movemap::mm_enable, jumpselect );
			}
		}
	}
	if ( tag->hasOption("chi") ) {
		bool const value( tag->getOption<bool>("chi") );
		if ( ! mmf_ ) mmf_ = core::select::movemap::MoveMapFactoryOP( new core::select::movemap::MoveMapFactory );
		mmf_->all_chi(value);
	}
	if ( tag->hasOption("bb") ) {
		bool const value( tag->getOption<bool>("bb") );
		if ( ! mmf_ ) mmf_ = core::select::movemap::MoveMapFactoryOP( new core::select::movemap::MoveMapFactory );
		mmf_->all_bb(value);
	}
	if ( tag->hasOption("bondangle") ) {
		bool const value( tag->getOption<bool>("bondangle") );
		if ( ! mmf_ ) mmf_ = core::select::movemap::MoveMapFactoryOP( new core::select::movemap::MoveMapFactory );
		mmf_->all_bondangles( value );
	}
	if ( tag->hasOption("bondlength") ) {
		bool const value( tag->getOption<bool>("bondlength") );
		if ( ! mmf_ ) mmf_ = core::select::movemap::MoveMapFactoryOP( new core::select::movemap::MoveMapFactory );
		mmf_->all_bondlengths( value );
	}

	if ( tag->hasOption("weight") ) {
		weight_ = tag->getOption<core::Real>("weight");
		autoset_wt_ = false;
	}

	if ( tag->hasOption("weight_min") ) {
		weight_min_ = tag->getOption<core::Real>("weight_min");
	}
	weight_scale_ = tag->getOption<core::Real>("weight_scale", 1.0);
	cartesian_ = tag->getOption<bool>("cartesian", false);
}

std::string SetCrystWeightMover::get_name() const {
	return mover_name();
}

std::string SetCrystWeightMover::mover_name() {
	return "SetCrystWeight";
}

std::string SetCrystWeightMover_subelement_ct_name( std::string const & name ) {
	return "SetCrystWeightMover_subelement_" + name + "Type";
}

void SetCrystWeightMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );

	attlist + XMLSchemaAttribute( "scorefxn_ref", xs_string, "String name of scorefunction defined in SCOREFUNCTION section; used to copy weights from; defaults to scorefunction defined by this element's scorefunction subelement)" );

	//empty SubelementList for MoveMap util function
	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( & SetCrystWeightMover_subelement_ct_name );
	rosetta_scripts::append_subelement_for_parse_movemap_factory_legacy(xsd, subelements);

	//argument legality checker for "jump":
	XMLSchemaRestriction jump_type;
	jump_type.name("jump_type");
	jump_type.base_type( xs_string );
	jump_type.add_restriction( xsr_pattern, "ALL|All|all|[0-9]+(,[0-9]+)*");
	xsd.add_top_level_element( jump_type );

	attlist + XMLSchemaAttribute( "jump", "jump_type", "if 'ALL', 'All', or 'all'; all jumps mobile in MoveMap.  else, interpreted as a comma-separated list of integers labeling jumps" );

	attlist + XMLSchemaAttribute( "chi", xsct_rosetta_bool, "chi mobile in MoveMap" );
	attlist + XMLSchemaAttribute( "bb", xsct_rosetta_bool, "backbone mobile in MoveMap" );
	attlist + XMLSchemaAttribute( "bondangle", xsct_rosetta_bool, "bond angles mobile in MoveMap" );
	attlist + XMLSchemaAttribute( "bondlength", xsct_rosetta_bool, "bond lengths mobile in MoveMap" );

	attlist + XMLSchemaAttribute( "weight", xsct_real, "Always set to this weight; overriding heuristics" );
	attlist + XMLSchemaAttribute( "weight_min", xsct_real, "Keep weight above this minimum" );

	attlist + XMLSchemaAttribute::attribute_w_default( "weight_scale", xsct_real, "scale weight by this factor (applies to automatic AND overriden weights)", "1.0");

	attlist + XMLSchemaAttribute::attribute_w_default( "cartesian", xsct_rosetta_bool, "Use cartesian ??? when automatically determining weights", "false" );

	moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "PHENIX crystallographic refinement interface code. Converts weights from one score function into weights for the xtal_ml term for some other scorefunction; MoveMap helps determine optimal weight", attlist, subelements);


	// XMLSchemaComplexTypeGenerator complex_type_generator;
	// complex_type_generator
	//  .element_name( mover_name() )
	//  .description( "PHENIX crystallographic refinement interface code. Converts weights from one score function into weights for the xtal_ml term for some other scorefunction; MoveMap helps determine optimal weight" )
	//  .complex_type_naming_func( & moves::complex_type_name_for_mover )
	//  .add_attributes( attlist )
	//  .set_subelements_repeatable( subelements ) //unclear that MoveMap is repeatable; not worth figuring out now
	//  .write_complex_type_to_schema( xsd );

	//protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string SetCrystWeightMoverCreator::keyname() const {
	return SetCrystWeightMover::mover_name();
}

protocols::moves::MoverOP
SetCrystWeightMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetCrystWeightMover );
}

void SetCrystWeightMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SetCrystWeightMover::provide_xml_schema( xsd );
}


////

void RecomputeDensityMapMover::apply( core::pose::Pose & pose )
{
	using namespace core::scoring::electron_density;
	std::string newMap = core::scoring::cryst::getPhenixInterface().calculateDensityMap( pose, !keep_sidechains_ );
	/*ElectronDensity& densityMap =*/ getDensityMap(newMap, true);
}

void RecomputeDensityMapMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	filters::Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/ )
{
	keep_sidechains_ = tag->getOption<bool>("sidechains", true);
}

std::string RecomputeDensityMapMover::get_name() const {
	return mover_name();
}

std::string RecomputeDensityMapMover::mover_name() {
	return "RecomputeDensityMap";
}

void RecomputeDensityMapMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default( "sidechains", xsct_rosetta_bool, "does nothing", "true" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "PHENIX crystallographic refinement interface code. Uses PHENIX to recompute the density map", attlist );
}

std::string RecomputeDensityMapMoverCreator::keyname() const {
	return RecomputeDensityMapMover::mover_name();
}

protocols::moves::MoverOP
RecomputeDensityMapMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RecomputeDensityMapMover );
}

void RecomputeDensityMapMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RecomputeDensityMapMover::provide_xml_schema( xsd );
}


////
void LoadDensityMapMover::apply( core::pose::Pose & /*pose*/ ) {
	using namespace core::scoring::electron_density;

	ElectronDensity& densityMap = getDensityMap(mapfile_, true);

	//fpd make these selectable?
	densityMap.setWindow(window_);
	densityMap.setSCscaling(sc_scale_);
}

void LoadDensityMapMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	filters::Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/ )
{
	mapfile_ = tag->getOption<std::string>("mapfile");
	sc_scale_ = tag->getOption<core::Real>("sc_scale", 1.0);
	window_ = tag->getOption<core::Size>("window", 3);
}

std::string LoadDensityMapMover::get_name() const {
	return mover_name();
}

std::string LoadDensityMapMover::mover_name() {
	return "LoadDensityMap";
}

void LoadDensityMapMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute( "mapfile", xs_string, "path to density map" );
	attlist + XMLSchemaAttribute::attribute_w_default( "sc_scale", xsct_real, "set sidechain scaling in density map to this factor", "1.0");
	attlist + XMLSchemaAttribute::attribute_w_default( "window", xsct_non_negative_integer, "set window in the density map to this value", "3");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "PHENIX crystallographic refinement interface code. Loads an electron density map", attlist );
}

std::string LoadDensityMapMoverCreator::keyname() const {
	return LoadDensityMapMover::mover_name();
}

protocols::moves::MoverOP
LoadDensityMapMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoadDensityMapMover );
}

void LoadDensityMapMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LoadDensityMapMover::provide_xml_schema( xsd );
}


////

void FitBfactorsMover::apply( core::pose::Pose & pose ) {
	if ( adp_strategy_ == "randomize" ) {
		randomize_bs( pose );
	} else {
		core::scoring::cryst::getPhenixInterface().set_adp_strategy( adp_strategy_ );
		core::scoring::cryst::getPhenixInterface().fitBfactors( pose );
	}

	// reinitialize fmodel object by rescoring pose
	core::scoring::cryst::getPhenixInterface().getScore( pose );
}

void FitBfactorsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	filters::Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/ )
{
	std::string adp_strat = tag->getOption<std::string>("adp_strategy", "individual");
	runtime_assert( adp_strat=="individual" || adp_strat=="group" || adp_strat=="randomize" );

	if ( adp_strat=="group" ) {
		auto group_adp_strat = tag->getOption<Size>("group_adp_mode", 1);
		runtime_assert( group_adp_strat==1 || group_adp_strat==2 );
		if ( group_adp_strat==1 ) {
			adp_strategy_ = "group1";
		} else {
			adp_strategy_ = "group2";
		}
	} else if ( adp_strat=="randomize" ) {
		adp_strategy_ = "randomize";
		b_min_ = tag->getOption<core::Real>("b_min", 5.0);
		b_max_ = tag->getOption<core::Real>("b_max", 5.0);
	} else {
		adp_strategy_ = "individual";
	}
}

void FitBfactorsMover::randomize_bs( core::pose::Pose & pose ) {
	for ( core::Size resid=1; resid<=pose.size(); ++resid ) {
		if ( pose.residue(resid).aa() == core::chemical::aa_vrt ) continue;
		core::conformation::Residue const &rsd_i = pose.residue(resid);
		for ( core::Size atmid=1; atmid<=rsd_i.natoms(); ++atmid ) {
			core::Real B = numeric::random::uniform()*(b_max_-b_min_) + b_min_;
			pose.pdb_info()->temperature( resid, atmid, B);
		}
	}
}

std::string FitBfactorsMover::get_name() const {
	return mover_name();
}

std::string FitBfactorsMover::mover_name() {
	return "FitBfactors";
}

void FitBfactorsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//argument legality checker for "adp_strategy":
	XMLSchemaRestriction adp_strategy_type;
	adp_strategy_type.name("adp_strategy_type");
	adp_strategy_type.base_type( xs_string );
	adp_strategy_type.add_restriction( xsr_enumeration, "individual");
	adp_strategy_type.add_restriction( xsr_enumeration, "group");
	adp_strategy_type.add_restriction( xsr_enumeration, "randomize");
	xsd.add_top_level_element( adp_strategy_type );

	attlist + XMLSchemaAttribute::attribute_w_default( "adp_strategy", "adp_strategy_type", "ADP (?) strategy.  Must be one of ('individual', 'group', 'randomize')", "individual" );


	//argument legality checker for "adp_strategy":
	XMLSchemaRestriction group_adp_mode;
	group_adp_mode.name("group_adp_mode_restricted_string_type");
	group_adp_mode.base_type( xsct_non_negative_integer );
	group_adp_mode.add_restriction( xsr_enumeration, "1");
	group_adp_mode.add_restriction( xsr_enumeration, "2");
	xsd.add_top_level_element( group_adp_mode );

	//XML XRW TODO: this option only relevant if adp_strategy == group; does not crash if not true but unchecked
	attlist + XMLSchemaAttribute::attribute_w_default( "group_adp_mode", "group_adp_mode_restricted_string_type", "group_adp_mode strategy, legal arguments are 1 or 2", "1" );

	//XML XRW TODO: these options are only checked if adp_strategy == random; does not crash if not true but unchecked
	attlist + XMLSchemaAttribute::attribute_w_default( "b_min", xsct_real, "b value min; used only with randomize strategy", "5.0" );
	attlist + XMLSchemaAttribute::attribute_w_default( "b_max", xsct_real, "b value max; used only with randomize strategy", "5.0" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "PHENIX crystallographic refinement interface code.  Fits b factors using an adp strategy", attlist );
}

std::string FitBfactorsMoverCreator::keyname() const {
	return FitBfactorsMover::mover_name();
}

protocols::moves::MoverOP
FitBfactorsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new FitBfactorsMover );
}

void FitBfactorsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FitBfactorsMover::provide_xml_schema( xsd );
}


////

void UpdateSolventMover::apply( core::pose::Pose & pose ) {
	// make sure fmodel is initialized by scoring the pose
	core::scoring::cryst::getPhenixInterface().getScore( pose );

	//fcalc
	if ( update_fcalc_ ) {
		core::scoring::cryst::getPhenixInterface().updateFcalc();
	}

	//mask
	if ( optimize_mask_ && optimize_params_ ) {
		core::scoring::cryst::getPhenixInterface().optimizeSolvParamsAndMask();
	} else if ( optimize_mask_ ) {
		core::scoring::cryst::getPhenixInterface().optimizeSolventMask();
	} else if ( optimize_params_ ) {
		core::scoring::cryst::getPhenixInterface().optimizeSolvParams();
	} else if ( update_mask_ ) {
		core::scoring::cryst::getPhenixInterface().updateSolventMask();
	}
}

void UpdateSolventMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	filters::Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/ )
{
	update_mask_ = tag->getOption<bool>("update_mask", true);
	update_fcalc_ = tag->getOption<bool>("update_fcalc", true);
	optimize_mask_ = tag->getOption<bool>("optimize_mask", false);
	optimize_params_ = tag->getOption<bool>("optimize_params", false);
}

std::string UpdateSolventMover::get_name() const {
	return mover_name();
}

std::string UpdateSolventMover::mover_name() {
	return "UpdateSolvent";
}

void UpdateSolventMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//XML XRW TODO: encode this dependency?
	attlist + XMLSchemaAttribute::attribute_w_default( "update_mask", xsct_rosetta_bool, "call PHENIX to update solvent mask. Only valid if optimize_mask and optimize_params are both false.", "true" );

	attlist + XMLSchemaAttribute::attribute_w_default( "update_fcalc", xsct_rosetta_bool, "call PHENIX to update fcalc", "true" );

	attlist + XMLSchemaAttribute::attribute_w_default( "optimize_mask", xsct_rosetta_bool, "call PHENIX to optimize solvent mask; can be combined with optimize_params", "false" );

	attlist + XMLSchemaAttribute::attribute_w_default( "optimize_params", xsct_rosetta_bool, "call PHENIX to optimize solvent mask; can be combined with optimize_mask", "false" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "PHENIX crystallographic refinement interface code.  Calls PHENIX to update and optimize solvent mask.", attlist );
}

std::string UpdateSolventMoverCreator::keyname() const {
	return UpdateSolventMover::mover_name();
}

protocols::moves::MoverOP
UpdateSolventMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new UpdateSolventMover );
}

void UpdateSolventMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	UpdateSolventMover::provide_xml_schema( xsd );
}


////

// helper function
void get_corresponding_CAs(
	core::pose::Pose const & model,
	core::pose::Pose const & native,
	core::sequence::SequenceAlignment const & aln,
	ObjexxFCL::FArray2D< core::Real > & p1a,
	ObjexxFCL::FArray2D< core::Real > & p2a
) {
	core::id::SequenceMapping mapping( aln.sequence_mapping(1,2) );

	core::Size natoms = 0;
	for ( core::Size ii = 1; ii <= model.size(); ++ii ) {
		Size const native_ii( mapping[ii] );
		if ( native_ii != 0 && native_ii <= native.size() && model.residue(ii).is_protein() ) {
			natoms++;
		}
	}
	p1a.dimension(3,natoms); p2a.dimension(3,natoms);

	Size n_gap(0);
	for ( core::Size ii = 1; ii <= model.size(); ++ii ) {
		Size const native_ii(mapping[ii]);
		if ( native_ii != 0 && native_ii <= native.size() && model.residue(ii).is_protein() ) {
			numeric::xyzVector< core::Real > model_xyz ( model.residue(ii).xyz("CA") );
			numeric::xyzVector< core::Real > native_xyz( native.residue(native_ii).xyz("CA") );
			for ( core::Size jj = 1; jj <= 3; ++jj ) {
				p1a(jj,ii - n_gap) = native_xyz[jj-1];
				p2a(jj,ii - n_gap) = model_xyz [jj-1];
			}
		} else {
			n_gap++;
		}
	}
}


void TagPoseWithRefinementStatsMover::apply( core::pose::Pose & pose ) {
	// make sure fmodel is initialized by scoring the pose
	/*core::Real xraytgt =*/ core::scoring::cryst::getPhenixInterface().getScore( pose );

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	//fpd if the pose is symmetric use a symmetric scorefunction
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *scorefxn );
	}
	core::Real score = (*scorefxn)(pose);
	core::io::RemarkInfo remark;

	std::ostringstream oss;

	oss << "[ " << tag_ << " ] R=" << core::scoring::cryst::getPhenixInterface().getInfoLine();

	//oss.str("");
	//oss << "[ " << tag_ << " ] xray_tgt = " << xraytgt << "   score = " << score;
	oss << " sc=" << score;

	core::pose::Pose pose_asu;
	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() || dump_pose_ ) {
		core::pose::symmetry::extract_asymmetric_unit(pose, pose_asu);
	}

	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		// rms to native
		if ( !get_native_pose() ) protocols::jd2::set_native_in_mover(*this);
		core::pose::PoseCOP native = get_native_pose();

		core::sequence::SequenceOP model_seq( new core::sequence::Sequence( pose_asu.sequence(),  "model",  1 ) );
		core::sequence::SequenceOP native_seq( new core::sequence::Sequence( native->sequence(), "native", 1 ) );
		core::sequence::SequenceAlignment aln = align_naive(model_seq,native_seq);
		ObjexxFCL::FArray2D< core::Real > p1a, p2a;
		get_corresponding_CAs( pose_asu, *native, aln, p1a, p2a );
		core::Real rms = numeric::model_quality::rms_wrapper( p1a.u2(), p1a, p2a );
		oss << " rms=" << rms;
	}

	if ( report_grads_ ) {
		ReportGradientsMover RGG(scorefxn);
		core::Real gradsum = RGG.compute(pose);
		oss << " grad=" << gradsum;
	}

	TR << oss.str() << std::endl;
	remark.num = 1; remark.value = oss.str();
	pose.pdb_info()->remarks().push_back( remark );

	if ( dump_pose_ ) {
		std::string out_name = protocols::jd2::current_output_name() + "_" + tag_ + ".pdb";
		pose_asu.dump_pdb( out_name );
	}
}


void TagPoseWithRefinementStatsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	filters::Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/ )
{
	tag_ = tag->getOption<std::string>("tag", "");
	dump_pose_ = tag->getOption<bool>("dump", false);
	report_grads_ = tag->getOption<bool>("report_grads", false);
}

std::string TagPoseWithRefinementStatsMover::get_name() const {
	return mover_name();
}

std::string TagPoseWithRefinementStatsMover::mover_name() {
	return "TagPoseWithRefinementStats";
}

void TagPoseWithRefinementStatsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute( "tag", xs_string, "tags Mover's output with this string" );
	attlist + XMLSchemaAttribute::attribute_w_default( "dump", xsct_rosetta_bool, "dump a copy of the PDB", "false" );
	attlist + XMLSchemaAttribute::attribute_w_default( "report_grads", xsct_rosetta_bool, "report gradients sum", "false" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "PHENIX crystallographic refinement interface code.  Modifies Pose with refinement statistics and can dump the pdb", attlist );
}

std::string TagPoseWithRefinementStatsMoverCreator::keyname() const {
	return TagPoseWithRefinementStatsMover::mover_name();
}

protocols::moves::MoverOP
TagPoseWithRefinementStatsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new TagPoseWithRefinementStatsMover );
}

void TagPoseWithRefinementStatsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TagPoseWithRefinementStatsMover::provide_xml_schema( xsd );
}


////
void SetRefinementOptionsMover::apply( core::pose::Pose & /*pose*/ ) {
	TR << "Setting options:" << std::endl;

	if ( res_high_ != 0 || res_low_ != 0 ) {
		TR << "  res_high:" << res_high_ << std::endl;
		TR << "  res_low:" << res_low_ << std::endl;
		core::scoring::cryst::getPhenixInterface().setResLimits(res_high_, res_low_);
	}
	if ( twin_law_.length() != 0 ) {
		TR << "  twin_law:" << twin_law_ << std::endl;
		core::scoring::cryst::getPhenixInterface().setTwinLaw(twin_law_);
	}
	if ( algo_.length() != 0 ) {
		TR << "  algorithm:" << algo_ << std::endl;
		core::scoring::cryst::getPhenixInterface().setAlgorithm(algo_);
	}
	if ( target_.length() != 0 ) {
		TR << "  target:" << std::endl;
		core::scoring::cryst::getPhenixInterface().set_target_function(target_);
	}

	if ( setmap_type_ ) {
		TR << "  map_type:" << map_type_ << std::endl;
		core::scoring::cryst::getPhenixInterface().set_map_type ( map_type_ );
	}

	if ( cif_files_.size() > 0 ) {
		TR << "Passing cif files to phenix refine:" << std::endl;
		for ( core::Size i=1; i<=cif_files_.size(); ++i ) {
			TR << "  " << cif_files_[i] << std::endl;
		}
		core::scoring::cryst::getPhenixInterface().set_cif_files ( cif_files_ );
	}

	if ( sharpen_b_ != 0.0 ) {
		core::scoring::cryst::getPhenixInterface().set_sharpen_b ( sharpen_b_ );
	}
}


void SetRefinementOptionsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	filters::Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/ )
{
	res_high_ = tag->getOption<core::Real>("res_high", 0.0);
	res_low_ = tag->getOption<core::Real>("res_low", 0.0);
	twin_law_ = tag->getOption<std::string>("twin_law", "");
	algo_ = tag->getOption<std::string>("algorithm", "");
	target_ = tag->getOption<std::string>("target", "");
	setmap_type_ = tag->hasOption("map_type");
	map_type_ = tag->getOption<std::string>("map_type", "");
	sharpen_b_ = tag->getOption<core::Real>("sharpen_b", 0.0);

	std::string allcifs = tag->getOption<std::string>("cifs", "");
	if ( allcifs != "" ) {
		cif_files_ = utility::string_split( allcifs, ',' );
	}
}

std::string SetRefinementOptionsMover::get_name() const {
	return mover_name();
}

std::string SetRefinementOptionsMover::mover_name() {
	return "SetRefinementOptions";
}

void SetRefinementOptionsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default( "res_high", xsct_real, "sets high ResLimit for PHENIX refinement", "0.0" );
	attlist + XMLSchemaAttribute::attribute_w_default( "res_low", xsct_real, "sets low ResLimit for PHENIX refinement", "0.0" );
	attlist + XMLSchemaAttribute( "twin_law", xs_string, "sets twin law for PHENIX refinement" );
	attlist + XMLSchemaAttribute( "algorithm", xs_string, "sets algorithm for PHENIX refinement" );
	attlist + XMLSchemaAttribute( "target", xs_string, "sets target function for PHENIX refinement" );
	attlist + XMLSchemaAttribute( "map_type", xs_string, "sets map type for PHENIX refinement" );
	attlist + XMLSchemaAttribute::attribute_w_default( "sharpen_b", xsct_real, "sets sharpen_b for PHENIX refinement", "0.0" );

	//XML XRW NOTE: cannot enforce "comma separated strings" with XSD restrictions
	attlist + XMLSchemaAttribute( "cifs", xs_string, "cif files passed along to PHENIX refinement.  Legal argument is cifs split by commas: 'cif1.cif,cif2.cif,cif3.cif" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "PHENIX crystallographic refinement interface code.  Sets refinement options to pass along through the PHENIX interface.", attlist );
}

std::string SetRefinementOptionsMoverCreator::keyname() const {
	return SetRefinementOptionsMover::mover_name();
}

protocols::moves::MoverOP
SetRefinementOptionsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetRefinementOptionsMover );
}

void SetRefinementOptionsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SetRefinementOptionsMover::provide_xml_schema( xsd );
}


}
}
