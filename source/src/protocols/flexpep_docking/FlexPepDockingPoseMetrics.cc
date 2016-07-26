// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   FlexPepDockingPoseMetrics.cc
///
/// @brief metrics calculations specific for FlexPepDock (at least for now)
/// @date March 29th, 2009
/// @author Barak Raveh / Nir London

#include <protocols/flexpep_docking/FlexPepDockingFlags.hh>
#include <protocols/flexpep_docking/FlexPepDockingPoseMetrics.hh>

#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <core/types.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>

// Utility headers
#include <utility/file/FileName.hh>
#include <utility/assert.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.hh>
#include <limits>
#include <string>

#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;
using core::pose::Pose;
using namespace core;

static THREAD_LOCAL basic::Tracer TR( "FlexPepDockingPoseMetrics" );

namespace protocols {
namespace flexpep_docking {

//@brief fraction of native contacts modeled in final pose (relative to native)
//
//@param  native    - native pose
//@param  final     - final model
//@param  threshold - threshold in Angstroms for native contacts (two
//                    (residues define a native contact if the distance
//                     between two of their atoms is lower than threshold)
core::Real
FlexPepDockingPoseMetrics::calc_frac_native_contacts(
	Pose const& native,
	Pose const& final,
	core::Real threashold
) const
{
	core::Size numOfNativeContacts = 0;
	core::Size subsetOfPredictedContacts = 0;

	if ( flags_->pep_fold_only ) {
		TR << "WARNING: calc_frac_native_contacts() should have not been invoked when -pep_fold_only flag (= no receptor) is active";
		return 0.0;
	}
	//iterate over the peptide
	for ( int i=flags_->peptide_first_res(); i<=flags_->peptide_last_res(); ++i ) {
		//iterate over the protein
		for ( int j=flags_->receptor_first_res(); j<=flags_->receptor_last_res(); ++j ) {
			if ( isInContact(native.residue(i),native.residue(j),threashold) ) {
				numOfNativeContacts++;
				if ( isInContact(final.residue(i),final.residue(j),threashold) ) {
					subsetOfPredictedContacts++;
				}
			}
		}
	}
	TR << "nat " <<  numOfNativeContacts << std::endl;
	TR << "rec " <<  subsetOfPredictedContacts << std::endl;
	core::Real fnat = ((core::Real)subsetOfPredictedContacts/(core::Real)numOfNativeContacts);
	return fnat;
}

bool FlexPepDockingPoseMetrics::isInContact(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2,
	core::Real threashold
) const
{
	core::Size natoms_res1 = res1.natoms();
	core::Size natoms_res2 = res2.natoms();
	for ( core::Size i = 1; i <= natoms_res1; ++i ) {
		for ( core::Size j = 1; j <= natoms_res2; ++j ) {
			core::Distance dist = res1.xyz(i).distance( res2.xyz(j) );
			if ( dist <= threashold ) {
				return true;
			}
		}
	}
	// in case no two atoms are below the threashold
	return false;
}


core::Real
FlexPepDockingPoseMetrics::calc_frac_atoms_kA_to_native(
	Pose const& pose1,
	Pose const & pose2,
	ObjexxFCL::FArray1D_bool const & res_subset,
	t_predicate_func predicate,
	double k,
	core::Size & ngood
) const
{
	using namespace core::scoring;

	core::Size const nres1 = pose1.total_residue();
	ASSERT_ONLY(core::Size const nres2 = pose2.total_residue();)
		assert( nres1 == nres2 );

	ngood = 0;
	core::Size natoms_total( 0 );
	for ( core::Size i = 1; i <= nres1; ++i ) {
		if ( !res_subset(i) ) {
			continue;
		}
		// Calculate phi psi RMSd only for the positions that actually have these properties defined
		core::chemical::ResidueType residue_type = pose1.conformation().residue_type(i);
		if ( !residue_type.is_protein() && !residue_type.is_peptoid() && !residue_type.is_carbohydrate() ) {
			continue;
		}
		core::Size natoms_res ( pose1.residue(i).natoms() );
		if ( predicate == &is_ligand_heavyatom ||
				predicate == &is_polymer_heavyatom ||
				predicate == &is_heavyatom ) {
			assert( pose1.residue(i).natoms() == pose2.residue(i).natoms() );
		} else if ( natoms_res > pose2.residue(i).natoms() ) {
			natoms_res = pose2.residue(i).natoms();
		}
		for ( core::Size j = 1; j <= natoms_res; ++j ) {
			if ( predicate( pose1, pose2, i, j ) ) {
				core::Distance dist = pose1.residue(i).xyz(j).distance( pose2.residue(i).xyz(j) );
				if ( dist <= k ) {
					ngood++;
				}
				natoms_total += 1;
			}
		}
	}
	return (ngood * 1.0) / natoms_total;
}


// check all sequential Kmers in the peptide, and output the best RMS Kmer
core::Real
FlexPepDockingPoseMetrics::best_Kmer_rms
( Pose const& pose1,
	Pose const& pose2,
	t_predicate_func predicate,
	core::Size k) const
{
	using namespace core;
	assert(pose1.total_residue() == pose2.total_residue());
	// NOTE: purposely an inefficient but simpler construction...
	// It would be more efficient to calc the RMS of each position separately
	// and make the calculation incrementaly, but why bother (Barak)
	Size nres = pose1.total_residue();
	Size res1_first = flags_->peptide_first_res()  ;
	Size res1_last = nres - k + 1;
	Real best_rms = std::numeric_limits<Real>::infinity();
	for ( Size res1 = res1_first ; res1 <= res1_last ; ++res1 ) {
		ObjexxFCL::FArray1D_bool res_subset( nres, false );
		for ( Size resi = res1; resi < res1 + k; ++resi ) {
			res_subset(resi) = true;
		}
		Real cur_rms = scoring::rmsd_no_super_subset( pose1, pose2, res_subset, *predicate );
		if ( cur_rms < best_rms ) {
			best_rms = cur_rms;
		}
	}
	return best_rms;
}


/////////////////////////////////////////////////////////////////////////////
///
//  @details
//  Calculate the RMSD in phi.psi angle over specified residues
//
//  @param
//  pose1 - the first structure to be assessed
//  pose2 - the second structure to be assessed
//  res_subset - an array of size (nres) indicating the residue
//               subset for the computation
//
// @return
// phi/psi torsion-RMSD between peptide backbones
////////////////////////////////////////////////////////////////////////////
core::Real
FlexPepDockingPoseMetrics::calc_phipsi_RMSD
( Pose const& pose1,
	Pose const& pose2,
	ObjexxFCL::FArray1D_bool const & res_subset) const
{
	using namespace core;
	using namespace basic;
	assert(pose1.total_residue() == pose2.total_residue());
	Size nres = pose1.total_residue();
	Real sumSqr = 0.0; // MSD = sum square deviation
	Size n = 0;
	for ( Size i = 1 ; i <= nres ; ++i ) {
		if ( !res_subset(i) ) {
			continue;
		}
		Real phidiff = subtract_degree_angles(pose1.phi(i), pose2.phi(i));
		Real psidiff = subtract_degree_angles(pose1.psi(i), pose2.psi(i));
		if ( !pose1.residue_type(i).is_lower_terminus() ) {
			sumSqr += pow(phidiff,2);
			n++;
		}
		if ( !pose1.residue_type(i).is_upper_terminus() ) {
			sumSqr += pow(psidiff,2);
			n++;
		}
	}
	return sqrt( sumSqr / n);

}


/////////////////////////////////////////
// Calculate pose metrics for interface
// and cache them to the score map
// calculates i_sc as well
/////////////////////////////////////////
std::map < std::string, core::Real >
FlexPepDockingPoseMetrics::calc_interface_metrics( core::pose::Pose & pose, Size rb_jump, core::scoring::ScoreFunctionOP scorefxn ) {
	using namespace core;
	using namespace core::pose::metrics;
	using namespace core::scoring;

	std::map < std::string, core::Real > if_metrics;

	core::scoring::ScoreFunctionOP scorefxn_no_cst = scorefxn->clone();
	scorefxn_no_cst->set_weight(coordinate_constraint, 0.0);
	scorefxn_no_cst->set_weight(atom_pair_constraint, 0.0);
	scorefxn_no_cst->set_weight(angle_constraint, 0.0);
	scorefxn_no_cst->set_weight(dihedral_constraint, 0.0);

	// calculate score, this is essential just because of bug in h-bond metrics in pose - TODO: fix this bug //
	(*scorefxn_no_cst)(pose);

	//////// calculate interface score ////////
	core::pose::Pose unbound_pose = pose;
	float trans_magnitude = 1000;
	rigid::RigidBodyTransMoverOP translate_away( new rigid::RigidBodyTransMover( unbound_pose, rb_jump ) );
	translate_away->step_size( trans_magnitude );

	float bound_energy = (*scorefxn_no_cst)( unbound_pose );
	translate_away->apply( unbound_pose );
	float unbound_energy = (*scorefxn_no_cst)( unbound_pose );

	float Isc = (bound_energy - unbound_energy);
	TR << "Isc: " << Isc << std::endl;
	///////////////////////////////////////////

	std::string sasa_calc_name = "sasa";
	std::string hbond_calc_name = "hbond";
	std::string packstat_calc_name = "packstat";
	std::string burunsat_calc_name = "burunsat";
	////std::string nonlocalcontacts_calc_name = "conts";

	// Register calculators (if they do not exist)

	// sasa
	if ( !CalculatorFactory::Instance().check_calculator_exists( sasa_calc_name ) ) {
		PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy );
		CalculatorFactory::Instance().register_calculator
			( sasa_calc_name, sasa_calculator );
	}
	// hbonds
	if ( !CalculatorFactory::Instance().check_calculator_exists( hbond_calc_name ) ) {
		core::pose::metrics::PoseMetricCalculatorOP hb_calc( new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator() );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator
			( hbond_calc_name, hb_calc );
	}
	// packstats
	if ( !CalculatorFactory::Instance().check_calculator_exists( packstat_calc_name ) ) {
		core::pose::metrics::PoseMetricCalculatorOP packstat_calc( new protocols::toolbox::pose_metric_calculators::PackstatCalculator() );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator
			( packstat_calc_name, packstat_calc );
	}
	// burried unsatisfied polar
	if ( !CalculatorFactory::Instance().check_calculator_exists( burunsat_calc_name ) ) {
		core::pose::metrics::PoseMetricCalculatorOP burunsat_calc( new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator
			(sasa_calc_name, hbond_calc_name) );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator
			( burunsat_calc_name, burunsat_calc );
	}

	// define containers for metrics for total complex
	basic::MetricValue<Real> tot_sasa_mval;
	basic::MetricValue<Size> tot_hb_mval;
	basic::MetricValue<Real> tot_packstat_mval;
	basic::MetricValue<Size> tot_unsat_mval;

	// define containers for metrics per peptide residue
	basic::MetricValue< utility::vector1< core::Real > >bound_sasa_per_res_mval;
	basic::MetricValue< utility::vector1< core::Size > >bound_hb_per_res_mval;
	basic::MetricValue< utility::vector1< core::Real > >bound_packstat_per_res_mval;
	basic::MetricValue< utility::vector1< core::Size > >bound_unsat_per_res_mval;

	basic::MetricValue< utility::vector1< core::Real > >unbound_sasa_per_res_mval;
	basic::MetricValue< utility::vector1< core::Size > >unbound_hb_per_res_mval;
	basic::MetricValue< utility::vector1< core::Real > >unbound_packstat_per_res_mval;
	basic::MetricValue< utility::vector1< core::Size > >unbound_unsat_per_res_mval;

	// calculate and store total metrics for bound and unbound poses
	float bound_sasa/* = 0.0*/, unbound_sasa;// = 0.0;
	Size  bound_hb/* = 0*/,   unbound_hb;// = 0;
	float bound_packstat = 0.0, unbound_packstat = 0.0;
	Size  bound_unsat/* = 0*/, unbound_unsat;// = 0;

	//delta sasa calculation
	pose.metric(sasa_calc_name,"total_sasa",tot_sasa_mval);
	bound_sasa = tot_sasa_mval.value();
	unbound_pose.metric(sasa_calc_name,"total_sasa",tot_sasa_mval);
	unbound_sasa = tot_sasa_mval.value();
	TR << "Total BSA  is: " << unbound_sasa - bound_sasa << std::endl;

	//interface hb calculation
	pose.metric(hbond_calc_name,"all_Hbonds", tot_hb_mval);
	bound_hb = tot_hb_mval.value();
	unbound_pose.metric(hbond_calc_name,"all_Hbonds", tot_hb_mval);
	unbound_hb = tot_hb_mval.value();
	TR << "Interface HB #: " << bound_hb - unbound_hb << std::endl;

	if ( !flags_->is_ligand_present(pose) ) {
		//packstat calculation
		pose.metric(packstat_calc_name,"total_packstat", tot_packstat_mval);
		bound_packstat = tot_packstat_mval.value();
		unbound_pose.metric(packstat_calc_name,"total_packstat", tot_packstat_mval);
		unbound_packstat = tot_packstat_mval.value();
		TR << "Total packstats: " << bound_packstat - unbound_packstat << std::endl;
	}

	//unsat polar calculation
	pose.metric(burunsat_calc_name,"all_bur_unsat_polars", tot_unsat_mval);
	bound_unsat = tot_unsat_mval.value();
	unbound_pose.metric(burunsat_calc_name,"all_bur_unsat_polars", tot_unsat_mval);
	unbound_unsat = tot_unsat_mval.value();
	TR << "Interface Unsat polar groups: " << bound_unsat - unbound_unsat << std::endl;

	//update answer:
	if_metrics["I_sc"] = Isc;
	if_metrics["I_bsa"] = unbound_sasa - bound_sasa;
	if_metrics["I_hb"] = bound_hb - unbound_hb;
	if_metrics["I_pack"] = bound_packstat - unbound_packstat;
	if_metrics["I_unsat"] = bound_unsat - unbound_unsat;

	///////////////// Per peptide residue calculations ///////////////////
	pose.metric(sasa_calc_name, "residue_sasa"  ,bound_sasa_per_res_mval);
	pose.metric(hbond_calc_name,"residue_Hbonds",bound_hb_per_res_mval);
	if ( !flags_->is_ligand_present(pose) ) {
		pose.metric(packstat_calc_name,"residue_packstat",bound_packstat_per_res_mval);
	}
	pose.metric(burunsat_calc_name,"residue_bur_unsat_polars",bound_unsat_per_res_mval);

	unbound_pose.metric(sasa_calc_name, "residue_sasa"  ,unbound_sasa_per_res_mval);
	unbound_pose.metric(hbond_calc_name,"residue_Hbonds",unbound_hb_per_res_mval);
	unbound_pose.metric(burunsat_calc_name,"residue_bur_unsat_polars",unbound_unsat_per_res_mval);

	for ( int i=flags_->peptide_first_res();
			i < flags_->peptide_first_res() + flags_->peptide_nres(); i++ ) {
		TR << i
			<< " bsa: " << unbound_sasa_per_res_mval.value()[i] - bound_sasa_per_res_mval.value()[i]
			<< " HB: " << bound_hb_per_res_mval.value()[i] - unbound_hb_per_res_mval.value()[i];
		if ( !flags_->is_ligand_present(pose) ) {
			TR << " pack: " << bound_packstat_per_res_mval.value()[i]; //- unbound_packstat_per_res_mval.value()[i]
		}
		TR << " unsat: " << bound_unsat_per_res_mval.value()[i] - unbound_unsat_per_res_mval.value()[i]
			<< std::endl;
	}
	return if_metrics;
}

/////////////////////////////////////////
// Calculate peptide score with and w/o
// fa_ref sequence reference energy
/////////////////////////////////////////
void
FlexPepDockingPoseMetrics::calc_pep_scores
( core::pose::Pose const & pose, Real& pepScore, Real& pepScore_noref ) const {

	pepScore = 0.0;
	pepScore_noref = 0.0;

	for ( int i=flags_->peptide_first_res(); i <= flags_->peptide_last_res(); ++i ) {
		using namespace core::scoring;
		Real ienergy = pose.energies().residue_total_energy(i);
		Real ifa_ref = pose.energies().residue_total_energies(i)[ref];
		pepScore += ienergy;
		pepScore_noref += ienergy - ifa_ref;
	}
}


} // end namespace flexpep_docking
} // end namespace protocols
