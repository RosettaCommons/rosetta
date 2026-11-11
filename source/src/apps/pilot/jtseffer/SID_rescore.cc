// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/jtseffer/SID_rescore
/// @brief Application to rescore poses from RosettaDock for use with SID.
/// @author Justin Seffernick (seffernick.9@osu.edu)

// In Surface-Induced Dissociation (SID), protein complexes are ionized in the gas phase and
// acelerated towards a rigid surface. Depending on the voltage applied (and thus acceleration
// energy of the complex), complexes break apart differently (upon surface collision), allowing
// the determination of appearance energy (AE) using mass spectrometry. AE is defined as the
// acceleration energy needed for 10% subcomplex formation. See Zhou, M.; Wysocki, V.H., Surface
// induced dissociation: Dissecting noncovalent protein complexes in the gas phase.
// Acc. Chem. Res. 2014; 47(4), 1010-1018

// To use this application, input structures from RosettaDock and SID AE.

#include <basic/options/option.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <utility/io/ozstream.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.fwd.hh>
#include <protocols/pose_metric_calculators/SaltBridgeCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>
#include <protocols/relax/FastRelax.fwd.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/scoring/rms_util.hh>

static basic::Tracer TR( "apps.pilot.jtseffer.SID_rescore" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

//local options
basic::options::StringOptionKey const interface( "interface" );
basic::options::RealOptionKey const AE( "AE" );
basic::options::RealOptionKey const n_ints( "n_ints" );
basic::options::BooleanOptionKey const skip_relax( "skip_relax" );

void read_in_pdbs( utility::vector1<core::pose::PoseOP> &poses, core::Size &num_chains ) {
	core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();
	//Import each pose
	while ( input.has_another_pose() ) {
		core::pose::PoseOP mypose( new core::pose::Pose );
		input.fill_pose( *mypose );
		poses.push_back(mypose);
	}
	num_chains = poses[1]->conformation().num_chains();
	TR << "Number of poses entered: " << num_chains << std::endl;
}

void calc_RF( const utility::vector1<core::pose::PoseOP> &pose_by_chain, core::Real &RF) {
	//Register sb metric calculator
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sb_metric", core::pose::metrics::PoseMetricCalculatorOP(new protocols::pose_metric_calculators::SaltBridgeCalculator));

	RF = 0.0;

	for ( core::Size i=1; i<=pose_by_chain.size(); i++ ) {

		core::pose::PoseOP pose = pose_by_chain[i];

		core::scoring::hbonds::HBondSet myhbset = core::scoring::hbonds::HBondSet(*pose);

		core::Real num_res = pose->size();

		//calculate salt bridges
		std::string SB_string = pose->print_metric("sb_metric", "salt_bridge");
		core::Real SB;
		std::stringstream convert_ss( SB_string );
		convert_ss >> SB;
		core::Real SB_norm = SB / num_res;

		//calculate hydrogen bonds
		core::Size HB = myhbset.nhbonds();
		core::Real HB_norm = HB / num_res;

		//calculate disulfides
		core::Size DS = 0;
		for ( core::Size i=1; i <= num_res; ++i ) {
			for ( core::Size j=1; j <= num_res; ++j ) {
				if ( pose->residue(i).is_bonded(j) && ! pose->residue(i).is_polymer_bonded(j) ) {
					DS++;
				}
			}
		}
		DS = DS / 2.0;
		core::Real DS_norm = DS / num_res;

		core::Real E_intra = 2.5*SB_norm + HB_norm + 60*DS_norm;

		if ( E_intra <= 0.41 ) {
			RF +=  0.0;
		} else if ( E_intra < 1.09 ) {
			RF += 1.47*E_intra - 0.603;
		} else {
			RF += 1.00;
		}
	}

	core::Real num_chains = pose_by_chain.size();

	RF = RF / num_chains;

}

void predict_AE( const core::pose::PoseOP &pose, core::Real &AE_pred, const core::pose::DockingPartners &interface, const core::Real &HSA_weight, const core::Real &NR_weight, const core::Real &RF_weight, const core::Real &RF, const core::Real &n_ints ) {
	protocols::analysis::InterfaceAnalyzerMoverOP IAM = utility::pointer::make_shared< protocols::analysis::InterfaceAnalyzerMover >(interface, true, core::scoring::ScoreFunctionFactory::create_score_function("ref2015")/*scorefxn_*/, false/*compute_packstat_*/, false/*pack_together_*/, false/*pack_separated_*/);

	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();

	bool do_relax = true;

	if ( option[ skip_relax ].user() ) {
		do_relax = false;
	}

	if ( do_relax ) {
		protocols::relax::FastRelaxOP relax( new protocols::relax::FastRelax( sfxn ) );
		relax->apply(*pose);
	}

	IAM->apply(*pose);

	core::Size NR = IAM->get_num_interface_residues();
	protocols::analysis::InterfaceData int_data = IAM->get_all_data();
	utility::vector1< core::Real > dhSASA_vec = int_data.dhSASA;
	core::Real dhSASA = dhSASA_vec[1];

	AE_pred = HSA_weight*dhSASA/n_ints + NR_weight*NR/n_ints - RF_weight*RF;
}

void calc_SID_score(core::Real &SID_score, const core::Real &AE_pred, const core::Real &AE_exp) {
	core::Real diff = std::abs(AE_pred - AE_exp);
	core::Real ub = 1750;
	core::Real lb = 100;
	if ( diff <=lb ) {
		SID_score = 0;
	} else if ( diff < ub ) {
		core::Real b = -(diff - ub) / (ub-lb);
		SID_score = 100*( 2*pow(b,3) -3*pow(b,2) +1 );
	} else {
		SID_score = 100.0;
	}
}

void output_results ( const utility::vector1<core::Real> &AE_preds, const utility::vector1<core::Real> &SID_scores, const utility::vector1<core::Real> &Rosetta_scores, const utility::vector1<core::Real> &Rosetta_SID_scores, const utility::vector1<core::Real> &rmsds ) {
	std::ostringstream out;
	out << "Pose_number\tAE_pred\t\tRosetta_score\tSID_score\tRosetta_SID_score\tRMSD" << std::endl;
	for ( core::Size i=1; i<=AE_preds.size(); i++ ) {
		out << i << "\t\t" << AE_preds[i] << "\t" << Rosetta_scores[i] << "\t" << SID_scores[i] << "\t" << Rosetta_SID_scores[i] << "\t" << rmsds[i] << std::endl;
	}
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::string outfile = "SID_rescore_default.out";
	if ( option[ out::file::o ].user() ) outfile = option[ out::file::o ]();
	utility::io::ozstream outz( outfile.c_str() );
	outz << out.str();
	outz.close();
	outz.clear();
}

int
main( int argc, char * argv [] )
{

	try{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::scoring;
		option.add( interface, "interface definition" ).def("A_B");
		option.add( AE, "Experimental AE from SID (eV)" ).def(0.0);
		option.add( n_ints, "Number of interfaes" ).def(1.0);
		option.add( skip_relax, "Skip relax step" ).def(false);
		devel::init( argc, argv );

		//Read in variables from command line
		core::Real AE_exp( option[ AE ].value() );
		core::Real ints( option[ n_ints ].value() );
		AE_exp = AE_exp / ints;
		core::pose::DockingPartners intf = core::pose::DockingPartners::docking_partners_from_string( option[ interface ].value() );

		//SID weights
		core::Real HSA_weight = 0.121422027;
		core::Real NR_weight = 5.15236;
		core::Real RF_weight = 208.741;

		//Import poses
		core::Size num_chains;
		utility::vector1<core::pose::PoseOP> poses;
		read_in_pdbs( poses, num_chains );
		TR << "number of chains: " << num_chains << std::endl;

		//Import native
		core::pose::Pose native;
		if ( option[ in::file::native ].user() ) {
			core::import_pose::pose_from_file( native, option[ basic::options::OptionKeys::in::file::native ]);
		}


		//Calculate Rigidity Factor (RF) of first pose
		//split first pose by chains
		utility::vector1<core::pose::PoseOP> pose1_by_chain;
		for ( core::Size i=1; i<=num_chains; i++ ) {
			pose1_by_chain.push_back(poses[1]->split_by_chain(i));
		}
		core::Real RF;
		calc_RF(pose1_by_chain, RF);
		TR << "RF: " << RF << std::endl;

		//Calculate Scores for each pose
		utility::vector1<core::Real> AE_preds;
		utility::vector1<core::Real> SID_scores;
		utility::vector1<core::Real> Rosetta_scores;
		utility::vector1<core::Real> Rosetta_SID_scores;
		utility::vector1<core::Real> rmsds;
		for ( core::Size i=1; i<=poses.size(); i++ ) {
			//Calculate Rosetta score
			core::scoring::ScoreFunctionOP sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "docking.wts" );
			core::Real Rosetta_score = (*sfxn)(*poses[i]);
			Rosetta_score = Rosetta_score/ints;
			Rosetta_scores.push_back(Rosetta_score);
			//Calculate predicted SID AE
			core::Real AE_pred;
			predict_AE(poses[i], AE_pred, intf, HSA_weight, NR_weight, RF_weight, RF, ints);
			AE_preds.push_back(AE_pred);
			//Calculate SID score
			core::Real SID_score;
			calc_SID_score(SID_score, AE_pred, AE_exp);
			SID_scores.push_back(SID_score);
			//Calculate Combined Score
			core::Real Rosetta_SID_score = Rosetta_score + 6.0*SID_score;
			Rosetta_SID_scores.push_back(Rosetta_SID_score);
			//Calculate RMSD
			core::Real rmsd;
			if ( option[ in::file::native ].user() ) {
				rmsd = core::scoring::CA_rmsd( *poses[i], native );
			} else {
				rmsd = 0.0;
			}
			rmsds.push_back(rmsd);

			//print results
			TR << "Pose number: " << i << std::endl;
			TR << "Predicted AE: " << AE_pred << std::endl;
			TR << "SID score: " << SID_score << std::endl;
			TR << "Rosetta score: " << Rosetta_score << std::endl;
			TR << "Rosetta/SID combined score: " << Rosetta_SID_score << std::endl;
			TR << "RMSD: " << rmsd << std::endl;
			output_results(AE_preds, SID_scores, Rosetta_scores, Rosetta_SID_scores, rmsds);
		}

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
