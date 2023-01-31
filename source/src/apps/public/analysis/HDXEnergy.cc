// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    HDXEnergy.cc
/// @brief   Scores Ab Initio Models using HDX-NMR rate Data
/// @author  DMarzolf (marzolf.4@osu.edu)

// App headers
#include <devel/init.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/prof.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/sasa/SasaCalc.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>

// Package headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

// utility
#include <utility/string_util.hh>
#include <core/pose/util.hh>
#include <core/scoring/sasa/util.hh>
#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/vector1.hh>
#include <numeric/NumericTraits.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace pose;

OPT_KEY( Boolean, StrongRes )
OPT_KEY( Boolean, ResPF )

static basic::Tracer TR( "apps.pilot.marzolf-daniel.HDXEnergy" );


//////////////////////////////////////////////////////////////////////

void register_options() {
	option.add_relevant(in::file::HDX);
}

//Output Results
void output_results( std::string const & outf ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::string outfile = "HDXRescore.txt";
	if ( option[ out::file::o ].user() ) {
		outfile = option[ out::file::o ]();
	}
	std::ofstream outy;
	outy.open( outfile, std::ofstream::out | std::ofstream::app);
	outy << outf;
	TR << "Scoring Finished, find results at " << outfile << std::endl;
}
core::Real count_neighbors( std::string neighbor_oxygen, core::pose::Pose pose, core::Size res_count_target, std::string target_atom_vector_start,  core::Real const distance_internal, numeric::xyzVector<core::Real> target_vector, core::Size res_count_neighbor ) {
	numeric::xyzVector<core::Real> neighbor_O_vector = pose.residue(res_count_neighbor).xyz(neighbor_oxygen) - pose.residue(res_count_target).xyz(target_atom_vector_start);
	core::Real const distance_to_O = pose.residue(res_count_neighbor).xyz(neighbor_oxygen).distance(pose.residue(res_count_target).xyz(target_atom_vector_start));
	numeric::xyzVector<core::Real> norm_target_vector = target_vector/distance_internal;
	numeric::xyzVector<core::Real> norm_O_vector = neighbor_O_vector/distance_to_O;
	core::Real angleO = std::acos(norm_target_vector.dot(norm_O_vector));
	core::Real resi_nc = 1.0/(1.0 + std::exp((distance_to_O-9.0)))*1.0/(1.0 + std::exp(6.28318530718*(angleO-1.57079632679)));
	return resi_nc;
}

void calculate_average_residue_scores (core::pose::Pose pose, utility::vector1<core::Real> &average_res_scores, core::Size num_res) {
	for ( core::Size j=1; j <= num_res; j++ ) {
		//Score of the current residue
		core::Real res_score = pose.energies().residue_total_energy(j);
		//Append residue scores to the average score vector
		average_res_scores.push_back(res_score);
	}
}


void calculate_order_scores( utility::vector1<core::Real> &ORDER_RES_Score, const utility::vector1<core::Real> &average_res_scores, core::Size num_res ) {
	core::Size ORDER = 5;
	core::Real Divide;
	//Loop through each residue
	for ( core::Size i=1; i<=num_res; i++ ) {
		core::Real TOTAL_Sum = 0.0;
		if ( i<=ORDER ) {
			//Case 1: within first 5 residues
			Divide = 0;
			for ( core::Size j=1; j<=(i+ORDER); j++ ) {
				TOTAL_Sum = TOTAL_Sum + average_res_scores[j];
				Divide++;
			}
		} else if ( (i+ORDER) > num_res ) {
			//Case 2: within 5 residues of end
			Divide = 0;
			for ( core::Size j=(i-ORDER); j<=num_res; j++ ) {
				TOTAL_Sum = TOTAL_Sum + average_res_scores[j];
				Divide++;
			}
		} else {
			//Case 3: not within 5 residues of either terminus
			Divide = 0;
			for ( core::Size j=(i-ORDER); j<=(i+ORDER); j++ ) {
				TOTAL_Sum = TOTAL_Sum + average_res_scores[j];
				Divide++;
			}
		}
		core::Real AVG = TOTAL_Sum / Divide;
		ORDER_RES_Score.push_back(AVG);
	}
}

core::Real order_score_score( utility::vector1<core::Real> &ORDER_RES_Score, core::Size strong_resi, core::Real os_pf, bool cate, bool num ){
	core::Real oscore_score (0.0);
	if ( cate == true ) {
		core::Real osupper (-1.98194512109);
		core::Real oslower (-2.79566984295);
		core::Real osmax (1.4848);
		core::Real osmin (-4.92147);
		if ( ORDER_RES_Score[strong_resi] >= osmax ) {
			oscore_score += 5.0;
		}
		if ( ORDER_RES_Score[strong_resi] < osmax and ORDER_RES_Score[strong_resi] > osupper ) {
			oscore_score += (2*pow((-((ORDER_RES_Score[strong_resi] - osupper) - (osmax - osupper))/(osmax - osupper)),3) - 3*pow((-((ORDER_RES_Score[strong_resi] - osupper) - (osmax - osupper))/(osmax - osupper)),2) + 1)  * 5.0;
		}
		if ( ORDER_RES_Score[strong_resi] > osmin and ORDER_RES_Score[strong_resi] < oslower ) {
			oscore_score += (2*pow((-((oslower - ORDER_RES_Score[strong_resi]) - (oslower - osmin))/(oslower - osmin)),3) - 3*pow((-((oslower - ORDER_RES_Score[strong_resi]) - (oslower - osmin))/(oslower - osmin)),2) + 1) * -5.0;
		}
		if ( ORDER_RES_Score[strong_resi] <= osmin ) {
			oscore_score += -5.0;
		}
	}
	if ( num == true ) {
		core::Real predicted_os = os_pf * -0.078144 - 2.382219;
		core::Real os_diff = fabs(predicted_os - ORDER_RES_Score[strong_resi]);
		if ( os_diff <= 0.2 ) {
			oscore_score += -4.0;
		}
		if ( os_diff > 0.2 and os_diff <= 0.6 ) {
			oscore_score += ((os_diff/2.0) - (0.7)) *4.0;
		}
	}
	return oscore_score;
}

core::Real nc_score( core::pose::Pose pose, core::Size nc_res, core::Real nc_pf, bool cate, bool num ){
	core::Real neighbor_count (0.0);
	core::Real neighbor_score (0.0);
	std::string target_atom_vector_start ("N");
	std::string target_atom ("H");
	if ( nc_res == 1 ) {
		target_atom = "2H";
	}
	if ( pose.residue(nc_res).type().name1() == 'P' ) {
		target_atom = "CB";
	}
	numeric::xyzVector<core::Real> target_vector = pose.residue(nc_res).xyz(target_atom) - pose.residue(nc_res).xyz(target_atom_vector_start);
	core::Real const distance_internal = pose.residue(nc_res).xyz(target_atom).distance(pose.residue(nc_res).xyz(target_atom_vector_start));
	for ( core::Size nc_res_target = 1; nc_res_target <= nres_protein( pose ); ++nc_res_target ) {
		if ( nc_res != nc_res_target ) {
			std::string neighbor_O ("O");
			neighbor_count += count_neighbors(neighbor_O, pose, nc_res, target_atom_vector_start, distance_internal, target_vector, nc_res_target);
			if ( pose.residue(nc_res_target).type().name1() == 'N' ) {
				std::string neighbor_OD1 ("OD1") ;
				neighbor_count += count_neighbors(neighbor_OD1, pose, nc_res, target_atom_vector_start, distance_internal, target_vector, nc_res_target);
			}
			if ( pose.residue(nc_res_target).type().name1() == 'Y' ) {
				std::string neighbor_OH ("OH") ;
				neighbor_count += count_neighbors(neighbor_OH, pose, nc_res, target_atom_vector_start, distance_internal, target_vector, nc_res_target);
			}
			if ( pose.residue(nc_res_target).type().name1() == 'D' ) {
				std::string neighbor_OD1 ("OD1") ;
				neighbor_count += count_neighbors(neighbor_OD1, pose, nc_res, target_atom_vector_start, distance_internal, target_vector, nc_res_target);
				std::string neighbor_OD2 ("OD2") ;
				neighbor_count += count_neighbors(neighbor_OD2, pose, nc_res, target_atom_vector_start, distance_internal, target_vector, nc_res_target);
			}
			if ( pose.residue(nc_res_target).type().name1() == 'S' ) {
				std::string neighbor_OG ("OG") ;
				neighbor_count += count_neighbors(neighbor_OG, pose, nc_res, target_atom_vector_start, distance_internal, target_vector, nc_res_target);
			}
			if ( pose.residue(nc_res_target).type().name1() == 'T' ) {
				std::string neighbor_OG1 ("OG1") ;
				neighbor_count += count_neighbors(neighbor_OG1, pose, nc_res, target_atom_vector_start, distance_internal, target_vector, nc_res_target);
			}
			if ( pose.residue(nc_res_target).type().name1() == 'E' ) {
				std::string neighbor_OE1 ("OE1") ;
				neighbor_count += count_neighbors(neighbor_OE1, pose, nc_res, target_atom_vector_start, distance_internal, target_vector, nc_res_target);
				std::string neighbor_OE2 ("OE2") ;
				neighbor_count += count_neighbors(neighbor_OE2, pose, nc_res, target_atom_vector_start, distance_internal, target_vector, nc_res_target);
			}
			if ( pose.residue(nc_res_target).type().name1() == 'Q' ) {
				std::string neighbor_OE1 ("OE1") ;
				neighbor_count += count_neighbors(neighbor_OE1, pose, nc_res, target_atom_vector_start, distance_internal, target_vector, nc_res_target);
			}
		}
	}
	if ( cate == true ) {
		core::Real ncupper (13.531024077);
		core::Real nclower (9.28908266135);
		core::Real ncmax (17.2781);
		core::Real ncmin (2.23526);
		if ( neighbor_count >= ncmax ) {
			neighbor_score += -1.0;
		}
		if ( neighbor_count < ncmax and neighbor_count > ncupper ) {
			neighbor_score += -1*(2 * pow((-((neighbor_count - ncupper) - (ncmax - ncupper))/(ncmax - ncupper)),3) - 3 * pow((-((neighbor_count - ncupper) - (ncmax - ncupper))/(ncmax - ncupper)),2) + 1);
		}
		if ( neighbor_count < nclower and neighbor_count > ncmin ) {
			neighbor_score += 2 * pow((-((nclower - neighbor_count) - (nclower-ncmin))/(nclower - ncmin)),3) - 3 * pow((-((nclower - neighbor_count)-(nclower - ncmin))/(nclower - ncmin)),2) + 1;
		}
		if ( neighbor_count <= ncmin ) {
			neighbor_score += 1;
		}
	}
	if ( num == true ) {
		core::Real predicted_nc = nc_pf * 0.505848 + 8.877848;
		core::Real nc_diff = fabs(predicted_nc - neighbor_count);
		if ( nc_diff <= 1.0 ) {
			neighbor_score += -2.0;
		}
		if ( nc_diff > 1.0 and nc_diff <= 2.0 ) {
			neighbor_score += ((nc_diff/2.0) - (1.25)) *2.0;
		}
	}
	return neighbor_score;
}

core::Real rel_sasa_score( core::pose::Pose pose, std::string sapseq, std::map<char,core::Real> left, std::map<char,core::Real> right, utility::vector1< Real > rel_sasa, core::Size sasa_res, core::Real sasa_pf, bool cate, bool num ){
	core::Real sasa_score (0.0);
	if ( cate == true ) {
		core::Real sasamax (1.00057);
		core::Real sasamin (0.0);
		core::Real sasaupper (0.30456718439);
		core::Real sasalower (0.16960955369);
		if ( sasa_res == 1 ) {
			sasaupper -= (0.2 * (left[sapseq[sasa_res-1]] + right['Z']));
			sasalower -= (0.2 * (left[sapseq[sasa_res-1]] + right['Z']));
		}
		if ( sasa_res == nres_protein( pose ) ) {
			sasaupper -= (0.2 * (left['Z'] + right[sapseq[sasa_res - 2]]));
			sasalower -= (0.2 * (left['Z'] + right[sapseq[sasa_res - 2]]));
		}
		if ( sasa_res > 1 and sasa_res < nres_protein( pose ) ) {
			sasaupper -= (0.2 * (left[sapseq[sasa_res-1]] + right[sapseq[sasa_res - 2]]));
			sasalower -= (0.2 * (left[sapseq[sasa_res-1]] + right[sapseq[sasa_res - 2]]));
		}
		if ( rel_sasa[sasa_res] >= sasamax ) {
			sasa_score += 5.0;
		}
		if ( rel_sasa[sasa_res] <= sasamax and rel_sasa[sasa_res] > sasaupper ) {
			sasa_score += (2*pow((-((rel_sasa[sasa_res] - sasaupper) - (sasamax - sasaupper))/(sasamax - sasaupper)),3) - 3*pow((-((rel_sasa[sasa_res] - sasaupper) - (sasamax - sasaupper))/(sasamax - sasaupper)),2) + 1) * 5.0;
		}
		if ( rel_sasa[sasa_res] < sasalower and rel_sasa[sasa_res] > sasamin ) {
			sasa_score += (2*pow((-((sasalower - rel_sasa[sasa_res])-(sasalower - sasamin))/(sasalower - sasamin)),3) - 3*pow((-((sasalower - rel_sasa[sasa_res])-(sasalower - sasamin))/(sasalower - sasamin)),2) + 1) * -5.0;
		}
		if ( rel_sasa[sasa_res] <= sasamin ) {
			sasa_score += -5.0;
		}
	}
	if ( num == true ) {
		core::Real predicted_sasa = sasa_pf * -0.036882 + 0.412818;
		core::Real sasa_diff = fabs(predicted_sasa - rel_sasa[sasa_res]);
		if ( sasa_diff <= 0.1 ) {
			sasa_score += -3.0;
		}
		if ( sasa_diff > 0.1 and sasa_diff <= 0.3 ) {
			sasa_score += ((sasa_diff/2.0) - (0.6)) *3.0;
		}
	}
	return sasa_score;
}

core::Real hb_score(std::string sapseq, std::list<core::Size> strong_residues_cat, std::vector<core::Size> strong_residues_pf, std::vector<core::Real> pf_list, core::pose::Pose pose, std::map<char,core::Real> left, std::map<char,core::Real> right, core::scoring::ScoreFunctionOP scorefunc, bool cate, bool num ){
	core::Real hbond_score (0.0);
	core::scoring::hbonds::HBondSet hbs1;
	(*scorefunc)(pose);
	core::scoring::hbonds::fill_hbond_set(pose, false, hbs1);
	if ( cate == true ) {
		for ( Size hb_val=1; hb_val<= Size(hbs1.nhbonds()); ++hb_val ) {
			core::scoring::hbonds::HBond const & hb( hbs1.hbond(hb_val) );
			core::conformation::Residue const & donor = pose.residue( hb.don_res() );
			if ( std::find (strong_residues_cat.begin(), strong_residues_cat.end(), pose.pdb_info()->number(donor.seqpos())) != strong_residues_cat.end() and hb.don_hatm_is_backbone() == 1 ) {
				core::Real hbenergy (0.0);
				core::Real hbupper (-1.19258593628);
				core::Real hblower (-1.53016315462);
				core::Real hbmax (-0.048);
				core::Real hbmin (-4.039);
				hbenergy += hb.energy();
				if ( pose.pdb_info()->number(donor.seqpos()) == 1 ) {
					hbupper -= (3.0 * (left[sapseq[pose.pdb_info()->number(donor.seqpos()) - 1]] + right['Z']));
					hblower -= (3.0 * (left[sapseq[pose.pdb_info()->number(donor.seqpos()) - 1]] + right['Z']));
				}
				if ( (unsigned)pose.pdb_info()->number(donor.seqpos()) == nres_protein( pose ) ) {
					hbupper -= (3.0 * (left['Z'] + right[sapseq[pose.pdb_info()->number(donor.seqpos()) - 2]]));
					hblower -= (3.0 * (left['Z'] + right[sapseq[pose.pdb_info()->number(donor.seqpos()) - 2]]));
				}
				if ( pose.pdb_info()->number(donor.seqpos()) > 1 and (unsigned)pose.pdb_info()->number(donor.seqpos()) < nres_protein( pose ) ) {
					hbupper -= (3.0 * (left[sapseq[pose.pdb_info()->number(donor.seqpos())- 1]] + right[sapseq[pose.pdb_info()->number(donor.seqpos()) - 2]]));
					hblower -= (3.0 * (left[sapseq[pose.pdb_info()->number(donor.seqpos()) - 1]] + right[sapseq[pose.pdb_info()->number(donor.seqpos()) - 2]]));
				}
				if ( hbenergy >= hbmax ) {
					hbond_score += 30.0;
				}
				if ( hbenergy <= hbmax and hbenergy > hbupper ) {
					hbond_score += (2*pow((-((hbenergy - hbupper) - (hbmax - hbupper))/(hbmax - hbupper)),3) - 3*pow((-((hbenergy - hbupper) - (hbmax - hbupper))/(hbmax - hbupper)),2) + 1) * 30.0;
				}
				if ( hbenergy < hblower and hbenergy > hbmin ) {
					hbond_score += (2*pow((-((hblower - hbenergy)-(hblower - hbmin))/(hblower - hbmin)),3) - 3*pow((-((hblower - hbenergy)-(hblower - hbmin))/(hblower - hbmin)),2) + 1) * -30.0;
				}
				if ( hbenergy <= hbmin ) {
					hbond_score += -30.0;
				}
			}
		}
	}
	if ( num == true ) {
		for ( Size hb_val=1; hb_val<= Size(hbs1.nhbonds()); ++hb_val ) {
			core::scoring::hbonds::HBond const & hb( hbs1.hbond(hb_val) );
			core::conformation::Residue const & donor = pose.residue( hb.don_res() );
			core::Size hb_count (0);
			core::Real hbenergy (0.0);
			hbenergy += hb.energy();
			if ( std::find(strong_residues_pf.begin(), strong_residues_pf.end(), pose.pdb_info()->number(donor.seqpos())) != strong_residues_pf.end() and hb.don_hatm_is_backbone() == 1 ) {
				while ( (unsigned)pose.pdb_info()->number(donor.seqpos()) != strong_residues_pf[hb_count] ) {
					hb_count += 1;
				}
				if ( (unsigned)pose.pdb_info()->number(donor.seqpos()) == strong_residues_pf[hb_count] ) {
					core::Real predicted_hb = pf_list[hb_count] * -0.11427 - 0.806135;
					core::Real hb_diff = fabs(predicted_hb - hbenergy);
					if ( hb_diff <= 0.2 ) {
						hbond_score += -3.0;
					}
					if ( hb_diff > 0.2 and hb_diff <= 0.8 ) {
						hbond_score += ((hb_diff/2.0) - (0.75)) *3.0;
					}
				}
			}
		}
	}
	return hbond_score;
}

int
main( int argc, char * argv [] )
{
	try {

		using namespace basic::options;
		using namespace utility;
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring::sasa;
		NEW_OPT( StrongRes, "Used if strong residue categorical list is provided", false);
		NEW_OPT( ResPF, "Used if protection factors are included for each residue", false);

		// initialize options, RNG, and factory-registrators
		devel::init(argc, argv);
		register_options();

		// check if PDB given
		if ( ! option[OptionKeys::in::file::s].user() ) {
			throw CREATE_EXCEPTION(utility::excn::Exception, "Please provide PDB file!");
		}

		bool cat = option[ StrongRes ];
		bool numer = option[ ResPF ];

		if ( cat == false and numer == false ) {
			throw CREATE_EXCEPTION(utility::excn::Exception, "Please choose -numerical or -categorical flags");
		}

		if ( cat == true and numer == true ) {
			throw CREATE_EXCEPTION(utility::excn::Exception, "Please choose -numerical or -categorical flags, not both");
		}

		// read in PDB to Pose
		Pose pose;
		core::import_pose::pose_from_file( pose, option[OptionKeys::in::file::s].value_string() , core::import_pose::PDB_file);
		core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
		//Score of the entire pose
		core::Real rosettascore = (*sfxn)(pose);
		//Initializing Pose HDX Score
		core::Real pose_score (0.0);
		pose_score += rosettascore;
		//Initializing residue list
		std::list<core::Size> input_strong_residues_cat = {};
		std::vector<core::Size> input_strong_residues_pf = {};
		std::vector<core::Real> input_strong_pf_pf = {};
		core::Size sasa_count (0);
		core::Size nc_count (0);
		core::Size os_count (0);

		utility::io::izstream input(option[OptionKeys::in::file::HDX]);
		std::string line;
		while ( getline(input,line) ) {
			std::istringstream ss(line);
			core::Size resi;
			core::Real PF;
			ss >> resi >> PF;
			if ( cat == true ) {
				input_strong_residues_cat.push_back( resi );
			}
			if ( numer == true ) {
				input_strong_residues_pf.push_back( resi );
				input_strong_pf_pf.push_back(PF);
			}
		}
		if ( cat==true ) {
			//Initializing SAP values and sequence
			std::map<char,core::Real> leftres = {{'V',-0.74}, {'I',-0.91}, {'L',-0.57}, {'E',-0.90}, {'Q',-0.47},{'D',0.90}, {'N',-0.58}, {'H',-0.80}, {'W',-0.40},
				{'F',-0.52}, {'Y',-0.41}, {'R',-0.59}, {'K',-0.56}, {'S',-0.44}, {'T',-0.79}, {'M',-0.64}, {'A',0.0},{'G',-0.22}, {'P',0.0}, {'C',-0.54} ,{'Z',0.96}};
			std::map<char,core::Real> rightres = {{'V',-0.30}, {'I',-0.59}, {'L',-0.13}, {'E',0.31}, {'Q',-0.27},{'D',0.58}, {'N',-0.13}, {'H',-0.51}, {'W',-0.44},
				{'F',-0.43}, {'Y',-0.37}, {'R',-0.32}, {'K',-0.29}, {'S',-0.39}, {'T',-0.47}, {'M',-0.28}, {'A',0.0},{'G',0.22}, {'P',-0.19}, {'C',-0.46}, {'Z',-1.32}};
			std::string sequence = pose.sequence();

			//Calculate OS
			utility::vector1<core::Real> average_res_scores;
			//Calculate average residue scores
			calculate_average_residue_scores( pose, average_res_scores, nres_protein( pose ));
			//Calculate the order scores (i.e. window average score of each residue)
			utility::vector1<core::Real> ORDER_RES_Score;
			calculate_order_scores( ORDER_RES_Score, average_res_scores, nres_protein( pose ) );
			for ( core::Size os_num = 1; os_num <= nres_protein( pose ); ++os_num ) {
				if ( std::find (input_strong_residues_cat.begin(), input_strong_residues_cat.end(), os_num) != input_strong_residues_cat.end() ) {
					pose_score += order_score_score( ORDER_RES_Score, os_num, os_count, cat, numer);
				}
			}
			//Calculating HB
			pose_score += hb_score( sequence, input_strong_residues_cat, input_strong_residues_pf, input_strong_pf_pf, pose, leftres, rightres, sfxn, cat, numer );

			//Calculating rel_sasa for the pose
			utility::vector1< Real > rel_sasa = core::scoring::sasa::rel_per_res_sc_sasa( pose );
			for ( core::Size sasa_val = 1; sasa_val <= nres_protein( pose ); ++sasa_val ) {
				if ( std::find (input_strong_residues_cat.begin(), input_strong_residues_cat.end(), sasa_val) != input_strong_residues_cat.end() ) {
					pose_score += rel_sasa_score( pose, sequence, leftres, rightres, rel_sasa, sasa_val, sasa_count, cat, numer);
				}
			}

			//Neighbor Count Portion (Finished calc, need scoring part)
			for ( core::Size nc_res =1; nc_res <= nres_protein( pose ); ++nc_res ) {
				if ( std::find (input_strong_residues_cat.begin(), input_strong_residues_cat.end(), nc_res) != input_strong_residues_cat.end() ) {
					pose_score += nc_score( pose, nc_res, nc_count, cat, numer);
				}
			}
		}

		if ( numer==true ) {
			//Initializing SAP values and sequence
			std::map<char,core::Real> leftres = {{'V',-0.74}, {'I',-0.91}, {'L',-0.57}, {'E',-0.90}, {'Q',-0.47},{'D',0.90}, {'N',-0.58}, {'H',-0.80}, {'W',-0.40},
				{'F',-0.52}, {'Y',-0.41}, {'R',-0.59}, {'K',-0.56}, {'S',-0.44}, {'T',-0.79}, {'M',-0.64}, {'A',0.0},{'G',-0.22}, {'P',0.0}, {'C',-0.54} ,{'Z',0.96}};
			std::map<char,core::Real> rightres = {{'V',-0.30}, {'I',-0.59}, {'L',-0.13}, {'E',0.31}, {'Q',-0.27},{'D',0.58}, {'N',-0.13}, {'H',-0.51}, {'W',-0.44},
				{'F',-0.43}, {'Y',-0.37}, {'R',-0.32}, {'K',-0.29}, {'S',-0.39}, {'T',-0.47}, {'M',-0.28}, {'A',0.0},{'G',0.22}, {'P',-0.19}, {'C',-0.46}, {'Z',-1.32}};
			std::string sequence = pose.sequence();

			//Calculate OS
			utility::vector1<core::Real> average_res_scores;
			//Calculate average residue scores
			calculate_average_residue_scores( pose, average_res_scores, nres_protein( pose ));
			//Calculate the order scores (i.e. window average score of each residue)
			utility::vector1<core::Real> ORDER_RES_Score;
			calculate_order_scores( ORDER_RES_Score, average_res_scores, nres_protein( pose ) );
			for ( core::Size os_num = 1; os_num <= nres_protein( pose ); ++os_num ) {
				if ( std::find (input_strong_residues_pf.begin(), input_strong_residues_pf.end(), os_num) != input_strong_residues_pf.end() ) {
					pose_score += order_score_score( ORDER_RES_Score, os_num, input_strong_pf_pf[os_count], cat, numer);
					os_count += 1;
				}
			}
			//Calculating HB

			pose_score += hb_score( sequence, input_strong_residues_cat, input_strong_residues_pf, input_strong_pf_pf, pose, leftres, rightres, sfxn, cat, numer );

			//Calculating rel_sasa for the pose

			utility::vector1< Real > rel_sasa = core::scoring::sasa::rel_per_res_sc_sasa( pose );
			for ( core::Size sasa_val = 1; sasa_val <= nres_protein( pose ); ++sasa_val ) {
				if ( std::find (input_strong_residues_pf.begin(), input_strong_residues_pf.end(), sasa_val) != input_strong_residues_pf.end() ) {
					pose_score += rel_sasa_score( pose, sequence, leftres, rightres, rel_sasa, sasa_val, input_strong_pf_pf[sasa_count], cat, numer);
					sasa_count += 1;
				}
			}

			//Neighbor Count Portion (Finished calc, need scoring part)
			for ( core::Size nc_res =1; nc_res <= nres_protein( pose ); ++nc_res ) {
				if ( std::find (input_strong_residues_pf.begin(), input_strong_residues_pf.end(), nc_res) != input_strong_residues_pf.end() ) {
					pose_score += nc_score( pose, nc_res, input_strong_pf_pf[nc_count], cat, numer);
					nc_count +=1;
				}
			}
		}


		std::ostringstream output;
		output << pose_score << std::endl;
		output_results(output.str());
	}

catch (utility::excn::Exception const & e ) {
	e.display();
	return -1;
}

	return 0;
}
