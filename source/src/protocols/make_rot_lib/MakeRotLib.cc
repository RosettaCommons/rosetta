// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// unit headers
#include <protocols/make_rot_lib/RotData.hh>
#include <protocols/make_rot_lib/MakeRotLib.hh>

// core headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>

// utility headers
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/exit.hh>
#include <utility/vector1.functions.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// numeric headers
#include <numeric/angle.functions.hh>

// protocol headers
#include <protocols/simple_moves/MinMover.hh>

// external headers
#include <boost/lexical_cast.hpp>

// C++ headers
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>


static THREAD_LOCAL basic::Tracer TR("protocols.make_rot_lib.MakeRotLib");

namespace protocols {
namespace make_rot_lib {

using namespace core;
using namespace utility;


// helper function returns the difference between to angles
Real
angle_diff( Real a1, Real a2 )
{
	Real t1( fabs(a1 - a2) );
	Real t2( 360 - t1 );
	return( t1 <= t2 ? t1 : t2 );
}

// hack hack hack
void
asp_corrections( RotVec & rotamers )
{
	for ( Size i = 1; i <= rotamers.size(); ++i ) {
		Real temp_chi( rotamers[i].get_min_chi(2) );
		if ( temp_chi >= 90 && temp_chi <= 180 ) {
			rotamers[i].set_min_chi( temp_chi+180, 2 );
		} else if ( temp_chi > 180 && temp_chi <= 270 ) {
			rotamers[i].set_min_chi( temp_chi-180, 2 );
		}
	}
}

void
glu_corrections( RotVec & rotamers )
{
	for ( Size i = 1; i <= rotamers.size(); ++i ) {
		Real temp_chi( rotamers[i].get_min_chi(3) );
		if ( temp_chi >= 90 && temp_chi <= 180 ) {
			rotamers[i].set_min_chi( temp_chi+180, 3 );
		} else if ( temp_chi > 180 && temp_chi <= 270 ) {
			rotamers[i].set_min_chi( temp_chi-180, 3 );
		}
	}
}

void
phe_tyr_corrections( RotVec & rotamers )
{
	for ( Size i = 1; i <= rotamers.size(); ++i ) {
		Real temp_chi( rotamers[i].get_min_chi(2) );
		if ( temp_chi >= 135 && temp_chi <= 180 ) {
			rotamers[i].set_min_chi( temp_chi+180, 2 );
		} else if ( temp_chi > 180 && temp_chi <= 315 ) {
			rotamers[i].set_min_chi( temp_chi-180, 2 );
		}
	}
}


void
min_rotamers( RotVec & rotamers,  core::scoring::ScoreFunctionOP scrfxn, std::string aa_name )
{
	using namespace core;
	using namespace pose;
	using namespace scoring;
	using namespace chemical;
	using namespace conformation;

	TR << "In min_rotamers()..." << std::flush << std::endl;

	// iterate over rotamers array
	Size const nrot( rotamers.size() );
	for ( Size i = 1; i <= nrot; ++i ) {

		TR << "Working on rotamer " << i << " of " << nrot << std::flush << std::endl;

		// get number of chi
		Size const nchi( rotamers[i].get_num_chi() );

		// make a pose, residue, and append residue to pose
		Pose pose;
		ResidueTypeSetCOP RTS( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		ResidueType const & RT( RTS->name_map( aa_name ) );
		Residue R( RT, true );
		pose.append_residue_by_jump( R, 1 );

		bool is_peptoid( RT.is_peptoid() );

		// create filenames
		/*std::stringstream start_name;
		if( is_peptoid ) {
		start_name << "tripeptide_" << rotamers[i].get_omg() << "_" << rotamers[i].get_phi() << "_" << rotamers[i].get_psi() << "_";
		} else {
		start_name << "tripeptide_" << rotamers[i].get_phi() << "_" << rotamers[i].get_psi() << "_";
		}
		for(Size j = 1; j <= nchi; ++j ) {
		start_name << rotamers[i].get_inp_chi( j ) << "_";
		}
		start_name << "start.pdb";
		//std::string start_filename( start_name.str() );

		std::stringstream end_name;
		if( is_peptoid ) {
		end_name << "tripeptide_" << rotamers[i].get_omg() << "_" << rotamers[i].get_phi() << "_" << rotamers[i].get_psi() << "_";
		} else {
		end_name << "tripeptide_" << rotamers[i].get_phi() << "_" << rotamers[i].get_psi() << "_";
		}
		for(Size j = 1; j <= nchi; ++j ) {
		end_name << rotamers[i].get_inp_chi( j ) << "_";
		}
		end_name << "end.pdb";
		*/
		//std::string end_filename( end_name.str() );

		// score the pose
		//Real debug_ener( (*scrfxn)( pose ) );
		//TR << "DEBUG ENER: " << debug_ener << std::endl;

		// set phi, psi, omg, eps and chi(s)
		id::TorsionID bb1( 1, id::BB, 1 );
		id::TorsionID bb2( 1, id::BB, 2 );
		id::TorsionID bb3( 1, id::BB, 3 );
		id::TorsionID bb4( 1, id::BB, 4 );

		pose.set_torsion( bb1, rotamers[i].get_omg() );
		pose.set_torsion( bb2, rotamers[i].get_phi() );
		pose.set_torsion( bb3, rotamers[i].get_psi() );
		pose.set_torsion( bb4, rotamers[i].get_eps() );

		for ( Size j = 1; j <= nchi; ++j ) {
			pose.set_chi( j, 1, rotamers[i].get_inp_chi( j ) );
		}

		// output starting pose

		// score the pose
		Real orig_ener( (*scrfxn)( pose ) );

		TR << "ORIG ENER: " << orig_ener << "\t" << pose.energies().total_energy() << std::endl;

		// minimize the pose
		kinematics::MoveMapOP mvmp( new kinematics::MoveMap );
		mvmp->set_chi( 1, true );
		mvmp->set( bb4, true );
		// keep omg fixed if we are a peptoid
		if ( !is_peptoid ) {
			mvmp->set( bb1, true );
		}

		protocols::simple_moves::MinMover mnmvr( mvmp, scrfxn, "linmin", 0.0001, true );


		Real current_ener( orig_ener), previous_ener( orig_ener );
		Size iter(0);
		do {
			previous_ener = current_ener;
			mnmvr.apply( pose );
			current_ener = pose.energies().total_energy();
			iter++;
			if ( iter > 100 ) break;
		} while ( current_ener + 0.001 < previous_ener );

		// score the pose
		Real min_ener( pose.energies().total_energy() );

		// output ending pose

		TR << "MIN ENER: " << min_ener << " found in " << iter << " steps" << std::endl;

		EnergyMap em( pose.energies().residue_total_energies(1) );
		rotamers[i].set_twist( em[ mm_twist ] );
		rotamers[i].set_inter_rep( em[ mm_lj_inter_rep ] );
		rotamers[i].set_inter_atr( em[ mm_lj_inter_atr ] );
		rotamers[i].set_intra_rep( em[ mm_lj_intra_rep ] );
		rotamers[i].set_intra_atr( em[ mm_lj_intra_atr ] );
		rotamers[i].set_solvation( em[ fa_sol ] );

		// record min chi, energy
		for ( Size j = 1; j <= nchi; ++j ) {
			rotamers[i].set_min_chi( pose.chi( j, 1 ), j );
		}
		rotamers[i].set_energy( min_ener );
		rotamers[i].set_min_omg( numeric::nonnegative_principal_angle_degrees( pose.torsion( bb1 ) ) );
		rotamers[i].set_min_eps( numeric::nonnegative_principal_angle_degrees( pose.torsion( bb4 ) ) );

	}

}

void
init_rotamers_centroids
( RotVec & rotamers,
	RotVec & centroids,
	Size & ncluster,
	std::string options_filename,
	std::string & aa_name,
	bool is_peptoid,
	Real omg_start_val,
	Real eps_start_val
)
{
	//init stuff to hold lines
	utility::vector1< std::string > lines;
	std::string line;

	// check to see if options file exists
	if ( ! utility::file::file_exists( options_filename.c_str() ) ) {
		utility_exit_with_message("Cannot find options file"+options_filename);
	}

	// open file
	std::ifstream options( options_filename.c_str() );

	// read in each line
	while ( getline( options, line ) ) {
		//std::istringstream l( line );
		if ( line.size() < 1 || line[0] == '#' ) continue;
		lines.push_back( line );
	}

	// close options file
	options.close();

	// get number of lines
	Size const nlines( lines.size() );

	// go through file once to get number of chi angles, number of clusters,
	// and aa name before parsing the rest of the file
	Size nchi( 0 );
	std::string base_aa_name;
	for ( Size i=1; i<= nlines; ++i ) {
		std::string const & line( lines[i] );
		std::istringstream l( line );
		std::string tag;

		l >> tag;

		// get number of chi
		if ( tag == "NUM_CHI" ) {
			l >> nchi;
		} else if ( tag == "CENTROID" ) {
			// get number of clusters
			++ncluster;
		} else if ( tag == "AA_NAME" ) {
			// get amino acid name
			l >> base_aa_name;
		}

	}

	// create full_aa_name
	std::string patch;
	if ( is_peptoid ) { // use this patch if we want peptoids
		patch = ":AcetylatedNtermDimethylatedCtermPeptoidFull";
		TR << "Making a PEPTOID rotamer library..." << std::endl;
	} else { // use this patch if protein
		patch = ":MethylatedCtermProteinFull:AcetylatedNtermProteinFull";
		TR << "Making a PROTEIN rotamer library..." << std::endl;
	}

	aa_name = base_aa_name + patch;


	// parse the rest of the file
	int phi_lower(0), phi_upper(0), phi_increment(0);
	int psi_lower(0), psi_upper(0), psi_increment(0);
	int omg_lower(0), omg_upper(0), omg_increment(0);
	utility::vector1<int> chi_lower, chi_upper, chi_increment;
	chi_lower.resize(nchi, 0); chi_upper.resize(nchi, 0); chi_increment.resize(nchi, 0);

	for ( Size i=1; i<= nlines; ++i ) {
		std::string const & line( lines[i] );
		std::istringstream l( line );
		std::string tag;

		l >> tag;

		// get phi range
		if ( tag == "PHI_RANGE" ) {
			l >> phi_lower >> phi_upper >> phi_increment;
		} else if ( tag == "PSI_RANGE" ) {
			// get psi range
			l >> psi_lower >> psi_upper >> psi_increment;
		} else if ( tag == "OMG_RANGE" ) {
			// get omg range
			l >> omg_lower >> omg_upper >> omg_increment;
		} else if ( tag == "CHI_RANGE" ) {
			// get chi range(s)
			Size chi_num(0); int temp_lower, temp_upper, temp_increment;
			l >> chi_num >> temp_lower >> temp_upper >> temp_increment;
			assert( chi_num <= nchi );
			chi_lower[chi_num] = temp_lower;
			chi_upper[chi_num] = temp_upper;
			chi_increment[chi_num] = temp_increment;
		} else if ( tag == "CENTROID" ) {
			// get centoids
			RotData temp( nchi, ncluster );
			for ( Size i = 1; i<= nchi; ++i ) {
				Real chi_val(0); int lib_chi_val(0);
				l >> chi_val >> lib_chi_val;
				temp.set_inp_chi( chi_val, i );
				temp.set_min_chi( chi_val, i );
				temp.set_lib_chi_val( lib_chi_val, i );
			}
			centroids.push_back( temp );
		}
	}

	// create intital based on rotamer range data
	utility::vector1<Size> total_chi_num_bins;
	total_chi_num_bins.resize(nchi, 0);
	Size nrotamers(1);
	for ( Size i = 1; i <= nchi; ++i ) {
		if ( (chi_upper[i] - chi_lower[i])%chi_increment[i] != 0 ) {
			utility_exit_with_message("Number of chi bins not integer for chi: " +
				boost::lexical_cast<std::string>(i));
		}
		total_chi_num_bins[i] = ((chi_upper[i] - chi_lower[i])/chi_increment[i])+1;
		nrotamers *= total_chi_num_bins[i];
	}

	// init rotamer array
	RotData temp( nchi, ncluster );
	rotamers.resize(nrotamers, temp);

	for ( Size i = 1; i <= nchi; ++i ) {

		// calc all chi(i) valuse based on upper, lower, inc
		utility::vector1<Real> chi_values;
		for ( int k = chi_lower[i]; k <= chi_upper[i]; k+=chi_increment[i] ) {
			chi_values.push_back( k );
		}

		// iterator to keep track of where in rotamers we are
		Size rot_iter(1);

		// calc number of rotmater combos for larger chi numbers
		Size count_up(1);
		if ( i >= nchi ) {
			count_up = 1;
		} else {
			for ( Size m = i+1; m <= nchi; ++m ) {
				count_up *= total_chi_num_bins[m];
			}
		}
		// calc number of rotmater combos for smaller chi numbers
		Size count_down(1);
		if ( i <= 1 ) {
			count_down = 1;
		} else {
			for ( Size m = 1; m < i; ++m ) {
				count_down *= total_chi_num_bins[m];
			}
		}

		// assign chi values
		for ( Size n = 1; n <= count_down; ++n ) {
			for ( Size j = 1; j <= chi_values.size(); ++j ) {
				for ( Size l = 1; l <= count_up; ++l ) {
					rotamers[rot_iter].set_inp_chi(chi_values[j], i);
					++rot_iter;
				}
			}
		}
	}

	// set phi and psi
	// ultimatly this will be done higher up in the code to initialize EVERYTHING
	// but i don't think that we will run that way at first
	for ( Size i = 1; i <= nrotamers; ++i ) {
		rotamers[i].set_phi( phi_lower );
		rotamers[i].set_psi( psi_lower );
		if ( is_peptoid ) {
			rotamers[i].set_omg( omg_lower );
		} else {
			rotamers[i].set_omg( omg_start_val );
		}
		rotamers[i].set_eps( eps_start_val );
	}
}

bool
calc_rotamer_clusters( RotVec & rotamers )
{
	// for all rotamers find closest centroid and assign cluster_num
	// compare pt to its old cluster see if it changed locations
	bool clust_change = false;
	for ( Size i=1; i<=rotamers.size(); ++i ) {
		Size oldclust = rotamers[i].get_cluster_num();
		Size clust = rotamers[i].get_min_cent_dist(); //to get closest centroid
		rotamers[i].set_cluster_num(clust);
		if ( oldclust != clust ) {
			clust_change = true;
		}
	}

	return clust_change;
}

bool
calc_centroids(RotVec & rotamers, RotVec & centroids)
{
	// for cluster i calc_dist divide by # in clust
	// find mean of all pts in cluster
	// set centroid to this mean

	typedef vector1<Real> rvec;

	Size ncluster( centroids.size() );
	Size nchi( rotamers[1].get_num_chi() );

	// init num_rots
	vector1<Size> num_rots;
	num_rots.resize(ncluster,0);

	// init total_chi_angle
	vector1< rvec > total_chi_angle;
	total_chi_angle.resize(ncluster);
	for ( Size i = 1; i <= total_chi_angle.size(); ++i ) {
		total_chi_angle[i].resize(nchi, 0);
	}

	// sum up chi angle totals
	for ( Size i = 1; i <= rotamers.size(); ++i ) {
		Size clusternum = rotamers[i].get_cluster_num();
		for ( Size j = 1; j <= nchi; ++j ) {
			total_chi_angle[clusternum][j] += rotamers[i].get_min_chi(j);
		}
		num_rots[clusternum]++;
	}

	//get average angles and assign new centroid chi's
	bool centroid_change(false);
	for ( Size i=1; i <= ncluster; ++i ) {
		for ( Size j=1; j <= nchi; ++j ) {
			Real old_angle = centroids[i].get_min_chi( j );
			Real avg_angle = total_chi_angle[i][j] / num_rots[i];
			centroids[i].set_min_chi( avg_angle, j );
			if ( old_angle != avg_angle ) {
				centroid_change = true;
			}
		}
	}

	return centroid_change;
}

Real
calc_dist( RotData & point1, RotData & point2 )
{
	Size nchi(point1.get_num_chi());
	Real sd(0);

	for ( Size i = 1; i <= nchi; ++i ) {
		sd += pow( angle_diff(point1.get_min_chi(i), point2.get_min_chi(i) ), 2 );
	}

	return sqrt( sd/nchi );
}

Real
avg_cluster_cen_dist(RotVec & rotamers, Size & ncluster)
{
	// calc_dist from centroid to rotamers in cluster

	vector1<Real> total_dist;
	vector1<Size> num_rots;

	total_dist.resize(ncluster,0);
	num_rots.resize(ncluster,0);


	for ( Size i = 1; i <= rotamers.size(); ++i ) {
		Size clusternum = rotamers[i].get_cluster_num();
		total_dist[clusternum] += rotamers[i].get_cen_dist(clusternum);
		num_rots[clusternum]++;
	}

	//print out how many rots in each cluster
	for ( Size i=1; i <= ncluster; ++i ) {
		TR << "Cluster: "<< i << " contains "<< num_rots[i] << " of "<< rotamers.size()<< "   or %=" << static_cast<Real>(num_rots[i]) / static_cast<Real>(rotamers.size())  <<std::endl;
	}

	TR <<"AVG_CLUST_CENT_DST:   " <<std::endl;
	Real avgdist(0);
	for ( Size i = 1; i <= ncluster; ++i ) {
		TR  << total_dist[i]/num_rots[i] << "\t";
		avgdist += (total_dist[i]/num_rots[i]);
	}
	TR << std::endl;
	TR << "AVG_CENT_DST:"  << "\t";
	return avgdist / total_dist.size();
}

void
calc_all_dist(RotVec & rotamers, RotVec & centroids)
{
	for ( Size i=1; i <= rotamers.size(); ++i ) {
		for ( Size j=1; j<= centroids.size(); ++j ) {
			Real dist = calc_dist (rotamers[i], centroids[j]);
			rotamers[i].set_cen_dist(dist, j);
		}
	}
}

// pull out best rots from rotamers and add them to final_rotamers
void
get_final_rots(RotVec & rotamers , RotVec & final_rotamers, Size & nclusters ){
	RotData temp(rotamers[1].get_num_chi(), nclusters);
	final_rotamers.resize(nclusters, temp);

	//sets high energies for later minimization (probably better way to do this)
	for ( Size j=1; j<= nclusters ; ++j ) {
		final_rotamers[j].set_energy(6000);
	}
	//finds lowest E rotamer in each cluster
	for ( Size i=1; i<= rotamers.size(); ++i ) {
		Size cluster_num( rotamers[i].get_cluster_num() );
		if ( rotamers[i].get_energy() <= final_rotamers[cluster_num].get_energy() ) {
			final_rotamers[cluster_num] = rotamers[i];
		}
	}
}


// calc probabilites for final rots & normilize
void
get_final_rot_probs( RotVec & final_rotamers)
{
	vector1<Real> prob_temp;
	prob_temp.resize(final_rotamers.size(),0);
	vector1<Real> normalized_ener;
	normalized_ener.resize(final_rotamers.size(),0);
	Real KbT(0.5961); //constant
	Real total_prob(0);

	// get the minimum energy
	Real min_ener(final_rotamers[1].get_energy()); // seed
	for ( Size i=1; i<=final_rotamers.size(); ++i ) {
		if ( final_rotamers[i].get_energy() < min_ener ) {
			min_ener = final_rotamers[i].get_energy();
		}
	}

	// calc the normalized energies
	for ( Size i=1; i<=final_rotamers.size(); ++i ) {
		normalized_ener[i] = final_rotamers[i].get_energy() - min_ener;
	}

	// finds probabilities from energy
	for ( Size i=1; i<=final_rotamers.size(); ++i ) {
		prob_temp[i] = exp( ( -normalized_ener[i] ) / KbT ) ;
		total_prob += prob_temp[i];
		TR <<"Cluster:  " << i << "   Probability=" << prob_temp[i] << std::endl;
	}
	TR <<"Total Probability= " << total_prob << std::endl;
	TR <<"Kb T value used: " << KbT << std::endl;

	//normalizes and sets probabilities of final_rotamers
	for ( Size i=1; i<=final_rotamers.size(); ++i ) {
		final_rotamers[i].set_probability(prob_temp[i] / total_prob);
	}
}

//find standard deviation of chi
void
calc_std_dev (RotVec & final_rotamers, core::scoring::ScoreFunctionOP scrfxn, std::string aa_name)
{
	//copy from min_rotamers function
	using namespace pose;
	using namespace scoring;
	using namespace chemical;
	using namespace conformation;

	// make a pose, residue, and append residue to pos
	Pose pose;
	ResidueTypeSetCOP RTS( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	ResidueType const & RT( RTS->name_map( aa_name ) );
	Residue R( RT, true );
	pose.append_residue_by_jump( R, 1 );

	bool is_peptoid( RT.is_peptoid() );

	// itterate over final rotamers
	for ( Size i = 1; i<= final_rotamers.size(); ++i ) {
		Size const nchi (final_rotamers[i].get_num_chi() );

		// set phi, psi, chi
		id::TorsionID bb1( 1, id::BB, 1 );
		id::TorsionID bb2( 1, id::BB, 2 );
		id::TorsionID bb3( 1, id::BB, 3 );
		id::TorsionID bb4( 1, id::BB, 4 );
		pose.set_torsion( bb1, final_rotamers[i].get_omg() );
		pose.set_torsion( bb2, final_rotamers[i].get_phi() );
		pose.set_torsion( bb3, final_rotamers[i].get_psi() );
		pose.set_torsion( bb4, final_rotamers[i].get_eps() );

		for ( Size l = 1; l <= nchi; ++l ) {
			pose.set_chi( l, 1, final_rotamers[i].get_inp_chi( l ) );
		}

		// minimize the pose
		kinematics::MoveMapOP mvmp( new kinematics::MoveMap );
		mvmp->set_chi( 1, true );
		mvmp->set( bb4, true );
		if ( !is_peptoid ) {
			mvmp->set( bb1, true );
		}
		protocols::simple_moves::MinMover mnmvr( mvmp, scrfxn, "linmin", 0.0001, true );
		for ( Size j=1; j<= 25; ++j ) {
			//TR << " ----- " << j << " ----- " << std::flush << std::endl;
			mnmvr.apply( pose );
		}

		// dump coords to a file
		/*std::stringstream final_name;
		if( is_peptoid ) {
		final_name << pose.residue(1).type().name3() << "_" <<  final_rotamers[i].get_omg() << "_" << final_rotamers[i].get_phi() << "_" << final_rotamers[i].get_psi() << "_";
		} else {
		final_name << pose.residue(1).type().name3() << "_" <<  final_rotamers[i].get_phi() << "_" << final_rotamers[i].get_psi() << "_";
		}
		final_name << "INP_CHI_";
		for(Size j = 1; j <= nchi; ++j ) {
		final_name << final_rotamers[i].get_inp_chi( j ) << "_";
		}
		final_name << "MIN_CHI_";
		for(Size j = 1; j <= nchi; ++j ) {
		final_name << final_rotamers[i].get_min_chi( j ) << "_";
		}
		final_name << "final.pdb";
		*/
		//std::string final_filename( final_name.str() );

		//set energy barrier to stop at
		Real cut_off_ener (final_rotamers[i].get_energy() + 0.5);

		//test to see if getting right Energy for this rotamer
		Real test_ener ( (*scrfxn)(pose) );
		Real found_ener (final_rotamers[i].get_energy());
		if ( test_ener != found_ener ) {
			TR << "--------WARNING---------"<< std::endl;
			TR << "For Final_Rot: " << i << "  Scored E: "<<  test_ener  << "  but Min Rot E: "<< found_ener <<std::endl;
		}

		//set step size of chi
		Real chi_step (0.1);

		// search around minimized chi's
		//TR <<"pre-for loop" <<std::endl;
		for ( Size j = 1; j <= nchi; ++j ) {
			Size k(0);
			Real new_ener (0);
			//TR <<"pre-while loop" <<std::endl;
			while ( new_ener < cut_off_ener && (k * chi_step) <= 30 ) {
				k=k+1;
				pose.set_chi( j, 1, numeric::nonnegative_principal_angle_degrees( final_rotamers[i].get_min_chi( j ) + k * chi_step ));
				Real plus_ener( (*scrfxn)( pose ) );
				pose.set_chi( j, 1, numeric::nonnegative_principal_angle_degrees( final_rotamers[i].get_min_chi( j ) - k * chi_step ));
				Real neg_ener( (*scrfxn)( pose ) );
				// set new_ener to highest energy in either direction
				new_ener = (plus_ener >= neg_ener ? plus_ener : neg_ener);

				//TR << "While Loop run: " << k <<" for chi: "<< j
				//<<" for final_rot: " << i <<  "   Current E= " << new_ener << "  Cutoff E=" << cut_off_ener << std::endl;
			}
			final_rotamers[i].set_std_dev(k*chi_step, j);

			//reset pose
			pose.set_chi( j, 1, final_rotamers[i].get_min_chi( j ) );
		}
	}
}


void
pretty_print_rd( RotData & rot )
{
	using namespace std;

	cout << "PHI: " << setw(4) << setprecision(2) << fixed << rot.get_phi() << " "
		<< "PSI: " << setw(4) << setprecision(2) << fixed << rot.get_psi() << " "
		<< "OMG: " << setw(8) << setprecision(4) << fixed << rot.get_min_omg() << " "
		<< "ESL: " << setw(8) << setprecision(4) << fixed << rot.get_min_eps() << " "
		<< "PRB: " << setw(4) << setprecision(4) << fixed << rot.get_probability() << " "
		<< "ENR: " << setw(8) << setprecision(4) << fixed << rot.get_energy() << " "
		<< "TWS: " << setw(8) << setprecision(4) << fixed << rot.get_twist() << " "
		<< "RAR: " << setw(8) << setprecision(4) << fixed << rot.get_intra_rep() << " "
		<< "RAA: " << setw(8) << setprecision(4) << fixed << rot.get_intra_atr() << " "
		//<< "ERR: " << setw(8) << setprecision(4) << fixed << rot.get_inter_rep() << " "
		// << "ERA: " << setw(8) << setprecision(4) << fixed << rot.get_inter_atr() << " "
		// << "SOL: " << setw(8) << setprecision(4) << fixed << rot.get_solvation() << " "
		<< "NCH: " << setw(4) << fixed << rot.get_num_chi() << " "
		<< "CLN: " << setw(4) << fixed << rot.get_cluster_num() << " ";

	// print inp_chi_ vector
	cout << "ICHI:";
	for ( Size i=1; i<=rot.get_num_chi(); ++i ) {
		cout << setw(7) << setprecision(2) << fixed << rot.get_inp_chi(i) << " ";
	}
	// print min_chi_ vector
	cout << "MCHI:";
	for ( Size i=1; i<=rot.get_num_chi(); ++i ) {
		cout <<  setw(7) <<  setprecision(2) << fixed << rot.get_min_chi(i) << " ";
	}
	// print std_dev_ vector
	cout << "STDD:";
	for ( Size i=1; i<=rot.get_num_chi(); ++i ) {
		cout <<  setw(7) <<  setprecision(2) << fixed << rot.get_std_dev(i) << " ";
	}
	// print cen_dst_ vector
	cout << "CDST:";
	for ( Size i=1; i<=rot.get_num_clusters(); ++i ) {
		cout <<  setw(7) <<  setprecision(2) << fixed << rot.get_cen_dist(i) << " ";
	}

	cout << endl;
}

/*
Dunbrack file format (unfortunatly limited to 4 chi angles) so some of this code is ugly and repetitive
PHE   -60  -40  2936    2 1 0 0  0.617431   179.7    78.9     0.0     0.0      10.0    12.3     0.0     0.0
PHE   -60  -40  2936    3 1 0 0  0.222414   -78.1   101.2     0.0     0.0      10.3    19.5     0.0     0.0
PHE   -60  -40  2936    3 2 0 0  0.117451   -69.2   -15.4     0.0     0.0       9.3    19.5     0.0     0.0
PHE   -60  -40  2936    2 2 0 0  0.028434  -173.6    23.5     0.0     0.0       9.3    22.3     0.0     0.0
PHE   -60  -40  2936    1 1 0 0  0.014263    76.9    94.0     0.0     0.0       9.0     8.7     0.0     0.0
PHE   -60  -40  2936    1 2 0 0  0.000007    61.2    -7.4     0.0     0.0      12.6    29.3     0.0     0.0
AAAXXBBBBXCCCCXXDDDDXXXXEEEEEEEXXFFFFFFFFXXGGGGGGXXHHHHHHXXIIIIIIXXJJJJJJXXXXKKKKKKXXLLLLLLXXMMMMMMXXNNNNNN
*/
void
dunbrack_print( RotVec  & final_rotamers, RotVec & centroids, std::string aa_name_full )
{
	using namespace std;

	// sanity check
	assert( final_rotamers.size() >= 1 && centroids.size() >= 1 && final_rotamers.size() == centroids.size() );

	// get nchi and num rots
	Size nchi( final_rotamers[1].get_num_chi() );
	Size nrot( final_rotamers.size() );
	Real phi( final_rotamers[1].get_phi() );
	Real psi( final_rotamers[1].get_psi() );
	Real omg( final_rotamers[1].get_omg() );

	// get aa name
	std::string aa_name3( aa_name_full.substr( 0, 3 ) );

	//counts (imaginary and I believe unused)
	Size count(9999);

	// sort final_rotamers based on probability adjust centroids too
	for ( Size i = 1; i <= nrot; ++i ) {
		Size max(i);
		for ( Size j = i+1; j <= nrot; ++j ) {
			if ( final_rotamers[j].get_probability() > final_rotamers[max].get_probability() ) {
				max = j;
			}
		}
		RotData temp_rot( final_rotamers[i] );
		RotData temp_cen( centroids[i] );
		final_rotamers[i] = final_rotamers[max];
		centroids[i] = centroids[max];
		final_rotamers[max] = temp_rot;
		centroids[max] = temp_cen;
	}

	// get "chi bin" number assignments from centroids
	vector1<Size> chi1_bin_nums; vector1<Size> chi2_bin_nums; vector1<Size> chi3_bin_nums; vector1<Size> chi4_bin_nums;
	chi1_bin_nums.resize( nrot, 0 ); chi4_bin_nums.resize( nrot, 0 ); chi3_bin_nums.resize( nrot, 0 ); chi2_bin_nums.resize( nrot, 0 );

	if ( nchi >= 1 ) {
		for ( Size i = 1; i <= nrot; ++i ) {
			chi1_bin_nums[i] = Size( centroids[i].get_lib_chi_val(1) );
		}
	}

	if ( nchi >= 2 ) {
		for ( Size i = 1; i <= nrot; ++i ) {
			chi2_bin_nums[i] = Size( centroids[i].get_lib_chi_val(2) );
		}
	}

	if ( nchi >= 3 ) {
		for ( Size i = 1; i <= nrot; ++i ) {
			chi3_bin_nums[i] = Size( centroids[i].get_lib_chi_val(3) );
		}
	}

	if ( nchi >= 4 ) {
		for ( Size i = 1; i <= nrot; ++i ) {
			chi4_bin_nums[i] = Size( centroids[i].get_lib_chi_val(4) );
		}
	}

	// create filename with string stream
	std::stringstream lib_name;
	lib_name << aa_name3 << "_" << omg << "_" << phi << "_" << psi << "_bbdep.rotlib";

	// open file for writing
	std::ofstream rotlib;
	rotlib.open( lib_name.str().c_str() );

	// now print out lines
	for ( Size i = 1; i <= nrot; ++i ) {
		rotlib << setw(3) << aa_name3
			<< "  "
			<< setw(4) << setprecision(0) << fixed << phi
			<< " "
			<< setw(4) << setprecision(0) << fixed << psi
			<< "  "
			<< setw(4) << setprecision(0) << fixed << omg
			<< "  "
			<< setw(4) << count
			<< "    "
			<< chi1_bin_nums[i] << " " << chi2_bin_nums[i] << " " << chi3_bin_nums[i] << " " << chi4_bin_nums[i]
			<< "  "
			<< setw(8) << setprecision(6) << fixed << final_rotamers[i].get_probability()
			<< "  ";

		for ( Size j = 1; j <= 4; ++j ) {
			if ( j <= nchi ) rotlib << setw(6) << setprecision(1) << fixed << numeric::principal_angle_degrees( final_rotamers[i].get_min_chi(j) );
			else rotlib << setw(6) << setprecision(1) << fixed << "0.0";
			rotlib << "  ";
		}

		rotlib << "  ";

		for ( Size j = 1; j <= 4; ++j ) {
			if ( j <= nchi ) rotlib << setw(6) << setprecision(1) << fixed << final_rotamers[i].get_std_dev(j);
			else rotlib << setw(6) << setprecision(1) << fixed << "0.0";
			rotlib << "  ";
		}

		rotlib << endl;

	}

	rotlib.close();

}


} // namespace make_rot_lib
} // namespace protocols
