// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/metrics.cc
/// @brief Routines to measure antibodies
/// @author Jeffrey J. Gray
/// @author Oana Lungu
/// @author Nick Marze
/// @author Jared Adolf-Bryfogle

// Project Headers
#include <protocols/antibody/util.hh>
#include <protocols/antibody/metrics.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/AntibodyEnum.hh>

// Basic
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>

// Core
//#include <core/pose/metrics/simple_calculators/SasaCalculator2.hh>
#include <core/id/SequenceMapping.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/sasa/SasaCalc.hh>
#include <core/scoring/rms_util.hh>


// Numeric Headers
#include <numeric/PCA.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

// Utility Headers
#include <utility/vector1.hh>


static basic::Tracer TR( "antibody.metrics" );
using namespace core;
using utility::vector1;

namespace protocols {
namespace antibody {


vector1< Real >
vl_vh_orientation_coords ( const pose::Pose & pose_in, const protocols::antibody::AntibodyInfo & ab_info ) {
	vector1< Real > angle_set(4, 0.0);
	vector1< Size > vl_vh_residues = ab_info.get_PackingAngleResidues();

	//Check if camelid antibody.
	if ( ab_info.is_camelid() ) {
		TR << "Could not obtain Vl Vh coords as antibody is a camelid antibody.  Setting coords to zero" << std::endl;
		return angle_set;
	}



	vector1< numeric::xyzVector< Real > > vl_coord_set;
	for ( Size i=1; i<=8; ++i ) {
		vl_coord_set.push_back( pose_in.residue( vl_vh_residues[i] ).xyz( "CA" ) );
	}

	vector1< numeric::xyzVector< Real > > vh_coord_set;
	for ( Size i=9; i<=16; ++i ) {
		vh_coord_set.push_back( pose_in.residue( vl_vh_residues[i] ).xyz( "CA" ) );
	}

	Size vl_n_res = vl_coord_set.size();
	numeric::xyzVector< Real > vl_centroid(0.0);
	for ( Size i = 1; i <= vl_n_res; ++i ) {
		vl_centroid += vl_coord_set[i];
	}
	vl_centroid /= vl_n_res;

	Size vh_n_res = vh_coord_set.size();
	numeric::xyzVector< Real > vh_centroid(0.0);
	for ( Size i = 1; i <= vh_n_res; ++i ) {
		vh_centroid += vh_coord_set[i];
	}
	vh_centroid /= vh_n_res;

	numeric::xyzVector< Real > vl_first_principal_component = numeric::first_principal_component( vl_coord_set );
	numeric::xyzVector< Real > vh_first_principal_component = numeric::first_principal_component( vh_coord_set );

	numeric::xyzVector< Real > point_one = vl_centroid + vl_first_principal_component;
	numeric::xyzVector< Real > point_one_prime = vl_centroid - vl_first_principal_component;
	numeric::xyzVector< Real > point_four = vh_centroid + vh_first_principal_component;
	numeric::xyzVector< Real > point_four_prime = vh_centroid - vh_first_principal_component;

	numeric::xyzVector< Real > top_point_one;
	numeric::xyzVector< Real > top_point_four;

	if ( point_one.distance( vl_coord_set[1] ) > point_one_prime.distance( vl_coord_set[1] ) ) {
		top_point_one = point_one_prime;
	} else {
		top_point_one = point_one;
	}

	if ( point_four.distance( vh_coord_set[1] ) > point_four_prime.distance( vh_coord_set[1] ) ) {
		top_point_four = point_four_prime;
	} else {
		top_point_four = point_four;
	}

	//TR << "Point 1: " << top_point_one << std::endl;
	//TR << "Point 2: " << vl_centroid << std::endl;
	//TR << "Point 3: " << vh_centroid << std::endl;
	//TR << "Point 4: " << top_point_four << std::endl;

	Real packing_angle = numeric::dihedral_degrees( top_point_one, vl_centroid, vh_centroid, top_point_four );
	Real opening_angle = numeric::angle_degrees( vl_centroid, vh_centroid, top_point_four );
	Real opposite_angle = numeric::angle_degrees( top_point_one, vl_centroid, vh_centroid );
	Real vl_vh_distance = vl_centroid.distance( vh_centroid );


	angle_set[ 1 ] = vl_vh_distance ;
	angle_set[ 2 ] = opening_angle ;
	angle_set[ 3 ] = opposite_angle ;
	angle_set[ 4 ] =  packing_angle ;

	return angle_set;
}


////////////////// Kink Measures ////////////////

std::pair<core::Real,core::Real>
kink_dihedral(const core::pose::Pose & pose, const protocols::antibody::AntibodyInfo & ab_info, bool debug) {

	Size kb = ab_info.kink_begin(pose);
	core::conformation::Residue kr0 = pose.residue(kb);
	core::conformation::Residue kr1 = pose.residue(kb+1);
	core::conformation::Residue kr2 = pose.residue(kb+2);
	core::conformation::Residue kr3 = pose.residue(kb+3);
	if ( debug ) {
		//  std::string kseq = std::string(kr0.name1()) + std::string(kr1.name1()) +
		//                 std::string(kr2.name1()) + std::string(kr3.name1());
		TR << "Kink is defined from pose residues " << kb << "-" << kb+3 << ": " << kr0.name1() << std::endl;
	}

	Vector CA0=kr0.xyz("CA");
	Vector CA1=kr1.xyz("CA");
	Vector CA2=kr2.xyz("CA");
	Vector CA3=kr3.xyz("CA");

	Real q = (CA0 - CA3).magnitude();
	Real qbase = dihedral(CA0,CA1,CA2,CA3);

	std::pair<core::Real,core::Real> qout(q,qbase);

	return qout;
}

core::Real
kink_bb_Hbond(const core::pose::Pose & pose, const protocols::antibody::AntibodyInfo & ab_info) {

	Size Di = ab_info.kink_anion_residue(pose);  // N-1
	core::conformation::Residue D  = pose.residue(Di);
	Vector DN = D.xyz("N");

	Size Ri = ab_info.kink_cation_residue(pose);  // 0
	core::conformation::Residue R  = pose.residue(Ri);
	Vector RO = R.xyz("O");

	TR << "H3_DN   (" << Di << "): " << D.name3() << " - " << DN << std::endl;
	TR << "H3_RO   (" << Ri << "): " << R.name3() << " - " << RO << std::endl;

	core::Real bbHBdist = ( DN - RO ).norm();

	return bbHBdist;
}


core::Real
kink_RD_Hbond(const core::pose::Pose & pose, const protocols::antibody::AntibodyInfo & ab_info) {

	std::vector<Vector> Aatoms = ab_info.kink_anion_atoms(pose); // N-1
	std::vector<Vector> Catoms = ab_info.kink_cation_atoms(pose); // 0

	Real HBdist = 100.0;
	for ( auto & Aatom : Aatoms ) {
		for ( auto & Catom : Catoms ) {
			HBdist = std::min( HBdist, (Aatom - Catom).norm());
		}
	}
	if ( HBdist == 100.0 ) {
		HBdist = 0;
	}
	return HBdist;
}


core::Real
kink_Trp_Hbond(const core::pose::Pose & pose, const protocols::antibody::AntibodyInfo & ab_info) {

	Size Wi = ab_info.kink_trp(pose); // N+1
	core::conformation::Residue W  = pose.residue(Wi);
	if ( W.name3() != "TRP" ) return 0.0;
	Vector W_NE1 = W.xyz("NE1");

	Size kb1 = ab_info.kink_begin(pose); // N-2
	core::conformation::Residue kb = pose.residue(kb1);
	Vector kb1_O = kb.xyz("O");

	TR << "H3_W/KB4 (" << Wi  << "): " << W.name3()  << " - " << W_NE1 << std::endl;
	TR << "H3_KB1   (" << kb1 << "): " << kb.name3() << " - " << kb1_O << std::endl;

	core::Real WHBdist = ( W_NE1 - kb1_O ).norm();

	return WHBdist;
}


//////////////////// SASA /////////////////


std::pair<ParatopeMetric< core::Real >, ParatopeMetric< core::Real> >
paratope_sasa( const core::pose::Pose & pose, const protocols::antibody::AntibodyInfo & ab_info, bool include_de_loop /* false */){
	using utility::vector1;

	ParatopeMetric< core::Real > sasa_results;
	ParatopeMetric< core::Real > hsasa_results;

	sasa_results.cdr.resize(ab_info.get_total_num_CDRs(include_de_loop), 0.0);
	hsasa_results.cdr.resize(ab_info.get_total_num_CDRs(include_de_loop), 0.0);

	///Create a basic new sasa calc, and get needed values.
	scoring::sasa::SasaCalcOP sasa_calc( new scoring::sasa::SasaCalc() );

	Real total_sasa = sasa_calc->calculate(pose);
	Real total_hsasa = sasa_calc->get_total_hsasa();

	vector1< Real > res_sasa = sasa_calc->get_residue_sasa();
	vector1< Real > res_hsasa = sasa_calc->get_residue_hsasa();

	Real paratope_sasa = 0;
	Real hydrophobic_sasa = 0;
	Real total_polar_sasa = 0;
	Real paratope_polar_sasa = 0;

	/////////////////////////////////////total SASA////////////////////////////////////////
	TR << "Total SASA is: " << total_sasa << std::endl;
	TR << "Total hydrophobic SASA is: " << total_hsasa << std::endl;
	total_polar_sasa = total_sasa - total_hsasa;
	TR << "Total polar SASA is: " << total_polar_sasa << std::endl;

	////////iterate to define antibody paratope/////////////////////////////////////////////
	for ( auto const & loop : ab_info.get_all_cdrs_present(include_de_loop) ) {
		Size loop_start = ab_info.get_CDR_start(loop, pose);
		Size loop_end = ab_info.get_CDR_end(loop, pose);
		TR.Debug << loop << " loop_start " << loop_start << std::endl;
		TR.Debug << loop << " loop_end " << loop_end << std::endl;

		///////loop SASA, hydrophobic SASA over paratope/////////////////////////////////////
		Real hydrop_loop_sasa = 0;
		Real loop_sasa = 0;
		Real polar_loop_sasa = 0;
		for ( Size ii=loop_start; ii<=loop_end; ++ii ) {
			core::conformation::Residue const & irsd( pose.residue( ii ) );
			loop_sasa += res_sasa[ii];
			paratope_sasa += res_sasa[ii];
			hydrop_loop_sasa += res_hsasa[ii];
			hydrophobic_sasa += res_hsasa[ii];
			TR.Debug << "residue " << irsd.name3() << ii << std::endl;
		}
		polar_loop_sasa = loop_sasa - hydrop_loop_sasa;

		sasa_results.cdr[loop] = loop_sasa;
		hsasa_results.cdr[loop] = hydrop_loop_sasa;

		/////////////////////out CDR values//////////////////////////////////////////////
		TR << "Loop " << ab_info.get_CDR_name(loop) << ": CDR  sasa " << loop_sasa

			<< "\tCDR hydrophobic sasa: " << hydrop_loop_sasa
			<< "\tCDR polar sasa:  " << polar_loop_sasa << std::endl;
		///////////////////////////////out paratope values////////////////////////////////////
	}
	TR << "Paratope_sasa " << paratope_sasa << std::endl;
	TR << "Paratope hydrophobic sasa " << hydrophobic_sasa << std::endl;
	paratope_polar_sasa = paratope_sasa - hydrophobic_sasa;
	TR << "Paratope polar sasa " << paratope_polar_sasa << std::endl;

	sasa_results.paratope = paratope_sasa;
	hsasa_results.paratope = hydrophobic_sasa;

	std::pair<ParatopeMetric<core::Real>, ParatopeMetric<core::Real> > sasa_out(sasa_results , hsasa_results);
	return sasa_out;
}

core::SSize
pose_charge( core::pose::Pose const & pose ) {  // not really an antibody fn
	SSize pose_net_charge(0);
	for ( Size i(1); i<=pose.size(); ++i ) {
		std::string name3 = pose.residue(i).name3();
		if ( name3 == "ARG" || name3 == "LYS" ) {
			//TR << name3 << std::endl;
			pose_net_charge++;
		} else if ( name3 == "ASP" || name3 == "GLU" ) {
			pose_net_charge--;
			//TR << name3 << std::endl;
		}
	}
	return pose_net_charge;
}

ParatopeMetric<core::SSize>
paratope_charge( core::pose::Pose const & pose, const protocols::antibody::AntibodyInfo & ab_info, bool include_de_loop /* false */ ) {
	using namespace core;

	ParatopeMetric<core::SSize> charge_results;
	charge_results.cdr.resize(ab_info.get_total_num_CDRs( include_de_loop ), 0);
	core::SSize paratope_charge = 0;

	for ( auto const & loop : ab_info.get_all_cdrs_present( include_de_loop ) ) {
		Size loop_start = ab_info.get_CDR_start(loop, pose);
		Size loop_end = ab_info.get_CDR_end(loop, pose);

		SSize plus_charge(0);
		SSize minus_charge(0);
		SSize loop_charge(0);

		for ( Size ib=loop_start; ib<=loop_end; ++ib ) {
			std::string name3 = pose.residue(ib).name3();
			if ( name3 == "ARG" || name3 == "LYS" ) {
				plus_charge++;
			} else if ( name3 == "ASP" || name3 == "GLU" ) {
				minus_charge--;
			}
		}
		loop_charge = plus_charge + minus_charge;
		paratope_charge += loop_charge;
		charge_results.cdr[loop] = loop_charge;

		TR << "Loop " << ab_info.get_CDR_name(loop) << ": "
			<< std::setw(3) << minus_charge << " anions, "
			<< std::setw(3) << plus_charge  << " cations = "
			<< std::setw(3) << loop_charge  << " total charge" << std::endl;
	}
	charge_results.paratope = paratope_charge;
	return charge_results;
}

core::Real
cdr_energy(core::pose::Pose const & pose, AntibodyInfoCOP ab_info, core::scoring::ScoreFunctionCOP scorefxn, CDRNameEnum const & cdr){

	core::pose::Pose new_pose = pose;

	//Allow Hbonding scores in pose
	core::scoring::ScoreFunctionOP new_scorefxn = scorefxn->clone();
	core::scoring::methods::EnergyMethodOptions options(new_scorefxn->energy_method_options());
	options.hbond_options().decompose_bb_hb_into_pair_energies(true);
	new_scorefxn->set_energy_method_options(options);
	(*new_scorefxn)(new_pose);

	core::Size cdr_start = ab_info->get_CDR_start(cdr, pose);
	core::Size cdr_end  = ab_info->get_CDR_end(cdr, pose);
	core::Real energy = 0.0;
	for ( core::Size i = cdr_start; i <= cdr_end; ++i ) {
		energy = energy + new_pose.energies().residue_total_energy(i);
	}

	return energy;
}

core::Real
cdr_CN_anchor_distance(core::pose::Pose const & pose, AntibodyInfoCOP ab_info, CDRNameEnum const & cdr) {
	core::Size n_term_res = ab_info->get_CDR_start(cdr, pose) - 1;
	core::Size c_term_res = ab_info->get_CDR_end(cdr, pose) +1;

	numeric::xyzVector<core::Real> n_pos= pose.residue(n_term_res).xyz("C");
	numeric::xyzVector<core::Real> c_pos = pose.residue(c_term_res).xyz("N");

	return n_pos.distance(c_pos);
}

vector1< Real >
cdr_backbone_rmsds(pose::Pose & p1, pose::Pose const & p2, AntibodyInfoOP const a1, AntibodyInfoOP const a2, Size aln_thresh /*10*/) {
	// vector to store results h1, h2, h3, l1, l2, l3, ocd
	vector1< Real > results_set(9, 0.0);

	// calculate OCD before moving anything
	vector1<Real> p1_OCDs = vl_vh_orientation_coords(p1, *a1);
	vector1<Real> p2_OCDs = vl_vh_orientation_coords(p2, *a2);
	// calculate OCD
	/*
	Calculated as a z-score. Mean/Stdev are from Nick's paper:
	academic.oup.com/peds/article/29/10/409/2462315

	Figure 3 is ths source of the below values:
	distance:    14.6 +/- 0.32
	H open:      97.2 +/- 2.55
	L open:      99.4 +/- 1.93
	pack angle: -52.3 +/- 3.83
	*/
	// is this the best practice? we need the st. devs calculated from database!
	//utility::vector1<core::Real> means = {14.6, 97.2, 99.4, -52.3};
	vector1<Real> sds = {0.32, 2.55, 1.93, 3.83};
	Real ocd = 0.0;

	//TR << "OCD X: mobile fixed" << std::endl;
	for ( Size i=1; i<=4; ++i ) {
		//TR << "OCD " << i << ": " << p1_OCDs[i] << " " << p2_OCDs[i] << std::endl;
		ocd += ( (p1_OCDs[i] - p2_OCDs[i]) * (p1_OCDs[i] - p2_OCDs[i]) / sds[i]*sds[i] );
	}
	// ocd is sum of squared z-scores, so return square root
	// Nick didn't do this in his paper, so a direct comparison is not possible
	results_set[1] = std::sqrt(ocd);

	// now for light/heavy RMSDs: p2 is fixed... p1 will move
	TR << "Pose1: " << p1.sequence() << std::endl;
	TR << "Pose2: " << p2.sequence() << std::endl;

	// for the alignment generate sequence maps, then use the maps to map the frameworks
	// this allows for residue mismatches here and there, but warn if there are too many
	id::SequenceMapping seq_map = sequence::map_seq1_seq2( utility::pointer::make_shared< core::sequence::Sequence >(p1), utility::pointer::make_shared< sequence::Sequence >(p2) );

	// construct atom maps of frh/frl/cdr residues from sequence mapping
	std::map< CDRNameEnum, std::map< id::AtomID, id::AtomID > > cdr_maps;
	std::map< char, id::AtomID_Map< id::AtomID > > fr_maps;

	for ( auto const cdr_enum : {h1, h2, h3, l1, l2, l3} ) {
		std::map< id::AtomID, id::AtomID > atomid_map;
		cdr_maps[cdr_enum] = atomid_map;
	}

	for ( char const c : {'H', 'L'} ) {
		id::AtomID_Map< id::AtomID > atomid_map;
		initialize_atomid_map( atomid_map, p1, id::AtomID::BOGUS_ATOM_ID() );
		fr_maps[c] = atomid_map;
	}

	Size n_zero = 0;

	// define conserved residues for alignment in Chothia numbering
	for ( auto i : get_conserved_residue_list('H') ) {

		Size k = p1.pdb_info()->pdb2pose('H', i);
		if ( k == 0 ) { // res doesnt exist in p1
			n_zero += 1;
			continue;
		}

		Size j = seq_map[k]; // map p1 to p2
		if ( j==0 ) { // res doesn't exist in p2
			n_zero += 1;
			continue;
		}

		TR.Debug << "Mapping H: " << i << " to " << p2.pdb_info()->pose2pdb(j) << std::endl;
		// map N, CA, C, O
		for ( std::string const atom_name : {"N", "CA", "C", "O"} ) {
			id::AtomID const id1( p1.residue(k).atom_index(atom_name), k );
			id::AtomID const id2( p2.residue(j).atom_index(atom_name), j );
			fr_maps['H'].set(id1, id2);
		}
	}

	for ( auto i : get_conserved_residue_list('L') ) {

		Size k = p1.pdb_info()->pdb2pose('L', i);
		if ( k == 0 ) { // res doesnt exist in p1
			n_zero += 1;
			continue;
		}

		Size j = seq_map[k]; // map p1 to p2
		if ( j==0 ) {
			n_zero += 1;
			continue;
		}

		TR.Debug << "Mapping L:" << i << " to " << p2.pdb_info()->pose2pdb(j) << std::endl;
		// map N, CA, C, O
		for ( std::string const atom_name : {"N", "CA", "C", "O"} ) {
			id::AtomID const id1( p1.residue(k).atom_index(atom_name), k );
			id::AtomID const id2( p2.residue(j).atom_index(atom_name), j );
			fr_maps['L'].set(id1, id2);
		}
	}

	// loop over pose residues, check which region, generate atom maps for region
	for ( Size i=1; i<=p1.size(); ++i ) {
		Size j = seq_map[i]; // map p1 to p2
		if ( a1->get_region_of_residue(p1, i) == cdr_region ) { // in CDR
			// we must always map or something is very wrong
			if ( j==0 ) {
				throw CREATE_EXCEPTION(utility::excn::BadInput, "Could not map the first pose to the second pose in the CDRs, something is very wrong.");
			}
			CDRNameEnum region = a1->get_CDRNameEnum_of_residue(p1, i);
			// map N, CA, C, O
			TR.Debug << "CDR : " << region << std::endl;
			TR.Debug << p1.pdb_info()->pose2pdb(i) << std::endl;
			TR.Debug << p2.pdb_info()->pose2pdb(j) << std::endl;
			for ( std::string const atom_name : {"N", "CA", "C", "O"} ) {
				id::AtomID const id1( p1.residue(i).atom_index(atom_name), i );
				id::AtomID const id2( p2.residue(j).atom_index(atom_name), j );
				cdr_maps[region][id1] = id2;
			}
		}
	}
	//TR << cdr_maps << std::endl;

	if ( n_zero > aln_thresh ) {
		// or do we throw an exception?
		TR << "With an alignment threshold of " << aln_thresh << ", there were too many (" << n_zero << ") mismatches." << std::endl;
		//throw CREATE_EXCEPTION(utility::excn::BadInput, "Too many mismatches between the input poses. Alignment would be poor");
	}

	// align heavy, calc CDR H1, H2, H3 rmsd
	results_set[2] = scoring::superimpose_pose( p1, p2, fr_maps['H'] );
	results_set[3] = scoring::rms_at_corresponding_atoms_no_super( p1, p2, cdr_maps[h1] );
	results_set[4] = scoring::rms_at_corresponding_atoms_no_super( p1, p2, cdr_maps[h2] );
	results_set[5] = scoring::rms_at_corresponding_atoms_no_super( p1, p2, cdr_maps[h3] );

	// align light, calc CDR L1, L2, L3 rmsd
	if ( a1->is_camelid() == false && a2->is_camelid() == false ) {
		results_set[6] = core::scoring::superimpose_pose( p1, p2, fr_maps['L'] );
		results_set[7] = scoring::rms_at_corresponding_atoms_no_super( p1, p2, cdr_maps[l1] );
		results_set[8] = scoring::rms_at_corresponding_atoms_no_super( p1, p2, cdr_maps[l2] );
		results_set[9] = scoring::rms_at_corresponding_atoms_no_super( p1, p2, cdr_maps[l3] );
	}

	return results_set;
}

}
}
