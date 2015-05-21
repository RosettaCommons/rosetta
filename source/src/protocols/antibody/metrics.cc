// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/metrics.cc
/// @brief Routines to measure antibodies
/// @author Jeffrey J. Gray
/// @author Oana Lungu
/// @author Nick Marze
/// @author Jared Adolf-Bryfogle

// Project Headers
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
#include <core/scoring/sasa/SasaCalc.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/Energies.hh>


// Numeric Headers
#include <numeric/PCA.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

// Utility Headers
#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "antibody.metrics" );
using namespace core;

namespace protocols {
namespace antibody {


vector1< Real >
vl_vh_orientation_coords ( const pose::Pose & pose_in, const protocols::antibody::AntibodyInfo & ab_info ) {

	vector1< Size > vl_vh_residues = ab_info.get_PackingAngleResidues();

	vector1< numeric::xyzVector< Real > > vl_coord_set;
	for (Size i=1; i<=8; ++i) {
		vl_coord_set.push_back( pose_in.residue( vl_vh_residues[i] ).xyz( "CA" ) );
	}

	vector1< numeric::xyzVector< Real > > vh_coord_set;
	for (Size i=9; i<=16; ++i) {
		vh_coord_set.push_back( pose_in.residue( vl_vh_residues[i] ).xyz( "CA" ) );
	}

	Size vl_n_res = vl_coord_set.size();
	numeric::xyzVector< Real > vl_centroid(0.0);
	for (Size i = 1; i <= vl_n_res; ++i) {
		vl_centroid += vl_coord_set[i];
	}
	vl_centroid /= vl_n_res;

	Size vh_n_res = vh_coord_set.size();
	numeric::xyzVector< Real > vh_centroid(0.0);
	for (Size i = 1; i <= vh_n_res; ++i) {
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
	
	vector1< Real > angle_set;
	angle_set.push_back( vl_vh_distance );
	angle_set.push_back( opening_angle );
	angle_set.push_back( opposite_angle );
	angle_set.push_back( packing_angle );

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
	if (debug) {
		//		std::string kseq = std::string(kr0.name1()) + std::string(kr1.name1()) +
		//                 std::string(kr2.name1()) +	std::string(kr3.name1());
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
	for (std::vector<Vector>::iterator Aa = Aatoms.begin(); Aa != Aatoms.end(); ++Aa) {
		for (std::vector<Vector>::iterator Ca = Catoms.begin(); Ca != Catoms.end(); ++Ca) {
			HBdist = std::min( HBdist, (*Aa - *Ca).norm());
		}
	}
	if (HBdist == 100.0) {
		HBdist = 0;
	}
	return HBdist;
}


core::Real
kink_Trp_Hbond(const core::pose::Pose & pose, const protocols::antibody::AntibodyInfo & ab_info) {

	Size Wi = ab_info.kink_trp(pose); // N+1
	core::conformation::Residue W  = pose.residue(Wi);
	if (W.name3() != "TRP") return 0.0;
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
paratope_sasa( const core::pose::Pose & pose, const protocols::antibody::AntibodyInfo & ab_info){
	using utility::vector1;
	
	ParatopeMetric< core::Real > sasa_results;
	ParatopeMetric< core::Real > hsasa_results;
	
	sasa_results.cdr.resize(6, 0.0);
	hsasa_results.cdr.resize(6, 0.0);
	
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
	for (core::Size i=1; i<=core::Size(ab_info.get_total_num_CDRs()); ++i){
		CDRNameEnum loop = static_cast<CDRNameEnum>(i);
		Size loop_start = ab_info.get_CDR_loop(loop).start();
		Size loop_end = ab_info.get_CDR_loop(loop).stop();
		TR.Debug << loop << " loop_start " << loop_start << std::endl;
		TR.Debug << loop << " loop_end " << loop_end << std::endl;

		///////loop SASA, hydrophobic SASA over paratope/////////////////////////////////////
		Real hydrop_loop_sasa = 0;
		Real loop_sasa = 0;
		Real polar_loop_sasa = 0;
		for (Size ii=loop_start; ii<=loop_end; ++ii){
			core::conformation::Residue const & irsd( pose.residue( ii ) );
			loop_sasa += res_sasa[ii];
			paratope_sasa += res_sasa[ii];
			hydrop_loop_sasa += res_hsasa[ii];
			hydrophobic_sasa += res_hsasa[ii];
			TR.Debug << "residue " << irsd.name3() << ii << std::endl;
		}
		polar_loop_sasa = loop_sasa - hydrop_loop_sasa;
		
		sasa_results.cdr[i] = loop_sasa;
		hsasa_results.cdr[i] = hydrop_loop_sasa;
		
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
	for(Size i(1); i<=pose.total_residue(); ++i) {
		std::string name3 = pose.residue(i).name3();
		if(name3 == "ARG" || name3 == "LYS") {
			//TR << name3 << std::endl;
			pose_net_charge++;
		}
		else if(name3 == "ASP" || name3 == "GLU") {
			pose_net_charge--;
			//TR << name3 << std::endl;
		}
	}
	return pose_net_charge;
}

ParatopeMetric<core::SSize>
paratope_charge( core::pose::Pose const & pose, const protocols::antibody::AntibodyInfo & ab_info ) {
	using namespace core;
	
	ParatopeMetric<core::SSize> charge_results;
	charge_results.cdr.resize(6, 0);
	core::SSize paratope_charge = 0;
	
	for(core::Size ia=1; ia<=core::Size(ab_info.get_total_num_CDRs()); ++ia){
		CDRNameEnum loop = static_cast<CDRNameEnum>(ia);
		Size loop_start = ab_info.get_CDR_loop(loop).start();
		Size loop_end = ab_info.get_CDR_loop(loop).stop();

		SSize plus_charge(0);
		SSize minus_charge(0);
		SSize loop_charge(0);

		for(Size ib=loop_start; ib<=loop_end; ++ib) {
			std::string name3 = pose.residue(ib).name3();
			if(name3 == "ARG" || name3 == "LYS") {
				plus_charge++;
			}
			else if(name3 == "ASP" || name3 == "GLU") {
				minus_charge--;
			}
		}
		loop_charge = plus_charge + minus_charge;
		paratope_charge += loop_charge;
		charge_results.cdr[ia] = loop_charge;
		
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
	for( core::Size i = cdr_start; i <= cdr_end; ++i){
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

}
}
