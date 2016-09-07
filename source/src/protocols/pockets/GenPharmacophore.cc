// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pockets/Fingerprint.cc
/// @brief  protocols::pockets::Fingerprint functions
/// @author Ragul Gowthaman

// Protocol Headers

#include <core/pose/Pose.hh>
#include <basic/MetricValue.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/fingerprint.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/after_opts.hh>

#include <protocols/simple_moves/ScoreMover.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/rna/RNA_ResidueType.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <basic/options/keys/gen_pharmacophore.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <protocols/pockets/GenPharmacophore.hh>
#include <core/chemical/rna/util.hh>
//#include <core/pose/rna/RNA_Util.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/io/ozstream.hh>

// C++ Headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

//Auto Headers

#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.fwd.hh>

// Boost headers
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


using namespace core;
using namespace std;


namespace protocols {
namespace pockets {

SmallMol::SmallMol(const SmallMol &other) {
	*this = other;
	parent = this;
}

SmallMol::~SmallMol() = default;

void SmallMol::add_atom(string line) {
	pdbContent.append(line + "\n");
	vector<double> newAtom(3);

	double x = atof(line.substr(30, 8).c_str());
	double y = atof(line.substr(38, 8).c_str());
	double z = atof(line.substr(46, 8).c_str());
	newAtom[0] = x;
	newAtom[1] = y;
	newAtom[2] = z;

	coordinates.push_back(newAtom);
	newAtom.clear();
}

void SmallMol::update_center() {
	int size = (int)coordinates.size();
	core::Real total_x = 0.0, total_y = 0.0, total_z = 0.0;
	for ( int i = 0; i < size; i++ ) {
		vector<core::Real> &current_atom = coordinates[i];
		total_x += current_atom[0];
		total_y += current_atom[1];
		total_z += current_atom[2];
	}

	cen[0] = total_x / (core::Real)size;
	cen[1] = total_y / (core::Real)size;
	cen[2] = total_z / (core::Real)size;
}

core::Real calDist(vector<core::Real> const &a1, vector<core::Real> const &a2) {
	if ( a1.size() != 3 || a2.size() != 3 ) {
		cout << "wrong size in vector";
		exit(-1);
	}
	core::Real total = 0.0;
	//TODO: vectorization
	for ( int i = 0; i < 3; i++ ) {
		total += (a1[i] - a2[i]) * (a1[i] - a2[i]);
	}
	return total;
}

core::Real SmallMol::calRMSD(SmallMol &mol1, SmallMol &mol2) {
	vector< vector<core::Real> > coord1 = mol1.coordinates;
	vector< vector<core::Real> > coord2 = mol2.coordinates;
	core::Real sumRMSD = 0.0;

	for (auto & i : coord1) {
		core::Real minDist = 999999999.0;
		for (auto & j : coord2) {
			core::Real currentDist = calDist(i, j);
			if ( currentDist < minDist ) {
				minDist = currentDist;
			}
		}
		sumRMSD += minDist;
	}
	mol2.rmsd = sumRMSD;
	return sumRMSD;
}

void SmallMol::printCoordinates() const {
	cout << molName << endl;
	for (auto currentAtom : coordinates) {
			for (double k : currentAtom) {
			cout << k << "\t";
		}
		cout << endl;
	}
}

void SmallMol::printContent() const {
	cout << pdbContent;
}

string SmallMol::getContent() const {
	return pdbContent;
}

core::Real SmallMol::get_center(int c) {
	return cen[c];
}

core::Real SmallMol::cal_distance(SmallMol *other) {
	// TODO: vectorization
	core::Real d_x = cen[0] - other->get_center(0);
	core::Real d_y = cen[1] - other->get_center(1);
	core::Real d_z = cen[2] - other->get_center(2);
	return (d_x * d_x + d_y * d_y + d_z * d_z) / 3;
}

core::Real SmallMol::cal_min_dist(SmallMol *other) {
	vector< vector<core::Real> > const &coord1 = coordinates;
	vector< vector<core::Real> > const &coord2 = other->get_coordinates();
	core::Real min_dist = 9999999999.9;

	for (const auto & i : coord1) {
		for (const auto & j : coord2) {
			core::Real currentDist = calDist(i, j);
			if ( currentDist < min_dist ) {
				min_dist = currentDist;
			}
		}
	}

	return min_dist;
}

SmallMol * SmallMol::findRoot() {
	SmallMol *u = this;
	while ( u != u->get_parent() )
			u = u->get_parent();

	return u;
}

bool SmallMol::connected(SmallMol *m) {
	SmallMol *r_m = m->findRoot();
	SmallMol *r_t = findRoot();
	return r_m == r_t;
}

void SmallMol::connect(SmallMol *m) {
	SmallMol *r_m = m->findRoot();
	SmallMol *r_t = findRoot();
	int size_m = r_m->get_size();
	int size_t = r_t->get_size();

	if ( size_m >= size_t ) {
		r_t->set_parent(r_m);
		r_m->set_size(size_m + size_t);
	} else {
		r_m->set_parent(r_t);
		r_t->set_size(size_m + size_t);
	}
}

bool GenPharmacophore::is_buried_ring(core::conformation::Residue const & rsd, core::Real const & ring_sasa, core::Real const & sasa_cutoff){
	if ( rsd.name3() == "  A" ) return ring_sasa <= sasa_cutoff; //46.81 ;
	else if ( rsd.name3() == "  C" )  return ring_sasa <= sasa_cutoff; //31.09 ;
	else if ( rsd.name3() == "  G" )  return ring_sasa <= sasa_cutoff; //45.06 ;
	else if ( rsd.name3() == "  U" )  return ring_sasa <= sasa_cutoff; //52.66 ;
	else return false;
}

core::Real GenPharmacophore::get_RNAring_sasa( core::conformation::Residue const & rsd, int const & rsdno, core::id::AtomID_Map<core::Real> const & pose_atom_sasa ) {
	if ( !rsd.is_RNA() ) {
		std::cout<<"Residue is not an RNA base. Cannot calculate RNA base SASA."<<std::endl;
		exit (1);
	}
	core::chemical::rna::RNA_ResidueType const & rsd_type = rsd.RNA_type();
	core::Real curr_ring_sasa = 0;
	for ( Size rj = 1, rj_end = rsd.nheavyatoms(); rj <= rj_end; ++rj ) {
		if ( !rsd_type.is_RNA_base_atom(rj) ) continue;
		id::AtomID const aid( rj, rsdno);
		curr_ring_sasa += pose_atom_sasa[aid];
	}
	return(curr_ring_sasa);
}

void GenPharmacophore::get_ideal_hydrogenBond_atoms(core::pose::Pose const & protein_pose){
	using namespace basic::options;
	core::Real const clash_dist = option[ OptionKeys::gen_pharmacophore::clash_distance_cutoff ];

	//ideal H-bond distance
	core::Real const opt_distance( 2.75 );
	//core::Real const distance( 3.5 );
	using namespace core::chemical;
	using namespace core::kinematics;

	// We only want to find donors/acceptors that are solvent accessible

	core::pose::metrics::PoseMetricCalculatorOP res_sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy );
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasaone", res_sasa_calculator );
	basic::MetricValue< utility::vector1< core::Real > > ressasa;
	protein_pose.metric( "sasaone", "residue_sasa", ressasa );
	core::pose::metrics::PoseMetricCalculatorOP atm_sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy );
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasatwo", atm_sasa_calculator );
	basic::MetricValue< core::id::AtomID_Map< core::Real> > atmsasa;
	protein_pose.metric( "sasatwo", "atom_sasa", atmsasa );
	core::id::AtomID_Map< core::Real > atom_sasas = atmsasa.value();

	std::list< numeric::xyzVector<core::Real> > dnr_coord_list;
	std::list< numeric::xyzVector<core::Real> > acp_coord_list;
	for ( Size j = 1, resnum = protein_pose.total_residue(); j <= resnum; ++j ) {
		core::conformation::Residue const & rsd( protein_pose.conformation().residue(j) );

		//int target=0;

		//fill in points that are ideal for a hydrogen acceptor with an O
		for ( auto hnum  = rsd.Hpos_polar().begin(), hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
			Size const hatm( *hnum );
			// Skip buried residues
			if ( atom_sasas(j, hatm) < 0.1 && atom_sasas(j, rsd.atom_base(hatm)) < 0.1 ) {
				//std::cout<<rsd.seqpos()+offset<<" Donor "<<rsd.name()<<" "<<rsd.atom_name(rsd.atom_base(hatm))<<" H SASA "<<atom_sasas(j, hatm)<<" Base SASA "<<atom_sasas(j, rsd.atom_base(hatm))<<" being ignored"<<std::endl;
				continue;
			}
			//std::cout<<rsd.seqpos()+offset<<" Donor "<<rsd.name()<<" "<<rsd.atom_name(rsd.atom_base(hatm))<<" H SASA "<<atom_sasas(j, hatm)<<" Base SASA "<<atom_sasas(j, rsd.atom_base(hatm))<<std::endl;

			numeric::xyzVector<core::Real> const & hatm_xyz( rsd.xyz( hatm ) );
			numeric::xyzVector<core::Real> const & datm_xyz( rsd.xyz( rsd.atom_base( hatm ) ) );
			for ( int step = 0; step<4; step++ ) {
				numeric::xyzVector<core::Real> const ro1(datm_xyz + opt_distance * ( hatm_xyz - datm_xyz ).normalized());
				std::cout<<"Optimal acceptor at "<<ro1.x()<<" "<<ro1.y()<<" "<<ro1.z()<<std::endl;
				acp_coord_list.push_back(ro1);
			}
		}

		//fill in points that are ideal for a hydrogen donor with an N
		for ( auto anum  = rsd.accpt_pos().begin(), anume = rsd.accpt_pos().end(); anum != anume; ++anum ) {
			Size const aatm( *anum );
			// Skip buried residues
			if ( atom_sasas(j, aatm) < 0.1 ) {
				//std::cout<<rsd.seqpos()+offset<<" Acceptor "<<rsd.name()<<" "<<rsd.atom_name(aatm)<<" SASA "<<atom_sasas(j, aatm)<<" being ignored"<<std::endl;
				continue;
			}
			//std::cout<<rsd.seqpos()+offset<<" Acceptor "<<rsd.name()<<" "<<rsd.atom_name(aatm)<<" SASA "<<atom_sasas(j, aatm)<<std::endl;

			numeric::xyzVector<core::Real> const & aatm_xyz( rsd.xyz( aatm ) );
			numeric::xyzVector<core::Real> aatm_base_xyz( rsd.xyz( rsd.atom_base( aatm ) ) );
			numeric::xyzVector<core::Real> const & aatm_base2_xyz( rsd.xyz( rsd.abase2( aatm ) ) );
			Hybridization const & hybrid( rsd.atom_type(aatm).hybridization() );

			core::Real theta(0.0);
			utility::vector1< core::Real > phi_list, phi_steps;
			phi_steps.push_back(  0 );
			switch( hybrid ) {
			case SP2_HYBRID :
				theta = 180.0 - 120.0;
				phi_list.push_back(   0.0 );
				phi_list.push_back( 180.0 );
				break;
			case SP3_HYBRID :
				theta = 180.0 - 109.0;
				phi_list.push_back( 120.0 );
				phi_list.push_back( 240.0 );
				break;
			case RING_HYBRID :
				{
				numeric::xyzVector<core::Real> const & avg_base_xyz (0.5 * ( aatm_base_xyz + aatm_base2_xyz ));
				aatm_base_xyz(1)=avg_base_xyz(1);
				aatm_base_xyz(2)=avg_base_xyz(2);
				aatm_base_xyz(3)=avg_base_xyz(3);
				theta = 0.0;
				phi_steps.clear();
				phi_steps.push_back( 0.0 );
				phi_list.push_back( 0.0 ); // doesnt matter
				break;
			}
			default :
				std::cerr << "Bad hybridization type for acceptor " << hybrid << '\n';
				exit(1000);
			}
			Stub stub( aatm_xyz, aatm_base_xyz, aatm_base2_xyz );
			for ( Size i=1; i<= phi_list.size(); ++i ) {
				numeric::xyzVector<core::Real> const ro1(stub.spherical( numeric::conversions::radians( phi_list[i]), numeric::conversions::radians( theta ), opt_distance));
				std::cout<<"Optimal donor at "<<ro1.x()<<" "<<ro1.y()<<" "<<ro1.z()<<std::endl;
				dnr_coord_list.push_back(ro1);
			}
		}
	}
	numeric::xyzVector<core::Real> protein_atom_coord(0.);
	for ( int j = 1, resnum = protein_pose.total_residue(); j <= resnum; ++j ) {
		core::conformation::Residue const & curr_rsd = protein_pose.residue(j);
		for ( Size i = 1, i_end = curr_rsd.natoms(); i <= i_end; ++i ) {
			protein_atom_coord.x() = curr_rsd.atom(i).xyz()(1);
			protein_atom_coord.y() = curr_rsd.atom(i).xyz()(2);
			protein_atom_coord.z() = curr_rsd.atom(i).xyz()(3);
			for ( auto aa = acp_coord_list.begin(); aa != acp_coord_list.end(); /* ++aa */ ) {
				if ( protein_atom_coord.distance(*aa) <= clash_dist ) {
					aa = acp_coord_list.erase(aa);
				} else {
					++aa;
				}
			}
			for ( auto bb = dnr_coord_list.begin(); bb != dnr_coord_list.end(); /* ++bb */ ) {
				if ( protein_atom_coord.distance(*bb) <= clash_dist ) {
					bb = dnr_coord_list.erase(bb);
				} else {
					++bb;
				}
			}
		}
	}
	utility::io::ozstream outPDB_stream;
	outPDB_stream.open("NPHR.pdb", std::ios::out);

	for (auto & aa : acp_coord_list) {
		outPDB_stream
			<<std::setw(6)<<"ATOM  "
			<<std::setw(5)<<" 1  "
			<<std::setw(5)<<"  O  "
			//<<std::setw(4)<<chemical::aa_from_oneletter_code(chemical::oneletter_code_from_aa(curr_rsd.aa()))
			<<std::setw(4)<<"ACP"
			<<" "
			<<std::setw(1)<<"A"
			<<std::setw(4)<<" 1  "
			<<"    "
			<<std::setw(8)<<std::fixed<<std::setprecision(3)<<aa.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<aa.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<aa.z()<<std::endl;
	}
	for (auto & bb : dnr_coord_list) {
		outPDB_stream
			<<std::setw(6)<<"ATOM  "
			<<std::setw(5)<<" 1  "
			<<std::setw(5)<<"  N  "
			//<<std::setw(4)<<chemical::aa_from_oneletter_code(chemical::oneletter_code_from_aa(curr_rsd.aa()))
			<<std::setw(4)<<"DNR"
			<<" "
			<<std::setw(1)<<"B"
			<<std::setw(4)<<" 2  "
			<<"    "
			<<std::setw(8)<<std::fixed<<std::setprecision(3)<<bb.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<bb.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<bb.z()<<std::endl;
	}

	outPDB_stream.close();

	return;
}

std::string  GenPharmacophore::extract_rna_rings_from_protein_rna_complex(core::pose::Pose const & protein_pose, core::pose::Pose const & rna_pose){

	std::ostringstream outRINGS_sstream;
	outRINGS_sstream.str("");

	//we do not want include water molecules for calculating rna base sasa
	//so strip any atoms that are not protein atoms
	core::pose::Pose old_pose = protein_pose;
	core::pose::Pose new_pose;
	new_pose.clear();
	for ( int ir = 1, ir_end = protein_pose.total_residue(); ir <= ir_end; ir++ ) {
		core::conformation::Residue const & curr_rsd = protein_pose.residue(ir);
		if ( !curr_rsd.is_protein() ) continue;
		new_pose.append_residue_by_jump(protein_pose.residue(ir), new_pose.total_residue(),"", "",  false);
	}

	using namespace basic::options;
	int const sasa_cutoff = option[ OptionKeys::gen_pharmacophore::ring_sasa_cutoff ];

	for ( int ir = 1, ir_end = rna_pose.total_residue(); ir <= ir_end; ir++ ) {

		core::conformation::Residue const & curr_rna_rsd = rna_pose.residue(ir);
		Size seq_pos = curr_rna_rsd.seqpos();
		char rna_chain_id = rna_pose.pdb_info()->chain(seq_pos);

		if ( !curr_rna_rsd.is_RNA() ) continue;
		core::pose::Pose temp_protein_rnabase_pose = new_pose;
		temp_protein_rnabase_pose.append_residue_by_jump(rna_pose.residue(ir), new_pose.total_residue(),"", "",  true);

		//Set-up atomID for SASA calculations by atom
		utility::vector1< core::Real > complex_rsd_sasa( temp_protein_rnabase_pose.total_residue(), 0.0 );
		core::id::AtomID_Map<core::Real> complex_atom_sasa;
		core::id::AtomID_Map<bool> complex_atom_subset;
		core::pose::initialize_atomid_map( complex_atom_sasa, temp_protein_rnabase_pose, 0.0 );
		core::pose::initialize_atomid_map( complex_atom_subset, temp_protein_rnabase_pose, true );

		//Calculate per-atom sasa
		core::Real probe_radius = basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius];
		core::scoring::calc_per_atom_sasa( temp_protein_rnabase_pose, complex_atom_sasa, complex_rsd_sasa, probe_radius, false, complex_atom_subset );

		for ( int ic = 1, ic_end = temp_protein_rnabase_pose.total_residue(); ic<=ic_end; ++ic ) {
			core::conformation::Residue const & curr_rsd = temp_protein_rnabase_pose.residue(ic);
			if ( !curr_rsd.is_RNA() ) continue;
			core::chemical::rna::RNA_ResidueType const & curr_rsd_type = curr_rsd.RNA_type();
			core::Real curr_ring_sasa = 0;
			for ( Size jc = 1, jc_end = curr_rsd.nheavyatoms(); jc <= jc_end; ++jc ) {
				if ( !curr_rsd_type.is_RNA_base_atom(jc) ) continue;
				id::AtomID const aid( jc, ic);
				curr_ring_sasa += complex_atom_sasa[aid];
			}

			// Filter RINGS based on base-sasa-cutoff
			if ( is_buried_ring(curr_rsd, curr_ring_sasa, sasa_cutoff) ) {
				for ( Size jc = 1, jc_end = curr_rsd.nheavyatoms(); jc <= jc_end; ++jc ) {
					if ( !curr_rsd_type.is_RNA_base_atom(jc) ) continue;//include only rings from rna base

					std::string str1 = curr_rsd.atom_name( jc );
					str1.erase( remove( str1.begin(), str1.end(), ' ' ), str1.end() ); //remove any blank space

					if ( curr_rsd.is_purine() ) {
						if ( str1.compare("N2") != 0 ) {
							if ( str1.compare("O6") != 0 ) {
								if ( str1.compare("N6") != 0 ) {
									outRINGS_sstream
										<<std::setw(6)<<"ATOM  "
										<<std::setw(5)<<jc
										<<std::setw(5)<< curr_rsd.atom_name( jc )
										<<std::setw(4)<<"RNG"
										<<" "
										<<std::setw(1)<<rna_chain_id
										<<std::setw(4)<<seq_pos
										<<"    "
										<<std::setw(8)<<std::fixed<<std::setprecision(3)<<curr_rsd.atom(jc).xyz()(1)<<std::setw(8)<<std::fixed<<std::setprecision(3)<<curr_rsd.atom(jc).xyz()(2)<<std::setw(8)<<std::fixed<<std::setprecision(3)<<curr_rsd.atom(jc).xyz()(3)<<std::endl;
								}
							}
						}
					} else if ( !curr_rsd.is_purine() ) {
						if ( str1.compare("O2") != 0 ) {
							if ( str1.compare("O4") != 0 ) {
								if ( str1.compare("N4") != 0 ) {
									outRINGS_sstream
										<<std::setw(6)<<"ATOM  "
										<<std::setw(5)<<jc
										<<std::setw(5)<< curr_rsd.atom_name( jc )
										<<std::setw(4)<<"RNG"
										<<" "
										<<std::setw(1)<<rna_chain_id
										<<std::setw(4)<<seq_pos
										<<"    "
										<<std::setw(8)<<std::fixed<<std::setprecision(3)<<curr_rsd.atom(jc).xyz()(1)<<std::setw(8)<<std::fixed<<std::setprecision(3)<<curr_rsd.atom(jc).xyz()(2)<<std::setw(8)<<std::fixed<<std::setprecision(3)<<curr_rsd.atom(jc).xyz()(3)<<std::endl;
								}
							}
						}
					}
				}
			}
		}
	}
	return outRINGS_sstream.str();

}


std::string GenPharmacophore::extract_Hbond_atoms_from_protein_rna_complex(core::pose::Pose const & protein_pose, core::pose::Pose const & rna_pose){

	std::ostringstream outHBOND_sstream;
	outHBOND_sstream.str("");

	//combine protein and rna pose
	pose::Pose prot_rna_complex_pose = protein_pose;
	for ( int ir = 1, ir_end = rna_pose.total_residue(); ir <= ir_end; ir++ ) {
		core::conformation::Residue const & curr_rna_rsd = rna_pose.residue(ir);
		if ( !curr_rna_rsd.is_RNA() ) continue;
		//int seqpos = curr_rna_rsd.seqpos();
		prot_rna_complex_pose.append_residue_by_jump(rna_pose.residue(ir), prot_rna_complex_pose.total_residue(),"", "",  true);
	}

	//find Hbond interactions and include it to the pharmacophore list
	core::scoring::hbonds::HBondSet rna_hb_set;
	scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	(*scorefxn)(prot_rna_complex_pose);

	assert( prot_rna_complex_pose.energies().residue_neighbors_updated() );
	core::scoring::EnergyGraph const & energy_graph( prot_rna_complex_pose.energies().energy_graph() );
	core::scoring::TenANeighborGraph const & tenA_neighbor_graph( prot_rna_complex_pose.energies().tenA_neighbor_graph() );
	core::scoring::hbonds::HBondDatabaseCOP hb_database = core::scoring::hbonds::HBondDatabase::get_database(core::scoring::hbonds::HBondOptions().params_database_tag());

	rna_hb_set.clear();

	for ( Size ic = 1, ic_end = prot_rna_complex_pose.total_residue(); ic<=ic_end; ic++ ) {
		core::conformation::Residue const & curr_rsd = prot_rna_complex_pose.residue(ic);
		if ( !curr_rsd.is_RNA() ) continue;

		int const nnrna = tenA_neighbor_graph.get_node( ic )->num_neighbors_counting_self_static();

		for ( core::graph::Graph::EdgeListConstIter
				nit = energy_graph.get_node(ic)->const_edge_list_begin(),
				nite = energy_graph.get_node(ic)->const_edge_list_end();
				nit != nite; ++nit ) {

			Size const pin( (*nit)->get_other_ind(ic) );

			core::conformation::Residue const& rn( prot_rna_complex_pose.residue( pin ) );
			int const nnn = tenA_neighbor_graph.get_node( pin )->num_neighbors_counting_self_static();

			if ( !rn.is_RNA() ) {

				identify_hbonds_1way( *hb_database, curr_rsd, rn, nnrna, nnn, false,
					false, false, false, false, rna_hb_set);

				identify_hbonds_1way( *hb_database, rn, curr_rsd, nnn, nnrna, false,
					false, false, false, false, rna_hb_set);
			}
		}
	}
	//TR << "\n\nHBOND SET:\n";
	for ( Size i=1; i<= Size(rna_hb_set.nhbonds()); ++i ) {

		core::scoring::hbonds::HBond const & hb( rna_hb_set.hbond(i) );
		core::conformation::Residue const & donor = prot_rna_complex_pose.residue( hb.don_res() );
		core::conformation::Residue const & accep = prot_rna_complex_pose.residue( hb.acc_res() );
		Size const donor_hatm_num = hb.don_hatm();
		Size const donor_base_atom_num = donor.atom_base( donor_hatm_num );
		Size const accep_atom_num = hb.acc_atm();

		/*TR << i << ":" <<
		chemical::oneletter_code_from_aa(donor.aa()) <<
		temp_protein_rnabase_pose.pdb_info()->number(donor.seqpos()) << temp_protein_rnabase_pose.pdb_info()->chain(donor.seqpos()) << ' ' <<
		'(' << donor.seqpos() << ')' <<
		donor.atom_name( donor_base_atom_num ) << " --- " <<
		chemical::oneletter_code_from_aa(accep.aa()) <<
		temp_protein_rnabase_pose.pdb_info()->number(accep.seqpos()) << temp_protein_rnabase_pose.pdb_info()->chain(accep.seqpos()) << ' ' <<
		'(' << accep.seqpos() << ')' <<
		accep.atom_name( accep_atom_num ) << "\n";
		*/
		//print PDB
		if ( donor.is_RNA() ) {
			//<<std::setw(4)<<chemical::aa_from_oneletter_code(chemical::oneletter_code_from_aa(donor.aa()))
			outHBOND_sstream
				<<std::setw(6)<<"ATOM  "
				<<std::setw(5)<<i
				<<std::setw(5)<< donor.atom_name( donor_base_atom_num )
				<<std::setw(4)<<"DNR"
				<<" "
				<<std::setw(1)<<prot_rna_complex_pose.pdb_info()->chain(donor.seqpos())
				<<std::setw(4)<<donor.seqpos()
				<<"    "
				<<std::setw(8)<<std::fixed<<std::setprecision(3)<<donor.atom(donor_base_atom_num).xyz()(1)<<std::setw(8)<<std::fixed<<std::setprecision(3)<<donor.atom(donor_base_atom_num).xyz()(2)<<std::setw(8)<<std::fixed<<std::setprecision(3)<<donor.atom(donor_base_atom_num).xyz()(3)<<std::endl;
		}

		if ( accep.is_RNA() ) {
			//<<std::setw(4)<<chemical::aa_from_oneletter_code(chemical::oneletter_code_from_aa(accep.aa()))
			outHBOND_sstream
				<<std::setw(6)<<"ATOM  "
				<<std::setw(5)<<i
				<<std::setw(5)<< accep.atom_name( accep_atom_num )
				<<std::setw(4)<<"ACP"
				<<" "
				<<std::setw(1)<<prot_rna_complex_pose.pdb_info()->chain(accep.seqpos())
				<<std::setw(4)<<accep.seqpos()
				<<"    "
				<<std::setw(8)<<std::fixed<<std::setprecision(3)<<accep.atom(accep_atom_num).xyz()(1)<<std::setw(8)<<std::fixed<<std::setprecision(3)<<accep.atom(accep_atom_num).xyz()(2)<<std::setw(8)<<std::fixed<<std::setprecision(3)<<accep.atom(accep_atom_num).xyz()(3)<<std::endl;
		}
	}

	return outHBOND_sstream.str();
}

std::string  GenPharmacophore::make_compatible_with_ROCS_custom_ForceField(std::string & input_filename){

	std::istringstream inPDB_sstream(input_filename);
	std::ostringstream outROCS_sstream;
	outROCS_sstream.str("");
	std::string line;
	while ( getline(inPDB_sstream,line) ) {
		std::string str1 = line.substr (17,3);
		str1.erase( remove( str1.begin(), str1.end(), ' ' ), str1.end() ); //remove any blank space
		if ( str1.compare("ACP") == 0 ) {
			line.replace(12,4," Ne ");
		} else if ( str1.compare("DNR") == 0 ) {
			line.replace(12,4," Be ");
		} else if ( str1.compare("RNG") == 0 ) {
			line.replace(11,5, ReplaceString(line.substr (11,5), "N", "C"));
			//line = ReplaceString(line, "N", "C");
			//std::string ATOMNAME = line.substr (11,5);
			//std::replace (line.begin(), line.end(), 'N',  'C');
			//std::cout<<"ATOMNAME "<<ATOMNAME<<std::endl;
			//line.boost::replace_all(ATOMNAME, "O", "C");
			//boost::replace_all(line.substr(11,5), "O", "C");
			//boost::replace_all_copy(line, line.substr(11,5), "N", "C");
		}
		outROCS_sstream<<line<<std::endl;
	}

	return outROCS_sstream.str();

}

std::string GenPharmacophore::ReplaceString(std::string subject, const std::string& search, const std::string& replace) {
	size_t pos = 0;
	while ( (pos = subject.find(search, pos)) != std::string::npos ) {
		subject.replace(pos, search.length(), replace);
		pos += replace.length();
	}
	return subject;
}

void GenPharmacophore::print_string_to_PDBfile(std::string const & input_filename, std::string const & output_filename) const {

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(output_filename, std::ios::out);
	std::istringstream inPDB_sstream(input_filename);
	std::string line;
	while ( getline(inPDB_sstream,line) ) {
		outPDB_stream<< line <<std::endl;
	}
	return;
}

void GenPharmacophore::cluster_KeyFeatures(std::string const & input_filename, std::string const & output_filename) const {

	using namespace basic::options;
	Size const num_ring = option[ OptionKeys::gen_pharmacophore::min_num_ring ];
	core::Real const ring_ring_dist = option[ OptionKeys::gen_pharmacophore::ring_ring_dist_cutoff ];
	core::Real const ring_atm_dist = option[ OptionKeys::gen_pharmacophore::ring_atm_dist_cutoff ];
	std::istringstream inPDB_sstream(input_filename);

	if ( !inPDB_sstream.good() ) {
		cout << "Cannot open file" << endl;
		exit(-1);
	}

	// string phr_name(input_filename);
	// int pdb_delim = phr_name.find('_');
	// string pdb_code = phr_name.substr(0, pdb_delim);

	typedef pair<string, string> mID;
	map<mID, SmallMol> moieties;
	vector<SmallMol> rings;
	vector<SmallMol> dnrAcp;
	vector<UnionEdge> edges;

	string line;

	while ( !inPDB_sstream.eof() ) {
		getline(inPDB_sstream, line);
		boost::trim(line);
		if ( line.find("ATOM") == string::npos ) {
			continue;
		}

		string chainID = line.substr(21, 1);
		string resSeq  = line.substr(22, 4);
		string mType   = line.substr(17, 3);
		mID newAtm(chainID, resSeq);

		if ( mType == "ACP" || mType == "DNR" ) {
			SmallMol newMol = SmallMol();
			newMol.set_name(chainID + " " + resSeq + " " + mType);
			newMol.add_atom(line);
			dnrAcp.push_back(newMol);
		} else if ( mType == "RNG" ) {
			moieties[newAtm].set_name(chainID + " " + resSeq + " " + mType);
			moieties[newAtm].add_atom(line);
		}
	}

	map<mID, SmallMol>::iterator it, end ;
	for ( it = moieties.begin(), end = moieties.end(); it != end; ++it ) {
		rings.push_back(it->second);
	}

	cout << "RNG molecules: " << rings.size() << "\t" << "DNR-ACP atoms: " << dnrAcp.size() << endl;

	// build ring - Dnr_Acp pairs
	// double max_dist = 5.0;
	vector<SmallMol>::iterator it1, it2, end1, end2;
	for ( it1 = dnrAcp.begin(), end1 = dnrAcp.end(); it1 != end1; ++it1 ) {
		for ( it2 = rings.begin(), end2 = rings.end(); it2 != end2; ++it2 ) {
			double dist = it1->cal_min_dist(&(*it2));
			if ( dist <= ring_atm_dist * ring_atm_dist ) {
				UnionEdge newEdge(&(*it1), &(*it2));
				edges.push_back(newEdge);
			}
		}
	}

	sort(edges.begin(), edges.end());

	vector<UnionEdge>::size_type l;
	for ( l = 0; l < edges.size(); l++ ) {
		SmallMol *a = edges[l].get_a();
		SmallMol *b = edges[l].get_b();

		if ( a->get_visited() == true ) continue;
		if ( !a->connected(b) ) {
			a->set_visited(true);
			a->connect(b);
			// int old_size = b->findRoot()->get_size();
			// b->findRoot()->set_size(old_size - 1);
		}
	}
	edges.clear();

	// max_dist = 5.0;
	for ( it1 = rings.begin(); it1 != end2; ++it1 ) {
		for ( it2 = it1 + 1; it2 != end2; ++it2 ) {
			double dist = it1->cal_min_dist(&(*it2));
			if ( dist <= ring_ring_dist * ring_ring_dist ) {
				UnionEdge newEdge(&(*it1), &(*it2));
				edges.push_back(newEdge);
			}
		}
	}

	sort(edges.begin(), edges.end());

	// connect only rings
	for ( l = 0; l < edges.size(); l++ ) {
		SmallMol *a = edges[l].get_a();
		SmallMol *b = edges[l].get_b();

		if ( !a->connected(b) ) {
			a->connect(b);
		}
	}


	// use UnionFind information to print out clusters
	map<SmallMol *, vector<SmallMol *> > clusters;
	for ( it1 = rings.begin(); it1 != end2; ++it1 ) {
		SmallMol *r = it1->get_root();
		vector<SmallMol *> &c = clusters[r];
		c.push_back(&(*it1));
	}
	for ( it2 = dnrAcp.begin(); it2 != end1; ++it2 ) {
		SmallMol *r = it2->get_root();
		vector<SmallMol *> &c = clusters[r];
		c.push_back(&(*it2));
	}

	map<SmallMol *, vector<SmallMol *> >::iterator it3, end3;
	int counter = 0;
	for ( it3 = clusters.begin(), end3 = clusters.end(); it3 != end3; ++it3 ) {
		vector<SmallMol *> &current = it3->second;
		Size numRng = 0;
		vector<SmallMol *>::size_type s;
		for ( s = 0; s < current.size(); s++ ) {
			if ( current[s]->numberOfAtoms() > 1 ) {
				numRng++;
			}
		}

		if ( numRng >= num_ring ) {
			ofstream output((output_filename + "_" + boost::lexical_cast<string>(counter++) + ".pdb").c_str());
			for ( s = 0; s < current.size(); s++ ) {
				output << current[s]->getContent();
			}
			output.close();
		}
	}

	moieties.clear();
	clusters.clear();
	edges.clear();
	rings.clear();
	dnrAcp.clear();

	return;
}

} // Pockets
} // protocols
