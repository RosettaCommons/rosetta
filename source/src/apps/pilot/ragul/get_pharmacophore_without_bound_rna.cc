// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author Ragul Gowthaman

// Project Headers
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
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
#include <core/chemical/AtomType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/rna/RNA_Info.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>

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
#include <core/io/pdb/pdb_writer.hh>
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

using namespace core;
using namespace core::pose::datacache;
using namespace core::optimization;
using namespace core::pose::metrics;
using namespace core::scoring;
using namespace core::scoring::constraints;
using namespace core::id;
using namespace core::chemical;
using namespace core::chemical::rna;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::conformation;


OPT_KEY( String, input_protein )
OPT_KEY( Integer, rna_base_sasa_cutoff )
OPT_KEY( Real, clash_dist_cutoff )


static THREAD_LOCAL basic::Tracer TR( "apps.pilot.ragul.rna_phr.get_pharmacophore_without_bound_rna" );

// Returns TRUE if SASA for respective base is <= histogram cutoff
// cutoff value choosen based on the PDB scan of protein-RNA complxes
bool
is_buried_ring(core::conformation::Residue const & rsd, core::Real ring_sasa, core::Real sasa_cutoff){
	if ( rsd.name3() == "  A" ) return ring_sasa <= sasa_cutoff;//46.81 ;
	else if ( rsd.name3() == "  C" )  return ring_sasa <= sasa_cutoff;//31.09 ;
	else if ( rsd.name3() == "  G" )  return ring_sasa <= sasa_cutoff;//45.06 ;
	else if ( rsd.name3() == "  U" )  return ring_sasa <= sasa_cutoff;//52.66 ;
	else return false;
}

numeric::xyzVector<core::Real> rotatePoint(core::Real x, core::Real y, core::Real z){
	numeric::xyzVector<core::Real> coord(x,y,z);
	numeric::xyzMatrix<core::Real>  rot_mat = numeric::x_rotation_matrix_degrees((core::Real)0);
	coord = rot_mat * coord;
	return coord;
}

core::Real get_RNAring_sasa( core::conformation::Residue const & rsd, int rsdno,
	core::id::AtomID_Map<core::Real> & pose_atom_sasa ){
	if ( !rsd.is_RNA() ) {
		std::cout<<"Residue is not an RNA base. Cannot calculate RNA base SASA."<<std::endl;
		exit (1);
	}
	core::chemical::rna::RNA_Info const & rsd_type = rsd.RNA_info();
	core::Real curr_ring_sasa = 0;
	for ( Size rj = 1, rj_end = rsd.nheavyatoms(); rj <= rj_end; ++rj ) {
		if ( !rsd_type.is_RNA_base_atom(rj) ) continue;
		id::AtomID const aid( rj, rsdno);
		curr_ring_sasa += pose_atom_sasa[aid];
	}
	return(curr_ring_sasa);
}

int main( int argc, char * argv [] ){

	try{

		NEW_OPT( input_protein, "rna protein name", "protein.pdb" );
		NEW_OPT( rna_base_sasa_cutoff, "rna_base_sasa_cutoff", 25 );
		NEW_OPT(  clash_dist_cutoff, " clash_dist_cutoff", 1.0 );

		devel::init(argc, argv);

		std::string const input_protein_pose = option[ input_protein ];
		//int const sasa_cutoff = option[ rna_base_sasa_cutoff ];
		core::Real const clash_dist = option[  clash_dist_cutoff ];


		pose::Pose protein_pose;
		core::import_pose::pose_from_file( protein_pose, input_protein_pose , core::import_pose::PDB_file);

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

		for ( Size j = 1, resnum = protein_pose.size(); j <= resnum; ++j ) {
			core::conformation::Residue const & rsd( protein_pose.conformation().residue(j) );

			/*int offset = 9;
			if (rsd.seqpos() < 72) offset=3;
			offset=0;*/
			//int target=0;
			/*core::Size total_atoms(0);
			using namespace basic::options;
			if (option[ OptionKeys::fingerprint::include_hydrogens ]()){
			total_atoms = rsd.natoms();
			} else {
			total_atoms = rsd.nheavyatoms();
			}*/

			//fill in points that are ideal for a hydrogen acceptor with an O
			for ( core::chemical::AtomIndices::const_iterator hnum  = rsd.Hpos_polar().begin(), hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
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
					numeric::xyzVector<core::Real> rrpoint = rotatePoint(ro1.x(),ro1.y(),ro1.z());
					std::cout<<"Optimal acceptor at "<<ro1.x()<<" "<<ro1.y()<<" "<<ro1.z()<<std::endl;
					acp_coord_list.push_back(rrpoint);
				}
			}


			//fill in points that are ideal for a hydrogen donor with an N
			for ( core::chemical::AtomIndices::const_iterator anum  = rsd.accpt_pos().begin(), anume = rsd.accpt_pos().end(); anum != anume; ++anum ) {
				Size const aatm( *anum );
				// Skip buried residues
				// int offset = 9;
				// if (rsd.seqpos() < 72) offset=3;
				// offset=0;
				if ( atom_sasas(j, aatm) < 0.1 ) {
					//std::cout<<rsd.seqpos()+offset<<" Acceptor "<<rsd.name()<<" "<<rsd.atom_name(aatm)<<" SASA "<<atom_sasas(j, aatm)<<" being ignored"<<std::endl;
					continue;
				}
				//std::cout<<rsd.seqpos()+offset<<" Acceptor "<<rsd.name()<<" "<<rsd.atom_name(aatm)<<" SASA "<<atom_sasas(j, aatm)<<std::endl;

				numeric::xyzVector<core::Real> const & aatm_xyz( rsd.xyz( aatm ) );
				numeric::xyzVector<core::Real> aatm_base_xyz( rsd.xyz( rsd.atom_base( aatm ) ) );
				numeric::xyzVector<core::Real> const & aatm_base2_xyz( rsd.xyz( rsd.abase2( aatm ) ) );
				Hybridization const & hybrid( rsd.atom_type(aatm).hybridization() );

				core::Real theta(0.0);// step_size(0.0);
				utility::vector1< core::Real > phi_list, phi_steps;
				phi_steps.push_back(  0 );
				switch( hybrid ) {
				case SP2_HYBRID :
					theta = 180.0 - 120.0;
					//step_size = 15.0;
					phi_list.push_back(   0.0 );
					phi_list.push_back( 180.0 );
					break;
				case SP3_HYBRID :
					theta = 180.0 - 109.0;
					//step_size = 10.0;
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
					//step_size = 0.0; // doesnt matter
					break;
				}
				default :
					std::cerr << "Bad hybridization type for acceptor " << hybrid << '\n';
					exit(1000);
				}
				Stub stub( aatm_xyz, aatm_base_xyz, aatm_base2_xyz );
				for ( Size i=1; i<= phi_list.size(); ++i ) {
					numeric::xyzVector<core::Real> const ro1(stub.spherical( numeric::conversions::radians( phi_list[i]), numeric::conversions::radians( theta ), opt_distance));
					numeric::xyzVector<core::Real> rrpoint = rotatePoint(ro1.x(),ro1.y(),ro1.z());
					std::cout<<"Optimal donor at "<<ro1.x()<<" "<<ro1.y()<<" "<<ro1.z()<<std::endl;
					dnr_coord_list.push_back(rrpoint);
				}
			}
		}

		numeric::xyzVector<core::Real> protein_atom_coord(0.);
		for ( int j = 1, resnum = protein_pose.size(); j <= resnum; ++j ) {
			core::conformation::Residue const & curr_rsd = protein_pose.residue(j);
			for ( Size i = 1, i_end = curr_rsd.natoms(); i <= i_end; ++i ) {
				protein_atom_coord.x() = curr_rsd.atom(i).xyz()(1);
				protein_atom_coord.y() = curr_rsd.atom(i).xyz()(2);
				protein_atom_coord.z() = curr_rsd.atom(i).xyz()(3);
				for ( std::list< numeric::xyzVector<core::Real> >::iterator aa = acp_coord_list.begin(); aa != acp_coord_list.end(); /* ++aa */ ) {
					if ( protein_atom_coord.distance(*aa) <= clash_dist ) {
						aa = acp_coord_list.erase(aa);
					} else {
						++aa;
					}
				}
				for ( std::list< numeric::xyzVector<core::Real> >::iterator bb = dnr_coord_list.begin(); bb != dnr_coord_list.end(); /* ++bb */ ) {
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

		for ( std::list< numeric::xyzVector<core::Real> >::iterator aa = acp_coord_list.begin(); aa != acp_coord_list.end(); ++aa ) {
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
				<<std::setw(8)<<std::fixed<<std::setprecision(3)<<aa->x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<aa->y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<aa->z()<<std::endl;
		}

		for ( std::list< numeric::xyzVector<core::Real> >::iterator bb = dnr_coord_list.begin(); bb != dnr_coord_list.end(); ++bb ) {
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
				<<std::setw(8)<<std::fixed<<std::setprecision(3)<<bb->x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<bb->y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<bb->z()<<std::endl;
		}

		outPDB_stream.close();

	}
catch ( utility::excn::EXCN_Base const & e ) {
	std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
}

	return 0;

}

