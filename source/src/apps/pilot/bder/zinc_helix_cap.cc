// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/bder/hydrophobic_sasa.cc
/// @brief Prints a list of buried unsatisfied polar atoms with a pymol selection for visualization.  This has code duplication but it is for convenient use.
/// @details
/// @author Bryan Der

#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/MetricValue.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <utility/file/FileName.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <set>


//tracers
using basic::Error;
using basic::Warning;
static basic::Tracer TR("apps.pilot.bder.zinc_helix_cap");


using namespace core;

namespace local {
basic::options::RealOptionKey const hsasa_residue_cutoff("hsasa_residue_cutoff");
}//local

/// @brief
class zinc_helix_cap : public protocols::moves::Mover {
public:
	zinc_helix_cap()
	{
	}
	virtual ~zinc_helix_cap(){};



	virtual
	void
	apply( core::pose::Pose & pose ){

		utility::file::FileName input_pdbname( pose.pdb_info()->name() );


		using namespace core::scoring;
		core::scoring::ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( PRE_TALARIS_2013_STANDARD_WTS, SCORE12_PATCH );
		scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn->energy_method_options() );
		energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
		scorefxn->set_energy_method_options( energymethodoptions );
		scorefxn->score( pose );


		using namespace core::pose::metrics;

		TR << "Registering SASA Calculator" << std::endl;
		if ( !CalculatorFactory::Instance().check_calculator_exists( "sasa_calc_name" ) ) {
			CalculatorFactory::Instance().register_calculator( "sasa_calc_name", new pose::metrics::simple_calculators::SasaCalculator() );
		}
		TR << "Registering num_Hbond Calculator" << std::endl;
		if ( !CalculatorFactory::Instance().check_calculator_exists( "num_hbonds_calc_name" ) ) {
			CalculatorFactory::Instance().register_calculator( "num_hbonds_calc_name", new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator() );
		}
		TR << "Registering BurUnsatPolar Calculator" << std::endl;
		if ( !CalculatorFactory::Instance().check_calculator_exists( "bur_unsat_calc_name" ) ) {
			CalculatorFactory::Instance().register_calculator( "bur_unsat_calc_name", new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator( "sasa_calc_name", "num_hbonds_calc_name", 0.01 /*sasa=0 is buried*/) );
		}

		basic::MetricValue<id::AtomID_Map< bool > > bur_unsat_atomid_map;
		pose.metric( "bur_unsat_calc_name", "atom_bur_unsat", bur_unsat_atomid_map );

		basic::MetricValue< core::id::AtomID_Map< core::Size > > atom_hbonds;
		pose.metric( "num_hbonds_calc_name", "atom_Hbonds", atom_hbonds );

		basic::MetricValue< core::id::AtomID_Map< core::Real > > atom_sasa;
		pose.metric( "sasa_calc_name", "atom_sasa", atom_sasa);






		std::set<Size> loop_positions;

		TR << "Now via DSSP" << std::endl;
		TR << "select resi ";
		for ( Size res = 5; res <= pose.size() - 5; res++ ) {


			core::scoring::dssp::Dssp dssp( pose );
			char dssp_i = dssp.get_dssp_secstruct(res);
			if ( dssp_i != 'H' ) {
				//TR << res << "+";
				loop_positions.insert(res - 3);
				loop_positions.insert(res - 2);
				loop_positions.insert(res - 1);
				loop_positions.insert(res);
				loop_positions.insert(res + 1);
				loop_positions.insert(res + 2);
				loop_positions.insert(res + 3);
			}

		}

		utility::vector1< core::id::AtomID > buried_unsat_atom_ids;


		for ( std::set< Size >::const_iterator it(loop_positions.begin()), end(loop_positions.end()); it!=end; ++it ) {
			//TR << *it << "+";

			//check for buried unsat of NH and O

			conformation::Residue const & rsd = pose.residue( *it );

			for ( Size at = 1; at <= 4; ++at ) { //mainchain atoms are 1 through 4
				core::id::AtomID atid( at, *it );


				if ( bur_unsat_atomid_map.value()[atid] /*bool*/ ) {

					buried_unsat_atom_ids.push_back( atid );
					TR << pose.pdb_info()->name() << " buried_unsat " << atid << "  SASA: " << atom_sasa.value()[atid] << "  HBOND: " << atom_hbonds.value()[atid] << std::endl;
				}

			}


		}
		TR << std::endl;
		TR << "Total number of buried unsat: " << pose.pdb_info()->name() << "  " << buried_unsat_atom_ids.size() << std::endl;
		core::Real norm_total = ((core::Real)buried_unsat_atom_ids.size())/loop_positions.size();
		TR << "Total number of normalized buried unsat: " << pose.pdb_info()->name() << "  " << norm_total << std::endl;

		std::string pdbname_base = input_pdbname.base();
		std::string pymol_bunsat_selection = "bunsat_" + pdbname_base;
		TR << "PYMOL_SELECTION: select " << pymol_bunsat_selection << ",";

		for ( Size i(1); i <= buried_unsat_atom_ids.size(); ++i ) {
			Size rosetta_resnum = buried_unsat_atom_ids[i].rsd();
			std::string atom_name = pose.residue(rosetta_resnum).atom_name( buried_unsat_atom_ids[i].atomno() );

			Size pdb_resnum = pose.pdb_info()->number( rosetta_resnum );
			char chain = pose.pdb_info()->chain( rosetta_resnum );

			if ( i != buried_unsat_atom_ids.size() ) {
				TR << " resi " << pdb_resnum << " and name " << atom_name << " and chain " << chain << " in " << pdbname_base << " +";
			} else {
				TR << " resi " << pdb_resnum << " and name " << atom_name << " and chain " << chain << " in " << pdbname_base;
			}


		}
		TR << std::endl;
		TR << "PYMOL_SELECTION: alter " << pymol_bunsat_selection << ", vdw=0.5" << std::endl;



		return;
	}


	virtual
	std::string
	get_name() const { return "zinc_helix_cap"; }



private:


};

typedef utility::pointer::owning_ptr< zinc_helix_cap > zinc_helix_capOP;

int main( int argc, char* argv[] )
{
	using basic::options::option;
	option.add( local::hsasa_residue_cutoff, "hydrophobic sasa per residue").def(50);


	devel::init(argc, argv);
	protocols::jd2::JobDistributor::get_instance()->go(new zinc_helix_cap);

	TR << "Recommended option:   -sasa_calculator_probe_radius 1.2" << std::endl;

	TR << "************************d**o**n**e**************************************" << std::endl;

	return 0;
}



// //compute set of hbonds
// core::scoring::hbonds::HBondSet hbond_set;
// hbond_set.setup_for_residue_pair_energies( pose, false, false );



// TR << "select resi ";


// for (Size res = 1; res <= pose.size(); res++) {
//  bool is_helix( false );

//   for (Size i = 1; i<= hbond_set.nhbonds(); i++) {
//    core::scoring::hbonds::HBondCOP hbond(hbond_set.hbond(i));

//    if(hbond->energy() > -0.5) { // arbitrary value for hbond strength cutoff
//     continue;
//    }

//    if( hbond->don_res() == res + 4 && hbond->don_hatm_is_backbone() ) {
//     //this is a helical residue
//     //TR << "res " << res << " don_res: " << hbond->don_res() << std::endl;
//     is_helix = true;
//     break;
//    }

//    else if( hbond->acc_res() == res - 4 && hbond->acc_atm_is_protein_backbone() ) {
//     //this is a helical residue
//     //TR << "res " << res << " acc_res: " << hbond->acc_res() << std::endl;
//     is_helix = true;
//     break;
//    }

//   }

//   if( is_helix ) {
//    //TR << "Residue " << res << " is helical" << std::endl;
//   }
//   else {
//    TR << res << "+";
//   }

// }
// TR << std::endl;
