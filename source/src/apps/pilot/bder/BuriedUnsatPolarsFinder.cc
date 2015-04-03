// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/bder/BuriedUnsatPolarsFinder.cc
/// @brief
/// @details
/// @author Bryan Der

#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
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


//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static thread_local basic::Tracer TR( "apps.pilot.bder.BuriedUnsatPolarsFinder" );


using namespace core;


/// @brief
class BuriedUnsatPolarsFinder : public protocols::moves::Mover {
public:
  BuriedUnsatPolarsFinder()
  {
  }
  virtual ~BuriedUnsatPolarsFinder(){};


  virtual
  void
  apply( core::pose::Pose & pose ){

		utility::file::FileName filename( pose.pdb_info()->name() );
		std::string pdbname_base = filename.base();
		std::string pymol_bunsat_selection = "bunsat_" + pdbname_base;

		using namespace core::scoring;
		core::scoring::ScoreFunctionOP scorefxn = get_score_function();
		scorefxn->score( pose );


		using namespace core::pose::metrics;

		TR << "Registering SASA Calculator" << std::endl;
		if( !CalculatorFactory::Instance().check_calculator_exists( "sasa_calc_name" ) ){
			CalculatorFactory::Instance().register_calculator( "sasa_calc_name", new pose::metrics::simple_calculators::SasaCalculatorLegacy() );
		}
		TR << "Registering num_Hbond Calculator" << std::endl;
		if( !CalculatorFactory::Instance().check_calculator_exists( "num_hbonds_calc_name" ) ){
			CalculatorFactory::Instance().register_calculator( "num_hbonds_calc_name", new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator() );
		}
		TR << "Registering BurUnsatPolar Calculator" << std::endl;
		if( !CalculatorFactory::Instance().check_calculator_exists( "bur_unsat_calc_name" ) ){
			CalculatorFactory::Instance().register_calculator( "bur_unsat_calc_name", new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator( "sasa_calc_name", "num_hbonds_calc_name", 0.01 /*sasa=0 is buried*/) );
		}

		basic::MetricValue<id::AtomID_Map< bool > > bur_unsat_atomid_map;
		pose.metric( "bur_unsat_calc_name", "atom_bur_unsat", bur_unsat_atomid_map );

		std::string num_unsathbonds_string = pose.print_metric("bur_unsat_calc_name", "all_bur_unsat_polars");
		TR << pose.pdb_info()->name() << " has " << num_unsathbonds_string << " unsat buried polar atoms." << std::endl;


		utility::vector1< core::id::AtomID > buried_unsat_atom_ids;

		for( Size i = 1; i <= pose.total_residue(); ++i){
			conformation::Residue const & rsd = pose.residue( i );
			for( Size at = 1; at <= rsd.nheavyatoms(); ++at){
				core::id::AtomID atid( at, i );
				if( bur_unsat_atomid_map.value()[atid] /*bool*/) {
					buried_unsat_atom_ids.push_back( atid );
				}
			}
		}
		Size total_number_buried_unsats = buried_unsat_atom_ids.size();


		TR << "PYMOL_SELECTION: select " << pymol_bunsat_selection << ",";

		for(Size i(1); i <= buried_unsat_atom_ids.size(); ++i) {
			Size rosetta_resnum = buried_unsat_atom_ids[i].rsd();
			std::string atom_name = pose.residue(rosetta_resnum).atom_name( buried_unsat_atom_ids[i].atomno() );

			Size pdb_resnum = pose.pdb_info()->number( rosetta_resnum );
			char chain = pose.pdb_info()->chain( rosetta_resnum );

			if( i != buried_unsat_atom_ids.size() ) {
				TR << " resi " << pdb_resnum << " and name " << atom_name << " in " << pdbname_base << " and chain " << chain << " +";
			}
			else {
				TR << " resi " << pdb_resnum << " and name " << atom_name << " in " << pdbname_base << " and chain " << chain;
			}


		}
		TR << std::endl;
		TR << "PYMOL_SELECTION: alter " << pymol_bunsat_selection << ", vdw=0.5" << std::endl;


    return;
  }


	virtual
	std::string
	get_name() const { return "BuriedUnsatPolarsFinder"; }


private:

	// basic::MetricValue< id::AtomID_Map< Real > > atom_sasa_;
	// basic::MetricValue< id::AtomID_Map< Size > > atom_hbonds_;

};

typedef utility::pointer::owning_ptr< BuriedUnsatPolarsFinder > BuriedUnsatPolarsFinderOP;

int main( int argc, char* argv[] )
{
	try {
	using basic::options::option;

  devel::init(argc, argv);
  protocols::jd2::JobDistributor::get_instance()->go(new BuriedUnsatPolarsFinder);

	TR << "Recommended option:   -sasa_calculator_probe_radius 1.2" << std::endl;

  TR << "************************d**o**n**e**************************************" << std::endl;

	} catch (utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

  return 0;
}

