// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/kevin/BuriedUnsatPolarsFinder3.cc
/// @brief
/// @details
/// @author Kevin Houlihan
/// @author Bryan Der

#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <devel/buns/BuriedUnsatisfiedPolarsCalculator2.hh>
#include <basic/options/keys/bunsat_calc2.OptionKeys.gen.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>

#include <core/pose/metrics/simple_calculators/SasaCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>
//#include <core/scoring/hbonds/hbonds.hh>
//#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/MetricValue.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>
//#include <devel/vardist_solaccess/LoadVarSolDistSasaCalculatorMover.hh>
#include <numeric/xyzVector.hh>

#include <core/scoring/sasa.hh>
#include <numeric/conversions.hh> //degrees-radians

//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static thread_local basic::Tracer TR( "apps.pilot.kevin.BuriedUnsatPolarsFinder3" );


using namespace basic::options;
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace utility;
using namespace core::pose::metrics;
using namespace protocols::toolbox::pose_metric_calculators;
using namespace basic::options;

typedef numeric::xyzVector< core::Real > Vector;

/*
basic::options::BooleanOptionKey const use_varsoldist_sasa_calc("use_varsoldist_sasa_calc");
basic::options::BooleanOptionKey const use_regular_sasa_calc("use_regular_sasa_calc");
//basic::options::BooleanOptionKey const add_probe_radius("add_probe_radius");

basic::options::BooleanOptionKey const generous_hbond_settings("generous_hbond_settings");
basic::options::RealOptionKey const sasa_cutoff("bunsat_sasa_cutoff");
basic::options::RealOptionKey const AHD_cutoff("bunsat_AHD_cutoff");
basic::options::RealOptionKey const dist_cutoff("bunsat_dist_cutoff");
basic::options::RealOptionKey const hydroxyl_dist_cutoff("bunsat_hydroxyl_dist_cutoff");
basic::options::RealOptionKey const sulphur_dist_cutoff("bunsat_sulphur_dist_cutoff");
basic::options::RealOptionKey const metal_dist_cutoff("bunsat_metal_dist_cutoff");
 */

///@brief
class BuriedUnsatPolarsFinder : public protocols::moves::Mover {
public:
	BuriedUnsatPolarsFinder() {}
	
	virtual ~BuriedUnsatPolarsFinder(){};
	
	virtual
	void
	apply(pose::Pose & pose)
	{
		
		TR << "PDB: " << pose.pdb_info()->name() << std::endl;
		
		/*
		assert(option[use_varsoldist_sasa_calc] ||
			   option[use_regular_sasa_calc]);
		 */
		
		utility::file::FileName filename(pose.pdb_info()->name());
		std::string pdbname_base = filename.base();
		//std::string pymol_bunsat_selection  = "bunsat_"  + pdbname_base;
		//std::string pymol_bunsat_selection2 = "bunsat2_" + pdbname_base;
		
		using namespace core::scoring;
		
		core::scoring::ScoreFunctionOP scorefxn = get_score_function();
		scorefxn->score(pose);

		// TODO continue from here
		// make calculators for buried unsat,
		// make buried unsat for buried unsat 2
		
		/*
		TR << "Registering num_Hbond Calculator" << std::endl;
		std::string num_hbonds_calc_name("num_hbonds_calc_name_");
		CalculatorFactory::Instance().remove_calculator(num_hbonds_calc_name);
		CalculatorFactory::Instance().register_calculator(num_hbonds_calc_name,
				new NumberHBondsCalculator());
		
		std::string sasa_calc_name("sasa_calc_name_");
		std::string bunsat_calc_name("bunsat_calc_name");
		Real sasa_buried = option[OptionKeys::bunsat_calc2::sasa_burial_cutoff]; //sasa_cutoff
		
		CalculatorFactory::Instance().remove_calculator(sasa_calc_name);
		CalculatorFactory::Instance().remove_calculator(bunsat_calc_name);
		
		if(basic::options::option[OptionKeys::bunsat_calc2::layered_sasa]) {
			TR << "Registering VarSolDist SASA Calculator" << std::endl;
			CalculatorFactory::Instance().register_calculator(sasa_calc_name,
					new devel::vardist_solaccess::VarSolDistSasaCalculator());
		}
		else {
			TR << "Registering SASA Calculator" << std::endl;
			CalculatorFactory::Instance().register_calculator(sasa_calc_name,
					new pose::metrics::simple_calculators::SasaCalculator());
		}
				
		CalculatorFactory::Instance().register_calculator(bunsat_calc_name,
				new BuriedUnsatisfiedPolarsCalculator(
				sasa_calc_name, num_hbonds_calc_name, sasa_buried));
		
		//basic::MetricValue< id::AtomID_Map< Real > > sasa_atomid_map;
		basic::MetricValue< id::AtomID_Map< bool > > bunsat_atomid_map;
		
		//pose.metric(sasa_calc_name, "atom_sasa", sasa_atomid_map);

		 */

		std::string bunsat_calc_name("default");
		basic::MetricValue< id::AtomID_Map< bool > > bunsat_atomid_map;
	
		std::string bunsat_calc2_name("bunsat_calc2_name");
		if( !CalculatorFactory::Instance().check_calculator_exists( bunsat_calc2_name ) ){
			CalculatorFactory::Instance().register_calculator(bunsat_calc2_name,
				new devel::buns::BuriedUnsatisfiedPolarsCalculator2( bunsat_calc_name ));
		}
		pose.metric(bunsat_calc2_name, "atom_bur_unsat", bunsat_atomid_map);
		
		
		// TODO Remove all this, for debugging
		//bunsat_atomid_map.print();

		/*
		bool nobunsat = true;
		for( Size i = 1; i <= pose.total_residue(); ++i){ //convert atomid-map to vector of atom ids
			conformation::Residue const & rsd = pose.residue( i );
			for( Size at = 1; at <= rsd.nheavyatoms(); ++at){
				core::id::AtomID atid( at, i );
				if( bunsat_atomid_map.value()[atid] ) {
					nobunsat = false;
					TR << "bunsat at residue " << i << ", atom " << at << std::endl;
				}
			}
		}
		if (nobunsat) TR << "no bunsat" << std::endl;
		 */
		
		
		return;
	}
		
	virtual
	std::string
	get_name() const { return "BuriedUnsatPolarsFinder"; }
	
	
	
private:
	
	// basic::MetricValue< id::AtomID_Map< Real > > atom_sasa_;
	// basic::MetricValue< id::AtomID_Map< Size > > atom_hbonds_;
	
};

typedef utility::pointer::owning_ptr<BuriedUnsatPolarsFinder> BuriedUnsatPolarsFinderOP;

int main(int argc, char* argv[])
{
	try {
		using basic::options::option;
		/*
		option.add(use_varsoldist_sasa_calc, "var sol dist sasa calculator").def(true);
		option.add(use_regular_sasa_calc, "regular sasa calculator").def(false);
		//option.add(add_probe_radius, "include_probe_radius_in_atom_radii").def(false);
		
		option.add(generous_hbond_settings, "decompose, allow sc/bb in secstruct, no env_dep").def(true);
		
		option.add(sasa_cutoff, "SASA cutoff to be considered buried").def(0.01);
		option.add(AHD_cutoff, "AHD angle above which hbonds are considered "
				   "for geometry-only based H-bond checking").def(120);
		option.add(dist_cutoff, "acceptor-hydrogen distance below which hbonds "
				   "are considered for geometry-only based H-bond checking").def(3.0);
		option.add(hydroxyl_dist_cutoff, "acceptor to hydroxyl oxygen ditance cutoff"
				   "for detecting h-bonds by heavyatom positions; meant to "
				   "account for ambiguity in hydrogen placement").def(3.5);
		option.add(sulphur_dist_cutoff, "sulphur to hydrogen distance cutoff below "
				   "which an h-bond is detected").def(3.3);
		option.add(metal_dist_cutoff, "metal atom to hydrogen distance cutoff below "
				   "which an h-bond is detected").def(2.7);
		 */
		
		
		devel::init(argc, argv);
		protocols::jd2::JobDistributor::get_instance()->go(new BuriedUnsatPolarsFinder);
		
		TR << "************************d**o**n**e**************************************" << std::endl;
		
	} catch (utility::excn::EXCN_Base const & e) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	
	return 0;
}
