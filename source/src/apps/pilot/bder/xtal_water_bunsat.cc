// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/bder/xtal_water_bunsat.cc
/// @brief
/// @details
/// @author Bryan Der

#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculator.hh>
#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>
#include <devel/vardist_solaccess/LoadVarSolDistSasaCalculatorMover.hh>
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
#include <numeric/xyzVector.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>

typedef numeric::xyzVector<core::Real> point;

//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static basic::Tracer TR("apps.pilot.bder.xtal_water_bunsat");


using namespace core;

//basic::options::FileOptionKey const dry_protein( "dry_protein" );
basic::options::BooleanOptionKey const use_varsoldist_sasa_calc( "use_varsoldist_sasa_calc" );


///@brief
class xtal_water_bunsat : public protocols::moves::Mover {
public:
	xtal_water_bunsat()
	{
	}
	virtual ~xtal_water_bunsat(){};



  	virtual
  	void
  	apply( core::pose::Pose & pose ){
		
		pose::Pose dry_pose( pose );

		for(Size i(1); i<=dry_pose.total_residue(); ++i) {
			if(! dry_pose.residue(i).is_protein() ) {
				dry_pose.conformation().delete_residue_slow(i);
				--i;
			}
		}
		dry_pose_ = dry_pose;
		
		dry_pose_.dump_pdb("dry_pose.pdb");
		
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
		scorefxn->score(dry_pose_);
		
		using namespace core::pose::metrics;
		
		if(basic::options::option[use_varsoldist_sasa_calc]) {
			TR << "Registering Andrew's SASA Calculator" << std::endl;
			if( !CalculatorFactory::Instance().check_calculator_exists( "sasa_calc_name" ) ){
				CalculatorFactory::Instance().register_calculator( "sasa_calc_name", new devel::vardist_solaccess::VarSolDistSasaCalculator() );
			}
		}
		else {
			TR << "Registering SASA Calculator" << std::endl;
			if( !CalculatorFactory::Instance().check_calculator_exists( "sasa_calc_name" ) ){
				CalculatorFactory::Instance().register_calculator( "sasa_calc_name", new pose::metrics::simple_calculators::SasaCalculator() );
			}
		}
		

		
		
		dry_pose_.metric( "sasa_calc_name", "atom_sasa", atom_sasa_);
		
		
		
		
		
	  
		water_points_.clear();
		for(Size i(1); i<=pose.total_residue(); ++i) {
			if(pose.residue(i).name() == "TP3") {
				water_points_.push_back(pose.residue(i).atom(1).xyz());
			}
		}
		TR << "waters: " << water_points_.size() << std::endl;

		//check NH
		for(Size i(2); i<= dry_pose_.total_residue(); ++i) {
			if(dry_pose_.residue(i).is_protein() && dry_pose_.residue(i).name3() != "PRO") {
				id::AtomID this_atomid = id::AtomID( dry_pose_.residue(i).atom_index("H"), i);
				if(atomid_is_contacting_crystallographic_water(this_atomid)) {
					Real this_sasa(check_sasa( this_atomid ));
					TR << " SASA: " << this_sasa << "   A2" << std::endl;
				}
			}
		}
		
		//check O
		for(Size i(1); i<= dry_pose_.total_residue() - 1; ++i) {
			if(dry_pose_.residue(i).is_protein()) {
				id::AtomID this_atomid = id::AtomID( 4, i );
				if(atomid_is_contacting_crystallographic_water(this_atomid)){
					Real this_sasa(check_sasa( this_atomid ));
					TR << " SASA: " << this_sasa << "   A2" << std::endl;
				}
			}
		}
		
		



		return;
	}
	
	
	virtual
	bool
	atomid_is_contacting_crystallographic_water( id::AtomID atid ) {
		Real min_water_dist(9999.9);
		Real distance_cutoff(3.0);
		if(dry_pose_.residue(atid.rsd()).atom_name(atid.atomno()) == " O  ") {
			distance_cutoff = 3.1;
		}
		else if(dry_pose_.residue(atid.rsd()).atom_name(atid.atomno()) == " H  ") {
			distance_cutoff = 2.3;
		}
		
		
		point atid_point = dry_pose_.residue(atid.rsd()).atom(atid.atomno()).xyz();
		
		
		for(Size i(1); i<=water_points_.size(); ++i) {
			Real water_dist = atid_point.distance(water_points_[i]);
			if(water_dist < min_water_dist) {
				min_water_dist = water_dist;
			}
		}
		
		if(min_water_dist < distance_cutoff) {
			TR << dry_pose_.pdb_info()->name() << " " << atid << " contacting water " << dry_pose_.residue(atid.rsd()).atom_name(atid.atomno()) << " " << min_water_dist;
			return true;
		}
		else {
			return false;
		}
	
	}
	
	virtual
	Real
	check_sasa( id::AtomID atid ) {
		
		Real this_sasa = atom_sasa_.value()[atid];
		
		
		if(dry_pose_.residue(atid.rsd()).atom_name(atid.atomno()) == " H  ") {
			id::AtomID bbN_atid( 1, atid.rsd() );
			this_sasa += atom_sasa_.value()[bbN_atid];
		}
		
		return this_sasa;
		
		
	}
	


	virtual
	std::string
	get_name() const { return "xtal_water_bunsat"; }



private:

	basic::MetricValue< core::id::AtomID_Map< core::Real > > atom_sasa_;

	utility::vector1<point> water_points_;
	
	core::pose::Pose dry_pose_;
	
	// basic::MetricValue< id::AtomID_Map< Real > > atom_sasa_;
	// basic::MetricValue< id::AtomID_Map< Size > > atom_hbonds_;

};

typedef utility::pointer::owning_ptr< xtal_water_bunsat > xtal_water_bunsatOP;

int main( int argc, char* argv[] )
{
	try {
	using basic::options::option;
	option.add( use_varsoldist_sasa_calc, "var sol d sasa calculator" ).def(true);


  devel::init(argc, argv);
  protocols::jd2::JobDistributor::get_instance()->go(new xtal_water_bunsat);

	TR << "Recommended option:   -sasa_calculator_probe_radius 1.2" << std::endl;

  TR << "************************d**o**n**e**************************************" << std::endl;

	} catch (utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

  return 0;
}

