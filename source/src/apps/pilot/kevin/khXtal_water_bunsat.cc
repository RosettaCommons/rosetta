// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps.pilot.kevin.khxtal_water_bunsat.cc
/// @brief
/// @details
/// @author Bryan Der
/// @author Kevin Houlihan

#include <iomanip>

#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
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
static THREAD_LOCAL basic::Tracer TR( "apps.pilot.kevin.khxtal_water_bunsat" );


using namespace core;

//basic::options::FileOptionKey const dry_protein( "dry_protein" );
basic::options::BooleanOptionKey const use_varsoldist_sasa_calc( "use_varsoldist_sasa_calc" );
basic::options::BooleanOptionKey const water_dist_H( "water_dist_H_cutoff" );
basic::options::BooleanOptionKey const water_dist_O( "water_dist_O_cutoff" );

/// @brief
class xtal_water_bunsat : public protocols::moves::Mover {
public:
	xtal_water_bunsat() {}
	virtual ~xtal_water_bunsat(){};

	virtual
	void
	apply( core::pose::Pose & pose ){

		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		scorefxn->score(pose);

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
				CalculatorFactory::Instance().register_calculator( "sasa_calc_name", new pose::metrics::simple_calculators::SasaCalculatorLegacy() );
			}
		}

		basic::MetricValue< core::id::AtomID_Map< core::Real > > atom_sasa;
		pose.metric( "sasa_calc_name", "atom_sasa", atom_sasa);

		utility::vector1<point> water_points;
		water_points.clear();
		//for(Size i(1); i<=wet_pose.size(); ++i) {
		//	if(wet_pose.residue(i).name() == "TP3") {
		//		water_points.push_back(wet_pose.residue(i).atom(1).xyz());
		//	}
		//}

        //utility::vector1< numeric::xyzVector<core::Real> > water_points;

        //for (auto w : pose.pdb_info()->get_unrecognized_atoms()) {
		const utility::vector1< core::pose::UnrecognizedAtomRecord >& ua = pose.pdb_info()->get_unrecognized_atoms();
        for (Size i=1, end=ua.size(); i<=end; i++) {
                const std::string& res_name = ua[i].res_name();
                const std::string& atom_name = ua[i].atom_name();
                if ( (res_name  == "HOH" ||
                      res_name  == "DOD" ||
                      res_name  == "SOL" ||
                      res_name  == "H2O" ||
                      res_name  == "D2O"   ) &&
                      atom_name == "O"       &&
                      ua[i].temp() <= 30
                ) {
                    water_points.push_back(ua[i].coords());
                }
        }
		TR << "waters: " << water_points.size() << std::endl;

		core::conformation::Conformation const dry_conf = pose.conformation();
		core::Size num_chains = dry_conf.num_chains();
		utility::vector1< core::Size > nter_resnums;
		utility::vector1< core::Size > cter_resnums;
		nter_resnums.resize(num_chains);
		cter_resnums.resize(num_chains);
		for (Size i=1; i<=num_chains; i++) {
			nter_resnums[i] = dry_conf.chain_begin(i);
			cter_resnums[i] = dry_conf.chain_end(i);
		}

		core::Real water_H_dist_cutoff = basic::options::option[water_dist_H]();
		water_H_dist_cutoff = 2.3;
		TR << "water_H_dist_cutoff: " << water_H_dist_cutoff << std::endl;
		//check NH
		for(Size i=1; i<=pose.size(); ++i) {
			if (!pose.residue(i).is_protein()) continue;
			if (nter_resnums.has_value(i)) {
				report_sasa_if_contacting_water(pose, i, "1H  ", atom_sasa,
						water_points, water_H_dist_cutoff);
				report_sasa_if_contacting_water(pose, i, "2H  ", atom_sasa,
						water_points, water_H_dist_cutoff);
				if(pose.residue(i).name3() != "PRO")
					report_sasa_if_contacting_water(pose, i, "3H  ", atom_sasa,
							water_points, water_H_dist_cutoff);
			}
			else if(pose.residue(i).name3() != "PRO") {
				report_sasa_if_contacting_water(pose, i, " H  ", atom_sasa,
						water_points, water_H_dist_cutoff);
			}
		}

		core::Real water_O_dist_cutoff = basic::options::option[water_dist_O];
		water_O_dist_cutoff = 3.1;
		TR << "water_O_dist_cutoff: " << water_O_dist_cutoff << std::endl;
		//check O
		for(Size i=1; i<=pose.size(); ++i) {
			if (!pose.residue(i).is_protein()) continue;
			report_sasa_if_contacting_water(pose, i, " O  ", atom_sasa,
					water_points, water_O_dist_cutoff);
			if (cter_resnums.has_value(i)) {
				report_sasa_if_contacting_water(pose, i, " OXT", atom_sasa,
						water_points, water_O_dist_cutoff);
			}
		}

		return;
	}

	virtual
	core::Real
	closest_crystallographic_water_dist(
		core::pose::Pose const & dry_pose,
		id::AtomID const atid,
		utility::vector1<point> const & water_points
	) {
		Real min_water_dist(9999.9);
		point atid_point = dry_pose.residue(atid.rsd()).atom(atid.atomno()).xyz();

		for(Size i(1); i<=water_points.size(); ++i) {
			Real water_dist = atid_point.distance(water_points[i]);
			if(water_dist < min_water_dist) {
				min_water_dist = water_dist;
			}
		}
		return min_water_dist;

	}

	virtual
	Real
	check_sasa(
		core::pose::Pose const &,
		id::AtomID const atid,
		basic::MetricValue< core::id::AtomID_Map< core::Real > > const & atom_sasa
	) {

		Real this_sasa = atom_sasa.value()[atid];
		return this_sasa;

	}
	
	virtual
	std::string
	get_name() const { return "xtal_water_bunsat"; }


private:
	virtual
	void
	report_sasa_if_contacting_water(
		core::pose::Pose const & dry_pose,
		core::Size const resi,
		std::string const & atom_name,
		basic::MetricValue< core::id::AtomID_Map< core::Real > > const & atom_sasa,
		utility::vector1<point> const & water_points,
		core::Real const max_dist_contact
	) {
		id::AtomID atid = id::AtomID(dry_pose.residue(resi).atom_index(atom_name), resi);
		core::Real closest_xtal_water_dist =
			closest_crystallographic_water_dist(dry_pose, atid, water_points);

		if(closest_xtal_water_dist <= max_dist_contact) {
			Real this_sasa = check_sasa(dry_pose, atid, atom_sasa);

			std::ios  state(NULL);
			state.copyfmt(TR);
	
			TR << std::setprecision(3) << std::fixed;
			TR << dry_pose.pdb_info()->name() << " res " << std::setw(6)
				<< dry_pose.pdb_info()->pose2pdb(resi) << ", atom \""
				<< dry_pose.residue(atid.rsd()).atom_name(atid.atomno())
				<< "\" contacting water at dist " << std::setw(7) << closest_xtal_water_dist
				<< " A, with SASA " << std::setw(7) << this_sasa << " A2" << std::endl;
	
			TR.copyfmt(state);
		}
		else return;

	}

};

typedef utility::pointer::owning_ptr< xtal_water_bunsat > xtal_water_bunsatOP;

int main( int argc, char* argv[] ) {
	try {
	using basic::options::option;
	option.add( use_varsoldist_sasa_calc, "var sol d sasa calculator" ).def(true);
	option.add( water_dist_H, "water_dist_H" ).def(core::Real(2.3));
	option.add( water_dist_O, "water_dist_O" ).def(core::Real(3.1));

	devel::init(argc, argv);
	protocols::jd2::JobDistributor::get_instance()->go(new xtal_water_bunsat);

	TR << "Recommended option:   -sasa_calculator_probe_radius 1.2" << std::endl;

	TR << "************************d**o**n**e**************************************" << std::endl;

	} catch (utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}

