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
/// @author Kevin Houlihan

//#include <iomanip>
#include <algorithm>
#include <cctype>
#include <limits>
#include <sstream>

#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>

//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static THREAD_LOCAL basic::Tracer TR( "apps.pilot.kevin.wat_distances" );

using namespace core;
using namespace core::pose::metrics;

basic::options::BooleanOptionKey const use_varsoldist_sasa_calc( "use_varsoldist_sasa_calc" );
basic::options::BooleanOptionKey const water_dist_H( "water_dist_H_cutoff" );
basic::options::BooleanOptionKey const water_dist_O( "water_dist_O_cutoff" );

/// @brief
class xtal_water_bunsat : public protocols::moves::Mover {
public:
	xtal_water_bunsat() {}
	virtual ~xtal_water_bunsat() {};

	virtual
	std::string
	get_name() const { return "xtal_water_bunsat"; }

	virtual
	void
	apply( core::pose::Pose & pose ){

		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		scorefxn->score(pose);

		utility::vector1< numeric::xyzVector<core::Real> > water_coords;

		const utility::vector1< pose::UnrecognizedAtomRecord >& unrec_atoms = pose.pdb_info()->get_unrecognized_atoms();

		for (Size i = 1, end = unrec_atoms.size(); i <= end; i++) {
			const pose::UnrecognizedAtomRecord& unrec_atom = unrec_atoms[i];
			const std::string& res_name = unrec_atom.res_name();
			const std::string& atom_name = unrec_atom.atom_name();
			if ( (res_name  == "HOH" ||
				  res_name  == "DOD" ||
				  res_name  == "SOL" ||
				  res_name  == "H2O" ||
				  res_name  == "D2O"   ) &&
				  atom_name == "O"       &&
				  unrec_atom.temp() <= 30
			) {

			}
		}

		utility::vector1<Real> probe_radii;
		for (Real probe_radius = 0.4; probe_radius <= 1.6; probe_radius += 0.1) {
			probe_radii.push_back(probe_radius);
		}

		utility::vector1<std::string> sasa_calc_names;

		utility::vector1< id::AtomID_Map< core::Real > > calc_atom_sasas;

		for (Size i = 1; i <= probe_radii.size(); i++) {
			std::stringstream sstream;
			sstream << "sasa_pr" << probe_radii[i];
			std::string calc_name(sstream.str());

			if (!CalculatorFactory::Instance().check_calculator_exists( calc_name ) ) {
				CalculatorFactory::Instance().register_calculator(
					calc_name, new pose::metrics::simple_calculators::SasaCalculatorLegacy(probe_radii[i]) );
			}
			sasa_calc_names.push_back(calc_name);
		}

		//basic::MetricValue< core::id::AtomID_Map< core::Real > > atom_sasa;
		//pose.metric(calc_name, "atom_sasa", atom_sasa);
		//calc_atom_sasas.push_back(atom_sasa.value());

		//std::cout << "WRITING SASA: " << calc_name << std::endl;

		std::string calc_name("vsasa");
		if( !CalculatorFactory::Instance().check_calculator_exists(calc_name) ){
			CalculatorFactory::Instance().register_calculator(calc_name, new devel::vardist_solaccess::VarSolDistSasaCalculator() );
		}
		//basic::MetricValue< core::id::AtomID_Map< core::Real > > atom_sasa;

		//pose.metric(calc_name, "atom_sasa", atom_sasa);
		//sasa_calc_names.push_back(calc_name);
		//calc_atom_sasas.push_back(atom_sasa.value());

		//std::cout << "WRITING SASA: " << calc_name << std::endl;

		//std::cout << "WTF????" << std::endl;

		const core::pose::PDBInfo& pdb_info = *(pose.pdb_info());
		for (Size resNum=1; resNum<=pose.total_residue(); ++resNum) {
			const core::conformation::Residue& res = pose.residue(resNum);
			const core::chemical::ResidueType& res_type = res.type();
			for (Size atomNum=1; atomNum<=res.natoms(); ++atomNum) {
				if (pdb_info.temperature(resNum, atomNum) > 30) continue;
				if (pdb_info.occupancy(resNum, atomNum) < 1) continue;
				//const core::conformation::Atom& atom = res.atom(atomNum);
				const core::chemical::AtomType& atom_type = res_type.atom_type(atomNum);
				utility::file::FileName filename(pdb_info.modeltag());
				platform::Real mindist = std::numeric_limits<platform::Real>::max();

				id::AtomID at(atomNum, resNum);
				for (size_t w=0, n=water_coords.size(); w!=n; w++) {
					platform::Real dist = res.atom(atomNum).xyz().distance(water_coords[w]);
					if (dist < mindist) mindist = dist;
				}
				if (mindist > 5.0) continue;
				const char& raw_icode = pdb_info.icode(resNum);
				//char icode = (isalnum(raw_icode) != 0) ? raw_icode : ' ';
				std::cout << filename.base() << ","
					<< pdb_info.chain(resNum)                             << ","
					<< pdb_info.number(resNum)                            << ","
					<< raw_icode                                          << ","
					<< res_type.name3()                                   << ","
					<< atom_type.name()                                   << ","
					<< atom_type.element()                                << ","
					<< mindist                                            ;


				for (Size i = 1; i <= sasa_calc_names.size(); i++) {
					std::cout << "READING SASA: " << sasa_calc_names[i] << std::endl;
					basic::MetricValue< core::id::AtomID_Map< core::Real > > atom_sasa;
					pose.metric(sasa_calc_names[i], "atom_sasa", atom_sasa);
					id::AtomID at(atomNum, resNum);
					std::cout << ", " << sasa_calc_names[i] << ": " << calc_atom_sasas[i][at];
				}
				std::cout << "\n";

					//<< (atom_type.is_polar_hydrogen() ? "true" : "false") << ","
					//<< (atom_type.is_acceptor() ? "true" : "false"      ) << ","
					//<< (atom_type.is_aromatic() ? "true" : "false"      ) << std::endl;
			}
		}

		std::cout.flush();

		return;
	}


private:

};

typedef utility::pointer::owning_ptr< xtal_water_bunsat > xtal_water_bunsatOP;

int main( int argc, char* argv[] ) {
	try {
	basic::options::option.add( water_dist_H, "water_dist_H" ).def(core::Real(2.3));
	basic::options::option.add( water_dist_O, "water_dist_O" ).def(core::Real(3.1));

    //if(!basic::options::option[basic::options::OptionKeys::in::remember_unrecognized_res]()
	//		|| !basic::options::option[basic::options::OptionKeys::in::remember_unrecognized_water]()){
	//	        TR.Warning << "Use -in:remember_unrecognized_res and -in:remember_unrecognized_water to locate unrecognized atoms." << std::endl;
	//}

	devel::init(argc, argv);
	protocols::jd2::JobDistributor::get_instance()->go(new xtal_water_bunsat);

	TR << "************************d**o**n**e**************************************" << std::endl;

	} catch (utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}
