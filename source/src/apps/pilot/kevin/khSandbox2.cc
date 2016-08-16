// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps.pilot.kevin.khSandbox2.cc
/// @brief
/// @details
/// @author Kevin Houlihan

#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculator.hh>
#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>
#include <devel/vardist_solaccess/LoadVarSolDistSasaCalculatorMover.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
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

//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static THREAD_LOCAL basic::Tracer TR( "apps.pilot.kevin.khSandbox2" );


using namespace core;

//basic::options::BooleanOptionKey const use_varsoldist_sasa_calc( "use_varsoldist_sasa_calc" );

/// @brief
class khSandbox2 : public protocols::moves::Mover {
public:
	khSandbox2()
	{
	}
	virtual ~khSandbox2(){};


	virtual
	void
	apply( core::pose::Pose & pose ){
		return;
	}

	virtual
	std::string
	get_name() const { return "khSandbox2"; }

};

typedef utility::pointer::owning_ptr< khSandbox2 > khSandbox2OP;

int main( int argc, char* argv[] )
{
	try {
//	using basic::options::option;
//	option.add( use_varsoldist_sasa_calc, "var sol d sasa calculator" ).def(true);

	devel::init(argc, argv);
	protocols::jd2::JobDistributor::get_instance()->go(new khSandbox2);

	TR << "Recommended option:   -sasa_calculator_probe_radius 1.2" << std::endl;

	TR << "************************d**o**n**e**************************************" << std::endl;

	} catch (utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}

