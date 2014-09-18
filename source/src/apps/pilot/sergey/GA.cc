// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test file for GA minimizer algorithm
///
/// @brief
/// @author Sergey Lyskov

//#include "GA_Minimizer.cc"  // debug in Pilot Apps
#include <core/optimization/GA_Minimizer.hh>


#include <core/optimization/Multifunc.hh>
#include <core/optimization/types.hh>

#include <numeric/random/random.hh>

#include <basic/Tracer.hh>

#include <core/types.hh>
#include <devel/init.hh>



#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


#include <cmath>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <numeric/NumericTraits.hh>

using core::Real;
using core::Size;

static thread_local basic::Tracer TR( "GA" );


class CheeseFunction : public core::optimization::Multifunc
{
public:
	CheeseFunction() {}

	Real operator ()( core::optimization::Multivec const & phipsi ) const {
		// "Fi(x) = 10 + abs(x)/10 - 10.*cos(2.*PI*x)";
		const utility::vector1< Real > &v(phipsi);

		Real x0 = 20.;
		Real res = 0.;
		const Real PI = numeric::NumericTraits<Real>::pi();

		for(Size i=1; i<=v.size(); i++) {
			Real x = v[i] - x0;
			res += 10+x*x - 10.*cos(2.*PI*x);
			//res += 10 + abs(x)/10 - 10.*cos(2.*PI*x);

			//double a = 10. + fabs(x)/10. - 10.*cos(2.*PI*x);  res += a*a;
			//double a = 10. + fabs(x)/10. - 10.*cos(2.*PI*x);  res += a;
		}
		return res;
    }

	void dfunc( core::optimization::Multivec const &, core::optimization::Multivec &) const {}
};

void test_GA(void)
{
	core::optimization::Multivec v;
	for(int i=0; i<200; i++) v.push_back(1000.1);

	CheeseFunction chf;
	core::optimization::MinimizerOptions minopt("GA", .1, false);
	minopt.minimize_tolerance(1.),
	minopt.ga_mutation_probability(.5);

	core::optimization::GA_Minimizer gam(chf, minopt);
	Real res = gam.run(v, 1000);

	TR << "Res: " << res << " " << v << std::endl;
}

void real_test_GA()
{
	using namespace core;

	pose::PoseOP start_pose;
	kinematics::MoveMap mm;
	core::scoring::ScoreFunctionOP scorefxn;
	core::optimization::AtomTreeMinimizer minimizer;

	start_pose = new pose::Pose();
	core::import_pose::pose_from_pdb(*start_pose, "test_in.pdb");

	scorefxn = core::scoring::get_score_function_legacy( scoring::PRE_TALARIS_2013_STANDARD_WTS );

	//kinematics::MoveMap mm;
	for ( int i=30; i<= 35; ++i ) {
		mm.set_bb ( i, true );
		mm.set_chi( i, true );
	}

	Real en = (*scorefxn)( *start_pose ); // to triger dunbrack loading/calcualtion
	TR << "Start minimization... En=" << en << " -------------------------------" << std::endl;

	std::string stype = "unknow";
	stype = "dfpmin";
	//stype = "GA";
	//stype = "dfpmin_armijo";
	//stype = "dfpmin_armijo_nonmonotone";

	core::optimization::MinimizerOptions options( stype/*"dfpmin"*/, 1e-4, true, true );
	options.minimize_tolerance(.001),
	options.ga_mutation_probability(.5);
	options.max_iter(1000);

	core::pose::Pose pose;
	pose = *start_pose;


	minimizer.run( pose, mm, *scorefxn, options );
	en = (*scorefxn)( pose );
	TR << "End minimization... En=" << en << " -------------------------------" << std::endl;
}

int main( int argc, char * argv [] )
{

	try {

	devel::init(argc, argv);

	//test_GA();
	real_test_GA();

	TR << "GA test ended. --------------------------------" << std::endl;
	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

