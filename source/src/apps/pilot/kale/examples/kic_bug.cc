// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh>
#include <devel/init.hh>
#include <numeric/random/random.hh>
#include <iomanip>

using namespace std;
using namespace core;
using namespace devel;
using namespace protocols::loops::loop_closure::kinematic_closure;

typedef utility::vector1<Real> ParameterList;
typedef utility::vector1<ParameterList> ParameterMatrix;

class PathologicalKinematicPerturber : public KinematicPerturber {

public:

	PathologicalKinematicPerturber() {};
	~PathologicalKinematicPerturber() {};
	string perturber_type() const { return "PathologicalKinematicPerturber"; }

	void perturb_chain(
			pose::Pose const &pose,
			ParameterList &torsion_angles,
			ParameterList &bond_angles,
			ParameterList &bond_lengths) {

		/* Internal Coordinate Data
		 * ------------------------
		 * torsion_angles     29.2461    171.427    176.577     147.79    12.8668    178.858   -80.5171    -143.65    178.923     14.635   -28.4338    178.845    39.2256   -6.25637    181.292    48.1442   -83.4468    179.689   -169.399   -42.8457    177.923    99.1484   -83.8437    177.051    157.905    93.0474    176.571    105.486    150.429    179.204   -127.218   -59.8358     176.83   -13.9848    74.3853    227.887    318.961    64.3249    323.456 
		 * bond_angles        163.122    110.949    111.247    114.693    110.981    114.701    120.399    110.977     113.26    120.354    110.976    115.985     123.07    110.942    113.978    120.437    110.995    113.837    120.808    110.998     116.09    122.525    110.966    114.171     120.44     110.96    114.442    120.881    110.995    115.729    121.958    111.023    114.901    121.872    111.011    121.441    142.219    110.958    126.139 
		 * bond_lengths       1.46079    1.51089    1.38058    1.45975    1.50969    1.35064    1.45991    1.50969    1.32866    1.45993    1.51056    1.33521    1.45958    1.51033    1.32498    1.45978    1.50988    1.33911    1.45966    1.50935     1.3387    1.45966     1.5105    1.34395    1.46054    1.50979     1.3149    1.45971    1.50918    1.33328    1.45981    1.50992    1.30402    1.45931    1.51015    1.26892     1.4601    1.51075    3.07789 
		 * origin               1.628      0.206      0.303 
		 * frame             0.403206 -0.0773553  -0.911834   0.631593  -0.697521    0.33846  -0.662205  -0.712377  -0.232387 
		 * solution        6 of 6
		 */

		torsion_angles[ 1] = 29.2461;
		torsion_angles[ 2] = 171.427;
		torsion_angles[ 3] = 176.577;
												 
		torsion_angles[ 4] = 147.79;
		torsion_angles[ 5] = 12.8668;
		torsion_angles[ 6] = 178.858;
												 
		torsion_angles[ 7] = -80.5171;
		torsion_angles[ 8] = -143.65;
		torsion_angles[ 9] = 178.923;
												 
		torsion_angles[10] = 14.635;
		torsion_angles[11] = -28.4338;
		torsion_angles[12] = 178.845;
												 
		torsion_angles[13] = 39.2256;
		torsion_angles[14] = -6.25637;
		torsion_angles[15] = 181.292;
												 
		torsion_angles[16] = 48.1442;
		torsion_angles[17] = -83.4468;
		torsion_angles[18] = 179.689;
												 
		torsion_angles[19] = -169.399;
		torsion_angles[20] = -42.8457;
		torsion_angles[21] = 177.923;
												 
		torsion_angles[22] = 99.1484;
		torsion_angles[23] = -83.8437;
		torsion_angles[24] = 177.051;
												 
		torsion_angles[25] = 157.905;
		torsion_angles[26] = 93.0474;
		torsion_angles[27] = 176.571;
												 
		torsion_angles[28] = 105.486;
		torsion_angles[29] = 150.429;
		torsion_angles[30] = 179.204;
												 
		torsion_angles[31] = -127.218;
		torsion_angles[32] = -59.8358;
		torsion_angles[33] = 176.83;
												 
		torsion_angles[34] = -13.9848;
		torsion_angles[35] = 74.3853;
		torsion_angles[36] = 227.887;
												 
		torsion_angles[37] = 318.961;
		torsion_angles[38] = 64.3249;
		torsion_angles[39] = 323.456;
	}
};

int main(int argc, char** argv) {
	devel::init(argc, argv);
	numeric::random::RG.set_seed(2);

	pose::Pose pose;
	string input_pdb = "cycle.pdb";

	KinematicMoverOP mover = new KinematicMover();
	KinematicPerturberOP perturber = new PathologicalKinematicPerturber();

	import_pose::pose_from_pdb(pose, input_pdb);
	pose.dump_pdb("problem.pdb");

	mover->set_pivots(2, 7, 12);
	mover->set_perturber(perturber);
	mover->apply(pose);

	pose.dump_pdb("solution.pdb");
}
