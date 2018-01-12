/****************************************************************************************************
	PCA_internal_test.cc

	A test of Rosetta's numeric library and its ability to do principal component analysis in
	more than 3 dimensions.
****************************************************************************************************/

#include <protocols/simple_moves/ScoreMover.hh>
#include <core/scoring/ScoreFunction.hh>
//#include <core/pose/util.hh>
//#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/chemical/ChemicalManager.hh>
//#include <core/import_pose/import_pose.hh>
//#include <protocols/cluster/cluster.hh>
//#include <protocols/loops/Loops.hh>
//#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
//#include <core/import_pose/pose_stream/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <devel/init.hh>
#include <iostream>
#include <string>
#include <deque>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/cluster.OptionKeys.gen.hh>
//#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <utility/vector1.hh>
//#include <core/scoring/rms_util.hh>
//#include <core/scoring/rms_util.tmpl.hh>
//#include <core/scoring/Energies.hh>
//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>
#include <ObjexxFCL/format.hh>
//#include <ObjexxFCL/FArray2D.hh>
//#include <protocols/minimization_packing/RepackSidechainsMover.hh>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <protocols/simple_moves/MutateResidue.hh>
//#include <protocols/relax/FastRelax.hh>
//#include <numeric/model_quality/rms.hh>
//#include <core/scoring/func/HarmonicFunc.hh>
//#include <core/scoring/func/CircularHarmonicFunc.hh>
//#include <core/scoring/constraints/AtomPairConstraint.hh>
//#include <core/scoring/constraints/AngleConstraint.hh>
//#include <core/scoring/constraints/DihedralConstraint.hh>
//#include <core/scoring/constraints/ConstraintSet.hh>
//#include <core/pose/annotated_sequence.hh>

//#include <numeric/EulerAngles.hh>

//For silent file output:
//#include <core/io/silent/SilentStructFactory.hh>
//#include <core/io/silent/SilentFileData.hh>
//#include <core/io/silent/SilentStruct.hh>

//#include <protocols/simple_moves/ConstraintSetMover.hh>

#include <numeric/PCA.hh>

#define PI 3.1415926535897932384626433832795
#define CNCa_ANGLE 121.7
#define CNH_ANGLE 119.15
#define CaCN_ANGLE 116.2
#define OCN_ANGLE 123.01

using ObjexxFCL::format::F;
//using namespace protocols::cluster;

//OPT_KEY( StringVector, v_extra_rms_atoms) //Extra atoms to use in the RMSD calculation.

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//utility::vector1<core::Size> empty_vector;

	//NEW_OPT ( v_extra_rms_atoms, "A list of additional atoms to use in the RMSD calculation, in the format \"residue:atomname residue:atomname residue:atomname\".  For example, \"-v_extra_rms_atoms 7:SG 12:CG 12:CD 12:CE 12:NZ 14:OG\".  Default empty list.", "");
}

using namespace core;
using namespace ObjexxFCL;
//using namespace core::pose;
using namespace protocols;
using namespace basic::options;

//MAIN
int main( int argc, char * argv [] ) {

	//using namespace protocols;
	//using namespace protocols::moves;
	//using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//using namespace utility::file;
	//using namespace protocols::cluster;
	//using namespace basic::options::OptionKeys::cluster;
	using namespace std;

	printf("Starting PCA_internal_test.cc\nFile created 17 June 2014 by Vikram K. Mulligan\n"); fflush(stdout);

	register_options();
	devel::init(argc, argv);

	utility::vector1< utility::vector1< core::Real > > inputmatrix;
	utility::vector1< core::Real > v1(5);
	utility::vector1< core::Real > v2(5);
	utility::vector1< core::Real > v3(5);
	v1[1]=5.3;
	v1[2]=3.3;
	v1[3]=2.1;
	v1[4]=1.2;
	v1[5]=0.5;
	v2[1]=10.1;
	v2[2]=6.1;
	v2[3]=4.4;
	v2[4]=2.0;
	v2[5]=1.1;
	v3[1]=15.8;
	v3[2]=9.3;
	v3[3]=6.3;
	v3[4]=3.9;
	v3[5]=1.6;
	inputmatrix.push_back(v1);
	inputmatrix.push_back(v2);
	inputmatrix.push_back(v3);

	printf("Input matrix:\n");
	for(core::Size j=1; j<=5; ++j){
		for(core::Size i=1; i<=inputmatrix.size(); ++i) {
			printf("%.2f\t", inputmatrix[i][j]);
		}
		printf("\n");
	}
	printf("\n"); fflush(stdout);

	std::pair<utility::vector1< utility::vector1< Real > >, utility::vector1< Real > > result = numeric::principal_components_and_eigenvalues_ndimensions( inputmatrix, true );

	utility::vector1 <utility::vector1 <core::Real> > &vectors (result.first);
	utility::vector1 <core::Real> &variances (result.second);

	printf("Output matrix:\n");
	for(core::Size j=1; j<=5; ++j){
		for(core::Size i=1; i<=5; ++i) {
			printf("%.4f\t", vectors[i][j]);
		}
		printf("\n");
	}
	printf("\n"); fflush(stdout);

	printf("Output variances:\n");
	for(core::Size i=1; i<=5; ++i) printf("%.1f\t", variances[i]);
	printf("\n"); fflush(stdout);
	
	printf("\n*****JOB COMPLETED*****\n"); fflush(stdout);
	return 0;
}

