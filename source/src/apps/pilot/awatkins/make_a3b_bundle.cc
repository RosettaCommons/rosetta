#include <protocols/simple_moves/ScoreMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/util.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/cluster/cluster.hh>
#include <protocols/loops/Loops.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <devel/init.hh>
#include <iostream>
#include <string>
#include <deque>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/cluster.OptionKeys.gen.hh>
//#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <utility/vector1.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <ObjexxFCL/format.hh>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <protocols/relax/FastRelax.hh>
#include <numeric/random/random.hh>
#include <numeric/random/uniform.hh>
#include <numeric/constants.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/TorsionID.hh>
#include <core/scoring/Energies.hh>
#include <protocols/simple_moves/RepackSidechainsMover.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <protocols/simple_moves/MutateResidue.hh>
//#include <core/chemical/PatchOperation.hh>
#include <core/chemical/ResidueType.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicWrapper.hh>

#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

//using ObjexxFCL::format::F;
using namespace protocols::cluster;

//#define PI 3.14159265358979323846264338327950288

OPT_KEY (Real, v_Dz_N)
OPT_KEY (Real, v_Dz_C)
OPT_KEY (Real, v_Dz_CM)
OPT_KEY (Real, v_phi1_N_denom)
OPT_KEY (Real, v_R0)
OPT_KEY (Real, v_R1)

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	utility::vector1<core::Size> empty_vector;
	NEW_OPT ( v_Dz_N, "", -1.088888);
	NEW_OPT ( v_Dz_C, "", 0.8836);
	NEW_OPT ( v_Dz_CM, "", -0.5);
	NEW_OPT ( v_phi1_N_denom, "", 2.9);
	NEW_OPT ( v_R0, "", 6.0);
	NEW_OPT ( v_R1, "", 3.5);
}


using namespace core;
using namespace ObjexxFCL;
using namespace core::pose;
using namespace protocols;
using namespace basic::options;

core::Real Xpos (
	core::Size t,
	core::Real R0,
	core::Real omega0,
	core::Real phi0,
	core::Real R1,
	core::Real omega1,
	core::Real phi1,
	core::Real alpha,
	core::Real Dz
) {
	core::Real phi0prime=phi0+Dz*tan(alpha)/R0;
	core::Real c0 = cos(omega0*(double)t+phi0prime);
	core::Real c1 = cos(omega1*(double)t+phi1);
	core::Real s0 = sin(omega0*(double)t+phi0prime);
	core::Real s1 = sin(omega1*(double)t+phi1);

	return R0*c0+R1*c0*c1-R1*cos(alpha)*s0*s1;
}

core::Real Ypos (
	core::Size t,
	core::Real R0,
	core::Real omega0,
	core::Real phi0,
	core::Real R1,
	core::Real omega1,
	core::Real phi1,
	core::Real alpha,
	core::Real Dz
) {
	core::Real phi0prime=phi0+Dz*tan(alpha)/R0;
	core::Real c0 = cos(omega0*(double)t+phi0prime);
	core::Real c1 = cos(omega1*(double)t+phi1);
	core::Real s0 = sin(omega0*(double)t+phi0prime);
	core::Real s1 = sin(omega1*(double)t+phi1);

	return R0*s0+R1*s0*c1+R1*cos(alpha)*c0*s1;
}

core::Real Zpos (
	core::Size t,
	core::Real R0,
	core::Real omega0,
	core::Real R1,
	core::Real omega1,
	core::Real phi1,
	core::Real alpha,
	core::Real Dz
) {

	return omega0*R0/tan(alpha)*(double)t-R1*sin(alpha)*sin(omega1*(double)t+phi1)+Dz;
}

//Functions to set the backbone dihedral angles of beta-amino acid peptides:
void betapeptide_setphi (
	core::pose::Pose &pose,
	core::Size resnumber,
	core::Real angle
) {
	using namespace core;
	using namespace core::id;
	pose.set_torsion( TorsionID( resnumber, id::BB, 1 ), angle );
	return;
}

void betapeptide_settheta (
	core::pose::Pose &pose,
	core::Size resnumber,
	core::Real angle
) {
	using namespace core;
	using namespace core::id;
	pose.set_torsion( TorsionID( resnumber, id::BB, 2 ), angle );
	return;
}

void betapeptide_setpsi (
	core::pose::Pose &pose,
	core::Size resnumber,
	core::Real angle
) {
	using namespace core;
	using namespace core::id;
	pose.set_torsion( TorsionID( resnumber, id::BB, 3 ), angle );
	return;
}

void betapeptide_setomega (
	core::pose::Pose &pose,
	core::Size resnumber,
	core::Real angle
) {
	using namespace core;
	using namespace core::id;
	pose.set_torsion( TorsionID( resnumber, id::BB, 4 ), angle );
	return;
}

void normalizevector (
	numeric::xyzVector<core::Real> &v
) {
	double normval = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	v[0] = v[0] / normval;
	v[1] = v[1] / normval;
	v[2] = v[2] / normval;

	return;
}

core::Real dotpdt (
	const numeric::xyzVector<core::Real> &v1,
	const numeric::xyzVector<core::Real> &v2
) {
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

bool use_in_rmsd(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	core::Size resno,
	core::Size atomno
) {
	if ( pose1.residue(resno).has( "N") && pose2.residue(resno).has( "N") && pose1.residue(resno).atom_index( "N")==atomno ) return true;
	if ( pose1.residue(resno).has("CA") && pose2.residue(resno).has("CA") && pose1.residue(resno).atom_index("CA")==atomno ) return true;
	if ( pose1.residue(resno).has("CM") && pose2.residue(resno).has("CM") && pose1.residue(resno).atom_index("CM")==atomno ) return true;
	if ( pose1.residue(resno).has( "C") && pose2.residue(resno).has( "C") && pose1.residue(resno).atom_index( "C")==atomno ) return true;
	if ( pose1.residue(resno).has( "O") && pose2.residue(resno).has( "O") && pose1.residue(resno).atom_index( "O")==atomno ) return true;

	return false;
}

void align_poses(
	core::pose::Pose &mypose,
	core::pose::Pose const &target
) {
	core::id::AtomID_Map< core::id::AtomID > amap;
	core::pose::initialize_atomid_map(amap,mypose,core::id::AtomID::BOGUS_ATOM_ID());
	for ( int ir = 1; ir <= (int)mypose.size(); ++ir ) {
		for ( int ia = 1; ia <= (int)mypose.residue(ir).nheavyatoms(); ++ia ) {
			if ( use_in_rmsd(mypose,target,ir,ia) ) {
				amap[core::id::AtomID(ia,ir)] = core::id::AtomID(ia,ir);
			}
		}
	}
	core::scoring::superimpose_pose( mypose, target, amap );
	return;
}

core::Real get_distance_measure(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2
) {
	std::map< core::id::AtomID, core::id::AtomID > amap;
	for ( int ir = 1; ir <= (int)pose1.size(); ++ir ) {
		for ( int ia = 1; ia <= (int)pose1.residue(ir).nheavyatoms(); ++ia ) {
			if ( use_in_rmsd(pose1,pose2,ir,ia) ) {
				amap.insert(std::pair< core::id::AtomID, core::id::AtomID> (core::id::AtomID(ia,ir), core::id::AtomID(ia, ir)));
			}
		}
	}
	return core::scoring::rms_at_all_corresponding_atoms(pose1,pose2,amap);
}

int main( int argc, char * argv [] ) {
	try {
		using namespace protocols;
		using namespace protocols::moves;
		using namespace core::scoring;
		using namespace core::optimization;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace utility::file;
		using namespace protocols::cluster;
		using namespace std;
		using namespace chemical;
		using namespace conformation;
		using namespace core::pose;
		//using namespace scoring::constraints;
		using namespace core::id;
		using namespace numeric;

		printf("Starting make_betahelicalbundle.cc\nFile created 22 September 2014 by Vikram K. Mulligan, Baker Laboratory.\n\n");
		fflush(stdout);

		//numeric::random::RandomGenerator RG( 923749 ); //Random generator and seed

		register_options();
		devel::init(argc, argv);
		//core::scoring::ScoreFunctionOP sfxn;
		//sfxn = core::scoring::getScoreFunction();

		core::Real R1;// = static_cast<core::Real>(option[v_R1]());
		core::Real omega1;// = 14.0*numeric::constants::d::pi/17.8;

		core::Real R0 = static_cast<core::Real>(option[v_R0]());
		core::Real omega0;// = 3.231/360.0*2.0*numeric::constants::d::pi;
		core::Real phi0 = 0;
		core::Real phi1 = 0;
		core::Real alpha;// = 12.0/360.0*2.0*numeric::constants::d::pi;
		core::Real Dz = 0;
		core::Real Dz_N = option[v_Dz_N]();
		core::Real Dz_C = option[v_Dz_C]();
		core::Real Dz_CM = option[v_Dz_CM]();
		core::Real phi1_N_denom = option[v_phi1_N_denom]();

		core::pose::Pose mypose;
		core::chemical::ResidueTypeSetCAP  fa_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		//core::pose::make_pose_from_sequence(mypose, "A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]", *fa_residue_set, true);

		const string sequence = "AAAA[B3A]AAAA[B3A]AAAA[B3A]AAAA[B3A]AAAA[B3A]AAAA[B3A]AAAA[B3A]AAAA[B3A]AAA";

		// grab residue types
		chemical::ResidueTypeCOPs requested_types = core::pose::residue_types_from_sequence( sequence, *( chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )), true );
		mypose.clear();

		for ( core::Size i = 1; i <= requested_types.size(); ++i ) {
			// grab the new residue
			chemical::ResidueType const & rsd_type = *requested_types[ i ];
			core::conformation::ResidueOP new_rsd( NULL );
			new_rsd = conformation::ResidueFactory::create_residue ( rsd_type );
			if ( i>1 ) mypose.append_residue_by_bond( *new_rsd, true );
			else {
				mypose.append_residue_by_jump(*new_rsd, 1, "", "", true);
			}
		}

		core::pose::Pose mypose2(mypose);

		numeric::xyzVector<core::Real> curatompos;
		numeric::xyzVector<core::Real> T1, T2, v1, v2;
		core::Real curphi, curpsi, curomega, curtheta;
		for ( core::Size i=1; i<90; i+=1 ) {
			omega0 = (-static_cast<core::Real>(i)*0.1-3.2)/360*2.0*numeric::constants::d::pi/R0*5.0;
			omega1 = 12.0*numeric::constants::d::pi/21.0 - omega0;
			alpha = asin(omega0)*3.6*R0/5.0;
			//alpha = omega0*4.0;
			//alpha = -(double)i/360.0*2.0*numeric::constants::d::pi;
			//R0 = (double)i*0.5;
			for ( core::Size ir=1; ir<=mypose.size(); ir++ ) {
				core::conformation::Residue rsd(mypose.residue(ir));

				Dz=0.0; R1=2.27; phi1=0.0;
				curatompos[0]=Xpos(ir,R0,omega0,phi0,R1,omega1,phi1,alpha,Dz);
				curatompos[1]=Ypos(ir,R0,omega0,phi0,R1,omega1,phi1,alpha,Dz);
				curatompos[2]=Zpos(ir,R0,omega0,R1,omega1,phi1,alpha,Dz);
				rsd.set_xyz("CA", curatompos);

				Dz=Dz_N; R1=1.61; phi1=-omega1/phi1_N_denom;
				curatompos[0]=Xpos(ir,R0,omega0,phi0,R1,omega1,phi1,alpha,Dz);
				curatompos[1]=Ypos(ir,R0,omega0,phi0,R1,omega1,phi1,alpha,Dz);
				curatompos[2]=Zpos(ir,R0,omega0,R1,omega1,phi1,alpha,Dz);
				rsd.set_xyz("N", curatompos);

				Dz=Dz_C; R1=1.68; phi1=omega1/4.0;
				curatompos[0]=Xpos(ir,R0,omega0,phi0,R1,omega1,phi1,alpha,Dz);
				curatompos[1]=Ypos(ir,R0,omega0,phi0,R1,omega1,phi1,alpha,Dz);
				curatompos[2]=Zpos(ir,R0,omega0,R1,omega1,phi1,alpha,Dz);
				rsd.set_xyz("C", curatompos);
				if ( mypose.residue( ir ).type().is_beta_aa() ) {
					Dz=Dz_CM; R1=1.68; phi1=omega1/2.0;
					curatompos[0]=Xpos(ir,R0,omega0,phi0,R1,omega1,phi1,alpha,Dz);
					curatompos[1]=Ypos(ir,R0,omega0,phi0,R1,omega1,phi1,alpha,Dz);
					curatompos[2]=Zpos(ir,R0,omega0,R1,omega1,phi1,alpha,Dz);
					rsd.set_xyz("CM", curatompos);
				}

				mypose.replace_residue(ir, rsd, false);
			}
			mypose.update_residue_neighbors();

			for ( core::Size ir=1; ir<=mypose.size(); ir++ ) {

				if ( mypose2.residue( ir ).type().is_beta_aa() ) {

					if ( ir>1 ) {
						dihedral_degrees( mypose.residue(ir-1).xyz("C"),
							mypose.residue(ir).xyz("N"),
							mypose.residue(ir).xyz("CA"),
							mypose.residue(ir).xyz("CM"),
							curphi);
						betapeptide_setphi(mypose2, ir, curphi);
					}
					dihedral_degrees( mypose.residue(ir).xyz("N"),
						mypose.residue(ir).xyz("CA"),
						mypose.residue(ir).xyz("CM"),
						mypose.residue(ir).xyz("C"),
						curtheta);
					betapeptide_settheta(mypose2, ir, curtheta);
					if ( ir<mypose.size() ) {
						dihedral_degrees( mypose.residue(ir).xyz("CA"),
							mypose.residue(ir).xyz("CM"),
							mypose.residue(ir).xyz("C"),
							mypose.residue(ir+1).xyz("N"),
							curpsi);
						betapeptide_setpsi(mypose2, ir, curpsi);
						dihedral_degrees( mypose.residue(ir).xyz("CM"),
							mypose.residue(ir).xyz("C"),
							mypose.residue(ir+1).xyz("N"),
							mypose.residue(ir+1).xyz("CA"),
							curomega);
						betapeptide_setomega(mypose2, ir, curomega);
					}
				} else {
					if ( ir>1 ) {
						dihedral_degrees( mypose.residue(ir-1).xyz("C"),
							mypose.residue(ir).xyz("N"),
							mypose.residue(ir).xyz("CA"),
							mypose.residue(ir).xyz("C"),
							curphi);
						betapeptide_setphi(mypose2, ir, curphi);
					}
					if ( ir<mypose.size() ) {
						dihedral_degrees( mypose.residue(ir).xyz("N"),
							mypose.residue(ir).xyz("CA"),
							mypose.residue(ir).xyz("C"),
							mypose.residue(ir+1).xyz("N"),
							curpsi);
						betapeptide_settheta(mypose2, ir, curpsi);
						dihedral_degrees( mypose.residue(ir).xyz("CA"),
							mypose.residue(ir).xyz("C"),
							mypose.residue(ir+1).xyz("N"),
							mypose.residue(ir+1).xyz("CA"),
							curomega);
						betapeptide_setpsi(mypose2, ir, curomega);
					}

				}
			}
			mypose2.update_residue_neighbors();
			printf("RMSD=%.8f\n", (double)get_distance_measure(mypose, mypose2));
			align_poses(mypose2, mypose);

			char outfile[25];
			sprintf(outfile, "cart%02lu.pdb", i);
			mypose.dump_pdb(outfile);
			sprintf(outfile, "tors%02lu.pdb", i);
			mypose2.dump_pdb(outfile);
		}
		printf("JOB COMPLETED.\n"); fflush(stdout);
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

