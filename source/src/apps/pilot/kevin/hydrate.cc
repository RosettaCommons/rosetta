#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <numeric/xyz.functions.hh>
#include <core/kinematics/Stub.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/Edge.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/HomogeneousTransform.hh>
#include <core/chemical/ResidueType.hh>

/*
using namespace core;
using namespace core::scoring;
using namespace core::pack;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
*/

static basic::Tracer TR( "hydrate" );

void place_water_acceptor(
	core::conformation::Residue const & rsd2,
	core::Size const accpt_atm,
	core::Real distance,
	core::Real LP_deg,
	core::Real tp5_about_acceptor_deg,
	core::Real tp5_about_hbond_deg,
	core::pose::Pose & pose
)
{
	// ? "base" ? one step up atomtree?
	core::Size const base1_accpt_atm( rsd2.atom_base( accpt_atm ) );
	// is abase2 base(base) ? 2 steps up atom tree?
	core::Size const base2_accpt_atm( rsd2.abase2( accpt_atm ) );

	// use this later
	// core::chemical::Hybridization const & hybrid( (rsd2.atom_type(accpt_atm)).hybridization() );

	// stub from these 3 atoms may rotate wrt acceptor LPs
	// ? get LP angle consistent with rest of pose?

	core::Vector const & accpt_xyz( rsd2.xyz( accpt_atm ) );
	core::Vector const & base1_accpt_xyz( rsd2.xyz( base1_accpt_atm ) );
	core::Vector const & base2_accpt_xyz( rsd2.xyz( base2_accpt_atm ) );

	// e1 points from carbonyl C along point to 0
	// if primary amide guarantees planarity then e2 in that plane
	core::kinematics::Stub stub( accpt_xyz, base1_accpt_xyz, accpt_xyz, base2_accpt_xyz );
	// for sp2, donating HOH:
	// 1) rotate about e3  by theta
	// 2) translate along e1, place H
	// 3) translate along e1, place O
	// 4) rotate about e1*cos(theta) by phi, place H
	// for placed HOH, one internal axis of rotation is free, change with homo trans

	// pivot on oxygen about axis normal to amide plane, turn toward a lone pair
	// LP position?? ridiculous VSEPR bunny ears for now
	// degree wrt C=O vector tilted away from N
	core::Real LP1_deg(LP_deg);
	numeric::xyzMatrix< core::Real > zrotmat = numeric::z_rotation_matrix_degrees( LP1_deg );

	// distance carbonyl oxygen acceptor to water oxygen
	// I have no idea, why not 2.5?
	core::Real dist_Osp2_Owat(distance);
	numeric::xyzVector< core::Real > crbnyl_O_to_wat_O( -1 * dist_Osp2_Owat, 0, 0);
	numeric::xyzVector< core::Real > wat_O_xyz = numeric::product(zrotmat, crbnyl_O_to_wat_O);

	// get a homo trans
	numeric::HomogeneousTransform< core::Real > HT_stub_to_wat(zrotmat, wat_O_xyz);

	// use this to rotate water about H-bond axis
	//core::Real tp5_about_hbond_deg(180);
	numeric::xyzMatrix< core::Real > tp5_about_hbond_rot_mat = numeric::x_rotation_matrix_degrees( tp5_about_hbond_deg );

	//core::Real tp5_about_acceptor_deg(180);
	numeric::xyzMatrix< core::Real > tp5_about_acceptor_rot_mat = numeric::x_rotation_matrix_degrees( tp5_about_acceptor_deg );

	core::pose::Pose tp5_pose;
	core::conformation::ResidueOP tp5_res( core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map("TP5") ) );
	tp5_pose.append_residue_by_jump( *tp5_res, 1 );
	// Play with TP5

	// place_tp5(xyzVector, HT, tp5)
	for ( core::Size ii = 1; ii <= tp5_res->natoms(); ii++ ) {
		std::string atm_name( tp5_res->type().atom_name( ii ) );
		core::Vector tp5_atm_coords(tp5_res->type().atom( ii ).ideal_xyz());
		tp5_res->set_xyz( ii, stub.local2global( numeric::product(tp5_about_acceptor_rot_mat, HT_stub_to_wat * numeric::product(tp5_about_hbond_rot_mat, tp5_res->type().atom( ii ).ideal_xyz() ) ) ) );
	}

	core::Size cutpoint ( pose.size() );
	pose.append_residue_by_jump( *tp5_res, cutpoint );
}


void place_water_donor(
	core::conformation::Residue const & rsd2,
	core::Size const don_atm,
	core::Real distance,
	core::Real tp5_about_don_deg,
	core::pose::Pose & pose
)
{
	// ? "base" ? one step up atomtree?
	core::Size const base1_don_atm( rsd2.atom_base( don_atm ) );
	// is abase2 base(base) ? 2 steps up atom tree?
	core::Size const base2_don_atm( rsd2.abase2( don_atm ) );

	// use this later
	// core::chemical::Hybridization const & hybrid( (rsd2.atom_type(accpt_atm)).hybridization() );

	// stub from these 3 atoms may rotate wrt acceptor LPs
	// ? get LP angle consistent with rest of pose?

	core::Vector const & don_xyz( rsd2.xyz( don_atm ) );
	core::Vector const & base1_don_xyz( rsd2.xyz( base1_don_atm ) );
	// why does this not exist?
	// core::Vector const & base2_don_xyz( rsd2.xyz( base2_don_atm ) );
	core::Vector const origin(0, 0, 0);

	// e1 points from carbonyl C along point to 0
	// if primary amide guarantees planarity then e2 in that plane
	//core::kinematics::Stub stub( don_xyz, base1_don_xyz, don_xyz, (origin - don_xyz).normalize() );
	core::kinematics::Stub stub( don_xyz, base1_don_xyz, (origin - don_xyz).normalize() );

	//numeric::xyzVector< core::Real > wat_O_xyz( -1 * distance, 0, 0);
	numeric::xyzVector< core::Real > wat_O_xyz( 1 * distance, 0, 0);

	// get a homo trans
	// numeric::HomogeneousTransform< core::Real > HT_stub_to_wat(zrotmat, wat_O_xyz);

	//core::Real tp5_about_acceptor_deg(180);
	numeric::xyzMatrix< core::Real > tp5_about_don_rot_mat = numeric::x_rotation_matrix_degrees( tp5_about_don_deg );

	core::pose::Pose tp5_pose;
	core::conformation::ResidueOP tp5_res( core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map("TP5") ) );
	core::kinematics::Stub tp5_stub( tp5_res->xyz("O"), tp5_res->xyz("EP1"), tp5_res->xyz("H1"));
	tp5_pose.append_residue_by_jump( *tp5_res, 1 );
	// Play with TP5
	for ( core::Size ii = 1; ii <= tp5_res->natoms(); ii++ ) {
		std::string atm_name( tp5_res->type().atom_name( ii ) );
		core::Vector tp5_atm_coords(tp5_res->type().atom( ii ).ideal_xyz());
		//tp5_res->set_xyz( ii, stub.local2global( numeric::product(tp5_about_don_rot_mat, tp5_stub.global2local( tp5_res->xyz(ii) ) ) ) + wat_O_xyz );
		tp5_res->set_xyz( ii, stub.local2global( numeric::product(tp5_about_don_rot_mat, tp5_stub.global2local( tp5_res->xyz(ii) ) + wat_O_xyz ) ) );
	}

	core::Size cutpoint ( pose.size() );
	pose.append_residue_by_jump( *tp5_res, cutpoint );
}


int main( int argc, char * argv [] ) {
	devel::init(argc, argv);

	// initialize task factory
	core::pack::task::TaskFactoryOP task_factory = new core::pack::task::TaskFactory;
	task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );
	// resfile
	if ( basic::options::option[ basic::options::OptionKeys::packing::resfile ].user() ) {
		task_factory->push_back( new core::pack::task::operation::ReadResfile );
	} else {
		core::pack::task::operation::RestrictToRepackingOP rtrop = new core::pack::task::operation::RestrictToRepacking;
		task_factory->push_back( rtrop );
	}
	if ( !basic::options::option[ basic::options::OptionKeys::in::file::s ].user() ) {
		// this function gets macro'd, argh
		utility_exit_with_message("Input PDB file not found");
	}

	// pose pose
	std::string pdb_file_name = basic::options::option[ basic::options::OptionKeys::in::file::s ]()[1];
	core::pose::Pose pose;
	core::import_pose::pose_from_file( pose, pdb_file_name , core::import_pose::PDB_file);
	pose.update_residue_neighbors();



	// go to residue 2, which happens to be a solvent exposed Q in 1ubq, add TP5
	core::Size const wet_residue_i(2);
	core::conformation::Residue const & rsd2( pose.residue(wet_residue_i) );

	// access pointer behavior?

	// at some point this will iterate over sites/atoms/residues
	core::chemical::AtomIndices::const_iterator accpt1_it = rsd2.accpt_pos().begin();
	// wander up to the sidechain carbonyl
	accpt1_it++;

	place_water_acceptor(rsd2, *accpt1_it, 2.5, 60, 0, 20, pose);
	// place_water_acceptor(rsd2, *accpt1_it, 2.5, 30, 60, 0, pose);
	place_water_acceptor(rsd2, *accpt1_it, 2.5, 60, 180, 0, pose);

	//core::chemical::AtomIndices::const_iterator don1_it = rsd2.accpt_pos().begin();
	// Why is this different!?
	core::chemical::AtomIndices::const_iterator don1_it = rsd2.Hpos_polar().begin();
	don1_it++;
	place_water_donor(rsd2, *don1_it, 3.0, 30, pose);
	don1_it++;
	place_water_donor(rsd2, *don1_it, 2.0, 0, pose);

	// write outputs
	//pose.dump_scored_pdb( "hydrate_"+pdb_file_name, *score_fxn);
	pose.dump_pdb( "hydrate_"+pdb_file_name);
	//pose.dump_pdb( "soak_"+water_file_name);

	TR << std::endl;
}
