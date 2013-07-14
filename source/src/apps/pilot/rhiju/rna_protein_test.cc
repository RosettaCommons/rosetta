// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueMatcher.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ResidueSelector.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

// AUTO-REMOVED #include <core/sequence/util.hh>

//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/rna/RNA_Util.hh>

//Mmmm.. constraints.
// AUTO-REMOVED #include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>

#include <core/kinematics/FoldTree.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
// AUTO-REMOVED #include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/rigid/RigidBodyMover.hh>

// AUTO-REMOVED #include <core/io/silent/RNA_SilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>
// AUTO-REMOVED #include <core/options/util.hh>

#include <core/options/option_macros.hh>

#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/util/basic.hh>

// AUTO-REMOVED #include <core/io/database/open.hh>


#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

#include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/conversions.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

//RNA stuff.
// AUTO-REMOVED #include <protocols/rna/RNA_FragmentsClasses.hh>
// AUTO-REMOVED #include <protocols/rna/RNA_DeNovoProtocol.hh>

//Job dsitributor
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef

#include <core/util/Tracer.hh>
using core::util::T;

// option key includes

#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/io/atom_tree_diffs/atom_tree_diff.hh>

#include <numeric/NumericTraits.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end



using core::util::Error;
using core::util::Warning;

using namespace core;
using namespace protocols;
using namespace core::options::OptionKeys;

using utility::vector1;

using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( Boolean, repack_test )
OPT_KEY( Boolean, prepack_test )
OPT_KEY( Boolean, pdbstats )
OPT_KEY( Boolean, rb_test )
OPT_KEY( Boolean, soft_rep )
OPT_KEY( Boolean, juke_sam )
OPT_KEY( Boolean, capri15 )
OPT_KEY( Boolean, no_repack )
OPT_KEY( Boolean, dump_the_pdb )

///////////////////////////////////////////////////////////////////////////////
void
rna_protein_repack_test()
{

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;


	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	pose::Pose pose;
	std::string pdb_file  = option[ in::file::s ][1];
	io::pdb::pose_from_pdb( pose, *rsd_set, pdb_file );

	dump_pdb( pose, "start.pdb");

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();
	for (Size i = 1; i <= pose.total_residue(); ++i) {
		task->nonconst_residue_task( i ).restrict_to_repacking();
	}

	ScoreFunctionOP scorefxn = getScoreFunctionLegacy( PRE_TALARIS_2013_STANDARD_WTS );
	scorefxn->set_weight( fa_pair, 0.0 );
	scorefxn->set_weight( fa_elec, 0.5 );
	//	scorefxn->energy_method_options().exclude_DNA_DNA( false );

	(*scorefxn)( pose );
	scorefxn->show( std::cout, pose );

	pack::pack_rotamers( pose, *scorefxn, task);
	dump_pdb( pose, "final.pdb" );
	scorefxn->show( std::cout, pose );


}

///////////////////////////////////////////////////////////////////////////////
void
rna_protein_prepack_test()
{

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;


	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	pose::Pose pose;
	std::string pdb_file  = option[ in::file::s ][1];
	io::pdb::pose_from_pdb( pose, *rsd_set, pdb_file );

	dump_pdb( pose, "start.pdb");

	//OK, phil.
	ScoreFunctionOP scorefxn = getScoreFunctionLegacy( PRE_TALARIS_2013_STANDARD_WTS );
	scorefxn->set_weight( fa_pair, 0.0 );
	scorefxn->set_weight( fa_elec, 0.5 );
	//	scorefxn->energy_method_options().exclude_DNA_DNA( false );



	//First push the protein and the RNA far apart.
	utility::vector1<Size> rna_residues, protein_residues;
	for (Size i = 1; i <= pose.total_residue(); ++i) {
		if (pose.residue(i).is_RNA() ) rna_residues.push_back( i );
		if (pose.residue(i).is_protein() ) protein_residues.push_back( i );
	}

	Size const num_protein = protein_residues.size() ;
	Size const upstream_residue = protein_residues[ num_protein/2 ];

	Size const num_rna = rna_residues.size() ;
	Size const downstream_residue = rna_residues[ num_rna/2 ];

	Size cutpoint( rna_residues[1] - 1 );
	if (protein_residues[1] > rna_residues[1] )  cutpoint = protein_residues[1]-1;

	kinematics::FoldTree f_original = pose.fold_tree();
	kinematics::FoldTree f_separate_protein_rna( pose.total_residue() );
	f_separate_protein_rna.new_jump( upstream_residue, downstream_residue, cutpoint );
	pose.fold_tree( f_separate_protein_rna );
	kinematics::Jump j( pose.jump( 1 ) );
	j.set_translation( Vector( 0.0, 0.0, 200.0 ) );
	pose.set_jump(1, j);

	(*scorefxn)(pose);

	dump_pdb( pose, "separated.pdb" );
	pose.fold_tree( f_original );

	////////////////////////////////////////////////////////////
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();
	for (Size i = 1; i <= pose.total_residue(); ++i) {
		task->nonconst_residue_task( i ).restrict_to_repacking();
	}
	(*scorefxn)( pose );
	scorefxn->show( std::cout, pose );

	pack::pack_rotamers( pose, *scorefxn, task);
	dump_pdb( pose, "prepack_"  +  pdb_file);
	scorefxn->show( std::cout, pose );


}




///////////////////////////////////////////////////////////////////////////////
void
pack_interface( pose::Pose & pose, scoring::ScoreFunction & scorefxn, ObjexxFCL::FArray1D<bool> & interface_residue ) {

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();

	Size const nres = pose.total_residue();

	//Switch to a simple fold_tree.
	kinematics::FoldTree f_save = pose.fold_tree();

	kinematics::FoldTree f ( nres );
	if ( options::option[ capri15 ] ) f.new_jump( 102+10,284+10,283+10 );
	pose.fold_tree( f );

	std::cout << "Packing residues: " ;
	for (Size i = 1; i <= nres; ++i) {
		task->nonconst_residue_task( i ).restrict_to_repacking();
		if ( interface_residue( i ) )  {
			std::cout << ' ' << i;
		} else {
			task->nonconst_residue_task( i ).prevent_repacking();
		}
	}

	std::cout << std::endl;

	pack::pack_rotamers( pose, scorefxn, task);

	scorefxn.show( std::cout, pose );

	pose.fold_tree( f_save );
}

///////////////////////////////////////////////////////////////////////////////
void
create_random_pose( pose::Pose & pose ){

	using namespace protocols::moves;

	Partner p( rigid::partner_downstream );
	if ( options::option[ capri15 ] ) p = partner_upstream; //Ooops, need RNA before protein/SAM.
	rigid::RigidBodyRandomizeMover mover( pose, 1 /* rb_jump_*/ ,  p );
	mover.apply( pose );

	//Translate to origin.
	kinematics::Jump j( pose.jump( 1 ) );
	j.set_translation( Vector( 0.0, 0.0, 0.0 ) );
	pose.set_jump( 1, j );

	// Keep jump atoms close, within some Gaussian sphere.
	Real trans_sigma = 4.0;
	if ( options::option[capri15] ) trans_sigma = 3.0;
	//	rigid::RigidBodyPerturbMoverOP rb_mover = new rigid::RigidBodyPerturbMover(
	//																															 pose, 1 /*jump_num*/, 0.0 /*rot*/, trans_sigma );
	//	rb_mover->apply( pose );
	using namespace numeric::random;
	Real const r = 2.5 + 1.0*uniform();
	Real const cos_theta = 2 * uniform() - 1;
	Real const sin_theta = std::sqrt( 1 - cos_theta*cos_theta );
	Real const phi = uniform() * 2 * numeric::NumericTraits<Real>::pi();
	Vector const rand_displacement = r * Vector( cos_theta, sin_theta*cos( phi ), sin_theta*sin( phi ) );
	j.set_translation( rand_displacement );
	pose.set_jump( 1, j );

}


///////////////////////////////////////////////////////////////////////////////
bool
check_protein_rna_clash( pose::Pose const & pose, ObjexxFCL::FArray1D< bool > & interface_residue ){

	static Real const dist_cutoff2 = 2.5 * 2.5;
	Size const nres = pose.total_residue();

	interface_residue.dimension( nres );
	interface_residue  = false;

	for (Size i=1; i <= nres; i++ ) {
		conformation::Residue const & rsd1 = pose.residue(i);
		if (!rsd1.is_RNA()) continue;

		//Quick check.
		Vector const & nbr_i( rsd1.xyz( rsd1.nbr_atom() ) );

		for (Size j=1; j <= pose.total_residue(); j++ ) {
			conformation::Residue const & rsd2 = pose.residue(j);
			if (rsd2.is_RNA()) continue;

			//Quick check.
			Vector const & nbr_j( rsd2.xyz( rsd2.nbr_atom() ) );

			Real const nbrcutoff = ( rsd1.nbr_radius() + rsd2.nbr_radius()) ;
			Real const nbrcutoff2 = nbrcutoff * nbrcutoff;

			if ( (nbr_i - nbr_j).length_squared() > nbrcutoff2 ) continue;

			interface_residue( i ) = true;
			interface_residue( j ) = true;

			//Now check for clashes between all RNA heavy atoms and protein backbone atoms.
			for (Size m = 1; m <= rsd1.nheavyatoms(); m++ ) {

				Vector const & atom_i( rsd1.xyz( m ) );

				// 5 atoms (up to C-alpha), except for glycine, I guess.
				Size atom_num_max =  std::min( Size( 5 ), rsd2.nheavyatoms() );

				//What if its the ligand? Look over all heavy atoms.
				if ( !rsd2.is_protein() ) atom_num_max = rsd2.nheavyatoms();

				for (Size n = 1; n <= atom_num_max; n++ ) {

					Vector const & atom_j( rsd2.xyz( n ) );

					if  ( ( atom_i - atom_j ).length_squared() <  dist_cutoff2 )  return true; //CLASH!!
				}
			}

		}
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////
void
rna_protein_rb_test(){

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;

	using namespace protocols::jobdist;

	pose::Pose pose;
	std::string pdb_file  = option[ in::file::s ][1];
	io::pdb::pose_from_pdb( pose, pdb_file );

	Size const nres = pose.total_residue();

	/////////////////////////////////////////////////////
	//Yea, this should *not* be hard-wired!
	// Instead could read in from a file, or at least
	// read in PDB's from two files and specify jump residues
	// from command line. Anyway...
	kinematics::FoldTree f( nres );
	if ( option[capri15] ) { //appropriate for rlm2 methyltransferase + sam + RNA.
		std::cout << "Setting up CAPRI15 fold tree " << nres << std::endl;
		assert( pose.residue( 5 ).aa() == na_rgu );
		f.new_jump( 5, 284+10, 10 );
		f.set_jump_atoms( 1, " N1 ", " CE ");
		f.new_jump( 102+10, 284+10, 283+10 );
	} else 	{
		std::cout << "Setting up 1zdi (MS2 coat protein/RNA) fold tree " << nres << std::endl;
		f.new_jump( 214, 267, 258 );
		f.set_jump_atoms( 1, " CB ", " O2*" );
	}
	pose.fold_tree( f );

	/////////////////////////////////////////////////////
	//Phil's suggestion for a fast score function
	ScoreFunctionOP scorefxn = getScoreFunctionLegacy( PRE_TALARIS_2013_STANDARD_WTS );

	if (option[soft_rep]) {
		scorefxn = ScoreFunctionFactory::create_score_function( SOFT_REP_WTS );
	}

	scorefxn->set_weight( fa_pair, 0.0 );
	scorefxn->set_weight( fa_elec, 0.5 );
	//	scorefxn->energy_method_options().exclude_DNA_DNA( false );

	dump_pdb( pose, "start.pdb");

	pose::Pose start_pose = pose;


	/////////////////////////////////////////////////////
	// Ian's magic silent file thingy.
	utility::vector1< BasicJobOP > input_jobs = load_s_and_l();
	std::string const atom_tree_diffs_file = option[ out::file::silent ] ;
	AtomTreeDiffJobDistributor jobdist( input_jobs, atom_tree_diffs_file );


	ObjexxFCL::FArray1D < bool > interface_residue( nres, false );

	////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////
	// Check score of start pose
	{
		pose = start_pose;
		std::cout << "Is start pose clashing? " << check_protein_rna_clash( pose, interface_residue ) << std::endl;
		//		pack_interface( pose, *scorefxn, interface_residue );
		std::string const tag( "START" );
		std::map< std::string, core::Real > scores;
		core::io::atom_tree_diffs::map_of_weighted_scores(pose, *scorefxn, scores);
		jobdist.dump_pose( tag, scores, start_pose, pose );
	}

	////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////
	// native pose
	//  hey, this probably goes well in the job distributor? anyway.
	pose::Pose native_pose;
	bool const use_native = option[in::file::native].active();
	if (use_native) 	io::pdb::pose_from_pdb( native_pose, option( in::file::native ) );

	////////////////////////////////////////////////////////////////////
	// MAIN LOOP
	Size const nstruct = option[ out::nstruct ];
	Size const num_tries = 50000;
	for ( Size n = 1; n <= nstruct; n++ ){

		std::cout << "Decoy time: " << n << std::endl;

		Size i( 1 );
		for ( i = 1; i <= num_tries; i++ ){
			pose = start_pose;
			create_random_pose( pose );
			if ( check_protein_rna_clash( pose, interface_residue ) ) continue;
			break;
		}

		std::cout << "Doing the repack after " << i << " orientations tested for clashes" <<  std::endl;
		if (!option[ no_repack] )	 pack_interface( pose, *scorefxn, interface_residue );


		//dump it.
		std::string const tag( "S_" + lead_zero_string_of( n, 5 ) );
		std::map< std::string, core::Real > scores;
		core::io::atom_tree_diffs::map_of_weighted_scores(pose, *scorefxn, scores);
		if ( use_native ){
			Real const rmsd = all_atom_rmsd( native_pose, pose );
			std::cout << "All atom rmsd: " << rmsd  << std::endl;
			scores[ "rms" ] =  rmsd;
		}
		jobdist.dump_pose( tag, scores, start_pose, pose );


		if (option[ dump_the_pdb] ) dump_pdb( pose, tag + ".pdb" );
	}

	//	rigid::RigidBodyPerturbMoverOP rb_mover = new rigid::RigidBodyPerturbMover(
	//			pose, 1 /*jump_num*/, 60.0 /*rot*/, 0.0 /*trans*/ );
	//	rb_mover->apply( pose );
	//	dump_pdb( pose, "perturb.pdb" );


	return;
}






///////////////////////////////////////////////////////////////////////////////
// Note that, in principle, the centroid could be a virtual atom,
// and I could let the atom-tree folding machinery figure out where the hell it is.
//
// Also, this copies some code from Phil's dna/base_geometry.cc
//
Vector
get_centroid( conformation::Residue const & rsd )
{

  Vector centroid( 0.0 );
  Size numatoms = 0;
  for ( Size i=rsd.first_sidechain_atom(); i<= rsd.nheavyatoms(); ++i ) {
    centroid += rsd.xyz(i);
    numatoms++;
  }
  if (numatoms > 0 ) {
		centroid /= static_cast< Real >( numatoms );
	} else { //Yo, is this a glycine?
		assert( rsd.aa() == chemical::aa_gly );
		centroid = rsd.xyz( "CA" );
	}

  return centroid;
}

///////////////////////////////////////////////////////////////////////////////
kinematics::Stub
get_base_coordinate_system( conformation::Residue const & rsd, Vector const & centroid )
{
  using namespace chemical;
  Size res_type = rsd.aa();

  assert( rsd.is_RNA() );

  Vector x,y,z;

  // Make an axis pointing from base centroid to Watson-Crick edge.
  std::string WC_atom;
  if ( res_type == na_rad ) WC_atom = " N1 ";
  if ( res_type == na_rcy ) WC_atom = " N3 ";
  if ( res_type == na_rgu ) WC_atom = " N1 ";
  if ( res_type == na_ura ) WC_atom = " N3 ";

  Vector const WC_coord (rsd.xyz( WC_atom ) );
  x = WC_coord - centroid;
  x.normalize();

  // Make a perpendicular axis pointing from centroid towards
  // Hoogstein edge (e.g., major groove in a double helix).
  std::string H_atom;
  if ( res_type == na_rad ) H_atom = "N7";
  if ( res_type == na_rcy ) H_atom = "C5";
  if ( res_type == na_rgu ) H_atom = "N7";
  if ( res_type == na_ura ) H_atom = "C5";

  Vector const H_coord (rsd.xyz( H_atom ) );
  y = H_coord - centroid; //not orthonormal yet...
  z = cross(x, y);
  z.normalize(); // Should poSize roughly 5' to 3' if in a double helix.

  y = cross(z, x);
  y.normalize(); //not necessary but doesn't hurt.

  //  std::cout << "WC : " << WC_coord << "   H : " << H_coord << "    centroid: " << centroid << std::endl;

  return kinematics::Stub( Matrix::cols( x, y, z ), centroid );
}


void
check_contact_and_output( conformation::Residue const & rsd1,
													conformation::Residue const & rsd2,
													Vector const & centroid_i,
													kinematics::Stub const & stub_i,
													Vector const & atom_j,
													Size const & n,
													utility::io::ozstream & base_out )
{

	static Real const Z_CUTOFF = 5.0;
	static Real const RHO_CUTOFF = 10.0;

	Matrix const & M_i( stub_i.M );
	Vector const & x_i = M_i.col_x();
	Vector const & y_i = M_i.col_y();
	Vector const & z_i = M_i.col_z();


	Vector dist_vector = atom_j - centroid_i;
	Real const x = dot( dist_vector, x_i   );
	Real const y = dot( dist_vector, y_i   );
	Real const z = dot( dist_vector, z_i   );
	Real const rho = std::sqrt( x*x + y*y );

	//Basic geometry cuts.
	if ( std::abs(z) > Z_CUTOFF ) return;
	if ( rho > RHO_CUTOFF ) return;

	////////////////////////////////
	// Do closeness check.
	////////////////////////////////
	// Now loop over atoms to see if there's a contact.
	bool there_is_a_contact = false;
	static Real const CONTACT_CUTOFF2 = 5.0*5.0;
	for (Size m = rsd1.first_sidechain_atom(); m <= rsd1.nheavyatoms(); m++ ) {
		Vector const & atom_i( rsd1.xyz( m ) );
		if ( (atom_i - atom_j).length_squared() < CONTACT_CUTOFF2 ){
			there_is_a_contact = true;
			break;
		}
	}
	if (!there_is_a_contact) return;


	base_out <<
		I(3, Size( rsd1.aa() )) << " " <<
		I(4, rsd1.seqpos()) <<  " " <<
		I(3, Size( rsd2.aa() ) ) << " " <<
		I(4, rsd2.seqpos()) <<  " " <<
		I(3, n) << " " <<
		F(8,4,x) << " " << F(8,4,y) << " " << F(8,4,z) <<
		std::endl;

}

///////////////////////////////////////////////////////////////////////////////
void
check_oxygen_contact_and_output( conformation::Residue const & rsd1,
																 conformation::Residue const & rsd2,
																 Vector const & atom_i, Vector const & atom_j,
																 Size const & m,
																 Size const & n,
																 utility::io::ozstream & oxygen_out )
{
	static Real const DIST_CUTOFF2 = 8.0;
	Distance dist = (atom_i-atom_j).length();
	if ( dist < DIST_CUTOFF2 ) {
	oxygen_out <<
		I(3, Size( rsd1.aa() )) << " " <<
		I(4, rsd1.seqpos()) <<  " " <<
		I(3, Size( rsd2.aa() ) ) << " " <<
		I(4, rsd2.seqpos()) <<  " " <<
		I(3, m) << " " <<
		I(3, n) << " " <<
		F(8,4,dist) <<
		std::endl;

	}

}


///////////////////////////////////////////////////////////////////////////////
void
rna_protein_pdbstats_test(){

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;

	//Should make this a loop over several pdb files.
	pose::Pose pose;
	std::string pdb_file  = option[ in::file::s ][1];
	io::pdb::pose_from_pdb( pose, pdb_file );

	Size const nres = pose.total_residue();

	if (nres < 1000 ) dump_pdb( pose, "start.pdb" );

	//////////////////////////
	// Silly guide for me:
	for (Size i = 1; i < num_aa_types; i++ ){
		std::cout << I(3,i) << " ==> " <<
			name_from_aa( chemical::AA(i) ) << std::endl;
	}

	utility::vector1 < std::string > RNA_backbone_oxygen_atoms;
	RNA_backbone_oxygen_atoms.push_back( " O1P");
	RNA_backbone_oxygen_atoms.push_back( " O2P");
	RNA_backbone_oxygen_atoms.push_back( " O5*");
	RNA_backbone_oxygen_atoms.push_back( " O4*");
	RNA_backbone_oxygen_atoms.push_back( " O3*");
	RNA_backbone_oxygen_atoms.push_back( " O2*");
	Size const num_RNA_backbone_oxygen_atoms = RNA_backbone_oxygen_atoms.size();

	//////////////////////////
	// Silly guide for me:
	std::cout << std::endl;
	for (Size i = 1; i <= num_RNA_backbone_oxygen_atoms; i++ ){
		std::cout << I(3,i) << " ==> " <<
			RNA_backbone_oxygen_atoms[i] << std::endl;
	}


	//Precalculate centroids.
	utility::vector1 < Vector> centroids;
	for (Size i=1; i <= nres; i++ ) {
		conformation::Residue const & rsd1 = pose.residue(i);
		centroids.push_back( get_centroid( rsd1 ) );
	}


	std::string const outfile = option[ out::file::o ];
	utility::io::ozstream base_out( outfile );
	utility::io::ozstream oxygen_out( "oxygen_"+outfile );

	for (Size i=1; i <= nres; i++ ) {
		conformation::Residue const & rsd1 = pose.residue(i);
		if (!rsd1.is_RNA()) continue;

		Vector const & centroid_i = centroids[i];
		kinematics::Stub const & stub_i = get_base_coordinate_system( rsd1, centroid_i );

		//Quick check.
		Vector const & nbr_i( rsd1.xyz( rsd1.nbr_atom() ) );

		for (Size j=1; j <= pose.total_residue(); j++ ) {
			conformation::Residue const & rsd2 = pose.residue(j);
			if (!rsd2.is_protein()) continue;

			//Quick check.
			Vector const & nbr_j( rsd2.xyz( rsd2.nbr_atom() ) );

			Real const nbrcutoff = ( rsd1.nbr_radius() + rsd2.nbr_radius()) ;
			Real const nbrcutoff2 = nbrcutoff * nbrcutoff;

			if ( (nbr_i - nbr_j).length_squared() > nbrcutoff2 ) continue;

			Vector const centroid_j = get_centroid( rsd2 );

			/////////////////////////////////////////////////////////////////////////////////
			// BASE stuff
			//Save info on x,y,z position of N, CA, CB, C, O, centroid --> to text file.
			// 5 atoms (up to C-alpha), except for glycine, I guess.
			for (Size n = 1; n <= std::min( Size( 5 ), rsd2.nheavyatoms() ); n++ ) {
				Vector const & atom_j( rsd2.xyz( n ) );
				check_contact_and_output( rsd1, rsd2, centroid_i, stub_i, atom_j, n, base_out );
			} // backbone atoms.

			check_contact_and_output( rsd1, rsd2, centroid_i, stub_i, centroid_j, 0, base_out );

			//////////////////////////////////////////////////////////////
			//Also distances to important RNA backbone oxygen atoms...
			//////////////////////////////////////////////////////////////
			for (Size m = 1 ; m <= num_RNA_backbone_oxygen_atoms; m++ ){
				Vector const & atom_i( rsd1.xyz( RNA_backbone_oxygen_atoms[m] ) );
				for (Size n = 1; n <= std::min( Size( 5 ), rsd2.nheavyatoms() ); n++ ) {
					Vector const & atom_j( rsd2.xyz( n ) );
					check_oxygen_contact_and_output( rsd1, rsd2, atom_i, atom_j, m, n, oxygen_out );
				} // backbone atoms.
				check_oxygen_contact_and_output( rsd1, rsd2, atom_i, centroid_j, m, 0, oxygen_out );
			}

		}// j
	}// i

	base_out.close();
	oxygen_out.close();


}

///////////////////////////////////////////////////////////////////////////////
void
juke_sam_test(){

	using namespace core::conformation;
	using namespace core::kinematics;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::id;
	using namespace core::optimization;

	//Should make this a loop over several pdb files.
	pose::Pose pose;
	utility::vector1 <std::string> pdb_files ( option[ in::file::s ]() );

	for (Size n = 1; n <= pdb_files.size(); n++ ){

		std::string const pdb_file = pdb_files[n];
		io::pdb::pose_from_pdb( pose, pdb_file );

		Size const nres = pose.total_residue();

		Size const sam = nres;
		std::cout <<"************************************************" << std::endl;
		std::cout << "HEY is this a SAM? " << sam << " ==> " << pose.residue( sam ).name3() << std::endl;
		std::cout <<"************************************************" << std::endl;

		//New fold tree -- connect sam to conserved glycine [102].
		FoldTree f( nres );
		f.new_jump( 102, sam, sam - 1 );
		pose.fold_tree( f );

		dump_pdb( pose, "start.pdb" );

		// Set up constraints
		ConstraintSetOP cst_set( new ConstraintSet() );
		Residue const  & sam_rsd = pose.residue( sam );
		cst_set->add_constraint( new AtomPairConstraint( AtomID( sam_rsd.atom_index( " C5'" ), sam ),
																										 AtomID( pose.residue(169).atom_index( " CG ") /*pro*/,  169 ) ,
																										 new HarmonicFunc( 3.5 /*anchor*/, 0.4 /*stdev*/) ) );
		cst_set->add_constraint( new AtomPairConstraint( AtomID( sam_rsd.atom_index( " N6 " ), sam ),
																										 AtomID( pose.residue(149).atom_index( " OD1") /*asp*/,  149 ) ,
																										 new HarmonicFunc( 3.01 /*anchor*/, 0.4 /*stdev*/) ) );
		cst_set->add_constraint( new AtomPairConstraint( AtomID( sam_rsd.atom_index( " N1 " ), sam ),
																										 AtomID( pose.residue(150).atom_index( " N  ") /*ile*/,  150 ) ,
																										 new HarmonicFunc( 3.01 /*anchor*/, 0.4 /*stdev*/) ) );
		cst_set->add_constraint( new AtomPairConstraint( AtomID( sam_rsd.atom_index( " C4 " ), sam ),
																										 AtomID( pose.residue(126).atom_index( " CG2") /*ile*/,  126 ) ,
																										 new HarmonicFunc( 3.44 /*anchor*/, 0.4 /*stdev*/) ) );
		cst_set->add_constraint( new AtomPairConstraint( AtomID( sam_rsd.atom_index( " N1 " ), sam ),
																										 AtomID( pose.residue(126).atom_index( " CD1") /*ile*/,  126 ) ,
																										 new HarmonicFunc( 4.34 /*anchor*/, 0.4 /*stdev*/) ) );
		cst_set->add_constraint( new AtomPairConstraint( AtomID( sam_rsd.atom_index( " O2'" ), sam ),
																										 AtomID( pose.residue(125).atom_index( " OD2") /*asp*/,  125 ) ,
																										 new HarmonicFunc( 2.44 /*anchor*/, 0.4 /*stdev*/) ) );
		cst_set->add_constraint( new AtomPairConstraint( AtomID( sam_rsd.atom_index( " C1'" ), sam ),
																										 AtomID( pose.residue(102).atom_index( " CA ") /*gly*/,  102 ) ,
																										 new HarmonicFunc( 3.79 /*anchor*/, 0.4 /*stdev*/) ) );
		cst_set->add_constraint( new AtomPairConstraint( AtomID( sam_rsd.atom_index( " O3'" ), sam ),
																										 AtomID( pose.residue(104).atom_index( " CA ") /*gly*/,  104 ) ,
																										 new HarmonicFunc( 2.93 /*anchor*/, 0.4 /*stdev*/) ) );
		cst_set->add_constraint( new AtomPairConstraint( AtomID( sam_rsd.atom_index( " O4'" ), sam ),
																										 AtomID( pose.residue(102).atom_index( " CA ") /*gly*/,  102 ) ,
																										 new HarmonicFunc( 3.11 /*anchor*/, 0.4 /*stdev*/) ) );
		cst_set->add_constraint( new AtomPairConstraint( AtomID( sam_rsd.atom_index( " O3'" ), sam ),
																										 AtomID( pose.residue(100).atom_index( " OD2") /*asp*/,  100 ) ,
																										 new HarmonicFunc( 3.64 /*anchor*/, 0.4 /*stdev*/) ) );
		cst_set->add_constraint( new AtomPairConstraint( AtomID( sam_rsd.atom_index( " SD " ), sam ),
																										 AtomID( pose.residue(100).atom_index( " O  ") /*phe*/,  167 ) ,
																										 new HarmonicFunc( 3.68 /*anchor*/, 0.4 /*stdev*/) ) );
		pose.constraint_set( cst_set );

		// Turn off rep and minimize.
		ScoreFunctionOP scorefxn_hard = getScoreFunctionLegacy( PRE_TALARIS_2013_STANDARD_WTS );
		scorefxn_hard->set_weight( fa_pair, 0.0 );
		scorefxn_hard->set_weight( fa_elec, 0.5 );
		scorefxn_hard->set_weight( atom_pair_constraint, 1.0 );

		ScoreFunctionOP scorefxn_soft = ScoreFunctionFactory::create_score_function( SOFT_REP_WTS );
		scorefxn_soft->set_weight( fa_pair, 0.0 );
		scorefxn_soft->set_weight( fa_elec, 0.5 );
		scorefxn_soft->set_weight( atom_pair_constraint, 1.0 );

		ScoreFunctionOP scorefxn_soft_norep = scorefxn_soft->clone();
		scorefxn_soft_norep->set_weight( fa_rep, 0.1 );

		kinematics::MoveMap mm;
		mm.set_bb( false );
		mm.set_chi( false );
		mm.set_jump( true );

		scorefxn_soft_norep->show( std::cout, pose );

		AtomTreeMinimizer minimizer;
		float const dummy_tol( 0.0000025);
		bool const use_nblist( true ), deriv_check( false );
		MinimizerOptions options( "dfpmin_armijo_nonmonotone", dummy_tol, use_nblist, deriv_check, deriv_check );
		options.nblist_auto_update( true );
		minimizer.run( pose, mm, *scorefxn_soft_norep, options );

		scorefxn_soft_norep->set_weight( fa_rep, 0.2 );
		minimizer.run( pose, mm, *scorefxn_soft_norep, options );

		minimizer.run( pose, mm, *scorefxn_soft, options );

		minimizer.run( pose, mm, *scorefxn_hard, options );


		dump_pdb( pose, "first_min.pdb" );

		scorefxn_hard->show( std::cout, pose );
		//	return;

		// soft repack
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			task->nonconst_residue_task( ii ).restrict_to_repacking();
		}
		pack::pack_rotamers( pose, *scorefxn_soft, task);

		dump_pdb( pose, "soft_repack.pdb" );

		// hard minimize
		mm.set_chi( true );
		minimizer.run( pose, mm, *scorefxn_hard, options );
		mm.set_bb( true );
		minimizer.run( pose, mm, *scorefxn_hard, options );
		mm.set_jump( true );
		minimizer.run( pose, mm, *scorefxn_hard, options );

		dump_pdb( pose, "hard_min.pdb" );
		scorefxn_hard->show( std::cout, pose );

		// soft repack
		pack::pack_rotamers( pose, *scorefxn_soft, task);
		dump_pdb( pose, "soft_repack2.pdb" );

		// turn off constraints, hard minimize
		minimizer.run( pose, mm, *scorefxn_hard, options );
		dump_pdb( pose, "hard_min2.pdb" );
		scorefxn_hard->show( std::cout, pose );

		// turn off constraints, hard minimize
		scorefxn_hard->set_weight( atom_pair_constraint, 0.0 );
		minimizer.run( pose, mm, *scorefxn_hard, options );
		dump_pdb( pose, "hard_min3.pdb" );
		scorefxn_hard->show( std::cout, pose );

		// Save the sucker.
		dump_pdb( pose, "final.pdb" );
		dump_pdb( pose, "sam_juke_"+pdb_file );
	}

}



///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace core::options;

	//Uh, options?
	NEW_OPT( repack_test, "Test repack with RNA/protein system", false );
	NEW_OPT( prepack_test, "Slide faraway and repack RNA/protein system", false );
	NEW_OPT( rb_test,     "Move RNA relative to protein", false );
	NEW_OPT( pdbstats,    "Get stats for RNA/protein interactions", false );
	NEW_OPT( juke_sam,    "Juke SAM for CAPRI15", false );
	NEW_OPT( soft_rep,    "Soft repulsive instead of hard rep", false );
	NEW_OPT( capri15,    "Simulation parameters for CAPRI15 Target T033", false );
	NEW_OPT( no_repack, "HACK OPTION", false );
	NEW_OPT( dump_the_pdb, "HACK OPTION", false );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	if ( option[ repack_test ] ) {
		rna_protein_repack_test();
	} else if ( option[ rb_test ] ) {
		rna_protein_rb_test();
	}	else if ( option[ prepack_test ] ) {
		rna_protein_prepack_test();
	} else if ( option[ pdbstats ] ) {
		rna_protein_pdbstats_test();
	} else if ( option[ juke_sam ] ) {
		juke_sam_test();
	}
	exit( 0 );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}
