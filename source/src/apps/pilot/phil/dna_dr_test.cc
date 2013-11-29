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
#include <devel/dna/protocols.hh>
#include <devel/dna/util.hh>
//#include <devel/dna/util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/frags/TorsionFragment.hh>
#include <utility/excn/Exceptions.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/rigid_body_moves.hh>

#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/dna/DNA_BasePotential.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/etable/Etable.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/AtomVDW.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/elec/FA_ElecEnergy.hh>
//#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
//#include <core/scoring/etable/count_pair/CountPair1BC4.hh>

#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/VariantType.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/AA.hh>

#include <core/conformation/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/rotamer_set/WaterPackingInfo.hh>

#include <core/kinematics/FoldTree.hh>
#include <protocols/viewer/visualize.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/id/AtomID_Map.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/options/util.hh>
//#include <basic/options/after_opts.hh>

#include <basic/prof.hh> // profiling
#include <basic/basic.hh>
#include <core/id/SequenceMapping.hh>

#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <set>
#include <cstdlib>
#include <sstream>
#include <math.h>

//silly using/typedef


#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/phil.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>



using basic::T;
using basic::Error;
using basic::Warning;

//static numeric::random::RandomGenerator RG(12323); // <- Magic number, do not change it!!!

using namespace core;
using namespace protocols;

using utility::vector1;
using std::string;
using std::cout;
using std::endl;
using io::pdb::dump_pdb;


static basic::Tracer tt( "demo.phil.dna_dr_test", basic::t_info );

///////////////////////////////////////////////////////////////////////////////

void
kono_sarai_stats()
{
	using utility::vector1;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace scoring;
	using namespace scoring::dna;
	using namespace conformation;
	using namespace chemical;
	using namespace optimization;
	using namespace pose;

	vector1< std::string > files( basic::options::start_files() );

	for ( Size nn=1; nn<= files.size(); ++nn ) {
		Pose pose;
		core::import_pose::pose_from_pdb( pose, files[nn] );
		if ( pose.empty() ) continue;

// 		set_base_partner( pose ); // fills base partner info
// 		BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );

// 		std::string const tag( lead_zero_string_of( nn, 4 )+".pdb" );

// 		io::pdb::dump_pdb( pose, tag );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {

			Residue const & rsd( pose.residue(i) );
			if ( !rsd.is_DNA() ) continue;

			// define coordinate frame
			bool const AG( rsd.aa() == na_ade || rsd.aa() == na_gua );

			Vector const origin( AG ? rsd.xyz("N9") : rsd.xyz("N1") );
			Vector p( AG ? rsd.xyz("C4") : rsd.xyz("C2") );

			Vector const x( ( p - origin ).normalized() );
			bool const old_way( true );
			Vector z,y;
			if ( old_way ) {
				z = ( get_z_axis( rsd, x ) );
				y = ( z.cross( x ) );
			} else {
				Vector q( AG ? rsd.xyz("N7") : rsd.xyz("C5") );
				z = ( ( x.cross( q - origin ) ).normalized() );
				y = ( z.cross( x ) );
			}

			for ( Size j=1; j<= pose.total_residue(); ++j ) {

				Residue const & rsd2( pose.residue(j) );
				if ( !rsd2.is_protein() ) continue;

				Vector const calpha( rsd2.xyz("CA") - origin );
				Real const xx( dot(calpha,x) );
				Real const yy( dot(calpha,y) );
				Real const zz( dot(calpha,z) );

				if ( ( xx >= -13.5 && xx <= 13.5  ) &&
						 ( yy >= -13.5 && yy <= 13.5  ) &&
						 ( zz >=  -6.0 && zz <=  6.0  ) ) {
					int aa_num = rsd2.aa() - 1;
					int na_num = rsd.aa() - chemical::first_DNA_aa;
					std::cout << "CONTACT " << F(9,3,xx) << F(9,3,yy) << F(9,3,zz) <<
						' ' << rsd.name1() << ' ' << na_num << ' ' << rsd2.name1() << ' ' << aa_num << ' ' << files[nn] << I(4,i) << I(4,j) << std::endl;
				}
			}
		}
	}
}

std::string
atom_line( std::string const & atom_name, std::string const & rsd_name, Vector const & xyz, int & counter )
{
	std::ostringstream os;
	assert( atom_name.size() == 4 && rsd_name.size() == 3 );
	++counter;
	os << "ATOM  " << I(6,counter) << atom_name << ' ' << rsd_name << " A" << I(4,counter) << "    " <<
		F(8,3,xyz(1)) <<F(8,3,xyz(2)) <<F(8,3,xyz(3)) << '\n';
	return os.str();
}

///////////////////////////////////////////////////////////////////////////////

void
kono_sarai_zscore()
{
	using utility::vector1;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace scoring;
	using namespace scoring::dna;
	using namespace conformation;
	using namespace chemical;
	using namespace optimization;
	using namespace pose;


	vector1< std::string > files( basic::options::start_files() );

	for ( Size nn=1; nn<= files.size(); ++nn ) {
		Pose pose;
		core::import_pose::pose_from_pdb( pose, files[nn] );
		if ( pose.empty() ) continue;
		Pose start_pose = pose;

 		set_base_partner( pose ); // fills base partner info
 		BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );

		pose.dump_pdb( "test.pdb");

// 		std::string const tag( lead_zero_string_of( nn, 4 )+".pdb" );

// 		io::pdb::dump_pdb( pose, tag );

		utility::vector1< int > motif_pos;
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			Residue const & rsd( pose.residue(i) );
			if ( rsd.is_DNA() && partner[i]>i ) {
				motif_pos.push_back(i);
			}
			if ( rsd.is_DNA() && partner[i] == 0 ) {
				std::cout << "unpaired base: " << files[nn] << ' ' << i << std::endl;
			}
		}

		for ( int k=0; k<500; k++ ) // number of random sequences (set as commandline option)
		{
			std::ofstream out( "zscore.tmp" );
			// recover the starting conformation
			pose = start_pose;

			set_base_partner( pose );

			// randomize the dna positions
			for ( Size i=1; i<= motif_pos.size(); ++i ) {
				devel::dna::randomize_motif_sequence( pose, motif_pos[i], 1 );
			}

			//std::stringstream ss;
			//std::string str;
			//ss << k;
			//ss >> str;
			//pose.dump_pdb( "test" + str + ".pdb" );
			//std::cout << k << std::endl;

			for ( Size i=1; i<= pose.total_residue(); ++i ) {

				Residue const & rsd( pose.residue(i) );
				if ( !rsd.is_DNA() ) continue;

				// define coordinate frame
				bool const AG( rsd.aa() == na_ade || rsd.aa() == na_gua );

				Vector const origin( AG ? rsd.xyz("N9") : rsd.xyz("N1") );
				Vector p( AG ? rsd.xyz("C4") : rsd.xyz("C2") );

				Vector const x( ( p - origin ).normalized() );
				bool const old_way( true );
				Vector z,y;
				if ( old_way ) {
					z = ( get_z_axis( rsd, x ) );
					y = ( z.cross( x ) );
				} else {
					Vector q( AG ? rsd.xyz("N7") : rsd.xyz("C5") );
					z = ( ( x.cross( q - origin ) ).normalized() );
					y = ( z.cross( x ) );
				}

				for ( Size j=1; j<= pose.total_residue(); ++j ) {

					Residue const & rsd2( pose.residue(j) );
					if ( !rsd2.is_protein() ) continue;

					Vector const calpha( rsd2.xyz("CA") - origin );
					Real const xx( dot(calpha,x) );
					Real const yy( dot(calpha,y) );
					Real const zz( dot(calpha,z) );

					// output contact info to tempfile

					if ( ( xx >= -13.5 && xx <= 13.5  ) &&
							 ( yy >= -13.5 && yy <= 13.5  ) &&
							 ( zz >=  -6.0 && zz <=  6.0  ) ) {

						out << "CONTACT " << F(9,3,xx) << F(9,3,yy) << F(9,3,zz) <<
							' ' << rsd.name1() << rsd2.name1() << ' ' << files[nn] << std::endl;
						//std::cout << "CONTACT " << F(9,3,xx) << F(9,3,yy) << F(9,3,zz) <<
						//	' ' << rsd.name1() << rsd2.name1() << ' ' << files[nn] << I(4,i) << I(4,j) << std::endl;
					}
				}
			}
			out.close();
			// correct place to call modified sarai_contacts python script on temp file
			system("python /mnt/cpb_home/pbat/python_amy/sarai_contact_zscores_rosetta.py");
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void
loop_modeling_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace scoring;
	using namespace scoring::constraints;
	using namespace conformation;
	using namespace chemical;
	using namespace optimization;
	using namespace id;
	using namespace io::pdb;
	using namespace pose;
	using namespace protocols::loops;
	using namespace protocols::frags;


	// read the pose
	Pose pose;
	core::import_pose::centroid_pose_from_pdb( pose, basic::options::start_file() );

	dump_pdb( pose, "premap.pdb" );

	{ // hacking

		for ( int na=first_DNA_aa; na <= last_DNA_aa; ++na ) {
			Real max_d( 0.0 );

			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				Residue const & rsd( pose.residue(i) );
				if ( rsd.aa() != AA(na) ) continue;

				Vector const & nbr_xyz( rsd.nbr_atom_xyz() );

				for ( Size j=1; j<= rsd.natoms(); ++j ) {
					Real const d( nbr_xyz.distance( rsd.xyz(j) ) );
					if ( d > max_d ) max_d = d;
				}
			}
			std::cout << "MAX_D: " << AA(na) << ' ' << max_d << std::endl;
		}
		//exit(0);
	}

	// read the alignment file
	id::SequenceMapping mapping;
	std::string source_seq, target_seq;
	read_alignment_file( option[ phil::align_file ], source_seq, target_seq, mapping );


	// may not represent the entire pose source sequence
	std::string const pose_seq( pose.sequence() );

	if ( source_seq != pose_seq ) {
		/// source sequence from align file does not cover entire pdb file sequence

		if ( pose_seq.find( source_seq ) == std::string::npos ) {
			utility_exit_with_message( "alignfile source sequence not contained in pose sequence" );
		}
		Size const offset( pose_seq.find( source_seq ) );

		std::string source_nterm_seq, target_nterm_seq;
		for ( Size i=1; i<= offset; ++i ) {
			// this is a residue in the input pdb that's not accounted for in the alignment
			// if its protein, unalign it
			// if its dna, align it
			Residue const & rsd( pose.residue( i ) );

			mapping.insert_residue( i );
			source_nterm_seq += rsd.name1();

			if ( rsd.is_DNA() ) {
				target_nterm_seq += rsd.name1();
				Size const target_pos( target_nterm_seq.size() );
				mapping.insert_target_residue( target_pos );
				mapping[i] = target_pos;
			} else {
				assert( rsd.is_protein() );
				// ligand etc could be tricky when reading an alignment from a file
			}
		}

		source_seq = source_nterm_seq + source_seq;
		target_seq = target_nterm_seq + target_seq;

		assert( pose_seq.find( source_seq ) == 0 );

		while( source_seq.size() < pose_seq.size() ) {
			assert( mapping.size1() == source_seq.size() );

			Size const pos( source_seq.size() + 1 );
			Residue const & rsd( pose.residue( pos ) );

			char const n1( rsd.name1() );
			source_seq += n1;
			mapping.push_back( 0 );

			if ( rsd.is_DNA() ) {
				target_seq += n1;
				mapping.insert_target_residue( target_seq.size() );
				mapping[pos] = target_seq.size();
			}
		}

		assert( pose_seq == source_seq );
	}

	//mapping.size2( target_seq.size() );

	std::cout << pose.sequence() << '\n' << target_seq << '\n';

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		std::cout << "M1: " << I(4,i) << ' ' << pose.residue(i).name() << I(4,mapping[i]);
		if ( mapping[i] ) {
			std::cout << ' ' << target_seq[ mapping[i]-1 ];
		}
		std::cout << std::endl;
	}

	protocols::loops::apply_sequence_mapping( pose, target_seq, mapping );


	dump_pdb( pose, "start.pdb" );



	Loops loops;

	{ // figure out where the loops are
		id::SequenceMapping m( mapping );
		m.reverse();

		// debug:
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			std::cout << "M2: " << I(4,i) << ' ' << pose.residue(i).name() << I(4,m[i]) << std::endl;
		}

		bool in_loop( false );
		Size loop_begin( 0 );

		for ( Size i=1; i<= m.size1(); ++i ) {
			if ( m[i] == 0 ) {
				if ( !in_loop ) {
					loop_begin = i;
					in_loop = true;
				}
			} else {
				if ( in_loop ) {
					loops.add_loop( loop_begin, i-1, i-1, 0.0, 0, true );
					in_loop = false;
				}
			}
		}

		//loops.add_loop(  82,  89,  84, 0, true );
		//loops.add_loop( 100, 107, 103, 0, true );
	}

	/**
		 fragment setup:
		 1. pick 3mer library from vall
		 2. derive 1mer libary from 3mers


	**/
	utility::vector1<int> frag_sizes; //( option[ OptionKeys::loops::frag_sizes ] );



	frag_sizes.push_back(3);
	frag_sizes.push_back(1);

	//FileVectorOption frag_files( option[OptionKeys::loops::frag_files] );
	//assert( frag_sizes.size() == frag_files.size() );

	std::map< Size, protocols::frags::TorsionFragmentLibraryOP > frag_libs;
	std::map< Size, bool > frag_libs_init;


	{ // scope
		using namespace protocols::frags;
		{ // 3mers
			Size const frag_size( 3 );
			protocols::frags::TorsionFragmentLibraryOP frag_lib_op( new protocols::frags::TorsionFragmentLibrary );
			frag_libs_init.insert( std::make_pair( frag_size , false ) );
			frag_libs.insert( std::make_pair(frag_size, frag_lib_op) );
			// gets 3mers from vall:
			setup_loops_fragment_libraries( pose.sequence(), loops, frag_libs, frag_libs_init );
		}

		{ // 1mers
			// gets 1mers from 3mers
			Size const frag_size = 1;
			Size const prev_size = 3;
			protocols::frags::TorsionFragmentLibraryOP frag_lib_op( new protocols::frags::TorsionFragmentLibrary );
			frag_libs_init.insert( std::make_pair( frag_size , false ) );
			frag_libs.insert( std::make_pair(frag_size, frag_lib_op) );

			protocols::frags::TorsionFragmentLibraryOP    prev_lib_op( frag_libs.find( prev_size)->second );

			frag_libs_init[ frag_size ] = frag_lib_op->derive_from_src_lib( frag_size, prev_size, prev_lib_op );
		}

		// confirm success
		std::cout << frag_libs_init[1] << frag_libs_init[3] << std::endl;
		assert( frag_libs_init[1] && frag_libs_init[3] );
	}

	// for graphics:
	protocols::viewer::add_conformation_viewer( pose.conformation(), "loops_pose" );

	basic::prof_reset();

	perturb_loops_with_ccd( pose, loops, frag_libs );

	basic::prof_show();

	{ // switch to fullatom mode, preserve DNA sidechain conformations
		Pose cen_pose;
		cen_pose = pose;

		core::util::switch_to_residue_type_set( pose, FA_STANDARD );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			Residue const &     rsd(     pose.residue(i) );
			Residue const & src_rsd( cen_pose.residue(i) );
			if ( rsd.is_DNA() ) {
				assert( rsd.natoms() == src_rsd.natoms() );
				for ( Size j=1; j<= rsd.natoms(); ++j ) {
					if ( rsd.atom_is_backbone(j) ) assert( rsd.xyz(j).distance_squared( src_rsd.xyz( rsd.atom_name(j) ) )<1e-2 );
					pose.set_xyz( AtomID(j,i), src_rsd.xyz( rsd.atom_name(j) ) );
				}
			}
		}
	}


	basic::prof_reset();
	refine_loops_with_ccd( pose, loops );

	basic::prof_show();

	dump_pdb( pose, "final.pdb" );

}


///////////////////////////////////////////////////////////////////////////////
void
loop_rebuilding_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace scoring;
	using namespace scoring::constraints;
	using namespace conformation;
	using namespace chemical;
	using namespace optimization;
	using namespace id;
	using namespace io::pdb;
	using namespace pose;
	using namespace protocols::loops;
	using namespace protocols::frags;


	// read the pose
	Pose pose;
	core::import_pose::centroid_pose_from_pdb( pose, basic::options::start_file() );

	// setup Chu's loops object
	Loops loops;

 	if (! loops.read_file( option[ OptionKeys::loops::loop_file ]().name() ) ) {

		basic::Error() << "ERROR -- loops_main: can not retrieve loops info\n"
											<< "exit ......\n";
 		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
 	}

 	std::cout << "************************************************\n"
						<< "loops to be modeled:\n" << loops
						<< "************************************************"
						<< '\n';


 	// need to set up fold_tree and allow_bb_move properly ?

	//////////////////////////////////////////////////////////////////////////////////////////////////
	/**
		 fragment setup:
		 1. pick 3mer library from vall
		 2. derive 1mer libary from 3mers
	**/
	utility::vector1<int> frag_sizes; //( option[ OptionKeys::loops::frag_sizes ] );

	frag_sizes.push_back(3);
	frag_sizes.push_back(1);

	std::map< Size, protocols::frags::TorsionFragmentLibraryOP > frag_libs;
	std::map< Size, bool > frag_libs_init;
	{ // scope
		using namespace protocols::frags;
		{ // 3mers
			Size const frag_size( 3 );
			protocols::frags::TorsionFragmentLibraryOP frag_lib_op( new protocols::frags::TorsionFragmentLibrary );
			frag_libs_init.insert( std::make_pair( frag_size , false ) );
			frag_libs.insert( std::make_pair(frag_size, frag_lib_op) );
			// gets 3mers from vall:
			setup_loops_fragment_libraries( pose.sequence(), loops, frag_libs, frag_libs_init );
		}

		{ // 1mers
			// gets 1mers from 3mers
			Size const frag_size = 1;
			Size const prev_size = 3;
			protocols::frags::TorsionFragmentLibraryOP frag_lib_op( new protocols::frags::TorsionFragmentLibrary );
			frag_libs_init.insert( std::make_pair( frag_size , false ) );
			frag_libs.insert( std::make_pair(frag_size, frag_lib_op) );

			protocols::frags::TorsionFragmentLibraryOP    prev_lib_op( frag_libs.find( prev_size)->second );

			frag_libs_init[ frag_size ] = frag_lib_op->derive_from_src_lib( frag_size, prev_size, prev_lib_op );
		}

		// confirm success
		std::cout << frag_libs_init[1] << frag_libs_init[3] << std::endl;
		assert( frag_libs_init[1] && frag_libs_init[3] );
	}

	// for graphics:
	//protocols::viewer::add_conformation_viewer( pose.conformation(), "loops_pose" );

	basic::prof_reset();

	// create scorefunction here


	//perturb_loops_with_ccd( pose, loops, frag_libs, scorefxn );

	basic::prof_show();

	dump_pdb( pose, "final.pdb" );

}


///////////////////////////////////////////////////////////////////////////////

void
dna_dr_test()
{
	using namespace pose;
	using namespace scoring;

	Pose pose;
	core::import_pose::pose_from_pdb( pose, "input/1a3q.pdb");
	pose.dump_pdb( "tmp.pdb" );

	ScoreFunction scorefxn;
	scorefxn.set_weight( dna_dr, 1.0 );

	tt << scorefxn(pose ) << '\n';


}

///////////////////////////////////////////////////////////////////////////////
void
dna_dr_loop_test()
{
  // > test.log 2> err&
  using namespace basic::options::OptionKeys::dna::specificity;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace scoring;
  using namespace scoring::constraints;
  using namespace conformation;
  using namespace chemical;
  using namespace optimization;
  using namespace id;
  using namespace io::pdb;
  using namespace pose;
  using namespace protocols::loops;
  using namespace protocols::frags;

  Pose pose;
	//core::import_pose::pose_from_pdb( pose, "input/1aay.pdb");
	core::import_pose::centroid_pose_from_pdb( pose, basic::options::start_file() );
  //pose.dump_pdb( "tmp.pdb" );

	// read the loop file
	// setup Chu's loops object
	Loops loops;

	if (! loops.read_file( option[ OptionKeys::loops::loop_file ]().name() ) ) {

		basic::Error() << "ERROR -- loops_main: can not retrieve loops info\n"
												<< "exit ......\n";
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}

	std::cout << "************************************************\n"
						<< "loops to be modeled:\n" << loops
						<< "************************************************"
						<< '\n';

	{ // setup a foldtree for the pose
		kinematics::FoldTree f( pose.fold_tree() );

		for( Loops::const_iterator it=loops.begin(), it_end=loops.end(); it != it_end; ++it ) {
			Size const loop_begin( it->start() );
			Size const loop_end  ( it->stop () );
			Size const cutpoint  ( it->cut  () );

			if ( loop_begin > 1 && loop_end < pose.total_residue() &&
					 !pose.residue( loop_begin ).is_terminus() &&
					 !pose.residue( loop_end   ).is_terminus() ) {

				f.new_jump( loop_begin-1, loop_end+1, cutpoint );
			}
		}

		Size dna_pos(1);
		for ( ; dna_pos <= pose.total_residue(); ++dna_pos ) if ( f.is_jump_point( dna_pos ) && pose.residue(dna_pos).is_DNA() ) break;
		assert( dna_pos <= pose.total_residue() );

		// this will tend to keep the DNA fixed in absolute orientation
		f.reorder( dna_pos );

		pose.fold_tree( f );
		std::cout << "starting foldtree: " << f << std::endl;

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////
	//fragment setup: 1. pick 3mer library from vall 2. derive 1mer libary from 3mers
	utility::vector1<int> frag_sizes; //( option[ OptionKeys::loops::frag_sizes ] );

	frag_sizes.push_back(3);
	frag_sizes.push_back(1);

	std::map< Size, protocols::frags::TorsionFragmentLibraryOP > frag_libs;
	std::map< Size, bool > frag_libs_init;
	{ // scope
		using namespace protocols::frags;
		{ // 3mers
			Size const frag_size( 3 );
			protocols::frags::TorsionFragmentLibraryOP frag_lib_op( new protocols::frags::TorsionFragmentLibrary );
			frag_libs_init.insert( std::make_pair( frag_size , false ) );
			frag_libs.insert( std::make_pair(frag_size, frag_lib_op) );
			// gets 3mers from vall:
			setup_loops_fragment_libraries( pose.sequence(), loops, frag_libs, frag_libs_init );
		}

		{ // 1mers
			// gets 1mers from 3mers
			Size const frag_size = 1;
			Size const prev_size = 3;
			protocols::frags::TorsionFragmentLibraryOP frag_lib_op( new protocols::frags::TorsionFragmentLibrary );
			frag_libs_init.insert( std::make_pair( frag_size , false ) );
			frag_libs.insert( std::make_pair(frag_size, frag_lib_op) );

			protocols::frags::TorsionFragmentLibraryOP    prev_lib_op( frag_libs.find( prev_size)->second );

			frag_libs_init[ frag_size ] = frag_lib_op->derive_from_src_lib( frag_size, prev_size, prev_lib_op );
		}

		// confirm success
		std::cout << frag_libs_init[1] << frag_libs_init[3] << std::endl;
		assert( frag_libs_init[1] && frag_libs_init[3] );
	}

	// for graphics:
	protocols::viewer::add_conformation_viewer( pose.conformation(), "loops_pose" );

  Pose const start_pose( pose );

	// setup scoring function
  //
	//utility::vector1< std::string > weights_files;
  //if ( option[ weights_tag ]() != "none" ) {
  //  weights_files.push_back( option[ weights_tag ] );
  //} else {
  //  read_list_file( option[ weights_tag_list ], weights_files );
	//}
  //for ( Size j=1; j<= weights_files.size(); ++j ) {

	for ( int k=-10; k<=4; k++ ) {
			pose = start_pose;
			// need to set up fold_tree and allow_bb_move properly ?
			basic::prof_reset();
			std::string function_tag("cen_std"), patch_tag("score4L");
			ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function(function_tag, patch_tag) );
			//ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( weights_files[j] ) );

			Real wt = pow ( 2, k );
			// dna specific mods
			scorefxn->set_weight( dna_dr, wt );
			// 			scorefxn->set_weight( fa_pair, 0.0 );
			// 			scorefxn->set_weight( fa_elec, option[ Wfa_elec ] );
			// 			scorefxn->set_weight( dna_bp, option[ Wdna_bp ] );
			// 			scorefxn->set_weight( dna_bs, option[ Wdna_bs ] );
			(*scorefxn)( pose );
			scorefxn->show( std::cout, pose );
			for ( int j=0; j< option[ out::nstruct ]; j++ ) // change to j=100 for wallaby, j=25 for running on cluster
			{
				pose = start_pose;
				perturb_loops_with_ccd( pose, loops, frag_libs, scorefxn );
				//perturb_loops_with_ccd( pose, loops, frag_libs );

				basic::prof_show();

				std::string const outfilename
					( option[ output_tag] + "_" +string_of( k )+"_"+lead_zero_string_of( j, 4 )+".pdb" );

				// rmsd calculation between initial pose and final pose
				// iterate over all loops, calculate and print rmsd
				for( Loops::const_iterator it=loops.begin(), it_end=loops.end(); it != it_end; ++it )
				{
					Real loop_rms( 0.0 );
					for ( Size ii=it->start(); ii<= it->stop(); ++ii )
					{
							Residue const & rsd1( start_pose.residue(ii) );
							Residue const & rsd2( pose.residue(ii) );
							loop_rms = loop_rms + rsd1.xyz("CA").distance_squared( rsd2.xyz("CA") );
							//cout << ii << " " << loop_rms << endl;
					}
					loop_rms = std::sqrt( loop_rms/ ( it->stop() - it->start() + 1 ) );
					cout << "**loop rmsd: " << k << ' ' << j << ' ' << outfilename << ' ' << (*scorefxn)(pose) << ' ' << loop_rms << ' ' <<
						pose.energies().total_energies()[ dna_dr ] << std::endl;
				}
				pose.dump_pdb( outfilename );
			}
			//std::cout << k << std::endl;
			// dump_pdb( pose, "final.pdb" );
			// add full atom stuff from loop_modeling_test() ?
			//tt << scorefxn(pose ) << '\n';
	}
}

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string const mode( option[ dna::specificity::mode ].value() );
	if ( mode == "kono_sarai_stats" ) {
		kono_sarai_stats();
		exit(0);
	} else if ( mode == "kono_sarai_zscore" ) {
		kono_sarai_zscore();
		exit(0);
	} else if ( mode == "dna_dr_test" ) {
		dna_dr_test();
		exit(0);
  } else if ( mode == "dna_dr_loop_test" ) {
    dna_dr_loop_test();
    exit(0);
	}


	exit(0); // add new mode strings

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
	// initialize option and random number system
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl; 
	} 
}

