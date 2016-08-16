// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/phil/test1.cc
/// @brief  Some simple examples of how to use basic functionality + some DNA functionality
/// @author Phil Bradley (pbradley@fhcrc.org)

// Project headers
#include <devel/init.hh>

#include <protocols/frags/VallData.hh>
#include <protocols/frags/TorsionFragment.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/OutputMovers.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/loop_closure/ccd/RamaCheck.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/viewer/visualize.hh>

#include <core/init/init.hh>
#include <core/id/AtomID_Map.hh>
#include <core/types.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <core/chemical/AtomType.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/GenBornPotential.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/methods/Methods.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/sequence/DerivedSequenceMapping.hh>
#include <core/sequence/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/datacache/CacheableDataType.hh> // profiling
#include <core/import_pose/import_pose.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/rtmin.hh>
#include <core/pack/min_pack.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

// Basic headers
#include <basic/options/util.hh>
#include <basic/prof.hh> // profiling
#include <basic/datacache/BasicDataCache.hh>
#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/phil.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// Utility headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

// External headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>


////////////////////////////////////////////////
// danger USING ////////////////////////////////
using namespace core;
using namespace protocols;
using namespace id;
using namespace pose;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace basic::options;
using namespace optimization;
using utility::vector1;
using std::string;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;
using core::import_pose::pose_from_file;
using io::pdb::dump_pdb; // deprecated though
static THREAD_LOCAL basic::Tracer tt( "demo.phil.test1", basic::t_trace );

///////////////////////////////////////////////////////////////////////////////
std::ostream &
operator<< ( std::ostream & out, Vector const & v ) {
	out << "( " << F(9,3,v(1)) << " , " << F(9,3,v(2)) << " , " << F(9,3,v(3)) << " )";
	return out;
}

///////////////////////////////////////////////////////////////////////////////
// some silly helper routines:
//
void
read_list_file( std::string const & filename, utility::vector1< std::string > & l )
{
	std::ifstream data( filename.c_str() );
	if ( !data.good() ) {
		utility_exit_with_message( "cant open list file: "+filename );
	}
	std::string line;
	while ( getline( data, line ) ) {
		std::cout << line << std::endl;
		l.push_back( line );
	}
}

///////////////////////////////////////////////////////////////////////////////
void
setup_atom_number(
	pose::Pose const & pose,
	id::AtomID_Map< int > & atom_number
)
{
	core::pose::initialize_atomid_map( atom_number, pose );

	int number(0);
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		for ( Size j=1; j<= Size(pose.residue(i).natoms()); ++j ) {
			atom_number[ id::AtomID( j, i ) ] = ++number;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
void
show_rasmol_hbonds(
	scoring::hbonds::HBondSet const & hbond_set,
	id::AtomID_Map< int > const & atom_number,
	std::ostream & out
)
{
	using namespace scoring::hbonds;

	for ( Size i=1; i<= Size(hbond_set.nhbonds()); ++i ) {
		if ( hbond_set.allow_hbond(i) ) {
			HBond const & hb( hbond_set.hbond(i) );
			id::AtomID const hatm( hb.don_hatm(), hb.don_res() );
			id::AtomID const aatm( hb.acc_atm() , hb.acc_res() );
			out << "monitor " << atom_number[ hatm ] << ' ' << atom_number[ aatm ] <<
				'\n';
		}
	}

}

///////////////////////////////////////////////////////////////////////////////
void
dump_hbond_pdb(
	pose::Pose & pose,
	std::string const & filename
)
{
	using namespace optimization;
	using id::AtomID;
	using id::DOF_ID;
	using id::PHI;
	using id::THETA;
	using id::D;


	scoring::hbonds::HBondSet hbond_set;
	pose.update_residue_neighbors();
	scoring::hbonds::fill_hbond_set( pose, false, hbond_set );

	// setup the global atom numbering that would be used for pdb output
	id::AtomID_Map< int > atom_number;
	setup_atom_number( pose, atom_number );

	std::ofstream out( filename.c_str() );
	out << "load pdb inline\n";
	show_rasmol_hbonds( hbond_set, atom_number, out );
	out << "exit\n";
	dump_pdb( pose, out );
	out.close();
}


///////////////////////////////////////////////////////////////////////////////
void
show_intrachain_energies(
	pose::Pose & pose,
	scoring::ScoreFunction const & scorefxn
)
{
	using namespace scoring;
	using namespace chemical;

	scorefxn( pose );

	//Size const nres( pose.total_residue() );
	Size const nchain( pose.conformation().num_chains() );

	FArray2D< EnergyMap > chain_energies( nchain, nchain );

	// cached energies object
	Energies & energies( pose.energies() );

	// the neighbor/energy links
	EnergyGraph & energy_graph( energies.energy_graph() );

	for ( Size i=1, i_end = pose.total_residue(); i<= i_end; ++i ) {
		conformation::Residue const & resl( pose.residue( i ) );
		for ( graph::Graph::EdgeListIter
				iru  = energy_graph.get_node(i)->upper_edge_list_begin(),
				irue = energy_graph.get_node(i)->upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge * edge( static_cast< EnergyEdge *> (*iru) );
			Size const j( edge->get_second_node_ind() );
			conformation::Residue const & resu( pose.residue( j ) );

			edge->add_to_energy_map( chain_energies( resl.chain(), resu.chain() ) );
		}
	}

	for ( Size i=1; i<= nchain; ++i ) {
		for ( Size j=i; j<= nchain; ++j ) {
			std::cout << "chainE:  " << i << ' ' << j;
			for ( Size k=1; k<= n_score_types; ++k ) {
				if ( scorefxn[ ScoreType(k) ] ) {
					std::cout << ' ' << ScoreType(k) << ' ' << F(9,3,chain_energies(i,j)[ ScoreType(k) ]);
				}
			}
			std::cout << std::endl;
		}
	}

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
test_dunbrack_io()
{
	//using namespace core;
	using namespace scoring;
	using namespace pack::dunbrack;

	RotamerLibrary const & rot_lib( *RotamerLibrary::get_instance() );

	pose::Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

	RotamerLibraryScratchSpace scratch;

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		rot_lib.rotamer_energy( pose.residue(i), scratch );
	}
}

///////////////////////////////////////////////////////////////////////////////
void
test_gb()
{
	//using namespace core;
	using namespace scoring;
	using namespace pose;
	using namespace conformation;
	using namespace basic::options;
	using namespace core::pose::datacache;

	pose::Pose pose;
	core::import_pose::pose_from_file( pose, start_file() , core::import_pose::PDB_file);

	ScoreFunction scorefxn;
	scorefxn.set_weight( gb_elec, 1.0 );

	scorefxn( pose );


	// show radii
	GenBornPoseInfo const & gb_info( static_cast< GenBornPoseInfo & >( pose.data().get( CacheableDataType::GEN_BORN_POSE_INFO ) ) );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		GenBornResidueInfo const & gb( gb_info.residue_info( i ) );
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			std::cout << "GBINFO: " << I(4,i) << ' ' << rsd.name3() << ' ' << rsd.atom_name(j) <<
				F(9,2,gb.born_radius(j) ) << F(9,2,gb.atomic_radius(j) ) << std::endl;
		}
	}

	LREnergyContainerCOP lrec = pose.energies().long_range_container( methods::gen_born_lr );
	assert ( !lrec->empty() );

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		std::cout << "GBPAIR: " << ii << ' ' << ii  << ' ' <<
			F(9,2,pose.energies().onebody_energies( ii )[ gb_elec ] ) << std::endl;

		for ( ResidueNeighborConstIteratorOP
				rni = lrec->const_upper_neighbor_iterator_begin( ii ),
				rniend = lrec->const_upper_neighbor_iterator_end( ii );
				(*rni) != (*rniend); ++(*rni) ) {
			Size const jj = rni->upper_neighbor_id();
			EnergyMap emap;
			rni->retrieve_energy( emap );
			std::cout << "GBPAIR: " << ii << ' ' << jj  << ' ' << F(9,2,emap[gb_elec]) << std::endl;
		}
	}


	dump_pdb( pose, "new_gb.pdb" );

	exit(0);

	/// now pack
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));

	task->initialize_from_command_line();
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( ii%2 == 0 ) task->nonconst_residue_task( ii ).restrict_to_repacking();
		else task->nonconst_residue_task( ii ).prevent_repacking();
	}

	pack::pack_rotamers( pose, scorefxn, task);

}

///////////////////////////////////////////////////////////////////////////////
void
test_rama()
{
	using namespace scoring;
	std::cout << "test rama" << std::endl;

	pose::Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

	ScoreFunction scorefxn;
	scorefxn.set_weight( rama , 0.2 );

	scorefxn( pose );

	std::cout << "rama succeeded" << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
void
test_scorefxn_io()
{
	using namespace scoring;
	pose::Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

	// standard packer wts
	ScoreFunctionOP scorefxn( get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS ) );
	scorefxn->set_weight( fa_rep, 0.0 );
	(*scorefxn)(pose);

	std::cout << *scorefxn;

	/// soft rep packer wts
	ScoreFunctionOP scorefxn2( ScoreFunctionFactory::create_score_function( SOFT_REP_WTS ) );
	(*scorefxn2)(pose);

	std::cout << *scorefxn2;

	/// scorefxn w/ std packer wts
	ScoreFunctionOP score12( get_score_function() );
	(*score12)(pose);

	std::cout << *score12;
}

///////////////////////////////////////////////////////////////////////////////
void
simple_rotamer_test()
{
	using namespace conformation;

	pose::Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

	ResidueOP ile_rsd( pose.residue(3).clone() );
	ResidueOP pro_rsd( pose.residue(70).clone() );

	assert( pro_rsd->name() == "PRO" && ile_rsd->name() == "ILE" );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		if ( rsd.is_terminus() ) continue;

		if ( rsd.nchi() >= 1 ) {
			Residue rot( rsd.type(), rsd, pose.conformation() );
			rot.set_chi( 1, 60.0 );

			Residue rot2( rsd.type(), *ile_rsd, pose.conformation() );
			rot2.set_chi( 1, 60.0 );

			// this fails b/c PRO has no backbone H
			Residue rot3( rsd.type(), *pro_rsd, pose.conformation() );
			rot3.set_chi( 1, 60.0 );
		}

		if ( true ) { //rsd.name() != "PRO" && rsd.name() != "GLY" ) {
			Residue rot( ile_rsd->type(), rsd, pose.conformation() );
			rot.set_chi(1,180.0);
			rot.set_chi(2,180.0);
			ResidueOP new_rsd( rot.create_residue() );

			pose.replace_residue( i, *new_rsd, false );

		}
	}

	dump_pdb( pose, "test_out.pdb" );

	//exit(0);
}

///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////


void
set_fullatom_flag_test()
{
	using namespace io::pdb;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace kinematics;
	using namespace basic::options;

	using namespace protocols::frags;
	using namespace protocols::moves;

	{ // if we dont already have a fullatom pose
		pose::Pose pose;
		core::import_pose::centroid_pose_from_pdb( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

		// now convert to fullatom:
		ResidueTypeSetCOP fullatom_residue_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			Residue const & rsd( pose.residue(i) );
			utility::vector1< std::string > const & variant_types( rsd.type().properties().get_list_of_variants() );
			for ( Size j=1; j<=variant_types.size(); ++j ) {
				std::cout << variant_types[j] << std::endl;
			}
			// get residue types with same AA and variant
			ResidueTypeCOP new_rsd_type( fullatom_residue_set->get_representative_type_aa( rsd.aa(), rsd.type().variant_types() ) );
			ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type, rsd, pose.conformation() ) );
			assert( new_rsd ); // really should print error msg

			pose.replace_residue( i, *new_rsd, false );
		}

		dump_pdb( pose, "fullatom_test1.pdb" );
	} // scope


	/// other way: if we already have a fullatom_pose
	{

		pose::Pose fullatom_pose;
		core::import_pose::pose_from_file( fullatom_pose, "input/test_in.pdb" , core::import_pose::PDB_file);

		pose::Pose pose;
		core::import_pose::centroid_pose_from_pdb( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

		// now convert to fullatom:
		// only tricky thing -- cutpoint residue types might not be preserved...
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			pose.replace_residue( i, fullatom_pose.residue(i), true /* orient_backbone */ );
		}

		dump_pdb( pose, "fullatom_test2.pdb" );

	} // scope

}

///////////////////////////////////////////////////////////////////////////////
void
delete_test()
{

	Pose start_pose;
	core::import_pose::pose_from_file( start_pose, "input/1aay.pdb" , core::import_pose::PDB_file);

	Size const nres( start_pose.total_residue() );

	Pose pose( start_pose );

	// standard packer wts
	ScoreFunctionOP scorefxn( get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS ) );

	// from beginning
	for ( Size i=1; i< nres; ++i ) {
		pose.conformation().delete_residue_slow( 1 );
		std::cout << "score after delete: " << F(9,3,(*scorefxn)(pose)) << ' ' << pose.fold_tree() << std::endl;
	}

	// from end
	pose = start_pose;
	for ( Size i=1; i< nres; ++i ) {
		pose.conformation().delete_residue_slow( pose.total_residue() );
		std::cout << "score after delete: " << F(9,3,(*scorefxn)(pose)) << ' ' << pose.fold_tree() << std::endl;
	}


	// from 2
	pose = start_pose;
	for ( Size i=1; i< nres; ++i ) {
		pose.conformation().delete_residue_slow( 2 );
		std::cout << "score after delete: " << F(9,3,(*scorefxn)(pose)) << ' ' << pose.fold_tree() << std::endl;
	}

	// residue_range
	pose = start_pose;
	pose.conformation().delete_residue_range_slow( 2, nres-1 );
	std::cout << "score after delete: " << F(9,3,(*scorefxn)(pose)) << ' ' << pose.fold_tree() << std::endl;
}


void
simple_loop_modeling_test()
{
	using namespace io::pdb;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace kinematics;
	using namespace id;
	using namespace basic::options;

	using namespace protocols::frags;
	using namespace protocols::moves;

	using core::Real;
	Real const init_phi  ( -150.0 );
	Real const init_psi  (  150.0 );
	Real const init_omega(  180.0 );

	// get rid of this crap when the MonteCarlo changes are reverted
	pose::PoseOP pose_ptr( new pose::Pose() );
	pose::Pose & pose( *pose_ptr );
	core::import_pose::centroid_pose_from_pdb( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

	std::cout << pose.sequence() << std::endl;

	core::sequence::DerivedSequenceMapping mapping;
	std::string source_seq, target_seq;
	read_alignment_file( start_file(), source_seq, target_seq, mapping );

	core::sequence::DerivedSequenceMapping start_mapping( mapping );

	assert( mapping.size1() == source_seq.size() && mapping.size2() == target_seq.size() );

	std::cout << "start mapping: " << std::endl;
	mapping.show();


	// alternate approach:

	// first delete all unaligned residues
	while ( !mapping.all_aligned() ) {
		for ( Size i=1; i<= mapping.size1(); ++i ) {
			if ( !mapping[i] ) {
				pose.conformation().delete_polymer_residue( i );
				mapping.delete_source_residue( i );
				break; // because the numbering is screwed up now, start the loop again
			}
		}
	}

	// now convert sequence of aligned positions
	ResidueTypeSet const & rsd_set( *pose.residue(1).residue_type_set() );
	{
		for ( Size i=1; i<= mapping.size1(); ++i ) {
			char const new_seq( target_seq[ mapping[i]-1 ] ); // strings are 0-indexed
			// Representative type should have no/minimal variants
			ResidueTypeCOP new_rsd_type( rsd_set.get_representative_type_name1( new_seq ) );
			ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type, pose.residue(i), pose.conformation() ) );
			pose.replace_residue( i, *new_rsd, false );
		}
	}

	// setup fold_tree
	{
		assert( pose.total_residue() == mapping.size1() );
		kinematics::FoldTree f( pose.fold_tree() );
		for ( Size i=1; i< pose.total_residue(); ++i ) {
			assert( mapping[i] );
			if ( mapping[i+1] != mapping[i]+1 ) {
				f.new_jump( i, i+1, i );
			}
		}
		pose.fold_tree(f);
		tt << "FOldtree: " << f << std::endl;
	}

	// add terminal residues
	// nterm
	while ( mapping[ 1 ] != 1 ) {
		int const aligned_pos( mapping[1] - 1 );
		char const new_seq( target_seq[ aligned_pos-1 ] ); // 0-indexed
		// Representative type should have no/minimal variants
		ResidueTypeCOP new_rsd_type( rsd_set.get_representative_type_name1( new_seq ) );
		ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type ) );
		pose.conformation().prepend_polymer_residue_before_seqpos( *new_rsd, 1, true );
		pose.set_omega( 1, 180.0 );
		mapping.insert_aligned_residue( 1, aligned_pos );
	}
	// cterm
	while ( mapping[ mapping.size1() ] != mapping.size2() ) {
		int const seqpos( mapping.size1() + 1 );
		int const aligned_pos( mapping[seqpos-1] + 1 );
		char const new_seq( target_seq[ aligned_pos-1 ] ); // 0-indexed
		// Representative type should have no/minimal variants
		ResidueTypeCOP new_rsd_type( rsd_set.get_representative_type_name1( new_seq ) );
		ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type ) );
		pose.conformation().append_polymer_residue_after_seqpos( *new_rsd, seqpos-1, true );
		pose.set_omega( seqpos-1, 180.0 );
		mapping.insert_aligned_residue( seqpos, aligned_pos );
	}

	// now fill in the breaks
	{
		for ( Size cut=1; cut<= Size(pose.fold_tree().num_cutpoint()); ++cut ) {
			while ( mapping[ pose.fold_tree().cutpoint( cut )+1] != Size(pose.fold_tree().cutpoint( cut )+1 ) ) {
				Size const cutpoint( pose.fold_tree().cutpoint( cut ) );
				assert( mapping[cutpoint] == cutpoint ); // we've fixed everything up til here

				// add at the nterm of the loop
				int const aligned_pos( cutpoint+1 );
				char const new_seq( target_seq[ aligned_pos - 1 ] ); // 0-indexed
				// Representative type should have no/minimal variants
				ResidueTypeCOP new_rsd_type( rsd_set.get_representative_type_name1( new_seq ) );
				ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type ) );
				pose.conformation().append_polymer_residue_after_seqpos( *new_rsd, cutpoint, true );
				pose.set_omega( cutpoint, 180.0 );
				mapping.insert_aligned_residue( cutpoint + 1, aligned_pos );
			} // while missing residues
			if ( true ) {
				int const cutpoint( pose.fold_tree().cutpoint( cut ) );
				core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, cutpoint   );
				core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, cutpoint+1 );
			}
		} // cut=1,ncutpoints
	} // scope

	assert( mapping.is_identity() );
	assert( pose.sequence() == target_seq );

	dump_pdb( pose, "loops_start.pdb" );

	// setup allow move
	Size const nres( pose.total_residue() );
	MoveMap mm;
	start_mapping.reverse();
	assert( start_mapping.size1() == nres );
	for ( Size i=1; i<= nres; ++i ) {
		if ( ( !start_mapping[i] ) || // unaligned
				( i>   1 && start_mapping[i] != start_mapping[i-1]+1 ) ||
				( i<nres && start_mapping[i] != start_mapping[i+1]-1 ) ) {
			std::cout << "moving: " << i << ' ' << start_mapping[i] << std::endl;
			mm.set_bb( i, true );
			if ( !start_mapping[i] ) {
				pose.set_phi  ( i, init_phi   );
				pose.set_psi  ( i, init_psi   );
				pose.set_omega( i, init_omega );
			}
		}
	}


	////////////////////////
	// pick fragments
	// read vall:
	VallData vall( option[ OptionKeys::phil::vall_file ] );

	TorsionFragmentLibrary lib;
	Size const frag_size(3);

	lib.resize( nres - frag_size + 1 );
	for ( Size i=1; i<= nres-frag_size+1; ++i ) {
		std::string const frag_seq( target_seq.substr(i-1,3) );
		utility::vector1< Size > homs_to_exclude;
		vall.get_frags( 200, frag_seq, "---", 1.0, 0.0, false, false, true, homs_to_exclude, lib[i] );
	}


	// setup scoring function
	// get rid of this crap when MonteCarlo changes are reverted
	core::scoring::ScoreFunctionOP scorefxn_ptr( new ScoreFunction() );
	ScoreFunction & scorefxn( *scorefxn_ptr );

	scorefxn.set_weight( scoring::env, 1.0 );
	scorefxn.set_weight( scoring::pair, 1.0 );
	scorefxn.set_weight( scoring::cbeta, 1.0 );
	scorefxn.set_weight( scoring::vdw, 1.0 );
	scorefxn.set_weight( scoring::chainbreak, 1.0 );

	for ( Size i=1; i<= nres; ++i ) {
		std::cout << "natoms " << i << ' ' << pose.residue(i).natoms() << std::endl;
	}

	scorefxn( pose );
	pose.energies().show( std::cout );

	{ // testing
		pose::Pose new_pose;
		new_pose = pose;
		Energies E( pose.energies() );

		std::ofstream out1( "out1" );
		std::ofstream out2( "out2" );
		std::ofstream out3( "out3" );

		pose.energies().show( out1 );
		new_pose.energies().show( out2 );
		E.show( out3 );

		out1.close();
		out2.close();
		out3.close();

	}

	{ // hacking
		FoldTree f( nres );
		f.new_jump( 1, 30, pose.fold_tree().cutpoint(1) );
		pose.fold_tree(f);
	}


	// fragment moves
	MonteCarlo mc( *pose_ptr, *scorefxn_ptr, 2.0 );

	protocols::viewer::add_conformation_viewer( pose.conformation() );
	protocols::viewer::add_monte_carlo_viewer( mc );

	basic::prof_reset();
	for ( int i=1; i<= 10; ++i ) {
		for ( int ii=1; ii<= 100; ++ii ) {
			// choose an insertion position
			int pos;
			while ( true ) {
				pos = static_cast< int > ( ( nres - frag_size + 1 ) * numeric::random::rg().uniform() + 1 );
				bool allowed( true );
				for ( Size k=0; k< frag_size; ++k ) {
					if ( !( mm.get( TorsionID( pos+k, BB,   phi_torsion ) ) &&
							mm.get( TorsionID( pos+k, BB,   psi_torsion ) ) &&
							mm.get( TorsionID( pos+k, BB, omega_torsion ) ) ) ) {
						allowed = false;
						break;
					}
				}
				if ( allowed ) break;
			}

			// choose a fragment
			Size const nfrags( lib[ pos ].size() );
			int const nn( static_cast< int >( nfrags * numeric::random::rg().uniform() + 1 ) );

			// make move
			lib[pos][nn].insert( pose, pos );

			std::cout << "score_before: " << scorefxn( pose ) << std::endl;

			mc.boltzmann( *pose_ptr );

			std::cout << "score_after: " << scorefxn( pose ) << ' ' << mc.last_accepted_score() << ' ' <<
				mc.lowest_score() << std::endl;

		}
		mc.recover_low( *pose_ptr );

		basic::prof_show();
	}

	basic::prof_show();
}

// for graphics
void*
simple_loop_modeling_test_wrapper( void* )
{
	simple_loop_modeling_test();
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
void
dna_deriv_test_old()
{
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace kinematics;
	using namespace optimization;
	using namespace basic::options;
	using namespace pose;

	using namespace io::pdb;
	using namespace protocols::frags;
	using namespace protocols::moves;
	using namespace scoring::dna;

	using core::Real;


	Pose pdb_pose;
	core::import_pose::pose_from_file( pdb_pose, start_file(), core::import_pose::PDB_file); // eg 1qn4
	Size const nres( pdb_pose.total_residue() );

	set_base_partner( pdb_pose );
	BasePartner const & partner( retrieve_base_partner_from_pose( pdb_pose ) );

	// create a mini-pose with just dna
	Pose pose, jump_pose;

	std::string const jump_anchor_atom( " C2 " );
	for ( Size chain1_begin=1; chain1_begin<= nres; ++chain1_begin ) {
		if ( partner[ chain1_begin ] ) {
			// found first dna residue
			int const chain1( pdb_pose.chain( chain1_begin ) );

			for ( Size j=chain1_begin; ( j<= nres && partner[j] && pdb_pose.chain(j) == chain1 ); ++j ) {
				Residue const &         rsd( pdb_pose.residue(            j ) );
				Residue const & partner_rsd( pdb_pose.residue( partner[ j ] ) );

				Size const rsd_seqpos( j - chain1_begin + 1 );
				if ( rsd_seqpos == 1 ) {
					pose.append_residue_by_bond( rsd );
					jump_pose.append_residue_by_bond( rsd );
				} else {
					pose.append_polymer_residue_after_seqpos( rsd, rsd_seqpos-1, false );
					jump_pose.insert_residue_by_jump( rsd, rsd_seqpos, rsd_seqpos-1,
						jump_anchor_atom, // anchor
						jump_anchor_atom ); // root_atomno
				}

				tt << "a " << j << pose.fold_tree() << std::endl;
				pose.insert_residue_by_jump( partner_rsd, rsd_seqpos + 1, rsd_seqpos, jump_anchor_atom,
					jump_anchor_atom );
				jump_pose.insert_residue_by_jump( partner_rsd, rsd_seqpos + 1, rsd_seqpos, jump_anchor_atom,
					jump_anchor_atom );
				tt << "b " << j << pose.fold_tree() << std::endl;
			}
			break;
		}
	}

	dump_pdb( pose, "mini_pose.pdb" );
	tt << "posefinale: " << pose.fold_tree() << std::endl;
	tt << "jump_posefinal: " << jump_pose.fold_tree() << std::endl;
	set_base_partner( pose );
	set_base_partner( jump_pose );

	Pose start_pose, start_jump_pose;
	start_pose = pose;
	start_jump_pose = jump_pose;

	MoveMap mm;
	mm.set_bb( true );
	mm.set_jump( true );


	{
		// setup the options
		MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.01, true /*use_nblist*/ );
		ScoreFunction scorefxn;
		scorefxn.set_weight( dna_bs, 0.5 );
		scorefxn.set_weight( dna_bp, 0.5 );


		AtomTreeMinimizer minimizer;
		minimizer.run( jump_pose, mm, scorefxn, options );
		dump_pdb( jump_pose, "jump_pose_bs_bp.pdb" );

		jump_pose = start_jump_pose;
	}


	{ // dna_bs minimize
		// setup the options
		MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.1, true /*use_nblist*/,
			true /*deriv_check*/, false /*no verbose-deriv-check, is default*/ );
		ScoreFunction scorefxn;
		scorefxn.set_weight( dna_bs, 0.5 );


		AtomTreeMinimizer minimizer;
		std::cout << "MINTEST: dna_bs" << std::endl;
		pose = start_pose;
		minimizer.run( pose, mm, scorefxn, options );
		dump_pdb( pose, "pose_bs.pdb" );

		std::cout << "MINTEST: dna_bs jump_pose" << std::endl;
		minimizer.run( jump_pose, mm, scorefxn, options );
		dump_pdb( jump_pose, "jump_pose_bs.pdb" );
	}
	{
		// dna_bs minimize
		MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.1, true /*use_nblist*/,
			true /*deriv_check*/, false /*no verbose-deriv-check, is default*/ );
		ScoreFunction scorefxn;
		scorefxn.set_weight( dna_bp, 0.5 );


		AtomTreeMinimizer minimizer;
		std::cout << "MINTEST: dna_bp" << std::endl;
		pose = start_pose;
		minimizer.run( pose, mm, scorefxn, options );
		dump_pdb( pose, "pose_bp.pdb" );
	}

}


///////////////////////////////////////////////////////////////////////////////
void
dna_deriv_test()
{
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace kinematics;
	using namespace optimization;
	using namespace basic::options;
	using namespace pose;

	using namespace io::pdb;
	using namespace protocols::frags;
	using namespace protocols::moves;
	using namespace scoring::dna;

	using core::Real;


	Pose pdb_pose;
	core::import_pose::pose_from_file( pdb_pose, "input/1aay.pdb" , core::import_pose::PDB_file); // use 1qn4 for massive dna bend
	Size const nres( pdb_pose.total_residue() );

	set_base_partner( pdb_pose );
	BasePartner const & partner( retrieve_base_partner_from_pose( pdb_pose ) );

	// create a mini-pose with just dna
	Pose pose, jump_pose;

	std::string const jump_anchor_atom( " C2 " );
	for ( Size chain1_begin=1; chain1_begin<= nres; ++chain1_begin ) {
		if ( partner[ chain1_begin ] ) {
			// found first dna residue
			int const chain1( pdb_pose.chain( chain1_begin ) );

			for ( Size j=chain1_begin; ( j<= nres && partner[j] && pdb_pose.chain(j) == chain1 ); ++j ) {
				Residue const &         rsd( pdb_pose.residue(            j ) );
				Residue const & partner_rsd( pdb_pose.residue( partner[ j ] ) );

				Size const rsd_seqpos( j - chain1_begin + 1 );
				if ( rsd_seqpos == 1 ) {
					pose.append_residue_by_bond( rsd );
					jump_pose.append_residue_by_bond( rsd );
				} else {
					pose.append_polymer_residue_after_seqpos( rsd, rsd_seqpos-1, false );
					jump_pose.insert_residue_by_jump( rsd, rsd_seqpos, rsd_seqpos-1,
						jump_anchor_atom, // anchor
						jump_anchor_atom ); // root_atomno
				}

				tt << "a " << j << pose.fold_tree() << std::endl;
				pose.insert_residue_by_jump( partner_rsd, rsd_seqpos + 1, rsd_seqpos, jump_anchor_atom,
					jump_anchor_atom );
				jump_pose.insert_residue_by_jump( partner_rsd, rsd_seqpos + 1, rsd_seqpos, jump_anchor_atom,
					jump_anchor_atom  );
				tt << "b " << j << pose.fold_tree() << std::endl;
			}
			break;
		}
	}

	tt << "posefinale: " << pose.fold_tree() << std::endl;
	tt << "jump_posefinal: " << jump_pose.fold_tree() << std::endl;
	set_base_partner( pose );
	set_base_partner( jump_pose );

	Pose start_pose, start_jump_pose;
	start_pose = pose;
	start_jump_pose = jump_pose;

	MoveMapOP mm( new MoveMap );
	mm->set_bb( true );
	mm->set_jump( true );

	// setup the options
	scoring::ScoreFunctionOP scorefxn( new scoring::ScoreFunction );
	protocols::simple_moves::MinMover min_mover( mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.01, true /*use_nblist*/ );

	{
		scorefxn->set_weight( dna_bs, 0.5 );
		scorefxn->set_weight( dna_bp, 0.5 );

		min_mover.apply( jump_pose );
		jump_pose = start_jump_pose;
	}


	{ // dna_bs minimize
		// setup the options
		min_mover.min_type( "linmin" );
		min_mover.tolerance( 0.1 );
		min_mover.deriv_check( true );

		scorefxn->reset();
		scorefxn->set_weight( dna_bs, 0.5 );

		std::cout << "MINTEST: dna_bs" << std::endl;
		pose = start_pose;
		min_mover.apply( pose );

		std::cout << "MINTEST: dna_bs jump_pose" << std::endl;
		min_mover.apply( jump_pose );
	}
	{
		// dna_bs minimize
		scorefxn->reset();
		scorefxn->set_weight( dna_bp, 0.5 );

		std::cout << "MINTEST: dna_bp" << std::endl;
		pose = start_pose;
		min_mover.apply( pose );
	}

}


///////////////////////////////////////////////////////////////////////////////
void
dna_io_test()
{
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace kinematics;
	using namespace basic::options;
	using namespace pose;

	using namespace io::pdb;
	using namespace protocols::frags;
	using namespace protocols::moves;

	using core::Real;

	// read the filenames
	utility::vector1< std::string > files( start_files() );


	// read the files
	basic::prof_reset();
	for ( Size n=1; n<= files.size(); ++n ) {
		Pose pose;
		std::cout << "Reading file: " << files[n] << std::endl;
		core::import_pose::pose_from_file( pose, files[n] , core::import_pose::PDB_file);
		for ( Size i=1; i<= pose.conformation().num_chains(); ++i ) {
			std::cout << "chain_sequence: " << files[n] << ' ' << i << ' ' << pose.chain_sequence(i) << std::endl;
		}
		basic::prof_show();
	}


}


///////////////////////////////////////////////////////////////////////////////
void
old_simple_conformation_test()
{


	pose::Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

	pose.set_phi( 3, pose.phi(3) );

	dump_pdb( pose, "output/test1.pdb" );

	pose.set_phi( 3, 0.0 );

	dump_pdb( pose, "output/test2.pdb" );

	pose.set_chi( 1, 3, 60.0 );

	dump_pdb( pose, "output/test3.pdb" );

	Size const nres( pose.total_residue() );
	kinematics::FoldTree f( nres );
	f.reorder( nres );

	pose.fold_tree( f );

	pose.set_phi( 3, 120.0 );

	dump_pdb( pose, "output/test4.pdb" );


}

///////////////////////////////////////////////////////////////////////////////
void
simple_centroid_test()
{
	using namespace chemical;

	// read centroid residue set
	ResidueTypeSetCOP rsd_set( ChemicalManager::get_instance()->residue_type_set( "centroid" ) );

	pose::Pose pose;
	//core::import_pose::pose_from_file( pose, *rsd_set, "input/test_in_cen.pdb" , core::import_pose::PDB_file);
	core::import_pose::pose_from_file( pose, *rsd_set, "input/test_in.pdb" , core::import_pose::PDB_file);

	scoring::ScoreFunction scorefxn;
	scorefxn.set_weight( scoring::env, 1.0 );
	scorefxn.set_weight( scoring::pair, 1.0 );
	scorefxn.set_weight( scoring::cbeta, 1.0 );
	scorefxn.set_weight( scoring::vdw, 1.0 );

	scorefxn( pose );

	/// Now handled automatically.  scorefxn.accumulate_residue_total_energies( pose );
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		//// THIS IS THE INCORRECT BEHAVIOR.  VDW IS NOT A ONE BODY ENERGY.
		//// THIS COODE SHOULD INSTEAD READ FROM THE RESIDUE TOTAL ENERGIES ARRAY.
		std::cout << "vdw: " << I(4,i) << F(9,2,pose.energies().onebody_energies(i)[ scoring::vdw ] ) << std::endl;
	}
	std::cout << "vdw_total: " << F(9,2,pose.energies().total_energies()[ scoring::vdw ] ) << std::endl;

	// if the input file was test_in.pdb, we can compare the placement of centroids in test_cen.pdb
	// with those in test_in_cen.pdb, which has rosetta++-placed centroids
	//
	// looks good!
	dump_pdb( pose, "test_cen.pdb" );

}

///////////////////////////////////////////////////////////////////////////////
void
simple_frag_test()
{
	using namespace protocols::frags;

	using namespace basic::options;


	// read vall:
	VallData vall( option[ OptionKeys::phil::vall_file ] );

	{
		SingleResidueTorsionFragmentLibrary lib;
		utility::vector1< Size > homs_to_exclude;
		vall.get_frags( 200, "LLL", "HHH", 1.0, 1.0, true, true, true, homs_to_exclude, lib );
	}


	// read centroid residue set
	chemical::ResidueTypeSetCOP rsd_set( chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" ) );

	pose::Pose pose;
	core::import_pose::pose_from_file( pose, *rsd_set, "input/test_in.pdb" , core::import_pose::PDB_file);

	std::string const sequence( pose.sequence() );

	TorsionFragmentLibrary lib;
	Size const nres( pose.total_residue() );
	Size const frag_size(3);

	lib.resize( nres - frag_size + 1 );
	for ( Size i=1; i<= nres-frag_size+1; ++i ) {
		std::string const frag_seq( sequence.substr(i-1,3) );
		utility::vector1< Size > homs_to_exclude;
		vall.get_frags( 200, frag_seq, "---", 1.0, 0.0, false, false, true, homs_to_exclude, lib[i] );
	}


	{ // try building an ideal peptide:
		using namespace pose;
		using namespace chemical;
		using namespace conformation;

		Pose pose;
		for ( Size i=1; i<= 20; ++i ) {
			utility::vector1< std::string > variants;
			if ( i == 1 ) {
				variants.push_back( "LOWER_TERMINUS_VARIANT" );
			}
			if ( i == 20 ) {
				variants.push_back( "UPPER_TERMINUS_VARIANT" );
			}
			ResidueTypeCOP rsd_type( rsd_set->get_representative_type_aa( static_cast<AA>(i) /*BAD*/, variants ) );
			ResidueOP new_rsd( ResidueFactory::create_residue( *rsd_type ) );
			pose.append_residue_by_bond( *new_rsd, true );
		}
		for ( Size i=1; i<= 20; ++i ) {
			pose.set_omega(i,180.0);
		}
		io::pdb::dump_pdb( pose, "test_ideal.pdb" );
	}
}


///////////////////////////////////////////////////////////////////////////////
void
dna_design_test_old( pose::Pose & pose )
{
	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace io::pdb;

	using namespace core::scoring::dna;

	Size const nres( pose.total_residue() );

	// for bp,bs stats:
	utility::vector1< Size > partner;
	find_basepairs( pose, partner );

	{ // test scoring
		set_base_partner( pose );

		ScoreFunction scorefxn;
		scorefxn.set_weight( dna_bp, 1.0 );
		scorefxn.set_weight( dna_bs, 1.0 );

		scorefxn( pose );
		return;
	}


	{ // test basepair params, basestep params
		utility::vector1< Real > params;
		for ( Size i=1; i<= nres; ++i ) {
			if ( partner[i] ) {
				get_base_pair_params( pose.residue(i), pose.residue( partner[i] ), params );
				if ( i<nres && i+1 != partner[i] && partner[i+1] && partner[i+1] == partner[i]-1 ) {
					assert( !pose.residue(i).is_upper_terminus() );
					get_base_step_params( pose.residue(i), pose.residue( i+1 ), params );
				}
			}
		}
	}

	return;


}

///////////////////////////////////////////////////////////////////////////////
void
dna_coupled_rotamer_design_test()
{

	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace io::pdb;

	using namespace core::scoring::dna;
	using namespace basic::options;

	Pose pose;
	core::import_pose::pose_from_file( pose, "input/1aay.pdb" , core::import_pose::PDB_file);
	Size const nres( pose.total_residue() );

	set_base_partner( pose );
	BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );


	ScoreFunction scorefxn;

	//if ( option[ OptionKeys::phil::soft_rep ] ) scorefxn.set_etable( FA_STANDARD_SOFT );

	scorefxn.set_weight( fa_atr, 0.80 ); // LJ attractive
	scorefxn.set_weight( fa_rep, 0.44 ); // LJ repulsize
	scorefxn.set_weight( fa_sol, 0.65 ); // LK solvation
	//scorefxn.set_weight( fa_pair, 0.49 ); // knowledge-based rsd-rsd electrostatics
	scorefxn.set_weight( fa_dun, 0.56 ); // Dunbrack rotamer energy
	scorefxn.set_weight( rama, 0.2 ); // ramachandran score
	scorefxn.set_weight( hbond_lr_bb, 1.17 ); // long-range backbone hbonds
	scorefxn.set_weight( hbond_sr_bb, 1.17 ); // short-range (helical) backbone hbonds
	scorefxn.set_weight( hbond_bb_sc, 1.17 ); // backbone-sidechain hbonds
	scorefxn.set_weight( hbond_sc   , 1.10 ); // sidechain-sidechain hbonds
	scorefxn.set_weight( fa_elec  , 1.00 ); // hack electrostatics

	show_intrachain_energies( pose, scorefxn );

	{ // try rotamer optimizations //////////////////////////////////////////////////////////////////////////////
		Energy score_orig = scorefxn( pose );

		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		// new interface for packer task: task starts out redesigning at all protein positions with all amino
		// acids, and starts out repacking at all other residues -- the new interface requires the protocol
		// make restrictions from here, providing commutativity of restrictions.
		task->initialize_from_command_line();
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( pose.residue(ii).is_protein() ) {
				task->nonconst_residue_task( ii ).restrict_to_repacking();
				//task->nonconst_residue_task( ii ).prevent_repacking();
				assert( task->pack_residue(ii) );
			} else {
				task->nonconst_residue_task( ii ).allow_aa( na_ade );
				task->nonconst_residue_task( ii ).allow_aa( na_thy );
				task->nonconst_residue_task( ii ).allow_aa( na_gua );
				task->nonconst_residue_task( ii ).allow_aa( na_cyt );
				assert( task->design_residue(ii) );
			}
		}

		{ // setup residue couplings
			using namespace pack::rotamer_set;
			RotamerCouplingsOP couplings( new RotamerCouplings() );
			couplings->resize( nres );
			for ( Size i=1; i<= nres; ++i ) {
				if ( partner[i] ) {
					(*couplings)[i].first = partner[i];
					(*couplings)[i].second = core::conformation::ResidueMatcherOP( new conformation::WatsonCrickResidueMatcher() );
				}
			}
			task->rotamer_couplings( couplings );
		}

		Energy rottrial_score;
		if ( true ) {
			pack::rotamer_trials( pose, scorefxn, task );

			rottrial_score = scorefxn( pose );

			std::cout << "Completed DNA-design rotamer_trials_test() with new score: " << rottrial_score << " vs orig: " <<
				score_orig << std::endl;

			{ // test for mismatches
				WatsonCrickResidueMatcher m;
				for ( Size i=1; i<= nres; ++i ) {
					if ( partner[i]>i ) {
						std::cout << i << pose.residue(i).aa() << ' ' << pose.residue(partner[i]).aa() << std::endl;
						assert( m( pose.residue(i), pose.residue(partner[i])) ); // fails if mismatch
					}
				}
			}
		} else {
			rottrial_score = scorefxn( pose );
		}

		// now try packing
		task->or_include_current( true ); // tmp hack
		pack::pack_rotamers( pose, scorefxn, task);
		Energy pack_score = scorefxn( pose );

		std::cout << "Completed DNA_design pack_rotamers_test() with new score: " << pack_score << " vs orig: " <<
			rottrial_score << std::endl;

		show_intrachain_energies( pose, scorefxn );

		{ // test for mismatches
			WatsonCrickResidueMatcher m;
			for ( Size i=1; i<= nres; ++i ) {
				if ( partner[i]>i ) {
					std::cout << i << pose.residue(i).aa() << ' ' << pose.residue(partner[i]).aa() << std::endl;
					assert( m( pose.residue(i), pose.residue(partner[i])));
				}
			}
		}

		// these are useful but cost a little time to get
		//scorefxn.accumulate_residue_total_energies( pose );
		//pose.energies().show( std::cout );

	}
}

///////////////////////////////////////////////////////////////////////////////
void
dna_design_test()
{
	using namespace io::pdb;

	// read the filenames
	utility::vector1< std::string > files( basic::options::start_files() );

	// read the files
	for ( Size n=1; n<= files.size(); ++n ) {
		pose::Pose pose;
		std::cout << "Reading file: " << files[n] << std::endl;
		core::import_pose::pose_from_file( pose, files[n] , core::import_pose::PDB_file);
		std::cout << "file total_residue: " << pose.total_residue() << std::endl;
		if ( pose.total_residue() ) {
			dna_design_test_old( pose );
		}
	}


}

///////////////////////////////////////////////////////////////////////////////

void
simple_dna_test()
{
	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace io::pdb;

	Pose pose;


	core::import_pose::pose_from_file( pose, /*residue_set, */ "input/1aay.pdb" , core::import_pose::PDB_file);

	std::cout << pose.sequence() << std::endl;

	ScoreFunction scorefxn;
	scorefxn.set_weight( fa_atr, 0.80 ); // LJ attractive
	scorefxn.set_weight( fa_rep, 0.44 ); // LJ repulsize
	scorefxn.set_weight( fa_sol, 0.65 ); // LK solvation
	//scorefxn.set_weight( fa_pair, 0.49 ); // knowledge-based rsd-rsd electrostatics
	scorefxn.set_weight( fa_dun, 0.56 ); // Dunbrack rotamer energy
	scorefxn.set_weight( rama, 0.2 ); // ramachandran score
	scorefxn.set_weight( hbond_lr_bb, 1.17 ); // long-range backbone hbonds
	scorefxn.set_weight( hbond_sr_bb, 1.17 ); // short-range (helical) backbone hbonds
	scorefxn.set_weight( hbond_bb_sc, 1.17 ); // backbone-sidechain hbonds
	scorefxn.set_weight( hbond_sc   , 1.10 ); // sidechain-sidechain hbonds

	{ // try rotamer optimizations //////////////////////////////////////////////////////////////////////////////
		Energy score_orig = scorefxn( pose );

		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		// new interface for packer task: task starts out redesigning at all protein positions with all amino
		// acids, and starts out repacking at all other residues -- the new interface requires the protocol
		// make restrictions from here, providing commutativity of restrictions.
		task->initialize_from_command_line();
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( pose.residue(ii).is_protein() ) {
				task->nonconst_residue_task( ii ).restrict_to_repacking();
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}

		pack::rotamer_trials( pose, scorefxn, task );

		Energy rottrial_score = scorefxn( pose );

		std::cout << "Completed SDT rotamer_trials_test() with new score: " << rottrial_score << " vs orig: " <<
			score_orig << std::endl;

		// now try packing
		task->or_include_current( true ); // tmp hack
		pack::pack_rotamers( pose, scorefxn, task);
		Energy pack_score = scorefxn( pose );

		std::cout << "Completed SDT pack_rotamers_test() with new score: " << pack_score << " vs orig: " <<
			rottrial_score << std::endl;

		// these are useful but cost a little time to get
		//scorefxn.accumulate_residue_total_energies( pose );
		//pose.energies().show( std::cout );

	}


}

///////////////////////////////////////////////////////////////////////////////
void
simple_copy_test()
{

	pose::Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

	pose::Pose new_pose;
	new_pose = pose;

	dump_pdb( new_pose, "test_copy.pdb" );
}


///////////////////////////////////////////////////////////////////////////////
void
simple_hbond_test()
{
	using namespace optimization;
	using id::AtomID;
	using id::DOF_ID;
	using id::PHI;
	using id::THETA;
	using id::D;

	pose::Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

	dump_hbond_pdb( pose, "test_out.hbonds.pdb" );

	// setup scorefxn
	scoring::ScoreFunctionOP scorefxn( new scoring::ScoreFunction );
	scorefxn->set_weight( scoring::hbond_lr_bb, 1.0 );
	scorefxn->set_weight( scoring::hbond_sr_bb, 1.0 );
	scorefxn->set_weight( scoring::hbond_bb_sc, 1.0 );
	scorefxn->set_weight( scoring::hbond_sc, 1.0 );

	(*scorefxn)( pose );

	//pose.energies().show( std::cout );

	// set the moving dofs
	kinematics::MoveMapOP mm1( new kinematics::MoveMap );
	kinematics::MoveMapOP mm2( new kinematics::MoveMap );
	kinematics::MoveMapOP mm3( new kinematics::MoveMap );
	kinematics::MoveMapOP mm4( new kinematics::MoveMap );

	// single backbone
	mm1->set_bb( 4, true );

	// all bb and chi
	mm2->set_bb( true );
	mm2->set_chi( true );

	// single dof
	mm3->set( DOF_ID( AtomID(1,4), id::PHI ), true );

	// everything!
	mm4->set( id::PHI, true );
	mm4->set( THETA, true );
	mm4->set( D, true );

	protocols::simple_moves::MinMover min_mover( mm1, scorefxn, "lbfgs_armijo_nonmonotone", 0.001, true /*use_nblist*/,
		true /*deriv_check*/ );

	dump_pdb( pose, "output/before.pdb" );

	//minimizer.run( pose, mm1, scorefxn, options );
	//dump_pdb( pose, "output/after1.pdb" );

	min_mover.movemap( mm2 );
	min_mover.apply( pose );
	dump_pdb( pose, "output/after2.pdb" );

	min_mover.movemap( mm3 );
	min_mover.apply( pose );
	dump_pdb( pose, "output/after3.pdb" );

	min_mover.movemap( mm4 );
	min_mover.apply( pose );
	dump_pdb( pose, "output/after4.pdb" );

}


///////////////////////////////////////////////////////////////////////////////
void
atom_tree_torsion_test()
{
	using id::AtomID;

	pose::Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

	pose.conformation().debug_residue_torsions();

	kinematics::AtomTree const & atom_tree( pose.atom_tree() );

	for ( Size i=2; i<= pose.total_residue() - 1; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		assert( Size( rsd.seqpos() ) == i );
		for ( Size chino=1; chino<= rsd.nchi(); ++chino ) {
			AtomID atom2( rsd.chi_atoms()[chino][2], i );
			AtomID atom3( rsd.chi_atoms()[chino][3], i );

			// loop over atoms bonded to atom2
			vector1< int > atom2_nbrs( rsd.bonded_neighbor( atom2.atomno() ) );
			vector1< int > atom3_nbrs( rsd.bonded_neighbor( atom3.atomno() ) );
			for ( Size ii=1; ii<= atom2_nbrs.size(); ++ii ) {
				for ( Size jj=1; jj<= atom3_nbrs.size(); ++jj ) {
					AtomID const atom1( atom2_nbrs[ii], i );
					AtomID const atom4( atom3_nbrs[jj], i );
					if ( atom1 == atom3 || atom2 == atom4 ) continue;
					Real offset;
					atom_tree.torsion_angle_dof_id( atom1, atom2, atom3, atom4, offset);
				}
			}
		}
	}


	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );

		for ( int r=1; r<= 2; ++r ) {
			id::TorsionType const type( r == 1 ? id::BB :
				id::CHI );
			Size const n( r == 1 ? rsd.mainchain_atoms().size() : rsd.nchi() );

			for ( Size j=1; j<= n; ++j ) {
				id::TorsionID const tor_id( i, type, j );
				//std::cout << "set_torsion: " << tor_id <<
				// " =============================================" << std::endl;
				pose.set_torsion( tor_id, 180.0 );
				pose.conformation().debug_residue_torsions();
				assert( std::abs( basic::subtract_degree_angles( pose.torsion( tor_id ),
					180.0 ) ) < 1e-3 );
			}
		}
	}
	dump_pdb( pose, "extended.pdb" );
}

///////////////////////////////////////////////////////////////////////////////
void
fa_scorefxn_test()
{
	// setup etable, fa scorefxn
	//using namespace core;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace pose;


	ScoreFunction scorefxn;
	scorefxn.set_weight( fa_atr, 0.80 );
	scorefxn.set_weight( fa_rep, 0.44 );
	scorefxn.set_weight( fa_sol, 0.65 );
	scorefxn.set_weight( fa_pair, 0.49 );
	scorefxn.set_weight( fa_dun, 0.49 );
	scorefxn.set_weight( rama, 0.2 );
	scorefxn.set_weight( hbond_lr_bb, 1.0 );
	scorefxn.set_weight( hbond_sr_bb, 1.0 );
	scorefxn.set_weight( hbond_bb_sc, 1.0 );
	scorefxn.set_weight( hbond_sc, 1.0 );

	Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

	std::cout << "pose-score: " << scorefxn( pose ) << std::endl;

	/// Now handled automatically.  scorefxn.accumulate_residue_total_energies( pose );

}

/// STOP REMOVING FUNCTIONS FROM TEST1
///////////////////////////////////////////////////////////////////////////////
void
rotamer_trials_test()
{
	//using namespace core;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace pose;

	ScoreFunction scorefxn;
	scorefxn.set_weight( fa_atr, 0.80 );
	scorefxn.set_weight( fa_rep, 0.44 );
	scorefxn.set_weight( fa_sol, 0.65 );
	scorefxn.set_weight( fa_pair, 0.49 );

	scorefxn.set_weight( hbond_sr_bb, 1.0 );
	scorefxn.set_weight( hbond_lr_bb, 1.0 );
	scorefxn.set_weight( hbond_bb_sc, 1.0 );
	scorefxn.set_weight( hbond_sc, 1.0 );

	std::cout << "atr-weight: " << scorefxn[fa_atr] << std::endl;
	std::cout << "rep-weight: " << scorefxn[fa_rep] << std::endl;


	Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);
	//  core::import_pose::pose_from_file( pose, "input/test_in_noPRO.pdb" , core::import_pose::PDB_file);
	Energy score_orig = scorefxn( pose );

	scorefxn.show( std::cout, pose );
	//for ( int jj = 1; jj <= pose.total_residue(); ++jj )
	//{
	pack::task::PackerTaskOP task
		( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().restrict_to_repacking();

	clock_t starttime = clock();
	pack::rotamer_trials( pose, scorefxn, task );

	clock_t stoptime = clock();
	std::cout << "TIMING rotamer_trials took: " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " seconds" << std::endl;

	Energy score = scorefxn( pose );
	std::cout << "Completed rotamer_trials_test() with new score: " << score << " vs orig: " << score_orig << std::endl;

	//pose.energies().show( std::cout );

	//std::string outname = "test_rotrials_";
	//std::stringstream ss;
	//ss << jj;
	//outname += ss.str();
	//outname += ".pdb";

	//std::cout << outname << std::endl;

	//dump_pdb( pose, outname.c_str() );
	dump_pdb( pose, "rotamer_trials_all.pdb");
	//}
}


///////////////////////////////////////////////////////////////////////////////
void
pack_rotamers_test()
{
	//using namespace core;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace pose;

	ScoreFunction scorefxn;
	scorefxn.set_weight( fa_atr, 0.80 );
	scorefxn.set_weight( fa_rep, 0.44 );
	scorefxn.set_weight( fa_sol, 0.65 );
	scorefxn.set_weight( fa_pair, 0.49 );

	scorefxn.set_weight( hbond_sr_bb, 1.0 );
	scorefxn.set_weight( hbond_lr_bb, 1.0 );
	scorefxn.set_weight( hbond_bb_sc, 1.0 );
	scorefxn.set_weight( hbond_sc, 1.0 );

	{ /* test 1 */
		Pose pose;
		core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);
		Energy score_orig = scorefxn( pose );

		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		task->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
		task->set_bump_check( true );

		clock_t starttime = clock();
		pack::pack_rotamers( pose, scorefxn, task);
		clock_t stoptime = clock();
		std::cout << "TIMING: pack_rotamers took " << ((double) stoptime - starttime)/CLOCKS_PER_SEC << std::endl;
		Energy score = scorefxn( pose );

		std::cout << "Completed pack_rotamers_test() #1 with new score: " << score << " vs orig: " << score_orig << std::endl;

		scorefxn.set_weight( ref, 1.0 );
		scorefxn.set_weight( p_aa_pp, 0.29 );
		score_orig = scorefxn( pose );

		dump_pdb( pose, "test_packrots.pdb" );
	} /* test1 */

	{ /* test2: design */
		Pose pose;
		core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);
		Energy score_orig = scorefxn( pose );
		pack::task::PackerTaskOP designtask( pack::task::TaskFactory::create_packer_task( pose ));
		designtask->initialize_from_command_line().or_include_current( true );
		designtask->set_bump_check( true );
		utility::vector1< bool > residues_to_repack( pose.total_residue(), false );
		for ( Size ii = 1; ii <= 20; ++ii ) residues_to_repack[ ii ] = true;
		designtask->restrict_to_residues( residues_to_repack );

		clock_t starttime = clock();
		pack::pack_rotamers( pose, scorefxn, designtask);
		clock_t stoptime = clock();
		std::cout << "TIMING pack_rotamers with design took " << ((double) stoptime - starttime)/CLOCKS_PER_SEC << std::endl;
		Energy design_score = scorefxn( pose );

		std::cout << "Completed pack_rotamers_test() #2 with new score: " << design_score << " vs orig: " << score_orig << std::endl;
		dump_pdb( pose, "test_design.pdb" );
	} /* test2 */

	{ /* test3: resfile */
		Pose pose;
		core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);
		Energy score_orig = scorefxn( pose );
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		task->initialize_from_command_line();
		core::pack::task::parse_resfile(pose, *task, "input/test_in.resfile");
		task->or_include_current( true );
		task->set_bump_check( true );

		clock_t starttime = clock();
		pack::pack_rotamers( pose, scorefxn, task);
		clock_t stoptime = clock();
		std::cout << "TIMING pack_rotamers after reading resfile took " << ((double) stoptime - starttime)/CLOCKS_PER_SEC << std::endl;
		Energy design_score = scorefxn( pose );

		std::cout << "Completed pack_rotamers_test() #3 with new score: " << design_score << " vs orig: " << score_orig << std::endl;
		dump_pdb( pose, "test_design.pdb" );
	} /* test3 */

}

///////////////////////////////////////////////////////////////////////////////
void
small_min_test()
{
	using namespace protocols::moves;
	using namespace scoring;

	using namespace basic::options;

	pose::Pose pose;
	core::import_pose::pose_from_file( pose, start_file() , core::import_pose::PDB_file);

	std::string const outfile_prefix( option[ OptionKeys::out::file::o ] );

	core::scoring::ScoreFunctionOP scorefxn( new ScoreFunction );

	// aiming for standard packer weights
	scorefxn->set_weight( fa_atr, 0.80 );
	scorefxn->set_weight( fa_rep, 0.44 );
	scorefxn->set_weight( fa_sol, 0.65 );
	scorefxn->set_weight( fa_pair, 0.49 );
	scorefxn->set_weight( fa_dun, 0.56 );
	scorefxn->set_weight( rama, 0.2 );
	scorefxn->set_weight( hbond_lr_bb, 1.17 );
	scorefxn->set_weight( hbond_sr_bb, 1.17 );
	scorefxn->set_weight( hbond_bb_sc, 1.17 );
	scorefxn->set_weight( hbond_sc   , 1.10 );

	// monte carlo object
	MonteCarloOP mc( new MonteCarlo( pose, *scorefxn, 0.8 /*temperature*/ ) );

	// the movable dof's
	kinematics::MoveMapOP mm( new kinematics::MoveMap );
	mm->set_bb ( true );
	mm->set_chi( true );

	// options for minimizer
	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.001, true /*use_nblist*/ ) );

	// packer options
	pack::task::PackerTaskOP task
		( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
	protocols::simple_moves::PackRotamersMoverOP pack_full_repack( new protocols::simple_moves::PackRotamersMover( scorefxn, task ) );
	/// @bug accumulate_residue_total_energies has to be called before anything can be done
	(*scorefxn)( pose );
	/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( pose );
	protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial( new protocols::simple_moves::EnergyCutRotamerTrialsMover( scorefxn, *task, mc, 0.01 /*energycut*/ ) );
	// pack_rottrial->setup_rottrial_task( pose, mc, 0.01 /*energycut*/ );

	// setup the move objects
	Size nmoves ( 5 );
	protocols::simple_moves::SmallMoverOP small_mover( new protocols::simple_moves::SmallMover( mm, 0.8/*temp*/, nmoves ) );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 2.0 );
	small_mover->angle_max( 'L', 3.0 );


	//dump_pdb( pose, "tmp_start.pdb" );
	{ // initial minimization
		TrialMoverOP min_trial( new TrialMover( min_mover, mc ) );
		min_trial->apply( pose );
	}

	pose::Pose start_pose;
	start_pose = pose;
	int const inner_cycle ( 30 );
	int const outer_cycle ( option[ OptionKeys::phil::nloop ] );

	basic::prof_reset();

	CycleMoverOP repack_cycle( new CycleMover );

	SequenceMoverOP main_min_seq( new SequenceMover );
	main_min_seq->add_mover( small_mover );
	main_min_seq->add_mover( pack_rottrial );
	main_min_seq->add_mover( min_mover );

	SequenceMoverOP pack_min_seq( new SequenceMover );
	pack_min_seq->add_mover( pack_full_repack );
	pack_min_seq->add_mover( pack_rottrial );
	pack_min_seq->add_mover( min_mover );

	TrialMoverOP main_min_trial( new TrialMover( main_min_seq, mc ) );
	TrialMoverOP pack_min_trail( new TrialMover( pack_min_seq, mc ) );
	RepeatMoverOP main_min_cycle( new RepeatMover( main_min_trial, inner_cycle ) );
	ProfilerMoverOP profiler( new ProfilerMover );

	SequenceMoverOP full_seq( new SequenceMover );
	full_seq->add_mover( main_min_cycle );
	full_seq->add_mover( pack_min_seq );
	full_seq->add_mover( profiler );

	RepeatMoverOP full_min_cycle( new RepeatMover( full_seq, outer_cycle ) );

	for ( int n=1; n<= option[ OptionKeys::out::nstruct ]; ++n ) {
		pose::Pose relax_pose;
		relax_pose = start_pose;

		full_min_cycle->apply( relax_pose );
		basic::prof_show();

		if ( outfile_prefix != "none" ) dump_pdb( relax_pose, outfile_prefix+string_of( n )+".pdb" );

	}

}

//-----------------------------------------------------------------------------------------A
//_________________________________________________________________________________________A
//
//    |
//    |
//    |
//    |
//
///////////////////////////////////////////////////////////////////////////////
void
rb_test()
{
	using namespace pose;
	using namespace moves;;
	using namespace scoring;
	using namespace kinematics;
	using namespace conformation;
	using namespace chemical;

	pose::Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

	int const nres( pose.total_residue() );
	assert( nres > 40 );

	{ // setup a foldtree
		FoldTree f( nres );
		f.new_jump( 8, 40, 15 );
		f.reorder( 8 );
		pose.fold_tree( f );
	}

	dump_pdb( pose, "tmp1.pdb" );

	//  kinematics::Stub stub1( pose.conformation().upstream_jump_stub(1) ),
	//   stub2( pose.conformation().downstream_jump_stub(1) );
	//  exit(0);

	{ // now add a pseudo residue at the end
		ResidueTypeSetCOP residue_set
			( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		ResidueOP rsd( ResidueFactory::create_residue( residue_set->name_map( "VRT" ) ) );
		pose.append_residue_by_jump( *rsd, 1 );
		dump_pdb( pose, "tmp2.pdb" );

		FoldTree f( pose.fold_tree() );
		f.reorder( f.nres() );
		pose.fold_tree( f );

		dump_pdb( pose, "tmp3.pdb" );

		pose.set_jump( 2, Jump() );

		dump_pdb( pose, "tmp4.pdb" );
	}

	// test scoring
	ScoreFunction scorefxn;

	// aiming for standard packer weights
	scorefxn.set_weight( fa_atr, 0.80 );
	scorefxn.set_weight( fa_rep, 0.44 );
	scorefxn.set_weight( fa_sol, 0.65 );
	scorefxn.set_weight( fa_pair, 0.49 );
	scorefxn.set_weight( fa_dun, 0.56 );
	scorefxn.set_weight( rama, 0.2 );
	if ( true ) {
		scorefxn.set_weight( hbond_lr_bb, 1.17 );
		scorefxn.set_weight( hbond_sr_bb, 1.17 );
		scorefxn.set_weight( hbond_bb_sc, 1.17 );
		scorefxn.set_weight( hbond_sc   , 1.10 );
	}

	scorefxn( pose );
	// pose.energies().show( std::cout );

	// now test the simple rigid-body move stuff
	MoveMap mm;
	mm.set_jump( 1, true );

	Pose start_pose;
	start_pose = pose;

	FoldTree f1( pose.fold_tree() ), f2( pose.fold_tree() );
	f1.reorder( 8 );
	f2.reorder( 40 );

	start_pose.fold_tree( f1 );

	// forward rotation
	rigid::RigidBodyPerturbMoverOP rb_mover( new rigid::RigidBodyPerturbMover(
		1 /*jump_num*/, 5.0 /*rot*/, 0.0 /*trans*/ ) );
	moves::PDBDumpMoverOP dumper( new PDBDumpMover( "tmp_fwd_rotation_" ) );
	moves::SequenceMoverOP sequencer( new SequenceMover );
	sequencer->add_mover( rb_mover );
	sequencer->add_mover( dumper );

	moves::RepeatMoverOP cycler( new RepeatMover( sequencer, 10 ) );
	pose = start_pose;
	cycler->apply( pose );

	// forward translation
	pose = start_pose;
	dumper->name( "tmp_fwd_translation_" );
	rb_mover->trans_magnitude( 1.0 );
	rb_mover->rot_magnitude( 0.0 );
	cycler->apply( pose );

	start_pose.fold_tree( f2 );

	// reverse rotation
	pose = start_pose;
	dumper->name( "tmp_rev_rotation_" );
	rb_mover->trans_magnitude( 0.0 );
	rb_mover->rot_magnitude( 5.0 );
	cycler->apply( pose );

	// reverse translation
	pose = start_pose;
	dumper->name( "tmp_rev_tranlation_" );
	rb_mover->trans_magnitude( 1.0 );
	rb_mover->rot_magnitude( 0.0 );
	cycler->apply( pose );
}


///////////////////////////////////////////////////////////////////////////////
void
ccd_test()
{
	using namespace pose;
	using namespace scoring;
	using namespace conformation;
	using namespace chemical;
	using namespace optimization;
	using namespace id;

	Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

	int const nres( pose.total_residue() );
	assert( nres > 40 );

	int const cutpoint( 18 );
	{ // setup a foldtree
		kinematics::FoldTree f( nres );
		f.new_jump( 8, 26, cutpoint );
		f.reorder( 8 );
		pose.fold_tree( f );

		// introduce new variant types
		if ( true ) { //false ) {
			core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, cutpoint   );
			core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, cutpoint+1 );
		}
	}

	pose.set_phi( 16, pose.phi(16) + 15.0 );
	pose.set_psi( 16, pose.psi(16) - 15.0 );

	dump_pdb( pose, "tmp_before.pdb" );

	Pose start_pose;
	start_pose = pose;

	kinematics::MoveMapOP mm( new kinematics::MoveMap );

	// setup moving dofs
	for ( int i=15; i<= 22; ++i ) {
		mm->set_bb ( i, true );
		mm->set( TorsionID( i, BB, 3 ), false ); // omega off
	}

	{ // ccd closure

		protocols::loops::Loop loop( 15, 22, cutpoint );
		protocols::loops::loop_closure::ccd::CCDLoopClosureMover ccd_loop_closure_mover( loop, mm );
		ccd_loop_closure_mover.tolerance( 0.05 );
		ccd_loop_closure_mover.rama()->max_rama_score_increase( 100 );
		ccd_loop_closure_mover.max_total_torsion_delta_per_residue( 100, 100, 100 );
		ccd_loop_closure_mover.apply( pose );

		dump_pdb( pose, "tmp_after_ccd.pdb" );

	}

	{ // minimization

		// setup the options
		MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.001, true /*use_nblist*/, true /*deriv_check*/ );

		AtomTreeMinimizer minimizer;
		scoring::ScoreFunction scorefxn;
		scorefxn.set_weight( scoring::chainbreak, 1.0 );

		pose = start_pose;
		minimizer.run( pose, *mm, scorefxn, options );

		dump_pdb( pose, "tmp_after_dfpmin.pdb" );
	}


}


///////////////////////////////////////////////////////////////////////////////
void
sasa_test()
{
	using namespace pose;
	using namespace scoring;
	using namespace conformation;

	Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_in_w_termini.pdb" , core::import_pose::PDB_file);

	id::AtomID_Map< Real > atom_sasa;
	utility::vector1< Real > rsd_sasa;

	Real const total_sasa( scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, 1.4 ) );

	std::cout << "total_sasa: " << total_sasa << std::endl;
	dump_pdb( pose, "dump_test_in_w_termini.pdb" );
}

///////////////////////////////////////////////////////////////////////////////
// should run to completion
//
void
simple_benchmark()
{
	// this should improve stability of the tests
	core::init::init_random_generators( 1000, "mt19937" );

	//ccd_test(  );
	//exit(0);

	// Simple_min_test(  );
	// Fri Oct 05 13:11:51 EDT 2007 @758 /Internet Time/
	/// Simple_min_test was moved to test/core/optimization/Minimizer.cxxtest.hh

	simple_rotamer_test(  );

	//for (Size ii = 1; ii <= 100; ++ii )
	rotamer_trials_test(  );

	// Wed Oct 10 10:03:07 EDT 2007 @627 /Internet Time/
	// rotamer_trials_test was moved to test/core/pack/RotamerTrials.cxxtest.hh
	// and put back in Oct 20 at 10:09 AM PDS

	pack_rotamers_test(  );

	fa_scorefxn_test(  );

	test_rama(  );

	// simple_conformation_test(  );
	// Thu Nov 08 08:23:43 EST 2007 @599 /Internet Time/
	// simple_conformation_test(  ); was moved to test/core/conformation/Conformation.cxxtest.hh

	atom_tree_torsion_test(  );

	ccd_test(  );

	sasa_test(  );

	simple_dna_test(  );

	dna_coupled_rotamer_design_test();

	dna_deriv_test();

	test_scorefxn_io();

	delete_test();
}


///////////////////////////////////////////////////////////////////////////////
void
mm_pack_test()
{
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace pose;

	// scrfxns
	ScoreFunctionOP score_12( get_score_function() );
	ScoreFunctionOP score_mod( get_score_function() );
	score_mod->set_weight( fa_dun, 0.00 );
	ScoreFunctionOP score_mm_only( new ScoreFunction ); score_mm_only->set_weight( mm_twist,   1.00 );

	// weighting
	Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

	Energy ener_12  = (*score_12)(pose);
	Energy ener_mod = (*score_mod)(pose);
	Energy ener_mm  = (*score_mm_only)(pose);

	float mm_weight = (ener_12 - ener_mod)/ener_mm;

	ScoreFunctionOP score_mm( get_score_function() );
	score_mm->set_weight( fa_dun, 0.00 ); score_mm->set_weight( mm_twist, mm_weight );

	// pack
	Pose pose_12;
	core::import_pose::pose_from_file( pose_12, "input/test_in.pdb" , core::import_pose::PDB_file);

	Pose pose_mm;
	core::import_pose::pose_from_file( pose_mm, "input/test_in.pdb" , core::import_pose::PDB_file);

	//  Pose pose_mod;
	//  core::import_pose::pose_from_file( pose_mod, "input/test_in.pdb" , core::import_pose::PDB_file);

	pack::task::PackerTaskOP task_12( pack::task::TaskFactory::create_packer_task( pose_12 ));
	task_12->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
	task_12->set_bump_check( true );

	pack::task::PackerTaskOP task_mm( pack::task::TaskFactory::create_packer_task( pose_mm ));
	task_mm->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
	task_mm->set_bump_check( true );

	//  pack::task::PackerTaskOP task_mod( pack::task::TaskFactory::create_packer_task( pose_mod ));
	//  task_mm->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
	// task_mm->set_bump_check( true );


	Energy orig_score_12 = (*score_12)(pose_12);
	clock_t starttime_12 = clock(); pack::pack_rotamers( pose_12, *score_12, task_12); clock_t stoptime_12 = clock();
	Energy final_score_12 = (*score_12)(pose_12);

	Energy orig_score_mm = (*score_mm)(pose_mm);
	clock_t starttime_mm = clock(); pack::pack_rotamers( pose_mm, *score_mm, task_mm); clock_t stoptime_mm = clock();
	Energy final_score_mm = (*score_mm)(pose_mm);

	//  Energy orig_score_mod = (*score_mod)(pose_mod);
	//  clock_t starttime_mod = clock(); pack::pack_rotamers( pose_mod, *score_mod, task_mod); clock_t stoptime_mod = clock();
	//  Energy final_score_mod = (*score_mod)(pose_mod);

	//output
	std::cout << "mm_pack_test score 12 orig: " << orig_score_12 << " final score: " << final_score_12 << "in "
		<<  ((double) stoptime_12 - starttime_12)/CLOCKS_PER_SEC << " seconds" << std::endl;
	std::cout << "mm_pack_test score mm orig: " << orig_score_mm << " final score: " << final_score_mm << "in "
		<<  ((double) stoptime_mm - starttime_mm)/CLOCKS_PER_SEC << " seconds" << std::endl;
	//  std::cout << "mm_pack_test score mod orig: " << orig_score_mod << " final score: " << final_score_mod << "in "
	//       <<  ((double) stoptime_mod - starttime_mod)/CLOCKS_PER_SEC << " seconds" << std::endl;
	std::cout << "MM Weight is: " << mm_weight << std::endl;

	dump_pdb(pose_12, "test_out_12.pdb");
	dump_pdb(pose_mm, "test_out_mm.pdb");
	// dump_pdb(pose_mod, "test_out_mod.pdb");

}

///////////////////////////////////////////////////////////////////////////////
void
start_file_test()
{
	utility::vector1< std::string > files( basic::options::start_files() );
	for ( Size i=1; i<= files.size(); ++i ) {
		std::cout << i << ' ' << files[i] << std::endl;
	}
	std::cout << "start_file(): " << basic::options::start_file() << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
// Phil's test function, a little expanded by Rhiju for centroid readin and
// testing multiple pdbs at once.
void
ss_test()
{

	using namespace basic::options;
	using namespace chemical;

	using namespace scoring;


	ScoreFunction sfxn;
	sfxn.set_weight( hs_pair, 1.0 );
	sfxn.set_weight( ss_pair, 1.0 );
	sfxn.set_weight( rsigma, 1.0 );

	/// Let's get a few PDBs to test
	utility::vector1< std::string > pdbnames( basic::options::start_files() );

	// Centroid! Affects the neighbor list.
	ResidueTypeSetCOP rsd_set( ChemicalManager::get_instance()->residue_type_set( "centroid" ) );

	for ( Size ii = 1; ii < pdbnames.size(); ++ii ) {
		pose::Pose pose;
		std::string const pdbname = pdbnames[ii];

		core::import_pose::pose_from_file( pose, *rsd_set, pdbname , core::import_pose::PDB_file);

		// protocols::loops::set_secstruct_from_psipred_ss2( pose ); // uses -in::file::psipred_ss2 <ss-filename>

		std::string dsspname = pdbname;
		dsspname.replace( dsspname.rfind(".pdb", dsspname.length() ), 4, ".dssp" );

		protocols::loops::set_secstruct_from_dssp( pose, dsspname );

		std::cout << sfxn( pose ) << std::endl;

		std::cout << pdbname;
		pose.energies().total_energies().show_if_nonzero_weight( std::cout, sfxn.weights() );
		std::cout << std::endl;
	}

}


///////////////////////////////////////////////////////////////////////////////
void
backrub_min_test()
{
	using namespace id;
	using namespace chemical;
	using namespace conformation;
	using namespace kinematics;

	pose::Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_short.pdb" , core::import_pose::PDB_file);

	utility::vector1< AtomID > mainchain;
	Size const seqpos( 3 );

	mainchain.push_back( AtomID( pose.residue( seqpos-1 ).atom_index("CA"), seqpos-1 ) );
	mainchain.push_back( AtomID( pose.residue( seqpos-1 ).atom_index( "C"), seqpos-1 ) );
	mainchain.push_back( AtomID( pose.residue( seqpos   ).atom_index( "N"), seqpos   ) );
	mainchain.push_back( AtomID( pose.residue( seqpos   ).atom_index("CA"), seqpos   ) );
	mainchain.push_back( AtomID( pose.residue( seqpos   ).atom_index( "C"), seqpos   ) );
	mainchain.push_back( AtomID( pose.residue( seqpos+1 ).atom_index( "N"), seqpos+1 ) );
	mainchain.push_back( AtomID( pose.residue( seqpos+1 ).atom_index("CA"), seqpos+1 ) );

	//AtomID const downstream_id( AtomID( pose.residue( seqpos+1 ).atom_index("C"), seqpos+1 ) );  // unused ~Labonte

	utility::vector1< std::pair< Size, Size > > edges;
	edges.push_back( std::make_pair( 1, 7 ) ); // CA of seqpos-1 --> CA of seqpos+1
	edges.push_back( std::make_pair( 7, 4 ) ); //       seqpos+1 -->       seqpos
	edges.push_back( std::make_pair( 4, 1 ) ); //       seqpos   -->       seqpos-1

	ResidueTypeSet const & rsd_set( *pose.residue(1).residue_type_set() );
	ResidueOP psd( ResidueFactory::create_residue( rsd_set.name_map( "VRT1" ) ) );
	pose.append_residue_by_jump( *psd, 1 );
	pose.append_residue_by_jump( *psd, 1 );
	pose.append_residue_by_jump( *psd, 1 );

	Size const first_new_pseudo_residue( pose.total_residue() - 2 );

	utility_exit_with_message("need to refactor backrub atomtree stuff");
	//pose.conformation().setup_backrub_segment( mainchain, downstream_id, edges, first_new_pseudo_residue );

	protocols::viewer::dump_pose_kinemage( "tmp2.kin", pose );


	pose::Pose start_pose;
	start_pose = pose;
	for ( Size i=1; i<= edges.size(); ++i ) {
		pose = start_pose;
		AtomID const id( 1, first_new_pseudo_residue + i - 1 );
		kinematics::tree::Atom const & atom( pose.atom_tree().atom( id ) );
		assert( atom.n_children() > 0 );

		DOF_ID const dof( atom.child(0)->id(), id::PHI );
		std::cout << "pseudo 1st child: " << i << ' ' <<dof << std::endl;
		Real const orig( pose.dof( dof ) );
		for ( int k=-5; k<=5; ++k ) {
			pose.set_dof( dof, orig+k*0.05 );
			pose.dump_pdb("backrub_move_"+string_of(i)+"_"+string_of(k)+".pdb" );
		}
	}

	pose = start_pose;

	{ // try minimizing those dofs
		using namespace optimization;
		using namespace scoring;

		MoveMap mm;
		for ( Size i=1; i<= edges.size(); ++i ) {
			AtomID const id( 1, first_new_pseudo_residue + i - 1 );
			kinematics::tree::Atom const & atom( pose.atom_tree().atom( id ) );
			assert( atom.n_children() > 0 );

			DOF_ID const dof( atom.child(0)->id(), id::PHI );
			mm.set( dof, true );
		}

		ScoreFunction scorefxn;
		scorefxn.set_weight( fa_atr, 0.80 );
		scorefxn.set_weight( fa_rep, 0.44 );
		scorefxn.set_weight( fa_sol, 0.65 );

		scorefxn.set_weight( hbond_lr_bb, 1.0 );
		scorefxn.set_weight( hbond_sr_bb, 1.0 );
		scorefxn.set_weight( hbond_bb_sc, 1.0 );
		scorefxn.set_weight( hbond_sc, 1.0 );

		pose.dump_pdb( "before.pdb");
		AtomTreeMinimizer().run( pose, mm, scorefxn, MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.001, true ) );
		pose.dump_pdb( "after.pdb");

	}

	//  Size const first_new_pseudo_residue( pose.total_residue() + 1 );

	//  kinematics::tree::Atom * root( pose.atom_tree().root()->clone(0) );
	//  AtomPointers old_atom_pointer;
	//  root->update_atom_pointer( old_atom_pointer );

	//  for ( Size i=1; i<= mainchain.size(); ++i ) {
	//   assert( mainchain[i] == old_atom_pointer[mainchain[i]]->id() );
	//   std::cout << mainchain[i] << std::endl;
	//  }

	//  kinematics::tree::Atom * new_root
	//   ( setup_backrub_atom_tree( mainchain, downstream_id, old_atom_pointer, edges, first_new_pseudo_residue ) );

	//  new_root->show();

	//  AtomTree at( new_root );

	//  dump_atomtree_kinemage( "tmp.kin", at, pose.conformation() );
}

///////////////////////////////////////////////////////////////////////////////
void
bk_test()
{
	pose::Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

	Size const nres( pose.total_residue() );
	pose.dump_pdb("test1.pdb");
	core::pose::remove_upper_terminus_type_from_pose_residue( pose, nres );
	pose.dump_pdb("test2.pdb");

	pose::Pose pose2;
	core::import_pose::pose_from_file( pose2, "input/test_in.pdb" , core::import_pose::PDB_file);
	core::pose::remove_lower_terminus_type_from_pose_residue( pose2, 1 );

	for ( Size i=1; i<= nres; ++i ) {
		pose.append_residue_by_bond( pose2.residue(i) );
	}
	pose.dump_pdb("test3.pdb");

	pose.conformation().insert_ideal_geometry_at_polymer_bond( nres );

	pose.dump_pdb("test4.pdb");

	pose.set_phi  ( nres  , -150.0 );
	pose.set_psi  ( nres  ,  150.0 );
	pose.set_omega( nres  ,  180.0 );
	pose.set_phi  ( nres+1, -150.0 );
	pose.set_psi  ( nres+1,  150.0 );
	pose.set_omega( nres+1,  180.0 );

	pose.dump_pdb("test5.pdb");

	exit(0);
}

///////////////////////////////////////////////////////////////////////////////
void
bk_test2()
{
	using numeric::conversions::radians;

	pose::Pose pose1, pose2;
	core::import_pose::pose_from_file( pose1, start_files()[1] , core::import_pose::PDB_file);
	core::import_pose::pose_from_file( pose2, start_files()[2] , core::import_pose::PDB_file);

	ResidueTypeSet const & rsd_set( *pose1.residue(1).residue_type_set() );
	Size const nres1( pose1.total_residue() );
	Size const nres2( pose2.total_residue() );

	Size cyx_pos(0);
	Size lyx_pos(0);
	{ // find a cysteine and lysine
		for ( Size i=2; i<= nres1-1; ++i ) {
			if ( !cyx_pos && pose1.residue(i).aa() == aa_cys ) {
				cyx_pos = i;
			}
			if ( !lyx_pos && pose1.residue(i).aa() == aa_lys ) {
				lyx_pos = i;
			}
		}
	}

	std::cout << "nres1= " << nres1 << " cyx_pos= " << cyx_pos << " lyx_pos= " << lyx_pos <<
		" nres2= " << nres2 << std::endl;

	ResidueType const & cyx_rsd_type( rsd_set.name_map("CYX") ); // have to be a tiny bit more careful if the cys has variants
	ResidueType const & lyx_rsd_type( rsd_set.name_map("LYX") ); // have to be a tiny bit more careful if the lys has variants
	pose1.dump_pdb("test1.pdb");
	core::pose::replace_pose_residue_copying_existing_coordinates( pose1, cyx_pos, cyx_rsd_type );
	core::pose::replace_pose_residue_copying_existing_coordinates( pose1, lyx_pos, lyx_rsd_type );
	pose1.dump_pdb("test2.pdb");


	core::pose::remove_variant_type_from_pose_residue( pose2, UPPER_TERMINUS_VARIANT, nres2 );
	core::pose::add_variant_type_to_pose_residue( pose2, CTERM_CONNECT, nres2 );


	// determine the desired connection points
	//
	//
	ResidueType const & ubq_rsd_type( pose2.residue_type( nres2 ) );
	Size const cyx_connid( 3 );
	Size const lyx_connid( 3 );
	Size const ubq_connid( 2 );

	assert( cyx_rsd_type.n_possible_residue_connections() == cyx_connid &&
		cyx_rsd_type.lower_connect_id() != cyx_connid &&
		cyx_rsd_type.upper_connect_id() != cyx_connid );

	assert( lyx_rsd_type.n_possible_residue_connections() == lyx_connid &&
		lyx_rsd_type.lower_connect_id() != lyx_connid &&
		lyx_rsd_type.upper_connect_id() != lyx_connid );

	assert( ubq_rsd_type.n_possible_residue_connections() == ubq_connid &&
		ubq_rsd_type.lower_connect_id() != ubq_connid &&
		ubq_rsd_type.is_upper_terminus() );


	tt << "before append_by_bond\n";
	pose1.conformation().show_residue_connections();

	{ /// now attach pose2 to the cysteine:
		pose1.append_residue_by_bond( pose2.residue( nres2 ), false, ubq_connid, cyx_pos, cyx_connid );

		Size const prepend_seqpos( pose1.total_residue() );

		tt << "after append_by_bond\n";
		pose1.conformation().show_residue_connections();

		for ( Size i=nres2-1; i>= 1; --i ) {
			pose1.prepend_polymer_residue_before_seqpos( pose2.residue(i), prepend_seqpos, false );
		}
		tt << "after all prepends: foldtree: " << pose1.fold_tree() << std::endl;
		pose1.conformation().show_residue_connections();
		pose1.dump_pdb("test3.pdb");

		pose1.conformation().insert_ideal_geometry_at_residue_connection( cyx_pos, cyx_connid );

		pose1.dump_pdb("test4.pdb");

		// set some torsions...
		// ca - cb - sg --- c - ca - n
		//
		Size const ubq_pos( pose1.total_residue() );
		AtomID const atom1( cyx_rsd_type.atom_index( "CA" ), cyx_pos );
		AtomID const atom2( cyx_rsd_type.atom_index( "CB" ), cyx_pos );
		AtomID const atom3( cyx_rsd_type.atom_index( "SG" ), cyx_pos );
		AtomID const atom4( ubq_rsd_type.atom_index( "C"  ), ubq_pos );
		AtomID const atom5( ubq_rsd_type.atom_index( "CA" ), ubq_pos );
		AtomID const atom6( ubq_rsd_type.atom_index( "N"  ), ubq_pos );

		// should convert this interface to use degrees...
		pose1.conformation().set_torsion_angle( atom1, atom2, atom3, atom4, radians(60.0) );
		pose1.conformation().set_torsion_angle( atom2, atom3, atom4, atom5, radians(60.0) );
		pose1.conformation().set_torsion_angle( atom3, atom4, atom5, atom6, radians(60.0) );

		pose1.dump_pdb("test5.pdb");

		// now try minimizing
		kinematics::MoveMap mm;
		mm.set_chi( cyx_pos, true );
		mm.set( pose1.atom_tree().torsion_angle_dof_id( atom1, atom2, atom3, atom4 ), true );
		mm.set( pose1.atom_tree().torsion_angle_dof_id( atom2, atom3, atom4, atom5 ), true );
		mm.set( pose1.atom_tree().torsion_angle_dof_id( atom3, atom4, atom5, atom6 ), true );
		mm.set_bb ( ubq_pos, true );

		ScoreFunction scorefxn;
		scorefxn.set_weight( fa_atr, 1.0 );
		scorefxn.set_weight( fa_rep, 1.0 );

		// this doen't work very well in my simple case since the cysteine is somewhat internal and the
		// clashes are so bad that it moves too far and violates the nblist
		//
		AtomTreeMinimizer().run( pose1, mm, scorefxn, MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.0001, true ) );

		pose1.dump_pdb("test6.pdb");
	}

	{ /// now attach pose2 to the lysine:
		pose1.append_residue_by_bond( pose2.residue( nres2 ), false, ubq_connid, lyx_pos, lyx_connid );

		tt << "after append_by_bond\n";
		pose1.conformation().show_residue_connections();

		Size const prepend_seqpos( pose1.total_residue() );

		for ( Size i=nres2-1; i>= 1; --i ) {
			pose1.prepend_polymer_residue_before_seqpos( pose2.residue(i), prepend_seqpos, false );
		}
		tt << "after all prepends: foldtree: " << pose1.fold_tree() << std::endl;
		pose1.conformation().show_residue_connections();
		pose1.dump_pdb("test7.pdb");

		pose1.conformation().insert_ideal_geometry_at_residue_connection( lyx_pos, lyx_connid );

		pose1.dump_pdb("test8.pdb");

		// set some torsions...


		pose1.dump_pdb("test9.pdb");
	}


	exit(0);
}

///////////////////////////////////////////////////////////////////////////////
void
ligrot_test()
{
	utility::vector1< std::string > files( basic::options::start_files() );

	std::ofstream out( "AKG_rotamers.pdb" );
	for ( Size n=1; n<= files.size(); ++n ) {
		pose::Pose pose;
		core::import_pose::pose_from_file( pose, files[n] , core::import_pose::PDB_file);

		conformation::Residue const & rsd( pose.residue(1) );
		for ( Size i=1; i<= rsd.nchi(); ++i ) {
			std::cout << "CHI " << I(4,i) << F(9,3,rsd.chi(i) ) << ' ' << files[i] << '\n';
		}

		pose.dump_pdb( out );
		//out << "TER\n";
	}
	out.close();
}


///////////////////////////////////////////////////////////////////////////////
void
proclose_test()
{
	using namespace optimization;

	// standard packer wts
	ScoreFunctionOP scorefxn( get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS ) );

	Pose pose;
	core::import_pose::pose_from_file( pose, start_file() , core::import_pose::PDB_file);

	Size const nres( pose.total_residue() );

	(*scorefxn)(pose);
	scorefxn->show( std::cout, pose );


	pose.dump_pdb( "test1.pdb" );

	// now try minimizing chi angles
	kinematics::MoveMap mm;
	for ( Size i=1; i<= nres; ++i ) {
		if ( pose.residue(i).aa() == aa_pro ) {
			mm.set_chi( i, true );
		}
	}
	AtomTreeMinimizer().run( pose, mm, *scorefxn, MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.001, true ) );
	(*scorefxn)(pose);
	scorefxn->show( std::cout, pose );


	pose.dump_pdb( "test2.pdb" );


}


/////////////////////////////////////////////////////////////////////////////////
//void
//hbond_plot_test()
//{
// scoring::hbonds::show_poly();
//}


///////////////////////////////////////////////////////////////////////////////
void
set_stub_transform_test()
{
	using kinematics::RT;
	using kinematics::Stub;

	pose::Pose pose;
	core::import_pose::pose_from_file( pose, "input/test_in.pdb" , core::import_pose::PDB_file);

	// set a new foldtree
	Size const nres( pose.total_residue() );
	kinematics::FoldTree f( nres );

	Size const pos1( 5 ), pos2( 45 ), cutpoint( 35 );
	f.new_jump( pos1, pos2, cutpoint );
	f.reorder( pos1 );
	pose.fold_tree( f );

	conformation::Residue const & rsd1( pose.residue( pos1 ) );
	conformation::Residue const & rsd2( pose.residue( pos2 ) );
	AtomID const a1( rsd1.atom_index("N" ), pos1 );
	AtomID const a2( rsd1.atom_index("CA"), pos1 );
	AtomID const a3( rsd1.atom_index("C" ), pos1 );
	AtomID const b1( rsd2.atom_index("N" ), pos2 );
	AtomID const b2( rsd2.atom_index("CA"), pos2 );
	AtomID const b3( rsd2.atom_index("C" ), pos2 );
	StubID const stub_id1( a1, a2, a3 );
	StubID const stub_id2( b1, b2, b3 );

	RT rt1( RT::Matrix::identity(), Vector( 50.0, 0.0, 0.0 ) );

	pose.conformation().set_stub_transform( stub_id1, stub_id2, rt1 );
	pose.dump_pdb( "test1.pdb" );

	// now transform so that pos2 is moved on top of pos3 (pos1,pos3 fixed, pos2 moving)
	Size const pos3( 25 );
	conformation::Residue const & rsd3( pose.residue( pos3 ) );
	AtomID const c1( rsd3.atom_index("N" ), pos3 );
	AtomID const c2( rsd3.atom_index("CA"), pos3 );
	AtomID const c3( rsd3.atom_index("C" ), pos3 );

	RT const rt2( Stub( pose.xyz(a1), pose.xyz(a2), pose.xyz(a3) ),
		Stub( pose.xyz(c1), pose.xyz(c2), pose.xyz(c3) ) );

	pose.conformation().set_stub_transform( stub_id1, stub_id2, rt2 );
	pose.dump_pdb( "test2.pdb" );

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
lk_ball_wtd_deriv_test()
{
	using namespace optimization;
	using std::cout;
	using std::endl;

	Pose pose;
	string const filename( start_file() );
	pose_from_file( pose, filename , core::import_pose::PDB_file);


	kinematics::MoveMap mm;
	mm.set_bb( false );
	mm.set_chi( true );
	mm.set_jump( false );

	/// scorefxn
	ScoreFunctionOP scorefxn( new ScoreFunction() );
	scorefxn->set_weight( fa_atr, 0.01 );
	scorefxn->set_weight( fa_rep, 0.01 );
	scorefxn->set_weight( lk_ball_wtd, 0.8 );

	Real const scorebefore( (*scorefxn)( pose ) );
	cout << "beforescore: " << F(9,3,scorebefore) << ' ' <<
		pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & irsd( pose.residue(i) );
		//scoring::methods::EnvElecResidueInfo const & irsd_info( retrieve_residue_info( pose, i ) );

		cout << "RSDE " << I(4,i) << ' ' << irsd.name1() <<
			" tot: " << F(9,3,pose.energies().residue_total_energies(i).dot( scorefxn->weights() ) ) <<
			' ' << pose.energies().residue_total_energies(i).weighted_string_of( scorefxn->weights() ) <<
			' ' << irsd.name() << ' ' << filename << endl;
	}


	pose.dump_pdb( "before.pdb" );


	{ // hacking
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ) );
		task->initialize_from_command_line();
		task->or_include_current( true );
		task->restrict_to_repacking();
		utility::vector1< bool > subset( pose.total_residue(), false );
		subset[10] = subset[11] = subset[12] = subset[13] = subset[14] = true;
		task->restrict_to_residues( subset );
		pack::RTMin rtmin;
		rtmin.rtmin( pose, *scorefxn, task );
		Real scoreafter( (*scorefxn)( pose ) );
		cout << "afterscore1: " << F(9,3,scoreafter) << ' ' <<
			pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;
		pack::min_pack( pose, *scorefxn, task );
		scoreafter = (*scorefxn)( pose );
		cout << "afterscore2: " << F(9,3,scoreafter) << ' ' <<
			pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

		//exit(0);
	}


	MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.1, true /*use_nblist*/,
		true /*deriv_check*/, false /*no verbose-deriv-check, is default*/ );
	AtomTreeMinimizer minimizer;

	minimizer.run( pose, mm, *scorefxn, options );

	pose.dump_pdb( "after.pdb" );

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
atom_types_test()
{
	using namespace std;
	using namespace core;
	using namespace chemical;


	ResidueTypeSet const & rsd_set( *ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	AtomTypeSet const & atom_set( *ChemicalManager::get_instance()->atom_type_set( FA_STANDARD ) );

	for ( Size i=1; i<= atom_set.n_atomtypes(); ++i ) {
		cout << "atom_set: " << I(4,i) << ' ' << atom_set[i].name() <<
			F(9,3,atom_set[i].lj_radius()) <<
			F(9,3,atom_set[i].lk_dgfree()) << std::endl;
	}

	{
		ResidueType const & rsd( rsd_set.name_map("ARG") );
		for ( Size i=1; i<= rsd.natoms(); ++i ) {
			cout << rsd.name() << I(3,i) << ' ' << rsd.atom_name(i) << ' ' <<
				rsd.atom_type(i).name() << endl;
		}
	}

	{
		ResidueType const & rsd( rsd_set.name_map("ADE") );
		for ( Size i=1; i<= rsd.natoms(); ++i ) {
			cout << rsd.name() << I(3,i) << ' ' << rsd.atom_name(i) << ' ' <<
				rsd.atom_type(i).name() << endl;
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
		using namespace basic::options;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// setup
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		devel::init(argc, argv);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// end of setup
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		/////////////////////////////////////
		// If you want to add a temporary
		// "run this and exit" behavior to
		// main, add it AFTER the following check
		if ( option[ OptionKeys::run::benchmark ] ) {
			simple_benchmark(  );
			exit(0);
		}
		/////////////////////////////////////

		// atom_types_test();
		// exit(0);

		lk_ball_wtd_deriv_test();
		exit(0);

		dna_deriv_test();
		exit(0);

		delete_test();
		exit(0);

		//hbond_plot_test();
		//exit(0);


		set_stub_transform_test();
		exit(0);

		bk_test2();
		exit(0);

		proclose_test();
		exit(0);

		ligrot_test();
		exit(0);

		bk_test();
		exit(0);

		backrub_min_test();
		exit(0);

		start_file_test();
		exit(0);

		ss_test();
		exit(0);

		pack_rotamers_test();
		exit( 0 );

		rb_test();
		exit(0);

		test_gb();
		exit(0);

		//dna_io_test();
		//exit(0);

		//list_dihedrals();
		//mm_library_test();
		//mm_score_test();
		mm_pack_test();
		exit(0);

		// rotamer_trials_test(  );
		// Wed Oct 10 10:03:07 EDT 2007 @627 /Internet Time/
		// rotamer_trials_test was moved to test/core/pack/RotamerTrials.cxxtest.hh
		pack_rotamers_test(  );
		exit(0);

		test_scorefxn_io();
		exit(0);

		dna_coupled_rotamer_design_test();
		exit(0);

		dna_deriv_test();
		exit(0);

		dna_design_test();
		exit(0);

		{ // loops graphics
			protocols::viewer::viewer_main( simple_loop_modeling_test_wrapper );
			exit( 0 );
		}

		simple_loop_modeling_test();
		exit(0);
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
