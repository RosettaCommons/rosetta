// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/Residue.cxxtest.hh
/// @brief  test suite for core::conformation::Residue
/// @author Christopher Miles (cmiles@uw.edu)

// Test Headers
#include <cxxtest/TestSuite.h>

// Unit Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/PseudoBond.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>

// Project headers
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/simple_moves/ModifyVariantTypeMover.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Basic Headers
#include <basic/Tracer.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.conformation.Residue.cxxtest" );

using core::pose::Pose;
using core::conformation::PseudoBond;
using core::conformation::PseudoBondCollection;
using core::conformation::PseudoBondCollectionOP;
using core::conformation::PseudoBondCollectionCOP;
using core::conformation::Residue;
using core::conformation::ResidueOP;

class ResidueTest : public CxxTest::TestSuite {
public:
	void setUp() {
		core_init_with_additional_options("-extra_res_fa core/chemical/params/U05.params");
	}

	void test_isDNA() {
		Pose pose;
		core::import_pose::pose_from_file(pose, "core/conformation/4gatA.pdb", core::import_pose::PDB_file);


		unsigned dna_start = 68;
		unsigned dna_end = 93;

		for ( unsigned i = 1; i < dna_start; ++i ) {
			TS_ASSERT(!pose.residue(i).is_DNA());
		}

		for ( unsigned i = dna_start; i <= dna_end; ++i ) {
			TS_ASSERT(pose.residue(i).is_DNA());
		}
	}

	void test_isLigand() {
		Pose pose;
		core::import_pose::pose_from_file(pose, "core/conformation/4gatA.pdb", core::import_pose::PDB_file);
		TS_ASSERT(pose.residue(67).is_ligand());
	}

	void test_recursive_identification_of_connection_dependencies() {
		//Note -- the determination is no longer recursive.  VKM -- 27 Oct 2017.
		Pose pose;
		core::import_pose::pose_from_file(pose, "core/conformation/4gatA.pdb", core::import_pose::PDB_file);

		protocols::simple_moves::ModifyVariantTypeMover add_nmethyl;
		add_nmethyl.set_additional_type_to_add("N_METHYLATION");
		core::select::residue_selector::ResidueIndexSelectorOP selector( new core::select::residue_selector::ResidueIndexSelector );
		selector->append_index( 5 );
		add_nmethyl.set_residue_selector(selector);
		add_nmethyl.apply(pose);

		TR << "\nATOM\tLOWER_DEP\tUPPER_DEP\n";
		for ( core::Size ia=1, iamax=pose.residue(5).natoms(); ia<=iamax; ++ia ) {
			std::string const atomname( pose.residue(5).atom_name(ia) );
			bool const lowerdep( pose.residue(5).atom_depends_on_lower(ia) );
			//bool const upperdep(false);
			bool const upperdep( pose.residue(5).atom_depends_on_upper(ia) );
			TR << atomname << "\t" << (lowerdep ? "TRUE" : "FALSE") << "\t" << (upperdep ? "TRUE" : "FALSE") << "\n";
			if ( atomname == " CN " || atomname == "1HN " || atomname == "2HN " || atomname == "3HN " ) {
				TS_ASSERT( lowerdep );
			} else {
				TS_ASSERT( !lowerdep );
			}
			if ( atomname == " O  " ) { TS_ASSERT(upperdep); }
			else { TS_ASSERT(!upperdep); }
		}
		TR << std::endl;
		TR.flush();
	}

	void test_recursive_identification_of_connection_dependencies_peptoid() {
		//Note -- the determination is no longer recursive.  VKM -- 27 Oct 2017.
		Pose pose;
		core::import_pose::pose_from_file(pose, "core/conformation/4gatA.pdb", core::import_pose::PDB_file);

		protocols::simple_moves::MutateResidue mut;
		mut.set_target(5);
		mut.set_res_name("601"); //A peptoid.
		mut.apply(pose);

		TR << "\nATOM\tLOWER_DEP\tUPPER_DEP\n";
		for ( core::Size ia=1, iamax=pose.residue(5).natoms(); ia<=iamax; ++ia ) {
			std::string const atomname( pose.residue(5).atom_name(ia) );
			bool const lowerdep( pose.residue(5).atom_depends_on_lower(ia) );
			//bool const upperdep(false);
			bool const upperdep( pose.residue(5).atom_depends_on_upper(ia) );
			TR << atomname << "\t" << (lowerdep ? "TRUE" : "FALSE") << "\t" << (upperdep ? "TRUE" : "FALSE") << "\n";
			if ( atomname == " CA1" || atomname == " CB1" || atomname == " CB2" || atomname == " CG1" || atomname == " CG2" ||
					atomname == " CD1" || atomname == " CD2" || atomname == " CE " || atomname == "1HA1" ||
					atomname == "1HB1" || atomname == "2HB1" || atomname == "3HB1" || atomname == "1HG1" ||
					atomname == "1HG2" || atomname == "1HD1" || atomname == "1HD2" || atomname == "1HE "
					) {
				TS_ASSERT( lowerdep );
			} else {
				TS_ASSERT( !lowerdep );
			}
			if ( atomname == " O  " ) { TS_ASSERT(upperdep); }
			else { TS_ASSERT(!upperdep); }
		}
		TR << std::endl;
		TR.flush();
	}

	void test_residue_serialization() {
		TS_ASSERT( true );
#ifdef SERIALIZATION
		Pose trpcage = create_trpcage_ideal_pose();
		ResidueOP trp6 = trpcage.residue( 6 ).clone();

		// Let's pretend that residue 7 is only a single atom large so we
		// can add some pseudobonds between residue 6 and residue 8, and
		// then test whether pseudobonds are properly serialized.
		PseudoBond pb68;
		pb68.lr( 6 ); pb68.ur( 8 );
		pb68.lr_conn_id( 2 ); // connecting at C on residue 6
		pb68.ur_conn_id( 1 ); // connecting at N on residue 8
		PseudoBondCollectionOP pbc( new PseudoBondCollection );
		pbc->push_back( pb68 );
		trp6->set_pseudobonds_to_residue( 8, pbc );

		TS_ASSERT( trp6->aa() == core::chemical::aa_trp ); // this should be trp in the input pose; not a test, so much as an assertion.

		// Now serialize the coordinates
		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arch( oss );
			arch( trp6 );
		}

		ResidueOP trp6_copy;
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arch( iss );
			arch( trp6_copy );
		}

		// trp6_copy should be exactly the same as trp6
		// and it should point to the same (global) instance of
		// the fullatom tryptophan residue type
		TS_ASSERT( & trp6->type() == & trp6_copy->type() );

		// The code below walks through the data members of class Residue in order

		// Coordinates and dihedrals should match exactly
		for ( core::Size ii = 1; ii <= trp6->natoms(); ++ii ) {
			for ( core::Size jj = 1; jj <= 3; ++jj ) {
				TS_ASSERT_EQUALS( trp6->xyz( ii )( jj ), trp6_copy->xyz( ii )( jj ) );
			}
			TS_ASSERT_EQUALS( trp6->atom( ii ).type(),    trp6_copy->atom( ii ).type() );
			TS_ASSERT_EQUALS( trp6->atom( ii ).mm_type(), trp6_copy->atom( ii ).mm_type() );
		}

		// Sadly, skipping the orbitals since I'm unsure how to activate them (I'm sure it's easy).
		// ADD CODE FOR ORBITALS

		TS_ASSERT_EQUALS( trp6->seqpos(), trp6_copy->seqpos() );
		TS_ASSERT_EQUALS( trp6->chain(),  trp6_copy->chain() );

		TS_ASSERT_EQUALS( trp6->chi().size(), trp6_copy->chi().size() );
		for ( core::Size ii = 1; ii <= trp6->chi().size(); ++ii ) {
			TS_ASSERT_EQUALS( trp6->chi()[ ii ], trp6_copy->chi()[ ii ] );
		}

		// now, this comparison isn't going to be all that useful because
		// trp doesn't have any nu dihedrals, but let's include the code
		// anyways so it could be compared in the future with some other
		// nu-containing residue
		TS_ASSERT_EQUALS( trp6->nus().size(), trp6_copy->nus().size() );
		for ( core::Size ii = 1; ii <= trp6->nus().size(); ++ii ) {
			TS_ASSERT_EQUALS( trp6->nus()[ ii ], trp6_copy->nus()[ ii ] );
		}

		TS_ASSERT_EQUALS( trp6->mainchain_torsions().size(), trp6_copy->mainchain_torsions().size() );
		for ( core::Size ii = 1; ii <= trp6->mainchain_torsions().size(); ++ii ) {
			TS_ASSERT_EQUALS( trp6->mainchain_torsions()[ ii ], trp6_copy->mainchain_torsions()[ ii ] );
		}

		for ( core::Size ii = 1; ii <= 3; ++ii ) {
			TS_ASSERT_EQUALS( trp6->actcoord()( ii ), trp6_copy->actcoord()( ii ) );
		}

		// ADD CODE HERE TO COMPARE DATA CACHE!

		// huh... can't directly access the nonstandard_polymer_ data member.

		// Make sure that the inter-residue connection information has been perfectly transfered
		// Connect map
		for ( core::Size ii = 1; ii <= trp6->n_possible_residue_connections(); ++ii ) {
			TS_ASSERT_EQUALS( trp6->actual_residue_connection( ii ).resid(), trp6_copy->actual_residue_connection( ii ).resid() );
			TS_ASSERT_EQUALS( trp6->actual_residue_connection( ii ).connid(), trp6_copy->actual_residue_connection( ii ).connid() );
		}

		TS_ASSERT_EQUALS( trp6->connections_to_residue( 5 ), trp6_copy->connections_to_residue( 5 ) );
		TS_ASSERT_EQUALS( trp6->connections_to_residue( 7 ), trp6_copy->connections_to_residue( 7 ) );
		for ( core::Size ii = 1; ii <= trpcage.size(); ++ii ) {
			TS_ASSERT_EQUALS( trp6_copy->is_bonded( ii ), ii == 5 || ii == 7 ); // the only inter-residue bonds are to 5 and 7
			TS_ASSERT_EQUALS( trp6_copy->is_pseudo_bonded( ii ), ii == 8 ); // the only pseudobond is to residue 8
		}

		PseudoBondCollectionCOP trp6_pbs      = trp6->get_pseudobonds_to_residue( 8 );
		PseudoBondCollectionCOP trp6_copy_pbs = trp6_copy->get_pseudobonds_to_residue( 8 );

		TS_ASSERT_EQUALS( trp6_pbs->size(), trp6_copy_pbs->size() );
		for ( PseudoBondCollection::PBIter trp6_pb = trp6_pbs->iter_begin(),
				trp6_copy_pb = trp6_copy_pbs->iter_begin();
				trp6_pb != trp6_pbs->iter_end(); ++trp6_pb, ++trp6_copy_pb ) {
			TS_ASSERT_EQUALS( trp6_pb->lr(), trp6_copy_pb->lr() );
			TS_ASSERT_EQUALS( trp6_pb->ur(), trp6_copy_pb->ur() );
			TS_ASSERT_EQUALS( trp6_pb->lr_conn_id(), trp6_copy_pb->lr_conn_id() );
			TS_ASSERT_EQUALS( trp6_pb->ur_conn_id(), trp6_copy_pb->ur_conn_id() );
		}

#endif // SERIALIZATION
	}

	void test_fill_missing_atoms() {
		core::chemical::ChemicalManager & chemman( *core::chemical::ChemicalManager::get_instance() );
		core::chemical::ResidueType const & restype( chemman.residue_type_set(core::chemical::FA_STANDARD)->name_map( "U05" ) );
		core::conformation::ResidueOP orig_res( new core::conformation::Residue(restype, true) );
		core::Size root_index( restype.atom_index( restype.root_atom() ) ); // C21
		core::Size stub2_index( restype.atom_index( "C13" ) );
		core::Size stub3_index( restype.atom_index( "C4" ) );


		// For a conformation
		core::pose::PoseOP pose( new core::pose::Pose );
		pose->append_residue_by_jump( *orig_res, 1 );
		//pose->dump_pdb( std::cout );

		TR << "All but stub rebuilding" << std::endl;
		// Should be able to rebuild everything given the three atom stub
		{
			core::conformation::Residue residue( pose->residue(1) );
			utility::vector1<bool> missing(residue.natoms(),true);
			missing[ root_index ] = false;
			missing[ stub2_index ] = false;
			missing[ stub3_index ] = false;
			for ( core::Size ii(1); ii <= orig_res->natoms(); ++ii ) {
				if ( ii != root_index && ii != stub2_index && ii != stub3_index ) {
					residue.set_xyz(ii, core::Vector( 0, 0, 0 ) );
				}
			}
			TSM_ASSERT_THROWS_NOTHING( "All rebuild", residue.fill_missing_atoms(missing, pose->conformation()) );
			for ( core::Size ii(1); ii <= orig_res->natoms(); ++ii ) {
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " X", residue.xyz(ii).x(), orig_res->xyz(ii).x(), 0.0001 );
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " Y", residue.xyz(ii).y(), orig_res->xyz(ii).y(), 0.0001 );
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " Z", residue.xyz(ii).z(), orig_res->xyz(ii).z(), 0.0001 );
			}
		}


		TR << "Single Atom rebuilding" << std::endl;
		// Can we rebuild each single atom? (This includes the root)
		for ( core::Size ii(1); ii <= orig_res->natoms(); ++ii ) {
			core::conformation::Residue residue( pose->residue(1) );
			utility::vector1<bool> missing(residue.natoms(),false);
			residue.set_xyz(ii, core::Vector( 0, 0, 0 ) );
			missing[ii] = true;
			TSM_ASSERT_THROWS_NOTHING( "Single missing: " + restype.atom_name(ii), residue.fill_missing_atoms(missing, pose->conformation()) );
			TSM_ASSERT_DELTA( restype.atom_name(ii) + " X", residue.xyz(ii).x(), orig_res->xyz(ii).x(), 0.0001 );
			TSM_ASSERT_DELTA( restype.atom_name(ii) + " Y", residue.xyz(ii).y(), orig_res->xyz(ii).y(), 0.0001 );
			TSM_ASSERT_DELTA( restype.atom_name(ii) + " Z", residue.xyz(ii).z(), orig_res->xyz(ii).z(), 0.0001 );
		}

		TR << "Atom and connection rebuilding" << std::endl;
		// Can we rebuild atoms when they and all their connecting atoms are rebuilt?
		for ( core::Size ii(1); ii <= orig_res->natoms(); ++ii ) {
			core::conformation::Residue residue( pose->residue(1) );
			utility::vector1<bool> missing(residue.natoms(),false);
			residue.set_xyz(ii, core::Vector( 0, 0, 0 ) );
			missing[ii] = true;
			core::chemical::AtomIndices nbrs( residue.bonded_neighbor( ii ) );
			for ( core::Size jj(1); jj <= nbrs.size(); ++jj ) {
				missing[nbrs[jj]] = true;
				residue.set_xyz(nbrs[jj], core::Vector( 0, 0, 0 ) );
			}
			TSM_ASSERT_THROWS_NOTHING( "Connects missing: " + restype.atom_name(ii), residue.fill_missing_atoms(missing, pose->conformation()) );
			TSM_ASSERT_DELTA( restype.atom_name(ii) + " X", residue.xyz(ii).x(), orig_res->xyz(ii).x(), 0.0001 );
			TSM_ASSERT_DELTA( restype.atom_name(ii) + " Y", residue.xyz(ii).y(), orig_res->xyz(ii).y(), 0.0001 );
			TSM_ASSERT_DELTA( restype.atom_name(ii) + " Z", residue.xyz(ii).z(), orig_res->xyz(ii).z(), 0.0001 );
		}

		TR << "Distal stub rebuild" << std::endl;
		// We have three connected atoms out in the distance - can we rebuild everything?
		{
			core::Size atom1( restype.atom_index( "C16" ) ), atom2( restype.atom_index( "C23" ) ), atom3( restype.atom_index( "N3" ) );
			core::conformation::Residue residue( pose->residue(1) );
			utility::vector1<bool> missing(residue.natoms(),true);
			missing[ atom1 ] = false;
			missing[ atom2 ] = false;
			missing[ atom3 ] = false;
			for ( core::Size ii(1); ii <= orig_res->natoms(); ++ii ) {
				if ( ii != atom1 && ii != atom2 && ii != atom3 ) {
					residue.set_xyz(ii, core::Vector( 0, 0, 0 ) );
				}
			}
			TSM_ASSERT_THROWS_NOTHING( "Distal rebuild", residue.fill_missing_atoms(missing, pose->conformation()) );
			for ( core::Size ii(1); ii <= orig_res->natoms(); ++ii ) {
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " X", residue.xyz(ii).x(), orig_res->xyz(ii).x(), 0.0001 );
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " Y", residue.xyz(ii).y(), orig_res->xyz(ii).y(), 0.0001 );
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " Z", residue.xyz(ii).z(), orig_res->xyz(ii).z(), 0.0001 );
			}
		}

		TR << "Distal stub non-icoord rebuild" << std::endl;
		// We have three atoms but not along icoord path - can we rebuild everything?
		{
			core::Size atom1( restype.atom_index( "C24" ) ), atom2( restype.atom_index( "C25" ) ), atom3( restype.atom_index( "C26" ) ); // Not in ICOORD path.
			core::conformation::Residue residue( pose->residue(1) );
			utility::vector1<bool> missing(residue.natoms(),true);
			missing[ atom1 ] = false;
			missing[ atom2 ] = false;
			missing[ atom3 ] = false;
			for ( core::Size ii(1); ii <= orig_res->natoms(); ++ii ) {
				if ( ii != atom1 && ii != atom2 && ii != atom3 ) {
					residue.set_xyz(ii, core::Vector( 0, 0, 0 ) );
				}
			}
			TSM_ASSERT_THROWS_NOTHING( "Distal rebuild", residue.fill_missing_atoms(missing, pose->conformation()) );
			for ( core::Size ii(1); ii <= orig_res->natoms(); ++ii ) {
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " X", residue.xyz(ii).x(), orig_res->xyz(ii).x(), 0.0001 );
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " Y", residue.xyz(ii).y(), orig_res->xyz(ii).y(), 0.0001 );
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " Z", residue.xyz(ii).z(), orig_res->xyz(ii).z(), 0.0001 );
			}
		}
	}

	void test_fill_missing_atoms_polymer() {
		core::chemical::ChemicalManager & chemman( *core::chemical::ChemicalManager::get_instance() );
		core::chemical::ResidueType const & restype( chemman.residue_type_set(core::chemical::FA_STANDARD)->name_map( "GLU" ) );
		core::conformation::ResidueOP orig_res( new core::conformation::Residue(restype, true) );

		// Form a conformation -- Note that the single residue is *not* a termini variant
		core::pose::PoseOP pose( new core::pose::Pose );
		pose->append_residue_by_jump( *orig_res, 1 );

		// Backbone rebuild --
		TR << "Single Residue backbone rebuild." << std::endl;
		{
			core::conformation::Residue residue( pose->residue(1) );
			utility::vector1<bool> missing(residue.natoms(),false);
			missing[ restype.atom_index( "N" ) ] = true;
			missing[ restype.atom_index( "H" ) ] = true;
			missing[ restype.atom_index( "CA" ) ] = true;
			missing[ restype.atom_index( "C" ) ] = true;
			missing[ restype.atom_index( "O" ) ] = true;
			TSM_ASSERT_THROWS_NOTHING( "Backbone rebuild", residue.fill_missing_atoms(missing, pose->conformation()) );
			for ( core::Size ii(1); ii <= orig_res->natoms(); ++ii ) {
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " X", residue.xyz(ii).x(), orig_res->xyz(ii).x(), 0.0001 );
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " Y", residue.xyz(ii).y(), orig_res->xyz(ii).y(), 0.0001 );
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " Z", residue.xyz(ii).z(), orig_res->xyz(ii).z(), 0.0001 );
			}
		}

		// Now make a 3 residue pose
		TR << "Polymeric Residue backbone rebuild." << std::endl;
		pose->append_residue_by_bond( *orig_res, true );
		pose->append_residue_by_bond( *orig_res, true );
		{
			core::conformation::Residue residue( pose->residue(2) ); // Copy, so as not to modify the pose
			utility::vector1<bool> missing(residue.natoms(),false);
			missing[ restype.atom_index( "N" ) ] = true;
			missing[ restype.atom_index( "H" ) ] = true;
			missing[ restype.atom_index( "CA" ) ] = true;
			missing[ restype.atom_index( "C" ) ] = true;
			missing[ restype.atom_index( "O" ) ] = true;
			TSM_ASSERT_THROWS_NOTHING( "Backbone rebuild", residue.fill_missing_atoms(missing, pose->conformation()) );
			for ( core::Size ii(1); ii <= residue.natoms(); ++ii ) {
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " X", residue.xyz(ii).x(), pose->residue(2).xyz(ii).x(), 0.0001 );
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " Y", residue.xyz(ii).y(), pose->residue(2).xyz(ii).y(), 0.0001 );
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " Z", residue.xyz(ii).z(), pose->residue(2).xyz(ii).z(), 0.0001 );
			}
		}

		// Now make a 3 "chain" pose with cutpoints between each residue (again, none are termini/cutpoint variants)
		core::kinematics::FoldTree ft( pose->fold_tree() );
		ft.new_jump( 1, 2, 1 ); // make chainbreak between 1 and 2
		ft.new_jump( 2, 3, 2 ); // make chainbreak between 2 and 3
		pose->fold_tree( ft );
		//pose->dump_pdb( std::cout );

		TS_ASSERT( pose->fold_tree().is_cutpoint(1) );
		TS_ASSERT( pose->fold_tree().is_cutpoint(2) );

		TR << "Polymeric Residue cutpoint backbone rebuild." << std::endl;
		pose->append_residue_by_bond( *orig_res, true );
		pose->append_residue_by_bond( *orig_res, true );
		{
			core::conformation::Residue residue( pose->residue(2) ); // Copy, so as not to modify the pose
			utility::vector1<bool> missing(residue.natoms(),false);
			missing[ restype.atom_index( "N" ) ] = true;
			missing[ restype.atom_index( "H" ) ] = true;
			missing[ restype.atom_index( "CA" ) ] = true;
			missing[ restype.atom_index( "C" ) ] = true;
			missing[ restype.atom_index( "O" ) ] = true;
			TSM_ASSERT_THROWS_NOTHING( "Backbone rebuild", residue.fill_missing_atoms(missing, pose->conformation()) );
			for ( core::Size ii(1); ii <= residue.natoms(); ++ii ) {
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " X", residue.xyz(ii).x(), pose->residue(2).xyz(ii).x(), 0.0001 );
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " Y", residue.xyz(ii).y(), pose->residue(2).xyz(ii).y(), 0.0001 );
				TSM_ASSERT_DELTA( restype.atom_name(ii) + " Z", residue.xyz(ii).z(), pose->residue(2).xyz(ii).z(), 0.0001 );
			}
		}

	}
};

