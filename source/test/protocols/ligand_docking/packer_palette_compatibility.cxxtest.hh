// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// @brief This file is to test that the behavior of PackerPalettes with non-standard ResidueTypes
// behaves as one would expect.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/PoseResidueTypeSet.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/residue_io.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pack/palette/CustomBaseTypePackerPalette.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/packer_neighbors.hh>

#include <protocols/ligand_docking/HighResDocker.hh>

#include <core/types.hh>

#include <utility/graph/Graph.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.ligand_docking.packer_palette_compatibility.cxxtest");

const core::Size AO_ROTAMERS = 14;
const core::Size AP_ROTAMERS = 12;
const core::Size BG_ROTAMERS = 30;
const core::Size BH_ROTAMERS = 28;

/// @brief A utility class to allow access to the make_packer_task of the HighResDocker
class PackFromDocker: public protocols::ligand_docking::HighResDocker {
public:
	using protocols::ligand_docking::HighResDocker::make_packer_task;
};

core::pose::PoseOP
make_pose( core::chemical::ResidueType const & rt1, core::chemical::ResidueType const & rt2 ) {
	core::pose::PoseOP pose = create_twores_1ubq_poseop();
	core::conformation::Residue res1( rt1, true );
	core::conformation::Residue res2( rt2, true );
	// make sure the residues are well separated from each other, so we don't have clash issues
	for ( core::Size ii(1); ii <= res1.natoms(); ++ii ) {
		res1.set_xyz( ii, res1.xyz(ii) + core::Vector(100,0,0) );
	}
	for ( core::Size ii(1); ii <= res2.natoms(); ++ii ) {
		res2.set_xyz( ii, res1.xyz(ii) + core::Vector(0,0,100) );
	}
	pose->append_residue_by_jump( res1, 1 );
	pose->append_residue_by_jump( res2, 1 );

	return pose;
}

core::pose::PoseOP
make_pose_cached( core::chemical::MutableResidueTypeOP mrt1, core::chemical::MutableResidueTypeOP mrt2 ) {
	using namespace core::chemical;
	core::pose::PoseOP pose = create_twores_1ubq_poseop();

	PoseResidueTypeSetOP pose_rts = pose->conformation().modifiable_residue_type_set_for_conf();
	pose_rts->add_unpatchable_residue_type( mrt1 ); // The  test_inpose_custom_palette & test_inpose_highresdocker tests don't work if this is a base type
	pose_rts->add_unpatchable_residue_type( mrt2 ); // ibid.
	pose->conformation().reset_residue_type_set_for_conf( pose_rts );

	core::chemical::ResidueType const & rt1 = pose->residue_type_set_for_pose()->name_map( mrt1->name() );
	core::chemical::ResidueType const & rt2 = pose->residue_type_set_for_pose()->name_map( mrt2->name() );
	core::conformation::Residue res1( rt1, true );
	core::conformation::Residue res2( rt2, true );
	// make sure the residues are well separated from each other, so we don't have clash issues
	for ( core::Size ii(1); ii <= res1.natoms(); ++ii ) {
		res1.set_xyz( ii, res1.xyz(ii) + core::Vector(100,0,0) );
	}
	for ( core::Size ii(1); ii <= res2.natoms(); ++ii ) {
		res2.set_xyz( ii, res1.xyz(ii) + core::Vector(0,0,100) );
	}
	pose->append_residue_by_jump( res1, 1 );
	pose->append_residue_by_jump( res2, 1 );

	return pose;
}

core::pack::rotamer_set::RotamerSetsOP
build_rotsets( core::pose::Pose & pose, core::pack::task::PackerTaskCOP task ) {
	using namespace core::pack;

	pose.update_residue_neighbors();

	core::scoring::ScoreFunction scfxn; // Default, empty scorefunction should be fine.
	scfxn.setup_for_packing( pose, task->repacking_residues(), task->designing_residues() );

	utility::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, scfxn, task );

	rotamer_set::RotamerSetsOP rotsets = utility::pointer::make_shared< rotamer_set::RotamerSets >();
	rotsets->set_task( task );
	rotsets->initialize_pose_for_rotsets_creation(pose);

	rotsets->build_rotamers( pose, scfxn, packer_neighbor_graph );

	return rotsets;
}

// Attempt to setup packing for two ligands in a pose.
// The residue types can either be with the same interchangability group or with two different ones
// The residue types can either be in the global RTS, in the Pose's RTS, or orphan (not in an RTS)
// Check to make sure that the identities and the rotamers that we get are appropriate.

class packer_palette_compatibility_global_tests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init_with_additional_options("-extra_res_fa protocols/ligand_docking/AO.params protocols/ligand_docking/AP.params protocols/ligand_docking/BG.params protocols/ligand_docking/BH.params");
	}

	void tearDown() {}

	void test_different_interchange_group() {
		using namespace core::chemical;
		using namespace core::pack;

		auto rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t );

		core::pose::PoseOP pose = make_pose( rts->name_map("AO"), rts->name_map("AP") );

		TS_ASSERT_EQUALS( pose->size(), 4 );

		task::PackerTaskOP task = task::TaskFactory::create_packer_task( *pose );
		rotamer_set::RotamerSetsOP rotsets = build_rotsets( *pose, task );

		TS_ASSERT_EQUALS( rotsets->nmoltenres(), pose->size() );
		TS_ASSERT_EQUALS( pose->residue_type(3).name(), "AO" );
		TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres(3), AO_ROTAMERS );
		TS_ASSERT_EQUALS( pose->residue_type(4).name(), "AP" );
		TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres(4), AP_ROTAMERS );
	}

	void test_same_interchange_group() {
		using namespace core::chemical;
		using namespace core::pack;

		auto rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t );

		core::pose::PoseOP pose = make_pose( rts->name_map("BG"), rts->name_map("BH") );

		TS_ASSERT_EQUALS( pose->size(), 4 );

		task::PackerTaskOP task = task::TaskFactory::create_packer_task( *pose );
		rotamer_set::RotamerSetsOP rotsets = build_rotsets( *pose, task );

		TS_ASSERT_EQUALS( rotsets->nmoltenres(), pose->size() );
		TS_ASSERT_EQUALS( pose->residue_type(3).name(), "BG" );
		TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres(3), BG_ROTAMERS );
		TS_ASSERT_EQUALS( pose->residue_type(4).name(), "BH" );
		TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres(4), BH_ROTAMERS );
	}

};

class packer_palette_compatibility_tests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_inpose_different_interchange_group() {
		using namespace core::chemical;
		using namespace core::pack;

		auto global_rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t );
		MutableResidueTypeOP res1 = read_topology_file( "protocols/ligand_docking/AO.params", global_rts );
		MutableResidueTypeOP res2 = read_topology_file( "protocols/ligand_docking/AP.params", global_rts );
		TS_ASSERT( ! global_rts->has_name("AO") );
		TS_ASSERT( ! global_rts->has_name("AP") );

		core::pose::PoseOP pose = make_pose_cached( res1, res2 );

		task::PackerTaskOP task = task::TaskFactory::create_packer_task( *pose );
		rotamer_set::RotamerSetsOP rotsets = build_rotsets( *pose, task );

		TS_ASSERT_EQUALS( rotsets->nmoltenres(), pose->size() );
		TS_ASSERT_EQUALS( pose->residue_type(3).name(), "AO" );
		TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres(3), AO_ROTAMERS );
		TS_ASSERT_EQUALS( pose->residue_type(4).name(), "AP" );
		TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres(4), AP_ROTAMERS );
	}

	void test_inpose_same_interchange_group() {
		using namespace core::chemical;
		using namespace core::pack;

		auto global_rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t );
		MutableResidueTypeOP res1 = read_topology_file( "protocols/ligand_docking/BH.params", global_rts );
		MutableResidueTypeOP res2 = read_topology_file( "protocols/ligand_docking/BG.params", global_rts );
		TS_ASSERT( ! global_rts->has_name("BH") );
		TS_ASSERT( ! global_rts->has_name("BG") );

		core::pose::PoseOP pose = make_pose_cached( res1, res2 );

		task::PackerTaskOP task = task::TaskFactory::create_packer_task( *pose );
		rotamer_set::RotamerSetsOP rotsets = build_rotsets( *pose, task );

		TS_ASSERT_EQUALS( rotsets->nmoltenres(), pose->size() );
		TS_ASSERT_EQUALS( pose->residue_type(3).name(), "BH" );
		TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres(3), BH_ROTAMERS );
		TS_ASSERT_EQUALS( pose->residue_type(4).name(), "BG" );
		TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres(4), BG_ROTAMERS );
	}

	void test_orphan_different_interchange_group() {
		using namespace core::chemical;
		using namespace core::pack;

		auto global_rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t );
		ResidueTypeCOP res1 = ResidueType::make( *read_topology_file( "protocols/ligand_docking/AO.params", global_rts ) );
		ResidueTypeCOP res2 = ResidueType::make( *read_topology_file( "protocols/ligand_docking/AP.params", global_rts ) );
		TS_ASSERT( ! global_rts->has_name("AO") );
		TS_ASSERT( ! global_rts->has_name("AP") );

		core::pose::PoseOP pose = make_pose( *res1, *res2 );

		TS_ASSERT_EQUALS( pose->size(), 4 );

		task::PackerTaskOP task = task::TaskFactory::create_packer_task( *pose );
		rotamer_set::RotamerSetsOP rotsets = build_rotsets( *pose, task );

		TS_ASSERT_EQUALS( rotsets->nmoltenres(), pose->size() );
		TS_ASSERT_EQUALS( pose->residue_type(3).name(), "AO" );
		TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres(3), AO_ROTAMERS );
		TS_ASSERT_EQUALS( pose->residue_type(4).name(), "AP" );
		TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres(4), AP_ROTAMERS );
	}

	void test_orphan_same_interchange_group() {
		using namespace core::chemical;
		using namespace core::pack;

		auto global_rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t );
		ResidueTypeCOP res1 = ResidueType::make( *read_topology_file( "protocols/ligand_docking/BG.params", global_rts ) );
		ResidueTypeCOP res2 = ResidueType::make( *read_topology_file( "protocols/ligand_docking/BH.params", global_rts ) );
		TS_ASSERT( ! global_rts->has_name("BG") );
		TS_ASSERT( ! global_rts->has_name("BH") );

		core::pose::PoseOP pose = make_pose( *res1, *res2 );

		TS_ASSERT_EQUALS( pose->size(), 4 );

		task::PackerTaskOP task = task::TaskFactory::create_packer_task( *pose );
		rotamer_set::RotamerSetsOP rotsets = build_rotsets( *pose, task );

		TS_ASSERT_EQUALS( rotsets->nmoltenres(), pose->size() );
		TS_ASSERT_EQUALS( pose->residue_type(3).name(), "BG" );
		TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres(3), BG_ROTAMERS );
		TS_ASSERT_EQUALS( pose->residue_type(4).name(), "BH" );
		TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres(4), BH_ROTAMERS );
	}

	void test_inpose_custom_palette() {
		// This should match what the protocols::ligand_docking::HighResDocker::make_packer_task_from_vector() is doing
		using namespace core::chemical;
		using namespace core::pack;

		auto global_rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t );
		MutableResidueTypeOP res1 = read_topology_file( "protocols/ligand_docking/BH.params", global_rts );
		MutableResidueTypeOP res2 = read_topology_file( "protocols/ligand_docking/BG.params", global_rts );
		TS_ASSERT( ! global_rts->has_name("BH") );
		TS_ASSERT( ! global_rts->has_name("BG") );

		core::pose::PoseOP pose = make_pose_cached( res1, res2 );

		palette::CustomBaseTypePackerPaletteOP palette( utility::pointer::make_shared< palette::CustomBaseTypePackerPalette >() );

		palette->add_type( res1->base_name() );
		palette->add_type( res2->base_name() );

		task::PackerTaskOP task = task::TaskFactory::create_packer_task( *pose, palette );

		// Allow design
		utility::vector1< std::string > types_to_allow = { res1->base_name(), res2->base_name() };
		task->nonconst_residue_task(3).restrict_restypes( types_to_allow );
		task->nonconst_residue_task(4).restrict_restypes( types_to_allow );

		rotamer_set::RotamerSetsOP rotsets = build_rotsets( *pose, task );

		TS_ASSERT_EQUALS( rotsets->nmoltenres(), pose->size() );
		TS_ASSERT_EQUALS( pose->residue_type(3).name(), "BH" );
		TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres(3), BH_ROTAMERS + BG_ROTAMERS );
		TS_ASSERT_EQUALS( pose->residue_type(4).name(), "BG" );
		TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres(4), BG_ROTAMERS + BH_ROTAMERS);
	}

	void test_inpose_highresdocker() {
		// Explicitly call protocols::ligand_docking::HighResDocker::make_packer_task()
		using namespace core::chemical;
		using namespace core::pack;

		auto global_rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t );
		MutableResidueTypeOP res1 = read_topology_file( "protocols/ligand_docking/BH.params", global_rts );
		MutableResidueTypeOP res2 = read_topology_file( "protocols/ligand_docking/BG.params", global_rts );
		TS_ASSERT( ! global_rts->has_name("BH") );
		TS_ASSERT( ! global_rts->has_name("BG") );

		core::pose::PoseOP pose = make_pose_cached( res1, res2 );

		PackFromDocker pfd;

		task::PackerTaskOP task = pfd.make_packer_task( *pose, true );

		rotamer_set::RotamerSetsOP rotsets = build_rotsets( *pose, task );

		TS_ASSERT_EQUALS( rotsets->nmoltenres(), pose->size() );
		TS_ASSERT_EQUALS( pose->residue_type(3).name(), "BH" );
		TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres(3), BH_ROTAMERS + BG_ROTAMERS );
		TS_ASSERT_EQUALS( pose->residue_type(4).name(), "BG" );
		TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres(4), BG_ROTAMERS + BH_ROTAMERS);
	}

};

