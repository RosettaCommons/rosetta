// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody/AntibodyInfo.cxxtest.hh
/// @brief  tests for the AntibodyInfo class
/// @author Jared Adolf-Bryfogle
/// @author Jeff Gray


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <utility/vector1.hh>

// Protocol Headers
#include <basic/Tracer.hh>
#include <boost/foreach.hpp>

#define BFE BOOST_FOREACH
using namespace protocols::antibody;
using utility::vector1;

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.AntibodyInfoTest");
class AntibodyInfoTests : public CxxTest::TestSuite {
	core::pose::Pose ab_pose_aho; //Full PDB
	core::pose::Pose ab_pose_chothia;
	core::pose::Pose ab_pose_aho_antigen;
	AntibodyInfoOP ab_info_north_aho;
	AntibodyInfoOP ab_info_chothia;
	AntibodyInfoOP ab_info_aroop;
	AntibodyInfoOP ab_info_aho_antigen;

	std::map< core::Size,  std::map< core::Size,  vector1< core::Size > > > numbering;//[scheme][start/end][cdr][pose number].  Used to test numbering and transforms.
	typedef std::map<CDRDefinitionEnum, std::pair< core::pose::Pose, AntibodyInfoOP> >  AbInfos;
	AbInfos infos;

public:

	void setUp(){

		core_init();
		core::import_pose::pose_from_file(ab_pose_aho, "protocols/antibody/1bln_AB_aho.pdb", core::import_pose::PDB_file);
		core::import_pose::pose_from_file(ab_pose_chothia, "protocols/antibody/1bln_AB_chothia.pdb", core::import_pose::PDB_file);
		core::import_pose::pose_from_file(ab_pose_aho_antigen, "protocols/antibody/aho_with_antigen.pdb", core::import_pose::PDB_file);
		ab_info_north_aho = AntibodyInfoOP(new AntibodyInfo(ab_pose_aho, AHO_Scheme, North));
		ab_info_chothia = AntibodyInfoOP( new AntibodyInfo(ab_pose_chothia, Chothia_Scheme, Chothia));
		ab_info_aroop = AntibodyInfoOP( new AntibodyInfo(ab_pose_chothia, Chothia_Scheme, Aroop) );
		ab_info_aho_antigen = AntibodyInfoOP( new AntibodyInfo(ab_pose_aho_antigen, AHO_Scheme, North));

		infos[North] = std::make_pair(ab_pose_aho, ab_info_north_aho);
		infos[Chothia] = std::make_pair(ab_pose_chothia, ab_info_chothia);
		infos[Aroop] = std::make_pair(ab_pose_chothia, ab_info_aroop);

		numbering[North][cdr_start].resize(6);
		numbering[North][cdr_end].resize(6);

		numbering[Chothia][cdr_start].resize(6);
		numbering[Chothia][cdr_end].resize(6);

		numbering[Aroop][cdr_start].resize(6);
		numbering[Aroop][cdr_end].resize(6);

		//North Aho
		numbering[North][cdr_start][l1] = 24;
		numbering[North][cdr_start][l2] = 54;
		numbering[North][cdr_start][l3] = 94;
		numbering[North][cdr_start][h1] = 135;
		numbering[North][cdr_start][h2] = 162;
		numbering[North][cdr_start][h3] = 209;
		numbering[North][cdr_end][l1] = 39;
		numbering[North][cdr_end][l2] = 61;
		numbering[North][cdr_end][l3] = 102;
		numbering[North][cdr_end][h1] = 147;
		numbering[North][cdr_end][h2] = 171;
		numbering[North][cdr_end][h3] = 220;

		//Chothia
		numbering[Chothia][cdr_start][l1] = 24;
		numbering[Chothia][cdr_start][l2] = 55;
		numbering[Chothia][cdr_start][l3] = 94;
		numbering[Chothia][cdr_start][h1] = 138;
		numbering[Chothia][cdr_start][h2] =164;
		numbering[Chothia][cdr_start][h3] =211 ;
		numbering[Chothia][cdr_end][l1] = 39;
		numbering[Chothia][cdr_end][l2] = 61;
		numbering[Chothia][cdr_end][l3] = 102;
		numbering[Chothia][cdr_end][h1] = 144;
		numbering[Chothia][cdr_end][h2] = 169;
		numbering[Chothia][cdr_end][h3] = 220;

		//Aroop
		numbering[Aroop][cdr_start][l1] = 24;
		numbering[Aroop][cdr_start][l2] = 55;
		numbering[Aroop][cdr_start][l3] = 94;
		numbering[Aroop][cdr_start][h1] = 138;
		numbering[Aroop][cdr_start][h2] =162;
		numbering[Aroop][cdr_start][h3] =211 ;
		numbering[Aroop][cdr_end][l1] = 39;
		numbering[Aroop][cdr_end][l2] = 61;
		numbering[Aroop][cdr_end][l3] = 102;
		numbering[Aroop][cdr_end][h1] = 147;
		numbering[Aroop][cdr_end][h2] = 178;
		numbering[Aroop][cdr_end][h3] = 220;

	}

	void tearDown(){
		ab_pose_aho.clear();
		ab_pose_chothia.clear();
	}

	void test_cdr_numbering_and_transform(){



		BFE(AbInfos::value_type outer_info, infos){
			CDRDefinitionEnum outer_definition = outer_info.first;

			std::pair<core::pose::Pose, AntibodyInfoOP > outer_pair = outer_info.second;
			AntibodyInfoOP outer_ab_info = outer_pair.second;
			std::string outer_definition_str = outer_ab_info->get_current_CDRDefinition();

			for ( core::Size i = 1; i<= 6; ++i ) {
				CDRNameEnum cdr = static_cast<CDRNameEnum>(i);

				//TR<< "Testing: "<< outer_scheme_str << std::endl;
				TS_ASSERT_EQUALS(numbering[outer_definition][cdr_start][cdr], outer_ab_info->get_CDR_loop(cdr).start());
				TS_ASSERT_EQUALS(numbering[outer_definition][cdr_end][cdr], outer_ab_info->get_CDR_loop(cdr).stop());

				//Test each numbering transform
				BFE(AbInfos::value_type inner_info, infos){
					CDRDefinitionEnum inner_definition = inner_info.first;
					std::pair< core::pose::Pose, AntibodyInfoOP > inner_pair = inner_info.second;
					AntibodyInfoOP inner_ab_info = inner_pair.second;
					std::string inner_definition_str = inner_ab_info->get_current_CDRDefinition();
					//TR<< "Testing: "<< outer_definition_str << " to " << inner_definition_str << std::endl;
					TS_ASSERT_EQUALS(numbering[inner_definition][cdr_start][cdr], outer_ab_info->get_CDR_loop(cdr, outer_pair.first, inner_definition).start());
					TS_ASSERT_EQUALS(numbering[inner_definition][cdr_end][cdr], outer_ab_info->get_CDR_loop(cdr, outer_pair.first, inner_definition).stop());

					//Test Full transform
					//outer_ab_info->set_transform_cdr_definition(inner_scheme);
					//TS_ASSERT_EQUALS(numbering[inner_scheme][cdr_start][cdr], outer_ab_info->get_CDR_loop(cdr).start());
					//TS_ASSERT_EQUALS(numbering[inner_scheme][cdr_end][cdr], outer_ab_info->get_CDR_loop(cdr).stop());
					//outer_ab_info->clear_transform_cdr_definition();
				}
			}
		}

		//Here, we test Numbering scheme/Landmark access.

		core::Size chothia_num = ab_pose_chothia.pdb_info()->pdb2pose('L', 43);
		core::Size aho_num = ab_info_north_aho->get_landmark_resnum(ab_pose_aho, Chothia_Scheme, 'L', 43);
		TS_ASSERT_EQUALS(chothia_num, aho_num);

	}

void test_info_functions(){
	TS_ASSERT_EQUALS("L1", ab_info_north_aho->get_CDR_name(l1));
	TS_ASSERT_EQUALS("AHO_Scheme", ab_info_north_aho->get_current_AntibodyNumberingScheme());
	TS_ASSERT_EQUALS("North", ab_info_north_aho->get_current_CDRDefinition());

	TS_ASSERT_EQUALS("Chothia", ab_info_chothia->get_current_CDRDefinition());
	TS_ASSERT_EQUALS("Chothia_Scheme", ab_info_aroop->get_current_AntibodyNumberingScheme());
	TS_ASSERT_EQUALS("Aroop", ab_info_aroop->get_current_CDRDefinition());

	TS_ASSERT_EQUALS(false, ab_info_north_aho->antigen_present());
	TS_ASSERT_EQUALS(16, ab_info_north_aho->get_CDR_length(l1));
	TS_ASSERT_EQUALS('H', ab_info_north_aho->get_CDR_chain(h3));
	TS_ASSERT(! ab_info_north_aho->antigen_present());
	TS_ASSERT_EQUALS(24, ab_info_north_aho->get_CDR_start_PDB_num(l1));
	TS_ASSERT_EQUALS(CDRNameEnum_total, ab_info_north_aho->get_total_num_CDRs());
	TS_ASSERT_EQUALS(0, ab_info_north_aho->get_antigen_chains().size());

	//Test Camelid
	TS_ASSERT(! ab_info_chothia->is_camelid());

	//Test Chothia numbering with North definitions in cluster
	TS_ASSERT_EQUALS("L1-16-1", ab_info_chothia->get_cluster_name(ab_info_chothia->get_CDR_cluster(l1)->cluster()));
	TS_ASSERT_EQUALS("L2-8-1", ab_info_chothia->get_cluster_name(ab_info_chothia->get_CDR_cluster(l2)->cluster()));
	TS_ASSERT_EQUALS("L3-9-2", ab_info_chothia->get_cluster_name(ab_info_chothia->get_CDR_cluster(l3)->cluster()));

	//Regional tests:
	TR << "Testing regions " << std::endl;
	//Test correct identification of CDR regions
	for ( core::Size i = 1; i <= 6; ++i ) {
		CDRNameEnum cdr = static_cast< CDRNameEnum >(i);
		for ( core::Size cdr_res = ab_info_north_aho->get_CDR_start( cdr, ab_pose_aho ); cdr_res <= ab_info_north_aho->get_CDR_end( cdr, ab_pose_aho ); ++ cdr_res ) {
			TS_ASSERT_EQUALS(cdr_region, ab_info_north_aho->get_region_of_residue( ab_pose_aho, cdr_res ) );
		}
	}
	//Test correct identification of framework regions
	assert(ab_info_north_aho->antigen_present() == false);
	for ( core::Size i = 1; i <= ab_pose_aho.total_residue(); ++i ) {
		if ( ab_info_north_aho->get_region_of_residue( ab_pose_aho, i ) != cdr_region ) {
			TS_ASSERT_EQUALS( framework_region, ab_info_north_aho->get_region_of_residue( ab_pose_aho, i ) );
		}
	}
	//Test correct identification of antigen regions
	utility::vector1<core::Size> ids = ab_info_aho_antigen->get_antigen_chain_ids(ab_pose_aho_antigen);

	for ( core::Size i = 1; i <= ids.size(); ++i ) {
		core::Size res_start = ab_pose_aho_antigen.conformation().chain_begin( ids[ i ] );
		core::Size res_end = ab_pose_aho_antigen.conformation().chain_end( ids[ i ] );
		for ( core::Size res = res_start; res <= res_end; ++res ) {
			TS_ASSERT_EQUALS(antigen_region, ab_info_aho_antigen->get_region_of_residue( ab_pose_aho_antigen, res ) );
		}
	}
}
void test_kink_functions(){
	// Aroop
	TS_ASSERT_EQUALS(218,ab_info_aroop->kink_begin(ab_pose_chothia));
	TS_ASSERT_EQUALS(221,ab_info_aroop->kink_end(ab_pose_chothia));
	TS_ASSERT_EQUALS(221,ab_info_aroop->kink_trp(ab_pose_chothia));
	core::conformation::Residue Wres = ab_pose_chothia.residue(ab_info_aroop->kink_trp(ab_pose_chothia));
	core::conformation::Residue Ares = ab_pose_chothia.residue(ab_info_aroop->kink_anion_residue(ab_pose_chothia));
	core::conformation::Residue Cres = ab_pose_chothia.residue(ab_info_aroop->kink_cation_residue(ab_pose_chothia));
	TS_ASSERT_EQUALS("TRP",Wres.name3());
	TS_ASSERT_EQUALS("ALA",Ares.name3());  // (not an anion in this PDB, but Ala is in the right position)
	TS_ASSERT_EQUALS("ARG",Cres.name3());


	// Modified_AHO - Also a test of transform as kinks are defined in chothia/aroop and must transform numbering from north_aho
	TS_ASSERT_EQUALS(218,ab_info_north_aho->kink_begin(ab_pose_aho));
	TS_ASSERT_EQUALS(221,ab_info_north_aho->kink_end(ab_pose_aho));
	TS_ASSERT_EQUALS(221,ab_info_north_aho->kink_trp(ab_pose_aho));
	core::conformation::Residue WresAHO = ab_pose_aho.residue(ab_info_north_aho->kink_trp(ab_pose_aho));
	core::conformation::Residue AresAHO = ab_pose_aho.residue(ab_info_north_aho->kink_anion_residue(ab_pose_aho));
	core::conformation::Residue CresAHO = ab_pose_aho.residue(ab_info_north_aho->kink_cation_residue(ab_pose_aho));
	TS_ASSERT_EQUALS("TRP",WresAHO.name3());
	TS_ASSERT_EQUALS("ALA",AresAHO.name3());  // (not an anion in this PDB, but Ala is in the right position)
	TS_ASSERT_EQUALS("ARG",CresAHO.name3());


	// TODO: add tests for kink_anion_atoms and kink_cation_atoms
	return;
}


};

