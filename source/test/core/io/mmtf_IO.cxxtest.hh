// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/mmtf.cxxtest.hh
/// @brief  test suite for basic mmtf reading/writing
/// @author Danny Farrell

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <core/chemical/AA.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <src/core/io/mmtf/mmtf_reader.hh>
#include <src/core/io/mmtf/mmtf_writer.hh>
#include <src/core/io/mmcif/cif_reader.hh>
#include <core/io/StructFileReaderOptions.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>
#include <core/types.hh>
#include <utility/string_util.hh>
#include <utility/version.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/extra_pose_info_util.hh>

#include <core/simple_metrics/util.hh>
#include <core/simple_metrics/SimpleMetricData.hh>
#include <core/simple_metrics/StringMetric.hh>
#include <core/simple_metrics/RealMetric.hh>

#include <core/simple_metrics/CompositeStringMetric.hh>
#include <core/simple_metrics/CompositeRealMetric.hh>

#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/PerResidueStringMetric.hh>

#include <core/simple_metrics/test_classes.hh>
#include <core/simple_metrics/test_classes.fwd.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>

#include <mmtf.hpp>
#include <cifparse/CifFile.h>
#include <cifparse/CifParserBase.h>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>

#include <utility/version.hh>

static basic::Tracer TR("core.io.mmtf_IO.cxxtest");

class mmtf_IO : public CxxTest::TestSuite
{
private:
	core::pose::PoseOP pdb_pose;
	utility::vector1< std::string > mmtf_test_files;

	// Helper function three letter code to one letter code (doesnt fail, returns X)
	char tlc_2_olc( std::string const & name )
	{
		if ( core::chemical::is_aa_name_unknown( name ) ) return 'X';
		return core::chemical::oneletter_code_from_aa( core::chemical::aa_from_one_or_three( name ) );
	}

	std::string atoms_to_sequence( core::io::ChainAtoms const & all_atoms ) {
		std::string sequence("");
		for ( auto const & atm : all_atoms ) sequence += tlc_2_olc( atm.resName );
		return sequence;
	}


public:
	mmtf_IO() {}

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options("-ignore_zero_occupancy false -ignore_unrecognized_res -packing::pack_missing_sidechains false");

		pdb_pose = core::import_pose::pose_from_file( "core/io/1QYS.pdb", false , core::import_pose::PDB_file);
		if ( !mmtf_test_files.empty() ) return;
		// "173D", "1BNA", // missing group name
		// "1IGT", // has BMA hetatm
		mmtf_test_files.push_back("1MSH");
		mmtf_test_files.push_back("1SKM");
		mmtf_test_files.push_back("4CK4");
		mmtf_test_files.push_back("4P3R");
		mmtf_test_files.push_back("5EMG"); // "empty-mmtfVersion99999999.mmtf",
		mmtf_test_files.push_back("1AA6");
		mmtf_test_files.push_back("1CAG");
		mmtf_test_files.push_back("1L2Q");
		mmtf_test_files.push_back("1O2F");
		mmtf_test_files.push_back("3NJW");
		mmtf_test_files.push_back("4CUP");
		// "4V5A", // too large for pdb file comparison
		mmtf_test_files.push_back("5ESW"); // "empty-numChains1.mmtf",
		mmtf_test_files.push_back("1AUY");
		mmtf_test_files.push_back("1HTQ");
		mmtf_test_files.push_back("1LPV");
		// mmtf_test_files.push_back("1R9V"); // failed to set abase2 for acceptor atom, it has no nbrs! (pdb fails)
		// mmtf_test_files.push_back("3ZYB"); // mmcif mmtf difference
		mmtf_test_files.push_back("4OPJ");
		mmtf_test_files.push_back("4Y60");
		// "empty-all0.mmtf", "empty-numModels1.mmtf"
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_mmtf_noncanonicals() {
		std::map< std::string, std::string > noncanonical_pdb_to_mmtf({
			{"core/io/mmtf/binselector_input_S_0002.pdb", "core/io/mmtf/binselector_input_S_0002.io.mmtf"},
			{"core/io/mmtf/crosslinkermover_octahedral_s2_symm_S_0008.pdb", "core/io/mmtf/crosslinkermover_octahedral_s2_symm_S_0008.io.mmtf"},
			{"core/io/mmtf/ncbb_packer_palette_cycpep_pre_stage_3.pdb", "core/io/mmtf/ncbb_packer_palette_cycpep_pre_stage_3.io.mmtf"},
			{"core/io/mmtf/oligourea_predict_native.pdb", "core/io/mmtf/oligourea_predict_native.io.mmtf"},
			});
		for ( auto const & pair : noncanonical_pdb_to_mmtf ) {
			std::string const & test_file(pair.first);
			std::string const & reference_test_file(pair.second);
			TR << "Working on noncanonical test file: " << test_file << std::endl;

			// 1. load pdb file
			core::pose::PoseOP og_pose = pose_from_file(test_file, false,
				core::import_pose::FileType::PDB_file);
			og_pose->dump_pdb(test_file + ".io.pdb");

			// 2. create the sfr
			core::io::StructFileRepOptionsOP options = utility::pointer::make_shared< core::io::StructFileRepOptions >();
			core::io::mmtf::set_mmtf_default_options( *options );
			core::io::pose_to_sfr::PoseToStructFileRepConverter converter( *options );
			converter.init_from_pose( *og_pose );
			core::io::StructFileRepOP sfr = converter.sfr();

			// 3. write mmtf via sfr
			utility::io::ozstream file(reference_test_file.c_str(), std::ios::out | std::ios::binary);
			core::io::mmtf::dump_mmtf(file, sfr, *options);
			file.close();

			// 4. load sfr from written mmtf
			core::io::StructFileReaderOptions opts;
			core::io::StructFileRepOP sfr_in(core::io::mmtf::create_sfr_from_mmtf_filename(
				noncanonical_pdb_to_mmtf[test_file], opts));

			// 5. compare 2 and 4
			TS_ASSERT_EQUALS(sfr->chains().size(), sfr_in->chains().size());
			for ( core::Size i=0; i<sfr->chains().size(); ++i ) {
				for ( core::Size j=0; j<sfr->chains()[i].size(); ++j ) {
					core::io::AtomInformation const & ogai(sfr->chains()[i][j]);
					core::io::AtomInformation const & cbai(sfr_in->chains()[i][j]);
					TS_ASSERT_EQUALS(ogai.isHet, cbai.isHet);
					TS_ASSERT_EQUALS(ogai.serial, cbai.serial);
					TS_ASSERT_EQUALS(utility::strip(ogai.name), utility::strip(cbai.name));
					TS_ASSERT_EQUALS(ogai.altLoc, cbai.altLoc);
					TS_ASSERT_EQUALS(ogai.resName, cbai.resName);
					TS_ASSERT_EQUALS(ogai.chainID, cbai.chainID);
					TS_ASSERT_EQUALS(ogai.resSeq, cbai.resSeq);
					TS_ASSERT_EQUALS(ogai.iCode, cbai.iCode);
					TS_ASSERT_DELTA(ogai.x, cbai.x, 0.001); // extra large beacause might add atoms etc
					TS_ASSERT_DELTA(ogai.y, cbai.y, 0.001);
					TS_ASSERT_DELTA(ogai.z, cbai.z, 0.001);
					TS_ASSERT_EQUALS(ogai.occupancy, cbai.occupancy);
					TS_ASSERT_EQUALS(ogai.temperature, cbai.temperature);
					TS_ASSERT_EQUALS(ogai.segmentID, cbai.segmentID);
					TS_ASSERT_EQUALS(utility::strip(ogai.element), cbai.element);
					TS_ASSERT_EQUALS(ogai.formalcharge, cbai.formalcharge);
					TS_ASSERT_EQUALS(ogai.terCount, cbai.terCount);
				}
			}
		}
	}

	void test_mmtf_model_io() {
		// core_init_with_additional_options("-ignore_zero_occupancy false -ignore_unrecognized_res -packing::pack_missing_sidechains false -include_sugars -auto_detect_glycan_connections -alternate_3_letter_codes pdb_sugar -ignore_unrecognized_res -load_PDB_components false -write_glycan_pdb_codes");
		std::string const pdb1_fn("core/io/5FYL.pdb");
		std::string const mmtf1_ofn("core/io/5FYL_m_io.mmtf");
		std::string const pdb1_ofn("core/io/5FYL_m_io.pdb");
		core::pose::PoseOP pose1 = pose_from_file(pdb1_fn, false,
			core::import_pose::FileType::PDB_file);
		pose1->dump_mmtf(mmtf1_ofn);
		std::string const pdb2_fn("core/io/1QYS.pdb");
		std::string const mmtf2_ofn("core/io/1QYS_m_io.mmtf");
		std::string const pdb2_ofn("core/io/1QYS_m_io.pdb");
		core::pose::PoseOP pose2 = core::import_pose::pose_from_file(pdb2_fn, false,
			core::import_pose::PDB_file);
		pose2->dump_mmtf(mmtf2_ofn);

		core::io::StructFileReaderOptions opts;
		utility::vector1< core::io::StructFileRepOP > sfrs;
		core::io::StructFileRepOptionsOP options =  core::io::StructFileRepOptionsOP( new core::io::StructFileRepOptions );
		core::io::mmtf::set_mmtf_default_options( *options );

		// 1. Instantiate our SFRs from mmtfs so we base our info only on the capabilities
		// of the mmtf_reader.  We are testing the ability to recreate 2 sfrs, after loading
		// them from a multi-model mmtf.
		core::io::StructFileRepOP sfr1(core::io::mmtf::create_sfr_from_mmtf_filename( mmtf1_ofn, opts ));
		core::io::StructFileRepOP sfr2(core::io::mmtf::create_sfr_from_mmtf_filename( mmtf2_ofn, opts ));
		sfrs.push_back(sfr1);
		sfrs.push_back(sfr2);

		// 2. Generate our multi-model mmtf from 2 sfrs
		std::string const out_mmtf_fn("core/io/5fyl_1qys.mmtf");
		utility::io::ozstream file(out_mmtf_fn.c_str(), std::ios::out | std::ios::binary);
		core::io::mmtf::dump_mmtf(file, sfrs, *options);
		file.close();

		// 3. load back our multi-model sfrs
		utility::vector1< core::io::StructFileRepOP > sfrs_in(core::io::mmtf::create_sfrs_from_mmtf_filename(
			out_mmtf_fn, opts, utility::vector1< core::Size >({0, 1})));
		// 4. compare input and output sfrs
		TS_ASSERT_EQUALS(sfrs.size(), sfrs_in.size());
		for ( core::Size i=1; i<=sfrs.size(); ++i ) {
			TS_ASSERT_EQUALS(sfrs[i]->chains().size(), sfrs_in[i]->chains().size());
			for ( core::Size j=0; j<sfrs[i]->chains().size(); ++j ) {
				for ( core::Size k=0; k<sfrs[i]->chains()[j].size(); ++k ) {
					core::io::AtomInformation const & ogai(sfrs[i]->chains()[j][k]);
					core::io::AtomInformation const & cbai(sfrs_in[i]->chains()[j][k]);
					TS_ASSERT_EQUALS(ogai.isHet, cbai.isHet);
					TS_ASSERT_EQUALS(ogai.serial, cbai.serial);
					TS_ASSERT_EQUALS(ogai.name, cbai.name);
					TS_ASSERT_EQUALS(ogai.altLoc, cbai.altLoc);
					TS_ASSERT_EQUALS(ogai.resName, cbai.resName);
					TS_ASSERT_EQUALS(ogai.chainID, cbai.chainID);
					TS_ASSERT_EQUALS(ogai.resSeq, cbai.resSeq);
					TS_ASSERT_EQUALS(ogai.iCode, cbai.iCode);
					TS_ASSERT_DELTA(ogai.x, cbai.x, 0.0003);
					TS_ASSERT_DELTA(ogai.y, cbai.y, 0.0003);
					TS_ASSERT_DELTA(ogai.z, cbai.z, 0.0003);
					TS_ASSERT_EQUALS(ogai.occupancy, cbai.occupancy);
					TS_ASSERT_EQUALS(ogai.temperature, cbai.temperature);
					TS_ASSERT_EQUALS(ogai.segmentID, cbai.segmentID);
					TS_ASSERT_EQUALS(ogai.element, cbai.element);
					TS_ASSERT_EQUALS(ogai.formalcharge, cbai.formalcharge);
					TS_ASSERT_EQUALS(ogai.terCount, cbai.terCount);
				}
			}
		}
	}

	void test_add_or_read_extra_data() {
		std::string const pdb_fn("core/io/5FYL.pdb");
		core::pose::PoseOP og_pose = pose_from_file(pdb_fn, false,
			core::import_pose::FileType::PDB_file);

		{ /// FOR residue_type_base_names ONLY
			core::io::StructFileRepOptionsOP options(core::io::StructFileRepOptionsOP( new core::io::StructFileRepOptions ));
			core::io::mmtf::set_mmtf_default_options( *options );
			core::io::pose_to_sfr::PoseToStructFileRepConverter converter = core::io::pose_to_sfr::PoseToStructFileRepConverter( *options );
			converter.init_from_pose( *og_pose );
			core::io::StructFileRepOP sfr(converter.sfr());

			utility::vector1< core::io::StructFileRepOP > sfrs;
			sfrs.push_back(sfr);
			{ // Stage 1: check if it works in a vacuum
				::mmtf::StructureData sd;
				core::io::mmtf::add_extra_data(sd, sfrs, *options);

				TS_ASSERT_EQUALS(sd.modelProperties.count("rosetta::residue_type_base_names"), 1);
				std::vector< std::map< std::string, std::pair< std::string, std::string > > > returned;
				::mmtf::MapDecoder const ep_MD(sd.modelProperties);

				ep_MD.decode("rosetta::residue_type_base_names", true, returned);
				TS_ASSERT_EQUALS(returned.size(), 1);
				for ( auto const & k_v : returned.at(0) ) {
					TS_ASSERT_EQUALS(k_v.second, sfr->residue_type_base_names().at(k_v.first));
				}
			}
			{ // Stage 2: test in context of dump_mmtf
				std::string const mmtf1_fn("core/io/5FYL_io.mmtf");
				core::io::mmtf::dump_mmtf(mmtf1_fn, sfr, *options);

				{ // simple
					::mmtf::StructureData sd;
					::mmtf::decodeFromFile(sd, mmtf1_fn);
					TS_ASSERT_EQUALS(sd.modelProperties.count("rosetta::residue_type_base_names"), 1);
				}

				core::io::StructFileReaderOptions opts;
				core::io::StructFileRepOP mmtfsfr(core::io::mmtf::create_sfr_from_mmtf_filename( mmtf1_fn, opts ));
				TS_ASSERT_EQUALS(mmtfsfr->residue_type_base_names().empty(), false);

				for ( auto const & k_v : mmtfsfr->residue_type_base_names() ) {
					TS_ASSERT_EQUALS(k_v.second, sfr->residue_type_base_names().at(k_v.first));
				}
				for ( auto const & k_v : sfr->residue_type_base_names() ) {
					TS_ASSERT_EQUALS(k_v.second, mmtfsfr->residue_type_base_names().at(k_v.first));
				}
			}
		} // END residue_type_base_names ONLY
		TR << "Tested residue_type_base_names" << std::endl;

		{ /// FOR heterogen_names ONLY
			core::io::StructFileRepOptionsOP options(core::io::StructFileRepOptionsOP( new core::io::StructFileRepOptions ));
			core::io::mmtf::set_mmtf_default_options( *options );
			options->set_use_pdb_format_HETNAM_records(true);
			core::io::pose_to_sfr::PoseToStructFileRepConverter converter = core::io::pose_to_sfr::PoseToStructFileRepConverter( *options );
			converter.init_from_pose( *og_pose );
			core::io::StructFileRepOP sfr(converter.sfr());
			utility::vector1< core::io::StructFileRepOP > sfrs;
			sfrs.push_back(sfr);
			// Stage 3: do this for heterogen_names
			{  // remove namespace
				// Also stolen from pdb io tests
				// These example lines are taken directly from the given examples at
				// http://www.wwpdb.org/documentation/file-format-content/format33/sect4.html
				// with the addition of a single "Rosetta-format" line.
				string const sample_pdb_lines(
					"HETNAM     NAG N-ACETYL-D-GLUCOSAMINE                                          \n"
					"HETNAM     SAD BETA-METHYLENE SELENAZOLE-4-CARBOXAMIDE ADENINE                 \n"
					"HETNAM  2  SAD DINUCLEOTIDE                                                    \n"
					"HETNAM     UDP URIDINE-5'-DIPHOSPHATE                                          \n"
					"HETNAM     UNX UNKNOWN ATOM OR ION                                             \n"
					"HETNAM     UNL UNKNOWN LIGAND                                                  \n"
					"HETNAM     B3P 2-[3-(2-HYDROXY-1,1-DIHYDROXYMETHYL-ETHYLAMINO)-                \n"
					"HETNAM   2 B3P  PROPYLAMINO]-2-HYDROXYMETHYL-PROPANE-1,3-DIOL                  \n"
					"HETNAM     Krp X  13Z Kryptonite, which will kill Superman -- Bwahaha!         \n" );

				utility::vector1< core::io::pdb::Record > const records( core::io::pdb::create_records_from_pdb_file_contents( sample_pdb_lines ) );

				core::Size const n_records( records.size() );
				TS_ASSERT_EQUALS( n_records, 9 );

				for ( core::Size i(1); i<=n_records; ++i ) {
					core::io::pdb::store_heterogen_name_record_in_sfr( records.at(i), *sfr );
				}
			}

			{ // Stage 4: check if it works in a vacuum
				std::string const mmtf2_fn("core/io/5FYL_io2.mmtf");
				::mmtf::StructureData sd;
				core::io::mmtf::add_extra_data(sd, sfrs, *options);
				TS_ASSERT_EQUALS(sd.modelProperties.count("rosetta::heterogen_names"), 1);

				std::vector< std::map< std::string, std::string > > returned;
				::mmtf::MapDecoder const ep_MD(sd.modelProperties);
				ep_MD.decode("rosetta::heterogen_names", true, returned);
				TS_ASSERT_EQUALS(returned.size(), 1);

				for ( auto const & k_v : returned.at(0) ) {
					TS_ASSERT_EQUALS(k_v.second, sfr->heterogen_names().at(k_v.first));
				}
				for ( auto const & k_v : sfr->heterogen_names() ) {
					TS_ASSERT_EQUALS(k_v.second, returned.at(0).at(k_v.first));
				}
				TS_ASSERT_EQUALS(returned.at(0).at("Krp"), "X  13Z Kryptonite, which will kill Superman -- Bwahaha!");
			}
			{ // test in context of dumping pose
				std::string const mmtf3_fn("core/io/5FYL_io3.mmtf");
				core::io::mmtf::dump_mmtf(mmtf3_fn, sfr, *options);

				core::io::StructFileReaderOptions opts;
				core::io::StructFileRepOP mmtfsfr(core::io::mmtf::create_sfr_from_mmtf_filename( mmtf3_fn, opts ));

				for ( auto const & k_v : mmtfsfr->heterogen_names() ) {
					TS_ASSERT_EQUALS(k_v.second, sfr->heterogen_names().at(k_v.first));
				}
				for ( auto const & k_v : sfr->heterogen_names() ) {
					TS_ASSERT_EQUALS(k_v.second, mmtfsfr->heterogen_names().at(k_v.first));
				}
			}
		} // END heterogen_names ONLY
		{ // TEST for versioning info
			core::io::StructFileRepOptionsOP options(core::io::StructFileRepOptionsOP( new core::io::StructFileRepOptions ));
			core::io::mmtf::set_mmtf_default_options( *options );
			core::io::pose_to_sfr::PoseToStructFileRepConverter converter = core::io::pose_to_sfr::PoseToStructFileRepConverter( *options );
			converter.init_from_pose( *og_pose );
			core::io::StructFileRepOP sfr(converter.sfr());
			std::string const mmtf_fn("core/io/5FYL_versioning_io.mmtf");
			core::io::mmtf::dump_mmtf(mmtf_fn, sfr, *options);
			::mmtf::StructureData sd;
			::mmtf::decodeFromFile(sd, mmtf_fn);
			::mmtf::MapDecoder const extraProperties_MD(sd.extraProperties);
			std::string io_rosetta_version, io_rosetta_commit, io_rosetta_date, io_rosetta_package, io_rosetta_revision;
			extraProperties_MD.decode("rosetta::version", true, io_rosetta_version);
			TS_ASSERT_EQUALS(utility::Version::version(), io_rosetta_version);
			extraProperties_MD.decode("rosetta::commit", true, io_rosetta_commit);
			TS_ASSERT_EQUALS(utility::Version::commit(), io_rosetta_commit);
			extraProperties_MD.decode("rosetta::date_created", true, io_rosetta_date);
			TS_ASSERT_EQUALS(utility::Version::date(), io_rosetta_date);
			// I'm not exactly sure how to test this, since it will change for a non-devel package i think right?
			if ( utility::Version::package() != "devel" ) {
				extraProperties_MD.decode("rosetta::package", true, io_rosetta_package);
				TS_ASSERT_EQUALS(utility::Version::package(), io_rosetta_package);
			}
			if ( utility::Version::revision() != "None" ) {
				extraProperties_MD.decode("rosetta::revision", true, io_rosetta_revision);
				TS_ASSERT_EQUALS(utility::Version::revision(), io_rosetta_revision);
				sd.extraProperties["rosetta::revision"] = msgpack::object(utility::Version::revision(), sd.msgpack_zone);
			}
		}
	}

	void test_bond_info_io() {
		std::string const pdb_fn("core/io/mmtf/1PEF.pdb");
		core::pose::PoseOP og_pose = pose_from_file(pdb_fn, false,
			core::import_pose::FileType::PDB_file);

		core::io::StructFileRepOptionsOP options(core::io::StructFileRepOptionsOP( new core::io::StructFileRepOptions ));
		core::io::mmtf::set_mmtf_default_options( *options );
		core::io::pose_to_sfr::PoseToStructFileRepConverter converter = core::io::pose_to_sfr::PoseToStructFileRepConverter( *options );
		converter.init_from_pose( *og_pose );
		core::io::StructFileRepOP sfr(converter.sfr());

		std::string const mmtf1_fn("core/io/mmtf/1PEF_1.mmtf");
		og_pose->dump_mmtf(mmtf1_fn);
		core::io::StructFileReaderOptions opts;
		core::io::StructFileRepOP mmtfsfr(core::io::mmtf::create_sfr_from_mmtf_filename( mmtf1_fn, opts ));
		TS_ASSERT_EQUALS(mmtfsfr->chains().size(), sfr->chains().size());
		for ( core::Size i=0; i<sfr->chains().size(); ++i ) {
			TS_ASSERT_EQUALS(mmtfsfr->chains()[i].size(), sfr->chains()[i].size());
			for ( core::Size j=0; j<sfr->chains()[i].size(); ++j ) {
				core::io::AtomInformation const & mmtf_ai(mmtfsfr->chains()[i][j]);
				core::io::AtomInformation const & pdb_ai(sfr->chains()[i][j]);
				TS_ASSERT_EQUALS(mmtf_ai.connected_indices.size(), pdb_ai.connected_indices.size());
				TS_ASSERT_EQUALS(mmtf_ai.connected_orders.size(), pdb_ai.connected_orders.size());
				std::map<core::Size, core::Size> mmtf_ids, pdb_ids;
				for ( core::Size k=1; k<=mmtf_ai.connected_indices.size(); ++k ) {
					mmtf_ids[mmtf_ai.connected_indices[k]] = mmtf_ai.connected_orders[k];
					pdb_ids[pdb_ai.connected_indices[k]] = pdb_ai.connected_orders[k];
				}
				for ( auto const & bond : pdb_ids ) {
					TS_ASSERT_EQUALS(mmtf_ids.count(bond.first), 1);
					if ( 1 <= bond.second && bond.second <= 3 ) {
						// For bond orders of 1-3, these should be the same.
						TS_ASSERT_EQUALS(mmtf_ids.at(bond.first), bond.second);
					} else {
						// MMTF smashes everything into an unknown bond order.
						TS_ASSERT_EQUALS( mmtf_ids.at(bond.first), 0 );
					}
				}
			}
		}
	}

	void test_add_link_and_ss_information() {
		::mmtf::StructureData sd;
		core::io::StructFileRepOP sfr( utility::pointer::make_shared< core::io::StructFileRep >() );
		std::vector< core::io::AtomInformation > all_AIs;
		// insert setup of sd.bondAtomList
		sd.bondAtomList.push_back(0); sd.bondAtomList.push_back(1);
		sd.bondAtomList.push_back(0); sd.bondAtomList.push_back(2);
		sd.bondAtomList.push_back(3); sd.bondAtomList.push_back(0);
		sd.bondAtomList.push_back(2); sd.bondAtomList.push_back(1);

		{
			core::io::AtomInformation ai0, ai1, ai2, ai3;
			ai0.name = "aa";
			ai0.serial = 0;
			ai0.resSeq = 1;
			ai0.chainID = 'A';

			ai1.name = "bb";
			ai1.serial = 1;
			ai1.resSeq = 2;
			ai1.chainID = 'A';

			ai2.name = "cc";
			ai2.serial = 2;
			ai2.resSeq = 3;
			ai2.chainID = 'A';

			ai3.name = "dd";
			ai3.serial = 2;
			ai3.resSeq = 4;
			ai3.chainID = 'B';

			all_AIs.push_back(ai0);
			all_AIs.push_back(ai1);
			all_AIs.push_back(ai2);
			all_AIs.push_back(ai3);
		}
		std::vector< core::Size > ai_to_model({0, 0, 0, 0});

		core::io::mmtf::add_link_and_ss_information(sd, *sfr, all_AIs, ai_to_model, 0);

		std::map< std::string, utility::vector1< core::io::LinkInformation > > const &
			link_map(sfr->link_map());
		TS_ASSERT(link_map.at("   1 A").size() == 3);
		TS_ASSERT(link_map.at("   2 A").size() == 1);

		{ // test first
			core::io::LinkInformation expected_link_01, expected_link_02, expected_link_03;
			expected_link_01.name1 = utility::pad_atom_name("aa");
			expected_link_01.name2 = utility::pad_atom_name("bb");
			expected_link_01.resSeq1 = 1;
			expected_link_01.resSeq2 = 2;
			expected_link_01.chainID1 = 'A';
			expected_link_01.chainID2 = 'A';
			expected_link_01.iCode1 = ' ';
			expected_link_01.iCode2 = ' ';
			expected_link_01.length = 0;
			{
				std::stringstream strstr;
				strstr << std::setw( 4 ) << std::right << expected_link_01.resSeq1 << expected_link_01.iCode1 << expected_link_01.chainID1;
				expected_link_01.resID1 = strstr.str();
			}
			{
				std::stringstream strstr;
				strstr << std::setw( 4 ) << std::right << expected_link_01.resSeq2 << expected_link_01.iCode2 << expected_link_01.chainID2;
				expected_link_01.resID2 = strstr.str();
			}
			TR << ">" << expected_link_01.resID1 << "<" << std::endl;

			expected_link_02.name1 = utility::pad_atom_name("aa");
			expected_link_02.name2 = utility::pad_atom_name("cc");
			expected_link_02.resSeq1 = 1;
			expected_link_02.resSeq2 = 3;
			expected_link_02.chainID1 = 'A';
			expected_link_02.chainID2 = 'A';
			expected_link_02.iCode1 = ' ';
			expected_link_02.iCode2 = ' ';
			{
				std::stringstream strstr;
				strstr << std::setw( 4 ) << std::right << expected_link_02.resSeq1 << expected_link_02.iCode1 << expected_link_02.chainID1;
				expected_link_02.resID1 = strstr.str();
			}
			{
				std::stringstream strstr;
				strstr << std::setw( 4 ) << std::right << expected_link_02.resSeq2 << expected_link_02.iCode2 << expected_link_02.chainID2;
				expected_link_02.resID2 = strstr.str();
			}
			expected_link_02.length = 0;

			expected_link_03.name1 = utility::pad_atom_name("aa");
			expected_link_03.name2 = utility::pad_atom_name("dd");
			expected_link_03.resSeq1 = 1;
			expected_link_03.resSeq2 = 4;
			expected_link_03.chainID1 = 'A';
			expected_link_03.chainID2 = 'B';
			expected_link_03.iCode1 = ' ';
			expected_link_03.iCode2 = ' ';
			{
				std::stringstream strstr;
				strstr << std::setw( 4 ) << std::right << expected_link_03.resSeq1 << expected_link_03.iCode1 << expected_link_03.chainID1;
				expected_link_03.resID1 = strstr.str();
			}
			{
				std::stringstream strstr;
				strstr << std::setw( 4 ) << std::right << expected_link_03.resSeq2 << expected_link_03.iCode2 << expected_link_03.chainID2;
				expected_link_03.resID2 = strstr.str();
			}
			expected_link_03.length = 0;

			utility::vector1< core::io::LinkInformation > expected_vector;
			expected_vector.push_back(expected_link_01);
			expected_vector.push_back(expected_link_02);
			expected_vector.push_back(expected_link_03);
			TR << expected_vector << std::endl;
			TR << expected_link_01.resID1 << std::endl;
			TR << link_map.at(expected_link_01.resID1) << std::endl;

			TS_ASSERT( link_map.at(expected_link_01.resID1) == expected_vector);
		}

		{ // test second
			core::io::LinkInformation expected_link_01, expected_link_02;
			expected_link_01.name1 = utility::pad_atom_name("aa");
			expected_link_01.name2 = utility::pad_atom_name("bb");
			expected_link_01.resSeq1 = 1;
			expected_link_01.resSeq2 = 2;
			expected_link_01.chainID1 = 'A';
			expected_link_01.chainID2 = 'A';
			expected_link_01.iCode1 = ' ';
			expected_link_01.iCode2 = ' ';
			{
				std::stringstream strstr;
				strstr << std::setw( 4 ) << std::right << expected_link_01.resSeq1 << expected_link_01.iCode1 << expected_link_01.chainID1;
				expected_link_01.resID1 = strstr.str();
			}
			{
				std::stringstream strstr;
				strstr << std::setw( 4 ) << std::right << expected_link_01.resSeq2 << expected_link_01.iCode2 << expected_link_01.chainID2;
				expected_link_01.resID2 = strstr.str();
			}

			expected_link_02.name1 = utility::pad_atom_name("aa");
			expected_link_02.name2 = utility::pad_atom_name("cc");
			expected_link_02.resSeq1 = 1;
			expected_link_02.resSeq2 = 3;
			expected_link_02.chainID1 = 'A';
			expected_link_02.chainID2 = 'A';
			expected_link_02.iCode1 = ' ';
			expected_link_02.iCode2 = ' ';
			{
				std::stringstream strstr;
				strstr << std::setw( 4 ) << std::right << expected_link_02.resSeq1 << expected_link_02.iCode1 << expected_link_02.chainID1;
				expected_link_02.resID1 = strstr.str();
			}
			{
				std::stringstream strstr;
				strstr << std::setw( 4 ) << std::right << expected_link_02.resSeq2 << expected_link_02.iCode2 << expected_link_02.chainID2;
				expected_link_02.resID2 = strstr.str();
			}
			utility::vector1< core::io::LinkInformation > expected_vector;
			expected_vector.push_back(expected_link_01);
			expected_vector.push_back(expected_link_02);
		}
	}

	// Read a pdb file, write it to a mmtf file, read that mmtf file, and then write it again, and then compare the two
	// mmtf files to make sure they're the same.
	void test_mmtf_roundtrip() {
		core::Size initial_size = pdb_pose->size();

		std::string mmtf_file_1 = "first_out.mmtf";
		pdb_pose->dump_mmtf(mmtf_file_1);

		core::pose::PoseOP mmtf_pose_1 = pose_from_file(mmtf_file_1, false,
			core::import_pose::FileType::MMTF_file);
		TS_ASSERT_EQUALS( initial_size, mmtf_pose_1->size() );

		std::string mmtf_file_2 = "second_out.mmtf";
		mmtf_pose_1->dump_mmtf(mmtf_file_2);

		core::pose::PoseOP mmtf_pose_2 = pose_from_file(mmtf_file_2, false,
			core::import_pose::FileType::MMTF_file);
		TS_ASSERT_EQUALS( initial_size, mmtf_pose_2->size() );

		::mmtf::StructureData sd1;
		::mmtf::decodeFromFile(sd1, mmtf_file_1);
		::mmtf::StructureData sd2;
		::mmtf::decodeFromFile(sd2, mmtf_file_2);
		TS_ASSERT_EQUALS( true, sd1 == sd2 );
	}

	/** we cannot check sd1 == sd2 because rosetta will add atoms
	*  so instead we simply compare a cif import to an mmtf import
	*  for each cif file in the MMTF test suite and see if they have
	*  the same # of residues, and that sequence is the same.
	*/
	void test_mmtf_testsuite() {
		for ( auto const & test_file : mmtf_test_files ) {
			TR << "Working on mmtf test file: " << test_file << std::endl;
			std::string mmtf_file = "core/io/mmtf/" + test_file + ".mmtf";
			std::string cif_file = "core/io/mmtf/" + test_file + ".cif";
			::mmtf::StructureData sd1;
			::mmtf::decodeFromFile(sd1, mmtf_file);
			TR << "Loaded files here" << std::endl;
			core::io::StructFileReaderOptions opts;

			core::io::StructFileRepOP mmtfsfr(core::io::mmtf::create_sfr_from_mmtf_filename( mmtf_file, opts ));
			TR << "Loaded mmtf file" << std::endl;

			std::string contents_of_file;
			utility::io::izstream file( cif_file );
			utility::slurp( file, contents_of_file );
			std::string diagnostics;
			CifFileOP cifFile( new CifFile );
			CifParserOP cifParser( new CifParser( cifFile.get() ) );
			cifParser->ParseString( contents_of_file, diagnostics );
			core::io::StructFileRepOP cifsfr( core::io::mmcif::create_sfr_from_cif_file_op( cifFile, opts ) );
			TR << "Loaded cif file" << std::endl;

			TS_ASSERT_EQUALS( cifsfr->chains().size(), mmtfsfr->chains().size() );
			for ( core::Size i=0; i < cifsfr->chains().size(); ++i ) {
				TS_ASSERT_EQUALS( cifsfr->chains()[i].size(), mmtfsfr->chains()[i].size() );
				utility::vector1< std::string > mmtfnames, cifnames;
				for ( auto const & x : mmtfsfr->chains()[i] ) mmtfnames.push_back(x.resName);
				for ( auto const & x : cifsfr->chains()[i] ) cifnames.push_back(x.resName);
				TS_ASSERT_EQUALS(mmtfnames, cifnames);
			}
		}
	}

	void test_aiPose_from_sfr() {
		core::io::StructFileRepOptionsOP options =  core::io::StructFileRepOptionsOP( new core::io::StructFileRepOptions );
		core::io::mmtf::set_mmtf_default_options( *options );
		core::io::pose_to_sfr::PoseToStructFileRepConverter converter = core::io::pose_to_sfr::PoseToStructFileRepConverter( *options );
		converter.init_from_pose( *pdb_pose );
		core::io::StructFileRepOP sfr =  converter.sfr();

		core::Size chain_num = 0, group_num = 1;
		auto ai_pose( core::io::mmtf::aiPose_from_sfr(*sfr) );

		TS_ASSERT_EQUALS(true, ai_pose[0][0][0].resName == pdb_pose->residue(1).name3() );

		for ( core::Size i=2; i<=pdb_pose->size(); ++i ) {
			char prev_chain_id( pdb_pose->pdb_info()->chain(i-1) );
			char chain_id( pdb_pose->pdb_info()->chain(i) );
			if ( prev_chain_id != chain_id ) {
				++chain_num;
				group_num = 0;
			}
			auto ai = ai_pose[chain_num][group_num][0];
			TS_ASSERT_EQUALS(true, ai.resName == pdb_pose->residue(i).name3() );
			TS_ASSERT_EQUALS(true, ai_pose[chain_num][group_num].size() == pdb_pose->residue(i).natoms());
			++group_num;
		}
	}

	void test_mmtf_string_float_map_IO() {
		using namespace core::pose;
		///TODO: Make sure option is on here.

		///Test Comment, and String/Float map
		add_comment(*pdb_pose, "test1", "comment1");
		add_comment(*pdb_pose, "test2", "comment2");

		setPoseExtraScore(*pdb_pose, "sscore1", "value1");
		setPoseExtraScore(*pdb_pose, "sscore2", "value2");

		setPoseExtraScore(*pdb_pose, "rscore1", 1);
		setPoseExtraScore(*pdb_pose, "rscore2", 2);

		std::string mmtf_file_1 = "extra_out.mmtf";
		pdb_pose->dump_mmtf(mmtf_file_1);

		core::pose::PoseOP mmtf_pose_1 = pose_from_file(mmtf_file_1, false,
			core::import_pose::FileType::MMTF_file);

		std::string c1; std::string c2;
		std::string v1; std::string v2;
		core::Real r1; core::Real r2;

		//Retrieval
		getPoseExtraScore(*mmtf_pose_1, "sscore1", v1);
		getPoseExtraScore(*mmtf_pose_1, "sscore2", v2);
		TS_ASSERT(v1 == "value1");
		TS_ASSERT(v2 == "value2");

		getPoseExtraScore(*mmtf_pose_1, "rscore1", r1);
		getPoseExtraScore(*mmtf_pose_1, "rscore2", r2);

		TS_ASSERT_DELTA(r1, 1.0, .001);
		TS_ASSERT_DELTA(r2, 2.0, .001);

		std::map< std::string, std::string > comments = get_all_comments(*mmtf_pose_1);
		TS_ASSERT(comments.size() > 0);

		get_comment(*mmtf_pose_1, "test1", c1);
		get_comment(*mmtf_pose_1, "test2", c2);
		TS_ASSERT(c1 == "comment1");
		TS_ASSERT(c2 == "comment2");

	}
	void test_mmtf_simple_metric_IO() {
		using namespace core::simple_metrics;

		//Setup SimpleMetrics
		TestStringMetric test_string = TestStringMetric();
		test_string.apply(*pdb_pose);

		TestRealMetric test_real = TestRealMetric();
		test_real.apply(*pdb_pose);

		TestCompositeStringMetric test_composite_string = TestCompositeStringMetric();
		test_composite_string.apply(*pdb_pose);

		TestCompositeRealMetric test_composite_real = TestCompositeRealMetric();
		test_composite_real.apply(*pdb_pose);


		TestPerResidueRealMetric test_per_residue_real = TestPerResidueRealMetric();
		test_per_residue_real.apply(*pdb_pose);


		TestPerResidueStringMetric test_per_residue_string = TestPerResidueStringMetric();
		test_per_residue_string.apply(*pdb_pose);

		//Ouptut mmTF and Read back in
		std::string mmtf_file_1 = "extra_out2.mmtf";
		pdb_pose->dump_mmtf(mmtf_file_1);

		core::pose::PoseOP new_pose = pose_from_file(mmtf_file_1, false,
			core::import_pose::FileType::MMTF_file);


		//Assert all values are present.
		TS_ASSERT( has_sm_data(*new_pose));

		//get_sm_data(*new_pose)->show();

		///StringMetric
		std::string svalue;
		bool present = get_sm_data(*new_pose)->get_value( "SomeString", svalue );
		TS_ASSERT( present );
		TS_ASSERT( svalue == "TESTING");

		///RealMetric
		core::Real rvalue;
		present = get_sm_data(*new_pose)->get_value( "SomeReal", rvalue );
		TS_ASSERT( present );
		TS_ASSERT_DELTA( rvalue, 1.0, .001 );

		///CompositeStringMetric
		std::map< std::string, std::string > csvalues;
		present = get_sm_data(*new_pose)->get_value("SomeCompositeString", csvalues );
		TS_ASSERT( present );
		if ( present ) {
			TS_ASSERT( csvalues["s_data1"] == "value1");
			TS_ASSERT( csvalues["s_data2"] == "value2");
		}

		///CompositeRealMetric
		std::map< std::string, core::Real > crvalues;
		present = get_sm_data( *new_pose )->get_value( "SomeCompositeReal", crvalues );
		TS_ASSERT( present );
		if ( present ) {
			TS_ASSERT_DELTA( crvalues["r_data1"], 1.0, .001);
			TS_ASSERT_DELTA( crvalues["r_data2"], 2.0, .001);
		}

		///PerResidueStringMetric
		std::map< core::Size, std::string > res_str_values;
		present = get_sm_data(*new_pose)->get_value( "SomePerResidueString", res_str_values);
		TS_ASSERT(present);
		if ( present ) {
			TS_ASSERT( res_str_values.at(1) == "value1");
			TS_ASSERT( res_str_values.at(2) == "value2");
		}

		///PerResidueRealMetric
		std::map< core::Size, core::Real > res_r_values;
		present = get_sm_data(*new_pose)->get_value( "SomePerResidueReal", res_r_values);
		TS_ASSERT(present);
		if ( present ) {
			TS_ASSERT_DELTA( res_r_values.at(1), 1.0, .001);
			TS_ASSERT_DELTA( res_r_values.at(2), 2.0, .001);
		}

		//Assert per-res output makes round-trip
		std::map< std::string, std::map< std::string, core::Real >> const & real_output = get_sm_data(*new_pose)->get_per_residue_real_metric_output();
		TS_ASSERT(real_output.size());
		TS_ASSERT(real_output.count("SomePerResidueReal"));


		std::map< std::string, std::map< std::string, std::string >> const & string_output = get_sm_data(*new_pose)->get_per_residue_string_metric_output();
		TS_ASSERT(string_output.size());
		TS_ASSERT(string_output.count("SomePerResidueString"));

	}
};
