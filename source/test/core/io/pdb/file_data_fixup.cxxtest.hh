// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/file_data_fixup.cxxtest.hh
/// @brief  test suite for functions associated with core/io/pdb/file_data_fixup.hh
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/io/pdb/file_data_fixup.hh>

// Program headers
#include <core/io/pdb/pdb_dynamic_reader.hh>
#include <core/import_pose/import_pose_options.hh>
#include <core/chemical/residue_io.cc>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers

// ObjexxFCL headers

// C++ headers
#include <string>

static basic::Tracer TR("core.io.pdb.file_data_fixup.cxxtest");

using namespace core;
using namespace core::io::pdb;

core::io::pdb::ResidueInformation create_ResidueInfo( std::string const & instring ) {
	using namespace core::io::pdb;
	import_pose::ImportPoseOptions options; //default options
	FileData fd = PDB_DReader::createFileData( instring, options );
	utility::vector1< ResidueInformation > rinfos;
	fd.create_working_data( rinfos, options );
	assert( rinfos.size() );
	if( rinfos.size() > 1 ) {
		utility_exit_with_message("Too many residues.");
	}
	return rinfos[1];
}

class geometric_rename_Tests : public CxxTest::TestSuite
{

	core::chemical::ResidueTypeCOP rsd_;
	core::chemical::ResidueTypeCOP gly_;
	ResidueInformation main_rinfo_;

public:
	// Shared initialization goes here.
	void setUp() {
		core_init();

		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		AtomTypeSetCAP atom_types = cm->atom_type_set(tag);
		ElementSetCAP element_types = cm->element_set("default");
		MMAtomTypeSetCAP mm_atom_types = cm->mm_atom_type_set(tag);
		orbitals::OrbitalTypeSetCAP orbital_types = cm->orbital_type_set(tag);

		ResidueTypeCOPs const & glycines( cm->residue_type_set(tag)->aa_map( aa_gly ) );
		assert( glycines.size() > 1 );
		gly_ = glycines[1];

		ResidueTypeSetOP rsd_types( new ResidueTypeSet );
		
		std::string filename("core/chemical/params/1aq1.mol2.params");
		rsd_ = read_topology_file(filename, atom_types, element_types, mm_atom_types, orbital_types, ResidueTypeSetCAP(rsd_types));

		main_rinfo_ = create_ResidueInfo(
"HETATM    1  O2  LG1 X 299      -0.058   0.605   0.915  1.00  0.00           O  \n"
"HETATM    2  C24 LG1 X 299       1.034  -0.296   1.023  1.00  0.00           C  \n"
"HETATM    3  C21 LG1 X 299       2.397   0.426   0.885  1.00  0.00           C  \n"
"HETATM    4  C22 LG1 X 299       2.168   1.910   0.585  1.00  0.00           C  \n"
"HETATM    5  C23 LG1 X 299       1.238   2.049  -0.655  1.00  0.00           C  \n"
"HETATM    6  C25 LG1 X 299      -0.153   1.537  -0.148  1.00  0.00           C  \n"
"HETATM    7  C26 LG1 X 299      -0.872   2.730   0.452  1.00  0.00           C  \n"
"HETATM    8  N3  LG1 X 299      -0.965   0.941  -1.242  1.00  0.00           N  \n"
"HETATM    9  C16 LG1 X 299      -1.883   1.537  -2.121  1.00  0.00           C  \n"
"HETATM   10  C8  LG1 X 299      -2.338   2.811  -2.270  1.00  0.00           C  \n"
"HETATM   11  C4  LG1 X 299      -3.258   3.103  -3.234  1.00  0.00           C  \n"
"HETATM   12  C2  LG1 X 299      -3.747   2.094  -4.117  1.00  0.00           C  \n"
"HETATM   13  C6  LG1 X 299      -3.266   0.816  -3.963  1.00  0.00           C  \n"
"HETATM   14  C10 LG1 X 299      -2.358   0.500  -2.981  1.00  0.00           C  \n"
"HETATM   15  C12 LG1 X 299      -1.717  -0.708  -2.627  1.00  0.00           C  \n"
"HETATM   16  C13 LG1 X 299      -1.767  -2.046  -3.084  1.00  0.00           C  \n"
"HETATM   17  C14 LG1 X 299      -1.044  -2.996  -2.549  1.00  0.00           C  \n"
"HETATM   18  C11 LG1 X 299      -0.120  -2.762  -1.426  1.00  0.00           C  \n"
"HETATM   19  C9  LG1 X 299       0.797  -3.520  -0.676  1.00  0.00           C  \n"
"HETATM   20  C5  LG1 X 299       1.182  -4.841  -0.625  1.00  0.00           C  \n"
"HETATM   21  C1  LG1 X 299       2.110  -5.264   0.279  1.00  0.00           C  \n"
"HETATM   22  C3  LG1 X 299       2.676  -4.369   1.218  1.00  0.00           C  \n"
"HETATM   23  C7  LG1 X 299       2.311  -3.043   1.224  1.00  0.00           C  \n"
"HETATM   24  C15 LG1 X 299       1.385  -2.627   0.266  1.00  0.00           C  \n"
"HETATM   25  N2  LG1 X 299       0.877  -1.390   0.070  1.00  0.00           N  \n"
"HETATM   26  C17 LG1 X 299      -0.099  -1.414  -0.954  1.00  0.00           C  \n"
"HETATM   27  C18 LG1 X 299      -0.838  -0.448  -1.507  1.00  0.00           C  \n"
"HETATM   28  C20 LG1 X 299      -1.268  -4.294  -3.262  1.00  0.00           C  \n"
"HETATM   29  N1  LG1 X 299      -2.169  -3.976  -4.202  1.00  0.00           N  \n"
"HETATM   30  C19 LG1 X 299      -2.600  -2.594  -4.258  1.00  0.00           C  \n"
"HETATM   31  O1  LG1 X 299      -0.769  -5.381  -3.040  1.00  0.00           O  \n"
"HETATM   32  O3  LG1 X 299       1.755   1.144  -1.649  1.00  0.00           O  \n"
"HETATM   33  C28 LG1 X 299       1.269   1.452  -2.941  1.00  0.00           C  \n"
"HETATM   34  N4  LG1 X 299       3.559   2.588   0.319  1.00  0.00           N  \n"
"HETATM   35  C27 LG1 X 299       4.002   3.352   1.518  1.00  0.00           C  \n"
"HETATM   36  H14 LG1 X 299       1.030  -0.730   2.034  1.00  0.00           H  \n"
"HETATM   37  H10 LG1 X 299       2.969  -0.029   0.062  1.00  0.00           H  \n"
"HETATM   38  H11 LG1 X 299       2.958   0.327   1.826  1.00  0.00           H  \n"
"HETATM   39  H12 LG1 X 299       1.684   2.407   1.438  1.00  0.00           H  \n"
"HETATM   40  H13 LG1 X 299       1.170   3.064  -1.074  1.00  0.00           H  \n"
"HETATM   41  H15 LG1 X 299      -0.350   3.656   0.166  1.00  0.00           H  \n"
"HETATM   42  H16 LG1 X 299      -1.906   2.762   0.077  1.00  0.00           H  \n"
"HETATM   43  H17 LG1 X 299      -0.881   2.639   1.548  1.00  0.00           H  \n"
"HETATM   44  H8  LG1 X 299      -1.970   3.593  -1.624  1.00  0.00           H  \n"
"HETATM   45  H4  LG1 X 299      -3.622   4.115  -3.330  1.00  0.00           H  \n"
"HETATM   46  H2  LG1 X 299      -4.471   2.328  -4.883  1.00  0.00           H  \n"
"HETATM   47  H6  LG1 X 299      -3.611   0.040  -4.630  1.00  0.00           H  \n"
"HETATM   48  H5  LG1 X 299       0.743  -5.551  -1.311  1.00  0.00           H  \n"
"HETATM   49  H1  LG1 X 299       2.418  -6.299   0.279  1.00  0.00           H  \n"
"HETATM   50  H3  LG1 X 299       3.400  -4.730   1.933  1.00  0.00           H  \n"
"HETATM   51  H7  LG1 X 299       2.725  -2.350   1.943  1.00  0.00           H  \n"
"HETATM   52  H9  LG1 X 299      -3.294  -2.098  -4.921  1.00  0.00           H  \n"
"HETATM   53  H21 LG1 X 299       0.693   0.600  -3.330  1.00  0.00           H  \n"
"HETATM   54  H22 LG1 X 299       0.621   2.339  -2.887  1.00  0.00           H  \n"
"HETATM   55  H23 LG1 X 299       2.117   1.658  -3.612  1.00  0.00           H  \n"
"HETATM   56  H24 LG1 X 299       3.476   3.214  -0.457  1.00  0.00           H  \n"
"HETATM   57  H25 LG1 X 299       4.234   1.878   0.118  1.00  0.00           H  \n"
"HETATM   58  H18 LG1 X 299       4.484   4.288   1.200  1.00  0.00           H  \n"
"HETATM   59  H19 LG1 X 299       3.130   3.584   2.147  1.00  0.00           H  \n"
"HETATM   60  H20 LG1 X 299       4.718   2.747   2.093  1.00  0.00           H  \n"
);
	}

	/// @brief Test the scoring scheme
	void test_scoring() {
		core::Size natoms = rsd_->natoms();
		core::Size natoms2 = natoms*natoms;
		core::Real score_delta = 0.1 * 1.0/ natoms2;

		NameBimap map;
		map.insert( NameBimap::value_type(" H1 ", " H3 ") );
		map.insert( NameBimap::value_type(" C1 ", " C3 ") );
		// Two atoms present, no name matching, no chirality
		TS_ASSERT_DELTA( score_mapping( map, main_rinfo_, *rsd_ ), 2.0/natoms, score_delta );
		map.insert( NameBimap::value_type(" C14", " C14") );
		map.insert( NameBimap::value_type(" C12", " C12") );
		// Four atoms present, two name matches, no chirality
		TS_ASSERT_DELTA( score_mapping( map, main_rinfo_, *rsd_ ), 4.0/natoms + 2.0/natoms2, score_delta );

		map.clear();
		map.insert( NameBimap::value_type(" C21", " C21") );
		map.insert( NameBimap::value_type(" C22", " C22") );
		map.insert( NameBimap::value_type(" C24", " C24") );
		map.insert( NameBimap::value_type(" H10", " H10") );
		// Four atoms present, four name matches, one three-member chirality
		TS_ASSERT_DELTA( score_mapping( map, main_rinfo_, *rsd_ ), 1.0 + 4.0/natoms + 4.0/natoms2, score_delta );

		map.clear();
		map.insert( NameBimap::value_type(" C21", " C21") );
		map.insert( NameBimap::value_type(" C22", " C22") );
		map.insert( NameBimap::value_type(" C24", " C24") );
		map.insert( NameBimap::value_type(" H10", " H11") );
		// Four atoms present, three name matches, mismatched three-member chirality
		TS_ASSERT_DELTA( score_mapping( map, main_rinfo_, *rsd_ ), 0.0 + 4.0/natoms + 3.0/natoms2, score_delta );
		map.insert( NameBimap::value_type(" H11", " H10") );
		// Five atoms present, three name matches, mismatched four-member chirality
		TS_ASSERT_DELTA( score_mapping( map, main_rinfo_, *rsd_ ), 0.0 + 5.0/natoms + 3.0/natoms2, score_delta );

		map.clear();
		map.insert( NameBimap::value_type(" C21", " C21") );
		map.insert( NameBimap::value_type(" C22", " C22") );
		map.insert( NameBimap::value_type(" C24", " C24") );
		map.insert( NameBimap::value_type(" H10", " H10") );
		map.insert( NameBimap::value_type(" H11", " H11") );
		// Five atoms present, five name matches, one four-member chirality
		TS_ASSERT_DELTA( score_mapping( map, main_rinfo_, *rsd_ ), 1.0 + 5.0/natoms + 5.0/natoms2, score_delta );

		map.clear();
		map.insert( NameBimap::value_type(" C26", " C27") );
		map.insert( NameBimap::value_type(" C25", " N4 ") );
		map.insert( NameBimap::value_type(" H15", " H18") );
		map.insert( NameBimap::value_type(" H16", " H19") );
		map.insert( NameBimap::value_type(" H17", " H20") );
		// Five atoms present, no name matches, one four-member chirality
		TS_ASSERT_DELTA( score_mapping( map, main_rinfo_, *rsd_ ), 1.0 + 5.0/natoms + 0.0/natoms2, score_delta );

		map.clear();
		map.insert( NameBimap::value_type(" C26", " C27") );
		map.insert( NameBimap::value_type(" C25", " N4 ") );
		map.insert( NameBimap::value_type(" H15", " H18") );
		map.insert( NameBimap::value_type(" H16", " H20") );
		map.insert( NameBimap::value_type(" H17", " H19") );
		// Five atoms present, no name matches, mismatched four-member chirality
		TS_ASSERT_DELTA( score_mapping( map, main_rinfo_, *rsd_ ), 0.0 + 5.0/natoms + 0.0/natoms2, score_delta );

		map.clear();
		map.insert( NameBimap::value_type(" C11", " C12") );
		map.insert( NameBimap::value_type(" C17", " C18") );
		map.insert( NameBimap::value_type(" C14", " C10") );
		map.insert( NameBimap::value_type(" C9 ", " C13") );
		// Four atoms present, no name matches, matched planar chirality
		TS_ASSERT_DELTA( score_mapping( map, main_rinfo_, *rsd_ ), 1.0 + 4.0/natoms + 0.0/natoms2, score_delta );

		map.clear();
		map.insert( NameBimap::value_type(" C11", " C12") );
		map.insert( NameBimap::value_type(" C17", " C18") );
		map.insert( NameBimap::value_type(" C14", " C13") );
		map.insert( NameBimap::value_type(" C9 ", " C10") );
		// Four atoms present, no name matches, matched planar chirality
		TS_ASSERT_DELTA( score_mapping( map, main_rinfo_, *rsd_ ), 1.0 + 4.0/natoms + 0.0/natoms2, score_delta );

		map.clear();

		map.insert( NameBimap::value_type(" C5 ", " C5 ") );
		map.insert( NameBimap::value_type(" C9 ", " C9 ") );
		map.insert( NameBimap::value_type(" C12", " C15") ); // C9-C15 has a bond, but C9-C12 does not.
		map.insert( NameBimap::value_type(" C10", " C7 ") );
		// Four atoms present, two name matches, one inappropriate bond.
		TS_ASSERT_DELTA( score_mapping( map, main_rinfo_, *rsd_ ), -2.0 + 4.0/natoms + 2.0/natoms2, score_delta );
	}

	/// @brief This is with the names of the atoms matched exactly.
	void test_identity_name_map() {

		TR << "Test identity" << std::endl;
		NameBimap name_map;
		remap_names_on_geometry(name_map, main_rinfo_, *rsd_);
		TS_ASSERT_EQUALS( name_map.left.size(), main_rinfo_.atoms.size());
		TS_ASSERT(name_map.left.count(" O2 "));
		TS_ASSERT(name_map.left.count(" N4 "));
		TS_ASSERT(name_map.left.count(" C27"));
		TS_ASSERT(name_map.left.count(" H16"));
		TS_ASSERT_EQUALS(name_map.left.find(" O2 ")->second," O2 ");
		TS_ASSERT_EQUALS(name_map.left.find(" N4 ")->second," N4 ");
		TS_ASSERT_EQUALS(name_map.left.find(" C27")->second," C27");
		TS_ASSERT_EQUALS(name_map.left.find(" H16")->second," H16");
	}

	/// @brief This is with all atoms present, but with name mismatches on all atoms.
	void test_all_mismatch_name_map() {

		TR << "Test all mismatch" << std::endl;
		ResidueInformation rinfo = create_ResidueInfo( // Name mismatches first digit gets incremented by two
"HETATM    1  O4  LG1 X 299      -0.058   0.605   0.915  1.00  0.00           O  \n"
"HETATM    2  C44 LG1 X 299       1.034  -0.296   1.023  1.00  0.00           C  \n"
"HETATM    3  C41 LG1 X 299       2.397   0.426   0.885  1.00  0.00           C  \n"
"HETATM    4  C42 LG1 X 299       2.168   1.910   0.585  1.00  0.00           C  \n"
"HETATM    5  C43 LG1 X 299       1.238   2.049  -0.655  1.00  0.00           C  \n"
"HETATM    6  C45 LG1 X 299      -0.153   1.537  -0.148  1.00  0.00           C  \n"
"HETATM    7  C46 LG1 X 299      -0.872   2.730   0.452  1.00  0.00           C  \n"
"HETATM    8  N5  LG1 X 299      -0.965   0.941  -1.242  1.00  0.00           N  \n"
"HETATM    9  C36 LG1 X 299      -1.883   1.537  -2.121  1.00  0.00           C  \n"
"HETATM   10  C1  LG1 X 299      -2.338   2.811  -2.270  1.00  0.00           C  \n"
"HETATM   11  C6  LG1 X 299      -3.258   3.103  -3.234  1.00  0.00           C  \n"
"HETATM   12  C4  LG1 X 299      -3.747   2.094  -4.117  1.00  0.00           C  \n"
"HETATM   13  C8  LG1 X 299      -3.266   0.816  -3.963  1.00  0.00           C  \n"
"HETATM   14  C30 LG1 X 299      -2.358   0.500  -2.981  1.00  0.00           C  \n"
"HETATM   15  C32 LG1 X 299      -1.717  -0.708  -2.627  1.00  0.00           C  \n"
"HETATM   16  C33 LG1 X 299      -1.767  -2.046  -3.084  1.00  0.00           C  \n"
"HETATM   17  C34 LG1 X 299      -1.044  -2.996  -2.549  1.00  0.00           C  \n"
"HETATM   18  C31 LG1 X 299      -0.120  -2.762  -1.426  1.00  0.00           C  \n"
"HETATM   19  C2  LG1 X 299       0.797  -3.520  -0.676  1.00  0.00           C  \n"
"HETATM   20  C7  LG1 X 299       1.182  -4.841  -0.625  1.00  0.00           C  \n"
"HETATM   21  C3  LG1 X 299       2.110  -5.264   0.279  1.00  0.00           C  \n"
"HETATM   22  C5  LG1 X 299       2.676  -4.369   1.218  1.00  0.00           C  \n"
"HETATM   23  C9  LG1 X 299       2.311  -3.043   1.224  1.00  0.00           C  \n"
"HETATM   24  C35 LG1 X 299       1.385  -2.627   0.266  1.00  0.00           C  \n"
"HETATM   25  N4  LG1 X 299       0.877  -1.390   0.070  1.00  0.00           N  \n"
"HETATM   26  C37 LG1 X 299      -0.099  -1.414  -0.954  1.00  0.00           C  \n"
"HETATM   27  C38 LG1 X 299      -0.838  -0.448  -1.507  1.00  0.00           C  \n"
"HETATM   28  C40 LG1 X 299      -1.268  -4.294  -3.262  1.00  0.00           C  \n"
"HETATM   29  N3  LG1 X 299      -2.169  -3.976  -4.202  1.00  0.00           N  \n"
"HETATM   30  C39 LG1 X 299      -2.600  -2.594  -4.258  1.00  0.00           C  \n"
"HETATM   31  O3  LG1 X 299      -0.769  -5.381  -3.040  1.00  0.00           O  \n"
"HETATM   32  O5  LG1 X 299       1.755   1.144  -1.649  1.00  0.00           O  \n"
"HETATM   33  C48 LG1 X 299       1.269   1.452  -2.941  1.00  0.00           C  \n"
"HETATM   34  N6  LG1 X 299       3.559   2.588   0.319  1.00  0.00           N  \n"
"HETATM   35  C47 LG1 X 299       4.002   3.352   1.518  1.00  0.00           C  \n"
"HETATM   36  H34 LG1 X 299       1.030  -0.730   2.034  1.00  0.00           H  \n"
"HETATM   37  H30 LG1 X 299       2.969  -0.029   0.062  1.00  0.00           H  \n"
"HETATM   38  H31 LG1 X 299       2.958   0.327   1.826  1.00  0.00           H  \n"
"HETATM   39  H32 LG1 X 299       1.684   2.407   1.438  1.00  0.00           H  \n"
"HETATM   40  H33 LG1 X 299       1.170   3.064  -1.074  1.00  0.00           H  \n"
"HETATM   41  H35 LG1 X 299      -0.350   3.656   0.166  1.00  0.00           H  \n"
"HETATM   42  H36 LG1 X 299      -1.906   2.762   0.077  1.00  0.00           H  \n"
"HETATM   43  H37 LG1 X 299      -0.881   2.639   1.548  1.00  0.00           H  \n"
"HETATM   44  H1  LG1 X 299      -1.970   3.593  -1.624  1.00  0.00           H  \n"
"HETATM   45  H6  LG1 X 299      -3.622   4.115  -3.330  1.00  0.00           H  \n"
"HETATM   46  H4  LG1 X 299      -4.471   2.328  -4.883  1.00  0.00           H  \n"
"HETATM   47  H8  LG1 X 299      -3.611   0.040  -4.630  1.00  0.00           H  \n"
"HETATM   48  H7  LG1 X 299       0.743  -5.551  -1.311  1.00  0.00           H  \n"
"HETATM   49  H3  LG1 X 299       2.418  -6.299   0.279  1.00  0.00           H  \n"
"HETATM   50  H5  LG1 X 299       3.400  -4.730   1.933  1.00  0.00           H  \n"
"HETATM   51  H9  LG1 X 299       2.725  -2.350   1.943  1.00  0.00           H  \n"
"HETATM   52  H2  LG1 X 299      -3.294  -2.098  -4.921  1.00  0.00           H  \n"
"HETATM   53  H41 LG1 X 299       0.693   0.600  -3.330  1.00  0.00           H  \n"
"HETATM   54  H42 LG1 X 299       0.621   2.339  -2.887  1.00  0.00           H  \n"
"HETATM   55  H43 LG1 X 299       2.117   1.658  -3.612  1.00  0.00           H  \n"
"HETATM   56  H44 LG1 X 299       3.476   3.214  -0.457  1.00  0.00           H  \n"
"HETATM   57  H45 LG1 X 299       4.234   1.878   0.118  1.00  0.00           H  \n"
"HETATM   58  H38 LG1 X 299       4.484   4.288   1.200  1.00  0.00           H  \n"
"HETATM   59  H39 LG1 X 299       3.130   3.584   2.147  1.00  0.00           H  \n"
"HETATM   60  H40 LG1 X 299       4.718   2.747   2.093  1.00  0.00           H  \n"
);
		NameBimap name_map;
		remap_names_on_geometry(name_map, rinfo, *rsd_);
		TS_ASSERT_EQUALS( name_map.left.size(), rinfo.atoms.size());
		TS_ASSERT(name_map.left.count(" O4 "));
		TS_ASSERT(name_map.left.count(" N6 "));
		TS_ASSERT(name_map.left.count(" C47"));
		TS_ASSERT(name_map.left.count(" H36"));
		TS_ASSERT_EQUALS(name_map.left.find(" O4 ")->second," O2 ");
		TS_ASSERT_EQUALS(name_map.left.find(" N6 ")->second," N4 ");
		TS_ASSERT_EQUALS(name_map.left.find(" C47")->second," C27");
		TS_ASSERT_EQUALS(name_map.left.find(" H36")->second," H16");
	}

	/// @brief Can we handle more atoms than we have in the residue?
	void test_too_many_atoms_name_map() {

		TR << "Test too many atoms" << std::endl;
		ResidueInformation rinfo = create_ResidueInfo( // Serine redone as glycine
"ATOM     19  N   GLY A 431      -3.764   5.641  23.700  1.00  0.00           N \n"
"ATOM     20  CA  GLY A 431      -2.495   6.302  23.526  1.00  0.00           C \n"
"ATOM     21  C   GLY A 431      -2.193   6.626  22.078  1.00  0.00           C \n"
"ATOM     22  O   GLY A 431      -2.749   6.028  21.148  1.00  0.00           O \n"
"ATOM     23  CB  GLY A 431      -1.389   5.433  24.093  1.00  0.00           C \n"
"ATOM     24  OG  GLY A 431      -1.544   5.260  25.474  1.00  0.00           O \n"
"ATOM     25  H   GLY A 431      -3.780   4.679  24.006  1.00  0.00           H \n"
"ATOM     26  HA  GLY A 431      -2.516   7.243  24.077  1.00  0.00           H \n"
"ATOM     27 1HB  GLY A 431      -1.400   4.462  23.598  1.00  0.00           H \n"
"ATOM     28 2HB  GLY A 431      -0.424   5.894  23.887  1.00  0.00           H \n"
"ATOM     29  HG  GLY A 431      -2.334   5.751  25.715  1.00  0.00           H \n"
);
		NameBimap name_map;
		remap_names_on_geometry(name_map, rinfo, *gly_);
		TS_ASSERT_EQUALS( name_map.left.size(), gly_->natoms() - 1);
		TS_ASSERT_EQUALS(name_map.left.find(" N  ")->second," N  ");
		TS_ASSERT_EQUALS(name_map.left.find(" CA ")->second," CA ");
		TS_ASSERT_EQUALS(name_map.left.find(" C  ")->second," C  ");
		TS_ASSERT_EQUALS(name_map.left.find(" HA ")->second,"1HA ");
	}

	/// @brief This tests if all hydrogens are missing.
	void test_no_hydro_name_map() {

		TR << "Test missing hydro" << std::endl;
		ResidueInformation rinfo = create_ResidueInfo( // Name mismatches first digit gets incremented by two
"HETATM    1  O4  LG1 X 299      -0.058   0.605   0.915  1.00  0.00           O  \n"
"HETATM    2  C44 LG1 X 299       1.034  -0.296   1.023  1.00  0.00           C  \n"
"HETATM    3  C41 LG1 X 299       2.397   0.426   0.885  1.00  0.00           C  \n"
"HETATM    4  C42 LG1 X 299       2.168   1.910   0.585  1.00  0.00           C  \n"
"HETATM    5  C43 LG1 X 299       1.238   2.049  -0.655  1.00  0.00           C  \n"
"HETATM    6  C45 LG1 X 299      -0.153   1.537  -0.148  1.00  0.00           C  \n"
"HETATM    7  C46 LG1 X 299      -0.872   2.730   0.452  1.00  0.00           C  \n"
"HETATM    8  N5  LG1 X 299      -0.965   0.941  -1.242  1.00  0.00           N  \n"
"HETATM    9  C36 LG1 X 299      -1.883   1.537  -2.121  1.00  0.00           C  \n"
"HETATM   10  C1  LG1 X 299      -2.338   2.811  -2.270  1.00  0.00           C  \n"
"HETATM   11  C6  LG1 X 299      -3.258   3.103  -3.234  1.00  0.00           C  \n"
"HETATM   12  C4  LG1 X 299      -3.747   2.094  -4.117  1.00  0.00           C  \n"
"HETATM   13  C8  LG1 X 299      -3.266   0.816  -3.963  1.00  0.00           C  \n"
"HETATM   14  C30 LG1 X 299      -2.358   0.500  -2.981  1.00  0.00           C  \n"
"HETATM   15  C32 LG1 X 299      -1.717  -0.708  -2.627  1.00  0.00           C  \n"
"HETATM   16  C33 LG1 X 299      -1.767  -2.046  -3.084  1.00  0.00           C  \n"
"HETATM   17  C34 LG1 X 299      -1.044  -2.996  -2.549  1.00  0.00           C  \n"
"HETATM   18  C31 LG1 X 299      -0.120  -2.762  -1.426  1.00  0.00           C  \n"
"HETATM   19  C2  LG1 X 299       0.797  -3.520  -0.676  1.00  0.00           C  \n"
"HETATM   20  C7  LG1 X 299       1.182  -4.841  -0.625  1.00  0.00           C  \n"
"HETATM   21  C3  LG1 X 299       2.110  -5.264   0.279  1.00  0.00           C  \n"
"HETATM   22  C5  LG1 X 299       2.676  -4.369   1.218  1.00  0.00           C  \n"
"HETATM   23  C9  LG1 X 299       2.311  -3.043   1.224  1.00  0.00           C  \n"
"HETATM   24  C35 LG1 X 299       1.385  -2.627   0.266  1.00  0.00           C  \n"
"HETATM   25  N4  LG1 X 299       0.877  -1.390   0.070  1.00  0.00           N  \n"
"HETATM   26  C37 LG1 X 299      -0.099  -1.414  -0.954  1.00  0.00           C  \n"
"HETATM   27  C38 LG1 X 299      -0.838  -0.448  -1.507  1.00  0.00           C  \n"
"HETATM   28  C40 LG1 X 299      -1.268  -4.294  -3.262  1.00  0.00           C  \n"
"HETATM   29  N3  LG1 X 299      -2.169  -3.976  -4.202  1.00  0.00           N  \n"
"HETATM   30  C39 LG1 X 299      -2.600  -2.594  -4.258  1.00  0.00           C  \n"
"HETATM   31  O3  LG1 X 299      -0.769  -5.381  -3.040  1.00  0.00           O  \n"
"HETATM   32  O5  LG1 X 299       1.755   1.144  -1.649  1.00  0.00           O  \n"
"HETATM   33  C48 LG1 X 299       1.269   1.452  -2.941  1.00  0.00           C  \n"
"HETATM   34  N6  LG1 X 299       3.559   2.588   0.319  1.00  0.00           N  \n"
"HETATM   35  C47 LG1 X 299       4.002   3.352   1.518  1.00  0.00           C  \n"
);
		NameBimap name_map;
		remap_names_on_geometry(name_map, rinfo, *rsd_);
		TS_ASSERT_EQUALS( name_map.left.size(), rinfo.atoms.size());
		TS_ASSERT(name_map.left.count(" O4 "));
		TS_ASSERT(name_map.left.count(" N6 "));
		TS_ASSERT(name_map.left.count(" C47"));
		TS_ASSERT(name_map.left.count(" H36") == 0);
		TS_ASSERT_EQUALS(name_map.left.find(" O4 ")->second," O2 ");
		TS_ASSERT_EQUALS(name_map.left.find(" N6 ")->second," N4 ");
		TS_ASSERT_EQUALS(name_map.left.find(" C47")->second," C27");
	}

	/// @brief Tests the ability of the matcher to make a symmetry item with a single matching name
	void test_sym_one_name_map() {

		TR << "Test symmetry and one name" << std::endl;
		ResidueInformation rinfo = create_ResidueInfo( // Name mismatches first digit gets incremented by two
"HETATM    7  C46 LG1 X 299      -0.872   2.730   0.452  1.00  0.00           C  \n"
"HETATM    8  N5  LG1 X 299      -0.965   0.941  -1.242  1.00  0.00           N  \n"
"HETATM    9  C36 LG1 X 299      -1.883   1.537  -2.121  1.00  0.00           C  \n"
"HETATM   10  C51 LG1 X 299      -2.338   2.811  -2.270  1.00  0.00           C  \n"
"HETATM   11  C56 LG1 X 299      -3.258   3.103  -3.234  1.00  0.00           C  \n"
"HETATM   12  C54 LG1 X 299      -3.747   2.094  -4.117  1.00  0.00           C  \n"
"HETATM   13  C58 LG1 X 299      -3.266   0.816  -3.963  1.00  0.00           C  \n"
"HETATM   14  C10 LG1 X 299      -2.358   0.500  -2.981  1.00  0.00           C  \n" // Correctly named item.
"HETATM   15  C32 LG1 X 299      -1.717  -0.708  -2.627  1.00  0.00           C  \n"
"HETATM   18  C31 LG1 X 299      -0.120  -2.762  -1.426  1.00  0.00           C  \n"
"HETATM   19  C52 LG1 X 299       0.797  -3.520  -0.676  1.00  0.00           C  \n"
"HETATM   20  C57 LG1 X 299       1.182  -4.841  -0.625  1.00  0.00           C  \n"
"HETATM   21  C53 LG1 X 299       2.110  -5.264   0.279  1.00  0.00           C  \n"
"HETATM   22  C55 LG1 X 299       2.676  -4.369   1.218  1.00  0.00           C  \n"
"HETATM   23  C59 LG1 X 299       2.311  -3.043   1.224  1.00  0.00           C  \n"
"HETATM   24  C35 LG1 X 299       1.385  -2.627   0.266  1.00  0.00           C  \n"
"HETATM   25  N4  LG1 X 299       0.877  -1.390   0.070  1.00  0.00           N  \n"
"HETATM   26  C37 LG1 X 299      -0.099  -1.414  -0.954  1.00  0.00           C  \n"
"HETATM   27  C38 LG1 X 299      -0.838  -0.448  -1.507  1.00  0.00           C  \n"
);
		NameBimap name_map;
		remap_names_on_geometry(name_map, rinfo, *rsd_);
		TS_ASSERT_EQUALS( name_map.left.size(), rinfo.atoms.size());
		TS_ASSERT(name_map.left.count(" C10"));
		TS_ASSERT(name_map.left.count(" C58"));
		TS_ASSERT(name_map.left.count(" C57"));
		TS_ASSERT(name_map.left.count(" N5 "));
		TS_ASSERT_EQUALS(name_map.left.find(" C10")->second," C10");
		TS_ASSERT_EQUALS(name_map.left.find(" C58")->second," C6 ");
		TS_ASSERT_EQUALS(name_map.left.find(" C57")->second," C5 ");
		TS_ASSERT_EQUALS(name_map.left.find(" N5 ")->second," N3 ");
	}

	/// @brief Tests the ability of the matcher to deal with disconnected fragments
	void test_disconnnected_name_map() {

		TR << "Test disconnected mapping" << std::endl;
		ResidueInformation rinfo = create_ResidueInfo( // Name mismatches first digit gets incremented by two
"HETATM    8  N5  LG1 X 299      -0.965   0.941  -1.242  1.00  0.00           N  \n"
"HETATM    9  C36 LG1 X 299      -1.883   1.537  -2.121  1.00  0.00           C  \n"
"HETATM   10  C51 LG1 X 299      -2.338   2.811  -2.270  1.00  0.00           C  \n"
"HETATM   11  C56 LG1 X 299      -3.258   3.103  -3.234  1.00  0.00           C  \n"
"HETATM   12  C54 LG1 X 299      -3.747   2.094  -4.117  1.00  0.00           C  \n"
"HETATM   13  C58 LG1 X 299      -3.266   0.816  -3.963  1.00  0.00           C  \n"
"HETATM   14  C10 LG1 X 299      -2.358   0.500  -2.981  1.00  0.00           C  \n" // Correctly named item.
"HETATM   19  C52 LG1 X 299       0.797  -3.520  -0.676  1.00  0.00           C  \n"
"HETATM   20  C57 LG1 X 299       1.182  -4.841  -0.625  1.00  0.00           C  \n"
"HETATM   21  C53 LG1 X 299       2.110  -5.264   0.279  1.00  0.00           C  \n"
"HETATM   22  C55 LG1 X 299       2.676  -4.369   1.218  1.00  0.00           C  \n"
"HETATM   23  C59 LG1 X 299       2.311  -3.043   1.224  1.00  0.00           C  \n"
"HETATM   24  C35 LG1 X 299       1.385  -2.627   0.266  1.00  0.00           C  \n"
"HETATM   25  N4  LG1 X 299       0.877  -1.390   0.070  1.00  0.00           N  \n"
"HETATM   50  H7  LG1 X 299       3.400  -4.730   1.933  1.00  0.00           H  \n" // Needed to disambiguate the second part
);
		NameBimap name_map;
		remap_names_on_geometry(name_map, rinfo, *rsd_);
		TS_ASSERT_EQUALS( name_map.left.size(), rinfo.atoms.size());
		TS_ASSERT(name_map.left.count(" C10"));
		TS_ASSERT(name_map.left.count(" C58"));
		TS_ASSERT(name_map.left.count(" C57"));
		TS_ASSERT(name_map.left.count(" N5 "));
		TS_ASSERT(name_map.left.count(" N4 "));
		TS_ASSERT(name_map.left.count(" H7 "));
		TS_ASSERT_EQUALS(name_map.left.find(" C10")->second," C10");
		TS_ASSERT_EQUALS(name_map.left.find(" C58")->second," C6 ");
		TS_ASSERT_EQUALS(name_map.left.find(" C57")->second," C5 ");
		TS_ASSERT_EQUALS(name_map.left.find(" N5 ")->second," N3 ");
		TS_ASSERT_EQUALS(name_map.left.find(" N4 ")->second," N2 ");
		TS_ASSERT_EQUALS(name_map.left.find(" H7 ")->second," H3 ");
	}

	/// @brief Appropriate chirality/connectivity should dominate
	void test_chirality_name_map() {

		TR << "Test chirality mapping" << std::endl;
		ResidueInformation rinfo = create_ResidueInfo(
"HETATM    8  N3  LG1 X 299      -0.965   0.941  -1.242  1.00  0.00           N  \n"
"HETATM    1  O2  LG1 X 299      -0.058   0.605   0.915  1.00  0.00           O  \n"
"HETATM    6  C25 LG1 X 299      -0.153   1.537  -0.148  1.00  0.00           C  \n"
"HETATM    5  C26 LG1 X 299       1.238   2.049  -0.655  1.00  0.00           C  \n"
"HETATM    7  C23 LG1 X 299      -0.872   2.730   0.452  1.00  0.00           C  \n"
"HETATM   42  H13 LG1 X 299      -1.906   2.762   0.077  1.00  0.00           H  \n"
"HETATM   40  H16 LG1 X 299       1.170   3.064  -1.074  1.00  0.00           H  \n"
);
		NameBimap name_map;
		remap_names_on_geometry(name_map, rinfo, *rsd_);
		TS_ASSERT_EQUALS( name_map.left.size(), 7);
		TS_ASSERT(name_map.left.count(" N3 "));
		TS_ASSERT(name_map.left.count(" O2 "));
		TS_ASSERT(name_map.left.count(" C25"));
		TS_ASSERT(name_map.left.count(" C26"));
		TS_ASSERT(name_map.left.count(" C23"));
		TS_ASSERT(name_map.left.count(" H13"));
		TS_ASSERT(name_map.left.count(" H16"));
		TS_ASSERT_EQUALS(name_map.left.find(" N3 ")->second," N3 ");
		TS_ASSERT_EQUALS(name_map.left.find(" O2 ")->second," O2 ");
		TS_ASSERT_EQUALS(name_map.left.find(" C25")->second," C25");
		TS_ASSERT_EQUALS(name_map.left.find(" C26")->second," C23");
		TS_ASSERT_EQUALS(name_map.left.find(" C23")->second," C26");
		//TS_ASSERT_EQUALS(name_map.left.find(" H13")->second," H16"); //Actually could any to any of the three hydrogens
		TS_ASSERT_EQUALS(name_map.left.find(" H16")->second," H13");
	}

	/// @brief Test if extra atoms cause issues
	void test_extra_atom_name_map() {

		TR << "Test extra atom mapping" << std::endl;
		ResidueInformation rinfo = create_ResidueInfo(
"HETATM    8  N3  LG1 X 299      -0.965   0.941  -1.242  1.00  0.00           N  \n"
"HETATM    1  O2  LG1 X 299      -0.058   0.605   0.915  1.00  0.00           O  \n"
"HETATM    6  C25 LG1 X 299      -0.153   1.537  -0.148  1.00  0.00           C  \n"
"HETATM    5  C26 LG1 X 299       1.238   2.049  -0.655  1.00  0.00           C  \n"
"HETATM    7  C23 LG1 X 299      -0.872   2.730   0.452  1.00  0.00           C  \n"
"HETATM   42  E13 LG1 X 299      -1.906   2.762   0.077  1.00  0.00           C  \n" // Extra atom - carbon instead of hydrogen
"HETATM   40  H16 LG1 X 299       1.170   3.064  -1.074  1.00  0.00           H  \n"
"HETATM   32  O3  LG1 X 299       1.755   1.144  -1.649  1.00  0.00           O  \n"
);
		NameBimap name_map;
		remap_names_on_geometry(name_map, rinfo, *rsd_);
		TS_ASSERT_EQUALS( name_map.left.size(), 7);
		TS_ASSERT(name_map.left.count(" N3 "));
		TS_ASSERT(name_map.left.count(" O2 "));
		TS_ASSERT(name_map.left.count(" C25"));
		TS_ASSERT(name_map.left.count(" C26"));
		TS_ASSERT(name_map.left.count(" C23"));
		TS_ASSERT(name_map.left.count(" E13") == 0);
		TS_ASSERT(name_map.left.count(" H16"));
		TS_ASSERT(name_map.left.count(" O3 "));
		TS_ASSERT_EQUALS(name_map.left.find(" N3 ")->second," N3 ");
		TS_ASSERT_EQUALS(name_map.left.find(" O2 ")->second," O2 ");
		TS_ASSERT_EQUALS(name_map.left.find(" C25")->second," C25");
		TS_ASSERT_EQUALS(name_map.left.find(" C26")->second," C23");
		TS_ASSERT_EQUALS(name_map.left.find(" C23")->second," C26");
		TS_ASSERT_EQUALS(name_map.left.find(" H16")->second," H13");
		TS_ASSERT_EQUALS(name_map.left.find(" O3 ")->second," O3 ");
	}


	/// @brief Can we determine elements from atom names? (First two letters.)
	void test_no_element_name_map() {

		TR << "Test no element mapping" << std::endl;
		ResidueInformation rinfo = create_ResidueInfo(
"HETATM    8  N3  LG1 X 299      -0.965   0.941  -1.242  1.00  0.00              \n"
"HETATM    1  O2  LG1 X 299      -0.058   0.605   0.915  1.00  0.00\n"
"HETATM    6  C25 LG1 X 299      -0.153   1.537  -0.148  1.00  0.00\n"
"HETATM    5  C26 LG1 X 299       1.238   2.049  -0.655  1.00  0.00\n"
"HETATM    7  C23 LG1 X 299      -0.872   2.730   0.452  1.00  0.00              \n"
"HETATM   42  H13 LG1 X 299      -1.906   2.762   0.077  1.00  0.00              \n"
"HETATM   40  H16 LG1 X 299       1.170   3.064  -1.074  1.00  0.00\n"
);
		NameBimap name_map;
		remap_names_on_geometry(name_map, rinfo, *rsd_);
		TS_ASSERT_EQUALS( name_map.left.size(), 7);
		TS_ASSERT(name_map.left.count(" N3 "));
		TS_ASSERT(name_map.left.count(" O2 "));
		TS_ASSERT(name_map.left.count(" C25"));
		TS_ASSERT(name_map.left.count(" C26"));
		TS_ASSERT(name_map.left.count(" C23"));
		TS_ASSERT(name_map.left.count(" H13"));
		TS_ASSERT(name_map.left.count(" H16"));
		TS_ASSERT_EQUALS(name_map.left.find(" N3 ")->second," N3 ");
		TS_ASSERT_EQUALS(name_map.left.find(" O2 ")->second," O2 ");
		TS_ASSERT_EQUALS(name_map.left.find(" C25")->second," C25");
		TS_ASSERT_EQUALS(name_map.left.find(" C26")->second," C23");
		TS_ASSERT_EQUALS(name_map.left.find(" C23")->second," C26");
		//TS_ASSERT_EQUALS(name_map.left.find(" H13")->second," H16"); //Actually could any to any of the three hydrogens
		TS_ASSERT_EQUALS(name_map.left.find(" H16")->second," H13");
	}
};

