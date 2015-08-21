// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/luki/sc.cc
/// @brief  Calculate shape complementarity of structures.
/// @author Luki Goldschmidt (luki@mbi.ucla.edu)

//core library
#include <basic/options/option.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>

//utilities
#include <devel/init.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


// local options
namespace basic { namespace options { namespace OptionKeys { namespace sc {
basic::options::StringOptionKey const molecule_1("sc:molecule_1");
basic::options::StringOptionKey const molecule_2("sc:molecule_2");
basic::options::BooleanOptionKey const verbose("sc:verbose");
basic::options::BooleanOptionKey const quick("sc:quick");
basic::options::RealOptionKey const density("sc:density");
basic::options::RealOptionKey const rp("sc:rp");
basic::options::RealOptionKey const sep("sc:sec");
basic::options::RealOptionKey const trim("sc:trim");
basic::options::RealOptionKey const weight("sc:weight");
}}}}//basic::options::OptionKeys::sc

using namespace std;

int process_pose(core::pose::Pose &pose, core::scoring::sc::ShapeComplementarityCalculator &sc, std::string &molecule_1, std::string &molecule_2, int verbose);

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// setup
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		option.add(basic::options::OptionKeys::sc::molecule_1, "Chains IDs for molecule 1").def("A");
		option.add(basic::options::OptionKeys::sc::molecule_2, "Chains IDs for molecule 2").def("B");
		option.add(basic::options::OptionKeys::sc::density, "Molecular suface dot density (dots/Angstrom^2)").def(15);
		option.add(basic::options::OptionKeys::sc::quick, "Quick mode (less accurate)").def(false);
		option.add(basic::options::OptionKeys::sc::rp, "Probe radius (Angstrom)").def(1.70);
		option.add(basic::options::OptionKeys::sc::trim, "Trimming distance for the peripheral shell").def(1.5);
		option.add(basic::options::OptionKeys::sc::sep, "Interface separation distance").def(8);
		option.add(basic::options::OptionKeys::sc::verbose, "Show verbose output").def(false);
		option.add(basic::options::OptionKeys::sc::weight, "Weight factor using in sc calculation").def(0.5);

		// initialize core
		devel::init(argc, argv);

		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << " Rosetta Tool:  sc - calculate the Lawrence & Coleman shape complementarity of a PDB file or silent files" << std::endl;
		std::cout << " Usage:" << std::endl;
		std::cout << "   PDB input:     -in:file:s *.pdb  or" << std::endl;
		std::cout << "                  -in:file:l list_of_pdbs" << std::endl;
		std::cout << "   Silent input:  -in:file:silent silent.out            Silent input filesname" << std::endl;
		std::cout << "                  -in:file:s                            Specify specific tags to be extracted, if left out all will be taken" << std::endl;
		std::cout << "                  -in:file:fullatom                     For full atom structures" << std::endl;
		std::cout << "                  -in:file:silent_struct_type <type>    Specify the input silent-file format" << std::endl;
		std::cout << "   Sc:            -sc:molecule_1 <string>               Chain IDs for molecule 1" << std::endl;
		std::cout << "                  -sc:molecule_2 <string>               Chain IDs for molecule 2" << std::endl;
		std::cout << "                  -sc:verbose                           Display verbose output" << std::endl;
		std::cout << "                  -sc:density <float>                   Molecular dot density (dots/A^2)" << std::endl;
		std::cout << "                  -sc:trim <float>                      Trimming distance for the peripheral shell (A)" << std::endl;
		std::cout << "                  -sc:rp <float>                        Probe radius (A)" << std::endl;
		std::cout << "                  -sc:sep <float>                       Interface separation threshold (A)" << std::endl;
		std::cout << "                  -sc:weight <float>                    Distance weighting factor" << std::endl;
		std::cout << " Example: " << std::endl;
		std::cout << "   sc -database ~/rosetta_database -sc:verbose -s 2g38.pdb" << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// setup up and run
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();
		core::Size count = 0;

		if ( !input.has_another_pose() ) {
			return 0;
		}

		core::scoring::sc::ShapeComplementarityCalculator sc;

		sc.settings.density = option[ basic::options::OptionKeys::sc::density ];
		sc.settings.rp = option[ basic::options::OptionKeys::sc::rp ];
		sc.settings.sep = option[ basic::options::OptionKeys::sc::sep ];
		sc.settings.band = option[ basic::options::OptionKeys::sc::trim ];
		sc.settings.weight = option[ basic::options::OptionKeys::sc::weight ];

		if ( option[ basic::options::OptionKeys::sc::quick ] ) {
			sc.settings.density = 5.0;
		}

		if ( !sc.Init() ) {
			return 1;
		}

		core::chemical::ResidueTypeSetCOP rsd_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" ) );
		std::string molecule_1 = option[ basic::options::OptionKeys::sc::molecule_1 ];
		std::string molecule_2 = option[ basic::options::OptionKeys::sc::molecule_2 ];

		while ( input.has_another_pose() ) {

			if ( input.has_another_pose() ) {
				std::cout << std::endl;
				std::cout << "Structure " << (count+1) << ":" << std::endl;
				std::cout << std::endl;
			}

			core::pose::Pose pose;
			input.fill_pose( pose, *rsd_set );

			if ( process_pose(pose, sc, molecule_1, molecule_2, option[ basic::options::OptionKeys::sc::verbose ]) ) {
				count++;
			}
		}

		if ( count != 1 ) {
			std::cout << count << " structures processed." << std::endl;
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

int process_pose(
	core::pose::Pose &pose,
	core::scoring::sc::ShapeComplementarityCalculator &sc,
	std::string &molecule_1,
	std::string &molecule_2,
	int verbose
)
{
	sc.Reset();
	core::pose::PDBInfoCOP pdb_info = pose.pdb_info();

	// Split PDB into two surfaces
	for ( core::Size i = 1; i <= pose.n_residue(); i++ ) {
		char chain = pdb_info->chain(i);
		if ( molecule_1.find(chain) != string::npos ) {
			sc.AddResidue(0, pose.residue(i));
		} else if ( molecule_2.find(chain) != string::npos ) {
			sc.AddResidue(1, pose.residue(i));
		}
	}

	if ( !sc.Calc() ) {
		return 0;
	}

	core::scoring::sc::RESULTS const &r = sc.GetResults();

	if ( verbose ) {

		// Verbose view
		std::cout << "==================================================" << endl;
		std::cout << endl;
		for ( int i = 0; i <= 2; i++ ) {
			if ( i < 2 ) {
				std::cout << "Molecule " << (i+1) << ":" << endl;
			} else {
				std::cout << "Total/Average for both molecules:" << endl;
			}

			std::cout << "\t  Total Atoms: " << r.surface[i].nAtoms << endl;
			std::cout << "\t Buried Atoms: " << r.surface[i].nBuriedAtoms << endl;
			std::cout << "\tBlocked Atoms: " << r.surface[i].nBlockedAtoms << endl;
			std::cout << "\t   Total Dots: " << r.surface[i].nAllDots << endl;
			//        std::cout << "   Buried Dots: " << r.surface[i].nBuriedDots << endl;
			//        std::cout << "      Accessible Dots: " << r.surface[i].nAccessibleDots << endl;
			std::cout << " Trimmed Surface Dots: " << r.surface[i].nTrimmedDots << endl;
			std::cout << "\t Trimmed Area: " << r.surface[i].trimmedArea << endl;
			std::cout << endl;
		}
		std::cout << endl;

		for ( int i = 0; i <= 2; i++ ) {
			if ( i < 2 ) {
				std::cout << "Molecule " << (i+1) << "->" << ((i+1)%2+1) << ": " << endl;
			} else {
				std::cout << "Average for both molecules:" << endl;
			}
			std::cout << "      Mean Separation: " << r.surface[i].d_mean << endl;
			std::cout << "    Median Separation: " << r.surface[i].d_median << endl;
			std::cout << "    Mean Shape Compl.: " << r.surface[i].s_mean << endl;
			std::cout << "  Median Shape Compl.: " << r.surface[i].s_median << endl;
			std::cout << endl;
		}
	}

	std::cout << "==================================================" << endl;
	std::cout << "Shape Complementarity:          " << r.sc << endl;
	std::cout << "Interface separation (A):       " << r.distance << endl;
	std::cout << "Area buried in interface (A^2): " << r.area << endl;
	std::cout << "==================================================" << endl;

	return 1;
}
