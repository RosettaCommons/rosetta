// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/r_frag_quality.cc
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange

// Project headers
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragID_Iterator.hh>
#include <core/fragment/BBTorsionSRFD.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>


#include <core/fragment/util.hh>
#include <protocols/jumping/JumpSample.hh>
#include <protocols/jumping/JumpSetup.hh>
#include <protocols/jumping/SecondaryStructure.hh>
#include <utility/excn/Exceptions.hh>

#include <core/scoring/rms_util.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <core/pose/Pose.hh>
#include <devel/init.hh>

#include <core/chemical/ChemicalManager.hh>


#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

// Utility Headers
#include <utility/io/izstream.hh>

#include <numeric/angle.functions.hh>
#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>


#include <map>


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>


using namespace core;
using namespace fragment;
using namespace pose;
using namespace kinematics;

using namespace protocols::jumping;


using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace ObjexxFCL::format;


static basic::Tracer tr( "main" );


class ThisApplication  {
public:
	ThisApplication();
	void run();
	void analyze_tree( kinematics::FoldTree f, SecondaryStructureOP );
	void show_jumps( std::ostream& );
	void show_cuts( std::ostream& );
	static void register_options();
private:
	typedef std::pair< Size, Size > Pairing;
	typedef std::map< Pairing, Size > Jumps;
	typedef std::map< Size, Size > Cuts;

	Jumps jumps_;
	Cuts cuts_;
};

ThisApplication::ThisApplication()
{}

OPT_KEY( File, fold_trees )
OPT_KEY( Boolean, regen_tree )
OPT_KEY( File, pymol_pic )

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( fold_trees, "file of FoldTrees to be analysed", "" );
	NEW_OPT( regen_tree, "tree (cuts) will be regenerated", false );
	NEW_OPT( pymol_pic, "creates a 'viol' file that shows the jumps in pymol-cnc plugin", "fold_tree.pymol.py");
}

void ThisApplication::analyze_tree( kinematics::FoldTree f, SecondaryStructureOP ss_def ) {
	JumpSample jumps( f );
	if ( option[ pymol_pic ].user() ) jumps.dump_pymol( option[ pymol_pic ] );
	bool bReSample = option[ regen_tree ];
	if ( bReSample ) {
		assert( ss_def );
		kinematics::FoldTree newf;
		Size fails( 0 );
		while ( (fails < 10) && !newf.random_tree_from_jump_points( f.nres(), jumps.size(), jumps.jumps(), ss_def->loop_fraction() ) ) {
			fails++;
		};
		if ( fails>=10 ) {
			tr.Debug<< "unable to re-build tree " << jumps <<  std::endl;
		}
		JumpSample new_jumps( newf );
		jumps = new_jumps;
	}

	bool special( false );
	for ( Size i=1; i<=jumps.size(); i ++ ) {
		jumps_[ Pairing( jumps.jumps() (1,i), jumps.jumps() (2,i) ) ]++;
		cuts_[ jumps.cuts() (i) ]++;
		special |= ( jumps.cuts() (i) ==43 );
		//    cuts_[ new_jumps.cuts() (i) ]++;
	}
	if ( special ) {
		std::cout << jumps << std::endl;
	}
}

void ThisApplication::show_jumps( std::ostream & os ) {
	for ( Jumps::iterator it = jumps_.begin(), eit = jumps_.end(); it!=eit; ++it ) {
		os << RJ(5, it->first.first) << RJ(5, it->first.second) << RJ(20, it->second) << std::endl;
	}
}


void ThisApplication::show_cuts( std::ostream & os ) {
	for ( Cuts::iterator it = cuts_.begin(), eit = cuts_.end(); it!=eit; ++it ) {
		os << RJ(5, it->first) << RJ(20, it->second) << std::endl;
	}
}

void ThisApplication::run() {
	utility::io::izstream data( option[ fold_trees ]() );
	tr.Info << "read Fold-Trees from " << option[ fold_trees ]()  << std::endl;
	if ( !data ) {
		std::cerr << "ERROR:: Unable to open file: "
			<< option[ fold_trees ]() << std::endl;
		std::exit( 1 );
	}

	SecondaryStructureOP ss_def( NULL );
	if ( option[ regen_tree ] ) {
		ConstantLengthFragSet fragset;
		fragset.read_fragment_file( option[ in::file::frag3 ] );
		ss_def = new SecondaryStructure( fragset, false /*no JustUseCentralResidue */ );
	}
	std::string line;
	Size ct( 1 );
	while ( getline(data,line) ) {

		// get rid of the tag at end of line
		Size found = line.find("  S_");
		if ( found != std::string::npos ) {
			line.erase( found );
		}

		std::istringstream line_stream( line );

		FoldTree f;
		line_stream >> f;
		//std::cout << f;
		analyze_tree( f, ss_def );
		if ( numeric::mod( (int) ct++, 100) == 0 ) {
			tr.Info << "read " << ct << std::endl;
		}
	}

	std::cout << "write JUMP Histogram to file jumps.dat" << std::endl;
	std::ofstream jout( "jumps.dat" );
	show_jumps( jout );

	std::ofstream ctsout( "cuts.dat" );
	std::cout << "CUT Histogram to file cuts.dat" << std::endl;
	show_cuts( ctsout );
}

int main( int argc, char** argv ) {
	try{
		ThisApplication::register_options();
		devel::init( argc, argv );

		ThisApplication app;
		app.run();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}


}
