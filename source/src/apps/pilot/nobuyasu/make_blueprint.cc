// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/nobuyasu/make_blueprint.cc
/// @brief makes blueprint file
/// @author Nobuyasu Koga ( 02/07/2010 )

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/import_pose/import_pose.hh>
#include <core/sequence/ABEGOManager.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/Sheet.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/BetaAlphaBetaMotif.hh>
#include <protocols/fldsgn/topology/util.hh>
#include <protocols/moves/DsspMover.hh>

#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <ObjexxFCL/format.hh>

#include <fstream>
#include <utility/excn/Exceptions.hh>

static THREAD_LOCAL basic::Tracer TR( "make_blueprint" );

typedef core::Size Size;
typedef std::string String;

class ThisApplication  {
public:
	ThisApplication(){};
	static void register_options();
};

OPT_KEY( File, blue )
OPT_KEY( Integer, abego )
OPT_KEY( Boolean, torsion )

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::s );
	NEW_OPT( blue, "output filename", "default.blueprint" );
	NEW_OPT( abego, "abego output level [ 1-3 ] ", 1 );
	NEW_OPT( torsion, "print torsions ", false );
}


int
main( int argc, char * argv [] )
{
	try{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using core::chemical::oneletter_code_from_aa;
		using protocols::fldsgn::topology::SS_Info2;
		using protocols::fldsgn::topology::SS_Info2_OP;
		using protocols::fldsgn::topology::SheetSet;
		using protocols::fldsgn::topology::SheetSetOP;
		using protocols::fldsgn::topology::StrandPairing;
		using protocols::fldsgn::topology::StrandPairings;
		using protocols::fldsgn::topology::StrandPairingSet;
		using protocols::fldsgn::topology::StrandPairingSetOP;
		using protocols::fldsgn::topology::BetaAlphaBetaMotifSet;

		ThisApplication::register_options();
		devel::init(argc, argv);

		// blueprint output file
		std::ofstream output;
		std::ostringstream filename;
		filename <<  option[ blue ]();
		output.open( filename.str().c_str() ,std::ios::out );

		// calc secondary structure info
		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, option[ in::file::s ].value().at( 1 ) );
		protocols::moves::DsspMover dsm;
		dsm.apply( pose );

		// out secondaray structure info
		SS_Info2_OP ssinfo( new SS_Info2( pose, pose.secstruct() ) );
		output << *ssinfo;

		// calc strand pairing set
		StrandPairingSetOP spairset( new StrandPairingSet( protocols::fldsgn::topology::calc_strand_pairing_set( pose, ssinfo ) ) );
		output << *spairset;
		output << "SSPAIR " << spairset->name() << std::endl;

		// calc sheet
		SheetSetOP sheet_set( new SheetSet( ssinfo, spairset ) );
		output << *sheet_set;

		// calc bab
		BetaAlphaBetaMotifSet bab( ssinfo, sheet_set );
		output << bab;


		core::sequence::ABEGOManager am;
		Size level( option[ abego ]() );
		utility::vector1< std::string > abego = core::sequence::get_abego( pose, level );

		if ( level <= am.alllevel() ) {
			for ( Size i=1; i<= pose.total_residue(); i++ ) {
				// Size abego = am.torsion2index( pose.phi( i ), pose.psi( i ), pose.omega( i ), level );
				output << i << " " << oneletter_code_from_aa( pose.aa( i ) ) << " " << pose.secstruct( i ) << abego[ i ] << " ." << std::endl;
			}
		} else {
			for ( Size i=1; i<= pose.total_residue(); i++ ) {
				output << i << " " << oneletter_code_from_aa( pose.aa( i ) ) << " " << pose.secstruct( i ) << " ." << std::endl;
			}
		}

		using namespace ObjexxFCL::format;
		if ( option[ torsion ]() ) {
			output << "## TORSION ANGLES" << std::endl;
			for ( Size ii=1; ii<=pose.total_residue(); ii++ ) {
				output << "# " << ii << " " << pose.secstruct( ii ) << " " << abego[ ii ] << " "
					<< F( 8, 3, pose.phi( ii ) ) << F( 8, 3, pose.psi( ii ) ) << F( 8, 3, pose.omega( ii ) )
					<< std::endl;
			}
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

