// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Nobuyasu Koga

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID_Mask.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>


// Utility Headers
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>

#include <map>
#include <fstream>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/exit.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/excn/Exceptions.hh>


static basic::Tracer TR( "pick_lowestscore_pose" );

typedef core::Size Size;
typedef core::Real Real;
typedef std::string String;

class ThisApplication  {
public:
	ThisApplication(){};
	static void register_options();
};

OPT_KEY( Integer, top )
OPT_KEY( File, output )
OPT_KEY( Real, bin )
OPT_KEY( Real, max )
OPT_KEY( String, name )

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::silent );
	OPT( in::file::native );
	NEW_OPT( output, "output file name", "lowestscore_pose.out" );
	NEW_OPT( top, "number of structures to pick at each rmsd", 1 );
	NEW_OPT( bin, "bin size of rmsd ", 0.5 );
	NEW_OPT( max, "max rmsd for picking ", 10 );
	NEW_OPT( name, "score name to pick ex, bk_tot, fa_atr....", "score" );
}


int
main( int argc, char * argv [] )
{
	try{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using core::io::silent::SilentFileData;
		using core::io::silent::ProteinSilentStruct;
		using core::scoring::calpha_superimpose_pose;
		using core::scoring::CA_rmsd;
		typedef core::io::silent::SilentFileData::iterator iterator;

		ThisApplication::register_options();
		devel::init(argc, argv);

		// parameters
		Size topX = option[ top ];
		Real rmsbin = option[ bin ];
		Real maxrms = option[ max ];

		// initialization
		Size maxbin = static_cast< Size > (maxrms/rmsbin);

		std::map< Size, std::multimap< Real, String>  > datamap;
		for ( Size i=1; i<=maxbin; i++ ) {
			std::multimap< Real, String > scoretag;
			datamap.insert( std::map< Size, std::multimap< Real, String> >::value_type( i, scoretag ) );
		}

		// read native structure for superimpose
		core::pose::Pose native_pose;
		utility::vector1< Size > positions;
		if ( option[ in::file::native ].user() ) {
			core::import_pose::pose_from_file( native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);
			for ( Size i=1; i<=native_pose.size(); i++ ) {
				positions.push_back( i );
			}
		}

		// read a silent file
		SilentFileData sfd;
		sfd.read_file( option[ in::file::silent ].value().at( 1 ) );

		// iterate thourgh data of silent file
		for ( iterator iter = sfd.begin(), it_end = sfd.end(); iter != it_end; ++iter ) {

			if ( !iter->has_energy( option[ name ]() ) || !iter->has_energy( "rms" ) ) {
				String msg( "Error: can't find score in SilentStuct!" );
				msg += "\nSilentStruct object has scores:\n";
				// little bit of indirection here to get the SilentStruct scores.
				std::ostringstream mystream;
				iter->print_score_header( mystream );
				mystream << std::endl;
				iter->print_scores( mystream );
				mystream << std::endl;
				msg += mystream.str();
				utility_exit_with_message( msg );
			}

			Real rms = iter->get_energy( "rms" );
			if ( rms < maxrms ) {
				Size nbin = static_cast< Size > ( rms/rmsbin + 1 );
				runtime_assert( nbin <= maxbin );

				Real score = iter->get_energy( option[ name ]() );
				String tag = iter->decoy_tag();
				datamap[ nbin ].insert( std::multimap< Real, String >::value_type( score, tag ) );
			}
		}

		// set output
		std::ofstream out2;
		std::ostringstream filename ;
		filename <<  option[ output ]();
		out2.open( filename.str().c_str() ,std::ios::out );

		out2 << "binsize=" << rmsbin << std::endl;
		out2 << "rms " << "decoynum " << "rms_lowestscore_pose " << "score " << "tag " << "output" << std::endl;


		for ( Size i=1; i<=maxbin; i++ ) {
			std::multimap< Real, String > data = datamap[ i ];
			std::multimap< Real, String >::const_iterator ite = data.begin();

			Real rms = ( i-1 )*rmsbin + rmsbin/2;

			Size count( 0 );
			while ( ite != data.end() && count < topX ) {

				out2 << rms << ' ' << data.size() << ' ';

				// set pdb name
				std::ostringstream out;
				out << "rms_" << rms << "_" << std::setw(3) << std::setfill('0') << ++count << ".pdb";

				Real score = ite->first;
				String tag = ite->second;

				core::pose::Pose pose;
				sfd.get_structure( tag ).fill_pose( pose );

				if ( option[ in::file::native ].user() ) {
					// superimpose pose onto native_pose
					calpha_superimpose_pose( pose, native_pose );
					Real rms_calc = CA_rmsd( pose, native_pose );
					// dump pdbs
					utility::io::ozstream file( out.str(), std::ios::out | std::ios::binary);
					file << "REMARK model_1: native, model_2: decoy, rms= " << rms_calc << std::endl;
					core::id::AtomID_Mask mask( true );
					core::pose::initialize_atomid_map( mask, native_pose );
					core::io::pdb::dump_pdb( native_pose, file, mask, "1" );
					core::pose::initialize_atomid_map( mask, pose );
					core::io::pdb::dump_pdb( pose, file, mask, "2" );
					file.close();
					// print data in out2
					out2 << rms_calc << ' ' << score << ' ' << tag << ' ' << out.str() << std::endl;
				} else {
					// dump pdb
					pose.dump_pdb( out.str() );
					// print data in out2
					out2 << rms << ' ' << score << ' ' << tag << ' ' << out.str() << std::endl;
				}
				ite++;
			}
		}
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

