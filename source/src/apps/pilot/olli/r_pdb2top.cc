// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange

#include <protocols/abinitio/Templates.hh>
#include <protocols/abinitio/PairingStatistics.hh>

#include <core/pose/Pose.hh>
#include <devel/init.hh>

#include <core/scoring/dssp/PairingsList.hh>
#include <core/scoring/dssp/StrandPairing.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/silent.fwd.hh>


// Utility headers
#include <basic/options/option_macros.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <utility/excn/Exceptions.hh>
// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>


static THREAD_LOCAL basic::Tracer tr( "main" );

using namespace core;
using namespace core::scoring::dssp;
using namespace protocols;
using namespace abinitio;
//using namespace jumping;
using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_1GRP_KEY( File, out, top )
OPT_1GRP_KEY( File, in, top )
OPT_1GRP_KEY( File, in, ref_top )
OPT_KEY( Integer, ntest )

int main( int argc, char** argv ) {
	try{
		Templates::register_options();
		OPT( in::file::s );
		NEW_OPT( out::top, "write topology info to this file", "");
		NEW_OPT( in::top, "read topology from this file for checking", "");
		NEW_OPT( in::ref_top, "read a reference topology for comparison", "" );
		NEW_OPT( ntest, "perform N cycles of topology selection for top2top scoring", 50 );
		devel::init( argc, argv );
		if ( !basic::options::option[ basic::options::OptionKeys::in::path::database ].user() ) {
			basic::options::option[ basic::options::OptionKeys::in::path::database ].def( "/work/olange/minirosetta_database");
		}

		if ( !basic::options::option[ basic::options::OptionKeys::jumps::max_strand_gap_allowed].user() ) {
			basic::options::option[ basic::options::OptionKeys::jumps::max_strand_gap_allowed].def( 10 );
		}

		if ( !basic::options::option[ basic::options::OptionKeys::jumps::contact_score ].user() ) {
			basic::options::option[ basic::options::OptionKeys::jumps::contact_score ].def( 0.2 );
		}

		if ( option[ in::file::s ].user() ) {
			Pose pose;
			core::import_pose::pose_from_file( pose, option[ in::file::s ]()[ 1 ] , core::import_pose::PDB_file);

			// get strand pairings
			StrandPairingSet strand_pairings( pose );
			//std::cout << strand_pairings << std::endl;

			PairingStatistics ps( strand_pairings );

			if ( option[ out::top ].user() ) {
				utility::io::ozstream file( option[ out::top ] );
				file << ps << std::endl;
			} else {
				std::cout << ps << std::endl;
			}

			if ( tr.Debug.visible() ) {
				PairingList pl;
				strand_pairings.get_beta_pairs( pl );
				tr.Debug << pl << std::endl;
			}
		}
		if ( option[ in::file::silent ].user() ) {

			utility::io::ozstream file( option[ out::top ] );

			//read silent file for input
			core::io::silent::SilentFileData sfd;
			sfd.read_file( *(option [ in::file::silent ]().begin()) );

			// run thru all structures
			Size ct ( 0 );

			PairingStatistics ps;
			PairingStatistics::ModelFreq model_freq;

			for ( core::io::silent::SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it, ++ct ) {
				Pose pose;
				std::string tag = it->decoy_tag();
				it->fill_pose( pose );

				// get strand pairings
				StrandPairingSet strand_pairings( pose );
				tr.Debug << strand_pairings << std::endl;
				tag = ObjexxFCL::right_string_of( ct, 4, '0');
				ps.add_topology( strand_pairings, tag );
				model_freq[ tag.substr(0, 4) ] += 1;
				//screen output
				if ( sfd.size() < 10 ) {
					tr.Info << tag << " " << std::endl;
				} else {
					if ( (ct % 50) == 0 ) {
						std::cout << ".";
						std::cout.flush();
					}
				}
			} // for decoys in silent file
			// set score terms
			ps.compute_model_weights( model_freq );
			file << ps << std::endl;
		}
		if ( option[ in::top ].user() ) {

			utility::io::izstream is( option[ in::top ] );
			PairingStatisticsOP ps( new PairingStatistics );
			is >> *ps;
			tr.Debug << *ps << std::endl;
			std::string model;
			StrandPairingSet const& sp( ps->suggest_topology( model ) );
			if ( tr.Debug.visible() ) {
				PairingList pl;
				sp.get_beta_pairs( pl );
				tr.Debug << pl << std::endl;
			}
			if ( option[ in::ref_top ].user() ) { ///Compute a compatibility score
				utility::io::izstream is( option[ in::ref_top ] );
				PairingStatisticsOP ref_ps( new PairingStatistics );
				is >> *ref_ps;
				std::string ref_model;
				PairingList ref_pl;
				StrandPairingSet const& ref_sp( ref_ps->suggest_topology( ref_model ) );
				ref_sp.get_beta_pairs( ref_pl );

				for ( Size ct=1; ct <= (Size) option[ ntest ]; ct++ ) {
					StrandPairingSet const& sp( ps->suggest_topology( model ) );
					core::Size positives( 0 );
					core::Size false_positives( 0 );
					core::Size ignored_strands( 0 );

					PairingList pl;
					sp.get_beta_pairs( pl );

					for ( PairingList::const_iterator it_ref = ref_pl.begin(); it_ref != ref_pl.end(); ++it_ref ) { //for all native pairings, which can be found in the enforced topology... ?
						PairingList::const_iterator it = std::find( pl.begin(), pl.end(), *it_ref );
						if ( it != pl.end() ) ++positives;
					}

					for ( PairingList::const_iterator it = pl.begin(); it != pl.end(); ++it ) { //for all enforced pairings, which can not be found in the native topology... ?
						PairingList::const_iterator it_ref = std::find( ref_pl.begin(), ref_pl.end(), *it );
						if ( it_ref == ref_pl.end() ) ++false_positives;
					}

					core::Size strand_positives( 0 );
					core::Size strand_false_positives( 0 );
					for ( StrandPairingSet::const_iterator it_ref = ref_sp.begin(); it_ref != ref_sp.end(); ++it_ref ) {
						bool found=false;
						for ( StrandPairingSet::const_iterator it= sp.begin(); it != sp.end() && !found; ++ it ) {
							//check for common pairing:
							found = it->has_common_pairing( *it_ref );
							tr.Trace << *it_ref << " and " << *it << " have " << ( found ? "    " : " no " ) << "common pairing" << std::endl;
						}
						if ( found ) ++strand_positives;
					}

					for ( StrandPairingSet::const_iterator it = sp.begin(); it != sp.end(); ++it ) {
						bool found=false;
						if ( ps->strand_weight( *it )/ps->weight( 1 ) ) {
							for ( StrandPairingSet::const_iterator it_ref= ref_sp.begin(); it_ref != ref_sp.end() && !found; ++ it_ref ) {
								found = it->has_common_pairing( *it_ref ) || it->mergeable( *it_ref );
							}
							if ( !found ) ++strand_false_positives;
						} else {
							++ignored_strands;
						}
					}

					using ObjexxFCL::format::F;
					tr.Info << ref_model << " vs. " << model << "  Positives: " << F( 5, 2, positives*1.0/ref_pl.size() ) << "  " << "FalsePositives: " << F( 5, 2, false_positives*1.0/pl.size() )
						<< " StrandPositives: " << F( 5,2, strand_positives*1.0/ref_sp.size() ) << " " << " StrandFalsePositives " << F( 5,2, strand_false_positives*1.0/(sp.size()-ignored_strands) ) << std::endl;
				}
			}
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
// finish with 0 when there's no error
	return 0;
}
