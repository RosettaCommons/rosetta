// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/doug/test_chain_sequence.cc
/// @brief Test the PDBInfo::chain_sequence functionality
/// @author P. Douglas Renfrew (doug.renfrew@gmail.com')

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/jd2/JobDistributor.hh>

#include <protocols/simple_moves/ScoreMover.hh>

// core headers
#include <core/import_pose/import_pose.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/util.hh>

#include <core/id/SequenceMapping.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/options/option_macros.hh>

static basic::Tracer TR("test_chain_sequence");


void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( run::preserve_header );
	option.add_relevant( in::file::s );
	option.add_relevant( in::file::l );
}
/// @brief Super simple function to convert three-letter codes to one-letter codeswith a little lbit of logic for stuff Rosetta auto-converts
std::string
convert_to_olc( utility::vector1< std::string > tlc_seq_vec ) {
	using namespace core::chemical;

	std::stringstream olc;

	for ( auto res_tlc : tlc_seq_vec ) {
		if ( res_tlc == "MSE" ) {
			res_tlc = "MET";
		}
		olc << oneletter_code_from_aa( aa_from_one_or_three( res_tlc ) );
	}

	return olc.str();
}

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::import_pose;
		using namespace core::pose;

		devel::init( argc, argv );
		register_options();

		Pose pose;

		if ( option[ in::file::s ].user() ) {
			utility::vector1< std::string > structs = option[ in::file::s ]();

			for ( core::Size i( 1 ); i <= structs.size(); ++i ) {

				core::import_pose::pose_from_file( pose, structs[ i ] );

				// get protein only sequnce of chain 'A'
				utility::vector1< std::string > pose_chain_tlc;
				for ( core::Size i(1); i <= pose.size() ; ++i ) {
					if ( pose.residue_type( i ).is_protein() && pose.pdb_info()->chain( i ) == 'A' ) {
						pose_chain_tlc.push_back( pose.residue_type( i ).name3() );
					}
				}

				std::string pose_seq_str( convert_to_olc( pose_chain_tlc ) );
				std::string seqres_seq_str( convert_to_olc( pose.pdb_info()->chain_sequences( 'A' ) ) );

				std::cout << pose_seq_str.size() << " " <<  pose_seq_str << std::endl;
				std::cout << seqres_seq_str.size() << " " << seqres_seq_str << std::endl;

				core::sequence::SequenceOP seqres_seq( new core::sequence::Sequence( seqres_seq_str, "SEQRES" ) );
				core::sequence::SequenceOP pose_seq( new core::sequence::Sequence( pose_seq_str, "RESOLVED" ) );

				core::id::SequenceMapping seq_map = core::sequence::map_seq1_seq2( seqres_seq, pose_seq );
				core::sequence::SequenceAlignment msa = core::sequence::mapping_to_alignment( seq_map, seqres_seq, pose_seq );

				std::cout << "msa1 " << msa.sequence( 1 )->sequence() << std::endl;
				std::cout << "msa2 " << msa.sequence( 2 )->sequence() << std::endl;
				std::cout << "seq len: " <<  msa.sequence( 1 )->sequence().size() << " " <<  msa.sequence( 2 )->sequence().size() << std::endl;

				std::cout << seq_map << std::endl;

				/*
				for ( core::Size i(1); i <= pose.size() ; ++i ) {
				if ( pose.residue_type( i ).is_protein() && pose.pdb_info()->chain( i ) == 'A' ) {
				std::cout << i << "\t" << pose.pdb_info()->number( i ) << "\t" << seq_map[i] << "\t" << seq_map.get_corresponding_residue_in_current( i ) << std::endl;
				}
				}
				*/
				std::cout << "Size1 Size 2: " << seq_map.size1() << " " << seq_map.size2() << std::endl;
				/*
				/// BACKWARDS
				std::cout << "BACKWARDS BACKWARDS BACKWARDS" << std::endl;

				core::id::SequenceMapping seq_map_bak = core::sequence::map_seq1_seq2( pose_seq, seqres_seq );
				core::sequence::SequenceAlignment msa_bak = core::sequence::mapping_to_alignment( seq_map_bak, pose_seq, seqres_seq );

				std::cout << "msa1 " << msa_bak.sequence( 1 )->sequence() << std::endl;
				std::cout << "msa2 " << msa_bak.sequence( 2 )->sequence() << std::endl;
				std::cout << "seq len: " <<  msa_bak.sequence( 1 )->sequence().size() << " " <<  msa_bak.sequence( 2 )->sequence().size() << std::endl;

				std::cout << seq_map_bak << std::endl;

				std::cout << "Size1 Size 2: " << seq_map_bak.size1() << " " << seq_map_bak.size2() << std::endl;

				//for (core::Size i(1); i <= )

				//auto cs = pose.pdb_info()->chain_sequences( );
				//std::cout << "DOUG cs.size(): " << cs.size() << std::endl;


				//for ( auto j : cs) {
				// std::cout << j << std::endl;
				//}
				*/
			}
		}


	} catch ( utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
