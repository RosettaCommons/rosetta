// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file score_aln.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <devel/init.hh>

#include <basic/options/option.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/ScoringScheme.fwd.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char* argv [] ) {
	try {
		using std::map;
		using std::string;
		using core::Size;
		using core::Real;
		using utility::vector1;
		using utility::file::FileName;
		using namespace basic;
		using namespace core::sequence;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		// options, random initialization
		devel::init( argc, argv );

		std::string output_fn( option[ out::file::alignment ]() );
		utility::io::ozstream output( output_fn );

		vector1< SequenceAlignment > aln;
		vector1< std::string > align_fns = option[ in::file::alignment ]();

		vector1< SequenceAlignment > aln_map;
		vector1< std::string > align_map_fns = option[ in::file::alignment2 ]();

		vector1< SequenceAlignment > final_aln;


		Size count = 100000;

		std::cout << "Reading alignments: " <<  align_fns[1] << std::endl;
		aln = core::sequence::read_aln( option[ cm::aln_format ](), align_fns[1] );

		final_aln = aln;  // keep original alignments

		std::cout << "Reading structural cross alignments: " << align_map_fns[1] << std::endl;
		aln_map = core::sequence::read_aln( option[ cm::aln_format ](), align_map_fns[1] );

		typedef vector1< SequenceAlignment >::iterator iter;
		typedef vector1< SequenceAlignment >::const_iterator const_iter;
		for ( iter it = aln.begin(), end = aln.end(); it != end; ++it ) {
			std::cout << "Original alignment (A->B): " << std::endl << *it << std::endl;

			SequenceOP seq1_copy = it->sequence(1)->clone();
			SequenceOP seq2_copy = it->sequence(2)->clone();
			seq1_copy->sequence( seq1_copy->ungapped_sequence() );
			seq2_copy->sequence( seq2_copy->ungapped_sequence() );

			std::string idA1 = it->sequence(1)->id();
			std::string idA2 = it->sequence(2)->id();

			core::id::SequenceMapping map1 = it->sequence_mapping ( 1, 2 );

			for ( iter jt = aln_map.begin(),jend = aln_map.end(); jt !=jend; ++jt ) {


				std::string idB1 = jt->sequence(1)->id();
				std::string idB2 = jt->sequence(2)->id();

				if ( idA2.substr(0,5) == idB1.substr(0,5) ) {
					SequenceOP bseq1_copy = jt->sequence(1)->clone();
					SequenceOP bseq2_copy = jt->sequence(2)->clone();
					bseq1_copy->sequence( bseq1_copy->ungapped_sequence() );
					bseq2_copy->sequence( bseq2_copy->ungapped_sequence() );
					core::id::SequenceMapping map2 = jt->sequence_mapping( 1, 2 );

					std::cout << idA1 << "  " << idA2 << "<-->" << idB1 << "   " << idB2 << std::endl;

					std::cout << "Crossmap alignment (B->C): " << std::endl << *jt << std::endl;

					//std::cout << map1 << std::endl;
					//std::cout << "------" << std::endl;
					//std::cout << map2 << std::endl;
					core::id::SequenceMapping map3 = transitive_map( map1, map2 );
					//std::cout << "------" << std::endl;
					//std::cout << map3 << std::endl;

					//std::cout << " seq1_copy:" <<  seq1_copy->sequence() << std::endl;
					//std::cout << " seq1_copy:" <<  seq1_copy->start() << std::endl;
					//std::cout << "bseq2_copy:" << bseq2_copy->sequence() << std::endl;
					//std::cout << " seq1_copy:" << bseq2_copy->start() << std::endl;

					bseq2_copy->id( idB2.substr(0,5) + "_" + utility::to_string( count ) );

					SequenceAlignment cross_aln = mapping_to_alignment(
						map3, seq1_copy, bseq2_copy
					);

					final_aln.push_back( cross_aln );
					count++;
				}
			}
		} // for alns


		Real filter_threshold = -1;
		if ( option[ cm::aln_length_filter_quantile ].user() ) {
			Real quantile = option[ cm::aln_length_filter_quantile ]();
			// make list of lengths
			//create jobs for each alignment
			std::vector <int> length_list;
			for ( const_iter align_it = final_aln.begin(); align_it != final_aln.end(); ++align_it ) {
				length_list.push_back(  align_it->length() - align_it->gapped_positions() );
				std::cout << "Len " <<  align_it->length() - align_it->gapped_positions() << std::endl;
			}
			std::vector<int>::iterator  i = length_list.begin();
			std::vector<int>::size_type m =(size_t)( length_list.size() * quantile );
			std::nth_element(i, i + m, length_list.end());
			filter_threshold = length_list.at(m);
			std::cout << "Quantile filter threshold = " << filter_threshold << std::endl;
		}


		for ( iter it = final_aln.begin(), end = final_aln.end(); it != end; ++it ) {

			Real const alignment_coverage( it->length() - it->gapped_positions() );

			if ( option[ cm::aln_length_filter ].user() ) {
				filter_threshold = option[ cm::aln_length_filter ]();
			}
			if ( ( filter_threshold > 0 ) && ( alignment_coverage <  filter_threshold ) ) {
				std::cout << "Skipping alignment " << it->sequence(1)->id() << "  " << it->sequence(2)->id() << ": length = " << int( alignment_coverage )
					<< "  threshold = " << int( filter_threshold ) << std::endl;
				continue;
			}
			// if alignment passes filter, print it to the output file
			it->printGrishinFormat(output);
		}

		output.close();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}


