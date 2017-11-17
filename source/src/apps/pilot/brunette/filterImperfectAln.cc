// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/brunette/filterImperfectAln.cc
///
/// @brief  This takes two sets of alignments. 1 perfect and 1 imperfect. The output is the residues of the imperfect alignment that are correctly aligned

/// @usage: -in:file:alignent <imperfect_alignment> -in:file:perfect_alignment <perfect_alignment> -out:file:alignment <output_alignment>

/// @author TJ Brunette

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/id/SequenceMapping.hh>
#include <core/sequence/util.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/Sequence.hh>

#include <devel/init.hh>
//utilities
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>


#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <map>

static basic::Tracer tr( "filterImperfectAln" );

using utility::vector1;
using core::Size;
using std::multimap;
using std::string;
using namespace core::sequence;
using namespace core::id;

namespace filterImperfectAln {
basic::options::FileVectorOptionKey perfect_alignment("filterImperfectAln:perfect_alignment");
}

multimap<string,SequenceAlignment> input_alignmentsMapped(vector1< std::string > align_fns){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	vector1< SequenceAlignment > tmp_alns =core::sequence::read_aln(
		option[ cm::aln_format ](), align_fns[1]);
	multimap<string,SequenceAlignment> alns;
	for ( Size jj = 1; jj <= tmp_alns.size(); ++jj ) {
		string aln_id = tmp_alns[jj].sequence(2)->id();
		string pdbid = aln_id.substr(0,5);
		alns.insert(std::pair<string,SequenceAlignment>(pdbid,tmp_alns[jj]));
	}
	return(alns);
}

void output_alignments(vector1 <SequenceAlignment> alns, std::ostream & out){
	for ( Size ii = 1; ii <= alns.size(); ++ii ) {
		alns[ii].printGrishinFormat(out);
	}
}

SequenceAlignment filterAln(SequenceAlignment perfect_aln, SequenceAlignment orig_aln){

	SequenceMapping perfectMapping( perfect_aln.sequence_mapping( 2, 1 ) );
	SequenceMapping origMapping( orig_aln.sequence_mapping( 2, 1 ) );
	const Size nres1 = orig_aln.sequence(1)->ungapped_sequence().size()+orig_aln.sequence(1)->start()-1;
	const Size nres2 = orig_aln.sequence(2)->ungapped_sequence().size()+orig_aln.sequence(2)->start()-1;
	SequenceMapping filt_mapping(nres1,nres2);
	for ( Size ii=1; ii<=nres2; ++ii ) {
		//create new mapping
		if ( (perfectMapping[ii] == origMapping[ii]) && (origMapping[ii] != 0) ) {
			filt_mapping[origMapping[ii]]=ii;
		}
	}
	SequenceOP orig_aln_sequence1_ = new Sequence(orig_aln.sequence(1)->ungapped_sequence(),orig_aln.sequence(1)->id(),orig_aln.sequence(1)->start());
	SequenceOP orig_aln_sequence2_ = new Sequence(orig_aln.sequence(2)->ungapped_sequence(),orig_aln.sequence(2)->id(),orig_aln.sequence(2)->start());
	SequenceAlignment filt_aln = mapping_to_alignment(filt_mapping,orig_aln_sequence1_,orig_aln_sequence2_);
	return filt_aln;
}

int main( int argc, char * argv [] ) {
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::sequence;
		using utility::file_basename;
		option.add (filterImperfectAln::perfect_alignment, "perfect alignments");
		devel::init(argc, argv);
		vector1< std::string > align_fns = option[ in::file::alignment ]();
		vector1< std::string > perfect_align_fns = option [filterImperfectAln::perfect_alignment]();
		vector1 <SequenceAlignment> out_alns;
		multimap<string,SequenceAlignment> align_map = input_alignmentsMapped(align_fns);
		multimap<string,SequenceAlignment> perfect_align_map = input_alignmentsMapped(perfect_align_fns);
		multimap<string,SequenceAlignment>::iterator itr,foundAln_itr;
		for ( itr = align_map.begin(); itr != align_map.end(); ++itr ) {
			foundAln_itr = perfect_align_map.find(itr->first);
			SequenceAlignment tmp_aln = filterAln(foundAln_itr->second,itr->second);
			out_alns.push_back(tmp_aln);
		}
		string out_filename = option[out::file::alignment ]();
		std::ofstream out_aln_stream( out_filename.c_str() );
		output_alignments(out_alns,out_aln_stream);

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

