// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

#include <devel/init.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>

#include <devel/init.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>

#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <core/kinematics/Jump.hh>

#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <protocols/hybridization/util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loops.hh>

// C++ headers
#include <iostream>
#include <fstream>
#include <string>

OPT_KEY( String, pdb )
OPT_KEY( IntegerVector, fragsizes )
OPT_KEY( IntegerVector, fragamounts )

using namespace ObjexxFCL::format;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::pose;

static basic::Tracer TR("grower_prep");

utility::vector1<protocols::loops::Loop>
get_unaligned( core::id::SequenceMapping const & sequencemap )
{
	utility::vector1<protocols::loops::Loop> loops;

	int startpos=0;
	bool inloop=false;
	for ( unsigned int currentpos = 1; currentpos <= sequencemap.size1(); ++currentpos ) {

		bool ismapped=false;
		if ( sequencemap[currentpos]==0 ) {
			ismapped = false;
		} else {
			ismapped = true;
		}
		if ( !inloop && !ismapped ) {
			startpos = currentpos;
		}

		if ( inloop && ismapped ) {
			loops.push_back(protocols::loops::Loop(startpos,currentpos-1));
		}
		inloop = !ismapped;
	}
	if ( inloop ) {
		loops.push_back(protocols::loops::Loop(startpos, sequencemap.size1()));
	}
	return loops;
}

void Prepare(){
	bool haspdb = false;
	bool hasfasta = false;
	core::pose::Pose pose;
	core::sequence::SequenceOP sequence;
	utility::vector1< core::Size > cbreaks;

	//Get the pose
	if ( option[pdb].user() ) {
		core::import_pose::pose_from_file( pose, option[ pdb ] );
		haspdb = true;
	}
	//Get the chainbreaks and the clean sequence
	if ( option[ in::file::fasta].user() ) {
		sequence = core::sequence::read_fasta_file( option[ in::file::fasta ]()[1])[1];

		// read chainbreaks
		std::string fulllength_clean;
		std::string seq = sequence->sequence();
		for ( int i=0; i<(int)seq.length(); ++i ) {
			if ( seq[i] == '/' ) {
				cbreaks.push_back(fulllength_clean.length());
			} else {
				fulllength_clean += seq[i];
			}
		}
		sequence = core::sequence::SequenceOP(new core::sequence::Sequence( fulllength_clean, "target" ));
		hasfasta = true;
	}
	//Get the Fragsizes and respective amounts
	utility::vector1< core::Size > frag_sizes;
	utility::vector1< core::Size > frag_amounts;
	if ( option[fragsizes].user() ) {
		frag_sizes = option[fragsizes]();
		if ( option[fragamounts].user() ) {
			frag_amounts = option[fragamounts]();
		} else {
			TR << " You specificed the size of the fragments you want to pick fragments but didn't list how many you want. Exiting the protocol" << std::endl;
			exit(0);
		}
		if ( frag_sizes.size() != frag_amounts.size() ) {
			TR << " The list of frag sizes doesn't equal the list of the respective amounts. Exiting the protocol" << std::endl;
			exit(0);
		}
	}
	//Pick the fragments
	if ( hasfasta and frag_sizes.size() > 0 ) {
		for ( core::uint i=1; i<=frag_sizes.size(); ++i ) {
			core::fragment::FragSetOP fragments = protocols::hybridization::create_fragment_set_no_ssbias(sequence->sequence(), frag_sizes[i], frag_amounts[i]);
			std::string totalfrags = utility::to_string(frag_amounts[i]);
			std::string fraglength = utility::to_string(frag_sizes[i]);
			std::string fragname = totalfrags+"."+fraglength+"mers";
			core::fragment::FragmentIO().write_data( fragname, *(fragments) );
		}
	}
	//do the alignment
	core::Size total_loops = 0;
	if ( haspdb && hasfasta ) {
		core::sequence::SequenceOP t_pdb_seq( new core::sequence::Sequence( pose.sequence(), "pose_seq" ));
		core::sequence::SWAligner sw_align;
		core::sequence::ScoringSchemeOP ss(  new core::sequence::SimpleScoringScheme(120, 0, -100, 0));
		core::sequence::SequenceAlignment fasta2template;

		fasta2template = sw_align.align(sequence, t_pdb_seq, ss);
		core::id::SequenceMapping sequencemap = fasta2template.sequence_mapping(1,2);

		TR << "Aligning:" << std::endl;
		TR << fasta2template << std::endl;
		for ( core::Size ii = 1; ii <= sequence->sequence().size(); ++ii ) {
			if ( sequencemap[ii] ) {
				TR << sequence->at(ii);
			} else {
				TR << "-";
			}
		}
		TR << std::endl;

		utility::vector1<protocols::loops::Loop> loops_from_aln = get_unaligned(sequencemap), loops_to_build;


		// split chainbreaks
		for ( int i=1; i<=(int)loops_from_aln.size() ; ++i ) {
			bool loop_was_split=false;
			for ( int j=1; j<=(int)cbreaks.size(); ++j ) {
				if ( loops_from_aln[i].start() <= cbreaks[j] && loops_from_aln[i].stop() > cbreaks[j] ) {
					runtime_assert( !loop_was_split );
					loop_was_split = true;
					loops_to_build.push_back(protocols::loops::Loop(loops_from_aln[i].start(), cbreaks[j]));
					loops_to_build.push_back(protocols::loops::Loop(cbreaks[j]+1, loops_from_aln[i].stop()));
				}
			}
			if ( !loop_was_split ) {
				loops_to_build.push_back(loops_from_aln[i]);
			}
		}
		total_loops = loops_to_build.size();
	}
	//report the total loops and fragments
	std::ofstream outfile;
	outfile.open("grower_prep.txt");
	if ( total_loops != 0 ) {
		outfile << "total loops = " << utility::to_string(total_loops) << std::endl;
	}
	for ( core::Size i=1; i<=frag_sizes.size(); i++ ) {
		outfile << "fragments = " << utility::to_string(frag_amounts[i]) << "." << utility::to_string(frag_sizes[i]) << "mers" << std::endl;
	}
	outfile.close();
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		NEW_OPT( pdb, "The input pdb you will be using for the grower", "");
		NEW_OPT( fragsizes, "The size of the fragments to pick", utility::vector1<int>());
		NEW_OPT( fragamounts, "How many fragments of each set do you want", utility::vector1<int>());

		devel::init( argc, argv );
		Prepare();
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
