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
/// @author Frank DiMaio
/// @author Srivatsan Raman
#include <protocols/rbsegment_relax/RBSegment.hh>

// Rosetta Headers
//#include <core/scoring/ScoringManager.hh>
#include <basic/Tracer.hh>

// Random number generator
#include <numeric/random/random.hh>


#include <string>
#include <fstream>

#include <core/id/SequenceMapping.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>
#include <map>


namespace protocols {
namespace rbsegment_relax {

using namespace core;

static THREAD_LOCAL basic::Tracer TR_seg( "RBSegment" );


/////////////////////
/// @brief helper function
inline core::Real square( core::Real X ) {
	return X*X;
}

//////////////////////////////////////////////////////////
/// @brief Helper function tokenizes a str
/////////////////////////////////////////////////////////
// don't use use string_utils.hh:split or string_split instead
void Tokenize(const std::string              &str,
	utility::vector1<std::string>  &tokens,
	const std::string              &delimiters=" ") {
	tokens.clear();
	std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

	while ( std::string::npos != pos || std::string::npos != lastPos ) {
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		lastPos = str.find_first_not_of(delimiters, pos);
		pos = str.find_first_of(delimiters, lastPos);
	}
}


void RBSegment::set_movement( core::Real sigAxisR, core::Real sigAxisT, core::Real sigOffAxisR, core::Real sigOffAxisT) {
	if ( sigOffAxisR == 0 && sigOffAxisT == 0 ) {
		TR_seg << "Setting movement parameters: [" << segments_[1].start() << "...] :"
			<< sigAxisR << "  ,  " << sigAxisT << std::endl;
	} else {
		TR_seg << "Setting movement parameters: [" << segments_[1].start() <<  "...] :"
			<< sigAxisR << "  ,  " << sigAxisT  << "  ,  " << sigOffAxisR  << "  ,  " << sigOffAxisT << std::endl;
	}
	sigAxisR_ = sigAxisR;
	sigAxisT_ = sigAxisT;
	sigOffAxisR_ = sigOffAxisR;
	sigOffAxisT_ = sigOffAxisT;
}

void RBSegment::get_movement( core::Real &sigAxisR, core::Real &sigAxisT, core::Real &sigOffAxisR, core::Real &sigOffAxisT) const {
	sigAxisR = sigAxisR_;
	sigAxisT = sigAxisT_;
	sigOffAxisR = sigOffAxisR_;
	sigOffAxisT = sigOffAxisT_;
}


//////////////////////////////////////////////////////////
/// @brief Parses an RB segment file into a vector of RBsegments
/////////////////////////////////////////////////////////
void read_RBSegment_file(
	utility::vector1< RBSegment > &rbsegs,
	protocols::loops::Loops &loops,
	std::string filename,
	bool autoGenerateLoops /* = false */,
	int nres /* =0 */,
	utility::vector1< core::Size > cutpoints /*=utility::vector1< core::Size >(0)*/
) {
	rbsegs.clear();
	loops.clear();

	std::ifstream infile( filename.c_str() );

	if ( !infile.good() ) {
		TR_seg << "[ERROR] Error opening RBSeg file '" << filename << "'" << std::endl;
		exit(1);
	}

	utility::vector1< utility::vector1< int > > locked_segments;
	utility::vector1< utility::vector1< int > > compound_segments;
	utility::vector1< std::pair< core::Real, core::Real > > locked_mvmt;
	utility::vector1< std::pair< core::Real, core::Real > > compound_mvmt;
	std::map < int, RBSegment > simple_segments;
	utility::vector1 < RBSegment > simple_seg_list;

	std::string line;
	utility::vector1< std::string > tokens;

	while ( getline(infile,line) ) {
		Tokenize( line, tokens ) ;

		if ( tokens.size() > 0 ) {
			/////////////////////////////
			//// segment (simple)
			if ( tokens[1] == "SEGMENT" || tokens[1] == "SEG" ) {
				if ( tokens.size() < 5 ) {
					TR_seg << "[ERROR] Error parsing line '" << line << "'" << std::endl;
					exit(1);
				}
				int id = atoi( tokens[2].c_str() );
				if ( simple_segments.find( id ) != simple_segments.end() ) {
					TR_seg << "[ERROR] Segment id " << id << " defined twice in RBSeg file!" << std::endl;
					exit(1);
				}
				int seg_start=atoi( tokens[3].c_str() ), seg_end = atoi( tokens[4].c_str() );
				char ss_key = tokens[5][0];
				simple_segments[id] = RBSegment( seg_start, seg_end, ss_key );

				if ( tokens.size() == 7 ) {
					simple_segments[id].set_movement( atof( tokens[6].c_str() ),
						atof( tokens[7].c_str() ) );
				} else if ( tokens.size() == 9 ) {
					simple_segments[id].set_movement( atof( tokens[6].c_str() ),
						atof( tokens[7].c_str() ),
						atof( tokens[8].c_str() ),
						atof( tokens[9].c_str() ) );
				}

				simple_seg_list.push_back( simple_segments[id] );
				/////////////////////////////
				//// lock
			} else if ( tokens[1] == "LOCKED" ) {
				if ( tokens.size() < 5 ) {
					TR_seg << "[ERROR] Error parsing line '" << line << "'" << std::endl;
					exit(1);
				}
				utility::vector1< int > thisLock;
				for ( core::Size i=4; i<=tokens.size(); ++i ) {
					thisLock.push_back(  atoi( tokens[i].c_str() ) );
				}
				locked_segments.push_back( thisLock );

				locked_mvmt.push_back(
					std::pair< core::Real, core::Real > ( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ) ) );
				/////////////////////////////
				//// compound
			} else if ( tokens[1] == "COMPOUND" ) {
				/*
				if ( tokens.size() < 5 ) {
				TR_seg << "[ERROR] Error parsing line '" << line << "'" << std::endl;
				exit(1);
				}
				utility::vector1< int > thisCompound;
				for (core::Size i=4; i<=tokens.size(); ++i)
				thisCompound.push_back(  atoi( tokens[i].c_str() ) );
				compound_segments.push_back( thisCompound );

				compound_mvmt.push_back(
				std::pair< core::Real, core::Real > ( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ) ) );
				*/
				TR_seg.Error << "[ERROR] compund segments unsupported ... ignoring line '" << line << "'" << std::endl;
			} else if ( tokens[1] == "LOOP" ) {
				if ( autoGenerateLoops ) {
					TR_seg.Warning << "[WARNING] LOOP specified in RBSegment file but autoGenerateLoops is true!"
						<< " Ignoring line." << std::endl;
					continue;
				}

				if ( tokens.size() < 3 ) {
					TR_seg << "[ERROR] Error parsing line '" << line << "'" << std::endl;
					exit(1);
				}
				core::Size start_res = (core::Size) atoi(tokens[2].c_str());
				core::Size end_res   = (core::Size) atoi(tokens[3].c_str());
				core::Size cutpt = 0;        // default - let LoopRebuild choose cutpoint
				core::Real skip_rate = 0.0;  // default - never skip
				std::string extend_loop_str;
				bool extend_loop = false;

				if ( tokens.size() > 3 ) {
					cutpt = (core::Size) atoi(tokens[4].c_str());
				}
				if ( tokens.size() > 4 ) {
					skip_rate = atof(tokens[5].c_str());
				}
				if ( tokens.size() > 5 ) {
					extend_loop_str = tokens[6];
				}

				if ( extend_loop_str.length() > 0 ) {
					extend_loop = true;
				}

				loops.push_back( protocols::loops::Loop(start_res, end_res, cutpt, skip_rate, extend_loop) );
			} else if ( tokens[1][0] != '#' ) {
				TR_seg << "[WARNING] Skipping line '" << line << "'" << std::endl;
			}
		}
	}

	// auto-gen loops (if necessary)
	if ( autoGenerateLoops && simple_seg_list.size() > 0 ) {
		std::sort( simple_seg_list.begin(), simple_seg_list.end(), RB_lt());
		int start_res=1, end_res=simple_seg_list[1][1].start()-1;
		//int cutpt = (start_res+end_res)/2;
		int nsegs = simple_seg_list.size();

		if ( end_res > start_res ) {
			loops.push_back( protocols::loops::Loop(start_res, end_res, 0, 0.0, false) );
		} else if ( end_res == start_res ) {
			loops.push_back( protocols::loops::Loop(start_res, end_res/*+1*/, 0, 0.0, false) );
		}
		for ( int i=1; i<nsegs; ++i ) {
			start_res = simple_seg_list[i][1].end()+1;
			end_res   = simple_seg_list[i+1][1].start()-1;
			//cutpt = (start_res+end_res)/2;  // set but never used ~Labonte
			if ( end_res > start_res ) {
				loops.push_back( protocols::loops::Loop(start_res, end_res, 0, 0.0, false) );
			} else if ( end_res == start_res ) {
				loops.push_back( protocols::loops::Loop(start_res/*-1*/, end_res/*+1*/, 0, 0.0, false) );
			} else {  // end_res < start_res
				loops.push_back( protocols::loops::Loop(end_res/*-1*/, start_res/*+1*/, 0, 0.0, false) );
			}
		}
		start_res = simple_seg_list[nsegs][1].end()+1;
		end_res   = nres;
		//cutpt = (start_res+end_res)/2;  // set but never used ~Labonte
		if ( end_res > start_res ) {
			loops.push_back( protocols::loops::Loop(start_res, end_res, 0, 0.0, false) );
		} else if ( end_res == start_res ) {
			loops.push_back( protocols::loops::Loop(start_res/*-1*/, end_res, 0, 0.0, false) );
		}

		// split loops on cutpoints
		std::sort( cutpoints.begin(), cutpoints.end() );
		for ( int i=1; i<=(int)cutpoints.size(); ++i ) {
			for ( int j=1; j<=(int)loops.size(); ++j ) {
				if ( cutpoints[i] > loops[j].start() && cutpoints[i] < loops[j].stop() ) {
					// make 2 new loops
					protocols::loops::Loop l1(loops[j].start(), cutpoints[i]), l2(cutpoints[i]+1, loops[j].stop());
					loops[j] = l1;          // replace the original loop
					loops.push_back( l2 );  // add second half to end
				}
			}
		}
	}

	// now assemble compound segments
	rbsegs.clear();

	// create _locked_ segments, removing constituent simple segments from map as compound seg is created
	for ( core::Size i=1; i<=locked_segments.size(); ++i ) {
		utility::vector1< RBSegment > thisLockSeg;
		for ( core::Size j=1; j<=locked_segments[i].size(); ++j ) {
			int id = locked_segments[i][j];
			if ( simple_segments.find( id ) == simple_segments.end() ) {
				TR_seg << "[ERROR] Segment id " << id << " undefined or used in multiple locked compound segments" << std::endl;
				exit(1);
			}
			thisLockSeg.push_back( simple_segments[ id ] );
			simple_segments.erase( id );
		}
		RBSegment toAdd( thisLockSeg );
		toAdd.set_movement( locked_mvmt[i].first, locked_mvmt[i].second );
		rbsegs.push_back( toAdd );
	}

	// create _fixed_ segments, keeping constituent simple segments in map
	for ( core::Size i=1; i<=compound_segments.size(); ++i ) {
		utility::vector1< RBSegment > thisCompoundSeg;
		for ( core::Size j=1; j<=compound_segments[i].size(); ++j ) {
			int id = compound_segments[i][j];
			if ( simple_segments.find( id ) == simple_segments.end() ) {
				TR_seg << "[ERROR] Segment id " << id << " undefined or used in multiple locked compound segments" << std::endl;
				exit(1);
			}
			thisCompoundSeg.push_back( simple_segments[ id ] );
		}
		RBSegment toAdd( thisCompoundSeg );
		toAdd.set_movement( compound_mvmt[i].first, compound_mvmt[i].second );
		rbsegs.push_back( toAdd );
	}

	// place remaining simple structs in vector
	// ids no longer matter
	for ( std::map < int, RBSegment >::iterator i=simple_segments.begin(); i!=simple_segments.end(); ++i ) {
		rbsegs.push_back( i->second );
	}

	// sort by start residue
	std::sort( rbsegs.begin(), rbsegs.end(), RB_lt());
	std::sort( loops.v_begin(), loops.v_end(), protocols::loops::Loop_lt());

	// debug ... print in the data structure
	for ( core::Size i=1; i<=rbsegs.size(); ++i ) {
		if ( rbsegs[i].isSimple() ) {
			TR_seg << "[debug] SIMPLE rb   " << rbsegs[i][1].start() << "--" << rbsegs[i][1].end() << "   [type " << (int)rbsegs[i][1].type() << "] " << "       r = " << rbsegs[i].getSigAxisR()
				<< " t = " << rbsegs[i].getSigAxisT()
				<< " r_o = " << rbsegs[i].getSigOffAxisR()
				<< " t_o = " << rbsegs[i].getSigOffAxisT() << std::endl;
		} else {
			TR_seg << "[debug] COMPOUND rb       r = " << rbsegs[i].getSigAxisR()
				<< " t = " << rbsegs[i].getSigAxisT()
				<< " r_o = " << rbsegs[i].getSigOffAxisR()
				<< " t_o = " << rbsegs[i].getSigOffAxisT() << std::endl;
			for ( core::Size j=1; j<=rbsegs[i].nContinuousSegments(); ++j ) {
				TR_seg << rbsegs[i][j].start() << "--" << rbsegs[i][j].end() << "   [type " << (int)rbsegs[i][j].type() << "] ";
			}
			TR_seg << std::endl;
		}
	}
	for ( core::Size i=1; i<=loops.size(); ++i ) {
		TR_seg << "[debug] LOOP  " << loops[i].start() << "--" << loops[i].stop()
			<< "  " << loops[i].cut() << "  [skip " << loops[i].skip_rate() << "]" << std::endl;
	}
}


//////////////////////////////////////////////////////////
/// @brief Select a single RB segment to perturb + attached loops
/////////////////////////////////////////////////////////
void select_RBsegments(
	utility::vector1< RBSegment > const &rbsegs_in,
	protocols::loops::Loops             const &loops_in,
	utility::vector1< RBSegment >       &rbsegs_selected,
	protocols::loops::Loops                   &loops_selected
)
{
	rbsegs_selected.clear();
	loops_selected.clear();

	int nRBSegs = (int)rbsegs_in.size();
	if ( nRBSegs == 0 ) {
		return;
	}

	int k;
	if ( nRBSegs == 1 ) {
		k=1;
	} else {
		k = numeric::random::rg().random_range(1, nRBSegs);
	}

	rbsegs_selected.push_back( rbsegs_in[k] );

	// now find _all_ attached loops
	for ( core::Size i=1; i<=loops_in.size(); ++i ) {
		for ( core::Size j=1; j<=rbsegs_in[k].nContinuousSegments(); ++j ) {
			bool adjLoopN = loops_in[i].stop() == rbsegs_in[k][j].start()-1;
			bool adjLoopC = loops_in[i].start() == rbsegs_in[k][j].end()+1;
			if ( adjLoopN || adjLoopC ) {
				protocols::loops::Loop bounding_loop = loops_in[i];

				loops_selected.push_back( bounding_loop );
			}
		}
	}

	std::cerr << rbsegs_selected.size() << " rigid body segments input\n";
	std::cerr << loops_selected.size() << " loops\n";
} // void LoopRebuild::select_loops


RBSegment RBSegment::remap( core::id::SequenceMapping const &mapping ) const {
	utility::vector1 < RBSegment > remappedSimpleSegs;
	for ( core::Size i=1 ; i<=segments_.size(); ++i ) {
		int oldS = (int)segments_[i].start(),
			oldE = (int)segments_[i].end();
		int newS, newE;

		if ( mapping[ oldS ] != 0 ) {
			newS = mapping[ oldS ];
		} else {
			TR_seg << "[ ERROR ] mapping(" << segments_[i].start() << ") not found!" << std::endl;
			exit(1);
		}

		if ( mapping[ oldE ] != 0 ) {
			newE = mapping[ oldE ];
		} else {
			TR_seg << "[ ERROR ] mapping(" << segments_[i].end() << ") not found!" << std::endl;
			exit(1);
		}
		TR_seg << "[ DEBUG ] remapping(" << oldS << "," << oldE << ") to (" << newS << "," <<  newE << ")" << std::endl;

		remappedSimpleSegs.push_back( RBSegment( newS, newE, segments_[i].type() ) );
	}
	return RBSegment( remappedSimpleSegs );
}


RBSegment::RBSegment ( utility::vector1 < RBSegment > const &segs_in ) {
	for ( core::Size i=1; i<=segs_in.size(); ++i ) {
		if ( !segs_in[i].isSimple() ) {
			TR_seg << "[ ERROR ]  Error parsing RBSeg file: trying to nest compound segments." << std::endl;
			exit(1);
		} else {
			segments_.push_back( segs_in[i][1] );
		}
	}
	sigAxisR_ = sigAxisT_ = sigOffAxisR_ = sigOffAxisT_ = 0.0;
}


}
}
