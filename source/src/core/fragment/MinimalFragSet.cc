// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragments/FragSet.cc
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @author James Thompson (tex@u.washington.edu)
/// @author David E. Kim (dekim@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007


// Unit Headers
#include <core/fragment/MinimalFragSet.hh>

// Package Headers
#include <core/fragment/MinimalFragSetIterator_.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/BBTorsionSRFD.fwd.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FragData.hh>

// Project Headers
#include <core/types.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/io/izstream.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>
#include <ostream>

#include <core/fragment/FrameIterator.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// C/C++ headers
#include <string>
#include <map>

#ifdef WIN32
#include <iterator>
#endif

namespace core {
namespace fragment {

using namespace kinematics;

static THREAD_LOCAL basic::Tracer tr( "core.fragments" );
// preliminary reader method --- reads classic rosetta++ frag files


MinimalFragSet::MinimalFragSet() {}
MinimalFragSet::~MinimalFragSet() {}

FragSetOP MinimalFragSet::clone() const {
	return FragSetOP( new MinimalFragSet( *this ) );
}
FragSetOP MinimalFragSet::empty_clone() const {
	return FragSetOP( new MinimalFragSet() );
}

/// @brief get fragments that start somewhere between start and end
Size MinimalFragSet::region(
	MoveMap const&,
	core::Size start,
	core::Size end, //not used
	core::Size, //min_overlap not used
	core::Size, //min_length not used
	FrameList &frame_list
) const {
	Size count( 0 );
	for ( Size pos=start; pos<=end; pos++ ) {
		count += frames( pos, frame_list );
	}
	return count;
}


/// @brief Accessor for the Frame at the specified insertion position. Returns false if
/// there is no frame at the specified position.
Size MinimalFragSet::frames( Size pos, FrameList &out_frames ) const
{
	FrameMap::const_iterator it = frames_.find(pos);
	if ( it == frames_.end() ) return 0;
	if ( it->second.begin() != it->second.end() ) {
		copy( it->second.begin(), it->second.end(), back_inserter( out_frames ) ); // should append frames
		return it->second.size();
	}

	return 0;
}

void MinimalFragSet::read_fragment_file( std::string filename, Size top25, Size ncopies ) {
	using std::cerr;
	using std::endl;
	using std::istringstream;
	using std::string;

	utility::io::izstream data( filename );
	if ( !data.good() ) {
		cerr << "Open failed for file: " << data.filename() << endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}
	// read torsions only vall
	string line;
	utility::vector1<torsions> vall_torsions;
	Size insertion_pos = 1;
	FrameOP frame;
	std::map<std::pair<Size,Size>, Size> frame_counts;
	Size n_frags( 0 );
	while ( getline( data, line ) ) {
		if ( line.substr(0,1) == "#" ) {
			std::string vallname;
			Size vall_start_line;
			Size vall_end_line;
			Size vall_last_residue_key;
			istringstream in( line );
			in >> vallname; // #
			in >> vallname; // index
			in >> vall_start_line;
			in >> vall_end_line;
			in >> vall_last_residue_key;
			in >> vallname;
			Size totalsize = vall_last_residue_key+vall_end_line;
			if ( totalsize > vall_torsions.size() ) {
				vall_torsions.resize(totalsize);
			}
			tr.Info << "Reading Vall: " << vallname << ".torsions";
			utility::io::izstream valldata( basic::database::full_name("sampling/" + vallname + ".torsions") );
			if ( !valldata.good() ) {
				cerr << "Open failed for vall torsions database file: " << valldata.filename() << endl;
				utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
			}
			Size valllinecnt = 0;
			string vallline;
			while ( getline( valldata, vallline ) ) {
				if ( vallline.substr(0,1) == "#" ) continue;
				valllinecnt++;
				if ( valllinecnt >= vall_start_line && valllinecnt <= vall_end_line ) {
					// Chunk residue key is the last vall residue key + the vall line number for the residue
					Size residueindex = vall_last_residue_key+valllinecnt; // See VallProvider for residue index
					torsions t;
					istringstream vin( vallline );
					vin >> t.phi >> t.psi >> t.omega;
					vall_torsions[residueindex] = t;
				}
			}
		} else {
			// read in fragments

			// skip blank lines
			if ( line == "" || line == " " ) {
				insertion_pos++;
				continue;
			}
			istringstream in( line );
			Size vall_residue_key;
			Size fraglen;
			FragDataOP current_fragment( NULL );
			current_fragment = FragDataOP( new FragData );
			in >> vall_residue_key >> fraglen;
			if ( vall_residue_key == 0 ) { // read torsions from fragment file
				for ( Size i = 1; i <= fraglen; i++ ) {
					getline( data, line );
					istringstream fin( line );
					Real phi;
					Real psi;
					Real omega;
					fin >> phi >> psi >> omega;
					BBTorsionSRFDOP res( new BBTorsionSRFD(3, 'L', 'G') );
					res->set_torsion(1, phi);
					res->set_torsion(2, psi);
					res->set_torsion(3, omega);
					current_fragment->add_residue(res);
				}
			} else { // read torsions from torsion only vall
				Size last_residue_key = vall_residue_key + fraglen - 1;
				for ( Size i = vall_residue_key; i <= last_residue_key; i++ ) {
					// set ss 'L' and aa 'G' just as placeholders to prevent run time error in
					// GunnCost moves (fragment_as_pose method)
					BBTorsionSRFDOP res( new BBTorsionSRFD(3, 'L', 'G') );
					res->set_torsion(1, vall_torsions[i].phi);
					res->set_torsion(2, vall_torsions[i].psi);
					res->set_torsion(3, vall_torsions[i].omega);
					current_fragment->add_residue(res);
				}
			}
			current_fragment->set_valid();
			std::pair<Size,Size> p(insertion_pos,fraglen);
			if ( !top25 || frame_counts[p] < top25*ncopies ) {
				FrameOP frame = FrameOP( new Frame( insertion_pos ) );
				if ( !frame->add_fragment( current_fragment ) ) {
					tr.Fatal << "Incompatible Fragment in file: " << data.filename() << endl;
					utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
				} else {
					frame_counts[p]++;
					for ( Size i = 2; i <= ncopies; i++ ) {
						frame->add_fragment( current_fragment );
						frame_counts[p]++;
					}
				}
				if ( frame && frame->is_valid() ) {
					add( frame );
				}
				n_frags = std::max( n_frags, frame_counts[p] );
			}
		}
	}

	tr.Info << "finished reading top " << n_frags << " "
		<< max_frag_length() << "mer fragments from file " << data.filename()
		<< endl;
}

ConstFrameIterator MinimalFragSet::begin() const {
	return ConstFrameIterator( FrameIteratorWorker_OP( new MinimalFragSetIterator_( frames_.begin(), frames_.end() ) ) );
}

ConstFrameIterator MinimalFragSet::end() const {
	return ConstFrameIterator( FrameIteratorWorker_OP( new MinimalFragSetIterator_( frames_.end(), frames_.end() ) ) );
}

FrameIterator MinimalFragSet::nonconst_begin() {
	return FrameIterator( FrameIteratorWorker_OP( new MinimalFragSetIterator_( frames_.begin(), frames_.end() ) ) );
}

FrameIterator MinimalFragSet::nonconst_end() {
	return FrameIterator( FrameIteratorWorker_OP( new MinimalFragSetIterator_( frames_.end(), frames_.end() ) ) );
}

bool MinimalFragSet::empty() const {
	return frames_.empty();
}


void MinimalFragSet::add_( FrameOP aframe )
{
	Size seqpos( aframe->start() );
	frames_[ seqpos ].push_back( aframe );
}

}//fragment
}// core
