// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SequenceCoupling.hh
/// @brief class definition for a sequence coupling profile that
//   represents a  probability distribution over the entire protein.
/// @author Hetu Kamisetty

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/SequenceCoupling.hh>

#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <utility/pointer/owning_ptr.hh>

#include <core/chemical/AA.hh>

#include <iostream>
#include <string>

#include <utility/vector1.hh>

namespace core {
namespace sequence {

static basic::Tracer tr( "core.sequence.SequenceCoupling" );

/*
void SequenceCoupling::profile(utility::vector1< utility::vector1< Real > > const & new_profile){
SequenceProfile::profile(new_profile);
}

utility::vector1< utility::vector1< Real > > const & SequenceCoupling::profile() const{
SequenceProfile::profile();
}
*/
void SequenceCoupling::delete_position( core::Size pos ){
	using utility::vector1;
	vector1< vector1< vector1< Real > > > newEdgePots;
	for ( core::Size i=1; i<=edgeList_.size(); i++ ) {
		vector1 < Real > vertices = edgeList_[i];
		if ( !(vertices[1]==pos || vertices[2]==pos) ) {
			newEdgePots.push_back(edgePots_[i]);
		}
		if ( vertices[1]>pos ) {
			vertices[1]--;
		}
		if ( vertices[2]>pos ) {
			vertices[2]--;
		}
	}
	edgePots(newEdgePots);
	numVerts_--;
	numEdges_ = edgePots_.size();

	SequenceProfile::delete_position(pos);
}
void SequenceCoupling::insert_char(core:: Size pos, char new_char){
	SequenceProfile::insert_char(pos,new_char);
	numVerts_++;
	for ( core::Size i=1; i<=edgeList_.size(); i++ ) {
		utility::vector1 < Real > vertices = edgeList_[i];
		if ( vertices[1]>=pos ) {
			vertices[1]++;
		}
		if ( vertices[2]>=pos ) {
			vertices[2]++;
		}
	}
}
void SequenceCoupling::read_from_file(
	utility::file::FileName const & fn
) {
	read_from_file(fn,1.0);
}
void SequenceCoupling::read_from_file(
	utility::file::FileName const & fn,
	core::Real temp)
{

	static utility::vector1< core::chemical::AA > order;
	order.resize( 20 );
	order[ 1] = core::chemical::aa_from_oneletter_code( 'A');
	order[ 2] = core::chemical::aa_from_oneletter_code( 'R');
	order[ 3] = core::chemical::aa_from_oneletter_code( 'N');
	order[ 4] = core::chemical::aa_from_oneletter_code( 'D');
	order[ 5] = core::chemical::aa_from_oneletter_code( 'C');
	order[ 6] = core::chemical::aa_from_oneletter_code( 'Q');
	order[ 7] = core::chemical::aa_from_oneletter_code( 'E');
	order[ 8] = core::chemical::aa_from_oneletter_code( 'G');
	order[ 9] = core::chemical::aa_from_oneletter_code( 'H');
	order[10] = core::chemical::aa_from_oneletter_code( 'I');
	order[11] = core::chemical::aa_from_oneletter_code( 'L');
	order[12] = core::chemical::aa_from_oneletter_code( 'K');
	order[13] = core::chemical::aa_from_oneletter_code( 'M');
	order[14] = core::chemical::aa_from_oneletter_code( 'F');
	order[15] = core::chemical::aa_from_oneletter_code( 'P');
	order[16] = core::chemical::aa_from_oneletter_code( 'S');
	order[17] = core::chemical::aa_from_oneletter_code( 'T');
	order[18] = core::chemical::aa_from_oneletter_code( 'W');
	order[19] = core::chemical::aa_from_oneletter_code( 'Y');
	order[20] = core::chemical::aa_from_oneletter_code( 'V');
	//order[21] = core::chemical::aa_from_oneletter_code( '-');
	utility::io::izstream input( fn );

	temp_ = temp;
	if ( !input ) {
		std::string msg(
			"ERROR: Unable to open file " +
			static_cast< std::string > (fn) +
			"!"
		);
		utility_exit_with_message( msg );
	}
	std::string line;

	tr.Debug << "reading from " << fn << std::endl;

	/*format*
	* GREMLIN - v1
	* <num_verts> <num_edges>
	* vert_pot1_A vert_pot1_R ...
	* vert_pot2_A vert_pot2_R ...
	* ..(num_vert_lines)
	* edgei edgej edge_pot1_AA edge_pot1_AR ...
	* edgei' edgej'  edge_pot2_AA ..
	* ..(num_edge lines)
	*/
	//one header line. ignored for now.
	getline( input, line );
	// initialize headers
	getline( input, line );
	std::istringstream line_stream( line );
	line_stream >> numVerts_;
	if ( line_stream.fail() ) {
		utility_exit_with_message( "ERROR: failed while reading numVerts!" );
	}
	line_stream >> numEdges_;
	if ( line_stream.fail() ) {
		utility_exit_with_message( "ERROR: failed while reading numEdges!" );
	}
	// AMW: cppcheck correctly notices that aa_seq is actually never used anywhere
	// but here, where it is gradually built up
	// if you would like to do something with it, uncomment here and on the line below
	//std::string aa_seq;
	core::Size lineCount=0;
	utility::vector1< utility::vector1< Real > > vertexPots;
	while ( lineCount<numVerts_ && getline(input,line) ) {
		utility::vector1< core::Real > vertPot;
		vertPot.resize( order.size() );
		std::istringstream ls( line );
		core::Real aa_pot;
		core::Real min_pot = 1e9;
		core::Size minPotAA=0;
		Size count=1;
		while ( ls>>aa_pot ) {
			vertPot[order[count]]=aa_pot;
			if ( min_pot>aa_pot ) {
				min_pot = aa_pot;
				minPotAA=count;
			}
			count++;
		}
		if ( minPotAA>order.size() ) {
			std::cout << " minpotAA" << minPotAA << "count " << count << std::endl;
			utility_exit_with_message("minPotAA was of wrong size. quitting");
		}
		//aa_seq +=core::chemical::oneletter_code_from_aa(order[minPotAA]);//dummy sequence set to most favored aa according to vertexPot
		vertexPots.push_back(vertPot);
		lineCount++;
	}
	//done vertexpots
	profile(vertexPots);
	core::Size edgeCount=0;
	while ( edgeCount<numEdges_ && getline(input,line) ) {
		utility::vector1< utility::vector1<core::Real> > edgePot;
		edgePot.resize(order.size());
		utility::vector1<core::Size>  edge;
		std::istringstream ls( line );
		core::Size verti, vertj;

		//read vertices
		ls>>verti;
		ls>>vertj;
		edge.push_back(verti);
		edge.push_back(vertj);
		edgeList_.push_back(edge);

		//read edgepots
		for ( core::Size i=1; i<=order.size(); i++ ) {
			utility::vector1 < core::Real > edgePotLine;
			edgePotLine.resize( order.size() );
			for ( core::Size j=1; j<=order.size(); j++ ) {
				core::Real edgeWt;
				ls>>edgeWt;
				//edgePotLine.push_back(edgeWt);
				edgePotLine[order[j]] = edgeWt;
			}
			//edgePot.push_back(edgePotLine);
			edgePot[order[i]] = edgePotLine;
		}
		edgePots_.push_back(edgePot);
		edgeCount++;
	}
}
core::Size SequenceCoupling::findEdgeId(Size vert1, Size vert2) const{
	Size i=0;
	for ( i=1; i<=edgeList_.size(); i++ ) {
		utility::vector1 < Real > vertices = edgeList_[i];
		if ( vertices[1]==vert1 && vertices[2] ==vert2 ) {
			return i;
		}
	}
	return 0;
}
utility::vector1< utility::vector1 < Real > > const & SequenceCoupling::edgePotBetween(Size vert1, Size vert2) const {//direction important. if (i,j) edge exists, will not return if (j,i) asked
	runtime_assert( vert1<= profile().size() && vert2 <= profile().size());
	Size edgeId = findEdgeId(vert1,vert2);
	return edgePotBetween(edgeId);
}
utility::vector1< utility::vector1 < Real > > const & SequenceCoupling::edgePotBetween(Size edgeId) const {//direction important. if (i,j) edge exists, will not return if (j,i) asked
	runtime_assert(edgeId>0);
	return edgePots_[edgeId];
}
}//ns sequence
}//ns core

