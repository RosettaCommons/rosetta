// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/SidechainContactDistCutoff.cc
/// @brief  Defines sidechain contact distance cutoffs.
/// @author David E. Kim (dekim@u.washington.edu)

// Unit headers
#include <protocols/frag_picker/SidechainContactDistCutoff.hh>


// Package headers

// Project headers
#include <core/types.hh>
#include <core/chemical/AA.hh>

#include <basic/database/open.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/io/izstream.hh>

// C/C++ headers
#include <iostream>
#include <string>


namespace protocols {
namespace frag_picker {

// @brief Auto-generated virtual destructor
SidechainContactDistCutoff::~SidechainContactDistCutoff() = default;

using namespace core;

static basic::Tracer TR( "protocols.frag_picker.SidechainContactDistCutoff" );

SidechainContactDistCutoff::SidechainContactDistCutoff() {
	SidechainContactDistCutoff(1.0);
}

SidechainContactDistCutoff::SidechainContactDistCutoff(core::Real scale_factor) {
	scale_factor_ = scale_factor;
	initialize();
}

void SidechainContactDistCutoff::initialize() {

	// set aa mapping
	// GLY   ALA   SER   CYS   VAL   THR   ILE   PRO   MET   ASP   ASN   LEU   LYS   GLU   GLN   ARG   HIS   PHE   TYR   TRP
	aa_map_.resize(20);
	aa_map_[1] = 'G';
	aa_map_[2] = 'A';
	aa_map_[3] = 'S';
	aa_map_[4] = 'C';
	aa_map_[5] = 'V';
	aa_map_[6] = 'T';
	aa_map_[7] = 'I';
	aa_map_[8] = 'P';
	aa_map_[9] = 'M';
	aa_map_[10] = 'D';
	aa_map_[11] = 'N';
	aa_map_[12] = 'L';
	aa_map_[13] = 'K';
	aa_map_[14] = 'E';
	aa_map_[15] = 'Q';
	aa_map_[16] = 'R';
	aa_map_[17] = 'H';
	aa_map_[18] = 'F';
	aa_map_[19] = 'Y';
	aa_map_[20] = 'W';

	// set aa to index mapping
	for ( core::Size i=1; i<=aa_map_.size(); ++i ) {
		aa_to_index_map_[aa_map_[i]] = i;
	}

	utility::vector1<utility::vector1<core::Real> >  mean_dist;
	utility::vector1<utility::vector1<core::Real> >  stddev_dist;

	// read distance tables
	// first table is the mean and the second is the standard deviation
	utility::io::izstream data(basic::database::full_name("sampling/sidechain_contact.txt"));
	TR.Info << "read sidechain contact data from sidechain_contact.txt" << std::endl;
	if ( !data ) {
		utility_exit_with_message("[ERROR] Unable to open file: sidechain_contact.txt");
	}

	std::string line;
	core::Size tablecnt = 0;
	while ( getline(data, line) ) {
		if ( line.size() > 3 && line.substr(0,3) == "GLY" ) { // tables start with GLY
			tablecnt++;
			core::Size rowcnt = 1;
			utility::vector1<core::Real> col(20);
			std::istringstream line_stream( line );
			std::string col1;
			line_stream >> col1;
			if ( core::chemical::oneletter_code_from_aa(core::chemical::aa_from_name(col1)) != aa_map_[rowcnt] ) {
				utility_exit_with_message( "Error in sidechain name mapping in sidechain_contact.txt!" );
			}
			for ( core::Size i=1; i<=col.size(); ++i ) line_stream >> col[i];
			if ( line_stream.fail() ) {
				utility_exit_with_message( "Error reading in SidechainContactDistCutoff()!" );
			}
			if ( tablecnt == 1 ) {
				mean_dist.push_back(col);
			} else if ( tablecnt == 2 ) {
				stddev_dist.push_back(col);
			}
			for ( core::Size j=1; j<=19; ++j ) { // read the remaining table rows
				getline(data, line);
				rowcnt++;
				std::istringstream line_stream_next(line);
				line_stream_next >> col1;
				if ( core::chemical::oneletter_code_from_aa(core::chemical::aa_from_name(col1)) != aa_map_[rowcnt] ) {
					utility_exit_with_message( "Error in sidechain name mapping in sidechain_contact.txt!" );
				}
				for ( core::Size i=1; i<=col.size(); ++i ) line_stream_next >> col[i];
				if ( line_stream_next.fail() ) {
					utility_exit_with_message( "Error reading in SidechainContactDistCutoff()!" );
				}
				if ( tablecnt == 1 ) {
					mean_dist.push_back(col);
				} else if ( tablecnt == 2 ) {
					stddev_dist.push_back(col);
				}
			}
		}
	}
	data.close();

	// calculate cutoffs
	// from I-TASSER paper
	cutoff_.resize(20);
	cutoff_squared_.resize(20);
	for ( core::Size i=1; i<=mean_dist.size(); ++i ) {
		for ( core::Size j=1; j<=mean_dist[i].size(); ++j ) {
			cutoff_[i].push_back( scale_factor_*(mean_dist[i][j] + 2.5*stddev_dist[i][j]) );
			cutoff_squared_[i].push_back( cutoff_[i][j]*cutoff_[i][j] );
			TR.Debug << "cutoff_squared: " << i << " " << j << " " << cutoff_squared_[i][j] << std::endl;
		}
	}
}

core::Real SidechainContactDistCutoff::get_cutoff( char aa_i, char aa_j ) {
	return cutoff_[aa_to_index_map_[aa_i]][aa_to_index_map_[aa_j]];
}

core::Real SidechainContactDistCutoff::get_cutoff_squared( char aa_i, char aa_j ) {
	return cutoff_squared_[aa_to_index_map_[aa_i]][aa_to_index_map_[aa_j]];
}

core::Real SidechainContactDistCutoff::scale_factor() {
	return scale_factor_;
}

} // frag_picker
} // protocols


