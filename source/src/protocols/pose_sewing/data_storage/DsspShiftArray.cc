// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/data_storage/DsspShiftArray.cc
/// @brief A container for an array of secondary structure boundaries
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/pose_sewing/data_storage/DsspShiftArray.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.pose_sewing.data_storage.DsspShiftArray" );


namespace protocols {
namespace pose_sewing {
namespace data_storage {

DsspShiftArray::DsspShiftArray():
	utility::VirtualBase()
{
	shift_array_.clear();
}

DsspShiftArray::DsspShiftArray( core::pose::Pose const & pose):
	utility::VirtualBase()
{
	shift_array_.clear();
	this->populate(pose);
}

DsspShiftArray::~DsspShiftArray(){}

DsspShiftArray::DsspShiftArray( DsspShiftArray const & ) = default;

void
DsspShiftArray::populate(core::pose::Pose const & pose){
	core::scoring::dssp::Dssp dssp(pose);
	char current_dssp = dssp.get_dssp_secstruct(1);

	if ( current_dssp == 'H' || current_dssp == 'G' || current_dssp == 'I' ) {
		current_dssp = 'H';
	} else if ( current_dssp == 'B' || current_dssp == 'E' ) {
		current_dssp = 'E';
	} else {
		current_dssp = 'L';
	}


	shift_array_.clear();
	shift_array_.push_back(1);
	for ( core::Size current_residue = 1; current_residue <= pose.size(); current_residue++ ) {
		char new_dssp = dssp.get_dssp_secstruct(current_residue);

		if ( new_dssp == 'H' || new_dssp == 'G' || new_dssp == 'I' ) {
			new_dssp = 'H';
		} else if ( new_dssp == 'B' || new_dssp == 'E' ) {
			new_dssp = 'E';
		} else {
			new_dssp = 'L';
		}

		if ( new_dssp!=current_dssp ) {
			current_dssp = new_dssp;
			shift_array_.push_back(current_residue);
		}
	}
	shift_array_.push_back(pose.size());
}
core::Size
DsspShiftArray::get_distance_to_nth_shift(core::Size N_res, core::Size shift_gap){
	//loop through residues in shift array. after interating 'shift_gap' number of times, return the distance
	core::Size seen_new_ss_count = 0;
	for ( utility::vector1<core::Size>::iterator it = shift_array_.begin(); it!=shift_array_.end(); ++it ) {
		if ( *it > N_res ) {
			++seen_new_ss_count;
			if ( seen_new_ss_count == shift_gap + 1 ) {
				return *it-N_res;
			}
		}
	}
	return 0;
}

DsspShiftArrayOP
DsspShiftArray::clone() const {
	return DsspShiftArrayOP( new DsspShiftArray( *this ) );
}


} //protocols
} //pose_sewing
} //data_storage






