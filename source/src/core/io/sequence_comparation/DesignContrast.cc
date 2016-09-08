// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// (C) 199x-2008 University of Washington
// (C) 199x-2008 University of California Santa Cruz
// (C) 199x-2008 University of California San Francisco
// (C) 199x-2008 Johns Hopkins University
// (C) 199x-2008 University of North Carolina, Chapel Hill
// (C) 199x-2008 Vanderbilt University

/// @file   core/io/sequence_comparation/DesignContrast.cc
///
/// @brief
/// @author Yi Liu


// C++ headers
#include <fstream>
#include <iostream>
#include <cstdlib>

// Unit Headers
#include <basic/Tracer.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>


// Project Headers
#include <core/conformation/Residue.hh>
#include <core/io/sequence_comparation/DesignContrast.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <basic/options/option.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>


namespace core {
namespace io {
namespace sequence_comparation {

using namespace core;
using namespace basic;

using namespace core::pose;

using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace ObjexxFCL::format;

using basic::T;
using basic::Error;
using basic::Warning;

using core::Size;

using utility::vector1;

/// Tracer instance for this file
static THREAD_LOCAL basic::Tracer TR( "core.io.sequence_comparation.DesignContrast" );

DesignContrast::DesignContrast(DesignContrast const & dc) {
	list_file_names_ = dc.getListNames();
	pdb_file_names_ = dc.getPdbNames();
	nneighbs_ = dc.getNeighbors();
	secstructs_ = dc.getSecStruct();
	pdb_codes_ = dc.getPdbCodes();
}

/// @details go through residues in the pose, count the number of neighbors of each residue.
/// currently using tenA neighbor graph to # neighbors within 10 Angstroms
void DesignContrast::setNeighbors(pose::Pose & pose) {
	pose.update_residue_neighbors();
	nneighbs_.clear();
	for ( Size i=1; i<= pose.size(); ++i ) {
		nneighbs_.push_back(pose.energies().tenA_neighbor_graph().get_node(i)->num_neighbors_counting_self());
	}
}

vector1<int> & DesignContrast::getNeighbors(){
	return nneighbs_;
}

vector1<int> const & DesignContrast::getNeighbors() const {
	return nneighbs_;
}


/// @details Set secondary structure for residues in the pose. Have to call
/// set_ss_from_phipsi to set up the secondary structure infomation for pose.
void DesignContrast::setSecStruct(pose::Pose & pose){
	secstructs_.clear();
	set_ss_from_phipsi(pose); // This function is defined in pose/util.cc
	std::string secstruct;
	for ( Size i=1; i<= pose.size(); ++i ) {
		secstruct = pose.secstruct(i);
		secstructs_.push_back(secstruct);
	}
}

vector1<std::string> & DesignContrast::getSecStruct(){
	return secstructs_;
}
vector1<std::string> const & DesignContrast::getSecStruct() const {
	return secstructs_;
}

/// @details Copied the major parts of this function from Andrew. Read the lines
/// of the list files and store the pdb names in pdb_file_names.
void DesignContrast::setNames () {
	list_file_names_.clear();
	pdb_file_names_.clear();
	if ( option[ in::file::s ].active() ) {
		pdb_file_names_ = option[ in::file::s ]().vector();
	}
	if ( option[ in::file::l ].active() ) {
		list_file_names_ = option[ in::file::l ]().vector(); // make the copy of list file
	}
	for ( vector1< FileName >::iterator i = list_file_names_.begin(), i_end = list_file_names_.end(); i != i_end; ++i ) {
		std::string listname( i->name() );
		utility::io::izstream data( listname.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open file: " + listname + '\n' );
		}
		std::string line;
		while ( getline(data, line) ) {
			pdb_file_names_.push_back( line );
		}
		data.close();
	}
}

vector1<FileName> & DesignContrast::getPdbNames(){
	return pdb_file_names_;
}

vector1<FileName> const & DesignContrast::getPdbNames() const {
	return pdb_file_names_;
}

vector1<FileName> & DesignContrast::getListNames() {
	return list_file_names_;
}

vector1<FileName> const & DesignContrast::getListNames() const {
	return list_file_names_;
}

void DesignContrast::setPdbCodes(){
	std::string code;
	for ( Size i=1; i<=pdb_file_names_.size(); ++i ) {
		code = std::string(pdb_file_names_[i]).substr(std::string(pdb_file_names_[i]).size()-8,4);
		TR << "in DesignContrast check the code: " << code <<std::endl;
		pdb_codes_.push_back(code);
	}
}

vector1<std::string> & DesignContrast::getPdbCodes(){
	return pdb_codes_;
}

vector1<std::string> const & DesignContrast::getPdbCodes() const {
	return pdb_codes_;
}

void DesignContrast::output_sqc_file (
	pose::Pose & native_pose,
	pose::Pose & decoy_pose,
	std::string const & single_code,
	std::ofstream & sqc
){

	TR << "in the start of output_sqc_file"<< std::endl;
	setNeighbors(native_pose);
	TR << "after setNeighbors "<< std::endl;
	setSecStruct(decoy_pose);
	TR << "after setSecStruct "<< std::endl;

	vector1<std::string> native_res_name;
	vector1<std::string> decoy_res_name;
	vector1<char> chain_ids;
	vector1<int> pdb_res_num;

	// this block is commented out due to deprecation of old PDBInfo interface
	// chain ids
	//    chain_ids = decoy_pose.pdb_chains();
	//  TR << "after get chains"<< std::endl;
	// sequence position for each residue.
	//    pdb_res_num = decoy_pose.pdb_numbering();
	//  TR << "after get pdb number "<< std::endl;

	for ( Size i=1; i<=native_pose.size(); ++i ) {
		native_res_name.push_back(native_pose.residue(i).name3());
	}
	for ( Size j=1 ; j<=decoy_pose.size(); ++j ) {
		decoy_res_name.push_back(decoy_pose.residue(j).name3());
	}
	TR << "total residues: " << decoy_pose.size() << std::endl;
	for ( Size n=1; n<=decoy_pose.size(); ++n ) {
		sqc
			<< single_code << "  "
			<< decoy_pose.pdb_info()->chain(n) << " "
			<< I( 4, decoy_pose.pdb_info()->number(n) ) << ' '
			<< A( 3, native_res_name[n] ) << ' '
			<< A( 3, decoy_res_name[n] ) << ' '
			<< I( 5, nneighbs_[n]) << ' '
			<< A( 2, secstructs_[n] ) << ' '
			<< '\n';
	}
	TR << "end of sqc output" << std::endl;
}

void DesignContrast::clear() {
	list_file_names_.clear();
	pdb_file_names_.clear();
	nneighbs_.clear();
	secstructs_.clear();
	pdb_codes_.clear();
}

} // namespace sequence_comparation
} // namespace io
} // namespace core
