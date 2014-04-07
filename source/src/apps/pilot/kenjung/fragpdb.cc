// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Ken Jung
/// @brief
//	Takes as input an index and offset and returns the corresponding frag in a pdb
//	OR a bin key and returns the pdbs of all the frags in the bin
//	should take only one loopsize

// libRosetta headers
//#include <basic/options/keys/in.OptionKeys.gen.hh>
//j#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <basic/Tracer.hh>
// for making poses from sequences
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/annotated_sequence.hh>

// C++ headers
#include <iostream>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <devel/init.hh>

#include <boost/lexical_cast.hpp>

using basic::T;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core;
using namespace protocols::loophash;
static basic::Tracer TR("main");


//forward declaration
void output_pdb( std::vector < BackboneSegment> & bs_vec_, core::Size loopsize, std::string out_path, std::string prefix );

int
main( int argc, char * argv [] )
{
    try {
	// initialize core
	devel::init(argc, argv);

  utility::vector1 < core::Size > loopsizes = option[lh::loopsizes]();
	utility::vector1 < core::Size > indexoffset = option[lh::fragpdb::indexoffset]();
	utility::vector1 < std::string > testkey = option[lh::fragpdb::bin]();
	utility::vector1 < boost::uint64_t > key;
	if ( testkey.size() != 0 ) {
			for (core::Size i = 1; i <= testkey.size(); i++) {
					key.push_back( boost::lexical_cast < boost::uint64_t > (testkey[i] ) );
			}
	}


	std::string out_path = option[lh::fragpdb::out_path];
	if (out_path != "")
					out_path = out_path + "/";

	// check if index/offset or key is null, no point in loading db if those are null
	if (indexoffset[1] == -1 && key.size() == 0 ) {
				 TR << "Need to set fragpdb:indexoffset or fragpdb:bin" << std::endl;
				 return 1;
	}
	if (indexoffset[1] != -1 && indexoffset.size() % 2 == 1 ) {
				TR << "Indexoffset must be paired" << std::endl;
				return 1;
	}


	if (loopsizes.size() > 1) {
				TR << "Too many loopsizes!" << std::endl;
				return 1;
	}
	core::Size loopsize = loopsizes[1];
	// load the whole library
	LoopHashLibraryOP loop_hash_library = new LoopHashLibrary( loopsizes, 1, 1);
	loop_hash_library->load_db();

	// if input is index and offset, go grab the angles
	std::vector < BackboneSegment > bs_vec_;
	const BackboneDB & bbdb_ = loop_hash_library->backbone_database();
	BackboneSegment backbone_;
	if (indexoffset[1] != -1 ){
					for ( core::Size i = 1; i <= indexoffset.size(); i=i+2 ) {
									bbdb_.get_backbone_segment( indexoffset[i], indexoffset[i+1], loopsize, backbone_ );
									bs_vec_.push_back(backbone_);
					}
					if( bs_vec_.size() == 0 ) {
								TR << "No frags matched" << std::endl;
					}
					else {
								output_pdb(bs_vec_, loopsize, out_path, "");
					}
	}
	else if (key.size() != 0 ){
					LoopHashMap &hashmap = loop_hash_library->gethash(loopsize);
					std::vector < core::Size > leap_index_list;
					for ( core::Size i = 1; i <= key.size(); i++ ) {
									hashmap.lookup_withkey( key[1], leap_index_list );
									for( std::vector < core::Size >::const_iterator itx = leap_index_list.begin();
											itx != leap_index_list.end(); ++itx ){
												core::Size bb_index = *itx;
												LeapIndex cp = hashmap.get_peptide( bb_index );
												bbdb_.get_backbone_segment( cp.index, cp.offset , loopsize , backbone_ );
												bs_vec_.push_back( backbone_ );
									}
									if( bs_vec_.size() == 0 ) {
												TR << "No frags matched" << std::endl;
									}
									else {
												output_pdb(bs_vec_, loopsize, out_path, utility::to_string(key[i]) + ".");
									}
					}

	}

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}

void output_pdb( std::vector < BackboneSegment> & bs_vec_, core::Size loopsize, std::string out_path, std::string prefix ) {
					core::pose::Pose frame;
					std::string sequence( loopsize+4, 'A' );
					core::pose::make_pose_from_sequence( frame, sequence,
						*( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ))  );
					for ( core::Size i = 0; i < bs_vec_.size(); i++) {
									for ( core::Size pos = 1; pos <= frame.total_residue(); pos++ ) {
										frame.set_phi( pos, -150 );
										frame.set_psi( pos, 150 );
										frame.set_omega( pos, 180 );
									}
									bs_vec_[i].apply_to_pose(frame, 3);
									frame.dump_pdb(out_path + prefix + utility::to_string(i) + ".pdb");
					}
}
