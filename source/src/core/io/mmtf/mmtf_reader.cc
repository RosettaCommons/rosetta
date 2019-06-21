// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/mmtf/mmtf_reader.hh
/// @author Daniel Farrell (danpf@uw.edu)


// Unit headers
#include <core/io/mmtf/mmtf_reader.hh>
#include <core/io/StructFileReaderOptions.hh>

// Package headers
#include <core/io/pdb/pdb_reader.hh>  // TODO: Pull out pseudo-duplicated code and move to sfr_storage.cc.

// When you move PDBReader and PoseUnbuilder, take these.
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/ResidueConnection.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/string_constants.hh>
#include <utility/io/izstream.hh>
#include <numeric/random/random.hh>

#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/Field.hh>
#include <core/io/HeaderInformation.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/io/StructFileRep.hh>

#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>

#include <core/io/Remarks.hh>

// Project headers
#include <core/types.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/string_constants.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_map.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <mmtf.hpp>

// C++ headers
#include <cstdlib>
#include <cstdio>
#include <algorithm>


static basic::Tracer TR( "core.io.mmtf.mmtf_reader" );

using basic::Error;
using basic::Warning;


core::io::StructFileRepOP
core::io::mmtf::create_sfr_from_mmtf_filename( std::string mmtf_filename,
	StructFileReaderOptions const & options ) {

	if ( options.read_pdb_header() ) {
		TR.Warning << "read_pdb_header set, but not used in reading mmtf_file" << std::endl;
		// TODO
	}

	::mmtf::StructureData sd;
	::mmtf::decodeFromFile(sd, mmtf_filename);

	core::io::StructFileRepOP sfr( new core::io::StructFileRep );

	std::map<char, core::io::ChainAtoms> atom_chain_map;
	std::vector< char > chain_list; // preserve order
	std::map<char, core::Size> chain_to_idx;

	core::Size modelIndex = 0;
	core::Size chainIndex = 0;
	core::Size groupIndex = 0;
	core::Size atomIndex = 0;
	core::Size atomSerial = 0;
	utility::vector1<char> known_chainIDs;

	// Ignore anything past first MODEL
	// for ( int i = 0; i < sd.numModels; i++, modelIndex++ ) {
	for ( int i = 0; i < 1; i++, modelIndex++ ) {
		// traverse chains
		for ( int j = 0; j < sd.chainsPerModel[modelIndex]; j++, chainIndex++ ) {
			// traverse groups
			for ( int k = 0; k < sd.groupsPerChain[chainIndex]; k++, groupIndex++ ) {
				const ::mmtf::GroupType& group =
					sd.groupList[sd.groupTypeList[groupIndex]];
				int groupAtomCount = group.atomNameList.size();
				for ( int l = 0; l < groupAtomCount; l++, atomIndex++ ) {
					++atomSerial;
					core::io::AtomInformation ai;

					if ( !sd.entityList.empty() ) {
						ai.isHet = ::mmtf::is_hetatm(chainIndex, sd.entityList, group);
					} else ai.isHet = ::mmtf::is_hetatm(group.chemCompType.c_str());

					if ( !::mmtf::isDefaultValue(sd.atomIdList) ) {
						ai.serial = sd.atomIdList[atomIndex];
					} else ai.serial = atomSerial;
					ai.name = group.atomNameList[l]; // can be up to 5 characters
					if ( !::mmtf::isDefaultValue(sd.altLocList) ) {
						ai.altLoc = sd.altLocList[atomIndex]==0x00 ? ' ' : sd.altLocList[atomIndex];
					} else ai.altLoc = ' ';
					ai.resName = group.groupName; // can be up to 5 characters

					// TODO we default to chainNameList because that's consistent with pdbs.
					// However, they are not required. so we always take chainName if it's available, but
					// otherwise use chainId
					if ( !::mmtf::isDefaultValue(sd.chainNameList) ) {
						ai.chainID = sd.chainNameList[chainIndex].at(0); // TODO can be up to 4 characters
					} else ai.chainID = sd.chainIdList[chainIndex].at(0); // TODO can be up to 4 characters

					if ( options.new_chain_order() ) {
						auto chainID_find = std::find(known_chainIDs.begin(), known_chainIDs.end(), ai.chainID);
						core::Size distance = 0;
						if ( chainID_find == known_chainIDs.end() ) {
							known_chainIDs.push_back(ai.chainID);
							distance = known_chainIDs.size();
						} else distance = std::distance(known_chainIDs.begin(), chainID_find);
						ai.chainID = utility::ALPHANUMERICS.at(distance);
					}

					ai.resSeq = sd.groupIdList[groupIndex];

					if ( !::mmtf::isDefaultValue(sd.insCodeList) ) {
						ai.iCode = sd.insCodeList[groupIndex]==0x00 ? ' ' : sd.insCodeList[groupIndex];
					} else ai.iCode = ' ';

					if ( !::mmtf::isDefaultValue(sd.occupancyList) ) {
						ai.occupancy = sd.occupancyList[atomIndex];
					} else ai.occupancy = 1.0;

					ai.x = sd.xCoordList[atomIndex];
					ai.y = sd.yCoordList[atomIndex];
					ai.z = sd.zCoordList[atomIndex];

					if ( !::mmtf::isDefaultValue(sd.bFactorList) ) {
						ai.temperature = sd.bFactorList[atomIndex];
					} else ai.temperature = 1.0;
					ai.segmentID = "    "; // what should this be?
					ai.element = group.elementList[l];
					ai.terCount = 0;  // what should this be? maybe each chain?
					atom_chain_map[ai.chainID].push_back(ai);
					if ( std::find( chain_list.begin(), chain_list.end(), ai.chainID ) == chain_list.end() ) {
						chain_list.push_back( ai.chainID );
					}
				}
			}
		}
	}

	for ( char i : chain_list ) { // std::vector
		sfr->chains().push_back( atom_chain_map.find( i )->second );
	}
	return sfr;
}
