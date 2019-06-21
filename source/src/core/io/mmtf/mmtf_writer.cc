// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/mmtf/mmtf_writer.cc
/// @brief Functions for MMCIF writing.
/// @author Daniel Farrell (danpf@uw.edu)


#include <core/io/StructFileRepOptions.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>

// When you move PDBReader and PoseUnbuilder, take these.
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/ResidueConnection.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <numeric/random/random.hh>

#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/Field.hh>
#include <core/io/HeaderInformation.hh>
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
#include <utility/vector1.hh>
#include <utility/tools/make_map.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <msgpack.hpp>
#include <mmtf.hpp>

// C++ headers
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <tuple>


#include <basic/Tracer.hh>

#include <core/io/mmtf/mmtf_writer.hh>

static basic::Tracer TR( "core.io.mmtf.mmtf_writer" );

using basic::Error;
using basic::Warning;
using utility::to_string;

namespace core {
namespace io {
namespace mmtf {

bool
dump_mmtf(
	core::pose::Pose const & pose,
	std::string const & file_name,
	core::io::StructFileRepOptionsOP options )
{
	set_mmtf_default_options(*options);
	core::io::pose_to_sfr::PoseToStructFileRepConverter converter = core::io::pose_to_sfr::PoseToStructFileRepConverter( *options );

	converter.init_from_pose( pose );
	core::io::StructFileRepOP sfr =  converter.sfr();

	utility::io::ozstream file(file_name.c_str(), std::ios::out | std::ios::binary);
	if ( !file ) {
		TR.Error << "mmtf_writer:dump_mmtf Unable to open file:" << file_name << " for writing!!!" << std::endl;
		return false;
	}
	dump_mmtf( file, sfr, *options);
	file.close();
	return true;
}


/// @brief Dump an MMTF from a pose to an ostream.
core::io::StructFileRepOP
dump_mmtf(
	core::pose::Pose const & pose,
	std::ostream & out,
	core::io::StructFileRepOptionsCOP options)
{

	core::io::pose_to_sfr::PoseToStructFileRepConverter converter = core::io::pose_to_sfr::PoseToStructFileRepConverter( *options );

	converter.init_from_pose( pose );
	core::io::StructFileRepOP sfr =  converter.sfr();


	dump_mmtf( out, sfr, *options);
	return sfr;
}

/// @brief Dump an MMTF from a pose, optionally extracting extra info.
core::io::StructFileRepOP
dump_mmtf(
	core::pose::Pose const & pose,
	std::string const & jd2_job_data,
	utility::io::ozstream & out)

{
	core::io::StructFileRepOptionsOP options =  core::io::StructFileRepOptionsOP( new core::io::StructFileRepOptions );
	core::io::pose_to_sfr::PoseToStructFileRepConverter converter = core::io::pose_to_sfr::PoseToStructFileRepConverter();

	converter.init_from_pose( pose );
	core::io::StructFileRepOP sfr =  converter.sfr();
	sfr->append_to_additional_string_output( jd2_job_data );

	dump_mmtf( out, sfr, *options);
	return sfr;
}


bool
dump_mmtf(
	std::string const & file_name,
	core::io::StructFileRepOP sfr,
	core::io::StructFileRepOptions const & options
) {

	utility::io::ozstream file(file_name.c_str(), std::ios::out | std::ios::binary);
	if ( !file ) {
		TR.Error << "mmtf_writer:dump_mmtf Unable to open file:" << file_name << " for writing!!!" << std::endl;
		return false;
	}

	dump_mmtf(file, sfr, options);
	file.close();
	return true;
}


void
set_mmtf_default_options(core::io::StructFileRepOptions & options) {
	// mmtfs require all bonded atoms, we override and set this always
	options.set_write_all_connect_info(true);
	options.set_write_pdb_link_records(true);
}


::mmtf::GroupType
make_current_group(aiGroup const & ai_group) {
	AtomInformation const & ai_1 = ai_group[0];

	// basic parts from first atom
	::mmtf::GroupType gt;
	gt.groupName = ai_1.resName;
	bool is_unknown = core::chemical::is_aa_name_unknown(
		utility::strip(ai_1.resName));

	if ( !is_unknown ) {
		gt.singleLetterCode = core::chemical::oneletter_code_from_aa(
			core::chemical::aa_from_name(utility::strip(ai_1.resName)));
	} else {
		gt.singleLetterCode = 'X'; // should be '?' for non-polymer types
	}
	gt.chemCompType = ai_1.chem_comp_type;

	for ( auto ai : ai_group ) {
		gt.atomNameList.push_back(utility::strip(ai.name));
		gt.elementList.push_back(utility::strip(ai.element));
		gt.formalChargeList.push_back(ai.formalcharge);
	}
	return gt;
}

bool
is_in_bondAtomList(std::vector<int32_t> const & bondAtomList, core::Size lower_atom, core::Size upper_atom) {
	int32_t lower_candidate = std::min(lower_atom, upper_atom);
	int32_t upper_candidate = std::max(lower_atom, upper_atom);
	for ( core::Size i=0; i<bondAtomList.size(); i += 2 ) {
		int32_t lower_known = std::min(bondAtomList[i], bondAtomList[i+1]);
		int32_t upper_known = std::max(bondAtomList[i], bondAtomList[i+1]);
		if ( lower_candidate == lower_known && upper_candidate == upper_known )  {
			return true;
		}
	}
	return false;
}

std::map<core::Size, sd_index>
make_atom_num_to_sd_map(::mmtf::StructureData const & sd) {
	int32_t modelIndex = 0;
	int32_t chainIndex = 0;
	int32_t groupIndex = 0;
	core::Size atomSerial = 0;

	std::map<core::Size, sd_index> ret_map;

	for ( int32_t i = 0; i < sd.numModels; i++, modelIndex++ ) {
		// traverse chains
		for ( int32_t j = 0; j < sd.chainsPerModel[modelIndex]; j++, chainIndex++ ) {
			// traverse groups
			for ( int32_t k = 0; k < sd.groupsPerChain[chainIndex]; k++, groupIndex++ ) {
				const ::mmtf::GroupType& group =
					sd.groupList[sd.groupTypeList[groupIndex]];
				int groupAtomCount = group.atomNameList.size();
				for ( int32_t l = 0; l < groupAtomCount; l++, atomSerial++ ) {
					sd_index this_sd_index(chainIndex, groupIndex, l);
					ret_map[atomSerial] = this_sd_index;
				}
			}
		}
	}
	return ret_map;
}


aiChain
make_chain(utility::vector0<AtomInformation> const & chain_atoms) {
	aiChain AIC_out;
	decltype(chain_atoms.end()) upper;
	for ( auto lower=chain_atoms.begin(); lower != chain_atoms.end(); lower = upper ) {
		upper = std::upper_bound(chain_atoms.begin(), chain_atoms.end(), *lower, ai_cmp());
		AIC_out.emplace_back(lower, upper);
	}
	return AIC_out;
}

/*
* Convert linear sfr to vec0[vec0[vec0[AtomInformation]]]
* aka                   pose -> chain -> group -> atoms
*/
aiPose
aiPose_from_sfr(core::io::StructFileRep const & sfr) {
	aiPose AIP;
	for ( Size i = 0; i < sfr.chains().size(); ++i ) {
		if ( sfr.chains()[i].size() == 0 ) continue;
		aiChain this_chain = make_chain(sfr.chains()[i]);
		AIP.push_back(this_chain);
	}
	return AIP;
}


void
add_bonds_to_sd(::mmtf::StructureData & sd,
	aiPose & AIP, std::map<core::Size, sd_index> & atom_num_to_sd_map) {
	int32_t groupIndex = 0;
	int32_t chainIndex = 0;
	unsigned int atomIndex = 0;
	for ( core::Size i=0; i<AIP.size(); ++i, ++chainIndex ) {  // for each chain
		for ( core::Size j=0; j<AIP[i].size(); ++j, ++groupIndex ) {  // for each group
			for ( core::Size k=0; k<AIP[i][j].size(); ++k, ++atomIndex ) {  // for each atom
				AtomInformation const & ai = AIP[i][j][k];
				sd_index upper_connect = atom_num_to_sd_map[atomIndex];
				for ( core::Size l=1; l<=ai.connected_indices.size(); ++l ) {
					core::Size link_buddy = ai.connected_indices[l];
					if ( link_buddy > atomIndex+1 ) continue;
					sd_index lower_connect = atom_num_to_sd_map[link_buddy-1];
					if ( upper_connect.group_index == lower_connect.group_index ) {  // same group
						::mmtf::GroupType& group = sd.groupList[sd.groupTypeList[groupIndex]];
						bool already_in_grouplist = is_in_bondAtomList(group.bondAtomList, lower_connect.group_atom_index,
							upper_connect.group_atom_index); // based on group atom number
						if ( !already_in_grouplist ) {
							group.bondAtomList.push_back(lower_connect.group_atom_index);
							group.bondAtomList.push_back(upper_connect.group_atom_index);
							if ( ai.connected_orders[l] == 0 ) {
								group.bondOrderList.push_back((int8_t)-1);
								group.bondResonanceList.push_back((int8_t)1);
							} else {
								group.bondOrderList.push_back(ai.connected_orders[l]);
								group.bondResonanceList.push_back((int8_t)0);
							}
						}
						++sd.numBonds;
					} else {  // Not same group
						bool already_in_bondAtomList = is_in_bondAtomList(sd.bondAtomList, atomIndex,
							link_buddy-1);  // based on overall atom number
						if ( !already_in_bondAtomList ) {
							sd.bondAtomList.push_back(atomIndex);
							sd.bondAtomList.push_back(link_buddy-1);
							// Rosetta only has single bonds inter-group
							// TODO should a peptide bond be single?
							sd.bondOrderList.push_back((int8_t)1);
							sd.bondResonanceList.push_back((int8_t)-1);
						}
						++sd.numBonds;
					}
				}
			}
		}
	}
}

void
dump_mmtf(
	std::ostream & out,
	core::io::StructFileRepOP sfr,
	core::io::StructFileRepOptions const & options
) {
	::mmtf::StructureData sd;
	aiPose AIP = aiPose_from_sfr(*sfr);

	// ATOM/HETATM
	// hard set: will not change
	sd.numModels = 1;

	// will change
	sd.numBonds = 0;
	sd.numAtoms = 0;
	sd.numGroups = 0;
	sd.numChains = 0;
	sd.chainsPerModel.push_back((int32_t)AIP.size());
	std::string current_chain = "";

	// Add basic easy info
	for ( core::Size i=0; i<AIP.size(); ++i ) {  // for each chain
		for ( core::Size j=0; j<AIP[i].size(); ++j ) {  // for each group
			::mmtf::GroupType gt = make_current_group(AIP[i][j]);
			auto it = std::find(sd.groupList.begin(), sd.groupList.end(), gt);
			if ( it == sd.groupList.end() ) {
				sd.groupList.push_back(gt);
				sd.groupTypeList.push_back(sd.groupList.size()-1);
			} else {
				auto index = std::distance(sd.groupList.begin(), it);
				sd.groupTypeList.push_back(index);
			}
			for ( core::Size k=0; k<AIP[i][j].size(); ++k ) {  // for each atom
				AtomInformation const & ai = AIP[i][j][k];
				current_chain = ai.chainID;
				sd.xCoordList.push_back(ai.x);
				sd.yCoordList.push_back(ai.y);
				sd.zCoordList.push_back(ai.z);
				++sd.numAtoms;
			}
			++sd.numGroups;
			sd.groupIdList.push_back(AIP[i][j][0].resSeq);
		}
		sd.groupsPerChain.push_back(AIP[i].size());
		sd.chainIdList.push_back(current_chain);
		sd.chainNameList.push_back(current_chain);
		++sd.numChains;
	}

	// Bond information
	// atom num to chain->group->atom map
	std::map<core::Size, sd_index> atom_num_to_sd_map =  make_atom_num_to_sd_map(sd);

	// Add bonds to sd
	add_bonds_to_sd(sd, AIP, atom_num_to_sd_map);
	// TODO
	if ( options.preserve_header() ) {
		TR << "mmtf preserve_header not implemented yet!" << std::endl;
	}

	::mmtf::encodeToStream(sd, out);
}

} // mmtf
} // io
} // core
