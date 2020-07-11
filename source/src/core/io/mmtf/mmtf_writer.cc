// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/mmtf/mmtf_writer.cc
/// @brief Functions for MMTF writing.
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
#include <utility/version.hh>
#include <numeric/random/random.hh>

#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/Field.hh>
#include <core/io/HeaderInformation.hh>
#include <core/io/StructFileRep.hh>
#include <core/io/mmtf/ExtraDataEnum.hh>
#include <core/io/mmtf/ExtraDataEnumManager.hh>
#include <core/simple_metrics/SimpleMetricStruct.hh>

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

namespace core {
namespace io {
namespace mmtf {

bool
dump_mmtf(
	core::pose::Pose const & pose,
	std::string const & file_name,
	StructFileRepOptionsCOP options )
{
	StructFileRepOptionsOP local_options = options->clone();
	set_mmtf_default_options(*local_options);
	pose_to_sfr::PoseToStructFileRepConverter converter = pose_to_sfr::PoseToStructFileRepConverter( *local_options );

	converter.init_from_pose( pose );
	core::io::StructFileRepOP sfr =  converter.sfr();

	utility::io::ozstream file(file_name.c_str(), std::ios::out | std::ios::binary);
	if ( !file ) {
		TR.Error << "mmtf_writer:dump_mmtf Unable to open file:" << file_name << " for writing!!!" << std::endl;
		return false;
	}
	dump_mmtf( file, sfr, *local_options);
	file.close();
	return true;
}


/// @brief Dump an MMTF from a pose to an ostream.
core::io::StructFileRepOP
dump_mmtf(
	core::pose::Pose const & pose,
	std::ostream & out,
	StructFileRepOptionsCOP options)
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
	std::ostream & out)

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
	gt.groupName = utility::strip(ai_1.resName);
	bool const is_unknown = core::chemical::is_aa_name_unknown(gt.groupName);

	if ( !is_unknown ) {
		core::chemical::AA const aa(core::chemical::aa_from_name(gt.groupName));
		if ( core::chemical::is_canonical_L_aa_or_gly(aa) ) {
			gt.singleLetterCode = core::chemical::oneletter_code_from_aa(aa);
		} else gt.singleLetterCode = 'X';  // how do you get DNA slc?
	} else gt.singleLetterCode = 'X';  // should be '?' for non-polymer types
	gt.chemCompType = ai_1.chem_comp_type;

	for ( auto ai : ai_group ) {
		gt.atomNameList.push_back(utility::strip(ai.name));
		gt.elementList.push_back(utility::strip(ai.element));
		gt.formalChargeList.push_back(ai.formalcharge);
	}
	return gt;
}

bool
is_in_bondAtomList(std::vector<int32_t> const & bondAtomList, core::Size const lower_atom, core::Size const upper_atom) {
	int32_t const lower_candidate = std::min(lower_atom, upper_atom);
	int32_t const upper_candidate = std::max(lower_atom, upper_atom);
	for ( core::Size i=0; i<bondAtomList.size(); i += 2 ) {
		int32_t const lower_known = std::min(bondAtomList[i], bondAtomList[i+1]);
		int32_t const upper_known = std::max(bondAtomList[i], bondAtomList[i+1]);
		if ( lower_candidate == lower_known && upper_candidate == upper_known )  {
			return true;
		}
	}
	return false;
}

///@warning groups (aka residues) are made by grouping together identical ai.{chain,resSeq,resName,iCode}
///         I don't think there's a better way than that.
aiChain
make_chain(utility::vector0<AtomInformation> const & chain_atoms)
{
	aiChain AIC_out;
	std::map< std::tuple< char, core::Size, std::string, char >, core::Size >rsd_map;

	core::Size current_index(0);
	for ( core::Size i=0; i<chain_atoms.size(); ++i ) {
		core::io::AtomInformation const & ai(chain_atoms[i]);
		std::tuple< char, core::Size, std::string, char > const current_id(ai.chainID, ai.resSeq, ai.resName, ai.iCode);
		auto const search(rsd_map.find(current_id));
		if ( search != rsd_map.end() ) AIC_out[search->second].push_back(chain_atoms[i]);
		else {
			rsd_map[current_id] = current_index;
			if ( current_index >= AIC_out.size() ) AIC_out.resize(AIC_out.size()+50);  // try to cut down on resize
			AIC_out[current_index].push_back(chain_atoms[i]);
			++current_index;
		}
	}
	// remove unused, will happen if things aren't in sequence.
	AIC_out.erase(std::remove_if(AIC_out.begin(), AIC_out.end(),
		[](utility::vector0<AtomInformation> const & grp) { return grp.empty(); }), AIC_out.end());
	return AIC_out;
}

/*
* Convert linear sfr to vec0[vec0[vec0[AtomInformation]]]
* aka                   pose -> chain -> group -> atoms
*/
aiPose
aiPose_from_sfr(core::io::StructFileRep const & sfr) {
	aiPose AIP;
	for ( core::Size i = 0; i < sfr.chains().size(); ++i ) {
		if ( sfr.chains()[i].size() == 0 ) continue;
		aiChain this_chain = make_chain(sfr.chains().at(i));
		AIP.push_back(this_chain);
	}
	return AIP;
}

aiModels
aiModels_from_sfrs(utility::vector1< core::io::StructFileRepOP > const & sfrs)
{
	aiModels AIM;
	for ( auto const & sfr : sfrs ) {
		AIM.push_back(aiPose_from_sfr(*sfr));
	}
	return AIM;
}

void
add_bonds_to_sd(::mmtf::StructureData & sd,
	aiModels const & AIM, std::map<core::Size, sd_index> const & atom_num_to_sd_map)
{
	int32_t groupIndex = 0;
	int32_t chainIndex = 0; // unused
	int32_t modelIndex = 0;
	unsigned int atomIndex = 0;
	int group_bonds(0), inter_bonds(0);
	std::vector<core::Size>type_check;
	// TODO this function sucks :(
	for ( core::Size i=0; i<AIM.size(); ++i, ++modelIndex ) {  // for each model
		aiPose const & AIP(AIM[i]);
		for ( core::Size j=0; j<AIP.size(); ++j, ++chainIndex ) {  // for each chain
			for ( core::Size k=0; k<AIP[j].size(); ++k, ++groupIndex ) {  // for each group
				for ( core::Size l=0; l<AIP[j][k].size(); ++l, ++atomIndex ) {  // for each atom
					AtomInformation const & ai = AIP[j][k][l];
					sd_index const upper_connect = atom_num_to_sd_map.at(atomIndex);
					for ( core::Size m=1; m<=ai.connected_indices.size(); ++m ) {
						core::Size const& link_buddy(ai.connected_indices[m]);
						if ( link_buddy > atomIndex+1 ) continue;
						sd_index const lower_connect = atom_num_to_sd_map.at(link_buddy-1);
						if ( upper_connect.group_index == lower_connect.group_index ) {  // same group
							::mmtf::GroupType& group = sd.groupList[sd.groupTypeList[groupIndex]];
							bool const already_in_grouplist = is_in_bondAtomList(group.bondAtomList, lower_connect.group_atom_index,
								upper_connect.group_atom_index); // based on group atom number
							if ( !already_in_grouplist ) {
								group.bondAtomList.push_back(lower_connect.group_atom_index);
								group.bondAtomList.push_back(upper_connect.group_atom_index);
								switch ( ai.connected_orders[m] ) {
								case core::chemical::SingleBond:
								case core::chemical::DoubleBond:
								case core::chemical::TripleBond :
									group.bondOrderList.push_back((int8_t)ai.connected_orders[m]);
									group.bondResonanceList.push_back((int8_t)0);
									break;
								case core::chemical::AromaticBond :
									group.bondOrderList.push_back((int8_t)-1);
									group.bondResonanceList.push_back((int8_t)1);
									break;
								default :
									group.bondOrderList.push_back((int8_t)-1); // There's a bond, we don't know what it is.
									group.bondResonanceList.push_back((int8_t)1); // Unfortunately, the MMTF requires unknown bond orders to be resonant
									break;
								}
							}
							++group_bonds;
						} else {  // Not same group
							bool const already_in_bondAtomList = is_in_bondAtomList(sd.bondAtomList, atomIndex,
								link_buddy-1);  // based on overall atom number
							if ( !already_in_bondAtomList ) {
								sd.bondAtomList.push_back(link_buddy-1);
								sd.bondAtomList.push_back(atomIndex);
								// Rosetta only has single bonds inter-group
								// TODO should a peptide bond be single/not resonate?
								sd.bondOrderList.push_back((int8_t)1);
								sd.bondResonanceList.push_back((int8_t)-1);
							}
							++inter_bonds;
						}
						++sd.numBonds;
					}
				}
			}
		}
	}
}

// currently unused, but might be useful
template< typename T >
void
add_if_not_empty(
	std::string const & given_name,
	T const & content,
	std::map<std::string, msgpack::object> & target_map,
	msgpack::zone & zone)
{
	if ( !content.empty() ) {
		target_map[given_name] = msgpack::object(content, zone);
	}
}

template< typename T >
void
resize_and_add_if_not_empty(
	utility::vector1< core::io::StructFileRepOP > const & sfrs,
	std::string const & key,
	std::map<std::string, msgpack::object> & xxxProperties,
	msgpack::zone & zone,
	T const & fn)
{
	std::vector< msgpack::object > outdata;
	for ( core::Size i=1; i<=sfrs.size(); ++i ) {
		core::io::StructFileRepOP const & sfr(sfrs[i]);
		core::Size c_idx(i-1);
		auto const & c_data(fn(sfr));
		if ( c_data.empty() ) continue;
		else if ( outdata.empty() ) outdata.resize(sfrs.size(), msgpack::object());
		outdata[c_idx] = msgpack::object(fn(sfr), zone);
	}
	if ( !outdata.empty() ) xxxProperties[key] = msgpack::object(outdata, zone);
}

/* insert function here */
/*T const & input_data,*/ // might not need
/* core::Size const index, */
void
add_extra_data(
	::mmtf::StructureData & sd,
	utility::vector1< core::io::StructFileRepOP > const & sfrs,
	core::io::StructFileRepOptions const & options )
{
	// extra sfr data
	std::vector< std::map< std::string, std::string > > model_heterogen_names;
	std::vector< std::map< std::string, std::pair< std::string, std::string > > > model_residue_type_base_names;

	// 1. IO type helpers
	if ( options.use_pdb_format_HETNAM_records() ) {
		resize_and_add_if_not_empty(sfrs, "rosetta::heterogen_names", sd.modelProperties, sd.msgpack_zone,
			[](core::io::StructFileRepOP const & sfr) {return sfr->heterogen_names();});
	} else if ( !options.write_glycan_pdb_codes() ) {
		resize_and_add_if_not_empty(sfrs, "rosetta::residue_type_base_names", sd.modelProperties, sd.msgpack_zone,
			[](core::io::StructFileRepOP const & sfr) {return sfr->residue_type_base_names();});
	}

	// 2. reproducibility / logging info
	if ( ! options.integration_test_mode() ) { // Don't include version info for integration tests
		sd.extraProperties["rosetta::version"] = msgpack::object(utility::Version::version(), sd.msgpack_zone);
		sd.extraProperties["rosetta::commit"] = msgpack::object(utility::Version::commit(), sd.msgpack_zone);
		sd.extraProperties["rosetta::date_created"] = msgpack::object(utility::Version::date(), sd.msgpack_zone);
		if ( utility::Version::package() != "devel" ) {
			sd.extraProperties["rosetta::package"] = msgpack::object(utility::Version::package(), sd.msgpack_zone);
		}
		if ( utility::Version::revision() != "None" ) {
			sd.extraProperties["rosetta::revision"] = msgpack::object(utility::Version::revision(), sd.msgpack_zone);
		}
	}

	/// 3. SimpleMetrics -- These are experimental and may change!
	if ( options.mmtf_extra_data_io() ) {
		TR.Trace << "Writing extra data" << std::endl;
		ExtraDataEnumManager manager = ExtraDataEnumManager();

		resize_and_add_if_not_empty(sfrs, manager.enum_to_string(pose_cache_string_data), sd.modelProperties, sd.msgpack_zone,
			[](core::io::StructFileRepOP const & sfr) {return sfr->pose_cache_string_data();});

		resize_and_add_if_not_empty(sfrs, manager.enum_to_string(pose_cache_real_data), sd.modelProperties, sd.msgpack_zone,
			[](core::io::StructFileRepOP const & sfr) {return sfr->pose_cache_real_data();});

		resize_and_add_if_not_empty(sfrs, manager.enum_to_string(pdb_comments), sd.modelProperties, sd.msgpack_zone,
			[](core::io::StructFileRepOP const & sfr) {return sfr->pdb_comments();});

		resize_and_add_if_not_empty(sfrs, manager.enum_to_string(pose_cache_string_data), sd.modelProperties, sd.msgpack_zone,
			[](core::io::StructFileRepOP const & sfr) {return sfr->pose_cache_string_data();});

		//SimpleMetric Data conversion
		resize_and_add_if_not_empty(sfrs, manager.enum_to_string(simple_metric_real_data), sd.modelProperties, sd.msgpack_zone,
			[](core::io::StructFileRepOP const & sfr) {return sfr->simple_metric_data().real_data_;});

		resize_and_add_if_not_empty(sfrs, manager.enum_to_string(simple_metric_string_data), sd.modelProperties, sd.msgpack_zone,
			[](core::io::StructFileRepOP const & sfr) {return sfr->simple_metric_data().string_data_;});

		resize_and_add_if_not_empty(sfrs, manager.enum_to_string(simple_metric_composite_real_data), sd.modelProperties, sd.msgpack_zone,
			[](core::io::StructFileRepOP const & sfr) {return sfr->simple_metric_data().composite_real_data_;});

		resize_and_add_if_not_empty(sfrs, manager.enum_to_string(simple_metric_composite_real_data), sd.modelProperties, sd.msgpack_zone,
			[](core::io::StructFileRepOP const & sfr) {return sfr->simple_metric_data().composite_real_data_;});

		resize_and_add_if_not_empty(sfrs, manager.enum_to_string(simple_metric_composite_string_data), sd.modelProperties, sd.msgpack_zone,
			[](core::io::StructFileRepOP const & sfr) {return sfr->simple_metric_data().composite_string_data_;});

		resize_and_add_if_not_empty(sfrs, manager.enum_to_string(simple_metric_per_residue_real_data), sd.modelProperties, sd.msgpack_zone,
			[](core::io::StructFileRepOP const & sfr) {return sfr->simple_metric_data().per_residue_real_data_;});

		resize_and_add_if_not_empty(sfrs, manager.enum_to_string(simple_metric_per_residue_string_data), sd.modelProperties, sd.msgpack_zone,
			[](core::io::StructFileRepOP const & sfr) {return sfr->simple_metric_data().per_residue_string_data_;});

		resize_and_add_if_not_empty(sfrs, manager.enum_to_string(simple_metric_per_residue_real_output), sd.modelProperties, sd.msgpack_zone,
			[](core::io::StructFileRepOP const & sfr) {return sfr->simple_metric_data().per_residue_real_output_;});

		resize_and_add_if_not_empty(sfrs, manager.enum_to_string(simple_metric_per_residue_string_data), sd.modelProperties, sd.msgpack_zone,
			[](core::io::StructFileRepOP const & sfr) {return sfr->simple_metric_data().per_residue_string_data_;});

		resize_and_add_if_not_empty(sfrs, manager.enum_to_string(simple_metric_per_residue_string_output), sd.modelProperties, sd.msgpack_zone,
			[](core::io::StructFileRepOP const & sfr) {return sfr->simple_metric_data().per_residue_string_output_;});
	}
}

std::map<std::tuple<core::Size, core::Size, core::Size, core::Size>, core::Size>
make_AIM_to_atom_num(aiModels const & AIM) {
	core::Size idx(0);
	std::map<std::tuple<core::Size, core::Size, core::Size, core::Size>, core::Size> AIM_to_atom_num;
	for ( core::Size i=0; i<AIM.size(); ++i ) {  // for each model
		aiPose const & AIP(AIM[i]);
		for ( core::Size j=0; j<AIP.size(); ++j ) {  // for each chain
			for ( core::Size k=0; k<AIP[j].size(); ++k ) {  // for each group
				for ( core::Size l=0; l<AIP[j][k].size(); ++l ) {  // for each atom
					std::tuple<core::Size, core::Size, core::Size, core::Size> const mytup(i, j, k, l);
					AIM_to_atom_num[mytup] = idx;
					++idx;
				}
			}
		}
	}
	return AIM_to_atom_num;
}

void
dump_mmtf(
	std::ostream & out,
	utility::vector1< core::io::StructFileRepOP > sfrs,
	core::io::StructFileRepOptions const & options)
{
	::mmtf::StructureData sd(sfrs_to_sd(sfrs, options));
	sd.hasConsistentData(true, 4);
	::mmtf::encodeToStream(sd, out);
}

void
dump_mmtf(
	std::ostream & out,
	core::io::StructFileRepOP sfr,
	core::io::StructFileRepOptions const & options)
{
	::mmtf::StructureData sd(sfr_to_sd(sfr, options));
	sd.hasConsistentData(true, 4);
	::mmtf::encodeToStream(sd, out);
}


::mmtf::StructureData
sfr_to_sd(
	core::io::StructFileRepOP sfr,
	core::io::StructFileRepOptions const & options
) {
	utility::vector1< core::io::StructFileRepOP > sfrs;
	sfrs.push_back(sfr);
	return sfrs_to_sd(sfrs, options);
}


::mmtf::StructureData
sfrs_to_sd(
	utility::vector1< core::io::StructFileRepOP > sfrs,
	core::io::StructFileRepOptions const & options)
{
	aiModels AIM(aiModels_from_sfrs(sfrs));
	::mmtf::StructureData sd;

	sd.numBonds = 0;
	sd.numAtoms = 0;
	sd.numGroups = 0;
	sd.numChains = 0;
	sd.numModels = 0;
	for ( aiPose const & AIP : AIM ) sd.chainsPerModel.push_back((int32_t)AIP.size());

	std::string current_chain = "";
	core::Size current_model_prefix(0);

	// Add basic easy info
	for ( core::Size i=0; i<AIM.size(); ++i ) {  // for each model
		current_model_prefix = sd.numAtoms;
		aiPose const & AIP(AIM[i]);
		for ( core::Size j=0; j<AIP.size(); ++j ) {  // for each chain
			for ( core::Size k=0; k<AIP[j].size(); ++k ) {  // for each group
				// It is possible to set bonds here too, but no matter how you set them
				// it will be complicated, so might as well do it later :D
				::mmtf::GroupType const gt = make_current_group(AIP[j][k]);
				auto it = std::find(sd.groupList.begin(), sd.groupList.end(), gt);
				if ( it == sd.groupList.end() ) {
					sd.groupList.push_back(gt);
					sd.groupTypeList.push_back(sd.groupList.size()-1);
				} else {
					auto const index = std::distance(sd.groupList.begin(), it);
					sd.groupTypeList.push_back(index);
				}

				for ( core::Size l=0; l<AIP[j][k].size(); ++l ) {  // for each atom
					{ // alter bonding parterners for models>0
						AtomInformation & ai = AIM[i][j][k][l];
						if ( current_model_prefix != 0 ) {
							for ( auto & bp : ai.connected_indices ) bp += current_model_prefix;
						}
					}
					AtomInformation const & ai = AIP[j][k][l];
					current_chain = ai.chainID;
					sd.xCoordList.push_back(ai.x);
					sd.yCoordList.push_back(ai.y);
					sd.zCoordList.push_back(ai.z);
					sd.bFactorList.push_back(ai.temperature);
					sd.occupancyList.push_back(ai.occupancy);
					++sd.numAtoms;

					// Add insertion code always (it's easily compressed to 2 ints if not used)
					if ( l == 0 ) {
						if ( ai.iCode == ' ' ) sd.insCodeList.push_back(0x00);
						else sd.insCodeList.push_back(ai.iCode);
					}
				}
				++sd.numGroups;
				sd.groupIdList.push_back(AIP[j][k][0].resSeq);
			}
			sd.groupsPerChain.push_back(AIP[j].size());
			sd.chainIdList.push_back(current_chain);
			sd.chainNameList.push_back(current_chain);
			++sd.numChains;
		}
		++sd.numModels;
	}
	// Bond information
	// atom num to chain->group->atom map
	std::map<core::Size, sd_index> const atom_num_to_sd_map(core::io::mmtf::make_atom_num_to_sd_map(sd));

	// Add bonds to sd
	add_bonds_to_sd(sd, AIM, atom_num_to_sd_map);

	// Add extra data!
	add_extra_data(sd, sfrs, options);
	return sd;
}


} // mmtf
} // io
} // core
