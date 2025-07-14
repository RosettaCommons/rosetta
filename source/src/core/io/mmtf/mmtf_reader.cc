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

// When you move PDBReader and PoseUnbuilder, take these.
#include <utility/string_util.hh>
#include <utility/string_constants.hh>

#include <core/io/AtomInformation.hh>
#include <core/io/StructFileRep.hh>
#include <core/io/mmtf/ExtraDataEnum.hh>
#include <core/io/mmtf/ExtraDataEnumManager.hh>
#include <core/simple_metrics/SimpleMetricStruct.hh>



// Project headers
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers

// External headers
#include <mmtf.hpp>
//#include <mmtf/map_decoder.hpp>

// C++ headers
#include <algorithm>


static basic::Tracer TR( "core.io.mmtf.mmtf_reader" );

namespace core {
namespace io {
namespace mmtf {

// NOTE: this is templated to work with AtomInformation or LinkInformation
// LinkInformation is similar enough to work, but if you ever want to use
// sym1 or sym2 it will have to be split from this function
template < typename T >
inline void
add_xbond_information(
	std::map< ResID, utility::vector1< T > >& xbond_map,
	core::io::AtomInformation const & atm_1,
	core::io::AtomInformation const & atm_2)
{
	numeric::xyzVector<core::Real> const atm_1_xyz(atm_1.x, atm_1.y, atm_1.z);
	numeric::xyzVector<core::Real> const atm_2_xyz(atm_2.x, atm_2.y, atm_2.z);

	T xbond;
	xbond.name1 = utility::pad_atom_name(atm_1.name);
	xbond.resName1 = atm_1.resName;
	xbond.chainID1 = atm_1.chainID;
	xbond.resSeq1 = atm_1.resSeq;
	xbond.iCode1 = atm_1.iCode;
	xbond.resID1 = ResID( xbond.resSeq1, xbond.iCode1, xbond.chainID1);

	xbond.name2 = utility::pad_atom_name(atm_2.name);
	xbond.resName2 = atm_2.resName;
	xbond.chainID2 = atm_2.chainID;
	xbond.resSeq2 = atm_2.resSeq;
	xbond.iCode2 = atm_2.iCode;
	xbond.resID2 = ResID( xbond.resSeq2, xbond.iCode2, xbond.chainID2 );
	xbond.length = atm_1_xyz.distance(atm_2_xyz);

	if ( xbond_map.count( xbond.resID1 ) ) {
		xbond_map[ xbond.resID1 ].push_back(xbond);
		//links.push_back();
	} else {
		utility::vector1< T > links;
		links.push_back( xbond );
		xbond_map[ xbond.resID1 ] = links;
	}
}

template < typename T >
bool
sort_xbond_func( T const & lhs, T const & rhs ) {
	return ( lhs.chainID2 < rhs.chainID2 ) || ( lhs.chainID2 == rhs.chainID2 && lhs.resSeq2 < rhs.resSeq2 );
}

void
add_link_and_ss_information(
	::mmtf::StructureData const & sd,
	core::io::StructFileRep & sfr,
	std::vector< core::io::AtomInformation > const & all_AIs,
	std::vector< core::Size > const & ai_to_model,
	core::Size const & current_model
) {
	auto ai_sort_func = []( core::io::AtomInformation const & lhs, core::io::AtomInformation const & rhs ) {
		return ( lhs.chainID < rhs.chainID ) || ( lhs.chainID == rhs.chainID && lhs.resSeq < rhs.resSeq );
	};
	for ( core::Size i=0; i < sd.bondAtomList.size(); i += 2 ) {
		if ( ai_to_model[sd.bondAtomList[i]] != current_model || ai_to_model[sd.bondAtomList[i+1]] != current_model ) continue;
		// Make sure atm_1 is less than atm_2
		bool const less_than(ai_sort_func(all_AIs[sd.bondAtomList[i]], all_AIs[sd.bondAtomList[i+1]]));
		core::io::AtomInformation const &atm_1((less_than) ? all_AIs[sd.bondAtomList[i]] : all_AIs[sd.bondAtomList[i+1]]);
		core::io::AtomInformation const &atm_2((less_than) ? all_AIs[sd.bondAtomList[i+1]] : all_AIs[sd.bondAtomList[i]]);
		// remove unnecessary between N and C
		if ( (atm_2.resSeq - atm_1.resSeq == 1) && (atm_1.name == " C  ") && (atm_2.name == " N  " ) ) continue;
		// Try SSBond first
		if ( (atm_1.name == " SG ") && (atm_2.name == " SG " ) ) core::io::mmtf::add_xbond_information(sfr.ssbond_map(), atm_1, atm_2);
		else core::io::mmtf::add_xbond_information(sfr.link_map(), atm_1, atm_2);
	}

	// sort all links
	for ( auto & map_result : sfr.link_map() ) {
		std::sort( map_result.second.begin(), map_result.second.end(), sort_xbond_func<core::io::LinkInformation> );
	}
	for ( auto & map_result : sfr.ssbond_map() ) {
		std::sort( map_result.second.begin(), map_result.second.end(), sort_xbond_func<core::io::SSBondInformation> );
	}
}


core::io::AtomInformation
make_atom_information(
	::mmtf::StructureData const &sd,
	::mmtf::GroupType const & group,
	int const groupAtomIndex,
	core::Size const atomIndex,
	core::Size const atomSerial,
	core::Size const groupIndex,
	core::Size const chainIndex,
	utility::vector1<std::string> & known_chainIDs,
	core::io::StructFileReaderOptions const & options )
{
	core::io::AtomInformation ai;
	if ( !sd.entityList.empty() ) {
		ai.isHet = ::mmtf::is_hetatm(chainIndex, sd.entityList, group);
	} else ai.isHet = ::mmtf::is_hetatm(group.chemCompType.c_str());
	ai.chem_comp_type = group.chemCompType;

	if ( !::mmtf::isDefaultValue(sd.atomIdList) ) {
		ai.serial = sd.atomIdList[atomIndex];
	} else ai.serial = atomSerial;

	if ( group.atomNameList[groupAtomIndex].size() > 4 ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput,
			"We do not support atom names with size > 4. Failed at >" + group.atomNameList[groupAtomIndex] + "<");
	}
	// ** TODO for some brave soul: The majority of Rosetta can handle the stripped atom_names. However
	// The glycan code seems unable to do that.
	ai.name = utility::pad_atom_name(group.atomNameList[groupAtomIndex]); // can be up to 5 characters
	if ( !::mmtf::isDefaultValue(sd.altLocList) ) {
		ai.altLoc = sd.altLocList[atomIndex] == 0x00 ? ' ' : sd.altLocList[atomIndex];
	} else ai.altLoc = ' ';
	ai.resName = group.groupName; // can be up to 5 characters

	// NOTE we default to chainNameList because that's consistent with pdbs.
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
		ai.iCode = sd.insCodeList[groupIndex] == 0x00 ? ' ' : sd.insCodeList[groupIndex];
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
	ai.element = utility::strip(group.elementList[groupAtomIndex]);
	ai.terCount = 0;  // this gets set later
	return ai;
}


std::vector< core::io::AtomInformation >
make_all_atom_information(::mmtf::StructureData const & sd,
	std::vector< core::Size > & model_indexes,
	std::map< core::Size, core::Size > & model_index_to_starting_index,
	StructFileReaderOptions const & options )
{
	core::Size modelIndex = 0;
	core::Size chainIndex = 0;
	core::Size groupIndex = 0;
	core::Size atomIndex = 0;
	core::Size atomSerial = 0;
	utility::vector1<char> known_chainIDs;

	// Just make the bare-bones atom information
	std::vector< core::io::AtomInformation > all_AIs;

	// traverse models
	for ( int i = 0; i < sd.numModels; i++, modelIndex++ ) {
		model_index_to_starting_index[i] = atomSerial;
		// traverse chains
		for ( int j = 0; j < sd.chainsPerModel[modelIndex]; j++, chainIndex++ ) {
			// traverse groups
			for ( int k = 0; k < sd.groupsPerChain[chainIndex]; k++, groupIndex++ ) {
				const ::mmtf::GroupType& group =
					sd.groupList[sd.groupTypeList[groupIndex]];
				int groupAtomCount = group.atomNameList.size();
				for ( int l = 0; l < groupAtomCount; l++, atomIndex++ ) {
					++atomSerial;
					model_indexes.push_back((core::Size)i);
					all_AIs.push_back(
						core::io::mmtf::make_atom_information(
						sd, group, l, atomIndex, atomSerial, groupIndex, chainIndex, known_chainIDs, options));
				}
			}
		}
	}
	return all_AIs;
}





/// A warning for this function:
/// mmtf has following bond order options
/// Bond-order   Resonance  Explanation
/// -1             1  kekulized form is unavailable, but resonance is known
/// 1(or 2,3,4)    1  kekulized form is known, and resonance is known and exists
/// 1(or 2,3,4)    0  kekulized form is known, but resonance is nonexistant
/// 1(or 2,3,4)   -1  kekulized form is known, but resonance is not known
/// -Warning pt 2- mmtf supports up to 4 bonds, rosetta  only up to 3
///       i set 4x bonds to be 0 (unk) in rosetta
void
add_bond_information(::mmtf::StructureData const & sd,
	std::vector< core::io::AtomInformation > & all_AIs,
	std::map<core::Size, sd_index> const & atom_num_to_sd_map)
{
	for ( core::Size i=0; i<all_AIs.size(); ++i ) {
		sd_index const & c_sd_index(atom_num_to_sd_map.at(i));
		::mmtf::GroupType const & group(
			sd.groupList[sd.groupTypeList[c_sd_index.group_index]]);
		for ( core::Size j=0,k=0; j<group.bondAtomList.size(); j+=2, ++k ) {
			if ( c_sd_index.group_atom_index == group.bondAtomList[j] ) {
				int const difference((group.bondAtomList[j+1]-group.bondAtomList[j]));
				// i = 0-indexed, so add 1
				all_AIs[i].connected_indices.push_back(i+difference+1);
				// bond orders are NOT required for mmtf
				// bond Resonance is also not required
				if ( !::mmtf::isDefaultValue(group.bondOrderList) ) {
					if ( group.bondOrderList[k] == -1 || group.bondOrderList[k] == 4 ) {
						all_AIs[i].connected_orders.push_back( core::chemical::UnknownBond );
					} else {
						all_AIs[i].connected_orders.push_back( core::chemical::BondName(group.bondOrderList[k]) );
					}
				}
			} else if ( c_sd_index.group_atom_index == group.bondAtomList[j+1] ) {
				int const difference((group.bondAtomList[j]-group.bondAtomList[j+1]));
				all_AIs[i].connected_indices.push_back(i+difference+1);
				// bond orders are NOT required for mmtf
				// bond Resonance is also not required
				if ( !::mmtf::isDefaultValue(group.bondOrderList) ) {
					if ( group.bondOrderList[k] == -1 || group.bondOrderList[k] == 4 ) {
						all_AIs[i].connected_orders.push_back( core::chemical::UnknownBond );
					} else {
						all_AIs[i].connected_orders.push_back( core::chemical::BondName(group.bondOrderList[k]) );
					}
				}
			}
		}
	}
	for ( core::Size i=0, j=0; i<sd.bondAtomList.size(); i+=2, ++j ) {
		core::io::AtomInformation & atm1(all_AIs[sd.bondAtomList[i]]);
		core::io::AtomInformation & atm2(all_AIs[sd.bondAtomList[i+1]]);
		atm1.connected_indices.push_back(sd.bondAtomList[i+1]+1);
		atm2.connected_indices.push_back(sd.bondAtomList[i]+1);
		if ( !::mmtf::isDefaultValue(sd.bondOrderList) && !::mmtf::isDefaultValue(sd.bondResonanceList) ) {
			if ( sd.bondOrderList[j] == -1 || sd.bondOrderList[j] == 4 ) {
				atm1.connected_orders.push_back( core::chemical::UnknownBond );
				atm2.connected_orders.push_back( core::chemical::UnknownBond );
			} else {
				atm1.connected_orders.push_back( core::chemical::BondName(sd.bondOrderList[j]) );
				atm2.connected_orders.push_back( core::chemical::BondName(sd.bondOrderList[j]) );
			}
		}
	}
}


// mmtf doesn't have TER... so instead look for bonds to the C term and
// use that to set terCount
// I'm sure there's a better way to do this, but thatll do pig... thatll do.
// TODO carbohydrates and DNA/RNA? ideally we could use a smarter graph system
//      to do this, since we can trace where polymers start and stop.
void
add_ters_via_bonds(std::vector< core::io::AtomInformation > & all_AIs,
	std::vector< core::Size > const & ai_to_model)
{
	// 1. Find all things bound to N and all Cs
	std::set< core::Size > all_things_bound_to_N, all_Cs, terminal_Cs;
	for ( core::Size i=0; i<all_AIs.size(); ++i ) {
		if ( all_AIs[i].name == " C  " ) all_Cs.insert(i);
		else if ( all_AIs[i].name == " N  " ) {
			for ( core::Size const & bound_idx : all_AIs[i].connected_indices ) {
				all_things_bound_to_N.insert(bound_idx-1);
			}
		}
	}
	// 2. find the Cs not bound to Ns via checking all things bound to N
	std::set_difference(
		all_Cs.begin(), all_Cs.end(),
		all_things_bound_to_N.begin(), all_things_bound_to_N.end(),
		std::inserter(terminal_Cs, terminal_Cs.end()));

	// 3. Use Terminal Cs to set AtomInformation terCount.  it accumulates
	// so you have to loop through all.  Note! reset the counter at 0 for new
	// models!
	int terCount(0), current_residue(0);
	std::string current_chain("^");
	core::Size current_model(0);
	for ( core::Size i=0; i<all_AIs.size(); ++i ) {
		if ( i > 0 && ai_to_model[i] != current_model ) {
			current_residue = 0;
			current_chain = "^";
			current_model = ai_to_model[i];
			terCount = 0;
		}
		if ( terminal_Cs.count(i) ) {
			current_residue = all_AIs[i].resSeq;
			current_chain = all_AIs[i].chainID;
		}
		if ( current_residue != 0 && (current_residue != all_AIs[i].resSeq || current_chain != all_AIs[i].chainID) ) {
			current_residue = 0;
			current_chain = "^";
			++terCount;
		}
		all_AIs[i].terCount = terCount;
	}
}



core::io::StructFileRepOP
create_sfr_from_mmtf_filename( std::string const & mmtf_filename,
	StructFileReaderOptions const & options )
{
	utility::vector1<core::Size> const first_only({0});
	utility::vector1< core::io::StructFileRepOP > sfrs(
		create_sfrs_from_mmtf_filename(mmtf_filename, options, first_only));
	return sfrs.front();
}


template < typename T >
void
set_model_index_if_not_empty(core::Size const model_index,
	std::string const & info_tag,
	std::vector< T > const & all_model_data, T & target_data)
{
	if ( all_model_data.empty() ) return;
	if ( model_index >= all_model_data.size() ) {
		std::stringstream ss;
		ss << "Unable to load model index: " << model_index <<
			" for data type: " << info_tag << " found model size to be: " <<
			all_model_data.size() << ". Remember it's std::vector, not vector1!" <<
			std::endl;
		throw CREATE_EXCEPTION(
			utility::excn::BadInput,
			ss.str());
	}
	T const & current_model_data = all_model_data[model_index];
	if ( !current_model_data.empty() ) {
		target_data = current_model_data;
	}
}


/// Warning. This *might* be inefficient! time will tell
/// We access and decode map_key for every model... that could be
/// expensive... but maybe it's not that bad! otherwise we have to
/// find a way to decode all at once....
template < typename T >
void
map_decode_from_model(
	::mmtf::MapDecoder const & md,
	std::string const & map_key,
	bool const required,
	core::Size const model_index,
	T & target)
{
	// We could T here, but I think it will be cheaper to keep as
	// msgpack::object.
	std::vector< msgpack::object > target_model_data;
	md.decode(map_key, required, target_model_data);

	if ( model_index+1 <= target_model_data.size() ) {
		target_model_data[model_index].convert_if_not_nil(target);
	} else if ( required ) {
		std::stringstream ss;
		ss << "Unable to load model index: " << model_index <<
			" for key: " << map_key << " and it was set to be required..." <<
			" found size: " << target_model_data.size() << std::endl;
		throw CREATE_EXCEPTION(
			utility::excn::BadInput,
			ss.str());
	}
}


void
read_extra_data(
	core::io::StructFileRepOP & sfr,
	::mmtf::StructureData const & sd,
	core::Size const model_index )
{
	// ::mmtf::MapDecoder const extraProperties_MD(sd.extraProperties); // to use at some point
	::mmtf::MapDecoder const modelProperties_MD(sd.modelProperties);

	///// 1. Rosetta IO
	map_decode_from_model(modelProperties_MD, "rosetta::heterogen_names", false,
		model_index, sfr->heterogen_names());

	std::map< std::string, std::pair< std::string, std::string > > base_names;
	map_decode_from_model(modelProperties_MD, "rosetta::residue_type_base_names", false,
		model_index, base_names);
	for (auto const & pairing: base_names ) {
		sfr->residue_type_base_names()[ resid_from_tag(pairing.first) ] = pairing.second;
	}

	///// 2. SimpleMetrics
	ExtraDataEnumManager manager = ExtraDataEnumManager();

	map_decode_from_model(modelProperties_MD, manager.enum_to_string(pose_cache_string_data),
		false /*required*/, model_index, sfr->pose_cache_string_data());

	map_decode_from_model(modelProperties_MD, manager.enum_to_string(pose_cache_real_data),
		false /*required*/, model_index, sfr->pose_cache_real_data());

	map_decode_from_model(modelProperties_MD, manager.enum_to_string(pdb_comments),
		false /*required*/, model_index, sfr->pdb_comments());

	core::simple_metrics::SimpleMetricStruct & sm_info = sfr->simple_metric_data();

	map_decode_from_model(modelProperties_MD, manager.enum_to_string(simple_metric_real_data),
		false /*required*/, model_index, sm_info.real_data_);

	map_decode_from_model(modelProperties_MD, manager.enum_to_string(simple_metric_string_data),
		false /*required*/, model_index, sm_info.string_data_);

	map_decode_from_model(modelProperties_MD, manager.enum_to_string(simple_metric_composite_real_data),
		false /*required*/, model_index, sm_info.composite_real_data_);

	map_decode_from_model(modelProperties_MD, manager.enum_to_string(simple_metric_composite_string_data),
		false /*required*/, model_index, sm_info.composite_string_data_);

	map_decode_from_model(modelProperties_MD, manager.enum_to_string(simple_metric_per_residue_real_data),
		false /*required*/, model_index, sm_info.per_residue_real_data_);

	map_decode_from_model(modelProperties_MD, manager.enum_to_string(simple_metric_per_residue_string_data),
		false /*required*/, model_index, sm_info.per_residue_string_data_);

	map_decode_from_model(modelProperties_MD, manager.enum_to_string(simple_metric_per_residue_real_output),
		false /*required*/, model_index, sm_info.per_residue_real_output_);

	map_decode_from_model(modelProperties_MD, manager.enum_to_string(simple_metric_per_residue_string_output),
		false /*required*/, model_index, sm_info.per_residue_string_output_);
}


std::map<core::Size, core::Size>
make_old_to_new_mapping(
	std::vector<std::vector<core::Size>> const atom_idx_chains,
	core::Size const starting_index)
{
	std::map<core::Size, core::Size>old_to_new_idx;
	core::Size count(0);
	for ( core::Size i=0; i<atom_idx_chains.size(); ++i ) {
		for ( core::Size j=0; j<atom_idx_chains[i].size(); ++j ) {
			old_to_new_idx[atom_idx_chains[i][j]-starting_index] = count;
			++count;
		}
	}
	return old_to_new_idx;
}

void
update_bond_indices(core::io::StructFileRep & sfr,
	std::map<core::Size, core::Size> const & old_to_new_idx)
{
	for ( core::Size i=0; i<sfr.chains().size(); ++i ) {
		for ( core::Size j=0; j<sfr.chains()[i].size(); ++j ) {
			for ( core::Size & connected_index : sfr.chains()[i][j].connected_indices ) {
				connected_index = old_to_new_idx.at(connected_index-1)+1;
			}
		}
	}
}



utility::vector1< core::io::StructFileRepOP >
create_sfrs_from_mmtf_filename(
	std::string const & mmtf_filename,
	core::io::StructFileReaderOptions const & options,
	utility::vector1< core::Size > const & model_indexes)
{
	if ( options.read_pdb_header() ) {
		TR.Warning << "read_pdb_header set, but not used in reading mmtf_file" << std::endl;
		// TODO
	}
	::mmtf::StructureData sd;
	::mmtf::decodeFromFile(sd, mmtf_filename);

	std::vector< core::Size > ai_to_model;
	std::map< core::Size, core::Size > model_index_to_starting_index;

	std::map<core::Size, sd_index> const atom_num_to_sd_map(make_atom_num_to_sd_map(sd));
	std::vector< core::io::AtomInformation > all_AIs(core::io::mmtf::make_all_atom_information(sd, ai_to_model, model_index_to_starting_index, options));

	// OK this is a little hairy, but I think we can handle it.
	// - mmtf has connected_indices (bonds) based on the index of the atom in the mmtf file + the current group
	// - Rosetta takes the atoms and sorts them into different chains thereby invalidating some bonds
	// - To bring the bonds back into a validated state we have to alter the connected_indices to the new indicies
	core::io::mmtf::add_bond_information(sd, all_AIs, atom_num_to_sd_map);
	core::io::mmtf::add_ters_via_bonds(all_AIs, ai_to_model);

	// alter correct bond_indices for models
	for ( core::Size i=0; i<all_AIs.size(); ++i ) {
		if ( ai_to_model[i] != 0 ) {
			core::Size const starting_index(model_index_to_starting_index[ai_to_model[i]]);
			for ( auto & bp : all_AIs[i].connected_indices ) bp -= starting_index;
			// only reset serial if atomIdList isn't set, because then we manually set
			if ( ::mmtf::isDefaultValue(sd.atomIdList) ) all_AIs[i].serial -= starting_index;
		}
	}

	utility::vector1< core::io::StructFileRepOP > sfrs;

	for ( core::Size const & model_index : model_indexes ) {
		core::io::StructFileRepOP sfr( utility::pointer::make_shared<core::io::StructFileRep>() );
		std::map<std::string, core::io::ChainAtoms> atom_chain_map;
		std::map<std::string, std::vector<core::Size>> atom_idx_chain_map;
		std::vector< std::string > chain_list; // preserve order
		for ( core::Size i=0; i<all_AIs.size(); ++i ) {
			core::io::AtomInformation const & ai(all_AIs[i]);
			if ( model_index != ai_to_model[i] ) continue;
			atom_chain_map[ai.chainID].push_back(ai);
			atom_idx_chain_map[ai.chainID].push_back(i);
			if ( std::find( chain_list.begin(), chain_list.end(), ai.chainID ) == chain_list.end() ) {
				chain_list.push_back( ai.chainID );
			}
		}
		add_link_and_ss_information(sd, *sfr, all_AIs, ai_to_model, model_index);

		std::vector<std::vector<core::Size>> atom_idx_chains;
		for ( std::string const & i : chain_list ) { // std::vector
			sfr->chains().push_back( atom_chain_map.find( i )->second );
			atom_idx_chains.push_back( atom_idx_chain_map.find( i )->second );
		}
		std::map<core::Size, core::Size> const old_to_new_idx(make_old_to_new_mapping(atom_idx_chains, model_index_to_starting_index[model_index]));
		update_bond_indices(*sfr, old_to_new_idx);
		sfr->filename() = mmtf_filename;
		read_extra_data(sfr, sd, model_index);
		sfrs.push_back(sfr);
	}
	return sfrs;
}

} //mmtf
} //io
} //core
