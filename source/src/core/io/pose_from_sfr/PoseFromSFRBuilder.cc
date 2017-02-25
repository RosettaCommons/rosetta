// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/pose_builder/PoseFromSFRBuilder.cc
/// @brief  Method definitions for PoseFromSFRBuilder class and related classes.
/// @author Sergey Lyskov
/// @author Andrew Leaver-Fay


// Note: AVOID ACCESSING THE OPTIONS SYSTEM DIRECTLY IN THIS FILE, ESPECIALLY FOR PDB INPUT!
// Doing so will mean that two different Poses cannot be constructed differently.
// Instead, modify StructFileRepOptions to include the option.


// Unit headers
#include <core/io/pose_from_sfr/PoseFromSFRBuilder.hh>
#include <core/io/pose_from_sfr/chirality_resolution.hh>

// Package headers
#include <core/io/StructFileRepOptions.hh>
#include <core/io/StructFileRep.hh>
#include <core/io/HeaderInformation.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/NomenclatureManager.hh>
#include <core/io/ResidueInformation.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedAtomID_Map.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>
#include <core/chemical/types.hh>
#include <core/chemical/rings/RingConformerSet.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/carbohydrates/GlycanTreeSetObserver.hh>
#include <core/pose/util.hh>
#include <core/pose/ncbb/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/cryst/util.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

// External headers
#include <ObjexxFCL/format.hh>
#include <boost/lexical_cast.hpp>

// C++ headers
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <utility>


namespace core {
namespace io {
namespace pose_from_sfr {

// Tracer instance for this file
static THREAD_LOCAL basic::Tracer TR( "core.io.pose_from_sfr.PoseFromSFRBuilder" );

using namespace ObjexxFCL;


PoseFromSFRBuilder::PoseFromSFRBuilder( chemical::ResidueTypeSetCOP rts, StructFileRepOptions const & options ) :
	residue_type_set_( rts ),
	options_( options ),
	coordinates_assigned_( false ),
	outputted_ignored_water_warning_( false )
{}

PoseFromSFRBuilder::~PoseFromSFRBuilder() {}

/// @details The process of building a Pose occurs in several phases.  First, in the setup() function,
/// the builder converts the SFR into a vector of ResidueInformation objects -- renaming some atoms along
/// the way. The next five phases involve five passes over this ResidueInformation vector.
void
PoseFromSFRBuilder::build_pose( StructFileRep const & sfr, pose::Pose & pose )
{
	pose.clear();

	setup( sfr );
	pass_1_merge_residues_as_necessary();
	pass_2_resolve_residue_types();
	pass_3_verify_sufficient_backbone_atoms();
	pass_4_redo_termini();
	pass_5_note_discarded_atoms();
	build_initial_pose( pose );
	refine_pose( pose );
	// build_pdb_info( pose );
}

id::AtomID_Mask const &
PoseFromSFRBuilder::missing_atoms() const
{
	return missing_;
}

bool
missing_O2prime( utility::vector1< core::io::AtomInformation > const & atoms )
{
	for ( Size n = 1; n <= atoms.size(); n++ ) {
		std::string const & name =  atoms[ n ].name;
		if ( name == " O2*" || name == " O2'" )  return false;
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////
//  introduced in 2013 with changes of RNA and DNA atom types to match PDB -- rhiju.
//  Could also include cleanup/handling of chirality of hydrogens (e.g., H5' <--> H5'') in here,
//  but actually easier to do it later when we actually are instantiating ResidueType and have Rosetta's
//  ideal coordinates.
void
PoseFromSFRBuilder::convert_nucleic_acid_residue_info_to_standard()
{
	// following is to show warnings or cap number.
	static Size nfix( 0 );
	static Size const max_fix( 2 );
	static bool const show_all_fixup( options_.show_all_fixes() );
	static bool showed_warning( false );

	for ( Size ii = 1; ii <= rinfos_.size(); ++ii ) {
		core::io::ResidueInformation & rinfo  = rinfos_[ ii ];
		std::string const original_name = rinfo.resName();

		// first establish if this is DNA or RNA (or something else)
		if ( !options_.guarantee_no_DNA()
				&& core::io::NomenclatureManager::is_old_DNA( rinfo.resName() )
				&& missing_O2prime( rinfo.atoms() ) )  {
			std::string new_name( original_name ); new_name.replace( 1, 1, "D" ); // A --> dA
			rinfo.resName( new_name );
			if ( ++nfix <= max_fix || show_all_fixup ) {
				TR << "Converting residue name " <<  original_name <<  " to " << rinfo.resName() << std::endl;
			}
		}

		if ( core::io::NomenclatureManager::is_old_RNA( rinfo.resName() ) ) {
			std::string new_name( original_name ); new_name.replace( 1, 1, " " ); // rA --> A
			rinfo.resName( new_name );
			if ( ++nfix <= max_fix || show_all_fixup ) {
				TR << "Converting residue name " <<  original_name <<  " to " << rinfo.resName() << std::endl;
			}
		}

		for ( Size jj = 1 ; jj <= rinfo.atoms().size(); jj++ ) {
			AtomInformation const & atom_info = rinfo.atoms()[jj];
			std::string const original_atom_name = atom_info.name;
			//  final stars (*)  are changed to primes (').
			// Don't assume atom names are a fixed length (mmCIF)
			if ( atom_info.name.size() >= 1 && atom_info.name[atom_info.name.size()-1] == '*' ) {
				std::string new_atom_name = atom_info.name.substr(0,atom_info.name.size()-1) + "\'";
				rinfo.rename_atom( original_atom_name, new_atom_name );
				if ( ++nfix <= max_fix || show_all_fixup ) {
					TR << "Converting atom name    " << original_atom_name << " to " << new_atom_name << std::endl;
				}
			}
		}
	}

	if ( nfix > max_fix && !show_all_fixup && !showed_warning ) {
		TR << "Number of nucleic acid residue fixups exceeds output limit. ";
		TR << "Rerun with -show_all_fixes to show everything." << std::endl;
		showed_warning = true;
	}
}

void
PoseFromSFRBuilder::setup( StructFileRep const & sfr ) {
	// clear any data that might have existed from previous build_pose executions
	missing_.resize( 0 );
	// deep copy of the input data, which we will be modifying.
	sfr_ = sfr;

	// Prune out LINK records that refer to metal residues here.
	// These are not appropriately represented by the variant types we handle
	// in the -in:auto_setup_metals code.

	// Likewise, prune LINK records that refer to saccharides unless the -include sugars flag is on.
	// And LINK records to or from unrecognized residues?
	using namespace core::chemical;
	typedef utility::vector1< LinkInformation > LinkInformationVect;
	typedef std::map< std::string, utility::vector1< LinkInformation > > LinkMap;

	// AMW: In CIF nomenclature, LINKs to be handled:
	// This is actually a "close contact" but it's O3' of adenine to C of FME, UPPER-UPPER
	// 1  1 "O3'" v A   76   ? ? C     v FME 77   ? ? 1.56

	LinkMap pruned_links;
	for ( LinkMap::const_iterator iter = sfr.link_map().begin(), iter_end = sfr.link_map().end();
			iter != iter_end; ++iter ) {
		std::string const & resID1 = iter->first;
		LinkInformationVect const & links = iter->second;
		for ( LinkInformationVect::const_iterator link_iter = links.begin(), link_iter_end = links.end();
				link_iter != link_iter_end; ++link_iter ) {
			LinkInformation const & link = *link_iter;

			bool const link_has_metal( NomenclatureManager::is_metal( link.resName1 ) ||
				NomenclatureManager::is_metal( link.resName2 ) );
			bool const link_has_sugar( NomenclatureManager::is_sugar( link.resName1 ) ||
				NomenclatureManager::is_sugar( link.resName2 ) );

			if ( link_has_metal ) {
				TR.Debug << "Omitting LINK record that uses a metal. These will be processed ";
				TR.Debug << "by -in:auto_setup_metals." << std::endl;
			} else if ( options_.ignore_sugars() && link_has_sugar ) {
				TR.Debug << "Omitting LINK record that uses a saccharide residue. ";
				TR.Debug << "Did you mean to use the -include_sugars flag?" << std::endl;
			} else if ( ! link_has_sugar && (  // PDB sugar codes will not be recognized by the ResidueTypeFinder!
					! ResidueTypeFinder( *residue_type_set_ ).name3( link.resName1 ).get_representative_type() ||
					! ResidueTypeFinder( *residue_type_set_ ).name3( link.resName2 ).get_representative_type() ) ) {
				// One or more residues in this LINK is not recognized.  Move on!
				TR.Debug << "Omitting LINK record that uses an unrecognized residue." << std::endl;
			} else if ( ( link.resSeq1 == link.resSeq2 - 1 && link.name1 == " O3'" && link.name2 == " P  " )
					||  ( link.resSeq1 == link.resSeq2 + 1 && link.name1 == " P  " && link.name2 == " O3'" ) ) {
				// We have a normal polymeric connection written as a LINK.
				TR.Debug << "Omitting LINK record that represents the canonical polymeric connectivity of a NCNT." << std::endl;
			} else if ( ( link.resSeq1 == link.resSeq2 - 1 && link.name1 == " C  " && link.name2 == " N  " )
					||  ( link.resSeq1 == link.resSeq2 + 1 && link.name1 == " N  " && link.name2 == " C  " ) ) {
				// We have a normal polymeric connection written as a LINK.
				TR.Debug << "Omitting LINK record that represents the canonical polymeric connectivity of a NCAA." << std::endl;
			} else if ( link.resName1 == "CYS" && link.resName2 == "CYS" && utility::strip(link.name1) == "SG" && utility::strip(link.name2) == "SG" ) {
				// We have an SSBOND redundantly specified as a LINK.
				TR.Debug << "Omitting LINK record that gives a SECOND specification of a disulfide bond." << std::endl;
			} else {
				if ( pruned_links.count( resID1 ) ) {
					pruned_links[ resID1 ].push_back( link );
				} else {
					pruned_links[ resID1 ] = LinkInformationVect( 1, link );
				}
			}
		}
	}
	sfr_.link_map() = pruned_links;
}

void
PoseFromSFRBuilder::pass_1_merge_residues_as_necessary()
{
	using namespace core::io::pdb;
	using namespace core::chemical;

	create_working_data( options_, sfr_, rinfos_ );
	convert_nucleic_acid_residue_info_to_standard();

	rosetta_residue_name3s_.resize(     rinfos_.size() );
	residue_types_.resize(              rinfos_.size() );
	is_lower_terminus_.resize(          rinfos_.size(), false );
	is_branch_lower_terminus_.resize(   rinfos_.size(), false );
	same_chain_prev_.resize(            rinfos_.size(), false );
	residue_was_recognized_.resize(     rinfos_.size(), true );
	remapped_atom_names_.resize(        rinfos_.size() );
	merge_behaviors_.resize(            rinfos_.size() );
	merge_atom_maps_.resize(            rinfos_.size() );

	for ( Size ii = 1; ii <= rinfos_.size(); ++ii ) {
		MergeBehaviorManager::ResidueMergeInstructions const & behavior_map_pair =
			residue_type_set_->merge_behavior_manager().merge_behavior_for_name3( rinfos_[ ii ].resName() );

		merge_behaviors_[ ii ] = behavior_map_pair.first;
		if ( merge_behaviors_[ ii ] != mrb_do_not_merge ) {
			merge_atom_maps_[ ii ] = behavior_map_pair.second;
		}

		if ( merge_behaviors_[ ii ] == mrb_merge_w_prev || merge_behaviors_[ ii ] == mrb_merge_w_next ) {
			// put the atoms from residue ii into the core::io::ResidueInformation for residue ii-1
			if ( ii == 1 && merge_behaviors_[ ii ] == mrb_merge_w_prev ) {
				utility_exit_with_message( "Residue 1, \"" + rinfos_[ ii ].resName() + "\" has been indicated"
					" to merge with the previous residue from the core::io::NomenclatureManager" );
			}

			if ( ii == rinfos_.size() && merge_behaviors_[ ii ] == mrb_merge_w_next ) {
				utility_exit_with_message( "The last residue, residue" + utility::to_string( ii ) + " named \"" + rinfos_[ ii ].resName() +
					"\" has been indicated to merge with the next residue from the core::io::NomenclatureManager" );
			}

			core::io::ResidueInformation & rmerged_into = rinfos_[ merge_behaviors_[ ii ] == mrb_merge_w_prev ? ii-1 : ii+ 1 ];

			chemical::MergeBehaviorManager::AtomRenamingMap nowhitespace_rename_map;
			for ( chemical::MergeBehaviorManager::AtomRenamingMap::const_iterator
					iter = merge_atom_maps_[ ii ].begin(), iter_end = merge_atom_maps_[ ii ].end();
					iter != iter_end; ++iter ) {
				nowhitespace_rename_map[ utility::strip( iter->first ) ] = iter->second;
			}

			for ( Size jj = 1; jj <= rinfos_[ ii ].atoms().size(); ++jj ) {
				AtomInformation jjatom = rinfos_[ ii ].atoms()[ jj ];
				if ( merge_atom_maps_[ ii ].count( jjatom.name ) ) {
					TR << "Renaming atom " << jjatom.name << " as " << merge_atom_maps_[ ii ][ jjatom.name ] << std::endl;
					jjatom.name = merge_atom_maps_[ ii ][ jjatom.name ];
				} else if ( nowhitespace_rename_map.count( jjatom.name ) ) {
					TR << "Renaming atom " << jjatom.name << " as " << nowhitespace_rename_map[ jjatom.name ] << std::endl;
					jjatom.name = nowhitespace_rename_map[ jjatom.name ];
				}
				rmerged_into.append_atom( jjatom );
			}
		}
	}
}

/// @details at the end of this pass, the residue_types_ array will have a pointer to the
/// ResidueType used for each entry in the rinfos_ array.  Any null pointer entry in this
/// array should be interpreted as an unrecognized, unreadable, or otherwise unread residue.
void
PoseFromSFRBuilder::pass_2_resolve_residue_types()
{
	using namespace core::io::pdb;
	using namespace core::chemical;

	Size const nres_pdb( rinfos_.size() );
	for ( Size ii = 1; ii <= rinfos_.size(); ++ii ) {
		if ( merge_behaviors_[ ii ] == mrb_merge_w_next || merge_behaviors_[ ii ] == mrb_merge_w_prev ) continue;

		core::io::ResidueInformation const & rinfo = rinfos_[ ii ];
		char chainID = rinfo.chainID();
		std::string const & pdb_name = rinfo.resName();
		std::string const & resid = rinfo.resid();

		// Convert PDB 3-letter code to Rosetta 3-letter code, if a list of alternative codes has been provided.
		std::pair< std::string, std::string > const & rosetta_names(
			NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( pdb_name ) );
		std::string const & name3( rosetta_names.first );
		rosetta_residue_name3s_[ ii ] = name3;
		if ( rosetta_names.second != "" ) {
			sfr_.residue_type_base_names()[ resid ] = std::make_pair( name3, rosetta_names.second );
		}

		bool const separate_chemical_entity = determine_separate_chemical_entity( chainID );
		bool same_chain_prev = determine_same_chain_prev( ii, separate_chemical_entity );
		bool same_chain_next = determine_same_chain_next( ii, separate_chemical_entity );
		bool const check_Ntermini_for_this_chain = determine_check_Ntermini_for_this_chain( chainID );
		bool const check_Ctermini_for_this_chain = determine_check_Ctermini_for_this_chain( chainID );

		bool is_branch_point( false ), is_branch_lower_terminus( false );
		utility::vector1< std::string > branch_points_on_this_residue;
		determine_residue_branching_info(
			rinfos_[ ii ].resid(), name3, same_chain_prev, same_chain_next, branch_points_on_this_residue, is_branch_point, is_branch_lower_terminus );

		bool const is_lower_terminus( ( ( ii == 1 || rinfos_.empty() || ( ! same_chain_prev && ! is_branch_lower_terminus) )
			&& check_Ntermini_for_this_chain ) );
		bool const is_upper_terminus( ( ii == nres_pdb || ! same_chain_next ) && check_Ctermini_for_this_chain );

		// Determine if this residue is a D-AA residue, an L-AA residue, or neither.
		StructFileRep::ResidueCoords const & xyz = rinfo.xyz();
		bool is_d_aa = NomenclatureManager::get_instance()->decide_is_d_aa( name3 );
		bool is_l_aa = NomenclatureManager::get_instance()->decide_is_l_aa( name3 );
		bool is_achiral = NomenclatureManager::get_instance()->decide_is_known_achiral( name3 );
		if ( ! (is_d_aa || is_l_aa || is_achiral ) ) {
			chemical::detect_ld_chirality_from_polymer_residue( xyz, name3, is_d_aa, is_l_aa );
		}

		// Get a list of ResidueTypes that could apply for this particular 3-letter PDB residue name.
		if ( ! is_residue_type_recognized( ii, name3 ) ) {
			residue_was_recognized_[ ii ] = false;
			continue;
		}

		TR.Debug << "Residue " << ii << std::endl;
		TR.Trace << "...same_chain_prev: " << same_chain_prev << std::endl;
		TR.Trace << "...same_chain_next: " << same_chain_next << std::endl;
		TR.Trace << "...is_lower_terminus: " << is_lower_terminus << std::endl;
		TR.Trace << "...check_Ntermini_for_this_chain: "<< check_Ntermini_for_this_chain << std::endl;
		TR.Trace << "...is_upper_terminus: " << is_upper_terminus << std::endl;
		TR.Trace << "...check_Ctermini_for_this_chain: "<< check_Ctermini_for_this_chain << std::endl;
		TR.Trace << "...is_branch_point: " << is_branch_point << std::endl;
		TR.Trace << "...is_branch_lower_terminus: "<< is_branch_lower_terminus << std::endl;
		TR.Trace << "...last_residue_was_recognized: " << last_residue_was_recognized( ii ) << std::endl;
		TR.Trace << "...is_d_aa: " << is_d_aa << std::endl;
		TR.Trace << "...is_l_aa: " << is_l_aa << std::endl;

		ResidueTypeCOP rsd_type_cop = get_rsd_type( name3, ii, branch_points_on_this_residue,
			resid, is_lower_terminus, is_upper_terminus, is_branch_point, is_branch_lower_terminus, is_d_aa, is_l_aa );

		if ( rsd_type_cop == 0 ) {
			std::string variant;
			if ( is_lower_terminus ) {
				variant += " lower-terminal";
			} else if ( is_branch_lower_terminus ) {
				variant += " branch-lower-terminal";
			}
			if ( is_upper_terminus ) {
				variant += " upper-terminal";
			}
			if ( is_branch_point ) {
				variant += " branch-point";
			}
			utility_exit_with_message( "No match found for unrecognized residue at position " +
				boost::lexical_cast< std::string >(ii) +
				"\nLooking for" + variant + " residue with 3-letter code: " + name3 );
		}

		residue_types_[ ii ] = rsd_type_cop;
		is_lower_terminus_[ ii ] = is_lower_terminus;
		same_chain_prev_[ ii ] = same_chain_prev;
		is_branch_lower_terminus_[ ii ] = is_branch_lower_terminus;
		fill_name_map( ii );
	}
}

/// @details Check for missing mainchain atoms; if there is no contiguous block of 3 mainchain
/// atoms that are not missing, then we will be unable to build an initial stub, or to fill in
/// the coordinates of those missing atoms, so, we will have to reject the residue.
void PoseFromSFRBuilder::pass_3_verify_sufficient_backbone_atoms()
{
	using namespace core::io::pdb;

	for ( core::Size ii = 1; ii <= rinfos_.size(); ++ii ) {
		if ( ! residue_types_[ ii ] ) continue;

		chemical::ResidueType const & iirestype( *residue_types_[ ii ] );
		if ( !iirestype.is_polymer() ) continue;

		StructFileRep::ResidueCoords const & iixyz( rinfos_[ ii ].xyz() );
		chemical::AtomIndices const & iimainchain( iirestype.mainchain_atoms() );
		NameBimap const & ii_atom_names = remapped_atom_names_[ ii ];
		Size const nbb( iimainchain.size() );
		if ( nbb >= 3 ) {
			bool mainchain_core_present( false );
			for ( Size jj=1; jj<= nbb-2; ++jj ) {
				std::string const & name1(iirestype.atom_name(iimainchain[jj  ]));
				std::string const & name2(iirestype.atom_name(iimainchain[jj+1]));
				std::string const & name3(iirestype.atom_name(iimainchain[jj+2]));
				if ( ! ii_atom_names.right.count(name1) ||
						!  ii_atom_names.right.count(name2) ||
						!  ii_atom_names.right.count(name3) ) {
					continue;
				}
				std::string const & rinfo_name1( ii_atom_names.right.find( name1 )->second );
				std::string const & rinfo_name2( ii_atom_names.right.find( name2 )->second );
				std::string const & rinfo_name3( ii_atom_names.right.find( name3 )->second );
				if ( iixyz.count( rinfo_name1 ) && iixyz.count( rinfo_name2 ) && iixyz.count( rinfo_name3 ) ) {
					mainchain_core_present = true;
					break;
				}
			}
			if ( !mainchain_core_present ) {
				TR.Warning << "[ WARNING ] skipping pdb residue b/c it's missing too many mainchain atoms: " <<
					rinfos_[ ii ].resid() << ' ' << rosetta_residue_name3s_[ ii ] << ' ' << iirestype.name() << std::endl;
				for ( Size jj=1; jj<= nbb; ++jj ) {
					std::string const & name(iirestype.atom_name(iimainchain[jj]));
					if ( ! ii_atom_names.right.count(name) ||
							!iixyz.count( ii_atom_names.right.find(name)->second ) ) {
						// Use of unmapped name deliberate
						TR << "missing: " << name << std::endl;
					}
				}
				if ( options_.exit_if_missing_heavy_atoms() == true ) {
					utility_exit_with_message("quitting due to missing heavy atoms");
				}
				// Designate that this residue should not be included in the Pose by setting the residue_type pointer to 0
				residue_types_[ ii ] = 0;
				residue_was_recognized_[ ii ] = false;
			}
		}
	}
}

/// @details Now take another pass over the residue types and for those residues that are now
/// the first residue of their chain (because the residue(s) before them had insufficient
/// backbone atoms) add the lower- or upper-termini variants to their residue types.
/// This behavior means that Rosetta will end up modeling chemistry that isn't there -- and
/// so this is more than a little questionable.  Alternatively, this code could and maybe
/// should add the truncation variants to the residue types.
void PoseFromSFRBuilder::pass_4_redo_termini()
{
	for ( core::Size ii = 1; ii <= rinfos_.size(); ++ii ) {
		if ( ! residue_types_[ ii ] ) { continue; }

		ResidueInformation const & rinfo = rinfos_[ ii ];
		char chainID = rinfo.chainID();

		if ( ! determine_check_Ntermini_for_this_chain( chainID ) ) { continue; }
		if ( ! determine_check_Ctermini_for_this_chain( chainID ) ) { continue; }

		Size const ii_prev = prev_residue_skipping_null_residue_types( ii );
		Size const ii_next = next_residue_skipping_null_residue_types( ii );

		bool type_changed( false );
		if ( ! residue_types_[ ii ]->is_polymer() ) { continue; }
		if ( ! residue_types_[ ii ]->is_lower_terminus() &&
				( ii == 1 || ii == ii_prev || ( residue_types_[ ii_prev ] && (
				!residue_types_[ ii_prev ]->is_polymer() ||
				(residue_types_[ ii_prev ]->is_upper_terminus() &&
				!residue_types_[ ii ]->is_branch_lower_terminus() ) ) ) ) ) {
			TR << "Adding undetected lower terminus type to residue " << ii << std::endl;
			residue_types_[ ii ] = residue_type_set_->get_residue_type_with_variant_added(
				*residue_types_[ ii ], chemical::LOWER_TERMINUS_VARIANT ).get_self_ptr();
			type_changed = true;
		}
		// AMW: new--don't add  an upper terminus type to protein residues following RNA
		// if that RNA is itself not upper-terminal
		if ( !residue_types_[ ii ]->is_upper_terminus() &&
				!( residue_types_[ ii ]->is_protein() && residue_types_[ ii_prev ]->is_RNA()
				&& !residue_types_[ ii_prev ]->is_upper_terminus() ) &&
				!( residue_types_[ ii ]->is_RNA() && residue_types_[ ii_next ]->is_protein()
				&& !residue_types_[ ii_next ]->is_upper_terminus() ) &&
				( ii == ii_next ||
				( residue_types_[ ii_next ] && (
				!residue_types_[ ii_next ]->is_polymer() ||
				residue_types_[  ii_next ]->is_lower_terminus() /* ||
				residue_types_[  ii_next ]->has_variant_type(BRANCH_LOWER_TERMINUS_VARIANT)*/ ) ) ) ) {
			TR << "Adding undetected upper terminus type to residue " << ii << std::endl;
			residue_types_[ ii ] = residue_type_set_->get_residue_type_with_variant_added(
				*residue_types_[ ii ], chemical::UPPER_TERMINUS_VARIANT ).get_self_ptr();
			type_changed = true;
		}

		if ( type_changed ) {
			fill_name_map( ii );
		}
	}
}

void PoseFromSFRBuilder::pass_5_note_discarded_atoms()
{
	for ( core::Size ii = 1; ii <= rinfos_.size(); ++ii ) {
		if ( ! residue_types_[ ii ] ) continue;

		debug_assert( residue_types_[ ii ]->natoms() >= remapped_atom_names_[ii].left.size() );
		core::Size missing_atoms( residue_types_[ ii ]->natoms() - remapped_atom_names_[ii].left.size() );
		if ( missing_atoms > 0 ) {
			TR.Debug << "Match: '" << residue_types_[ ii ]->name() << "'; missing " << missing_atoms << " coordinates" << std::endl;
		}

		debug_assert( rinfos_[ ii ].xyz().size() >= remapped_atom_names_[ii].left.size() );
		core::Size discarded_atoms( rinfos_[ ii ].xyz().size() - remapped_atom_names_[ ii ].left.size() );
		if ( is_lower_terminus_[ ii ] && rinfos_[ ii ].xyz().count(" H  ") && ! remapped_atom_names_[ ii ].left.count(" H  ") ) {
			// Don't worry about missing BB H if Nterm
			--discarded_atoms;
		}
		if ( discarded_atoms > 0 ) {
			TR.Warning << "[ WARNING ] discarding " << discarded_atoms
				<< " atoms at position " << ii << " in file " << sfr_.filename()
				<< ". Best match rsd_type:  " << residue_types_[ ii ]->name() << std::endl;
		}
	}
}


void PoseFromSFRBuilder::build_initial_pose( pose::Pose & pose )
{
	using namespace core::io::pdb;
	using namespace core::chemical;
	using namespace core::conformation;

	// TODO:
	// ERROR HANDLING HERE IF WE DID NOT RESOLVE ANY RESIDUE TYPES


	for ( Size ii = 1; ii <= residue_types_.size(); ++ii ) {
		if ( ! residue_types_[ ii ] ) { continue; }

		ResidueType const & ii_rsd_type( *residue_types_[ ii ] );
		TR.Trace << "ResidueType " << ii << ": " << ii_rsd_type.name() << std::endl;
		//TR.Trace << "went by, in the file context, " << rinfos_[ ii ].resName()
		// << " at " << rinfos_[ ii ].chainID() << rinfos_[ ii ].resSeq() << rinfos_[ ii ].iCode() << " " << rinfos_[ ii ].segmentID() << std::endl;

		ResidueOP ii_rsd( ResidueFactory::create_residue( ii_rsd_type ) );
		for ( StructFileRep::ResidueCoords::const_iterator iter = rinfos_[ ii ].xyz().begin(), iter_end = rinfos_[ ii ].xyz().end();
				iter != iter_end; ++iter ) {

			std::string const & rinfo_name( iter->first );
			if ( remapped_atom_names_[ ii ].left.count( rinfo_name ) ) {
				// offsetting all coordinates by a small constant prevents problems with atoms located
				// at position (0,0,0).
				// This is a bit of a dirty hack but it fixes the major problem of reading in rosetta
				// pdbs which usually start at 0,0,0. However the magnitude of this offset is so small
				// that the output pdbs should still match input pdbs. hopefully. yes. ahem.
				// RM: I'm not sure what the problem with having coordinates exactly at the origin is.
				// RM: If we do have a problem with that, it seems it should be fixed there and not here.
				// RM: (As you could imagine theoretically hitting (0,0,0) during minimization or packing.)
				double offset = 1e-250; // coordinates now double, so we can use _really_ small offset.
				std::string const & pose_name( remapped_atom_names_[ ii ].left.find( rinfo_name )->second );
				ii_rsd->atom( pose_name ).xyz( iter->second + offset );
				// +1 here as we haven't added the residue to the pose yet.
				id::NamedAtomID atom_id( pose_name, pose.size()+1 );
				coordinates_assigned_.set( atom_id, true);
			}
			//else runtime_assert( iter->first == " H  " && ii_rsd_type.is_terminus() ); // special casee
		}

		// These sister atoms are THUS FAR only found in RNA (it's to perfectly copy
		// native structure coordinates). So, at the moment, we only need to call
		// this function for RNA. Eventually, we may implement more proteinaceous use
		// cases, at which point you may want to move this condition into the
		// check_and_correct function.
		if ( ii_rsd->is_RNA() ) {
			check_and_correct_sister_atoms( ii_rsd );
		}

		Size const old_nres( pose.size() );

		TR.Trace << "...new residue is a polymer: " << ii_rsd->type().is_polymer() << std::endl;
		if ( old_nres >= 1 ) {
			TR.Trace << "...old residue is a polymer: " << pose.residue_type(old_nres).is_polymer() << std::endl;
		}

		// Add the first new residue to the pose.
		if ( !old_nres ) {
			TR.Trace << ii_rsd_type.name() << " " << ii << " is the start of a new pose." << std::endl;
			pose.append_residue_by_bond( *ii_rsd );

			// If this is a lower terminus, AND it's not about to get UPPER-UPPER bonded from prev
		} else if ( ( ( is_lower_terminus_[ ii ] && !( !pose.residue_type( old_nres ).is_upper_terminus() && pose.residue_type( old_nres ).is_RNA() && !ii_rsd_type.is_upper_terminus() ) )
				&& determine_check_Ntermini_for_this_chain( rinfos_[ ii ].chainID() )) ||
				! same_chain_prev_[ ii ] ||
				/* is_branch_lower_terminus || */
				pose.residue_type( old_nres ).has_variant_type( "C_METHYLAMIDATION" ) ||
				! ii_rsd->is_polymer() ||
				! pose.residue_type( old_nres ).is_polymer() ||
				! last_residue_was_recognized( ii ) ) {
			// A new chain because this is a lower terminus (see logic above for designation)
			// and if we're not checking it then it's a different chain from the previous

			core::Size root_index = 1;
			// connect metal ions by a jump to the closest metal-binding residue that is lower in sequence.
			if ( ii_rsd->is_metal() && basic::options::option[basic::options::OptionKeys::in::auto_setup_metals] ) {
				root_index = find_atom_tree_root_for_metal_ion( pose, ii_rsd );
			}
			if ( root_index>1 ) { TR << ii_rsd_type.name() << " " << ii << " was added by a jump, with base residue " << root_index << std::endl;}
			pose.append_residue_by_jump( *ii_rsd, root_index /*pose.size()*/ );

		} else { // Append residue to current chain dependent on bond length.
			if ( ! options_.missing_dens_as_jump() ) {
				TR.Trace << ii_rsd_type.name() << " " << ii << " is appended to chain " << rinfos_[ ii ].chainID() << std::endl;
				pose.append_residue_by_bond( *ii_rsd );
			} else {
				//fpd look for missing density in the input PDB
				//fpd if there is a bondlength > 3A
				//fpd we will consider this missing density
				// apl -- FIX THIS! Do not ask for the residue, as it forms an O(N^2) expense. Only ask for the residue type.
				Residue const & last_rsd( pose.residue( old_nres ) );
				core::Real bondlength = ( last_rsd.atom( last_rsd.upper_connect_atom() ).xyz() -
					ii_rsd->atom( ii_rsd->lower_connect_atom() ).xyz() ).length();

				if ( bondlength > 3.0 ) {
					TR << "[ WARNING ] missing density found at residue (rosetta number) " << old_nres << std::endl;
					pose.append_residue_by_jump( *ii_rsd, old_nres );

					if ( pose.residue_type(old_nres).is_protein() ) {
						if ( !pose.residue_type(old_nres).has_variant_type( UPPER_TERMINUS_VARIANT ) &&
								!pose.residue_type(old_nres).has_variant_type( UPPERTERM_TRUNC_VARIANT ) ) {
							core::pose::add_variant_type_to_pose_residue( pose, UPPERTERM_TRUNC_VARIANT, old_nres );
						}
					} else {
						if ( !pose.residue_type(old_nres).has_variant_type( UPPER_TERMINUS_VARIANT ) ) {
							core::pose::add_variant_type_to_pose_residue( pose, UPPER_TERMINUS_VARIANT, old_nres );
						}
					}

					if ( pose.residue_type(old_nres+1).is_protein() ) {
						core::pose::add_variant_type_to_pose_residue( pose, LOWERTERM_TRUNC_VARIANT, old_nres+1 );
					} else {
						core::pose::add_variant_type_to_pose_residue( pose, LOWER_TERMINUS_VARIANT, old_nres+1 );
					}
				} else {
					TR.Trace << ii_rsd_type.name() << " " << ii << " is appended to chain" << rinfos_[ ii ].chainID() << std::endl;
					pose.append_residue_by_bond( *ii_rsd );
				}
			}
		}

		// Report on residues with suspiciously bad atom coordinates.
		// (For now, this section just reports on bad rings, but other things could be added here to help users know
		// when an input PDB is probably bad.)
		Size const n_rings( ii_rsd_type.n_rings() );
		if ( n_rings ) {
			ii_rsd->update_nus();
		}
		for ( core::uint i( 1 ); i <= n_rings; ++i ) {
			if ( ii_rsd->ring_conformer( i ) != ii_rsd_type.ring_conformer_set( i )->get_lowest_energy_conformer() ) {
				TR.Warning << ii_rsd_type.name3() << ii << " has an unfavorable ring conformation; ";
				TR.Warning << "the coordinates for this input structure may have been poorly assigned. ";
				TR.Warning << "Hydrogen placement will likewise potentially be poor." << std::endl;
				TR.Debug << "  Measured: " << ii_rsd->ring_conformer( i ).specific_name << "  Expected: ";
				TR.Debug << ii_rsd_type.ring_conformer_set( i )->get_lowest_energy_conformer().specific_name << std::endl;

			}
		}

		// If newly added residue was a carbohydrate, set flag on conformation.
		if ( ii_rsd->is_carbohydrate() ) {
			pose.conformation().contains_carbohydrate_residues( true );
		}

		pose_to_rinfo_.push_back( ii );
		pose_temps_.push_back( rinfos_[ ii ].temps() );

		// Update the pose-internal chain label if necessary.
		if ( ( is_lower_terminus_[ ii ] || ! determine_check_Ntermini_for_this_chain( rinfos_[ ii ].chainID() ) ||
				is_branch_lower_terminus_[ ii ] ) && pose.size() > 1 ) {
			pose.conformation().insert_chain_ending( pose.size() - 1 );
		}
	}
}

void PoseFromSFRBuilder::refine_pose( pose::Pose & pose )
{
	// Sub-steps:
	// 1. Check termini and add terminal types if needed. Try to copy over any
	// additional matching atoms from the original Residue.
	// 2. Handle missing atoms.
	// 3. PDBInfo assembly
	// 3a. PDB-wide information (e.g. the file name)
	// 3b. Residue level information
	// 3c. Atom level information.
	// 4. Sanity check that DNA 5' phosphates are virtual
	// 5. Connectivity.

	using namespace core;
	using namespace pose::carbohydrates;
	using namespace io;
	using namespace pdb;
	using namespace chemical;
	using namespace conformation;


	//typedef std::map< std::string, double > ResidueTemps;

	// Note: _do not_ access Residue here. You do not need to. Doing so triggers
	// a refold and that is expensive. You probably can get away with
	// pose.residue_type( ii ). If you think you need Residue, ask someone else.

	// Step 2. Handle missing atoms.
	Size num_heavy_missing = 0;
	core::pose::initialize_atomid_map( missing_, pose ); // dimension the missing-atom mask

	// Poses with zero residues must exit here after establishing a zero
	// residue, all-atom-unrecognized PDBInfo.
	if ( pose.size() == 0 ) {
		// PDBInfo setup
		core::pose::PDBInfoOP pdb_info( new core::pose::PDBInfo( pose.size() ) );
		pdb_info->add_unrecognized_atoms( unrecognized_atoms_ );
		pose.pdb_info( pdb_info );
		return;
	}

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		Residue const & ii_rsd( pose.residue( ii ) );
		for ( Size jj = 1; jj <= ii_rsd.natoms(); ++jj ) {
			id::AtomID atom_id( jj, ii );
			id::NamedAtomID named_atom_id( ii_rsd.atom_name( jj ), ii );
			if ( ! coordinates_assigned_[ named_atom_id ] ) {
				missing_[ atom_id ] = true;
				if ( !ii_rsd.atom_is_hydrogen( jj ) ) num_heavy_missing++;
			}
		}
	}
	pose.conformation().fill_missing_atoms( missing_ );

	// Step 3. Build PDBInfo.
	build_pdb_info_1_everything_but_temps( pose );

	// Step 4. Sanity check that DNA 5' phosphates are virtual.

	// Most DNA structures lack 5' phosphate groups. 5' phosphates must be built
	// to serve as part of the backbone for atom/fold tree purposes, but they
	// must be made virtual so as not to affect physical calculations.
	for ( Size ii = 1, nres = pose.size(); ii <= nres; ++ii ) {
		Residue const & rsd = pose.residue( ii );
		if ( ! rsd.type().is_DNA() ) continue;

		for ( Size jj = 1, natoms = rsd.natoms(); jj <= natoms; ++jj ) {
			id::AtomID const id( jj, ii );
			if ( missing_[ id ] && rsd.atom_name( jj ) == " P  " ) {
				runtime_assert( rsd.has_variant_type( chemical::VIRTUAL_DNA_PHOSPHATE ) );
				break;
			}
		}
	}


	// Step 5. Reconcile connectivity data.

	// Look for and create any remaining non-mainchain (Edge::CHEMICAL) bonds
	// based on a specified radius from any unsatisfied residue connections.
	// This is used for such things as branched polymers, ubiquitination, or
	// covalent intermediates.
	// Note: The fold tree will remain with a jump between each such bond until
	// import_pose::set_reasonable_fold_tree() is called later, which actually
	// adds the CHEMICAL edges to fold tree; this method simply makes the bonds.
	pose.conformation().detect_bonds();

	// Also initialize other sidechain connectivity.
	core::pose::initialize_disulfide_bonds( pose, sfr_ );
	core::pose::ncbb::initialize_ncbbs( pose );

	// 1 residue fragments for ligand design.
	if ( pose.size() > 1 ) {
		pose.conformation().detect_pseudobonds();
	}

	// this is where pdb info used to get stored.
	build_pdb_info_2_temps( pose );

	//Create the GlycanTreeSet if carbohydrates are present.
	// This requires up-to-data connectivity information and is required to be in the pose if glycan residues are present.
	if ( pose.conformation().contains_carbohydrate_residues() ) {
		GlycanTreeSetObserverOP observer = GlycanTreeSetObserverOP( new GlycanTreeSetObserver( pose.conformation() ));
		pose.observer_cache().set( pose::datacache::GLYCAN_TREE_OBSERVER, observer, true /*auto_attach */ );
	}

	// Step 6. Addition of automatic constraints, bfactor repair, and comment reading.

	// Add constraints based on LINK records if desired
	if ( options_.constraints_from_link_records() ) {
		core::pose::get_constraints_from_link_records( pose, sfr_ );
	}

	// Fix the b factors of missing atoms using neighbors, and set H b factors
	// as 1.2x attached atom.
	if ( options_.preserve_crystinfo() ) {
		core::scoring::cryst::fix_bfactorsMissing( pose );
		core::scoring::cryst::fix_bfactorsH( pose );
	}

	if ( options_.pdb_comments() ) {
		std::map< std::string, std::string > const & pdb_comments( sfr_.pdb_comments() );
		for ( std::map< std::string, std::string >::const_iterator it = pdb_comments.begin(), end = pdb_comments.end();
				it != end; ++it ) {
			core::pose::add_comment( pose, it->first, it->second );
		}
	}
}

void
PoseFromSFRBuilder::build_pdb_info_1_everything_but_temps( pose::Pose & pose ) {

	using namespace core;
	using namespace io;
	using namespace pdb;
	using namespace chemical;
	using namespace conformation;

	// Step 3. PDBInfo assembly
	core::pose::PDBInfoOP pdb_info( new core::pose::PDBInfo( pose.size() ) );

	// Step 3a. Add PDB-wide information
	pdb_info->name( sfr_.filename() );
	pdb_info->modeltag( ( sfr_.modeltag() == "" ) ? sfr_.filename() : sfr_.modeltag() );
	if ( options_.preserve_header() ) {
		pdb_info->remarks( *( sfr_.remarks() ) );
		pdb_info->header_information( sfr_.header() );
	}

	// Step 3b. Collect residue level data.

	// Structure file residue indices can be negative: do not change to Size
	utility::vector1< int > pdb_numbering;
	utility::vector1< char > pdb_chains, insertion_codes;
	utility::vector1< std::string > segment_ids;
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		ResidueInformation const & rinfo = rinfos_[ pose_to_rinfo_[ ii ] ];
		pdb_numbering.push_back( rinfo.resSeq() );
		pdb_chains.push_back( rinfo.chainID() );
		insertion_codes.push_back( rinfo.iCode() );
		segment_ids.push_back( rinfo.segmentID() );
	}

	pdb_info->set_numbering( pdb_numbering );
	pdb_info->set_chains( pdb_chains );
	pdb_info->set_icodes( insertion_codes );
	pdb_info->set_segment_ids( segment_ids );
	if ( options_.preserve_crystinfo() ) {
		pdb_info->set_crystinfo( sfr_.crystinfo() );
	}

	// 3c. add atom level information.

	// Ensure enough space exists.
	pdb_info->resize_atom_records( pose );

	// Add unrecognized atoms.
	pdb_info->add_unrecognized_atoms( unrecognized_atoms_ );
	pose.pdb_info( pdb_info );

}

/// @details Now that the final residue types have been set for the Pose, it is safe to set the temperature data
/// in the PDBInfo.
void
PoseFromSFRBuilder::build_pdb_info_2_temps( pose::Pose & pose )
{

	using namespace core;
	using namespace pose;
	using namespace io;
	using namespace pdb;
	using namespace chemical;
	using namespace conformation;

	typedef std::map< std::string, double > ResidueTemps;

	PDBInfoOP pdb_info = pose.pdb_info();

	// Add temperatures.
	for ( core::Size ii = 1; ii <= pose.size(); ii++ ) {
		ResidueTemps const & res_temps( rinfos_[ pose_to_rinfo_[ ii ] ].temps() );
		NameBimap const & namemap( remapped_atom_names_[ pose_to_rinfo_[ ii ] ] );
		for ( ResidueTemps::const_iterator iter = res_temps.begin(); iter != res_temps.end(); ++iter ) {
			// namemap should only include atoms which have a presence in both rinfo and pose
			if ( namemap.left.count( iter->first ) ) {
				std::string const & pose_atom_name( namemap.left.find(iter->first)->second );
				if ( pose.residue( ii ).type().has( pose_atom_name ) ) { // There are issues with terminus patching which means atoms can sometimes disappear
					core::Size jj = pose.residue( ii ).type().atom_index( pose_atom_name );
					pdb_info->temperature( ii, jj, iter->second );
				}
			} else {
				if ( ( iter->first )[ 0 ] == 'H' || ( ( iter->first )[ 0 ] == ' ' && ( iter->first )[ 1 ] == 'H' ) ) {
					// Don't warn for missing temperatures on hydrogens.
				} else {
					TR << "[ WARNING ] can't find pose atom for file-residue " << ii << " atom " << iter->first << " (trying to store temperature in PDBInfo)" << std::endl;
				}
			}
		}
	}

	// Mark PDBInfo as okay and store in Pose.
	pdb_info->obsolete( false );
	// already stored -- pose.pdb_info( pdb_info );
}

// This function uses linkage information to determine main-chain and branch
// polymer connectivity.
/// @details  This function does two separate things that are related:
/// - It determines from linkage information whether a residue is a branch
///   lower terminus or branch point.
/// - It assigns main-chain connectivity to carbohydrate ResidueType base
///   names.
void
PoseFromSFRBuilder::determine_residue_branching_info(
	std::string const & resid,
	std::string const & name3,
	bool & same_chain_prev,
	bool & same_chain_next,
	utility::vector1< std::string > & branch_points_on_this_residue,
	bool & is_branch_point,
	bool & is_branch_lower_terminus )
{
	using namespace std;
	using namespace core::io::pdb;
	using namespace core::chemical;

	// Carbohydrate base names will have "->?)-" as a prefix if their main-
	// chain connectivity requires LINK records to determine.  Fortuitously,
	// position 2 is also the index for the atom name (such as " O2 ") that
	// provides the missing information.
	int const CARB_MAINCHAIN_CONN_POS = 2;

	bool unknown_main_chain_connectivity( false );
	if ( sfr_.residue_type_base_names().count( resid ) ) {
		TR.Trace << "Current residue has had its base name extracted from the PDB file: ";
		TR.Trace << sfr_.residue_type_base_names()[ resid ].second << endl;
		unknown_main_chain_connectivity =
			( sfr_.residue_type_base_names()[ resid ].second[ CARB_MAINCHAIN_CONN_POS ] == '?' );
	}

	//is_branch_point = false;
	TR.Trace << "Checking if resid " << resid << " is in the link map " << endl;
	if ( sfr_.link_map().count( resid ) ) {  // if found in the linkage map
		// The link map is keyed by resID of each branch point.
		TR.Trace << "Found resid " << resid << " in link map " << endl;
		Size const n_branches( sfr_.link_map()[ resid ].size() );
		for ( uint branch( 1 ); branch <= n_branches; ++branch ) {
			TR.Trace << "Examining branch " << branch << endl;
			LinkInformation const & link_info( sfr_.link_map()[ resid ][ branch ] );
			// LinkInformation contains both the branch point and branch lower
			// terminus of each linkage.
			if ( link_info.chainID1 == link_info.chainID2 && link_info.resSeq1 == ( link_info.resSeq2 - 1 )
					&& ! unknown_main_chain_connectivity ) {
				// If this occurs, the link is to the next residue on the same chain, so both residues are part of
				// the same main chain or branch, and this linkage information can be ignored, UNLESS this .pdb file
				// came from the PDB, in which case its 3-letter codes don't tell us the main chain, so we must get
				// the main chain from the LINK records.
				// Note that this all assumes insertion codes are not involved! It also assumes that the PDB file
				// makers did things reasonably.
				continue;
			}
			if ( unknown_main_chain_connectivity && branch == 1 ) {
				char const connectivity( link_info.name1[ CARB_MAINCHAIN_CONN_POS ] );
				TR.Trace << "Assigning main-chain connectivity to position " << connectivity;
				TR.Trace << " of this residue." << endl;
				sfr_.residue_type_base_names()[ resid ].second[ CARB_MAINCHAIN_CONN_POS ] = connectivity;
			} else {
				TR.Trace << "This residue is a branch point with corresponding branch lower terminus: ";
				TR.Trace << link_info.resID2 << endl;
				is_branch_point = true;
				branch_lower_termini_.push_back( link_info.resID2 );
				branch_points_on_this_residue.push_back( link_info.name1 );
			}
		}
	} else if ( unknown_main_chain_connectivity ) {
		// If .pdb 3-letter codes are being used, LINK records MUST be used to designate main-chain connectivity
		// for anything with a HETATM record.  So if we got here, it means that this must be an upper terminus.
		// First, assign an arbitrary main-chain connection and then reset this to be NOT of the same chain as
		// the next residue.
		TR.Trace << "Assigning main-chain connectivity arbitrarily to position 3 of this terminal residue." << endl;
		// TODO: In the future, we might want something else.
		sfr_.residue_type_base_names()[ resid ].second[ CARB_MAINCHAIN_CONN_POS ] = '3';
		same_chain_next = false;
	}
	is_branch_lower_terminus = branch_lower_termini_.contains( resid );
	if ( is_branch_lower_terminus && carbohydrates::CarbohydrateInfoManager::is_valid_sugar_code( name3 ) ) {
		same_chain_prev = false;
	}
}

/// @details Look at a list of potential residue types and finding the best one,
/// and record the set of unrecognized atoms.
bool
PoseFromSFRBuilder::is_residue_type_recognized(
	Size const pdb_residue_index,
	std::string const & rosetta_residue_name3,
	core::chemical::ResidueTypeCOPs const & rsd_type_list
){
	std::map< std::string, Vector > const & xyz( rinfos_[ pdb_residue_index ].xyz() );
	std::map< std::string, double > const & rtemp( rinfos_[ pdb_residue_index ].temps() );

	bool const is_HOH_to_ignore ( rosetta_residue_name3 == "HOH" && options_.ignore_waters() );

	if ( !rsd_type_list.empty() && !is_HOH_to_ignore ) {
		return true;
	}

	using namespace basic::options;
	if ( !(options_.ignore_unrecognized_res() ||
			options_.remember_unrecognized_res() ||
			is_HOH_to_ignore ) ) {
		// We should fail fast on unrecognized input rather than produce bad results!
		std::string message( "Unrecognized residue: " + rosetta_residue_name3 );
		if ( chemical::carbohydrates::CarbohydrateInfoManager::is_valid_sugar_code( rosetta_residue_name3 ) ) {
			message += "; did you intend to use the -include_sugars flag?";
		}
		utility_exit_with_message( message );
	}

	if ( !options_.remember_unrecognized_water() ) {
		// don't bother with water
		if ( rosetta_residue_name3 == "HOH" ) {
			return false;
		}
	}

	if ( options_.remember_unrecognized_res() ) {
		for ( std::map<std::string, Vector>::const_iterator iter=xyz.begin(), iter_end=xyz.end(); iter!= iter_end; ++iter ) {
			if ( unrecognized_atoms_.size() > 5000 ) {
				utility_exit_with_message("can't handle more than 5000 atoms worth of unknown residues\n");
			}
			TR << "remember unrecognized atom " << pdb_residue_index << " " << rosetta_residue_name3 << " " << stripped_whitespace(iter->first)
				<< " temp " << rtemp.find(iter->first)->second << std::endl;

			core::pose::UnrecognizedAtomRecord ua( pdb_residue_index, rosetta_residue_name3, stripped_whitespace(iter->first), iter->second, rtemp.find(iter->first)->second );

			unrecognized_atoms_.push_back( ua );
		}
	}

	if ( is_HOH_to_ignore ) output_ignore_water_warning_once();
	return false;
}

/// @brief Query the ResidueTypeSet using the residue-type-finder for a potential match based on the
/// name3 that has possibly been remapped by Rosetta through the Nomenclature manager.
bool
PoseFromSFRBuilder::is_residue_type_recognized(
	Size const pdb_residue_index,
	std::string const & rosetta_residue_name3
){
	using namespace core::chemical;

	// this residue list is only used to see if there are any residue_types with name3 at all:
	ResidueTypeCOPs rsd_type_list;
	ResidueTypeCOP  rsd_type = ResidueTypeFinder( *residue_type_set_ ).name3( rosetta_residue_name3 ).get_representative_type();
	if ( rsd_type != 0 ) rsd_type_list.push_back( rsd_type );
	return is_residue_type_recognized( pdb_residue_index, rosetta_residue_name3, rsd_type_list );
}



///////////////////////////////////////////////////////////////////////
// Use ResidueTypeFinder to efficiently figure out best match
//    residue_type to these PDB atom_names, name3, etc.
chemical::ResidueTypeCOP
PoseFromSFRBuilder::get_rsd_type(
	std::string const & rosetta_residue_name3,
	Size seqpos,
	utility::vector1< std::string > const & branch_points_on_this_residue,
	std::string const & resid,
	bool const is_lower_terminus,
	bool const is_upper_terminus,
	bool const is_branch_point,
	bool const is_branch_lower_terminus,
	bool const is_d_aa,
	bool const is_l_aa )
{
	// you can be neither but not both
	debug_assert( ! ( is_d_aa && is_l_aa ) );

	std::map< std::string, Vector > const & xyz( rinfos_[ seqpos ].xyz() );

	typedef std::map< std::string, Vector > ResidueCoords;
	using namespace core::chemical;
	using utility::tools::make_vector1;
	using utility::vector1;

	vector1< vector1< VariantType > > required_variants_in_sets;
	vector1< ResidueProperty > properties, disallow_properties;
	vector1< VariantType > disallow_variants;  // are variants different from properties?
	std::string residue_base_name( "" );  // used when we have more information then just a 3-letter code
	if ( sfr_.residue_type_base_names().count( resid ) ) {
		residue_base_name = sfr_.residue_type_base_names()[ resid ].second;
	}
	vector1< std::string > patch_names;
	bool disallow_carboxyl_conjugation_at_glu_asp( false );

	ResidueTypeCOP rsd_type = ResidueTypeFinder( *residue_type_set_ ).name3( rosetta_residue_name3 ).get_representative_type();
	if ( rsd_type->is_polymer() ) {
		if ( is_lower_terminus ) {
			required_variants_in_sets.push_back(
				make_vector1( LOWER_TERMINUS_VARIANT, LOWERTERM_TRUNC_VARIANT, METHYL_GLYCOSIDE ) );
		} else {
			disallow_variants.push_back( LOWER_TERMINUS_VARIANT );
			disallow_variants.push_back( LOWERTERM_TRUNC_VARIANT );
			disallow_variants.push_back( METHYL_GLYCOSIDE );
		}
		if ( is_upper_terminus ) {
			required_variants_in_sets.push_back( make_vector1( UPPER_TERMINUS_VARIANT, UPPERTERM_TRUNC_VARIANT ) );
		} else {
			disallow_variants.push_back( UPPER_TERMINUS_VARIANT );
			disallow_variants.push_back( UPPERTERM_TRUNC_VARIANT );
		}
		// note that carbohydrate branch points will be covered by patch_names below.
		if ( ! rsd_type->is_carbohydrate() ) {
			if ( is_branch_point ) {
				properties.push_back( BRANCH_POINT );
			} else {
				disallow_properties.push_back( BRANCH_POINT );
			}
		}
		if ( is_branch_lower_terminus ) {
			properties.push_back( BRANCH_LOWER_TERMINUS );
		} else {
			disallow_properties.push_back( BRANCH_LOWER_TERMINUS );
			disallow_variants.push_back( BRANCH_LOWER_TERMINUS_VARIANT );
		}
		if ( is_d_aa ) {
			properties.push_back( D_AA );
		}
		if ( is_l_aa ) {
			properties.push_back( L_AA );
		}
	}
	if ( rsd_type->aa() == aa_cys && rosetta_residue_name3 != "CYD" ) {
		disallow_variants.push_back( DISULFIDE );
	}
	if ( !is_branch_lower_terminus ) {
		disallow_carboxyl_conjugation_at_glu_asp = true;
	}

	if ( !options_.keep_input_protonation_state() ) {
		disallow_variants.push_back( PROTONATED );
		disallow_variants.push_back( DEPROTONATED );
	}
	if ( rsd_type->is_carbohydrate() ) {
		// The below assumes that ResidueTypes with fewer patches are selected 1st, that is, that an
		// :->2-branch ResidueType will be checked as a possible match before an :->2-branch:->6-branch
		// ResidueType.  If this were not the case, Rosetta could misassign an :->2-branch:->6-branch
		// ResidueType to a residue that actually only has a single branch at the 2 or 6 position.
		for ( core::uint k( 1 ); k <= branch_points_on_this_residue.size(); ++k ) {
			char const & branch_point = branch_points_on_this_residue[ k ][ 2 ];  // 3rd column (index 2) is the atom number.
			patch_names.push_back( "->" + std::string( 1, branch_point ) + ")-branch" );
		}
	}

	utility::vector1< std::string > xyz_atom_names;
	for ( ResidueCoords::const_iterator iter=xyz.begin(), iter_end=xyz.end(); iter!= iter_end; ++iter ) {
		std::string xyz_name = iter->first;
		xyz_atom_names.push_back( xyz_name );
	}

	rsd_type = ResidueTypeFinder( *residue_type_set_ )
		.name3( rosetta_residue_name3 )
		.residue_base_name( residue_base_name )
		.disallow_variants( disallow_variants )
		.variants_in_sets( required_variants_in_sets )
		.properties( properties )
		.disallow_properties( disallow_properties )
		.patch_names( patch_names )
		.ignore_atom_named_H( is_lower_terminus )
		.disallow_carboxyl_conjugation_at_glu_asp( disallow_carboxyl_conjugation_at_glu_asp )
		.check_nucleic_acid_virtual_phosphates( true )
		.get_best_match_residue_type_for_atom_names( xyz_atom_names );

	return rsd_type;
}

/// @brief Returns the input resid if it's either the first residue or it's
/// the second residue and the first residue has a mrb_merge_w_next behavior specified
Size PoseFromSFRBuilder::prev_residue_skipping_merges( Size resid ) const
{
	return resid > 1 ? ( merge_behaviors_[ resid-1 ] != core::chemical::mrb_merge_w_next ? resid-1 : ( resid-1 > 1 ? resid-2 : resid )) : resid;
}

/// @brief Returns the input resid if it's the first residue that has a non-null-pointing
/// entry in the residue_types_ array.
Size PoseFromSFRBuilder::prev_residue_skipping_null_residue_types( Size resid ) const
{
	if ( resid == 1 ) return resid;
	for ( Size ii = resid-1; ii >= 1; --ii ) {
		if ( residue_types_[ ii ] ) return ii;
	}
	return resid;
}

/// @brief Returns the input resid if it's either the last residue or it's
/// the second-to-last residue and the last residue has a mrb_merge_w_prev behavior specified
Size PoseFromSFRBuilder::next_residue_skipping_merges( Size resid ) const
{
	//return resid < rinfos_.size() ? ( merge_behaviors_[ resid+1 ] != core::chemical::mrb_merge_w_prev ? resid+1 : ( resid+1 < rinfos_.size() ? resid+2 : resid  )) : resid;
	if ( resid < rinfos_.size()  ) {
		if ( merge_behaviors_[ resid+1 ] != core::chemical::mrb_merge_w_prev ) {
			return resid+1;
		} else if ( resid+1 < rinfos_.size() ) {
			return resid+2;
		}
	}
	return resid;
}

/// @brief Returns the input resid if it's the last residue that has a non-null-pointing
/// entry in the residue_types_ array.
Size PoseFromSFRBuilder::next_residue_skipping_null_residue_types( Size resid ) const
{
	if ( resid == residue_types_.size() ) return resid;
	for ( Size ii = resid+1; ii <= residue_types_.size(); ++ii ) {
		if ( residue_types_[ ii ] ) return ii;
	}
	return resid;
}

bool PoseFromSFRBuilder::determine_separate_chemical_entity( char chainID ) const
{
	std::string const & sep_chains = options_.chains_whose_residues_are_separate_chemical_entities();
	return std::find( sep_chains.begin(), sep_chains.end(), chainID ) != sep_chains.end();
}

bool PoseFromSFRBuilder::determine_same_chain_prev( Size resid, bool separate_chemical_entity ) const
{
	Size res_prev = prev_residue_skipping_merges( resid );
	bool res_prev_unrecognized = ! residue_was_recognized_[ res_prev ];

	//Check to see if the unrecognized residue was at an Nterminus and should have started the new chain.
	if ( res_prev_unrecognized && res_prev != 1 ) {
		Size res_prev2 = prev_residue_skipping_merges( res_prev );
		bool res_prev_separate_chemical_entity =determine_separate_chemical_entity( rinfos_[ res_prev].chainID() );

		bool prev_res_same_chain = (resid != res_prev && rinfos_[ res_prev ].chainID() == rinfos_[ res_prev2 ].chainID() &&
			rinfos_[ res_prev ].terCount() == rinfos_[ res_prev2 ].terCount() && ! res_prev_separate_chemical_entity );

		if ( ! prev_res_same_chain ) {
			return false;
		}
	}

	return resid != res_prev && rinfos_[ resid ].chainID() == rinfos_[ res_prev ].chainID() &&
		rinfos_[ resid ].terCount() == rinfos_[ res_prev ].terCount() && ! separate_chemical_entity;
}

bool PoseFromSFRBuilder::determine_same_chain_next( Size resid, bool separate_chemical_entity ) const
{
	Size res_next = next_residue_skipping_merges( resid );
	return resid != res_next && rinfos_[resid].chainID() == rinfos_[ res_next ].chainID() &&
		rinfos_[resid].terCount() == rinfos_[res_next].terCount() && ! separate_chemical_entity;
}

bool PoseFromSFRBuilder::determine_check_Ntermini_for_this_chain( char chainID ) const
{
	std::string const & nterm_chains = options_.check_if_residues_are_Ntermini();
	return "ALL" == nterm_chains ||
		std::find( nterm_chains.begin(), nterm_chains.end(), chainID ) == nterm_chains.end();
}

bool PoseFromSFRBuilder::determine_check_Ctermini_for_this_chain( char chainID ) const
{
	std::string const & cterm_chains = options_.check_if_residues_are_Ctermini();
	return "ALL" == cterm_chains ||
		std::find( cterm_chains.begin(), cterm_chains.end(), chainID ) == cterm_chains.end();
}


/// @details modifies the remapped_atom_names_ bimap at position seqpos given the names that
/// are in residue_types_[ seqpos ]
void PoseFromSFRBuilder::fill_name_map( Size seqpos )
{
	assert( residue_types_[ seqpos ] );

	NameBimap & name_map( remapped_atom_names_[ seqpos ] );
	core::io::ResidueInformation const & rinfo( rinfos_[ seqpos ] );
	chemical::ResidueType const & rsd_type( *residue_types_[ seqpos ] );

	//Reset name map
	bool rename = rsd_type.remap_pdb_atom_names();
	for ( core::Size ii(1); ii <= options_.residues_for_atom_name_remapping().size(); ++ii ) {
		rename = rename || ( options_.residues_for_atom_name_remapping()[ ii ] == rsd_type.name3() );
	}
	if ( rename ) {
		// Remap names according to bonding pattern and elements
		// Will reset whatever is in name_map
		//AMW: move from file_data_fixup
		remap_names_on_geometry( name_map, rinfo, rsd_type );
	} else {
		name_map.clear();
		// Using names as-is, except for canonicalizing
		for ( utility::vector1< AtomInformation >::const_iterator
				iter=rinfo.atoms().begin(), iter_end=rinfo.atoms().end(); iter!= iter_end; ++iter ) {
			std::string const & name ( iter->name );
			if ( ! rinfo.xyz().count( name ) ) { // Only map atoms with coordinates.
				continue;
			}
			std::string strip_name( stripped_whitespace(name) );
			if ( rsd_type.has( strip_name ) ) {
				// We do the diversion through the index to canonicalize the atom name spacing to Rosetta standards.
				std::string canonical( rsd_type.atom_name( rsd_type.atom_index(strip_name) ) );
				name_map.insert( NameBimap::value_type( name, canonical ) );
			}
		}
	}
}

Size
PoseFromSFRBuilder::find_atom_tree_root_for_metal_ion(
	pose::Pose const & pose,
	conformation::ResidueCOP metal_rsd
)
{

	// If this is a metal ion and we're automatically setting up metals, search for the closest metal-binding residue
	// and make that the jump parent.  Otherwise, let the jump parent be the closest residue.
	numeric::xyzVector< core::Real > const metal_xyz = metal_rsd->xyz(1); //Atom 1 is always the metal of a residue representing a metal ion.  (There's a check for this in residue_io.cc).

	core::Size closest_metalbinding_residue=0;
	core::Size metalbinding_dist_sq = 0;
	core::Size closest_residue=0;
	core::Size closest_dist_sq = 0;

	for ( core::Size ii=1, nres=pose.size(); ii<=nres; ++ii ) { //Loop through all residues already added, looking for possible residues to root the metal onto.
		if ( !pose.residue(ii).is_protein() ) continue; //I'm not interested in tethering metals to non-protein residues.
		if ( !pose.residue(ii).has("CA") ) continue; //I'll be basing this on metal-alpha carbon distance, so anything without an alpha carbon won't get to be the root.

		numeric::xyzVector < core::Real > const residue_xyz = pose.residue(ii).xyz("CA");

		core::Real const current_dist_sq = residue_xyz.distance_squared(metal_xyz);

		if ( closest_residue==0 || current_dist_sq < closest_dist_sq ) {
			closest_residue = ii;
			closest_dist_sq = (core::Size)current_dist_sq;
		}
		if ( pose.residue(ii).is_metalbinding() &&
				(closest_metalbinding_residue==0 || current_dist_sq < metalbinding_dist_sq) ) {
			closest_metalbinding_residue = ii;
			metalbinding_dist_sq = (core::Size) current_dist_sq;
		}
	} //Inner loop through all residues

	if ( closest_metalbinding_residue!=0 ) return closest_metalbinding_residue; //If we found a metal-binding residue, it's the root; otherwise, the closest residue is.
	else if ( closest_residue!=0 ) return closest_residue;

	return 1;
}

bool PoseFromSFRBuilder::last_residue_was_recognized( Size seqpos ) const
{
	return residue_was_recognized_[ prev_residue_skipping_merges( seqpos )];
}

/// @brief Tell user about ancient Rosetta choice to ignore HOH waters if -ignore_unrecognized_res is set.
void
PoseFromSFRBuilder::output_ignore_water_warning_once()
{
	if ( outputted_ignored_water_warning_ ) return;
	if ( options_.ignore_unrecognized_res() && options_.ignore_waters() ) {
		TR << TR.Red << "For backwards compatibility, setting -ignore_unrecognized_res leads ALSO to -ignore_waters. You can set -ignore_waters false to get waters." << TR.Reset << std::endl;
		outputted_ignored_water_warning_ = true;
	}
}



/// for nucleic acids, slightly better mechanism is below (convert_nucleic_acid_residue_info_to_standard)
std::string
convert_atom_name( std::string const & res_name, std::string const & atom_name )
{
	if ( res_name == "MET" && ( atom_name == " S  " || atom_name == "S" ) ) {
		// " S  " for PDB format, "S" for mmCIF
		return " SD ";
	} else if ( res_name == "MET" && ( atom_name == "SE  " || atom_name == "SE" ) ) {
		// " SE " for PDB format, "SE" for mmCIF
		TR << "Reading Selenium SE from MSE as SD from MET" << std::endl;
		return " SD ";
	}
	return atom_name;
}

/// @details  Temporary hacky hack
/// Need better mechanism for this
/// for nucleic acids, slightly better mechanism is below (convert_nucleic_acid_residue_info_to_standard)
std::string
convert_res_name( std::string const & name )
{
	if ( name == "MSE" ) {
		TR << "Reading MSE as MET!" << std::endl;
		return "MET";
	}

	//If this is one of these metal ions and there is a 1 or 2 appended to the name (e.g. "CU2"), return just the first two characters as the name.
	// AMW: this is eventually a job for the nomenclature manager but for now I will settle for "not broken"
	if ( name == "CU1" ) return "CU ";
	if ( name == "ZN2" ) return "ZN ";

	return name;
}


void
create_working_data(
	StructFileRepOptions const & options,
	StructFileRep const & sfr,
	utility::vector1< core::io::ResidueInformation > & rinfos )
{
	rinfos.clear();

	for ( Size ch=0; ch < sfr.chains().size(); ++ch ) {
		for ( Size i=0; i < sfr.chains()[ch].size(); ++i ) {
			AtomInformation ai( sfr.chains()[ch][i] );
			// we should make a copy instead of taking a reference if "fixing" the names causes problems
			std::string const  res_name( convert_res_name( ai.resName ) );
			std::string const atom_name( convert_atom_name( res_name, ai.name ) );
			ai.resName = res_name;
			ai.name = atom_name;

			bool const ok = update_atom_information_based_on_occupancy( options, ai );
			if ( !ok ) continue;

			ResidueInformation new_res( ai );
			if ( rinfos.size() == 0 || rinfos.back() != new_res ) rinfos.push_back(new_res);
			rinfos.back().append_atom( ai );
		}
	}
}


// what to do if occupancy is 0.0?
//chu modify the logic how atoms are treated with zero or negative occupancy field.
bool
update_atom_information_based_on_occupancy(
	StructFileRepOptions const & options,
	AtomInformation & ai )
{

	if ( ai.occupancy == 0.0 ) {
		if ( options.randomize_missing_coords() ) {
			randomize_missing_coords( ai );
		} else if ( !options.ignore_zero_occupancy() ) {
			// do nothing and keep this atom as it is
		} else {
			//When flag default changes from true to false, change to TR.Debug and remove second line
			TR.Warning << "PDB reader is ignoring atom " << ai.name << " in residue " << ai.resSeq << ai.iCode << ai.chainID
				<< ".  Pass flag -ignore_zero_occupancy false to change this behavior" << std::endl;
			return false; // skip this atom with zero occ by default
		}
	} else if ( ai.occupancy < 0.0 ) { // always randomize coords for atoms with negative occ
		randomize_missing_coords( ai );
	} else {
		// do nothing for normal atoms with positive occ
	}
	return true;
}


/// @details The missing density regions in the input pdb should have 0.000 in the placeholders
/// this routine puts random coordinates wherever there is 0.000 for mainchain atoms.
/// tex - that's a stupid way of defining missing density, as atoms can be at the origin for other
/// reasons. This has been updated to check for occupancy to define missing density rather than atoms
/// located at the origin.
void
randomize_missing_coords( AtomInformation & ai ) {
	// if( ai.resSeq == 1 && ai.name == " N  ") return;//ignore first atom. Rosetta pdbs start with 0.000
	if ( ai.x == 0.000 && ai.y == 0.000 && ai.z == 0.000 && ai.occupancy <= 0.0 ) {
		TR << "Randomized: " << ai.name << " " << ai.resName << "  " << ai.resSeq << std::endl;
		//v  if ( ai.name == " N  " || ai.name == " CA " || ai.name == " C  " ||
		//v   ai.name == " O  " || ai.name == " CB " ) {
		ai.x = ai.x + 900.000 + numeric::random::rg().uniform()*100.000;
		ai.y = ai.y + 900.000 + numeric::random::rg().uniform()*100.000;
		ai.z = ai.z + 900.000 + numeric::random::rg().uniform()*100.000;
		//v  }
	}
	return;
}

} // namespace pose_from_sfr
} // namespace io
} // namespace core
