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
#include <core/io/NomenclatureManager.hh>
#include <core/io/ResidueInformation.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/io/util.hh>
#include <core/chemical/io/merge_and_split_behaviors_io.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedAtomID_Map.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/PoseResidueTypeSet.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>
#include <core/chemical/rings/RingConformerSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/ncbb/util.hh>
#include <core/pose/init_id_map.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/cryst/util.hh>
#include <core/simple_metrics/util.hh>
#include <core/simple_metrics/SimpleMetricData.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/random/random.hh>

// Utility headers
#include <utility>
#include <utility/string_constants.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh> // for nmers_of
#include <utility/string_util.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>

// External headers

// C++ headers
#include <sstream>
#include <utility>
#include <algorithm>    // std::sort std::find
#include <vector>


namespace core {
namespace io {
namespace pose_from_sfr {

// Tracer instance for this file
static basic::Tracer TR( "core.io.pose_from_sfr.PoseFromSFRBuilder" );

using namespace ObjexxFCL;
using utility::nmers_of;


PoseFromSFRBuilder::PoseFromSFRBuilder( chemical::ResidueTypeSetCOP rts, StructFileRepOptions const & options ) :
	residue_type_set_(std::move( rts )),
	options_( options ),
	coordinates_assigned_( false ),
	outputted_ignored_water_warning_( false ),
	nfix_(0),
	showed_nfix_warning_(false)
{}

PoseFromSFRBuilder::~PoseFromSFRBuilder() = default;

/// @details The process of building a Pose occurs in several phases.  First, in the setup() function,
/// the builder converts the SFR into a vector of ResidueInformation objects -- renaming some atoms along
/// the way. The next five phases involve five passes over this ResidueInformation vector.
void
PoseFromSFRBuilder::build_pose( StructFileRep const & sfr, pose::Pose & pose )
{
	pose.clear();
	/// Clearing cleared the RTS in the pose.
	/// If our ResidueTypeSet is a PoseRTS, we should stash a copy in the pose itself.
	auto pose_rts = utility::pointer::dynamic_pointer_cast< core::chemical::PoseResidueTypeSet const >( residue_type_set_ );
	if ( pose_rts != nullptr ) {
		pose.conformation().reset_residue_type_set_for_conf( pose_rts );
	}

	setup( sfr );
	pass_1_split_and_merge_residues_as_necessary();
	if ( options_.fast_restyping() ) {
		pass_2_quick_and_dirty_resolve_residue_types();
	} else {
		pass_2_resolve_residue_types();
	}
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
	static Size const max_fix( 2 );
	static bool const show_all_fixup( options_.show_all_fixes() );

	for ( Size ii = 1; ii <= rinfos_.size(); ++ii ) {
		core::io::ResidueInformation & rinfo  = rinfos_[ ii ];
		std::string const original_name = rinfo.resName();

		// first establish if this is DNA or RNA (or something else)
		if ( !options_.guarantee_no_DNA()
				&& NomenclatureManager::is_old_DNA( rinfo.resName() )
				&& missing_O2prime( rinfo.atoms() ) )  {
			std::string new_name( original_name ); new_name.replace( 1, 1, "D" ); // A --> dA
			rinfo.resName( new_name );
			if ( ++nfix_ <= max_fix || show_all_fixup ) {
				TR << "Converting residue name " <<  original_name <<  " to " << rinfo.resName() << std::endl;
			}
		}

		if ( NomenclatureManager::is_old_RNA( rinfo.resName() ) ) {
			std::string new_name( original_name ); new_name.replace( 1, 1, " " ); // rA --> A
			rinfo.resName( new_name );
			if ( ++nfix_ <= max_fix || show_all_fixup ) {
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
				if ( ++nfix_ <= max_fix || show_all_fixup ) {
					TR << "Converting atom name    " << original_atom_name << " to " << new_atom_name << std::endl;
				}
			}
		}
	}

	if ( nfix_ > max_fix && !show_all_fixup && !showed_nfix_warning_ ) {
		TR << "Number of nucleic acid residue fixups exceeds output limit. ";
		TR << "Rerun with -show_all_fixes to show everything." << std::endl;
		showed_nfix_warning_ = true;
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
	using LinkInformationVect = utility::vector1<LinkInformation>;
	typedef std::map< std::string, utility::vector1< LinkInformation > > LinkMap;

	// AMW: In CIF nomenclature, LINKs to be handled:
	// This is actually a "close contact" but it's O3' of adenine to C of FME, UPPER-UPPER
	// 1  1 "O3'" v A   76   ? ? C     v FME 77   ? ? 1.56

	LinkMap pruned_links;
	for ( const auto & iter : sfr.link_map() ) {
		std::string const & resID1 = iter.first;
		LinkInformationVect const & links = iter.second;
		for ( const auto & link : links ) {
			bool const link_has_metal( NomenclatureManager::is_metal( link.resName1 ) ||
				NomenclatureManager::is_metal( link.resName2 ) );
			bool const link_has_sugar( NomenclatureManager::is_sugar( link.resName1 ) ||
				NomenclatureManager::is_sugar( link.resName2 ) );

			if ( link.resName1 == "HOH" || link.resName2 == "HOH" ) {
				TR.Debug << "Omitting LINK record that implies a covalent bond to a nearby ";
				TR.Debug << "molecule of solvent. This is inappropriate even when solvent is ";
				TR.Debug << "read in and doubly so when it isn't." << std::endl;
			} else if ( link_has_metal ) {
				TR.Debug << "Omitting LINK record that uses a metal. These will be processed ";
				TR.Debug << "by -in:auto_setup_metals." << std::endl;
			} else if ( options_.ignore_sugars() && link_has_sugar ) {
				TR.Debug << "Omitting LINK record that uses a saccharide residue. ";
				TR.Debug << "Did you mean to use the -include_sugars flag?" << std::endl;
			} else if ( ( link.resSeq1 == link.resSeq2 - 1
					&& utility::strip(link.name1) == "O3'" && utility::strip(link.name2) == "P" )
					||  ( link.resSeq1 == link.resSeq2 + 1
					&& utility::strip(link.name1) == "P" && utility::strip(link.name2) == "O3'" ) ) {
				// Strip each string first to handle both PDB and CIF versions.
				// We have a normal polymeric connection written as a LINK.
				TR.Debug << "Omitting LINK record that represents the canonical polymeric connectivity of a NCNT." << std::endl;
			} else if ( ( link.resSeq1 == link.resSeq2 - 1
					&& utility::strip(link.name1) == "C" && utility::strip(link.name2) == "N" )
					||  ( link.resSeq1 == link.resSeq2 + 1
					&& utility::strip(link.name1) == "N" && utility::strip(link.name2) == "C" ) ) {
				// We have a normal polymeric connection written as a LINK.
				TR.Debug << "Omitting LINK record that represents the canonical polymeric connectivity of a NCAA." << std::endl;
			} else if ( link.resName1 == "CYS" && link.resName2 == "CYS" && utility::strip(link.name1) == "SG" && utility::strip(link.name2) == "SG" ) {
				// We have an SSBOND redundantly specified as a LINK.
				TR.Debug << "Omitting LINK record that gives a SECOND specification of a disulfide bond." << std::endl;
			} else if ( ! link_has_sugar && (  // PDB sugar codes will not be recognized by the ResidueTypeFinder!
					! ResidueTypeFinder( *residue_type_set_ ).name3( link.resName1 ).get_representative_type() ||
					! ResidueTypeFinder( *residue_type_set_ ).name3( link.resName2 ).get_representative_type() ) ) {
				// ResidueTypeFinder is a bit slow, so save this for the last check
				// One or more residues in this LINK is not recognized.  Move on!
				TR.Debug << "Omitting LINK record that uses an unrecognized residue." << std::endl;
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
PoseFromSFRBuilder::pass_1_split_and_merge_residues_as_necessary()
{
	using namespace core::io::pdb;
	using namespace core::chemical;
	using namespace core::chemical::io;

	create_working_data( options_, sfr_, rinfos_ );
	convert_nucleic_acid_residue_info_to_standard();

	// Split residues.
	{
		core::uint i( 1 );
		while ( i <= rinfos_.size() ) {
			SplitBehaviors const & residues_renamings_pair =
				residue_type_set_->merge_split_behavior_manager().split_behavior_for_name3( rinfos_[ i ].resName() );
			if ( ! residues_renamings_pair.first.empty() ) {
				ResidueInformation const rinfo( rinfos_[ i ] ); // make a copy before replacing.
				Size const n_res_into_which_to_split( residues_renamings_pair.first.size() );
				TR << "Splitting residue " << rinfo.resName() << " into " << n_res_into_which_to_split << " residues: ";
				rinfos_[ i ] = ResidueInformation();  // Clear the original residue.
				// Make space in the vector for the new residues.
				rinfos_.insert( rinfos_.begin() + i, n_res_into_which_to_split - 1, ResidueInformation() );
				for ( core::uint j( 0 ); j < n_res_into_which_to_split; ++j ) {
					ResidueInformation & res_frag( rinfos_[ i + j ] );
					std::string const & new_name( residues_renamings_pair.first[ j + 1 ].first );
					if ( j ) { TR << ", "; }
					TR << new_name;

					// Create the new residue.
					// Keep the same chain and sequence ID, but give it an insertion code.
					res_frag.resName( new_name );
					res_frag.chainID( rinfo.chainID() );
					res_frag.resSeq( rinfo.resSeq() );
					res_frag.iCode( utility::UPPERCASE_LETTERS[ j ] );
					res_frag.terCount( rinfo.terCount() );
					res_frag.segmentID( rinfo.segmentID() );

					sfr_.residue_type_base_names()[ res_frag.resid() ] = residues_renamings_pair.first[ j + 1 ];

					// Copy and rename the atoms.
					AtomRenamingMap const & renamings( residues_renamings_pair.second[ j + 1 ] );
					Size const n_atoms( rinfo.atoms().size() );
					for ( core::uint k( 1 ); k <= n_atoms; ++k ) {
						std::string const & old_name( rinfo.atoms()[ k ].name );
						if ( renamings.count( old_name ) ) {
							AtomInformation new_atom = rinfo.atoms()[ k ];
							new_atom.name = renamings.at( old_name );
							res_frag.append_atom( new_atom );
						}
					}
				}
				TR << std::endl;

				i += n_res_into_which_to_split;
			} else {
				++i;
			}
		}
	}


	residue_types_.resize(              rinfos_.size() );
	is_lower_terminus_.resize(          rinfos_.size(), false );
	same_chain_prev_.resize(            rinfos_.size(), true );
	residue_was_recognized_.resize(     rinfos_.size(), true );
	remapped_atom_names_.resize(        rinfos_.size() );
	merge_behaviors_.resize(            rinfos_.size() );
	merge_atom_maps_.resize(            rinfos_.size() );


	// Merge residues.
	for ( Size ii = 1; ii <= rinfos_.size(); ++ii ) {
		resid_to_index_[ rinfos_[ ii ].resid() ] = ii;
		ResidueMergeInstructions const & behavior_map_pair =
			residue_type_set_->merge_split_behavior_manager().merge_behavior_for_name3( rinfos_[ ii ].resName() );

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
				utility_exit_with_message( "The last residue, residue" + utility::to_string( ii ) +
					" named \"" +rinfos_[ ii ].resName() +
					"\" has been indicated to merge with the next residue from the core::io::NomenclatureManager" );
			}

			ResidueInformation & rmerged_into = rinfos_[ merge_behaviors_[ ii ] == mrb_merge_w_prev ? ii-1 : ii+ 1 ];

			AtomRenamingMap nowhitespace_rename_map;
			for ( AtomRenamingMap::const_iterator
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

	pre_process_residue_type_information();

	for ( Size ii = 1; ii <= rinfos_.size(); ++ii ) {
		if ( merge_behaviors_[ ii ] == chemical::io::mrb_merge_w_next ||
				merge_behaviors_[ ii ] == chemical::io::mrb_merge_w_prev ) {
			continue;
		}

		core::io::ResidueInformation const & rinfo = rinfos_[ ii ];
		char chainID = rinfo.chainID();
		std::string const & name3 = rinfo.rosetta_resName();
		std::string const & resid = rinfo.resid();

		utility::vector1< std::string > known_connect_atoms_on_this_residue;

		//----- the hairy stuff happens here -----------------------
		determine_residue_branching_info( ii, known_connect_atoms_on_this_residue, known_links_ );
		//----------------------------------------------------------

		bool const separate_chemical_entity = determine_separate_chemical_entity( chainID );
		bool same_chain_prev = determine_same_chain_prev( ii, separate_chemical_entity );
		bool same_chain_next = determine_same_chain_next( ii, separate_chemical_entity );
		bool const check_Ntermini_for_this_chain = determine_check_Ntermini_for_this_chain( chainID );
		bool const check_Ctermini_for_this_chain = determine_check_Ctermini_for_this_chain( chainID );

		bool is_lower_terminus( ( ii == 1 || rinfos_.empty() || ! same_chain_prev  )
			&& check_Ntermini_for_this_chain && ! lower_terminus_is_occupied_according_to_link_map( resid ) );
		bool is_upper_terminus( ( ii == nres_pdb || ! same_chain_next ) && check_Ctermini_for_this_chain  && ! upper_terminus_is_occupied_according_to_link_map( resid ) );

		// Determine if this residue is a D-AA residue, an L-AA residue, or neither.
		StructFileRep::ResidueCoords const & xyz = rinfo.xyz();
		bool is_d_aa = NomenclatureManager::get_instance()->decide_is_d_aa( name3 );
		bool is_l_aa = NomenclatureManager::get_instance()->decide_is_l_aa( name3 );
		bool is_achiral = NomenclatureManager::get_instance()->decide_is_known_achiral( name3 );
		bool is_chemical_component_ligand = false;

		// Get a list of ResidueTypes that could apply for this particular 3-letter PDB residue name.
		if ( ! is_residue_type_recognized( ii, name3, is_chemical_component_ligand ) ) {
			residue_was_recognized_[ ii ] = false;
			continue;
		}

		// Don't try ld chirality detection on known CCD ligands -- that way leads madness
		// (too many potential atom names, currently special cased, that might overlap)
		if (
				!is_chemical_component_ligand &&
				( ( !(is_d_aa || is_l_aa || is_achiral ) ) || ( d_l_threeletter_codes_are_same_for_aa(name3, false) /*Returns false for non-aa.*/ ) )
				) {
			chemical::detect_ld_chirality_from_polymer_residue( xyz, name3, is_d_aa, is_l_aa );
		}

		if ( is_chemical_component_ligand ) {
			// We cannot handle non-polymeric residues that are chemical component derived
			// that nonetheless have LINKs. These will come up in a variety of circumstances.
			// (The reason is because it confuses terminus assignment, sort of.)
			is_lower_terminus = true;
			is_upper_terminus = true;
			//known_connect_atoms_on_this_residue.clear();
		}

		TR.Trace << "Residue " << ii << "(PDB file numbering: " << resid << " )" << std::endl;
		TR.Trace << "...same_chain_prev: " << same_chain_prev << std::endl;
		TR.Trace << "...same_chain_next: " << same_chain_next << std::endl;
		TR.Trace << "...is_lower_terminus: " << is_lower_terminus << std::endl;
		TR.Trace << "...check_Ntermini_for_this_chain: "<< check_Ntermini_for_this_chain << std::endl;
		TR.Trace << "...is_upper_terminus: " << is_upper_terminus << std::endl;
		TR.Trace << "...check_Ctermini_for_this_chain: "<< check_Ctermini_for_this_chain << std::endl;
		TR.Trace << "...last_residue_was_recognized: " << last_residue_was_recognized( ii ) << std::endl;
		TR.Trace << "...known connects this residue: " << known_connect_atoms_on_this_residue << std::endl;
		TR.Trace << "...is_d_aa: " << is_d_aa << std::endl;
		TR.Trace << "...is_l_aa: " << is_l_aa << std::endl;

		// Try to get RT using all connect atoms.
		ResidueTypeCOP rsd_type_cop = get_rsd_type( name3, ii, known_connect_atoms_on_this_residue,
			resid, is_lower_terminus, is_upper_terminus, is_d_aa, is_l_aa );

		int kk = known_connect_atoms_on_this_residue.size() - 1; // for easier wraparound logic
		while ( !rsd_type_cop && kk >= 0 ) {
			// Include just kk connect atoms and see if we work
			// kk-membered combinations of connect atoms, vec<vec<string>>
			auto power_sets = nmers_of( known_connect_atoms_on_this_residue, kk );
			//TR << "power_sets: " << power_sets << std::endl;
			for ( auto const & power_set : power_sets ) {
				//rsd_type_cop = get_rsd_type( name3, ii, all_but( known_connect_atoms_on_this_residue, utility::vector1< Size >( jj ) ),
				//TR << "power_set " << power_set << std::endl;
				rsd_type_cop = get_rsd_type( name3, ii, power_set,
					resid, is_lower_terminus, is_upper_terminus, is_d_aa, is_l_aa );
				if ( rsd_type_cop ) break;
			}

			--kk;
		}

		if ( rsd_type_cop == nullptr ) {
			std::string variant;
			if ( is_lower_terminus ) {
				variant += " lower-terminal";
			}
			if ( is_upper_terminus ) {
				variant += " upper-terminal";
			}
			utility_exit_with_message( "No match found for unrecognized residue at position " +
				std::to_string(ii) + "( PDB ID: " + resid + " )" +
				"\nLooking for" + variant + " residue with 3-letter code: " + name3 +
				( known_connect_atoms_on_this_residue.empty() ? "" : "\nWith inter-residue connections to atom(s) " + utility::join(known_connect_atoms_on_this_residue,",") ) +
				"\nThis can be caused by wrong residue naming. E.g. a BMA (beta-mannose) is named MAN (alpha-mannose)");
		}

		TR.Debug << "Initial type of " << ii << " is " << rsd_type_cop->name() << std::endl;
		TR.Debug << " terminal? " << is_lower_terminus << " same chain prev? " << same_chain_prev << std::endl;
		residue_types_[ ii ] = rsd_type_cop;
		is_lower_terminus_[ ii ] = is_lower_terminus;
		same_chain_prev_[ ii ] = same_chain_prev;
		fill_name_map( ii );
	}
}


/// Utility function to help pass_2_quick_and_dirty_resolve_residue_types()
core::chemical::ResidueTypeCOP
find_restype(
	std::string const & name,
	std::string const & name3,
	core::chemical::ResidueTypeSet const & restypeset,
	std::set< std::string > & warned_on_ccd_confusion
) {
	core::chemical::ResidueTypeCOP rsd_type_cop = restypeset.name_mapOP( name );

	if ( rsd_type_cop != nullptr && restypeset.has_name( "pdb_" + name3 ) && warned_on_ccd_confusion.count( name3 ) == 0 ) { // Better check (e.g. without generating the PDB type.)
		TR.Warning << "The Rosetta ResidueType for code " << name3 << " is a different chemical entity than the PDB's version." << std::endl;
		TR.Warning << "    You may need to use `-fast_restyping false` to get the desired results." << std::endl;
		warned_on_ccd_confusion.insert( name3 );
	}

	return rsd_type_cop;
}

/// @details Like the pass_2_resolve_residue_types() function above.
/// However, we assume that if the three letter code corresponds to a (full!) name of the ResidueType,
/// unless we have the full name of the ResidueType specified in the HETNAM records
/// (this should be written by the PDB writer).
/// Failing that, we assume that the residue comes from the CCD.
void
PoseFromSFRBuilder::pass_2_quick_and_dirty_resolve_residue_types()
{
	using namespace core::io::pdb;
	using namespace core::chemical;

	pre_process_residue_type_information();

	Size const nres_pdb( rinfos_.size() );

	std::set< std::string > warned_on_ccd_confusion;

	for ( Size ii = 1; ii <= rinfos_.size(); ++ii ) {

		core::io::ResidueInformation const & rinfo = rinfos_[ ii ];
		std::string const & name3 = rinfo.rosetta_resName();
		std::string const & resid = rinfo.resid();

		/////////////////////////
		// Handle residue to skip
		if ( name3 == "UNL" ) {
			residue_was_recognized_[ ii ] = false;
			continue;
		}

		if ( name3 == "HOH" && options_.ignore_waters() ) {
			output_ignore_water_warning_once();

			if ( options_.remember_unrecognized_res() && options_.remember_unrecognized_water() ) {
				remember_unrecognized_res( ii );
			}

			residue_was_recognized_[ ii ] = false;
			continue;
		}

		/////////////////////////////////
		// Find the residue type needed.

		ResidueTypeCOP rsd_type_cop = nullptr;

		std::string residue_base_name = "";
		if ( sfr_.residue_type_base_names().count( resid ) ) {
			// We have annotated full name information
			residue_base_name = sfr_.residue_type_base_names()[ resid ].second;

			rsd_type_cop = residue_type_set_->name_mapOP( residue_base_name );
		}

		if ( rsd_type_cop == nullptr && residue_type_set_->has_name( name3 ) ) { // Full name lookup intentional!
			// We assume the three letter code corresponds to a Rosetta ResidueType full name
			rsd_type_cop = find_restype( name3, name3, *residue_type_set_, warned_on_ccd_confusion );
		}

		// The three letter code may have spaces, in which case we may need to strip it first.
		if ( rsd_type_cop == nullptr && (name3.front() == ' ' || name3.back() == ' ') ) {
			std::string stripped = utility::strip( name3 );
			// Special-case RNA residues
			if ( stripped == "A" || stripped == "rA" ) {
				stripped = "RAD";
			} else if ( stripped == "U" || stripped == "rU" ) {
				stripped = "URA";
			} else if ( stripped == "G" || stripped == "rG" ) {
				stripped = "RGU";
			} else if ( stripped == "C" || stripped == "rC" ) {
				stripped = "RCY";
			}
			rsd_type_cop = find_restype( stripped, name3, *residue_type_set_, warned_on_ccd_confusion );
		}

		// We attempt to special-case codes which have special AA designations (e.g. D-AA)
		if ( rsd_type_cop == nullptr && ! core::chemical::is_aa_name_unknown( name3 ) ) {
			core::chemical::AA aa = core::chemical::aa_from_name( name3 );
			std::string full_name = core::chemical::full_name_from_aa( aa );

			rsd_type_cop = find_restype( full_name, name3, *residue_type_set_, warned_on_ccd_confusion );
		}

		// Some additional Rosetta-specific Residues, which may be needed to round-trip
		if ( rsd_type_cop == nullptr && (name3 == "XXX" || name3 == "YYY" ) ) {
			if ( rinfo.xyz().count("ORIG") ) {
				if ( name3 == "XXX" ) {
					if ( residue_type_set_->mode() == FULL_ATOM_t ) {
						rsd_type_cop = residue_type_set_->name_mapOP( "GB_AA_PLACEHOLDER:NtermProteinFull" ); // Yup, probably should be VRT, but that's what we're actually doing with the normal way
					} else {
						rsd_type_cop = residue_type_set_->name_mapOP( "VRT" );
					}
				}
				if ( name3 == "YYY" ) {
					rsd_type_cop = residue_type_set_->name_mapOP( "INV_VRT" );
				}
			}
		}

		// At this point we can go directly to the PDB components
		if ( rsd_type_cop == nullptr && residue_type_set_->has_name( "pdb_" + name3 ) ) {
			rsd_type_cop = residue_type_set_->name_mapOP( "pdb_" + name3 );
		}

		//////////////////////////////////////
		// Handle residues we can't understand

		if ( rsd_type_cop == nullptr ) {
			if ( !options_.ignore_unrecognized_res() ) {
				std::string message = "No match found for unrecognized residue at position " +
					std::to_string(ii) + "( PDB ID: " + resid + " ) 3-letter code " + name3;
				if ( ! residue_base_name.empty() ) {
					message += " with base name " + residue_base_name;
				}
				if ( chemical::carbohydrates::CarbohydrateInfoManager::is_valid_sugar_code( name3 ) ) {
					message += "; did you intend to use the -include_sugars flag?";
				}
				message += "\nYou may have better success with the `-fast_restyping false` option.";
				utility_exit_with_message( message );
			}
			if ( options_.remember_unrecognized_res() ) {
				remember_unrecognized_res( ii );
			}
			residue_was_recognized_[ ii ] = false;
			continue;
		}

		///////////////////////////////////////
		// Some ugly special-casing
		// This isn't a license to expand this -- the HETNAM method should be the prefered approach.

		if ( rsd_type_cop->aa() == core::chemical::aa_his ) {
			// If we have the one tautomer proton (but not both!) we type as HIS_D.
			if ( rinfo.xyz().count(" HD1") && ! rinfo.xyz().count(" HE2") ) {
				rsd_type_cop = residue_type_set_->name_mapOP( "HIS_D" );
			}
		}

		/////////////////////////////////////
		// Deal with polymers and connections

		char chainID = rinfo.chainID();

		bool const separate_chemical_entity = determine_separate_chemical_entity( chainID );
		bool same_chain_prev = determine_same_chain_prev( ii, separate_chemical_entity );
		bool same_chain_next = determine_same_chain_next( ii, separate_chemical_entity );
		bool const check_Ntermini_for_this_chain = determine_check_Ntermini_for_this_chain( chainID );
		bool const check_Ctermini_for_this_chain = determine_check_Ctermini_for_this_chain( chainID );

		bool is_lower_terminus( ( ii == 1 || rinfos_.empty() || ! same_chain_prev  )
			&& check_Ntermini_for_this_chain && ! lower_terminus_is_occupied_according_to_link_map( resid ) );
		bool is_upper_terminus( ( ii == nres_pdb || ! same_chain_next ) && check_Ctermini_for_this_chain  && ! upper_terminus_is_occupied_according_to_link_map( resid ) );

		bool update_patches = false;

		// This is a quick & dirty typing method.
		// Don't add any checks which greatly increase computation.
		// (Prefer the HETNAM method for full restype specification)
		if ( is_lower_terminus && rsd_type_cop->lower_connect_id() != 0 ) {
			update_patches = true;
		}
		if ( is_upper_terminus && rsd_type_cop->upper_connect_id() != 0 ) {
			update_patches = true;
		}

		//----- the hairy stuff happens here -----------------------
		utility::vector1< std::string > known_connect_atoms_on_this_residue;
		determine_residue_branching_info( ii, known_connect_atoms_on_this_residue, known_links_ );
		//----------------------------------------------------------
		for ( std::string const & atm: known_connect_atoms_on_this_residue ) {
			if ( !rsd_type_cop->has(atm) || rsd_type_cop->n_residue_connections_for_atom( rsd_type_cop->atom_index(atm) ) < 1 ) {
				update_patches = true;
			}
		}

		TR.Debug << "Residue " << ii << "(PDB file numbering: " << resid << " ) name3 " << name3 << std::endl;
		//TR.Debug << "...same_chain_prev: " << same_chain_prev << std::endl;
		//TR.Debug << "...same_chain_next: " << same_chain_next << std::endl;
		TR.Debug << "...is_lower_terminus: " << is_lower_terminus << std::endl;
		//TR.Debug << "...check_Ntermini_for_this_chain: "<< check_Ntermini_for_this_chain << std::endl;
		TR.Debug << "...is_upper_terminus: " << is_upper_terminus << std::endl;
		//TR.Debug << "...check_Ctermini_for_this_chain: "<< check_Ctermini_for_this_chain << std::endl;
		//TR.Debug << "...last_residue_was_recognized: " << last_residue_was_recognized( ii ) << std::endl;
		TR.Debug << "...known connects this residue: " << known_connect_atoms_on_this_residue << std::endl;
		//TR.Debug << "...is_d_aa: " << is_d_aa << std::endl;
		//TR.Debug << "...is_l_aa: " << is_l_aa << std::endl;

		if ( update_patches ) {
			TR.Debug << "Need to update patches for " << rsd_type_cop->name() << " -- " << is_lower_terminus << " " << is_upper_terminus << " ";
			for ( std::string const & atm: known_connect_atoms_on_this_residue ) {
				TR.Debug << atm << " ";
			}
			TR.Debug << std::endl;

			// We fall back to the ResidueTypeFinder, but hopefully with a method which is quick.
			utility::vector1< ResidueProperty > properties;
			utility::vector1< VariantType > variants;

			if ( is_lower_terminus ) {
				if ( known_connect_atoms_on_this_residue.contains( "P" ) || known_connect_atoms_on_this_residue.contains( "N" ) ) {
					variants.push_back( CUTPOINT_UPPER );
				} else {
					properties.push_back( LOWER_TERMINUS );
				}
				// equivalent of check_nucleic_acid_virtual_phosphates:
				if ( rsd_type_cop->is_DNA() ) {
					variants.push_back( VIRTUAL_DNA_PHOSPHATE );
				}
				if ( rsd_type_cop->is_RNA() ) {
					variants.push_back( VIRTUAL_PHOSPHATE );
				}
			}
			if ( is_upper_terminus ) {
				if ( known_connect_atoms_on_this_residue.contains( "O3'" ) || known_connect_atoms_on_this_residue.contains( "C" ) ) {
					variants.push_back( CUTPOINT_LOWER );
				} else {
					properties.push_back( UPPER_TERMINUS );
				}
			}

			// TODO: Find a quicker way of patching for variant and property
			ResidueTypeCOP patched_rsd_type = ResidueTypeFinder( *residue_type_set_ )
				.name3( name3 )
				.residue_base_name( residue_base_name )
				.variants( variants )
				.properties( properties )
				.connect_atoms( known_connect_atoms_on_this_residue )
				.get_representative_type( false ); // Deliberately take the first and avoid metapatches for speed.

			if ( patched_rsd_type == nullptr ) {
				TR.Warning << "Attempted to find a ResidueType with three letter code " << name3;
				if ( !residue_base_name.empty() ) {
					TR.Warning << " (" << residue_base_name << ")";
				}
				TR.Warning << ( is_lower_terminus ? " lower-terminal" :"" );
				TR.Warning << ( is_upper_terminus ? " upper-terminal" :"" );
				if ( !known_connect_atoms_on_this_residue.empty() ) {
					TR.Warning << " with connections at";
					for ( std::string const & atm: known_connect_atoms_on_this_residue ) {
						TR.Warning << " " << atm << ",";
					}
				}
				TR.Warning << " -- Attempt was unsuccessful!" << std::endl;
			} else {
				TR.Debug << "Updated " << rsd_type_cop->name() << " to " << patched_rsd_type->name() << std::endl;
				rsd_type_cop = patched_rsd_type; // Successfull -- keep the patched varaint
			}
		}

		is_lower_terminus = is_lower_terminus || rsd_type_cop->lower_connect_id() == 0;

		TR.Debug << "Initial type of " << ii << " is " << rsd_type_cop->name() << std::endl;
		TR.Debug << " terminal? " << is_lower_terminus << " same chain prev? " << same_chain_prev << std::endl;
		residue_types_[ ii ] = rsd_type_cop;
		is_lower_terminus_[ ii ] = is_lower_terminus;
		same_chain_prev_[ ii ] = same_chain_prev;
		fill_name_map( ii );
	}
}

/// @details Common code from the pass_2 functions which do some (not computationally intensive) ResidueType annotation.
void PoseFromSFRBuilder::pre_process_residue_type_information() {

	// Convert PDB 3-letter code to Rosetta 3-letter code, if a list of alternative codes has been provided.
	for ( Size ii = 1; ii <= rinfos_.size(); ++ii ) {
		core::io::ResidueInformation const & rinfo = rinfos_[ ii ];

		std::pair< std::string, std::string > const & rosetta_names(
			NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( rinfo.resName() ) );
		std::string const & name3( rosetta_names.first );
		std::string const & base_name( rosetta_names.second );

		if ( base_name != "" ) {
			sfr_.residue_type_base_names()[ rinfo.resid() ] = std::make_pair( name3, base_name );
			sfr_.default_mainchain_connectivity()[ rinfo.resid() ] =
				NomenclatureManager::get_instance()->default_mainchain_connectivity_from_pdb_code( rinfo.resName() );
		}
		if ( core::chemical::carbohydrates::CarbohydrateInfoManager::is_valid_sugar_code( name3 ) ) {
			TR.Trace << "Identified glycan at position " << ii << std::endl;
			glycan_positions_.push_back( ii );
		}
	}

	known_links_ = core::io::explicit_links_from_sfr_linkage( sfr_.link_map(), rinfos_ );
	if ( options_.auto_detect_glycan_connections() ) {
		TR.Debug << "Auto-detecting glycan connections" << std::endl;
		//This clears the links before fix_glycan_order call which makes the clear in add_glycan_links_to_map redundant but I'm leaving it for now.
		if ( !options_.maintain_links() ) {
			TR.Debug << "Clearing linkages read from provided LINK records" << std::endl;
			known_links_.clear();
		} else {
			TR.Debug << "Keeping linkages read from provided LINK records" << std::endl;
		}
		utility::vector1< core::Size > chain_ends = core::io::fix_glycan_order( rinfos_, glycan_positions_, options_, known_links_ );
		// Reset correspondences for reordered glycans
		for ( core::Size pos : glycan_positions_ ) {
			resid_to_index_[ rinfos_[ pos ].resid() ] = pos;
		}
		for ( core::Size end : chain_ends ) {
			if ( end+1 <= same_chain_prev_.size() ) {
				TR << "Setting chain termination for " << end << std::endl;
				same_chain_prev_[ end + 1 ] = false;
			}
		}
		core::io::add_glycan_links_to_map( known_links_, core::io::determine_glycan_links( rinfos_, options_ ), rinfos_ );
		TR.Debug << "Finished auto-detecting glycan connections" << std::endl;
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
				core::uint const shadowed_atom_index( iirestype.atom_being_shadowed( iimainchain[ jj ] ) );  // required for rare short-mainchain cyclic residues
				std::string const & name1( (shadowed_atom_index) ? iirestype.atom_name( shadowed_atom_index ) : iirestype.atom_name( iimainchain[ jj ] ) );
				std::string const & name2( iirestype.atom_name( iimainchain[ jj + 1 ] ) );
				std::string const & name3( iirestype.atom_name( iimainchain[ jj + 2 ] ) );

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
				TR.Warning << "skipping pdb residue b/c it's missing too many mainchain atoms: " <<
					rinfos_[ ii ].resid() << ' ' << rinfos_[ii].rosetta_resName() << ' ' << iirestype.name() << std::endl;
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
				residue_types_[ ii ] = nullptr;
				residue_was_recognized_[ ii ] = false;
			}
		}
	}
}

bool PoseFromSFRBuilder::lower_terminus_is_occupied_according_to_link_map( std::string const & resid ) {

	if ( sfr_.link_map().count( resid ) ) {  // if found in the linkage map
		// The link map is keyed by resID of each branch point.
		Size const n_branches( sfr_.link_map()[ resid ].size() );
		for ( uint branch( 1 ); branch <= n_branches; ++branch ) {
			LinkInformation const & link_info( sfr_.link_map()[ resid ][ branch ] );

			// AMW Might want to filter for protein residue types here...
			if ( link_info.name1 == " N  " || link_info.name1 == " P  " ) {
				return true;
			}
		}
	}

	for ( auto const & elem : sfr_.link_map() ) {
		for ( LinkInformation const & link_info : elem.second ) {
			if ( resid == link_info.resID2 && ( link_info.name2 == " N  " || link_info.name2 == " P  " ) ) {
				return true;
			}
		}
	}

	return false;
}

bool PoseFromSFRBuilder::upper_terminus_is_occupied_according_to_link_map( std::string const & resid ) {
	if ( sfr_.link_map().count( resid ) ) {  // if found in the linkage map
		// The link map is keyed by resID of each branch point.
		Size const n_branches( sfr_.link_map()[ resid ].size() );
		for ( uint branch( 1 ); branch <= n_branches; ++branch ) {
			LinkInformation const & link_info( sfr_.link_map()[ resid ][ branch ] );

			// AMW Might want to filter for protein residue types here...
			if ( link_info.name1 == " C  " || link_info.name1 == " O3'" ) {
				return true;
			}
		}
	}

	for ( auto const & elem : sfr_.link_map() ) {
		for ( LinkInformation const & link_info : elem.second ) {
			if ( resid == link_info.resID2 && ( link_info.name2 == " C  " || link_info.name2 == " O3'" ) ) {
				return true;
			}
		}
	}

	return false;
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
		std::string const & resid = rinfo.resid();
		char chainID = rinfo.chainID();

		if ( ! determine_check_Ntermini_for_this_chain( chainID ) ) { continue; }
		if ( ! determine_check_Ctermini_for_this_chain( chainID ) ) { continue; }

		Size const ii_prev = prev_residue_skipping_null_residue_types( ii );
		Size const ii_next = next_residue_skipping_null_residue_types( ii );

		bool type_changed( false );
		if ( ! residue_types_[ ii ]->is_polymer() ) { continue; }
		if ( ! residue_types_[ ii ]->is_lower_terminus() && residue_types_[ ii ]->lower_connect_id() != 0 &&
				( ii == 1 || ii == ii_prev || ( residue_types_[ ii_prev ] && (
				( !residue_types_[ ii_prev ]->is_polymer() && !residue_types_[ ii]->residue_connection_is_polymeric( residue_types_[ ii]->lower_connect_id() ) ) ||
				residue_types_[ ii_prev ]->is_upper_terminus() ) ) ) ) {
			std::string lower_atom = residue_types_[ ii ]->atom_name( residue_types_[ ii ]->lower_connect_atom() );
			// Check to see if we have an explicit connection to the lower atom, and only adjust if we don't
			if ( ! ( known_links_.count( resid ) && known_links_[ resid ].count( lower_atom ) > 0 ) ) {
				// Do not adjust if it's an upper cutpoint. That variant type would have been assigned
				// for a good reason.
				if ( ! residue_types_[ ii ]->has_variant_type( chemical::CUTPOINT_UPPER ) ) {
					TR << "Adding undetected lower terminus type to residue " << ii << ", " << resid << std::endl;
					residue_types_[ ii ] = residue_type_set_->get_residue_type_with_variant_added(
						*residue_types_[ ii ], chemical::LOWER_TERMINUS_VARIANT ).get_self_ptr();
					type_changed = true;
				} else {
					TR.Debug << "Would have added an undetected lower terminus variant but for the " << std::endl;
					TR.Debug << "already-detected cutpoint-upper type!" << std::endl;
				}
			}
		}
		// AMW: new--don't add  an upper terminus type to protein residues following RNA
		// if that RNA is itself not upper-terminal
		if ( !residue_types_[ ii ]->is_upper_terminus() && residue_types_[ ii ]->upper_connect_id() != 0 &&
				// !( residue_types_[ ii ]->is_polymer() && residue_types_[ ii_prev ]->is_RNA()
				// && !residue_types_[ ii_prev ]->is_upper_terminus() ) &&
				// !( residue_types_[ ii ]->is_RNA() && residue_types_[ ii_next ]->is_polymer()
				// && !residue_types_[ ii_next ]->is_upper_terminus() ) &&
				( ii == ii_next ||
				( residue_types_[ ii_next ] && (
				!residue_types_[ ii_next ]->is_polymer() ||
				residue_types_[  ii_next ]->is_lower_terminus() ||
				chainID != rinfos_[ ii_next ].chainID() /* ||
				residue_types_[  ii_next ]->has_variant_type(BRANCH_LOWER_TERMINUS_VARIANT)*/ ) ) ) ) {

			// AMW TODO: print the residue type pdb_SEE and see what's wrong with it:
			// Why does it think that it can't make an upper terminus version of itself?
			// It wants it but can't have it -- one of those is Bad.
			std::string upper_atom = residue_types_[ ii ]->atom_name( residue_types_[ ii ]->upper_connect_atom() );
			if ( ! ( known_links_.count( resid ) && known_links_[ resid ].count( upper_atom ) > 0 ) ) {
				TR << "Adding undetected upper terminus type to residue " << ii << ", " << resid << std::endl;
				residue_types_[ ii ] = residue_type_set_->get_residue_type_with_variant_added(
					*residue_types_[ ii ], chemical::UPPER_TERMINUS_VARIANT ).get_self_ptr();
				type_changed = true;
			}
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
			TR.Warning << "discarding " << discarded_atoms
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


	// Setup the set of residues to use.
	// Note that we're deferring adding them to the pose, so we can't use the "current" state of the pose for any logic here

	utility::vector1< conformation::ResidueOP > residues;
	utility::vector1< core::Size > jump_connections; // Which residue to jump from, 0 for polymeric connection to previous

	utility::vector1< core::Size > chain_endings; // Which residues are the end of their respective chains?

	for ( Size ii = 1; ii <= residue_types_.size(); ++ii ) {
		if ( ! residue_types_[ ii ] ) { continue; }

		ResidueType const & ii_rsd_type( *residue_types_[ ii ] );
		TR.Trace << "ResidueType " << ii << ": " << ii_rsd_type.name() << std::endl;
		//TR.Trace << "went by, in the file context, " << rinfos_[ ii ].resName()
		// << " at " << rinfos_[ ii ].chainID() << rinfos_[ ii ].resSeq() << rinfos_[ ii ].iCode() << " " << rinfos_[ ii ].segmentID() << std::endl;

		ResidueOP ii_rsd( ResidueFactory::create_residue( ii_rsd_type ) );
		for ( auto iter = rinfos_[ ii ].xyz().begin(), iter_end = rinfos_[ ii ].xyz().end();
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
				// +1 here as we haven't added the residue yet.
				id::NamedAtomID atom_id( pose_name, residues.size()+1 );
				coordinates_assigned_.set( atom_id, true);
			}
			//else runtime_assert( iter->first == " H  " && ii_rsd_type.is_terminus() ); // special casee
		}

		// store residue name as in PDB file for interpretable error messages
		std::string const & PDB_resid = rinfos_[ ii ].resid();

		// These sister atoms are THUS FAR only found in RNA (it's to perfectly copy
		// native structure coordinates). So, at the moment, we only need to call
		// this function for RNA. Eventually, we may implement more proteinaceous use
		// cases, at which point you may want to move this condition into the
		// check_and_correct function.
		// AMW: This restriction means some odd behavior for automatically set up DNA residues
		// with "wrong" chirality phosphates. But changing this means that we make a lot
		// of swaps for DNA residues. We will have to see what's best, but not for now.
		if ( ii_rsd->is_RNA() ) {
			check_and_correct_sister_atoms( ii_rsd );
		}

		Size const old_nres( residues.size() );

		TR.Trace << "...new residue is a polymer: " << ii_rsd->type().is_polymer() << std::endl;
		if ( old_nres >= 1 ) {
			TR.Trace << "...old residue is a polymer: " << residues[ old_nres ]->is_polymer() << std::endl;
		}

		// Add the first new residue
		if ( !old_nres ) {
			TR.Trace << ii_rsd_type.name() << " " << ii << " is the start of a new pose." << std::endl;

			residues.push_back( ii_rsd );
			jump_connections.push_back( 0 );

			// If this is a lower terminus, AND it's not about to get UPPER-UPPER bonded from prev
			// Actually, building initial PDBs like this is all right: we will be
			// doing more complicated checks later.
		} else if (
				( is_lower_terminus_[ ii ] && determine_check_Ntermini_for_this_chain( rinfos_[ ii ].chainID() ) ) ||
				( !residues[ old_nres ]->type().is_upper_terminus() && residues[ old_nres ]->type().is_RNA() && !ii_rsd_type.is_upper_terminus() && ii_rsd_type.is_lower_terminus() ) ||
				! same_chain_prev_[ ii ] || // explicit new-chain signal
				residues[ old_nres ]->type().has_variant_type( "C_METHYLAMIDATION" ) || // C_METHYLAMIDATED can't bond, so jump AMW TODO VARIANT
				! ii_rsd->is_polymer() || // new residue isn't polymer, so jump
				// Don't connect to previous with an UPPER/LOWER polymeric bond if previous didn't have an upper.
				! residues[ old_nres ]->type().is_polymer() || // Nonpolymers can't bond-to-next automatically, so jump
				residues[ old_nres ]->type().is_upper_terminus() || // generalized from C_METHYLAMIDATION
				! last_residue_was_recognized( ii ) ) {
			// A new chain because this is a lower terminus (see logic above for designation)
			// and if we're not checking it then it's a different chain from the previous
			core::Size root_index = 1;
			// connect metal ions by a jump to the closest metal-binding residue that is lower in sequence.
			if ( ii_rsd->is_metal() && basic::options::option[basic::options::OptionKeys::in::auto_setup_metals] ) {
				root_index = find_atom_tree_root_for_metal_ion( residues, ii_rsd );
			}
			if ( root_index>1 ) {
				TR << ii_rsd_type.name() << " " << ii;
				TR << " was added by a jump, with base residue " << root_index << std::endl;
			}

			residues.push_back( ii_rsd );
			jump_connections.push_back( root_index );
			if ( residues.size() > 1 ) {
				chain_endings.push_back( residues.size() - 1 );
			}

		} else { // Append residue to current chain dependent on bond length.
			if ( ! options_.missing_dens_as_jump() ) {
				TR.Trace << ii_rsd_type.name() << " " << ii;
				TR.Trace << " (PDB residue: " << PDB_resid << ")";
				TR.Trace << " is appended to chain " << rinfos_[ ii ].chainID() << std::endl;

				residues.push_back( ii_rsd );
				jump_connections.push_back( 0 );

			} else {
				//fpd look for missing density in the input PDB
				//fpd if there is a bondlength > 3A
				//fpd we will consider this missing density
				ResidueOP last_rsd = residues[old_nres];
				core::Real bondlength = ( last_rsd->atom( last_rsd->upper_connect_atom() ).xyz() -
					ii_rsd->atom( ii_rsd->lower_connect_atom() ).xyz() ).length();
				if ( bondlength > 3.0 ) {
					TR.Warning << "missing density found at residue (rosetta number) " << old_nres << std::endl;

					// The pose here is only used for getting the appropriate ResidueTypeSet
					// (Missing atoms are filled in later.)

					if ( last_rsd->is_protein() ) {
						if ( !last_rsd->has_variant_type( UPPER_TERMINUS_VARIANT ) &&
								!last_rsd->has_variant_type( UPPERTERM_TRUNC_VARIANT ) ) {
							residues[ old_nres ] = this->add_variant_type_to_residue( last_rsd, UPPERTERM_TRUNC_VARIANT, pose );
						}
					} else {
						if ( !last_rsd->has_variant_type( UPPER_TERMINUS_VARIANT ) ) {
							residues[ old_nres ] = this->add_variant_type_to_residue( last_rsd, UPPER_TERMINUS_VARIANT, pose );
						}
						if ( !ii_rsd->has_variant_type( LOWER_TERMINUS_VARIANT) ) {
							TR << " add lower variant " << old_nres+1 << std::endl;
							ii_rsd = this->add_variant_type_to_residue( ii_rsd, LOWER_TERMINUS_VARIANT, pose );
						}
					}

					if ( ii_rsd->is_protein() ) {
						// If it doesn't have a lower term variant for some other reason
						if ( !ii_rsd->has_variant_type( LOWER_TERMINUS_VARIANT ) ) {
							ii_rsd = this->add_variant_type_to_residue( ii_rsd, LOWERTERM_TRUNC_VARIANT, pose );
						}
					} else if ( !ii_rsd->is_carbohydrate() ) {
						ii_rsd = this->add_variant_type_to_residue( ii_rsd, LOWER_TERMINUS_VARIANT, pose );
					}

					residues.push_back( ii_rsd );
					jump_connections.push_back( old_nres );
					// Keep this in the same chain -- no chain ending.


				} else {
					TR.Trace << ii_rsd_type.name() << " " << ii << " is appended to chain" << rinfos_[ ii ].chainID() << std::endl;

					residues.push_back( ii_rsd );
					jump_connections.push_back( 0 );

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
				TR.Warning << "the coordinates for this input structure may have been poorly assigned." << std::endl;
				TR.Debug << "  Measured: " << ii_rsd->ring_conformer( i ).specific_name << "  Expected: ";
				TR.Debug << ii_rsd_type.ring_conformer_set( i )->get_lowest_energy_conformer().specific_name << std::endl;

			}
		}

		pose_to_rinfo_.push_back( ii );
		pose_temps_.push_back( rinfos_[ ii ].temps() );

	}

	// Now actually append all the residues
	try {
		pose.append_residues( residues, jump_connections );
	} catch (utility::excn::Exception & e) {
		utility_exit_with_message( e.msg() );
	}

	// Make sure chain endings are up-to-date.
	for ( core::Size end: chain_endings ) {
		pose.conformation().insert_chain_ending( end );
	}

	TR.Trace << "Initial, pre-refined Pose built successfully." << std::endl;
}

core::conformation::ResidueOP
PoseFromSFRBuilder::add_variant_type_to_residue( core::conformation::ResidueOP old_rsd, chemical::VariantType const variant_type, core::pose::Pose const & pose )
{
	if ( old_rsd->has_variant_type( variant_type ) ) { return old_rsd; }

	// the type of the desired variant residue
	chemical::ResidueTypeSetCOP rsd_set( pose.residue_type_set_for_pose( old_rsd->type().mode() ) );
	chemical::ResidueType const & new_rsd_type( rsd_set->get_residue_type_with_variant_added( old_rsd->type(), variant_type ) );

	core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( new_rsd_type ) );

	// Copy coordinates
	for ( Size ii=1; ii <= new_rsd->natoms(); ++ii ) {
		std::string const & atom_name( new_rsd->atom_name(ii) );
		if ( old_rsd->has( atom_name ) ) {
			new_rsd->atom( ii ).xyz( old_rsd->atom( atom_name ).xyz() );
		} // We ignore missing atoms, as we'll fill them in later in the build process.
	}

	return new_rsd;
}

bool
is_connected( core::conformation::Conformation const & conf, core::Size ii, std::string const & ii_atm, core::Size jj, std::string const & jj_atm ) {
	core::conformation::Residue const & ii_res( conf.residue( ii ) );
	core::conformation::Residue const & jj_res( conf.residue( jj ) );
	if ( ! ii_res.is_bonded( jj ) || ! jj_res.is_bonded( ii ) ) { return false; }
	Size ii_conn = ii_res.connect_atom( jj_res );
	if ( ii_conn != ii_res.atom_index( ii_atm ) ) { return false; }
	Size jj_conn = jj_res.connect_atom( ii_res );
	if ( jj_conn != jj_res.atom_index( jj_atm ) ) { return false; }
	return true;
}

void
show_residue_connections( core::conformation::Conformation const & conf, core::Size i ) {
	core::conformation::Residue const &res( conf.residue(i) );
	Size const nconn(res.n_possible_residue_connections());
	TR << "RESCON: " << i << ' ' << res.name() << " n-conn= " << nconn <<
		" n-poly= " << res.n_polymeric_residue_connections() <<
		" n-nonpoly= " << res.n_non_polymeric_residue_connections();
	for ( Size j = 1; j <= nconn; ++j ) {
		TR << " conn# " << j << ' ' << res.residue_connect_atom_index( j )
			<< ' ' << res.connect_map(j).resid() << ' ' <<
			res.connect_map(j).connid();
	}
	TR << std::endl;
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
	using namespace io;
	using namespace pdb;
	using namespace chemical;
	using namespace conformation;


	//typedef std::map< std::string, core::Real > ResidueTemps;

	// Note: _do not_ access Residue here. You do not need to. Doing so triggers
	// a refold and that is expensive. You probably can get away with
	// pose.residue_type( ii ). If you think you need Residue, ask someone else.

	// Step 2. Handle missing atoms.
	//Size num_heavy_missing = 0;
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
		ResidueType const & ii_rsd( pose.residue_type( ii ) ); // Type to avoid refolding at this moment.
		for ( Size jj = 1; jj <= ii_rsd.natoms(); ++jj ) {
			id::AtomID atom_id( jj, ii );
			id::NamedAtomID named_atom_id( ii_rsd.atom_name( jj ), ii );
			if ( ! coordinates_assigned_.has( named_atom_id ) || ! coordinates_assigned_[ named_atom_id ] ) {
				missing_[ atom_id ] = true;
				//if ( !ii_rsd.atom_is_hydrogen( jj ) ) num_heavy_missing++;
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
		ResidueType const & rsd_type = pose.residue_type( ii );
		if ( ! rsd_type.is_DNA() ) continue;

		for ( Size jj = 1, natoms = rsd_type.natoms(); jj <= natoms; ++jj ) {
			id::AtomID const id( jj, ii );
			if ( missing_[ id ] && rsd_type.atom_name( jj ) == " P  " ) {
				runtime_assert( rsd_type.has_variant_type( chemical::VIRTUAL_DNA_PHOSPHATE ) );
				break;
			}
		}
	}

	// Step 5. Reconcile connectivity data.

	//Store carbohydrate atoms with ambiguous lower connect so that we can regenerate coordinates after lower connects are resolved. This should probably be changed to apply to all non proteins.
	utility::vector1< core::id::AtomID > lower_connect_atoms;
	for ( Size ii=1; ii<=pose.size(); ++ii ) {
		ResidueType const & restype = pose.residue_type( ii );

		if ( !restype.is_carbohydrate() ) continue;
		for ( Size jj=1; jj<=restype.natoms(); ++jj ) {
			if ( ( restype.icoor(jj).depends_on_polymer_lower() && !restype.is_lower_terminus() ) ||
					( restype.icoor(jj).depends_on_polymer_upper() && !restype.is_upper_terminus() ) ) {
				core::id::AtomID const
					stub_atom1( pose.residue(ii).icoor( jj ).stub_atom1().atom_id( pose.residue(ii), pose.conformation() ) ),
					stub_atom2( pose.residue(ii).icoor( jj ).stub_atom2().atom_id( pose.residue(ii), pose.conformation() ) ),
					stub_atom3( pose.residue(ii).icoor( jj ).stub_atom3().atom_id( pose.residue(ii), pose.conformation() ) );
				if ( stub_atom1 == id::AtomID::BOGUS_ATOM_ID() || stub_atom2 == id::AtomID::BOGUS_ATOM_ID() || stub_atom3 == id::AtomID::BOGUS_ATOM_ID() ) {
					lower_connect_atoms.push_back( core::id::AtomID(jj,ii) );
				}
			}
		}
	}

	// AMW LOOK HERE: you may need this area for ligand LINKs.
	// Add any links from the known link data (e.g. the link map)
	// that isn't already in the pose.
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		core::Size rinfo_ii( pose_to_rinfo_[ ii ] );
		std::string const & resid( rinfos_[ rinfo_ii ].resid() );
		if ( known_links_.count( resid ) ) {
			for ( auto const & link_pair : known_links_[resid] ) {
				//std::string const & my_atom( link_pair.first );
				std::string const & partner_resid( link_pair.second.first );
				core::Size partner_rinfo_ii( resid_to_index_[ partner_resid ] );
				core::Size partner( pose_to_rinfo_.index( partner_rinfo_ii ) );
				//std::string const & partner_atom( link_pair.second.second );

				// Just strip my_atom and partner_atom (CIF defense)
				std::string my_atom = stripped_whitespace( link_pair.first );
				std::string partner_atom = stripped_whitespace( link_pair.second.second );

				if ( partner == 0 ) {
					TR.Warning << "Cannot find " << partner_resid << " in Pose -- skipping connection" << std::endl;
					continue;
				}
				TR.Debug << "Dealing with link ii:" << ii << " resid: " << resid << " (" << pose.residue_type( ii ).name() << ") atom: " << my_atom << " partner: " << partner << " (" << pose.residue_type( partner ).name() << ") partner_resid: " << partner_resid << " partner_atom: " << partner_atom << std::endl;
				if ( ! pose.residue_type(ii).has( my_atom ) || ! pose.residue_type( partner ).has( partner_atom ) ) {
					TR.Warning << "Cannot form link between " << resid << " " << my_atom
						<< " and " << partner_resid << " " << partner_atom << " as atom(s) don't exist." << std::endl;
					continue;
				}
				if ( ! is_connected( pose.conformation(), ii, my_atom, partner, partner_atom ) ) {
					if ( pose.residue( ii ).has_incomplete_connection( pose.residue( ii ).atom_index( my_atom ) ) &&
							pose.residue( partner ).has_incomplete_connection( pose.residue( partner ).atom_index( partner_atom ) ) ) {
						TR.Debug << "Making a chemical connection between residue " << ii << " " << my_atom << " and residue " << partner << " " << partner_atom << std::endl;
						pose.conformation().declare_chemical_bond( ii, my_atom, partner, partner_atom );
					} else {
						// Try to find connectable variants.

						// First: if either atom is a hydrogen, continue.
						if ( pose.residue_type( ii ).atom_is_hydrogen( pose.residue_type( ii ).atom_index( my_atom ) ) ) { continue; }
						if ( pose.residue_type( partner ).atom_is_hydrogen( pose.residue_type( partner ).atom_index( partner_atom ) ) ) { continue; }

						// Likewise, if residue is a carbohydrate, continue.
						if ( pose.residue_type( ii ).is_carbohydrate() ) {
							TR.Warning << "Carbohydrate residue " << ii << " not properly linked. ";
							TR.Warning << "A patch file for branching from the parent might be missing." << std::endl;
							continue;
						}

						// the type of the desired variant residue
						chemical::ResidueTypeSetCOP rsd_set( pose.residue_type_set_for_pose( pose.residue_type( ii ).mode() ) );
						//chemical::ResidueType const & new_rsd_type(  );

						core::pose::replace_pose_residue_copying_existing_coordinates( pose, ii,
							rsd_set->name_map( pose.residue_type( ii ).name() + core::chemical::PATCH_LINKER + "MP-" + my_atom + "-connect" ) );

						core::pose::replace_pose_residue_copying_existing_coordinates( pose, partner,
							rsd_set->name_map( pose.residue_type( partner ).name() + core::chemical::PATCH_LINKER + "MP-" + partner_atom + "-connect" ) );

						if ( pose.residue( ii ).has_incomplete_connection( pose.residue( ii ).atom_index( my_atom ) ) &&
								pose.residue( partner ).has_incomplete_connection( pose.residue( partner ).atom_index( partner_atom ) ) ) {
							TR.Debug << "Making a metapatched chemical connection between residue " << ii << " " << my_atom << " and residue " << partner << " " << partner_atom << std::endl;
							pose.conformation().declare_chemical_bond( ii, my_atom, partner, partner_atom );
						} else {

							// The ResidueType selection code should find a connectable type if one is present.
							TR.Warning << "Explicit link between " << rinfos_[ ii ].resid() << " " << my_atom
								<<" and " << rinfos_[ partner ].resid() << " " << partner_atom
								<< " requested, but no availible connections are present!" << std::endl;
							TR.Warning << "Types are " << pose.residue_type(ii).name() << " and " << pose.residue_type(partner).name() << std::endl;
							show_residue_connections( pose.conformation(), ii );
							show_residue_connections( pose.conformation(), partner );
						}
					}
				} // else we're already connected.
			}
		}
	}


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
	if ( !options_.no_detect_pseudobonds() && pose.size() > 1 ) {
		pose.conformation().detect_pseudobonds();
	}

	// this is where pdb info used to get stored.
	build_pdb_info_2_temps( pose );

	//Create the GlycanTreeSet if carbohydrates are present.
	// This requires up-to-data connectivity information and is required to be in the pose if glycan residues are present.
	if ( pose.conformation().contains_carbohydrate_residues() ) {
		pose.conformation().setup_glycan_trees();
	}

	//rebuild the atoms that were placed with an ambiguous lower connect
	for ( core::Size ii=1; ii<=lower_connect_atoms.size(); ii++ ) {
		core::Size resno = lower_connect_atoms[ii].rsd();
		core::Size atomno = lower_connect_atoms[ii].atomno();
		core::chemical::AtomICoor icoor = pose.residue(resno).type().icoor(atomno);
		if ( icoor.stub_atom1().atom_id( pose.residue(resno), pose.conformation()) == id::AtomID::BOGUS_ATOM_ID() || icoor.stub_atom2().atom_id(pose.residue(resno), pose.conformation()) == id::AtomID::BOGUS_ATOM_ID()
				|| icoor.stub_atom3().atom_id( pose.residue(resno), pose.conformation()) == id::AtomID::BOGUS_ATOM_ID() ) continue;
		numeric::xyzVector<core::Real> icoor_xyz = icoor.build(pose.residue(resno), pose.conformation());
		pose.conformation().set_xyz( lower_connect_atoms[ii], icoor_xyz);
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

	if ( sfr_.pdb_comments().size() ) {
		std::map< std::string, std::string > const & pdb_comments( sfr_.pdb_comments() );
		for ( auto const & pdb_comment : pdb_comments ) {
			core::pose::add_comment( pose, pdb_comment.first, pdb_comment.second );
		}
	}

	if ( sfr_.pose_cache_string_data().size() ) {
		std::map< std::string, std::string > const & string_data( sfr_.pose_cache_string_data());
		for ( auto const & data : string_data ) {
			core::pose::setPoseExtraScore(pose, data.first, data.second);
		}
	}

	if ( sfr_.pose_cache_real_data().size() ) {
		std::map< std::string, core::Real > const & real_data( sfr_.pose_cache_real_data());
		for ( auto const & data: real_data ) {
			core::pose::setPoseExtraScore(pose, data.first, data.second);
		}
	}

	if ( sfr_.simple_metric_data().has_data() ) {
		core::simple_metrics::MetricKey mk;
		core::simple_metrics::SimpleMetricDataOP sm_data = core::simple_metrics::get_sm_data(pose);
		sm_data->set_all_data(mk, sfr_.simple_metric_data());
	}
	TR.Trace << "Pose refined successfully:" << std::endl;
	TR.Trace << pose << std::endl;
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
		pdb_info->chain_sequences( sfr_.chain_sequences() );
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

	typedef std::map< std::string, core::Real > ResidueTemps;

	PDBInfoOP pdb_info = pose.pdb_info();

	// Add temperatures.
	for ( core::Size ii = 1; ii <= pose.size(); ii++ ) {
		ResidueTemps const & res_temps( rinfos_[ pose_to_rinfo_[ ii ] ].temps() );
		NameBimap const & namemap( remapped_atom_names_[ pose_to_rinfo_[ ii ] ] );
		for ( const auto & res_temp : res_temps ) {
			// namemap should only include atoms which have a presence in both rinfo and pose
			if ( namemap.left.count( res_temp.first ) ) {
				std::string const & pose_atom_name( namemap.left.find(res_temp.first)->second );
				if ( pose.residue_type( ii ).has( pose_atom_name ) ) { // There are issues with terminus patching which means atoms can sometimes disappear
					core::Size jj = pose.residue_type( ii ).atom_index( pose_atom_name );
					pdb_info->temperature( ii, jj, res_temp.second );
				}
			} else {
				if ( ( res_temp.first )[ 0 ] == 'H' || ( ( res_temp.first )[ 0 ] == ' ' && ( res_temp.first )[ 1 ] == 'H' ) ) {
					// Don't warn for missing temperatures on hydrogens.
				} else {
					TR.Warning << "can't find pose atom for file-residue " << ii << " atom " << res_temp.first << " (trying to store temperature in PDBInfo)" << std::endl;
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
/// @details  This function does three separate things that are related:
/// - It determines which atoms on this residue are listed as chemically
///   connected to other residues
/// - It assigns main-chain connectivity to carbohydrate ResidueType base
///   names.
/// - It turns off implicit connections to the next residue
///   for carbohydrates where we don't know linkage information.
void
PoseFromSFRBuilder::determine_residue_branching_info(
	Size const seqpos,
	utility::vector1< std::string > & known_connect_atoms_on_this_residue,
	std::map< std::string, std::map< std::string, std::pair< std::string, std::string > > > const & explicit_link_mapping )
{
	using namespace std;
	using namespace core::io::pdb;
	using namespace core::chemical;

	std::string const & resid = rinfos_[ seqpos ].resid();
	std::string const & name3 = rinfos_[ seqpos ].resName(); // for debugging

	// Carbohydrate base names will have "->?)-" as a prefix if their main-
	// chain connectivity requires LINK records to determine.  Fortuitously,
	// position 2 is also the index for the atom name (such as " O2 ") that
	// provides the missing information.
	int const CARB_MAINCHAIN_CONN_POS = 2;

	bool unknown_main_chain_connectivity( false );
	if ( sfr_.residue_type_base_names().count( resid ) ) {
		TR.Trace << "Current residue '" << name3 <<
			"' has had its base name extracted from the PDB file: ";
		TR.Trace << sfr_.residue_type_base_names()[ resid ].second << endl;
		unknown_main_chain_connectivity =
			sfr_.residue_type_base_names()[ resid ].second.size() > CARB_MAINCHAIN_CONN_POS &&
			( sfr_.residue_type_base_names()[ resid ].second[ CARB_MAINCHAIN_CONN_POS ] == '?' );
	}

	TR.Trace << "Checking if resid " << resid << "(" << seqpos << ")" << " is in the link map " << endl;
	if ( explicit_link_mapping.count( resid ) ) {  // if found in the linkage map
		TR.Trace << "Found resid " << resid << " in link map " << endl;
		// We want to sort the linkages by partner number - std::set will allow us to do this
		std::set< std::tuple<core::Size, std::string, std::string> > connections;
		for ( auto const & elm_pair : explicit_link_mapping.at(resid) ) {
			std::string const & link_atom = elm_pair.first;
			std::string const & partner_resid = elm_pair.second.first;
			//TR.Trace << "Link atom is " << link_atom << std::endl;
			debug_assert( resid_to_index_.count( partner_resid ) );
			core::Size const & partner = resid_to_index_[ partner_resid ];
			connections.insert( std::make_tuple( partner, partner_resid, link_atom) );
		}
		for ( auto const & elm_tuple : connections ) {
			core::Size const & partner = std::get<0>( elm_tuple );;
			//std::string const & partner_resid = std::get<1>( elm_tuple );
			std::string const & link_atom = std::get<2>( elm_tuple );
			//TR.Trace << "Link atom is " << link_atom << std::endl;

			if ( unknown_main_chain_connectivity ) {
				char const connectivity( link_atom[ CARB_MAINCHAIN_CONN_POS ] );
				std::string const & rosetta_name = sfr_.residue_type_base_names()[ resid ].first;
				debug_assert( chemical::carbohydrates::CarbohydrateInfoManager::is_valid_sugar_code( rosetta_name ) );
				char const anomeric = chemical::carbohydrates::CarbohydrateInfoManager::anomeric_position_from_code( rosetta_name );
				if ( connectivity != anomeric ) {
					TR.Trace << "Assigning main-chain connectivity to position " << connectivity;
					TR.Trace << " of this residue." << endl;
					sfr_.residue_type_base_names()[ resid ].second[ CARB_MAINCHAIN_CONN_POS ] = connectivity;
					unknown_main_chain_connectivity = false;
				}
			}

			//bi_map[ link_info.resID1 ][ link_info.name1 ] = make_pair( link_info.resID2, link_info.name2 );
			//bi_map[ link_info.resID2 ][ link_info.name2 ] = make_pair( link_info.resID1, link_info.name1 );

			if ( rinfos_[ seqpos ].chainID() == rinfos_[ partner ].chainID() && ( // same nominal chain
					( seqpos == partner-1 && rinfos_[ seqpos ].resSeq() == rinfos_[ partner ].resSeq()-1 ) // next residue (both in Pose & PDB numbering)
					|| ( seqpos == partner+1 && rinfos_[ seqpos ].resSeq() == rinfos_[ partner ].resSeq()+1 ) ) // previous residue (both in Pose & PDB numbering)
					&& // AND doesn't include a nonstandard atom (RNA)
					( link_atom != " O2'" && explicit_link_mapping.at( resid ).at( link_atom ).second != " O2'" ) ) {
				// If this occurs, the link is to the next residue on the same chain, so both residues are part of
				// the same main chain or branch, and this linkage information can be ignored, UNLESS this .pdb file
				// came from the PDB, in which case its 3-letter codes don't tell us the main chain, so we must get
				// the main chain from the LINK records.
				// Note that this all assumes insertion codes are not involved! It also assumes that the PDB file
				// makers did things reasonably.
				continue;
			}

			known_connect_atoms_on_this_residue.push_back( utility::strip( link_atom ) );
		}
	}

	// Fallback, if we're still not satisfied with mainchain connections
	if ( unknown_main_chain_connectivity ) {
		char default_mainchain_connectivity( sfr_.default_mainchain_connectivity()[ resid ] );
		if ( ! default_mainchain_connectivity ) {
			default_mainchain_connectivity = '3';
		}
		TR.Trace << "Assigning main-chain connectivity arbitrarily to position "<< default_mainchain_connectivity;
		TR.Trace << " of this terminal residue." << endl;

		sfr_.residue_type_base_names()[ resid ].second[ CARB_MAINCHAIN_CONN_POS ] = default_mainchain_connectivity;
		if ( next_residue_skipping_merges( seqpos ) != seqpos ) { // Don't update for last residue
			same_chain_prev_[ next_residue_skipping_merges( seqpos ) ] = false; // Don't connect this residue with the next
		}
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
		remember_unrecognized_res( pdb_residue_index );
	}

	if ( is_HOH_to_ignore ) output_ignore_water_warning_once();
	return false;
}

void
PoseFromSFRBuilder::remember_unrecognized_res( core::Size rinfo_index ) {
	core::io::ResidueInformation const & rinfo = rinfos_[ rinfo_index ];
	std::map< std::string, Vector > const & xyz( rinfo.xyz() );
	std::map< std::string, core::Real > const & rtemp( rinfo.temps() );

	for ( const auto & iter : xyz ) {
		if ( unrecognized_atoms_.size() > 5000 ) {
			utility_exit_with_message("can't handle more than 5000 atoms worth of unknown residues\n");
		}
		TR << "remember unrecognized atom " << rinfo_index << " " << rinfo.rosetta_resName() << " " << stripped_whitespace(iter.first)
			<< " temp " << rtemp.find(iter.first)->second << std::endl;

		core::pose::UnrecognizedAtomRecord ua( rinfo_index, rinfo.rosetta_resName(), stripped_whitespace(iter.first), iter.second, rtemp.find(iter.first)->second );

		unrecognized_atoms_.push_back( ua );
	}
}


/// @brief Query the ResidueTypeSet using the residue-type-finder for a potential match based on the
/// name3 that has possibly been remapped by Rosetta through the Nomenclature manager.
bool
PoseFromSFRBuilder::is_residue_type_recognized(
	Size const pdb_residue_index,
	std::string const & rosetta_residue_name3,
	bool & is_chemical_component_ligand
){
	using namespace core::chemical;

	// UNL cannot be recognized -- as an "unknown ligand" it has zero atoms.
	// AMW TODO: if at some point a Rosetta type is developed that is supposed to
	// represent an unknown ligand -- some kind of repulsive ball, who knows! -- undo this.
	if ( rosetta_residue_name3 == "UNL" ) return false;

	// this residue list is only used to see if there are any residue_types with name3 at all:
	ResidueTypeCOPs rsd_type_list;
	ResidueTypeCOP  rsd_type = ResidueTypeFinder( *residue_type_set_ ).name3( rosetta_residue_name3 ).get_representative_type();
	if ( rsd_type != nullptr ) rsd_type_list.push_back( rsd_type );
	if ( rsd_type && rsd_type->name().substr(0,4) == "pdb_"  && !rsd_type->is_polymer() ) {
		is_chemical_component_ligand = true;
	}
	return is_residue_type_recognized( pdb_residue_index, rosetta_residue_name3, rsd_type_list );
}

/// @brief Given an amino acid three-letter code, retrive a
/// residue. If the residue is an alpha- or beta-amino acid, also
/// retrieve its mirrored type.  Return true if the three-letter
/// codes are the same, false if they differ.  If no resiude type
/// could be loaded or there is no mirror type, throw if
/// error_on_unrecognized is true, and return false otherwise.
/// @author Vikram K. Mulliga (vmulligan@flatironinstitute.org).
bool
PoseFromSFRBuilder::d_l_threeletter_codes_are_same_for_aa(
	std::string const & name3,
	bool const error_on_unrecognized
) const {
	std::string const errmsg( "Error in core::io::pose_from_sfr::PoseFromSFRBuilder::d_l_threeletter_codes_are_same_for_aa(): " );

	//Get the restype:
	core::chemical::ResidueTypeCOP restype(
		core::chemical::ResidueTypeFinder( *residue_type_set_ ).name3(name3).get_representative_type()
	);
	if ( restype == nullptr ) {
		if ( error_on_unrecognized ) {
			utility_exit_with_message( errmsg + "Did not recognize residue with three-letter code \"" + name3 + "\"." );
		} else {
			return false;
		}
	}

	if ( restype->is_achiral_backbone() ) {
		if ( error_on_unrecognized ) {
			utility_exit_with_message( errmsg + "Residue \"" + name3 + "\" is achiral.  Cannot get mirrored type." );
		} else {
			return false;
		}
	}
	if ( !( restype->is_alpha_aa() || restype->is_beta_aa() || restype->is_gamma_aa() ) ) {
		if ( error_on_unrecognized ) {
			utility_exit_with_message( errmsg + "Residue \"" + name3 + "\" is not an alpha amino acid, beta amino acid, or gamma amino acid.  Cannot perform check." );
		} else {
			return false;
		}
	}

	//Get the mirrored type:
	core::chemical::ResidueTypeCOP mirrored_restype( residue_type_set_->get_mirrored_type( restype ) );
	if ( mirrored_restype == nullptr ) {
		if ( error_on_unrecognized ) {
			utility_exit_with_message( errmsg + "Residue \"" + name3 + "\" seems not to have a mirrored type." );
		} else {
			return false;
		}
	}

	//Compare the names:
	return (restype->name3() == mirrored_restype->name3());
}

///////////////////////////////////////////////////////////////////////
// Use ResidueTypeFinder to efficiently figure out best match
//    residue_type to these PDB atom_names, name3, etc.
chemical::ResidueTypeCOP
PoseFromSFRBuilder::get_rsd_type(
	std::string const & rosetta_residue_name3,
	Size seqpos,
	utility::vector1< std::string > const & known_connect_atoms_on_this_residue,
	std::string const & resid,
	bool const is_lower_terminus,
	bool const is_upper_terminus,
	bool const is_d_aa,
	bool const is_l_aa )
{
	// AMW: Just changed so that protein residue types don't get patched for
	// links to their polymer main chain (i.e., deprecates the weird 'neutral'
	// patches that would make residues think they had no 'upper' but DID have a
	// connection to C -- by the way, this should eventually make us not need to
	// delete LINK records for NCAA main-chains...)

	// you can be neither but not both
	debug_assert( ! ( is_d_aa && is_l_aa ) );

	std::map< std::string, Vector > const & xyz( rinfos_[ seqpos ].xyz() );

	using namespace core::chemical;
	using utility::tools::make_vector1;
	using utility::vector1;

	vector1< ResidueProperty > preferred_properties, discouraged_properties;
	vector1< VariantType > variants, disallow_variants;

	std::string residue_base_name( "" );  // used when we have more information then just a 3-letter code
	if ( sfr_.residue_type_base_names().count( resid ) ) {
		residue_base_name = sfr_.residue_type_base_names()[ resid ].second;
	}
	vector1< std::string > patch_names;

	// NOTE: Previous versions attempted to get a representative type and base which properties/variants they want based on that
	// -- Don't do that, as the representative type might be substantially different than what's found for the full search
	// Let the ResidueTypeFinder do the search for you.

	if ( is_lower_terminus ) {
		if ( known_connect_atoms_on_this_residue.contains( "P" ) || known_connect_atoms_on_this_residue.contains( "N" ) ) {
			variants.push_back( CUTPOINT_UPPER );
		} else {
			preferred_properties.push_back( LOWER_TERMINUS );
		}
	} else {
		discouraged_properties.push_back( LOWER_TERMINUS );
	}
	if ( is_upper_terminus ) {
		if ( known_connect_atoms_on_this_residue.contains( "O3'" ) || known_connect_atoms_on_this_residue.contains( "C" ) ) {
			variants.push_back( CUTPOINT_LOWER );
		} else {
			preferred_properties.push_back( UPPER_TERMINUS );
		}
	} else {
		discouraged_properties.push_back( UPPER_TERMINUS );
	}
	if ( is_d_aa ) {
		preferred_properties.push_back( D_AA );
	}
	if ( is_l_aa ) {
		preferred_properties.push_back( L_AA );
	}

	if ( rosetta_residue_name3 != "CYD" ) {
		discouraged_properties.push_back( DISULFIDE_BONDED );
	}

	if ( !options_.keep_input_protonation_state() ) {
		disallow_variants.push_back( PROTONATED );
		disallow_variants.push_back( DEPROTONATED );
	}

	utility::vector1< std::string > xyz_atom_names;
	for ( auto const & xyz_elem : xyz ) {
		xyz_atom_names.push_back( xyz_elem.first );
	}

	ResidueTypeCOP rsd_type = ResidueTypeFinder( *residue_type_set_ )
		.name3( rosetta_residue_name3 )
		.residue_base_name( residue_base_name )
		.variants( variants )
		.disallow_variants( disallow_variants )
		.preferred_properties( preferred_properties )
		.discouraged_properties( discouraged_properties )
		.patch_names( patch_names )
		.ignore_atom_named_H( is_lower_terminus )
		.check_nucleic_acid_virtual_phosphates( true )
		.connect_atoms( known_connect_atoms_on_this_residue )
		.get_best_match_residue_type_for_atom_names( xyz_atom_names );

	return rsd_type;
}

/// @brief Returns the input resid if it's either the first residue or it's
/// the second residue and the first residue has a mrb_merge_w_next behavior specified
Size PoseFromSFRBuilder::prev_residue_skipping_merges( Size resid ) const
{
	return resid > 1 ? ( merge_behaviors_[ resid-1 ] != chemical::io::mrb_merge_w_next ? resid-1 : ( resid-1 > 1 ? resid-2 : resid )) : resid;
}

/// @brief Returns the input resid if it's the first residue that has a non-null-pointing
/// entry in the residue_types_ array.
Size PoseFromSFRBuilder::prev_residue_skipping_null_residue_types( Size resid ) const
{
	if ( resid == 1 ) { return resid; }
	for ( Size ii = resid-1; ii >= 1; --ii ) {
		if ( residue_types_[ ii ] ) { return ii; }
	}
	return resid;
}

/// @brief Returns the input resid if it's either the last residue or it's
/// the second-to-last residue and the last residue has a mrb_merge_w_prev behavior specified
Size PoseFromSFRBuilder::next_residue_skipping_merges( Size resid ) const
{
	if ( resid < rinfos_.size()  ) {
		if ( merge_behaviors_[ resid+1 ] != chemical::io::mrb_merge_w_prev ) {
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
	if ( resid == residue_types_.size() ) { return resid; }
	for ( Size ii = resid+1; ii <= residue_types_.size(); ++ii ) {
		if ( residue_types_[ ii ] ) { return ii; }
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
	if ( same_chain_prev_[ resid ] == false ) { // If we've explicitly set this to false (from default true) obey this
		return false;
	}
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
	if ( resid != res_next && same_chain_prev_[ res_next ] == false ) { return false; }
	return resid != res_next && rinfos_[resid].chainID() == rinfos_[ res_next ].chainID() &&
		rinfos_[resid].terCount() == rinfos_[res_next].terCount() && ! separate_chemical_entity;
}

bool PoseFromSFRBuilder::determine_check_Ntermini_for_this_chain( char chainID ) const
{
	std::string const & nterm_chains = options_.check_if_residues_are_Ntermini();
	return "ALL" == nterm_chains ||
		std::find( nterm_chains.begin(), nterm_chains.end(), chainID ) != nterm_chains.end();
}

bool PoseFromSFRBuilder::determine_check_Ctermini_for_this_chain( char chainID ) const
{
	std::string const & cterm_chains = options_.check_if_residues_are_Ctermini();
	return "ALL" == cterm_chains ||
		std::find( cterm_chains.begin(), cterm_chains.end(), chainID ) != cterm_chains.end();
}


/// @details modifies the remapped_atom_names_ bimap at position seqpos given the names that
/// are in residue_types_[ seqpos ]
void PoseFromSFRBuilder::fill_name_map( Size seqpos )
{
	debug_assert( residue_types_[ seqpos ] );

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
		for ( auto
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
	utility::vector1< core::conformation::ResidueOP > const & context,
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

	for ( core::Size ii=1, nres=context.size(); ii<=nres; ++ii ) { //Loop through all residues already added, looking for possible residues to root the metal onto.
		core::conformation::Residue const & context_res = *context[ii];
		if ( !context_res.is_protein() ) continue; //I'm not interested in tethering metals to non-protein residues.
		if ( !context_res.has("CA") ) continue; //I'll be basing this on metal-alpha carbon distance, so anything without an alpha carbon won't get to be the root.

		numeric::xyzVector < core::Real > const residue_xyz = context_res.xyz("CA");

		core::Real const current_dist_sq = residue_xyz.distance_squared(metal_xyz);

		if ( closest_residue==0 || current_dist_sq < closest_dist_sq ) {
			closest_residue = ii;
			closest_dist_sq = (core::Size)current_dist_sq;
		}
		if ( context_res.is_metalbinding() &&
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

			if ( rinfos.size() == 0 || !rinfos.back().matches( ai ) ) {
				rinfos.emplace_back(ai); // Creates a new item on back
			}
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
