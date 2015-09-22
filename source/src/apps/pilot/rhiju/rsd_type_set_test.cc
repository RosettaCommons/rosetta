// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/io/pdb/file_data.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/import_pose_options.hh>
#include <core/import_pose/import_pose_options.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <devel/init.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <basic/Tracer.hh>

//pose_frompdb headers
#include <core/io/pdb/Field.hh>
#include <core/io/pdb/HeaderInformation.hh>
#include <core/io/pdb/file_data.hh>

// Package headers
#include <core/io/pdb/pose_io.hh>
#include <core/io/pdb/file_data_options.hh>
#include <core/io/pdb/file_data_fixup.hh>
#include <core/io/pdb/pdb_dynamic_reader.hh>
#include <core/io/pdb/pdb_dynamic_reader_options.hh>
#include <core/io/pdb/NomenclatureManager.hh>

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
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/util.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/pose/PDBInfo.hh>

#include <core/pose/util.hh>
#include <core/pose/ncbb/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/cryst/util.hh>

// Basic headers
#include <basic/options/option.hh>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include <utility/io/izstream.hh>


#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

//silly using/typedef
// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>


using namespace core;
using namespace ObjexxFCL::format;
using namespace basic::options::OptionKeys;

using utility::vector1;

static THREAD_LOCAL basic::Tracer tr( "rsdtypeset_test" );

using ObjexxFCL::strip_whitespace;
using ObjexxFCL::stripped_whitespace;
using ObjexxFCL::rstripped_whitespace;
using namespace ObjexxFCL::format;
using namespace core::io::pdb;
using namespace core::import_pose;

using std::string;
using std::iostream;

using namespace core::chemical;

OPT_KEY( Boolean, sequence_test )
OPT_KEY( Boolean, original_rsd_type_set_map )
OPT_KEY( Boolean, original_rsd_type_set_map_only )

//////////////////////////////////////////////////////////////////////////////////
//
// sandbox for testing new rsd_type_set "on-the-fly" functions
//
// move this into core/chemical/ResidueTypeSet, etc. when ready.
//
// DELETE THIS FROM PILOT APPS after 2015, if ResidueTypeSet is finalized by then.
//
//   -- rhiju, june 2015
//
//////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Testing new functions for ResidueTypeSet
///////////////////////////////////////////////////////////////////////////////
void
residue_types_from_sequence_test()
{
	using namespace core::chemical;
	using namespace core::pose;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ResidueTypeSetCOP rsd_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

	vector1< std::string > sequences;
	if ( option[ in::file::seq ].user() ) {
		sequences = option[ in::file::seq ]();
	} else  {
		sequences.push_back( "K[LYS_p:NtermProteinFull]QFTKCTMFHTSGY[TYR_p:C_methylamidated]T[THR_p:N_acetylated]QAIVEYGLFQIS[SER_p:CtermProteinFull]" );
		sequences.push_back( "g[RGU:Virtual_Phosphate]ggcuu[URA:rna_cutpoint_lower]c[RCY:rna_cutpoint_upper]ggccu" );
	}

	bool auto_termini( true );
	for ( Size m = 1; m <= sequences.size(); m++ ) {
		std::string const sequence = sequences[ m ];

		// this shows up in pose_from_sequence
		ResidueTypeCOPs rsd_types  = residue_types_from_sequence ( sequence, *rsd_set, auto_termini );
		for ( Size n = 1; n <= rsd_types.size(); n++ ) {
			tr << rsd_types[ n ]->name() << std::endl;
		}
		std::cout << std::endl;
	}

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
// the original code (probably could keep this in file_data.cc, with a note to DELETE AFTER 2015).
ResidueTypeCOP
get_rsd_type_ORIGINAL(
	Size const i,
	std::string const & name3,
	std::map< std::string, Vector > const & xyz,
	std::map< std::string, double > const & rtemp,
	utility::vector1<Size> & UA_res_nums,
	utility::vector1<std::string> & UA_res_names,
	utility::vector1<std::string> & UA_atom_names,
	utility::vector1<numeric::xyzVector<Real> > & UA_coords,
	utility::vector1<core::Real> & UA_temps,
	FileDataOptions const & options,
	chemical::ResidueTypeSet const & residue_set,
	FileData & fd,
	utility::vector1< std::string > const &  branch_points_on_this_residue,
	std::string const & resid,
	bool const is_lower_terminus,
	bool const is_upper_terminus,
	bool const is_branch_point,
	bool const is_branch_lower_terminus,
	bool & last_residue_was_recognized )
{


	using namespace core::chemical;

		// Get a list of ResidueTypes that could apply for this particular 3-letter PDB residue name.
		ResidueTypeCOPs const & rsd_type_list( residue_set.name3_map_DO_NOT_USE( name3 ) );
		if ( ! is_residue_type_recognized(
				i, name3, rsd_type_list, xyz, rtemp,
				UA_res_nums, UA_res_names, UA_atom_names, UA_coords, UA_temps, options)) {
			last_residue_was_recognized = false;
			return 0;
		}

		if ( tr.Debug.visible() ) {
			tr.Debug << "Residue " << i << std::endl;
		}
		if ( tr.Trace.visible() ) {
			// tr.Trace << "...same_chain_prev: " << same_chain_prev << std::endl;
			// tr.Trace << "...same_chain_next: " << same_chain_next << std::endl;
			tr.Trace << "...is_lower_terminus: " << is_lower_terminus << std::endl;
			//tr.Trace << "...check_Ntermini_for_this_chain: "<< check_Ntermini_for_this_chain << std::endl;
			tr.Trace << "...is_upper_terminus: " << is_upper_terminus << std::endl;
			//tr.Trace << "...check_Ctermini_for_this_chain: "<< check_Ctermini_for_this_chain << std::endl;
			tr.Trace << "...is_branch_point: " << is_branch_point << std::endl;
			tr.Trace << "...is_branch_lower_terminus: "<< is_branch_lower_terminus << std::endl;
			tr.Trace << "...last_residue_was_recognized: " << last_residue_was_recognized << std::endl;
		}

		ResidueTypeCOPs filtered_rsd_type_list;

		for ( Size j=1; j<= rsd_type_list.size(); ++j ) {
			ResidueType const & rsd_type( *(rsd_type_list[j]) );
			bool const is_polymer( rsd_type.is_polymer() ); // need an example residue type, though this will
			// remain fixed for all residue_types with the same name3

			// only take the desired variants
			bool lower_term_type = rsd_type.has_variant_type( LOWER_TERMINUS_VARIANT ) ||
					rsd_type.has_variant_type( LOWERTERM_TRUNC_VARIANT );
			bool upper_term_type = rsd_type.has_variant_type( UPPER_TERMINUS_VARIANT ) ||
					rsd_type.has_variant_type( UPPERTERM_TRUNC_VARIANT );
			if ( is_polymer && (
				(is_lower_terminus != lower_term_type ) || (is_upper_terminus != upper_term_type ) ) ) {
				if ( tr.Trace.visible() ) {
					tr.Trace << "Discarding '" << rsd_type.name() << "' ResidueType" << std::endl;
					tr.Trace << "because of the terminus state" << std::endl;
				}
				continue;
			}
			if ( is_polymer && ( is_branch_point != rsd_type.is_branch_point() ) ) {
				if ( tr.Trace.visible() ) {
					tr.Trace << "Discarding '" << rsd_type.name() << "' ResidueType" << std::endl;
					tr.Trace << "because of the branch state" << std::endl;
				}
				continue;
			}
			if ( is_polymer && ( is_branch_lower_terminus != rsd_type.is_branch_lower_terminus() ) ) {
				if ( tr.Trace.visible() ) {
					tr.Trace << "Discarding '" << rsd_type.name() << "' ResidueType" << std::endl;
					tr.Trace << "because of the branch lower terminus state" << std::endl;
				}
				continue;
			}
			// Okay, this logic is NOT OBVIOUS, so I am going to add my explanation.
			// We DEFINITELY want to assign a disulfide type from the start if the PDB just up and says CYD. That's great!
			// But if we DO NOT see CYD, we do not want to assign a disulfide type. We want disulfide connections
			// to be inferred later on, in conformation's detect_disulfides as called in the pose-building process.
			// Commenting out this logic causes anything that lacks a HG (i.e. crystal structures) to be assigned as the
			// disulfide type instead of the CYS type--which MIGHT be right, but might be wrong and leads to wasteful
			// disulfide reversion.
			if ( rsd_type.aa() == aa_cys && rsd_type.has_variant_type( DISULFIDE ) && name3 != "CYD" ) {
				if ( tr.Trace.visible() ) {
					tr.Trace << "Discarding '" << rsd_type.name() << "' ResidueType" << std::endl;
					tr.Trace << "because of the disulfide state" << std::endl;
				}
				continue;
			}
			if ( ! options.keep_input_protonation_state() &&
					( rsd_type.has_variant_type( PROTONATED ) || rsd_type.has_variant_type( DEPROTONATED ) ) ) {
				if ( tr.Trace.visible() ) {
					tr.Trace << "Discarding '" << rsd_type.name() << "' ResidueType" << std::endl;
					tr.Trace << "because of the protonation state" << std::endl;
				}
				continue;
			}

			// special checks to ensure selecting the proper carbohydrate ResidueType
			if ( rsd_type.is_carbohydrate() &&
					 residue_type_base_name( rsd_type ) != fd.residue_type_base_names[ resid ] ) {
				if ( tr.Trace.visible() ) {
					tr.Trace << "Discarding '" << rsd_type.name() << "' ResidueType" << std::endl;
					tr.Trace << "because the residue is not a carbohydrate" << std::endl;
				}
				continue;
			}
			if ( rsd_type.is_carbohydrate() && rsd_type.is_branch_point() ) {
				// The below assumes that ResidueTypes with fewer patches are selected 1st, that is, that an
				// :->2-branch ResidueType will be checked as a possible match before an :->2-branch:->6-branch
				// ResidueType.  If this were not the case, Rosetta could misassign an :->2-branch:->6-branch
				// ResidueType to a residue that actually only has a single branch at the 2 or 6 position.
				char branch_point;
				bool branch_point_is_missing( false );
				Size const n_branch_points( branch_points_on_this_residue.size() );
				for ( core::uint k( 1 ); k <= n_branch_points; ++k ) {
					branch_point = branch_points_on_this_residue[ k ][ 2 ];  // 3rd column (index 2) is the atom number.
					if ( tr.Debug.visible() ) {
						tr.Debug << "Checking '" << rsd_type.name() <<
								"' for branch at position " << branch_point << std::endl;
					}
					if ( residue_type_all_patches_name( rsd_type ).find( string( 1, branch_point ) + ")-branch" ) ==
							string::npos ) {
						branch_point_is_missing = true;
						break;
					}
				}
				if ( branch_point_is_missing ) {
					if ( tr.Trace.visible() ) {
						tr.Trace << "Discarding '" << rsd_type.name() << "' ResidueType" << std::endl;
						tr.Trace << "because of a missing branch point" << std::endl;
					}
					continue;
				}
			}

			if ( tr.Debug.visible() ) {
				tr.Debug << "Trying '" << rsd_type.name() << "' ResidueType" << std::endl;
			}

			filtered_rsd_type_list.push_back( rsd_type_list[j] );
		}


		typedef std::map< std::string, Vector > ResidueCoords;
		utility::vector1< std::string > xyz_atom_names;
		for ( ResidueCoords::const_iterator iter=xyz.begin(), iter_end=xyz.end(); iter!= iter_end; ++iter ) {
			std::string xyz_name = iter->first;
			xyz_atom_names.push_back( strip_whitespace( xyz_name ) );
		}

		ResidueTypeCOP best_match_rsd_type = find_best_match( filtered_rsd_type_list, xyz_atom_names,
																													is_lower_terminus /* ignore_atom_named_H */ );


		return best_match_rsd_type;
}

////////////////////////////////////////////////////////////////////////
// This is what I want to introduce into build_pose_as_is1
////////////////////////////////////////////////////////////////////////
ResidueTypeCOP
get_rsd_type_TEST(
	Size const i,
	std::string const & name3,
	std::map< std::string, Vector > const & xyz,
	std::map< std::string, double > const & rtemp,
	utility::vector1<Size> & UA_res_nums,
	utility::vector1<std::string> & UA_res_names,
	utility::vector1<std::string> & UA_atom_names,
	utility::vector1<numeric::xyzVector<Real> > & UA_coords,
	utility::vector1<core::Real> & UA_temps,
	FileDataOptions const & options,
	chemical::ResidueTypeSet const & residue_set,
	FileData & fd,
	utility::vector1< std::string > const &  branch_points_on_this_residue,
	std::string const & resid,
	bool const is_lower_terminus,
	bool const is_upper_terminus,
	bool const is_branch_point,
	bool const is_branch_lower_terminus,
	bool & last_residue_was_recognized )
{
	typedef std::map< std::string, Vector > ResidueCoords;

	using namespace core::chemical;
	using utility::tools::make_vector1;

	vector1< vector1< VariantType > > required_variants_in_sets;
	vector1< ResidueProperty > properties;
	vector1< VariantType >     variants, disallow_variants;
	std::string                residue_base_name( "" ); // carbohydrates
	vector1< std::string >     patch_names;

	ResidueTypeCOP rsd_type = ResidueTypeFinder( residue_set ).name3( name3 ).get_representative_type();
	ResidueTypeCOPs rsd_types;
	if ( rsd_type != 0 ) rsd_types.push_back( rsd_type );

	// hold over -- probably should not include in this function, and refactor rsd_types --> found_residue_type
	// But -- need to be careful about HOH logic.
	// Should rename the function too.
	if ( ! is_residue_type_recognized( i, name3, rsd_types, xyz, rtemp,
																		 UA_res_nums, UA_res_names, UA_atom_names, UA_coords, UA_temps, options)) {
		last_residue_was_recognized = false;
		return 0;
	}

	if ( rsd_type->is_polymer() ) {
		if ( is_lower_terminus )     required_variants_in_sets.push_back( make_vector1( LOWER_TERMINUS_VARIANT, LOWERTERM_TRUNC_VARIANT ) );
		if ( is_upper_terminus )     required_variants_in_sets.push_back( make_vector1( UPPER_TERMINUS_VARIANT, UPPERTERM_TRUNC_VARIANT ) );
		if ( is_branch_point && !rsd_type->is_carbohydrate() )  properties.push_back( BRANCH_POINT ); // note that carbohydrates are covered by patch_names below.
		if ( is_branch_lower_terminus ) properties.push_back( BRANCH_LOWER_TERMINUS );
	}
	if ( rsd_type->aa() == aa_cys && name3 != "CYD" ) disallow_variants.push_back( DISULFIDE );
	if ( !options.keep_input_protonation_state() ) {
		disallow_variants.push_back( PROTONATED );
		disallow_variants.push_back( DEPROTONATED );
	}
	if ( rsd_type->is_carbohydrate() ) {
		runtime_assert( fd.residue_type_base_names[ resid ].substr(0,3) == name3 ); // is this true for carbohydrates? Or Mg(2+)?
		residue_base_name = fd.residue_type_base_names[ resid ];
		// The below assumes that ResidueTypes with fewer patches are selected 1st, that is, that an
		// :->2-branch ResidueType will be checked as a possible match before an :->2-branch:->6-branch
		// ResidueType.  If this were not the case, Rosetta could misassign an :->2-branch:->6-branch
		// ResidueType to a residue that actually only has a single branch at the 2 or 6 position.
		for ( core::uint k( 1 ); k <= branch_points_on_this_residue.size(); ++k ) {
			char const & branch_point = branch_points_on_this_residue[ k ][ 2 ];  // 3rd column (index 2) is the atom number.
			patch_names.push_back( "->" + string( 1, branch_point ) + ")-branch" );
		}
	}

	for ( Size k = 1; k <= required_variants_in_sets.size(); k++ ){
		vector1< VariantType > const & variant_set = required_variants_in_sets[ k ];
		for ( Size m = 1; m <= variant_set.size(); m++ ) {
			variants.push_back( variant_set[ m ] );
		}
	}

	utility::vector1< std::string > xyz_atom_names;
	for ( ResidueCoords::const_iterator iter=xyz.begin(), iter_end=xyz.end(); iter!= iter_end; ++iter ) {
		std::string xyz_name = iter->first;
		xyz_atom_names.push_back( strip_whitespace( xyz_name ) );
	}

	// following 'chaining' looks a lot like chemical::ResidueSelector -- use that instead? Rename to chemical::ResidueTypeSelector.
	rsd_type = ResidueTypeFinder( residue_set ).name3( name3 ).residue_base_name( residue_base_name ).disallow_variants( disallow_variants ).variants_in_sets( required_variants_in_sets ).properties( properties ).patch_names( patch_names ).ignore_atom_named_H( is_lower_terminus ).get_best_match_residue_type_for_atom_names( xyz_atom_names );

  return rsd_type;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fills the pose with the data from FileData
// huge chunk of code copied from file_data.cc to allow tests of new rsd_type identifier function.
void
build_pose_as_is1_TEST(
	io::pdb::FileData & fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	id::AtomID_Mask & missing,
	FileDataOptions const & options
)
{
	typedef std::map< std::string, double > ResidueTemps;
	typedef std::map< std::string, Vector > ResidueCoords;
	typedef utility::vector1< std::string > Strings;

	using namespace chemical;
	using namespace conformation;

	// reset current data
	pose.clear();

	utility::vector1< ResidueInformation > rinfos;
	id::NamedAtomID_Mask coordinates_assigned( false );
	// Map pose residue numbers to indices into rinfos.
	// Some residues in the input file may be discarded (missing atoms, unrecognized, etc.)
	utility::vector1< Size > pose_to_rinfo;
	fd.create_working_data( rinfos, options );
	fixup_rinfo_based_on_residue_type_set( rinfos, residue_set );
	utility::vector1<ResidueTemps> pose_temps;

	Strings branch_lower_termini;

	Size const nres_pdb( rinfos.size() );

	// Map rinfo atom names to Rosetta pose atom names (and pose->rinfo for the right map)
	utility::vector1< core::io::pdb::NameBimap > rinfo_name_map(nres_pdb);

	utility::vector1<Size> UA_res_nums;
	utility::vector1<std::string> UA_res_names, UA_atom_names;
	utility::vector1<numeric::xyzVector<Real> > UA_coords;
	utility::vector1<core::Real> UA_temps;

	std::string chains_whose_residues_are_separate_chemical_entities =
			options.chains_whose_residues_are_separate_chemical_entities();
	std::string::const_iterator const entities_begin = chains_whose_residues_are_separate_chemical_entities.begin();
	std::string::const_iterator const entities_end = chains_whose_residues_are_separate_chemical_entities.end();

	std::string chains_to_check_if_Ntermini= options.check_if_residues_are_Ntermini() ;
	std::string::const_iterator const check_Ntermini_begin = chains_to_check_if_Ntermini.begin();
	std::string::const_iterator const check_Ntermini_end = chains_to_check_if_Ntermini.end();

	std::string chains_to_check_if_Ctermini= options.check_if_residues_are_Ctermini() ;
	std::string::const_iterator const check_Ctermini_begin = chains_to_check_if_Ctermini.begin();
	std::string::const_iterator const check_Ctermini_end = chains_to_check_if_Ctermini.end();

	//mjo do not add residue by bond if the last residue was not recognized
	bool last_residue_was_recognized(true);

	// Loop over every residue in the FileData extracted from the PDB file, select appropriate ResidueTypes,
	// create Residues, and build the Pose.
	for ( Size i = 1; i <= nres_pdb; ++i ) {
		ResidueInformation const & rinfo = rinfos[ i ];
		char chainID = rinfo.chainID;
		std::string const & pdb_name = rinfo.resName;
		std::string const & resid = rinfo.resid();

		runtime_assert( resid.size() == 6 );

		// Convert PDB 3-letter code to Rosetta 3-letter code, if a list of alternative codes has been provided.
		std::pair< std::string, std::string > const & rosetta_names(
				NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( pdb_name ) );
		std::string const & name3( rosetta_names.first );
		if ( rosetta_names.second != "" ) {
			fd.residue_type_base_names[ resid ] = rosetta_names.second;
		}

		bool const separate_chemical_entity = find(entities_begin, entities_end, chainID ) !=  entities_end;
		bool const same_chain_prev = ( i > 1        && chainID == rinfos[i-1].chainID &&
				rinfo.terCount == rinfos[i-1].terCount && !separate_chemical_entity);
		bool const same_chain_next = ( i < nres_pdb && chainID == rinfos[i+1].chainID &&
				rinfo.terCount == rinfos[i+1].terCount && !separate_chemical_entity);
		bool const check_Ntermini_for_this_chain = ("ALL" == chains_to_check_if_Ntermini) ?
				true : find(check_Ntermini_begin, check_Ntermini_end, chainID ) ==  check_Ntermini_end;
		bool const check_Ctermini_for_this_chain = ("ALL" == chains_to_check_if_Ctermini) ?
				true : find(check_Ctermini_begin, check_Ctermini_end, chainID ) ==  check_Ctermini_end;

		// Determine polymer information: termini, branch points, etc.
		Strings branch_points_on_this_residue;
		bool is_branch_point( false );
		if ( fd.link_map.count( resid ) ) {  // if found in the linkage map
			// Find and store to access later:
			//     - associated 1st residue of all branches off this residue (determines branch lower termini)
			//     - positions of branch points
			for ( Size branch = 1, n_branches = fd.link_map[ resid ].size(); branch <= n_branches; ++branch ) {
				LinkInformation const & link_info( fd.link_map[ resid ][ branch ] );
				if ( link_info.chainID1 == link_info.chainID2 && link_info.resSeq1 == ( link_info.resSeq2 - 1 ) ) {
					// If this occurs, the link is to the next residue on the same chain, so both residues are part of
					// the same main chain or branch, and this linkage information can be ignored.
					// Note that this assumes insertion codes are not involved! It also assumes that the PDB file
					// makers did things reasonably.
					continue;
				}
				is_branch_point = true;
				branch_lower_termini.push_back( link_info.resID2 );
				branch_points_on_this_residue.push_back( link_info.name1 );
			}
		}
		bool const is_branch_lower_terminus = branch_lower_termini.contains(resid);
		bool const is_lower_terminus( ( i == 1 || rinfos.empty() || (!same_chain_prev && !is_branch_lower_terminus) )
				&& check_Ntermini_for_this_chain );
		bool const is_upper_terminus( ( i == nres_pdb || !same_chain_next ) && check_Ctermini_for_this_chain );


		ResidueCoords const & xyz = rinfo.xyz;
		ResidueTemps  const & rtemp = rinfo.temps;

		ResidueTypeCOP rsd_type_cop = get_rsd_type_TEST(     i, name3, xyz, rtemp, UA_res_nums, UA_res_names, UA_atom_names, UA_coords, UA_temps, options, residue_set,
																															fd, branch_points_on_this_residue, resid, is_lower_terminus, is_upper_terminus, is_branch_point,
																															is_branch_lower_terminus,  last_residue_was_recognized);

		if ( rsd_type_cop == 0 ) {
			if ( !last_residue_was_recognized ) continue;

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
																	 boost::lexical_cast<string>(i) +
																	 "\nLooking for" + variant + " residue with 3-letter code: " + name3 );
			}
		}

		ResidueType const & rsd_type = residue_set.name_map( rsd_type_cop->name() ); // ensure instantiated & saved in residue_set.

		if ( tr.Trace.visible() ) {
			// tr.Trace << "Naive match of " << rsd_type.name() << " with " << best_rsd_missing << " missing and "
			// 		<< best_xyz_missing << " discarded atoms." << std::endl;
		}


		// Map the atom names.
		fill_name_map( rinfo_name_map[i], rinfo, rsd_type, options );

		debug_assert( rsd_type.natoms() >= rinfo_name_map[i].left.size() );
		core::Size missing_atoms( rsd_type.natoms() - rinfo_name_map[i].left.size() );
		if ( missing_atoms > 0 ) {
			tr.Debug << "Match: '" << rsd_type.name() << "'; missing " << missing_atoms << " coordinates" << std::endl;
		}

		debug_assert( rinfo.xyz.size() >= rinfo_name_map[i].left.size() );
		core::Size discarded_atoms( rinfo.xyz.size() - rinfo_name_map[i].left.size() );
		if ( is_lower_terminus && rinfo.xyz.count(" H  ") && ! rinfo_name_map[i].left.count(" H  ") ) {
			// Don't worry about missing BB H if Nterm
			--discarded_atoms;
		}
		if ( discarded_atoms > 0 ) {
			tr.Warning << "[ WARNING ] discarding " << discarded_atoms
					<< " atoms at position " << i << " in file " << fd.filename
					<< ". Best match rsd_type:  " << rsd_type.name() << std::endl;
		}

		// check for missing mainchain atoms:
		if ( rsd_type.is_polymer() ) {
			AtomIndices const & mainchain( rsd_type.mainchain_atoms() );
			Size const nbb( mainchain.size() );
			if ( nbb >= 3 ) {
				bool mainchain_core_present( false );
				for ( Size k=1; k<= nbb-2; ++k ) {
					std::string const & name1(rsd_type.atom_name(mainchain[k  ]));
					std::string const & name2(rsd_type.atom_name(mainchain[k+1]));
					std::string const & name3(rsd_type.atom_name(mainchain[k+2]));
					if( !rinfo_name_map[i].right.count(name1) ||
							!rinfo_name_map[i].right.count(name2) ||
							!rinfo_name_map[i].right.count(name3) ){
						continue;
					}
					std::string const & rinfo_name1( rinfo_name_map[i].right.find( name1 )->second );
					std::string const & rinfo_name2( rinfo_name_map[i].right.find( name2 )->second );
					std::string const & rinfo_name3( rinfo_name_map[i].right.find( name3 )->second );
					if ( xyz.count( rinfo_name1 ) && xyz.count( rinfo_name2 ) && xyz.count( rinfo_name3 ) ) {
						mainchain_core_present = true;
						break;
					}
				}
				if ( !mainchain_core_present ) {
					tr.Warning << "[ WARNING ] skipping pdb residue b/c it's missing too many mainchain atoms: " <<
							resid << ' ' << name3 << ' ' << rsd_type.name() << std::endl;
					for ( Size k=1; k<= nbb; ++k ) {
						std::string const & name(rsd_type.atom_name(mainchain[k]));
						if( !rinfo_name_map[i].right.count(name) ||
								!xyz.count( rinfo_name_map[i].right.find(name)->second ) ) {
							// Use of unmapped name deliberate
							tr << "missing: " << name << std::endl;
						}
					}
					if( options.exit_if_missing_heavy_atoms() == true ) {
						utility_exit_with_message("quitting due to missing heavy atoms");
					}
					continue;
				}
			}
		}

		// found a match, create the residue...
		ResidueOP new_rsd( ResidueFactory::create_residue( rsd_type ) );

		// ...and now fill in the coords
		for ( ResidueCoords::const_iterator iter=xyz.begin(), iter_end=xyz.end(); iter!= iter_end; ++iter ) {
			std::string const & rinfo_name( iter->first );
			if ( rinfo_name_map[i].left.count( rinfo_name ) ) {
				// offsetting all coordinates by a small constant prevents problems with atoms located
				// at position (0,0,0).
				// This is a bit of a dirty hack but it fixes the major problem of reading in rosetta
				// pdbs which usually start at 0,0,0. However the magnitude of this offset is so small
				// that the output pdbs should still match input pdbs. hopefully. yes. aehm.
				// RM: I'm not sure what the problem with having coordinates exactly at the origin is.
				// RM: If we do have a problem with that, it seems it should be fixed there and not here.
				// RM: (As you could imagine theoretically hitting (0,0,0) during minimization or packing.)
				double offset = 1e-250; // coordinates now double, so we can use _really_ small offset.
				std::string const & pose_name( rinfo_name_map[i].left.find( rinfo_name )->second );
				new_rsd->atom( pose_name ).xyz( iter->second + offset );
				// +1 here as we haven't added the residue to the pose yet.
				id::NamedAtomID atom_id( pose_name, pose.total_residue()+1 );
				coordinates_assigned.set( atom_id, true);
			}
			//else runtime_assert( iter->first == " H  " && rsd_type.is_terminus() ); // special casee
		}

		check_and_correct_sister_atoms( new_rsd );

		Size const old_nres( pose.total_residue() );

		if ( tr.Trace.visible() ) {
			tr.Trace << "...new residue is a polymer: " << new_rsd->type().is_polymer() << std::endl;
			if ( old_nres >= 1 ) {
				tr.Trace << "The old residue is a polymer: " << pose.residue_type(old_nres).is_polymer() << std::endl;
			}
		}

		// Add the first new residue to the pose
		if ( !old_nres ) {
			if ( tr.Trace.visible() ) {
				tr.Trace << rsd_type.name() << " " << i << " is the start of a new pose" << std::endl;
			}
			pose.append_residue_by_bond( *new_rsd );

		} else if ( ( ( is_lower_terminus && check_Ntermini_for_this_chain ) || ! same_chain_prev )
						|| /* is_branch_lower_terminus || */
							pose.residue_type( old_nres ).has_variant_type( "C_METHYLAMIDATION" ) ||
						! new_rsd->is_polymer() ||
						! pose.residue_type( old_nres ).is_polymer() ||
						! last_residue_was_recognized ) {
			// A new chain because this is a lower terminus (see logic above for designation)
			// and if we're not checking it then it's a different chain from the previous

			core::Size rootindex=1;

			// Ensure that metal ions are connected by a jump to the closest metal-binding residue that is lower in sequence.
			if(new_rsd->is_metal() && basic::options::option[basic::options::OptionKeys::in::auto_setup_metals].user()) {
				// If this is a metal ion and we're automatically setting up metals, search for the closest metal-binding residue
				// and make that the jump parent.  Otherwise, let the jump parent be the closest residue.
				numeric::xyzVector < core::Real > const metal_xyz = new_rsd->xyz(1); //Atom 1 is always the metal of a residue representing a metal ion.  (There's a check for this in residue_io.cc).

				core::Size closest_metalbinding_residue=0;
				core::Size metalbinding_dist_sq = 0;
				core::Size closest_residue=0;
				core::Size closest_dist_sq = 0;

				for(core::Size jr=1, nres=pose.n_residue(); jr<=nres; ++jr) { //Loop through all residues already added, looking for possible residues to root the metal onto.
					if(!pose.residue(jr).is_protein()) continue; //I'm not interested in tethering metals to non-protein residues.
					if(!pose.residue(jr).has("CA")) continue; //I'll be basing this on metal-alpha carbon distance, so anything without an alpha carbon won't get to be the root.

					numeric::xyzVector < core::Real > const residue_xyz = pose.residue(jr).xyz("CA");

					core::Real const current_dist_sq = residue_xyz.distance_squared(metal_xyz);

					if(closest_residue==0 || current_dist_sq < closest_dist_sq) {
						closest_residue = jr;
						closest_dist_sq = (core::Size)current_dist_sq;
					}
					if(	pose.residue(jr).is_metalbinding() &&
								(closest_metalbinding_residue==0 || current_dist_sq < metalbinding_dist_sq)
						) {
							closest_metalbinding_residue = jr;
							metalbinding_dist_sq = (core::Size)current_dist_sq;
					}
				} //Inner loop through all residues

				if(closest_metalbinding_residue!=0) rootindex=closest_metalbinding_residue; //If we found a metal-binding residue, it's the root; otherwise, the closest residue is.
				else if(closest_residue!=0) rootindex=closest_residue;

			}	//If this is a metal

			if(rootindex>1) {tr << rsd_type.name() << " " << i << " was added by a jump, with base residue " << rootindex << std::endl;}

			pose.append_residue_by_jump( *new_rsd, rootindex /*pose.total_residue()*/ );

		} else { // Append residue to current chain dependent on bond length.
			if (!options.missing_dens_as_jump()) {
				if ( tr.Trace.visible() ) {
					tr.Trace << rsd_type.name() << " " << i << " is appended to chain " << chainID << std::endl;
				}
				pose.append_residue_by_bond( *new_rsd );
			} else {
				//fpd look for missing density in the input PDB
				//fpd if there is a bondlength > 3A
				//fpd we will consider this missing density
				Residue const &last_rsd( pose.residue( old_nres ) );
				core::Real bondlength = ( last_rsd.atom( last_rsd.upper_connect_atom() ).xyz() -
						new_rsd->atom( new_rsd->lower_connect_atom() ).xyz() ).length();

				if ( bondlength > 3.0 ) {
					tr << "[ WARNING ] missing density found at residue (rosetta number) " << old_nres << std::endl;
					pose.append_residue_by_jump( *new_rsd, old_nres );

					if (pose.residue_type(old_nres).is_protein()) {
						if (!pose.residue_type(old_nres).has_variant_type( UPPER_TERMINUS_VARIANT ) &&
						    !pose.residue_type(old_nres).has_variant_type( UPPERTERM_TRUNC_VARIANT ) )
							core::pose::add_variant_type_to_pose_residue( pose, chemical::UPPERTERM_TRUNC_VARIANT, old_nres );
					} else {
						if (!pose.residue_type(old_nres).has_variant_type( UPPER_TERMINUS_VARIANT ))
							core::pose::add_variant_type_to_pose_residue( pose, chemical::UPPER_TERMINUS_VARIANT, old_nres );
					}

					if (pose.residue_type(old_nres+1).is_protein()) {
						core::pose::add_variant_type_to_pose_residue( pose, chemical::LOWERTERM_TRUNC_VARIANT, old_nres+1 );
					} else {
						core::pose::add_variant_type_to_pose_residue( pose, chemical::LOWER_TERMINUS_VARIANT, old_nres+1 );
					}
				} else {
					if ( tr.Trace.visible() ) {
						tr.Trace << rsd_type.name() << " " << i << " is appended to chain" << chainID << std::endl;
					}
					pose.append_residue_by_bond( *new_rsd );
				}
			}
		}

		// If newly added residue was a carbohydrate, set flag on conformation.
		if ( new_rsd->is_carbohydrate() ) {
			pose.conformation().contains_carbohydrate_residues( true );
		}

		pose_to_rinfo.push_back( Size(i) );
		pose_temps.push_back( rinfo.temps );

		// Update the pose-internal chain label if necessary.
		if ( ( is_lower_terminus || ! check_Ntermini_for_this_chain || is_branch_lower_terminus ) &&
				pose.total_residue() > 1 ) {
			pose.conformation().insert_chain_ending( pose.total_residue() - 1 );
		}

		last_residue_was_recognized = true;
	} // i=1,nres_pdb


	// Check termini status of newly created pose residues.
	// RM: All considered, this is a poor place to do this - we should ideally be doing this back
	// when we're originally doing the typing - though some knowledge of downstream residues is necessary.
	Size const nres( pose.total_residue() );
	for ( Size i=1; i<= nres; ++i ) {
		// Need to map pose index to rinfo index, in case we're skipping residues
		ResidueInformation const & rinfo = rinfos[pose_to_rinfo[i]];
		char chainID = rinfo.chainID;

		bool const check_Ntermini_for_this_chain = ("ALL" == chains_to_check_if_Ntermini) ?
					true : find(check_Ntermini_begin, check_Ntermini_end, chainID ) ==  check_Ntermini_end;
		bool const check_Ctermini_for_this_chain = ("ALL" == chains_to_check_if_Ctermini) ?
					true : find(check_Ctermini_begin, check_Ctermini_end, chainID ) ==  check_Ctermini_end;

		if ( !check_Ntermini_for_this_chain ) { continue; }
		if ( !check_Ctermini_for_this_chain ) { continue; }

		//Residue const & rsd( pose.residue( i ) ); // THIS WAS A BAD BUG
		bool type_changed(false);
		if ( !pose.residue_type(i).is_polymer() ) { continue; }
		if ( !pose.residue_type(i).is_lower_terminus() &&
				( i == 1 ||
				!pose.residue_type( i-1 ).is_polymer() ||
				(pose.residue_type( i-1 ).is_upper_terminus() &&
						!pose.residue_type( i ).is_branch_lower_terminus() ) ) ) {
			tr << "Adding undetected lower terminus type to residue " << i << std::endl;
			core::pose::add_lower_terminus_type_to_pose_residue( pose, i );
			type_changed = true;
		}
		if ( !pose.residue_type(i).is_upper_terminus() &&
				( i == nres ||
				!pose.residue_type(i+1).is_polymer() ||
				pose.residue_type(i+1).is_lower_terminus() /*||
				pose.residue_type(i+1).has_variant_type(BRANCH_LOWER_TERMINUS_VARIANT)*/ ) ) {
			tr << "Adding undetected upper terminus type to residue " << i << std::endl;
			core::pose::add_upper_terminus_type_to_pose_residue( pose, i );
			type_changed = true;
		}
		if( type_changed ) {
			// add_terminus_type will copy coordinates for matching atoms - see if there's additional atoms we missed.
			for( core::Size ii(1); ii <= pose.residue_type(i).natoms(); ++ii ) {
				std::string const & name( pose.residue_type(i).atom_name(ii) );
				id::NamedAtomID atom_id( name, i );
				// Unfortunately we're doing only exact name matches here.
				if( ! coordinates_assigned[atom_id] && rinfo.xyz.count(name)  ) {
					pose.set_xyz( atom_id, rinfo.xyz.find(name)->second );
					coordinates_assigned.set(atom_id, true);
					rinfo_name_map[pose_to_rinfo[i]].insert( NameBimap::value_type( name, name ) );
					tr.Debug << "Setting coordinates for undetected atom " << name << " on residue " << i << std::endl;
				}
			}
		}
	}


	// now handle missing atoms
	//id::AtomID_Mask missing( false );
	Size num_heavy_missing = 0;

	core::pose::initialize_atomid_map( missing, pose ); // dimension the missing-atom mask
	if ( pose.total_residue() == 0 ) {
		// if unchecked it segfaults further down...

		// PDBInfo setup
		core::pose::PDBInfoOP pdb_info( new core::pose::PDBInfo( pose.total_residue() ) );
		for( Size i = 1; i <= UA_res_nums.size(); ++i ) {
			pdb_info->add_unrecognized_atom( UA_res_nums[i], UA_res_names[i], UA_atom_names[i], UA_coords[i], UA_temps[i] );
		}
		// store pdb info
		pose.pdb_info( pdb_info );
		return;
	}

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			id::AtomID atom_id( j, i );
			id::NamedAtomID named_atom_id( rsd.atom_name(j), i );
			if ( ! coordinates_assigned[named_atom_id] ) {
				missing[ atom_id ] = true;
				if( !rsd.atom_is_hydrogen(j) ) num_heavy_missing++;
			}
		}
	}

	pose.conformation().fill_missing_atoms( missing );

	//ja save the pdb residue indices in the Pose //well, PDBInfo
	//ja pdb residue indices can be negative
	utility::vector1< int > pdb_numbering;
	//sml chain char
	utility::vector1< char > pdb_chains, insertion_codes;
	//Size const nres( pose.total_residue() );
	for ( Size i(1); i <= nres; ++i ) {
		ResidueInformation const & rinfo = rinfos[pose_to_rinfo[i]];
		pdb_numbering.push_back( rinfo.resSeq );
		pdb_chains.push_back( rinfo.chainID );
		insertion_codes.push_back( rinfo.iCode );
	}

	// PDBInfo setup
	core::pose::PDBInfoOP pdb_info( new core::pose::PDBInfo( pose.total_residue() ) );

	// set pdb-wide information
	pdb_info->name( fd.filename );
	if(fd.modeltag=="") {
		pdb_info->modeltag( fd.filename );
	} else {
		pdb_info->modeltag( fd.modeltag );
	}

	if( options.preserve_header() == true ) {
		pdb_info->remarks( *fd.remarks );
		pdb_info->header_information( fd.header_information() );
	}

	// set residue level pdb information
	pdb_info->set_numbering( pdb_numbering );
	pdb_info->set_chains( pdb_chains );
	pdb_info->set_icodes( insertion_codes );
	if ( options.preserve_crystinfo() )
		pdb_info->set_crystinfo( fd.crystinfo );


	// most DNA structures lack 5' phosphate groups. 5' phosphates must be built to serve as part of the backbone for
	// atom/fold tree purposes. Here they are made virtual so as not to affect physical calculations.
	for ( core::uint seqpos(1), nres( pose.total_residue() ); seqpos <= nres; ++seqpos ) {
		Residue const & rsd( pose.residue( seqpos ) );
		if ( ! rsd.type().is_DNA() ) continue;
		for ( core::uint atomi(1), natoms( rsd.natoms() ); atomi <= natoms; ++atomi ) {
			id::AtomID const id( atomi, seqpos );
			if ( missing[ id ] && rsd.atom_name(atomi) == " P  " ) {
				tr << "Virtualizing missing phosphate that was built in at seqpos " << seqpos << std::endl;
				core::pose::add_variant_type_to_pose_residue( pose, chemical::VIRTUAL_DNA_PHOSPHATE, seqpos );
				break;
			}
		}
	}


	// Look for and create any remaining non-mainchain (Edge::CHEMICAL) bonds based on a specified radius from any
	// unsatisfied residue connections.  This is used for such things as branched polymers, ubiquitination, or covalent
	// intermediates.  Note: The fold tree will remain with a jump between each such bond until import_pose::
	// set_reasonable_fold_tree() is called later, which actually adds the CHEMICAL edges to fold tree; this method
	// simply makes the bonds.
	pose.conformation().detect_bonds();

	//mjo TODO: this can try to access pose->pdb_info() which is not yet
	//initialized. Moving it after the pose->pdb_info has been
	//initialized causes integration test changes
	core::pose::initialize_disulfide_bonds(pose);

	//kdrew: if detect_oops flag is set, initialize oops
	// This option should probably be moved to FileDataOptions. ~Labonte
	if ( basic::options::option[ basic::options::OptionKeys::in::detect_oops ].user() )
	{
		core::pose::ncbb::initialize_oops(pose);
	}

	if(pose.n_residue()>1){// 1 residue fragments for ligand design.
		pose.conformation().detect_pseudobonds();
	}

	// add contraints based on LINK records if desired
	if ( basic::options::option[ basic::options::OptionKeys::in::constraints_from_link_records ].value() ) {
		core::pose::get_constraints_from_link_records( pose, fd );
	}

	// ensure enough space for atom level pdb information
	pdb_info->resize_atom_records( pose );

	// add unrecognized atoms to PDBInfo
	for( Size i = 1; i <= UA_res_nums.size(); ++i ) {
		pdb_info->add_unrecognized_atom( UA_res_nums[i], UA_res_names[i], UA_atom_names[i], UA_coords[i], UA_temps[i] );
	}

	// add temps to PDBInfo
	for( core::Size ir = 1; ir <= pose.total_residue(); ir++ ) {
		// fill in b-factor from pdb file
		ResidueTemps & res_temps( rinfos[pose_to_rinfo[ir]].temps );
		NameBimap const & namemap( rinfo_name_map[pose_to_rinfo[ir]] );
		for( ResidueTemps::const_iterator iter=res_temps.begin(); iter != res_temps.end(); ++iter ) {
			//namemap should only include atoms which have a presence in both rinfo and pose
			if( namemap.left.count(iter->first) ) {
				// printf("setting temp: res %d atom %s temp %f\n",ir,iter->first.c_str(),iter->second);
				std::string const & pose_atom_name( namemap.left.find(iter->first)->second );
				if( pose.residue(ir).type().has( pose_atom_name ) ) { // There are issues with terminus patching which means atoms can sometimes disappear
					core::Size ia = pose.residue(ir).type().atom_index( pose_atom_name );
					pdb_info->temperature( ir, ia, iter->second );
				}
			} else {
				if( (iter->first)[0] == 'H' || ((iter->first)[0] == ' ' && (iter->first)[1] == 'H') ) {
					;// don't warn if H
				} else {
					tr << "[ WARNING ] can't find atom for res " << ir << " atom " << iter->first << " (trying to set temp)" << std::endl;
				}
			}
		}
	}

	// mark PDBInfo as ok and store in Pose
	pdb_info->obsolete( false );
	pose.pdb_info( pdb_info );

	//fpd fix bfactors of missing atoms using neighbors
	//fpd set hydrogen Bfactors as 1.2x attached atom
	if (options.preserve_crystinfo()) {
		core::scoring::cryst::fix_bfactorsMissing( pose );
		core::scoring::cryst::fix_bfactorsH( pose );
	}
	if ( basic::options::option[ basic::options::OptionKeys::out::file::pdb_comments]() ) {
		utility::io::izstream data(fd.filename);

		std::string line;
		while( getline( data, line ) ) {
			if( line != "##Begin comments##")
			continue;
			getline( data, line );
			while (line != "##End comments##") {
				//TR<<"Testing read comments! :"<<line<<std::endl;
				std::string const key;
				std::string const value;
				utility::vector1<std::string> comment_line(utility::string_split(line,' '));
				if (comment_line.size()<2) {
					getline( data, line );
					continue;
				}
				core::pose::add_comment(pose,comment_line[1],comment_line[2]);
				getline( data, line );
			}
		}
	}
}


// "super-simple" (C) by Phil
/// @brief Try to Build pose object from pdb 'as-is'. PDB file must be _really_ clean.
void build_pose_as_is_TEST(
	io::pdb::FileData & fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options
)
{
	id::AtomID_Mask missing( false );

	build_pose_as_is1_TEST( fd, pose, residue_set, missing, options ); // FOR TESTING.
	build_pose_as_is2( fd, pose, residue_set, missing, options ); // shared with file_data.cc
}


/// @brief Build Rosetta 3 Pose object from FileData.
void build_pose_TEST(
	io::pdb::FileData & fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options
)
{
	tr.Debug << "build_pose..." << std::endl;
	build_pose_as_is_TEST( fd, pose, residue_set, options);
	tr.Debug << "build_pose... Ok." << std::endl;
}


void
pose_from_pdb_TEST(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filenames_string
)
{

	ImportPoseOptions const options;
	bool read_fold_tree( false );

	utility::vector1<std::string> filenames = utility::split(filenames_string);

	std::string res;

	BOOST_FOREACH(std::string filename, filenames){
		utility::io::izstream file( filename );
		if (!file) {
			tr.Error << "PDB File:" << filename << " not found!" << std::endl;
			utility_exit_with_message( "Cannot open PDB file \"" + filename + "\"" );
		} else {
			tr.Debug << "read file: " << filename << std::endl;
		}
		utility::slurp( file, res );
	}

	//fpd If the conformation is not of type core::Conformation, reset it
	conformation::ConformationOP conformation_op( new conformation::Conformation() );
	if ( !pose.conformation().same_type_as_me( *conformation_op, true ) ) {
		pose.set_new_conformation( conformation_op );
	}

	io::pdb::FileData fd = PDB_DReader::createFileData(res, options);
	if ( fd.filename == "" ) {
		fd.filename = utility::join(filenames, "_");
	}
	build_pose_TEST(fd, pose, residue_set, options);

	// set secondary structure for centroid PDBs
	if ( residue_set.name() == core::chemical::CENTROID ) {
		core::pose::set_ss_from_phipsi( pose );
	}

	// check for foldtree info
	read_additional_pdb_data( res, pose, options, read_fold_tree );
}




////////////////////////////////////////////////////////////
void
pose_from_pdb_test() {

	using namespace core::chemical;
	using namespace core::pose;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->
		residue_type_set ( FA_STANDARD );

	Pose pose;
	import_pose::pose_from_pdb ( pose, *rsd_set, option[in::file::s]()[ 1 ] );
	tr << pose.annotated_sequence() << std::endl;

	Pose pose_TEST;
	pose_from_pdb_TEST ( pose_TEST, *rsd_set, option[in::file::s]()[ 1 ] );
	tr << pose_TEST.annotated_sequence() << std::endl;

	runtime_assert( pose.annotated_sequence() == pose_TEST.annotated_sequence() );

}



////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// oh ugly -- not up there to prevent namespace chemical:: changes that would confuse reintegration into core functions.
#include <basic/options/keys/chemical.OptionKeys.gen.hh>

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	try {

	NEW_OPT( sequence_test, "check residue_types_from_sequence()", false );
	NEW_OPT( original_rsd_type_set_map,  "compare to mode with figure out residue_types using original name3_map() code.", false );
	NEW_OPT( original_rsd_type_set_map_only,  "figure out residue_types using original name3_map() code only (for performance comparisons).", false );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);

	option[ OptionKeys::chemical::patch_selectors ].push_back( "PEPTIDE_CAP" ); // N_acetylated.txt and C_methylamidated.txt

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	if ( option[ sequence_test ]() ) {
		residue_types_from_sequence_test();
	} else {
		pose_from_pdb_test();
	}

	exit( 0 );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
