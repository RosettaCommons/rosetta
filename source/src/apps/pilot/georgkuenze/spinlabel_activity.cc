// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    spinlabel_activity.cc
/// @brief   This pilot app inserts a spin-label residue at a specified position
///          and tests all rotamers for this spin-label. Those romaters which remain
///          after the bump filter (with or without allowing packing of neighboring
///          residues) are written out as SDF files. Finally, the native residue is
///          replaced by the spin-label residue and the mutated pose is written out
///          as a PDB file.
/// @details last Modified: 08/18/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Project headers
#include <devel/init.hh>

#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/rotamers/SingleLigandRotamerLibrary.hh>
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResidueLevelTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/annealer/AnnealerFactory.hh>
#include <core/pack/annealer/SimAnnealerBase.hh>
#include <core/pack/interaction_graph/InteractionGraphFactory.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <core/pack/palette/CustomBaseTypePackerPalette.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/sdf/mol_writer.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>
#include <utility/graph/Graph.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// C++ headers
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

OPT_KEY ( Integer, label_position )
OPT_KEY ( Boolean, repacking )
OPT_KEY ( String, spinlabel_code )

int
main( int argc, char** argv )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::chemical;
		using namespace core::pack;
		using namespace core::pack::rotamers;
		using namespace core::pack::task;
		using namespace core::pack::palette;
		using namespace core::pack::annealer;
		using namespace core::pack::interaction_graph;
		using namespace core::pack::rotamer_set;
		using namespace core::conformation;
		using namespace core::pose;
		using namespace core::scoring;

		OPT(in::file::s);
		NEW_OPT( label_position, "Position of MTSL spinlabel", 1);
		NEW_OPT( repacking, "Repack side chains around spinlabel", true);
		NEW_OPT( spinlabel_code, "Three-letter code of spinlabel", "R1A" );

		// Initialize core
		// Get input
		devel::init( argc, argv );
		utility::vector1<std::string> filenames = option[ in::file::s ]();
		core::Size seqid = option[ label_position ];
		bool repack(option[ repacking ]);
		std::string sl_code = option[ spinlabel_code ];

		// Setup pose
		PoseOP myPose;
		if ( !(filenames.size() > 0) ) {
			utility_exit_with_message( "No input PDB file provided." );
		} else {
			myPose = core::import_pose::pose_from_file(filenames[1]);
		}

		// Setup scorefunction
		ScoreFunctionOP sfxn = core::scoring::get_score_function();

		////////////////////////////////////////////////////////////////////////////
		//////////                                                        //////////
		//////////      Part 1 - MTSL rotamer library and bump check      //////////
		//////////                                                        //////////
		////////////////////////////////////////////////////////////////////////////

		// Setup packer task
		PackerTaskOP pack1 = TaskFactory::create_packer_task( *myPose );
		utility::vector1<bool> residues_to_pack(myPose->total_residue(), false);
		if ( repack ) {
			std::fill(residues_to_pack.begin(), residues_to_pack.end(), true);
		}
		residues_to_pack[seqid] = true;
		pack1->restrict_to_residues(residues_to_pack);
		pack1->restrict_to_repacking();
		utility::graph::GraphOP packer_neighbor_graph = create_packer_graph( *myPose, *sfxn, pack1 );

		// Create residue type set
		ResidueTypeSetCOP residue_set( ChemicalManager::get_instance()->residue_type_set( "fa_standard" ) );

		// Get the MTSL spinlabel residue type
		ResidueTypeCOP sl_residue_type = residue_set->get_representative_type_name3(sl_code);
		std::cout << "Getting Spinlabel residue type: " << sl_residue_type->name3() << std::endl;

		// Get the MTSL spinlabel rotamer library
		SingleResidueRotamerLibraryFactory const & rotlib_fac( *SingleResidueRotamerLibraryFactory::get_instance() );
		SingleResidueRotamerLibraryCOP sl_rotlib(rotlib_fac.get( *sl_residue_type ));

		if ( sl_rotlib ) {
			SingleLigandRotamerLibraryCOP sl_lig_rotlib(utility::pointer::dynamic_pointer_cast< SingleLigandRotamerLibrary const > ( sl_rotlib ));
			if ( sl_lig_rotlib == 0 ) {
				utility_exit_with_message( "Failed to retrieve single ligand rotamer library for " + sl_residue_type->name3() );
			}

			utility::vector1<ResidueOP> rot_vector;

			// Fill the rotamer vector (we could also use build_base_rotamers()) and perform the bump_filter()
			//sl_lig_rotlib->build_base_rotamers( *sl_residue_type, rot_vector );
			//rot_vector.clear();
			bool buried = false;
			utility::vector1< utility::vector1< core::Real > > extra_chi_steps;
			sl_lig_rotlib->fill_rotamer_vector( *myPose, *sfxn, *pack1, packer_neighbor_graph, sl_residue_type, myPose->residue(seqid), extra_chi_steps, buried, rot_vector);
			sl_lig_rotlib->bump_filter(rot_vector, seqid, *sfxn, *myPose, *pack1, packer_neighbor_graph);
			std::cout << "Get MTSL rotamers and perform bump check." << std::endl;

			// For visualization output the MTSL rotamers as SDF file
			std::cout << "Write coordinates for accepted MTSL rotamers to file." << std::endl;
			core::chemical::sdf::MolWriter writer;
			for ( core::Size i = 1; i <= rot_vector.size(); ++i ) {
				std::ostringstream convert;
				convert << std::setw(4) << std::setfill('0') << i;
				// output name
				std::string output_filename = utility::strip(filenames[1], ".pdb");
				if ( repack ) { output_filename += "_" + sl_code + "_rotamer_with_repacking_" + convert.str() + ".mol"; }
				else { output_filename += "_" + sl_code + "_rotamer_no_repacking_" + convert.str() + ".mol"; }
				writer.output_residue(output_filename, *rot_vector[i]);
			}

		} else {
			utility_exit_with_message( "Failed to retrieve single residue rotamer library for " + sl_residue_type->name3() );
		}

		// Output original pdb too show that the native residue was not mutated.
		std::cout << "Write original protein PDB coordinates to file." << std::endl;
		std::string output = utility::strip(filenames[1], ".pdb");
		output += "_.pdb";
		myPose->dump_pdb( output );

		///////////////////////////////////////////////////////////////////////////////
		//////////                                                           //////////
		//////////      Part 2 - one MTSL rotamer at a time and packing      //////////
		//////////                                                           //////////
		///////////////////////////////////////////////////////////////////////////////

		// Create PackerPalette (VKM, Jan 2019):
		CustomBaseTypePackerPaletteOP palette( utility::pointer::make_shared< CustomBaseTypePackerPalette >() );
		palette->add_type( sl_code );

		// Setup packer task
		PackerTaskOP pack2 = TaskFactory::create_packer_task( *myPose, palette );
		if ( repack ) {
			std::fill(residues_to_pack.begin(), residues_to_pack.end(), true);
		} else {
			std::fill(residues_to_pack.begin(), residues_to_pack.end(), false);
		}
		residues_to_pack[seqid] = true;
		pack2->restrict_to_residues(residues_to_pack);
		for ( core::Size ii = 1; ii <= myPose->total_residue(); ++ii ) {
			if ( ii == seqid ) continue;
			pack2->nonconst_residue_task(ii).restrict_to_repacking();
		}

		// Turn on at the spinlabel position
		pack2->nonconst_residue_task(seqid).restrict_restypes( utility::vector1< std::string >( seqid ) );

		// Setup packing
		RotamerSetsOP rotsets = RotamerSetsOP(new rotamer_set::RotamerSets());
		pack_scorefxn_pose_handshake( *myPose, *sfxn);
		myPose->update_residue_neighbors();
		sfxn->setup_for_packing( *myPose, pack2->repacking_residues(), pack2->designing_residues() );
		packer_neighbor_graph = create_packer_graph( *myPose, *sfxn, pack2 );

		// Generate rotamer set
		rotsets->set_task( pack2 );
		rotsets->initialize_pose_for_rotsets_creation( *myPose );
		rotsets->build_rotamers( *myPose, *sfxn, packer_neighbor_graph );

		// Pull out the specific rotamer set for the MTSL spinlabel
		RotamerSetOP rotset_sl = rotsets->rotamer_set_for_residue( seqid );
		std::cout << "Generated " << rotset_sl->num_rotamers() << " rotamers of residue type " << rotset_sl->rotamer(1)->name3()
			<< " at position " << seqid << " with native residue " << myPose->residue(seqid).name3() << "." << std::endl;

		//std::cout << "Rotamer offsets before changing the spinlabel rotamer set: "
		//  << "Residue " << (seqid - 1) << ": " << rotsets->nrotamer_offset_for_moltenres(rotsets->resid_2_moltenres( seqid - 1 ) ) << "; "
		//  << "Residue " << (seqid    ) << ": " << rotsets->nrotamer_offset_for_moltenres(rotsets->resid_2_moltenres( seqid     ) ) << "; "
		//  << "Residue " << (seqid + 1) << ": " << rotsets->nrotamer_offset_for_moltenres(rotsets->resid_2_moltenres( seqid + 1 ) ) << std::endl;

		// Store the original spinlabel rotamers to iterate over them later
		core::Size num_rotamers(rotset_sl->num_rotamers());
		utility::vector1< ResidueOP > sl_rotamers(rotset_sl->num_rotamers());
		utility::vector1< ResidueOP >::const_iterator sl_rotamers_iter;
		core::Size ii(1);
		for ( sl_rotamers_iter = rotset_sl->begin(); sl_rotamers_iter != rotset_sl->end(); ++sl_rotamers_iter ) {
			sl_rotamers[ii] = *sl_rotamers_iter;
			++ii;
		}

		// Now fill the spinlabel rotamerset vector with only one type of rotamer
		// First, we add the one rotamer type to the end of the vector times the total number (N) of rotamers
		// Second, we delete the N number of rotamers at the beginning of the vector
		utility::vector1< bool > sl_rotamers_to_delete_mask(2*num_rotamers, false);
		std::fill(sl_rotamers_to_delete_mask.begin(), sl_rotamers_to_delete_mask.begin() + num_rotamers, true);
		for ( core::Size j = 1; j <= num_rotamers; ++j ) {
			for ( core::Size k = 1; k <= num_rotamers; ++k ) {
				rotset_sl->add_rotamer( *sl_rotamers[j] );
			}
			rotset_sl->drop_rotamers(sl_rotamers_to_delete_mask);

			// Did we change the offset of rotamer ids in the RotamerSets?
			// if (j == 1) {
			//  std::cout << "Rotamer offsets after  changing the spinlabel rotamer set: "
			//    << "Residue " << (seqid - 1) << ": " << rotsets->nrotamer_offset_for_moltenres(rotsets->resid_2_moltenres( seqid - 1 ) ) << "; "
			//    << "Residue " << (seqid    ) << ": " << rotsets->nrotamer_offset_for_moltenres(rotsets->resid_2_moltenres( seqid     ) ) << "; "
			//    << "Residue " << (seqid + 1) << ": " << rotsets->nrotamer_offset_for_moltenres(rotsets->resid_2_moltenres( seqid + 1 ) ) << std::endl;
			// }

			// Do packing
			rotsets->prepare_sets_for_packing( *myPose, *sfxn );
			AnnealableGraphBaseOP ig = InteractionGraphFactory::create_and_initialize_annealing_graph( *pack2, *rotsets, *myPose, *sfxn, packer_neighbor_graph );
			ObjexxFCL::FArray1D< int > bestrotamer_at_seqpos( myPose->total_residue() );
			core::PackerEnergy bestenergy( 0.0 );
			utility::vector0< int > rot_to_pack;
			// Packing with SA but no replacement of residue
			pack_rotamers_run( *myPose, pack2, rotsets, ig, rot_to_pack, bestrotamer_at_seqpos, bestenergy );

			std::cout << "Best spinlabel rotamer at resid " << seqid << ": ID " << rotsets->rotid_on_moltenresidue(bestrotamer_at_seqpos(seqid)) << std::endl;
			std::cout << "PackerEnergy " << bestenergy << std::endl;
		}

		///////////////////////////////////////////////////////////////////////////////
		//////////                                                           //////////
		//////////      Part 3 - mutation to MTSL and packing of others      //////////
		//////////                                                           //////////
		///////////////////////////////////////////////////////////////////////////////

		// Setup packer task
		PackerTaskOP pack3 = TaskFactory::create_packer_task( *myPose, palette ); //VKM, Jan 2019: Note that for now, we can reuse the same PackerPalette.  If this step used a *different* spin label, we'd have to change this out for a new one, here.
		if ( repack ) {
			std::fill(residues_to_pack.begin(), residues_to_pack.end(), true);
		} else {
			std::fill(residues_to_pack.begin(), residues_to_pack.end(), false);
		}
		residues_to_pack[seqid] = true;
		pack3->restrict_to_residues(residues_to_pack);
		for ( core::Size ii = 1; ii <= myPose->total_residue(); ++ii ) {
			if ( ii == seqid ) continue;
			pack3->nonconst_residue_task(ii).restrict_to_repacking();
		}
		core::chemical::ResidueTypeSetCOP residue_set_at_spinlabel_site = pack3->nonconst_residue_task(seqid).get_original_residue_set();
		std::string const & spinlabel_to_include(sl_code);

		// Does this ResidueTypeSet have ResidueTypes with the given interchangeability group?
		// Then add the Spinlabel residue type
		runtime_assert_string_msg ( residue_set_at_spinlabel_site->has_interchangeability_group(spinlabel_to_include),
			"Unable to add non-canonical amino acid(s) with interchangeability group " + spinlabel_to_include + " because there are no ResidueTypes with that interchangeability group in the ResidueTypeSet for residue " + std::to_string( seqid ) + "." );

		// Turn off canonical amino acids at the spinlabeled position
		pack3->nonconst_residue_task(seqid).restrict_restypes( utility::vector1< std::string >( { spinlabel_to_include } ) );

		// Perform packing of protein and mutation of spinlabeled residue to spinlabel residue
		pack_rotamers( *myPose, *sfxn, pack3);

		std::cout << "Write PDB file after spinlabel mutation and repacking of protein." << std::endl;
		std::string output_second = utility::strip(filenames[1], ".pdb");
		output_second += "_" + sl_code + "_labeled.pdb";
		myPose->dump_pdb( output_second );


	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}

