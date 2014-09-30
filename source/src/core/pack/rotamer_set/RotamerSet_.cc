// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/RotamerSet_.cc
/// @brief  amino acid rotamer set class implementation
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/pack/rotamer_set/RotamerSet_.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pack/rotamer_set/rotamer_building_functions.hh>
#include <core/pack/rotamer_set/rna_rotamer_building_functions.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/ChiSet.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/SingleResidueRotamerLibrary.hh>
#include <core/pack/interaction_graph/SurfacePotential.hh>

// Project Headers
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.functions.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/AbstractRotamerTrie.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
#include <core/scoring/methods/Methods.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <iostream>

namespace core {
namespace pack {
namespace rotamer_set {

static thread_local basic::Tracer tt( "core.pack.rotamer_set.RotamerSet_", basic::t_info );

RotamerSet_::RotamerSet_()
:
	n_residue_types_( 0 ),
	n_residue_groups_( 0 ),
	cached_tries_( scoring::methods::n_energy_methods, 0 ),
	id_for_current_rotamer_( 0 ),
	rotamer_offsets_require_update_( false )
{}

RotamerSet_::~RotamerSet_() {}

void
RotamerSet_::build_rotamers(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	task::PackerTask const & the_task,
	graph::GraphCOP packer_neighbor_graph,
	bool use_neighbor_context
)
{
	using namespace chemical;

	for ( task::ResidueLevelTask::ResidueTypeCOPListConstIter
			allowed_iter = the_task.residue_task( resid() ).allowed_residue_types_begin(),
			allowed_end = the_task.residue_task( resid() ).allowed_residue_types_end();
			allowed_iter != allowed_end; ++allowed_iter ) {
		build_rotamers_for_concrete_virt( pose, scorefxn, the_task, *allowed_iter, packer_neighbor_graph, use_neighbor_context );
		//std::cout << "Built rotamers for " << (*allowed_iter)->name() << " seqpos " << resid() << " " << n_residue_types_ << std::endl;
	}

	if ( num_rotamers() == 0 ) {
		tt << "[ WARNING ]  including current in order to get at least 1 rotamer !!!!!! " << resid() << ' ' <<
			pose.residue( resid() ).name() << '\n';
		id_for_current_rotamer_ = num_rotamers();
		prepare_for_new_residue_type( pose.residue_type( resid() ));
		ResidueOP rot = pose.residue( resid() ).create_rotamer();
		push_back_rotamer( rot );
	}

	for ( RotSetOperationListIterator
			rotsetop_iter = the_task.residue_task( resid() ).rotamer_set_operation_begin(),
			rotsetop_end = the_task.residue_task( resid() ).rotamer_set_operation_end();
			rotsetop_iter != rotsetop_end; ++rotsetop_iter ) {
		(*rotsetop_iter)->alter_rotamer_set( pose, scorefxn, the_task, packer_neighbor_graph, *this );
	}

	tt.flush();
	//std::cout << "Built " << num_rotamers() << " rotamers for residue " << resid() << " " << pose.residue(resid()).name() << std::endl;
}

void
RotamerSet_::add_rotamer(
	conformation::Residue const & rotamer
)
{
	//tt << "RotamerSet_::add_rotamer" << '\n';
	//++n_residue_types_; // TOTAL HACK -- possible to add 50 HIS's and declare there to be 50 different AA types.
	//tt << "TOTAL HACK BEING USED in " << __FILE__ << " line " << __LINE__ << ".  For efficiency, refactor!" << '\n';
	//residue_type_rotamers_begin_.push_back( num_rotamers() + 1);
	//n_rotamers_for_restype_.push_back( 1 );
	//rotamers_.push_back( rotamer.clone() );

	prepare_for_new_residue_type( rotamer.type() );
	push_back_rotamer( rotamer.clone() );
}


Size
RotamerSet_::get_n_residue_types() const
{
	//std::cout << "B get_n_residue_types: " << n_residue_types_ << std::endl;
	update_rotamer_offsets();
	//std::cout << "A get_n_residue_types: " << n_residue_types_ << std::endl;
	return n_residue_types_;
}

Size
RotamerSet_::get_n_residue_groups() const
{
	update_rotamer_offsets();
	return n_residue_groups_;
}

Size
RotamerSet_::get_residue_type_begin( Size which_restype ) const
{
	update_rotamer_offsets();
	assert( which_restype <= n_residue_types_ );
	return residue_type_rotamers_begin_[ which_restype ];
}

Size
RotamerSet_::get_residue_group_begin( Size which_resgroup ) const
{
	update_rotamer_offsets();
	assert( which_resgroup <= n_residue_groups_ );
	return residue_group_rotamers_begin_[ which_resgroup ];
}

Size
RotamerSet_::get_n_rotamers_for_residue_type( Size which_restype ) const
{
	update_rotamer_offsets();
	assert( which_restype <= n_residue_types_ );
	return n_rotamers_for_restype_[ which_restype ];
}

Size
RotamerSet_::get_n_rotamers_for_residue_group( Size which_resgroup ) const
{
	update_rotamer_offsets();
	assert( which_resgroup <= n_residue_groups_ );
	return n_rotamers_for_resgroup_[ which_resgroup ];
}

Size
RotamerSet_::get_residue_type_index_for_rotamer( Size which_rotamer ) const
{
	return residue_type_for_rotamers_[ which_rotamer ];
}

Size
RotamerSet_::get_residue_group_index_for_rotamer( Size which_rotamer ) const
{
	return residue_group_for_rotamers_[ which_rotamer ];
}


/// @details This will check for rotamer/background collisions if the PackerTask's
/// bump_check boolean is true -- it uses the bump_check_{sidechain/full} methods
/// defined by the energy methods contained in the scorefxn object
void
RotamerSet_::build_rotamers_for_concrete_virt(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	task::PackerTask const & task,
	chemical::ResidueTypeCOP concrete_residue,
	graph::GraphCOP packer_neighbor_graph,
	bool use_neighbor_context
)
{
	conformation::Residue const & existing_residue( pose.residue( resid() ));
	build_rotamers_for_concrete(
		pose, scorefxn, task, concrete_residue, existing_residue,
		packer_neighbor_graph, use_neighbor_context );
}

void
RotamerSet_::build_rotamers_for_concrete(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const & existing_residue,
		graph::GraphCOP packer_neighbor_graph,
		bool use_neighbor_context)
{
	using namespace conformation;
	using namespace pack::task;

	prepare_for_new_residue_type( *concrete_residue );

	if ( task.residue_task( resid() ).optimize_h() ) {
		build_optimize_H_rotamers( pose, task, concrete_residue, existing_residue, packer_neighbor_graph );
	// The behavior depends on the residue type.  This should be refactored -- at least into several separate methods
	// that this one switches between....
	} else if ( concrete_residue->is_DNA() ) { // DNA /////////////////////////////////////////////////////////////////
		utility::vector1< ResidueOP > new_rotamers;

		build_dna_rotamers( resid(), pose, concrete_residue, task, new_rotamers );

		if ( task.include_current( resid() ) && existing_residue.name() == concrete_residue->name() ) {
			ResidueOP rot = existing_residue.create_rotamer();
			new_rotamers.push_back( rot );
			id_for_current_rotamer_ = new_rotamers.size();
		}

		for ( Size ii=1; ii<= new_rotamers.size(); ++ii ) {
			assert( new_rotamers[ii]->seqpos() == resid() && new_rotamers[ii]->chain() == existing_residue.chain() );
			push_back_rotamer( new_rotamers[ii] );
		}

	} else if ( concrete_residue->is_RNA() ) { // RNA /////////////////////////////////////////////////////////////////
		utility::vector1< ResidueOP > new_rotamers;

		// sample chi, include_current, and proton chi expansion is inside here:
		build_rna_rotamers( resid(), pose, concrete_residue, task, new_rotamers, id_for_current_rotamer_ );

		for ( Size ii=1; ii<= new_rotamers.size(); ++ii ) {
			//assert( new_rotamers[ii]->seqpos() == resid() && new_rotamers[ii]->chain() == existing_residue.chain() );
			push_back_rotamer( new_rotamers[ii] );
		}

	} else if (concrete_residue->is_carbohydrate()) { // Carbohydrates ////////////////////////////////////////////////
		// For now, rotamer bins are stored in the params files themselves.  This code will generate all the rotamers
		// from the torsion angles listed in the params files.  (The params files do not contain PROTON_CHI records and
		// use CHI_ROTAMER records exclusively.)  All of this will likely change in the future.

		tt.Warning << "Residue " << resid() <<
				" is a carbohydrate; rotamers for carbohydrates are not yet fully implemented." << std::endl;

		utility::vector1<ResidueOP> new_rotamers;

		build_rotamers_from_rotamer_bins(existing_residue, new_rotamers);

		// Add current rotamer, if applicable.
		if (task.include_current(resid()) && existing_residue.name() == concrete_residue->name()) {
			ResidueOP rot = existing_residue.create_rotamer();
			new_rotamers.push_back(rot);
			id_for_current_rotamer_ = new_rotamers.size();
		}

		// Push back rotamers.
		Size n_rotamers = new_rotamers.size();
		for (uint i=1; i <= n_rotamers; ++i) {
			push_back_rotamer(new_rotamers[i]);
		}

	} else if ( concrete_residue->name() == "VRT1" ) { // Single-atom virtual residue /////////////////////////////////
		tt << "building VRT1 residue at " << resid() << ' ' << existing_residue.name() << "position\n";
		if ( existing_residue.name() == concrete_residue->name() ) {
			ResidueOP rot = existing_residue.clone();
			push_back_rotamer( rot );
			id_for_current_rotamer_ = num_rotamers();
		} else {
			ResidueOP rot = ResidueFactory::create_residue( *concrete_residue );
			rot->set_xyz( 1, existing_residue.nbr_atom_xyz() );
			rot->seqpos( existing_residue.seqpos() );
			rot->chain ( existing_residue.chain() );
			push_back_rotamer( rot );
		}

	} else if ( concrete_residue->name() == "TP3" ) { // TIP3 water /////////////////////////////////
		build_tp3_water_rotamers( pose, task, concrete_residue, existing_residue, packer_neighbor_graph );
	} else { // All other residues ///////////////////////////////////////////////////////////////////////

		utility::vector1< utility::vector1< Real > > extra_chi_steps( concrete_residue->nchi() );

		Size nneighbs(999);
		if ( use_neighbor_context ) {
			nneighbs = pose.energies().tenA_neighbor_graph().get_node( resid() )->num_neighbors_counting_self();
		}
		bool buried = ( nneighbs >= task.residue_task(resid()).extrachi_cutoff() || !use_neighbor_context );

		for ( Size ii = 1; ii <= concrete_residue->nchi(); ++ii ) {
			set_extra_samples(task, nneighbs, ii, concrete_residue, extra_chi_steps[ ii ] );
		}

		bump_selector_.reset();
		bump_selector_.set_max_rot_bumpenergy( task.max_rotbump_energy() );

		utility::vector1< ResidueOP > suggested_rotamers;
		dunbrack::SingleResidueRotamerLibraryCOP rotlib = dunbrack::RotamerLibrary::get_instance().get_rsd_library( *concrete_residue ); //For D-amino acids, returns the rotamer library for the corresponding L-amino acid
		if (rotlib) {
			/// DOUG DOUG DOUG DEBUG OUTPUT
			//std::cout << "EXTRA_CHI_STEPS::build_rotamers_for_concrete\t" << extra_chi_steps.size() << std::endl;
			//for ( Size i(1); i <= extra_chi_steps.size(); ++i ) std::cout << extra_chi_steps[i].size() << std::endl;

			//for ( Size i(1); i <= extra_chi_steps.size(); ++i ) {
			//	for ( Size j(1); j <= extra_chi_steps[i].size(); ++j ) {
			//		std::cout << i << "/" << j << ":\t" << extra_chi_steps[i][j] << "\t" << std::flush;
			//	}
			//	std::cout << std::endl;
			//}
			//std::cout << std::endl;

			rotlib->fill_rotamer_vector( pose, scorefxn, task, packer_neighbor_graph, concrete_residue, existing_residue, extra_chi_steps, buried, suggested_rotamers);
			if(core::chemical::is_canonical_D_aa( existing_residue.aa() ) && suggested_rotamers.size() > 0 ) { //If this is a D-amino acid, flip all the chi values in the suggested_rotamers vector
				for(core::Size i=1; i<=suggested_rotamers.size(); i++) {
					if(suggested_rotamers[i]->nchi() > 0) {
						for(core::Size j=1; j<=suggested_rotamers[i]->nchi(); j++) {
							suggested_rotamers[i]->set_chi(j, -1.0*suggested_rotamers[i]->chi(j));
						}
					}
				}
			}
		} else {
			// Add proton chi rotamers even if there's no rotamer library... this will disappear when
			// ligands get their own rotamer libraries.
			if ( concrete_residue->n_proton_chi() != 0 ) {
				utility::vector1< pack::dunbrack::ChiSetOP > proton_chi_chisets;
				proton_chi_chisets.push_back( dunbrack::ChiSetOP( new pack::dunbrack::ChiSet( concrete_residue->nchi() ) ) );
				for ( Size ii = 1; ii <= concrete_residue->n_proton_chi(); ++ii ) {
					pack::dunbrack::expand_proton_chi(
						task.residue_task( resid() ).extrachi_sample_level(
							( (Size) nneighbs >=  task.residue_task( resid() ).extrachi_cutoff() ),
							concrete_residue->proton_chi_2_chi( ii ),
							*concrete_residue ),
						concrete_residue,
						ii, proton_chi_chisets);
				}


				suggested_rotamers.reserve( proton_chi_chisets.size() );
				for ( Size ii = 1; ii <= proton_chi_chisets.size(); ++ii ) {
					suggested_rotamers.push_back( existing_residue.clone() );
					for ( Size jj = 1; jj <= concrete_residue->n_proton_chi(); ++jj ) {
						Size jj_protchi = concrete_residue->proton_chi_2_chi( jj );
						suggested_rotamers[ ii ]->set_chi(jj_protchi, proton_chi_chisets[ ii ]->chi[ jj_protchi ] );
					}
				}
			}
			tt.Debug << "Building rotamers for " << concrete_residue->name() << ": Total number suggested: " <<
					suggested_rotamers.size() << std::endl;
		}

		tt.Debug <<
			"Existing residue: " << existing_residue.name() <<
			" Concrete residue: " << concrete_residue->name() <<
			" Total number of rotamers suggested: " << suggested_rotamers.size() << std::endl;

		// Look through suggested rotamers and throw out those with obvious "bumps".
		for ( Size ii = 1; ii <= suggested_rotamers.size(); ++ii ) {
			ResidueOP rot = suggested_rotamers[ ii ];
			tt.Debug << " Suggested rotamer " << ii << ": ";
			for ( Size jj = 1; jj <= rot->nchi(); ++jj ) {
				tt.Debug << rot->chi()[jj] << ' ';
			}

			if ( task.bump_check() ) {
				core::PackerEnergy bumpenergy = bump_check( rot, scorefxn, pose, task, packer_neighbor_graph );
				tt.Debug << " Bump energy: " << bumpenergy;
				BumpSelectorDecision decision =  bump_selector_.iterate_bump_selector( bumpenergy );
				switch ( decision ) {
					case KEEP_ROTAMER :
						tt.Debug << " ... added" << std::endl;
						push_back_rotamer( rot );
						break;
					case DELETE_PREVIOUS_ROTAMER :
						tt.Debug << " ... replace previous" << std::endl;
						assert ( num_rotamers() > 0 );
						rotamers_[ num_rotamers() ] = rot;
						break;
					case DELETE_ROTAMER : // Do nothing.
						tt.Debug << " ... deleted" << std::endl;
						break;
				}
			} else {
				tt.Debug << " ... added unchecked" << std::endl;
				push_back_rotamer( rot );
			}
		}

		// Do we include the current residue?
		if ( task.include_current( resid() ) && existing_residue.name() == concrete_residue->name() ) {
			 tt.Debug << " Including current rotamer: " <<
					existing_residue.name() << ' ' << existing_residue.seqpos() << std::endl;
			ResidueOP rot = existing_residue.create_rotamer();
			push_back_rotamer( rot );
			id_for_current_rotamer_ = num_rotamers();

		// ...Or an "emergency rotamer"?
		} else if ( suggested_rotamers.size() == 0 && concrete_residue->nchi() == 0 ) {
			tt.Debug << " Using emergency rotamer" << std::endl;
			ResidueOP rot = ResidueFactory::create_residue( *concrete_residue, existing_residue, pose.conformation() );
			push_back_rotamer( rot );
		}

		// not ready for design yet.
		if ( task.residue_task( resid() ).include_virtual_side_chain() ) {
			if ( existing_residue.nchi() > 0 &&
					 existing_residue.aa() != chemical::aa_pro ) {
				dunbrack::RotamerLibraryScratchSpace scratch;
				Size n_min( 0 );
				Real fa_dun_min( 0.0 );
				for ( Size n = 1; n <= suggested_rotamers.size(); n++ ){
					Real const fa_dun = rotlib->rotamer_energy( *suggested_rotamers[ n ], scratch );
					if ( n_min == 0 || fa_dun < fa_dun_min ) {
						n_min = n; fa_dun_min = fa_dun;
					}
					// hack for now. PackerTask must be instantiated with a pose without virtual sidechains.
					runtime_assert( !suggested_rotamers[n]->has_variant_type( chemical::VIRTUAL_SIDE_CHAIN ) );
				}
				ResidueOP rot = core::pose::add_variant_type_to_residue( *suggested_rotamers[ n_min ], chemical::VIRTUAL_SIDE_CHAIN, pose);
				push_back_rotamer( rot );
			}
		}


	}
}

/// @details Creates a sets of rotamers for an "optimize H" repacking:
/// optimize H repositions hydroxyl hydrogens and resolve protonation
/// ambiguity and if the task's flip_HNQ flag is set, it will flip
/// amide groups ( ASN (N) and GLN (Q) ) and the ring orientation
/// in histidine (H).
void
RotamerSet_::build_optimize_H_rotamers(
	pose::Pose const & pose,
	task::PackerTask const & task,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const & existing_residue,
	graph::GraphCOP packer_neighbor_graph
)
{
	using namespace chemical;
	using namespace conformation;

	// Three cases:
	// 1) if we're flipping HNQ's and this is an HNQ, then push back two rotamers, an original
	// and a flipped rotamer.  If it's His and concrete.type != existing.type, then copy
	// the heavy coords and place the H's in their ideal conf.
	// 2) if the concrete residue is not the same as the existing residue, then copy
	// over the heavyatom coordinates and idealize the hydrogen coordinates. HIS case w/o flip.
	// 3) otherwise, create proton-chi-varying rotamers only, cloning the existing residue


	if ( task.residue_task( resid() ).flip_HNQ() && (
		  concrete_residue->aa() == aa_his ||
			concrete_residue->aa() == aa_asn ||
			concrete_residue->aa() == aa_gln ) ) {

		/// This is not the Richardson style HNQ flip.  Instead of preserving the bond geometry
		/// and rotating chi2/chi3 by 180, they swap the coordinates of the heavy atoms.
		/// e.g. ASN ND2 swaped with ASN OD1.
		/// rosetta++ simply rotates the bond angles.

		if ( concrete_residue->name() != pose.residue( resid() ).name() ) {
			/// HisE --> HisD or HisD --> HisE
			ResidueOP example_rotamer = ResidueFactory::create_residue( *concrete_residue );
			example_rotamer->seqpos( existing_residue.seqpos() );
			for ( Size ii = 1; ii <= existing_residue.nheavyatoms(); ++ii ) {
				if ( example_rotamer->has( existing_residue.atom_name( ii ) ) ) {
					example_rotamer->set_xyz(
						example_rotamer->atom_index( existing_residue.atom_name( ii ) ),
						existing_residue.xyz( ii ) );
				}
			}
			example_rotamer->chain( pose.residue( resid() ).chain() );
			example_rotamer->mainchain_torsions() = pose.residue( resid() ).mainchain_torsions();
			example_rotamer->copy_residue_connections( pose.residue( resid() ) );
			idealize_hydrogens( *example_rotamer, pose.conformation() );
			push_back_rotamer( example_rotamer );

			ResidueOP flipped_rotamer = example_rotamer->clone();
			Real flipped_chi2 = flipped_rotamer->chi(2) + 180;
			flipped_rotamer->set_chi( 2, flipped_chi2 );
			push_back_rotamer( flipped_rotamer );

		} else {
			push_back_rotamer( existing_residue.clone() );

			ResidueOP flipped_rotamer = existing_residue.clone();
			Size chi_to_flip( 0 );
			switch ( concrete_residue->aa() ) {
				case aa_his :
				case aa_asn :
					chi_to_flip = 2;
				break;
				case aa_gln :
					chi_to_flip = 3;
				break;
				default:
					utility_exit_with_message("Illegal case statement option.");
				break;
			}

			conformation::set_chi_according_to_coordinates( *flipped_rotamer );

			Real flipped_chi2 = flipped_rotamer->chi( chi_to_flip ) + 180;
			flipped_rotamer->set_chi( chi_to_flip, flipped_chi2 );
			push_back_rotamer( flipped_rotamer );
		}

	} else if ( concrete_residue->is_NA() ) {
			/// merely clone the input residue
			push_back_rotamer( existing_residue.clone() );
	} else if ( concrete_residue->name() != pose.residue( resid() ).name() ) {
		// in particular, HIS can be protonated on either ND1 or NE2
		// Note: there is an assumption here that there are no proton chi
		// in residues with alternate H placements.  This assumption is unneccessary
		// and can be changed by modifying the code below to expand proton chi.

		ResidueOP example_rotamer = ResidueFactory::create_residue( *concrete_residue );
		example_rotamer->seqpos( existing_residue.seqpos() );
		for ( Size ii = 1; ii <= existing_residue.nheavyatoms(); ++ii ) {
			if ( example_rotamer->has( existing_residue.atom_name( ii ) ) ) {
				example_rotamer->set_xyz(
					example_rotamer->atom_index( existing_residue.atom_name( ii ) ),
					existing_residue.xyz( ii ) );
			}
		}
		example_rotamer->chain( pose.residue( resid() ).chain() );
		example_rotamer->mainchain_torsions() = pose.residue( resid() ).mainchain_torsions();
		example_rotamer->copy_residue_connections( pose.residue( resid() ) );
		idealize_hydrogens( *example_rotamer, pose.conformation() );
		conformation::set_chi_according_to_coordinates( *example_rotamer );

		push_back_rotamer( example_rotamer );
	} else if ( concrete_residue->name() == "TP3" ) { // TIP3 water /////////////////////////////////
		build_tp3_water_rotamers( pose, task, concrete_residue, existing_residue, packer_neighbor_graph );
	} else {
		/// Rotatable proton chi
		utility::vector1< ResidueOP > suggested_rotamers;

		// REFACTOR!!!
		if ( concrete_residue->n_proton_chi() != 0 ) {
			utility::vector1< pack::dunbrack::ChiSetOP > proton_chi_chisets;
			proton_chi_chisets.push_back(
				dunbrack::ChiSetOP( new pack::dunbrack::ChiSet( concrete_residue->nchi() ) ) );
			for ( Size ii = 1; ii <= concrete_residue->n_proton_chi(); ++ii ) {
				pack::dunbrack::expand_proton_chi(
					task.residue_task( resid() ).extrachi_sample_level(
						true, // ignore buriedness when adding extra proton chi rotamers
						concrete_residue->proton_chi_2_chi( ii ),
						*concrete_residue ),
					concrete_residue,
					ii, proton_chi_chisets);
			}
			suggested_rotamers.reserve( proton_chi_chisets.size() );
			for ( Size ii = 1; ii <= proton_chi_chisets.size(); ++ii ) {
				suggested_rotamers.push_back( existing_residue.clone() );
				for ( Size jj = 1; jj <= concrete_residue->n_proton_chi(); ++jj ) {
					Size jj_protchi = concrete_residue->proton_chi_2_chi( jj );
					suggested_rotamers[ ii ]->set_chi(
						jj_protchi,
						proton_chi_chisets[ ii ]->chi[ jj_protchi ] );
				}
			}
			for ( Size ii = 1; ii <= suggested_rotamers.size(); ++ii ) {
				push_back_rotamer( suggested_rotamers[ ii ] );
			}
		} else {
			/// merely clone the input residue
			push_back_rotamer( existing_residue.clone() );
			id_for_current_rotamer_ = rotamers_.size();
		}

	}


}


void
RotamerSet_::set_extra_samples(
	task::PackerTask const & task,
	int num_10A_neighbors,
	int chi,
	chemical::ResidueTypeCOP concrete_residue,
	utility::vector1< Real > & extra_chi_steps
) const
{
	using namespace task;
	bool buried = ( num_10A_neighbors >= int(task.residue_task( resid()).extrachi_cutoff()) );
	switch ( task.residue_task( resid() ).extrachi_sample_level( buried, chi, *concrete_residue ) ) {
		case NO_EXTRA_CHI_SAMPLES :
		break;
		case EX_ONE_STDDEV :
			extra_chi_steps.push_back(1);
			extra_chi_steps.push_back(-1);
		break;
		case EX_ONE_HALF_STEP_STDDEV :
			extra_chi_steps.push_back(0.5);
			extra_chi_steps.push_back(-0.5);
		break;
		case EX_TWO_FULL_STEP_STDDEVS :
			extra_chi_steps.push_back(1);
			extra_chi_steps.push_back(2);
			extra_chi_steps.push_back(-1);
			extra_chi_steps.push_back(-2);
		break;
		case EX_TWO_HALF_STEP_STDDEVS :
			extra_chi_steps.push_back(0.5);
			extra_chi_steps.push_back(1);
			extra_chi_steps.push_back(-0.5);
			extra_chi_steps.push_back(-1);
		break;
		case EX_FOUR_HALF_STEP_STDDEVS :
			extra_chi_steps.push_back(0.5);
			extra_chi_steps.push_back(1);
			extra_chi_steps.push_back(1.5);
			extra_chi_steps.push_back(2.0);
			extra_chi_steps.push_back(-0.5);
			extra_chi_steps.push_back(-1);
			extra_chi_steps.push_back(-1.5);
			extra_chi_steps.push_back(-2);
		break;
		case EX_THREE_THIRD_STEP_STDDEVS :
			extra_chi_steps.push_back(0.33);
			extra_chi_steps.push_back(0.67);
			extra_chi_steps.push_back(1);
			extra_chi_steps.push_back(-0.33);
			extra_chi_steps.push_back(-0.67);
			extra_chi_steps.push_back(-1);
		break;
		case EX_SIX_QUARTER_STEP_STDDEVS :
			extra_chi_steps.push_back(0.25);
			extra_chi_steps.push_back(0.5);
			extra_chi_steps.push_back(0.75);
			extra_chi_steps.push_back(1);
			extra_chi_steps.push_back(1.25);
			extra_chi_steps.push_back(1.5);
			extra_chi_steps.push_back(-0.25);
			extra_chi_steps.push_back(-0.5);
			extra_chi_steps.push_back(-0.75);
			extra_chi_steps.push_back(-1);
			extra_chi_steps.push_back(-1.25);
			extra_chi_steps.push_back(-1.5);
		break;
		case ExtraRotSampleCardinality :
		default :
			std::cerr << "Error in RotamerSet_::set_extrachi_samples, invalid ExtraChiSample type" << '\n';
			utility_exit();
		break;
	}
}

/// @details The packer's 1-body is a combination of the rotamer internal energies (the
/// context dependent and independent one body energies), the intra-residue
/// energies defined by the two body energies, and the sum of the
/// two body energies with the background residues in this repacking.  The
/// PackerTask tells the RotamerSet_ what residues are part of the
/// background and which are being repacked.
void
RotamerSet_::compute_one_body_energies(
	pose::Pose const & pose,
	scoring::ScoreFunction const & sf,
	task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph,
	utility::vector1< core::PackerEnergy > & energies
) const
{
	using namespace conformation;
	using namespace scoring;

	std::fill( energies.begin(), energies.end(), core::PackerEnergy( 0.0 ) );

	int const nrotamers = num_rotamers(); // does not change in this function
	Size const theresid = resid();

	for ( int ii = 1; ii <= nrotamers; ++ii ) {
		EnergyMap emap;
		sf.eval_ci_1b( *rotamers_[ ii ], pose, emap );
		sf.eval_cd_1b( *rotamers_[ ii ], pose, emap );
		energies[ ii ] += static_cast< core::PackerEnergy > (sf.weights().dot( emap )); // precision loss here.
	}

	sf.evaluate_rotamer_intrares_energies( *this, pose, energies );

	for ( graph::Graph::EdgeListConstIter
			ir  = packer_neighbor_graph->get_node( theresid )->const_edge_list_begin(),
			ire = packer_neighbor_graph->get_node( theresid )->const_edge_list_end();
			ir != ire; ++ir ) {

		int const neighbor_id( (*ir)->get_other_ind( theresid ) );

		if ( task.pack_residue( neighbor_id ) ) continue;

		Residue const & neighbor( pose.residue( neighbor_id ) );
		sf.evaluate_rotamer_background_energies( *this, neighbor, pose, energies );

	}

	// long-range energy interactions with background
	// Iterate across the long range energy functions and use the iterators generated
	// by the LRnergy container object
	for ( ScoreFunction::LR_2B_MethodIterator
			lr_iter = sf.long_range_energies_begin(),
			lr_end  = sf.long_range_energies_end();
			lr_iter != lr_end; ++lr_iter ) {
		LREnergyContainerCOP lrec = pose.energies().long_range_container( (*lr_iter)->long_range_type() );
		if ( !lrec || lrec->empty() ) continue; // only score non-emtpy energies.

		// Potentially O(N) operation leading to O(N^2) behavior
		for ( ResidueNeighborConstIteratorOP
				rni = lrec->const_neighbor_iterator_begin( theresid ),
				rniend = lrec->const_neighbor_iterator_end( theresid );
				(*rni) != (*rniend); ++(*rni) ) {
			Size const neighbor_id = rni->neighbor_id();
			assert( neighbor_id != theresid );
			if ( task.pack_residue( neighbor_id ) ) continue;

			(*lr_iter)->evaluate_rotamer_background_energies(
				*this, pose.residue( neighbor_id ), pose, sf,
				sf.weights(), energies );

		} // (potentially) long-range neighbors of theresid [our resid()]
	} // long-range energy functions


}

/// @details Used in OptE.  Based on the function compute_one_body_energies().  OptE needs to store the energies for all score terms for each rotamer separately.  In this context there are only rotamers at a single position, with all other positions fixed.
void
RotamerSet_::compute_one_body_energy_maps(
	pose::Pose const & pose,
	scoring::ScoreFunction const & sf,
	task::PackerTask const & , // task
	graph::GraphCOP packer_neighbor_graph,
	utility::vector1< scoring::EnergyMap > & energies
) const
{
	using namespace conformation;
	using namespace scoring;

	EnergyMap default_map_of_zeroes;

	std::fill( energies.begin(), energies.end(), default_map_of_zeroes );

	int const nrotamers = num_rotamers(); // does not change in this function
	Size const theresid = resid();

	//ronj variables for surfaceE, see below for more info
	Real last_computed_surfaceE = 0.0;

	utility::vector1<Size> num_neighbors_;
	num_neighbors_.resize( pose.n_residue(), 0 );
	scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
	for ( Size id = 1; id <= pose.n_residue(); ++id ) {
		num_neighbors_[ id ] = tenA_neighbor_graph.get_node( id )->num_neighbors_counting_self();
	}

	for ( int ii = 1; ii <= nrotamers; ++ii ) {
		// one-body energies
		EnergyMap emap;
		sf.eval_ci_1b( *rotamers_[ ii ], pose, emap );
		sf.eval_cd_1b( *rotamers_[ ii ], pose, emap );

		// With apl's blessing, adding a special check for whether or not surface scoring is in use for the long-term
		// goal of trying to optimize the non-PD surface score together with the other EnergyMethods in the optE protocol.
		// Since this method is only used by optE, this ugly if statement here will not affect performance in regular runs. (ronj)

		// This calculation is inside the for loop above which goes through all the different possible rotamers at a sequence
		// position. To avoid making the expensive surface calculation, cache the last computed energy if the residue type
		// of the last rotamer is the same as that of the current rotamer since that will not change the surface score.
		// This kind of "state" can't be kept in the surface calculation function since that is not implemented as a class. (ronj)
		core::Real surface_weight( sf.get_weight( core::scoring::surface ) );
		if ( surface_weight ) {
			if ( ( ii > 1 ) && ( (*rotamers_[ii]).name() == (*rotamers_[ii-1]).name() ) ) {
				emap[ surface ] = last_computed_surfaceE;
			} else {
				pack::interaction_graph::SurfacePotential::get_instance()->compute_residue_surface_energy( *rotamers_[ii], pose, emap, theresid, num_neighbors_ );
				last_computed_surfaceE = emap[ surface ];
			}
			// emap[surface] should now have an energy in it for this rotamer
		}

		energies[ii] += emap;

		// add interactions for each rotamer with its (fixed) neighbors
		for ( graph::Graph::EdgeListConstIter
				ir  = packer_neighbor_graph->get_node( theresid )->const_edge_list_begin(),
				ire = packer_neighbor_graph->get_node( theresid )->const_edge_list_end();
				ir != ire; ++ir ) {
			int const neighbor_id( (*ir)->get_other_ind( theresid ) );
			Residue const & neighbor( pose.residue( neighbor_id ) );

			EnergyMap emap2b;

			emap2b.zero();
			sf.eval_ci_2b( neighbor, *rotamers_[ ii ], pose, emap2b );
			energies[ii] += emap2b;

			emap2b.zero();
			sf.eval_cd_2b( neighbor, *rotamers_[ ii ], pose, emap2b );
			energies[ii] += emap2b;
		}
	}

	// Intrares energies
	sf.evaluate_rotamer_intrares_energy_maps( *this, pose, energies );

	// Long Range Interactions
	for ( ScoreFunction::LR_2B_MethodIterator
			lr_iter = sf.long_range_energies_begin(),
			lr_end  = sf.long_range_energies_end();
			lr_iter != lr_end; ++lr_iter ) {

		LREnergyContainerCOP lrec = pose.energies().long_range_container( (*lr_iter)->long_range_type() );

		if ( !lrec || lrec->empty() ) continue; // only score non-empty energies.

		// Potentially O(N^2) operation...
		for ( ResidueNeighborConstIteratorOP rni = lrec->const_neighbor_iterator_begin( theresid ),
					rniend = lrec->const_neighbor_iterator_end( theresid );
					(*rni) != (*rniend); ++(*rni) ) {

			Size const neighbor_id = rni->neighbor_id();

			if( theresid == neighbor_id ) continue;

			(*lr_iter)->evaluate_rotamer_background_energy_maps(
				*this, pose.residue( neighbor_id ), pose, sf, sf.weights(), energies );

		} // (potentially) long-range neighbors of ii [our resid()]
	} // long-range energy functions
}


Size
RotamerSet_::num_rotamers() const
{
	return rotamers_.size();
}


Size
RotamerSet_::id_for_current_rotamer() const
{
	return id_for_current_rotamer_;
}


conformation::ResidueCOP
RotamerSet_::rotamer( Size rot_id ) const
{
	return rotamers_[ rot_id ];
}


conformation::Residue const &
RotamerSet_::rotamer_ref( Size rot_id ) const
{
	return *rotamers_[ rot_id ];
}


/// @details In handing out non-const data, the guarantee of rotamer-type contiguity
/// within the rotamers_ array, and the correspondence of the rotamer offset
/// data is lost.  Future access to rotamer offset data first requires an update
/// of the rotamer offset arrays.
conformation::ResidueOP
RotamerSet_::nonconst_rotamer( Size rot_id )
{
	rotamer_offsets_require_update_ = true;

	return rotamers_[ rot_id ];
}


void
RotamerSet_::store_trie(
	Size method_enum_id,
	conformation::AbstractRotamerTrieOP trie
)
{
	cached_tries_[ method_enum_id ] = trie;
}


conformation::AbstractRotamerTrieCOP
RotamerSet_::get_trie( Size method_enum_id ) const
{
	return cached_tries_[ method_enum_id ];
}

/// @details  O(n) operation; if you have a lot of rotamers you want to remove, use
/// drop_rotamers() instead.
void
RotamerSet_::drop_rotamer( Size rot_id )
{
	assert( rot_id <= rotamers_.size() );
	utility::vector1< conformation::ResidueOP > copy_rotamers( rotamers_.size() - 1, 0 );
	Size count_copy( 1 );
	for ( Size ii = 1; ii <= rotamers_.size(); ++ii ) {
		if ( ii != rot_id ) {
			copy_rotamers[ count_copy ] = rotamers_[ ii ];
			if ( ii == id_for_current_rotamer_ ) {
				id_for_current_rotamer_ = count_copy;
			}
			++count_copy;
		} else {
			if ( ii == id_for_current_rotamer_ ) {
				id_for_current_rotamer_ = 0;
			}
		}
	}
	copy_rotamers.swap( rotamers_ );
	rotamer_offsets_require_update_ = true;
	update_rotamer_offsets();

}

/// @brief rotamers_to_delete must be of size nrotmaers -- each position
/// in the array that's "true" is removed from the set of rotamers
void
RotamerSet_::drop_rotamers( utility::vector1< bool > const & rotamers_to_delete )
{
	assert( rotamers_to_delete.size() == rotamers_.size() );

	Size n_dropped = 0;
	for ( Size ii = 1; ii <= rotamers_.size(); ++ii ) {
		if ( rotamers_to_delete[ ii ] ) {
			/// if all rotamers end up dropped, then preserve the input rotamer.
			if ( ii == id_for_current_rotamer_ ) {
				current_rotamer_copy_ = rotamers_[ ii ];
			}
			rotamers_[ ii ] = 0;
			++n_dropped;
		}
	}
	if ( n_dropped == 0 ) return;

	if ( n_dropped == rotamers_.size() ) {
		if ( id_for_current_rotamer_ == 0 ) {
			utility_exit_with_message( "ERROR:: RotamerSet_::drop_rotamers attempted to remove all rotamers without available input_rotamer." );
		}
		// keep the input rotamer.
		rotamers_.resize( 1 );
		rotamers_[ 1 ] = current_rotamer_copy_;
		id_for_current_rotamer_ = 1;
		current_rotamer_copy_.reset();
	} else {
		utility::vector1< conformation::ResidueOP > new_rotamers( rotamers_to_delete.size() - n_dropped, 0 );
		Size count_new = 1;
		for ( Size ii = 1; ii <= rotamers_.size(); ++ii ) {
			if ( rotamers_[ ii ] != 0 ) {
				new_rotamers[ count_new ] = rotamers_[ ii ];
				if ( ii == id_for_current_rotamer_ ) {
					id_for_current_rotamer_ = count_new;
				}
				++count_new;
			}
		}
		new_rotamers.swap( rotamers_ );
	}
	rotamer_offsets_require_update_ = true;
	update_rotamer_offsets();
}

/// @brief deletes the rotamers in the list with the given indices.
/// The indices of these rotamers is presumed to be those before any delete operation.
/// e.g. if there are four rotamers, and rotamer_indices_to_delete includes 1 & 3,
/// then the rotamers that will remain are the rotamers originally indexed as 2 and 4,
/// even though their new indices will be 1 & 2.
void
RotamerSet_::drop_rotamers_by_index(
	utility::vector1< Size > const & rotamer_indices_to_delete
)
{
	utility::vector1< bool > rotamers_to_delete( rotamers_.size(), false );
	for ( Size ii = 1; ii <= rotamer_indices_to_delete.size(); ++ii ) {
		rotamers_to_delete[ rotamer_indices_to_delete[ ii ] ] = true;
	}
	drop_rotamers( rotamers_to_delete );
}

/// @details Bump check does not include long range energies,
/// though, maybe this should change.
core::PackerEnergy
RotamerSet_::bump_check(
	ResidueCOP rotamer,
	scoring::ScoreFunction const & sf,
	pose::Pose const & pose,
	task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph
) const
{
	using namespace scoring;
	using namespace conformation;

	EnergyMap emap;

	for ( graph::Graph::EdgeListConstIter
			ir  = packer_neighbor_graph->get_node( resid() )->const_edge_list_begin(),
			ire = packer_neighbor_graph->get_node( resid() )->const_edge_list_end();
			ir != ire; ++ir ) {
		int const neighbor_id( (*ir)->get_other_ind( resid() ) );
		Residue const & neighbor( pose.residue( neighbor_id ) );

		if ( ! task.pack_residue( neighbor_id ) ) {
			sf.bump_check_full( *rotamer, neighbor, pose, emap);
		} else {
			sf.bump_check_backbone( *rotamer, neighbor, pose, emap);
		}
	}
	return static_cast< core::PackerEnergy > (sf.weights().dot( emap ));
}

/// @details refactored into its own method so that it could be used in both regular packing
/// and also in optimizeH.
void
RotamerSet_::build_tp3_water_rotamers(
		pose::Pose const & pose,
		task::PackerTask const & task,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const & existing_residue,
		graph::GraphCOP packer_neighbor_graph
)
{

	// build rotamers for water
	utility::vector1< ResidueOP > new_rotamers;

	build_independent_water_rotamers( resid(), *concrete_residue, task, pose, packer_neighbor_graph, new_rotamers );

	for ( Size ii=1; ii<= new_rotamers.size(); ++ii ) {
		new_rotamers[ii]->seqpos( resid() );
		new_rotamers[ii]->chain( existing_residue.chain() );
		push_back_rotamer( new_rotamers[ii] );
	}

	if ( task.include_current( resid() ) && existing_residue.name() == concrete_residue->name() ) {
		ResidueOP rot = existing_residue.create_rotamer();
		push_back_rotamer( rot );
		id_for_current_rotamer_ = num_rotamers();
	}
}

void
RotamerSet_::prepare_for_new_residue_type( core::chemical::ResidueType const & restype )
{
	if ( n_residue_types_ == 0 ) {
		new_residue_type();
		new_residue_group();
		return;
	}
	if ( num_rotamers() == 0 ) {
		// n_residue_types_ and n_residue_groups_ is not zero -- how odd
		return;
	}

	if ( different_restype( rotamers_[ num_rotamers() ]->type(), restype )) {
		new_residue_type();
	}
	if (  different_resgroup( rotamers_[ num_rotamers() ]->type(), restype )) {
		new_residue_group();
	}
}


bool
RotamerSet_::different_restype( core::chemical::ResidueType const & rt1, core::chemical::ResidueType const & rt2 ) const
{
	return & rt1 != & rt2;
}

/// @details The logic to determine if two residue types should be classified as part of the same group.
/// The thinking is as follows.  Two residue types are in the same group if they have the same residue type.
/// They're in the same group if their residue types differ, but they have the same name3 (HIS vs HIS_D have
/// the same name3) and they have the same neighbor radius (SER and PhosphoSER should have different groups).
/// The goal is to organize residue types together which will be packed together (as happens in multistate design
/// with HIS and HISD) and that have the same reach (as is needed for the AANeighborSparseMatrix).
bool
RotamerSet_::different_resgroup( core::chemical::ResidueType const & rt1, core::chemical::ResidueType const & rt2 ) const
{
	return & rt1 != & rt2 && ( rt1.name3() != rt2.name3() || rt1.nbr_radius() != rt2.nbr_radius() );
}


void
RotamerSet_::new_residue_type()
{
	++n_residue_types_;
	residue_type_rotamers_begin_.push_back( num_rotamers() + 1);
	n_rotamers_for_restype_.push_back( 0 );
}


void
RotamerSet_::new_residue_group()
{
	++n_residue_groups_;
	residue_group_rotamers_begin_.push_back( num_rotamers() + 1 );
	n_rotamers_for_resgroup_.push_back( 0 );
}


void
RotamerSet_::push_back_rotamer( conformation::ResidueOP rotamer )
{
	rotamers_.push_back( rotamer );
	residue_type_for_rotamers_.push_back( n_residue_types_ );
	residue_group_for_rotamers_.push_back( n_residue_groups_ );
	++n_rotamers_for_restype_[ n_residue_types_ ];
	++n_rotamers_for_resgroup_[ n_residue_groups_ ];
}


void
RotamerSet_::build_dependent_rotamers(
	RotamerSets const & rotamer_sets,
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	task::PackerTask const & the_task,
	graph::GraphCOP packer_neighbor_graph
)
{
	using namespace chemical;
	conformation::Residue const & existing_residue( pose.residue( resid() ) );
	for ( task::ResidueLevelTask::ResidueTypeCOPListConstIter
			allowed_iter = the_task.residue_task( resid() ).allowed_residue_types_begin(),
			allowed_end = the_task.residue_task( resid() ).allowed_residue_types_end();
			allowed_iter != allowed_end; ++allowed_iter ) {
		build_dependent_rotamers_for_concrete( rotamer_sets, pose, scorefxn, the_task,
			existing_residue, *allowed_iter, packer_neighbor_graph);
	}
}


void
RotamerSet_::build_dependent_rotamers_for_concrete(
	RotamerSets const & rotamer_sets,
	pose::Pose const & pose,
	scoring::ScoreFunction const &,// scorefxn,
	task::PackerTask const & task,
	conformation::Residue const & existing_residue,
	chemical::ResidueTypeCOP concrete_residue,
	graph::GraphCOP packer_neighbor_graph
)
{
	using namespace conformation;
	using namespace pack::task;

	//if ( !task.residue_task( resid() ).build_dependent_rotamers() ) return;

	if ( concrete_residue->name() != "TP3" ) return; // logic only exists for this guy right now

	if ( concrete_residue->name() == "TP3" ) { // TIP3 water /////////////////////////////////

		// build rotamers for water
		utility::vector1< ResidueOP > new_rotamers;

		build_dependent_water_rotamers(
			rotamer_sets, resid(), *concrete_residue,
			task, pose, packer_neighbor_graph,
			new_rotamers );

		if ( new_rotamers.empty() ) return;

		prepare_for_new_residue_type( *concrete_residue );

		for ( Size ii=1; ii<= new_rotamers.size(); ++ii ) {
			new_rotamers[ii]->seqpos( resid() );
			new_rotamers[ii]->chain( existing_residue.chain() );
			push_back_rotamer( new_rotamers[ii] );
		}

	} else {
		utility_exit_with_message( "unsupported restype for dependent rotamer building: "+concrete_residue->name() );
	}
}


void
RotamerSet_::update_rotamer_offsets() const
{
	if ( ! rotamer_offsets_require_update_ ) return;

	if ( rotamers_.size() == 0 ) {
		n_residue_types_ = 0;
		n_residue_groups_ = 0;
		residue_type_for_rotamers_.resize( 0 );
		residue_group_for_rotamers_.resize( 0 );
		residue_type_rotamers_begin_.resize( 0 );
		residue_group_rotamers_begin_.resize( 0 );
		n_rotamers_for_restype_.resize( 0 );
		n_rotamers_for_resgroup_.resize( 0 );
		return;
	}

	/// From here forward, rotamers_.size() >= 1
	residue_type_for_rotamers_.resize( rotamers_.size() );
	residue_group_for_rotamers_.resize( rotamers_.size() );
	n_residue_types_ = 1;
	n_residue_groups_ = 1;
	residue_type_for_rotamers_[ 1 ] = n_residue_types_;
	residue_group_for_rotamers_[ 1 ] = n_residue_groups_;
	for ( Size ii = 2; ii <= rotamers_.size(); ++ii ) {
		// compare addresses of the two types
		// treat them as different amino acids only if they have different name3's
		// or if they have different radii
		//if ( & (rotamers_[ ii ]->type()) != & (rotamers_[ ii ]->type()) ) {
		if ( different_restype( rotamers_[ ii ]->type(), rotamers_[ ii-1 ]->type() ) ) {
			++n_residue_types_;
		}
		residue_type_for_rotamers_[ ii ] = n_residue_types_;

		if ( different_resgroup( rotamers_[ ii ]->type(), rotamers_[ ii-1 ]->type() ) ) {
			++n_residue_groups_;
		}
		residue_group_for_rotamers_[ ii ] = n_residue_groups_;
	}

	residue_type_rotamers_begin_.resize( n_residue_types_ );
	n_rotamers_for_restype_.resize( n_residue_types_ );
	std::fill( residue_type_rotamers_begin_.begin(), residue_type_rotamers_begin_.end(), 0 );
	std::fill( n_rotamers_for_restype_.begin(), n_rotamers_for_restype_.end(), 0 );

	residue_group_rotamers_begin_.resize( n_residue_groups_ );
	n_rotamers_for_resgroup_.resize( n_residue_groups_ );
	std::fill( residue_group_rotamers_begin_.begin(), residue_group_rotamers_begin_.end(), 0 );
	std::fill( n_rotamers_for_resgroup_.begin(), n_rotamers_for_resgroup_.end(), 0 );

	Size count_seen_residue_types( 1 );
	Size count_seen_residue_groups( 1 );
	n_rotamers_for_restype_[ count_seen_residue_types ] = 1;
	n_rotamers_for_resgroup_[ count_seen_residue_groups ] = 1;
	residue_type_rotamers_begin_[ count_seen_residue_types ] = 1;
	residue_group_rotamers_begin_[ count_seen_residue_groups ] = 1;

	for ( Size ii = 2; ii <= rotamers_.size(); ++ii ) {
		if ( residue_type_for_rotamers_[ ii ] != residue_type_for_rotamers_[ ii-1 ] ) {
			++count_seen_residue_types;
			residue_type_rotamers_begin_[ count_seen_residue_types ] = ii;
		}
		++n_rotamers_for_restype_[ count_seen_residue_types ];
		if ( residue_group_for_rotamers_[ ii ] != residue_group_for_rotamers_[ ii-1 ] ) {
			++count_seen_residue_groups;
			residue_group_rotamers_begin_[ count_seen_residue_groups ] = ii;
		}
		++n_rotamers_for_resgroup_[ count_seen_residue_groups ];
	}
	//std::cout << "nrestypes " << n_residue_types_ << std::endl;
	rotamer_offsets_require_update_ = false;
}

} // rotamer_set
} // pack
} // core
