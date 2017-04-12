// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/RotamerSet/ContinuousRotamerSet.cc
/// @brief  Implementation of the ContinuousRotamerSet and ContinuousRotamerSets classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/rotamer_set/ContinuousRotamerSet.hh>

// Package Headers
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/task/PackerTask.hh>

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/exit.hh>

#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace rotamer_set {

/// @details Auto-generated virtual destructor
ContinuousRotamerSets::~ContinuousRotamerSets() {}


ContinuousRotamerSet::ContinuousRotamerSet() :
	resid_( 0 ),
	input_rotamer_samplingrot_index_( 0 ),
	input_rotamer_rotblock_( 0 ),
	n_restypes_( 0 ),
	n_baserots_total_( 0 ),
	n_samplingrots_total_( 0 )
{}

ContinuousRotamerSet::~ContinuousRotamerSet() {}

void ContinuousRotamerSet::build_rotamers(
	pose::Pose const & pose,
	Size resid,
	task::PackerTask const & task
)
{
	resid_ = resid;

	/// Need to query the rotamer library for each available ResidueType that
	/// the packer task suggests should be available at this position.  For
	/// the ResidueTypes that have a dunbrack library, ask the library for the
	/// complete set of rotamer samples.  For the residue types that do not have
	/// a rotamer library (aka RotamerLibrary::get_rsd_library returns NULL),
	/// store that fact by keeping its vector in the samples_ array at size 0.
	/// NOTE: it's possible for there to be two rotamers at a position without
	/// a library being present: the first for the native rotamer, the second
	/// for the ideal rotamer.

	Size count_allowed( 0 );
	for ( task::ResidueLevelTask::ResidueTypeCOPListConstIter
			allowed_iter = task.residue_task( resid ).allowed_residue_types_begin(),
			allowed_end = task.residue_task( resid ).allowed_residue_types_end();
			allowed_iter != allowed_end; ++allowed_iter ) {
		++count_allowed;
	}
	if ( count_allowed == 0 ) {
		// save the input rotamer
	} else {
		n_restypes_ = count_allowed;
		aa_for_rotblock_.resize( n_restypes_ );
		restype_for_rotblock_.resize( n_restypes_ );

		n_baserotamers_for_rotblock_.resize( n_restypes_ );
		baserots_offsets_.resize( n_restypes_ );

		n_samplingrots_for_rotblock_.resize( n_restypes_ );
		samplingrot_offsets_.resize( n_restypes_ );

		samples_.resize( n_restypes_ );
	}

	Size count_restype_ind( 0 );
	for ( task::ResidueLevelTask::ResidueTypeCOPListConstIter
			allowed_iter = task.residue_task( resid ).allowed_residue_types_begin(),
			allowed_end = task.residue_task( resid ).allowed_residue_types_end();
			allowed_iter != allowed_end; ++allowed_iter ) {
		++count_restype_ind;
		aa_for_rotblock_[ count_restype_ind ] = (*allowed_iter)->aa();
		restype_for_rotblock_[ count_restype_ind ] = (*allowed_iter);
		rotamers::SingleResidueRotamerLibraryCOP rotlib =
			rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( **allowed_iter );

		if ( rotlib ) {
			// OK -- two options -- we're dealing with a Dunbrack library, in which case, we should
			// store DunbrackRotamerSampleData, or, we're dealing with some other kind of library,
			// in which case, we need to quit (TO DO: support ligand rotamers)
			dunbrack::SingleResidueDunbrackLibraryCOP dunlib =
				utility::pointer::dynamic_pointer_cast<dunbrack::SingleResidueDunbrackLibrary const > ( rotlib );
			if ( dunlib ) {
				// amw
				utility::fixedsizearray1< Real, 5 > bbs;
				bbs[ 1 ] = pose.phi( resid );
				bbs[ 2 ] = pose.psi( resid );
				samples_[ count_restype_ind ] = dunlib->get_all_rotamer_samples(bbs );// pose.phi( resid ), pose.psi( resid ) );
			} else {
				utility_exit_with_message("ContinuousRotamerSet cannot support non-Dunbrack rotamer libraries at this time" );
			}

			utility::fixedsizearray1< Real, 5 > bbs;
			// Explicit for alpha or beta.
			// Otherwise... we need a mainchain torsion sense that is
			// somehow not modified by terminal types!
			if ( pose.residue( resid ).type().is_beta_aa() ) {
				bbs[ 1 ] = pose.phi(   resid );
				bbs[ 2 ] = pose.theta( resid );
				bbs[ 3 ] = pose.psi(   resid );
			} else {
				bbs[ 1 ] = pose.phi( resid );
				bbs[ 2 ] = pose.psi( resid );
			}
			samples_[ count_restype_ind ] = dunlib->get_all_rotamer_samples( bbs );

		} else {
			// Ala or gly or something similar.
			// noop
		}

		if ( (*allowed_iter).get() == & pose.residue_type( resid ) && task.residue_task( resid ).include_current() ) {
			/// save coordinates for this residue
			input_rotamer_rotblock_ = count_restype_ind;
			conformation::Residue const & res = pose.residue( resid );
			input_rotamer_coords_.resize( res.natoms() );
			for ( Size ii = 1; ii <= res.natoms(); ++ii ) {
				input_rotamer_coords_[ ii ] = res.xyz( ii );
			}
		}
	}

	for ( Size ii = 1; ii <= n_restypes_; ++ii ) {
		determine_rotcounts_for_restype( ii );
		if ( ii != 1 ) {
			baserots_offsets_[ ii ] = baserots_offsets_[ ii - 1 ] + n_baserotamers_for_rotblock_[ ii - 1 ];
			samplingrot_offsets_[ ii ] = samplingrot_offsets_[ ii - 1 ] + n_samplingrots_for_rotblock_[ ii - 1 ];
		} else {
			baserots_offsets_[ ii ] = 0;
			samplingrot_offsets_[ ii ] = 0;
		}
		n_baserots_total_ += n_baserotamers_for_rotblock_[ ii ];
		n_samplingrots_total_ += n_samplingrots_for_rotblock_[ ii ];
		if ( ii == input_rotamer_rotblock_ ) {
			// mark the input rotamer as the last rotamer of the samplingrots for its restype.
			input_rotamer_samplingrot_index_ = n_samplingrots_for_rotblock_[ ii ] + samplingrot_offsets_[ ii ];
		}
	}

	aa_block_for_baserotamer_.resize( n_baserots_total_ );
	aa_block_for_samplingrotamer_.resize( n_samplingrots_total_ );
	for ( Size ii = 1; ii <= n_restypes_; ++ii ) {
		for ( Size jj = baserots_offsets_[ ii ] + 1; jj <= baserots_offsets_[ ii ] + n_baserotamers_for_rotblock_[ ii ]; ++jj ) {
			aa_block_for_baserotamer_[ jj ] = ii;
		}
		for ( Size jj = samplingrot_offsets_[ ii ] + 1; jj <= samplingrot_offsets_[ ii ] + n_samplingrots_for_rotblock_[ ii ]; ++jj ) {
			aa_block_for_samplingrotamer_[ jj ] = ii;
		}
	}
}

/// Note, its entirely possible for two residue types to be the "same amino acid",
/// e.g. HIS and HIS_D.
Size
ContinuousRotamerSet::get_n_residue_types() const
{
	return n_restypes_;
}


//Size
//ContinuousRotamerSet::get_residue_type_begin( Size which_restype ) const
//{
// return n_totalrots_offsets_[ which_restype ] + 1;
//}

Size
ContinuousRotamerSet::get_n_sampling_rotamers_for_rotblock( Size which_restype ) const
{
	return n_samplingrots_for_rotblock_[ which_restype ];
}

Size
ContinuousRotamerSet::get_n_baserotamers_for_rotblock( Size which_restype ) const
{
	return samples_[ which_restype ].size();
}

Size
ContinuousRotamerSet::get_rotblock_index_for_sampling_rotamer( Size which_rotamer ) const
{
	return aa_block_for_samplingrotamer_[ which_rotamer ];
}

Size
ContinuousRotamerSet::num_base_rotamers_total() const { return n_baserots_total_; }

Size
ContinuousRotamerSet::num_sampling_rotamers_total() const { return n_samplingrots_total_; }

Size
ContinuousRotamerSet::sampling_id_for_current_rotamer() const
{
	return input_rotamer_samplingrot_index_;
}

utility::vector1< Vector > const &
ContinuousRotamerSet::current_rotamer_coords() const
{
	return input_rotamer_coords_;
}

Size
ContinuousRotamerSet::resid() const { return resid_;}

chemical::ResidueTypeCOP
ContinuousRotamerSet::restype_for_rotblock( Size rotblock ) const
{
	return restype_for_rotblock_[ rotblock ];
}


dunbrack::DunbrackRotamerSampleData const &
ContinuousRotamerSet::baserotamer_data( Size rotblock_ind, Size rotid_for_aa ) const
{
	return samples_[ rotblock_ind ][ rotid_for_aa ];
}

Size
ContinuousRotamerSet::pick_baserotamer_from_rotblock( Size rotblock_ind, Real rand_btw_0_and_1 ) const
{
	debug_assert( rand_btw_0_and_1 >= 0.0 && rand_btw_0_and_1 <= 1.0 );
	Real accumulated_prob( 0.0 );
	Size const nsamples = samples_[ rotblock_ind ].size();
	for ( Size ii = 1; ii <= nsamples; ++ii ) {
		accumulated_prob += samples_[ rotblock_ind ][ ii ].probability();
		//debug_assert( accumulated_prob <= 1.0 + 1e-5 ); // should never exceed 1
		if ( accumulated_prob > rand_btw_0_and_1 ) return ii;
	}
	return nsamples;
}

void
ContinuousRotamerSet::determine_rotcounts_for_restype( Size restype_ind )
{
	bool restype_is_current = restype_ind == input_rotamer_rotblock_;
	if ( samples_[ restype_ind ].size() == 0 ) {
		// ala or gly
		n_baserotamers_for_rotblock_[ restype_ind ] = 1;
		n_samplingrots_for_rotblock_[ restype_ind ] = 1;
	} else {
		Real const prob_to_accumulate = 0.98;
		n_baserotamers_for_rotblock_[ restype_ind ] = samples_[ restype_ind ].size();
		Real cumulative_probability( 0.0 );
		for ( Size ii = 1; ii <= n_baserotamers_for_rotblock_[ restype_ind ]; ++ii ) {
			n_samplingrots_for_rotblock_[ restype_ind ] += 1;
			cumulative_probability += samples_[ restype_ind ][ ii ].probability();
			if ( cumulative_probability > prob_to_accumulate ) break;
		}
		n_samplingrots_for_rotblock_[ restype_ind ] = samples_[ restype_ind ].size();
	}

	/// Multiply the number of totalrots by the combinations of proton chi samples.
	/// this can be computed independently of whether we have a RotamerLibrary or not;
	/// though, neither alanine nor glycine have proton chi...
	for ( Size ii = 1; ii <= restype_for_rotblock_[ restype_ind ]->n_proton_chi(); ++ii ) {
		n_samplingrots_for_rotblock_[ restype_ind ] *=
			restype_for_rotblock_[ restype_ind ]->proton_chi_samples( ii ).size();
	}

	if ( restype_is_current ) {
		n_samplingrots_for_rotblock_[ restype_ind ] += 1;
	}

}

/////////////// ContinuousRotamerSets ////////////////////

ContinuousRotamerSets::ContinuousRotamerSets(
	pose::Pose const & pose,
	task::PackerTask const & task
)
{
	total_residue_ = pose.size();
	nmoltenres_ = task.num_to_be_packed();
	moltenresid_2_resid_.resize( nmoltenres_ );
	resid_2_moltenresid_.resize( total_residue_, 0 );
	rotamer_sets_.resize( nmoltenres_ );
	moltenres_sample_rot_offset_.resize( nmoltenres_ );
	Size count_moltenres( 0 );
	for ( Size ii = 1; ii <= total_residue_; ++ii ) {
		if ( ! task.being_packed( ii ) ) continue;
		++count_moltenres;
		moltenresid_2_resid_[ count_moltenres ] = ii;
		resid_2_moltenresid_[ ii ] = count_moltenres;
		rotamer_sets_[ count_moltenres ].build_rotamers( pose, ii, task );
	}

	n_sample_rotamers_ = 0;
	for ( Size ii = 1; ii <= nmoltenres_; ++ii ) {
		moltenres_sample_rot_offset_[ ii ] = n_sample_rotamers_;
		n_sample_rotamers_ += rotamer_sets_[ ii ].num_sampling_rotamers_total();
	}
	moltenres_for_sample_rot_.resize( n_sample_rotamers_, 0 );
	Size count_sample_rot = 1;
	for ( Size ii = 1; ii <= nmoltenres_; ++ii ) {
		Size const ii_nsamples = rotamer_sets_[ ii ].num_sampling_rotamers_total();
		for ( Size jj = 1; jj <= ii_nsamples; ++jj, ++count_sample_rot ) {
			moltenres_for_sample_rot_[ count_sample_rot ] = ii;
		}
	}

}

Size ContinuousRotamerSets::nmoltenres() const { return nmoltenres_; }
Size ContinuousRotamerSets::total_residue() const { return total_residue_; }

Size ContinuousRotamerSets::resid_2_moltenresid( Size resid ) const  { return resid_2_moltenresid_[ resid ]; }
Size ContinuousRotamerSets::moltenresid_2_resid( Size mresid ) const { return moltenresid_2_resid_[ mresid]; }

ContinuousRotamerSet const &
ContinuousRotamerSets::rotamer_set_for_moltenres( Size mresid ) const
{
	return rotamer_sets_[ mresid ];
}

ContinuousRotamerSet const &
ContinuousRotamerSets::rotamer_set_for_res( Size resid ) const
{
	return rotamer_sets_[ resid_2_moltenresid_[ resid ]];
}

Size
ContinuousRotamerSets::n_sample_rotamers() const
{
	return n_sample_rotamers_;
}

Size
ContinuousRotamerSets::moltenres_for_sample_rot( Size sample_rotno ) const
{
	return moltenres_for_sample_rot_[ sample_rotno ];
}

Size
ContinuousRotamerSets::full_sample_rot_index_2_moltenres_sample_rot_index( Size sample_rotno ) const
{
	return sample_rotno - moltenres_sample_rot_offset_[ moltenres_for_sample_rot_[ sample_rotno ] ];
}


} // namespace rotamer_set
} // namespace pack
} // namespace core


