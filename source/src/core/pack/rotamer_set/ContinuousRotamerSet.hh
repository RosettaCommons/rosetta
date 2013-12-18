// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/RotamerSet/ContinuousRotamerSet.hh
/// @brief  Declaration for the class which desribes rotamers statistically instead of as discrete samples
///         For use in stochastic packing (see Feng & Dokholyan, 2006)
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_rotamer_set_ContinuousRotamerSet_HH
#define INCLUDED_core_pack_rotamer_set_ContinuousRotamerSet_HH

// Unit Headers
#include <core/pack/rotamer_set/ContinuousRotamerSet.fwd.hh>

// Package Headers
// AUTO-REMOVED #include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/task/PackerTask.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
//#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#include <core/chemical/ResidueType.fwd.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>


namespace core {
namespace pack {
namespace rotamer_set {


class ContinuousRotamerSet : public utility::pointer::ReferenceCount
{
private: // typedefs
	typedef utility::pointer::ReferenceCount parent;
	typedef utility::vector1< dunbrack::DunbrackRotamerSampleData >  RotamersForAA;

public:
	ContinuousRotamerSet();
	virtual ~ContinuousRotamerSet();

	void build_rotamers(
		pose::Pose const & pose,
		Size resid,
		task::PackerTask const & task
	);

	Size
	get_n_residue_types() const;

	//Size
	//get_residue_type_begin( Size which_restype ) const;

	Size
	get_n_sampling_rotamers_for_rotblock( Size which_restype ) const;

	Size
	get_n_baserotamers_for_rotblock( Size which_restype ) const;

	/// @brief Rotamers i to i+j of all the same residue type are grouped together.
	/// This function returns the index of the residue type in a contiguous block
	/// of rotamers.  E.g. rotamers 100 to 120 might all be lysine rotamers, and might
	/// be the 8th residue type, with the first 7 residue types spanning rotamers 1 to 99.
	/// If new lysine rotamers are appended to the end of the rotamer set, they are
	/// considered to be in a separate residue type block.  Lysine rotamers 200 to 210 might
	/// be block 15 while lysine rotamers 100 to 120 are still block 7.
	Size
	get_rotblock_index_for_sampling_rotamer( Size which_rotamer ) const;

	Size
	num_base_rotamers_total() const;

	Size
	num_sampling_rotamers_total() const;

	Size
	sampling_id_for_current_rotamer() const;

	utility::vector1< Vector > const &
	current_rotamer_coords() const;
	
	Size resid() const;

	chemical::ResidueTypeCOP
	restype_for_rotblock( Size rotblock ) const;

	dunbrack::DunbrackRotamerSampleData const &
	baserotamer_data( Size aa_ind, Size rotid_for_aa ) const;

	Size
	pick_baserotamer_from_rotblock( Size aa_ind, Real rand_btw_0_and_1 ) const;

private:
	void determine_rotcounts_for_restype( Size );

private:

	Size resid_; //which residue is this?

	Size input_rotamer_samplingrot_index_;
	Size input_rotamer_rotblock_;
	utility::vector1< Vector > input_rotamer_coords_;

	Size n_restypes_; // aka n_rotblocks
	utility::vector1< chemical::AA >              aa_for_rotblock_;
	utility::vector1< chemical::ResidueTypeCOP >  restype_for_rotblock_;

	Size n_baserots_total_; // excluding proton chis, how many dihedral angle bins are there?
	utility::vector1< Size >          n_baserotamers_for_rotblock_;
	utility::vector1< Size >          baserots_offsets_;
	utility::vector1< Size >          aa_block_for_baserotamer_;

	Size n_samplingrots_total_; // // how many rotamers above the 98% probability threshold?
	utility::vector1< Size >          n_samplingrots_for_rotblock_; //top 98% by prob * nproton_chi_extra
	utility::vector1< Size >          samplingrot_offsets_;
	utility::vector1< Size >          aa_block_for_samplingrotamer_;

	utility::vector1< RotamersForAA > samples_;

};

class ContinuousRotamerSets : public utility::pointer::ReferenceCount
{
public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~ContinuousRotamerSets();
	ContinuousRotamerSets(
		pose::Pose const & pose,
		task::PackerTask const & task
	);

	Size nmoltenres() const;
	Size total_residue() const;

	Size resid_2_moltenresid( Size ) const;
	Size moltenresid_2_resid( Size ) const;

	ContinuousRotamerSet const &
	rotamer_set_for_moltenres( Size ) const;

	ContinuousRotamerSet const &
	rotamer_set_for_res( Size ) const;

	Size
	n_sample_rotamers() const;

	Size
	moltenres_for_sample_rot( Size sample_rotno ) const;

	Size
	full_sample_rot_index_2_moltenres_sample_rot_index( Size sample_rotno ) const;

private:
	Size total_residue_;
	Size nmoltenres_;
	utility::vector1< Size > moltenresid_2_resid_;
	utility::vector1< Size > resid_2_moltenresid_;
	utility::vector1< ContinuousRotamerSet > rotamer_sets_;
	Size n_sample_rotamers_;
	utility::vector1< Size > moltenres_for_sample_rot_;
	utility::vector1< Size > moltenres_sample_rot_offset_;
};


} // namespace rotamer_set
} // namespace pack
} // namespace core


#endif // INCLUDED_core_pack_RotamerSet_RotamerSet_HH
