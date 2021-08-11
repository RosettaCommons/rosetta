// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/data_storage/PoseWithTerminalSegmentsOfKnownDSSP.cc
/// @brief a region of a Pose with secondary structures of known DSSP at either teminus
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/pose_sewing/data_storage/PoseWithTerminalSegmentsOfKnownDSSP.hh>

#include <core/pose/subpose_manipulation_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.pose_sewing.data_storage.PoseWithTerminalSegmentsOfKnownDSSP" );


namespace protocols {
namespace pose_sewing {
namespace data_storage {

PoseWithTerminalSegmentsOfKnownDSSP::PoseWithTerminalSegmentsOfKnownDSSP():
	utility::VirtualBase()
{
	C_term_DSSP_ = 'X';
	N_term_DSSP_ = 'X';
}

PoseWithTerminalSegmentsOfKnownDSSP::~PoseWithTerminalSegmentsOfKnownDSSP(){}

PoseWithTerminalSegmentsOfKnownDSSP::PoseWithTerminalSegmentsOfKnownDSSP( PoseWithTerminalSegmentsOfKnownDSSP const &) = default;


PoseWithTerminalSegmentsOfKnownDSSPOP
PoseWithTerminalSegmentsOfKnownDSSP::clone() const {
	return PoseWithTerminalSegmentsOfKnownDSSPOP( new PoseWithTerminalSegmentsOfKnownDSSP( *this ) );
}
//getters and setters
void
PoseWithTerminalSegmentsOfKnownDSSP::set_N_term_DSSP (char new_nterm_dssp){
	N_term_DSSP_ = new_nterm_dssp;
}

char
PoseWithTerminalSegmentsOfKnownDSSP::get_N_term_DSSP() const {
	return N_term_DSSP_;
}

void
PoseWithTerminalSegmentsOfKnownDSSP::set_C_term_DSSP (char new_cterm_dssp){
	C_term_DSSP_ = new_cterm_dssp;
}

char
PoseWithTerminalSegmentsOfKnownDSSP::get_C_term_DSSP() const {
	return C_term_DSSP_;
}

void
PoseWithTerminalSegmentsOfKnownDSSP::set_N_term_length ( core::Size new_nterm_length){
	N_term_length_ = new_nterm_length;
}

core::Size
PoseWithTerminalSegmentsOfKnownDSSP::get_N_term_length () const {
	return N_term_length_;
}

void
PoseWithTerminalSegmentsOfKnownDSSP::set_C_term_length ( core::Size new_cterm_length){
	C_term_length_ = new_cterm_length;
}

core::Size
PoseWithTerminalSegmentsOfKnownDSSP::get_C_term_length () const {
	return C_term_length_;
}

void
PoseWithTerminalSegmentsOfKnownDSSP::set_filename(std::string const & new_filename){
	filename_ = new_filename;
}

void
PoseWithTerminalSegmentsOfKnownDSSP::set_segfile_path(const std::string &segfile_path){
	segfile_path_ = segfile_path;
}

std::string
PoseWithTerminalSegmentsOfKnownDSSP::get_segfile_path() const {
	return segfile_path_;
}

std::string
PoseWithTerminalSegmentsOfKnownDSSP::get_filename() const {
	return filename_;
}

void
PoseWithTerminalSegmentsOfKnownDSSP::set_secstruct(std::string const & new_secstruct){
	secstruct_ = new_secstruct;
}

std::string
PoseWithTerminalSegmentsOfKnownDSSP::get_secstruct() const {
	return secstruct_;
}

void
PoseWithTerminalSegmentsOfKnownDSSP::set_source_pose ( core::pose::PoseCOP new_source){
	source_pose_ = new_source;
}

core::pose::PoseOP
PoseWithTerminalSegmentsOfKnownDSSP::create_source_pose_for_segment(bool add_elements){
	core::pose::PoseOP out_pose = core::pose::PoseOP(new core::pose::Pose());
	core::chemical::ResidueTypeSetCOP res_type_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	core::pose::make_pose_from_sequence(*out_pose, source_pose_->sequence(),res_type_set,true,false);
	out_pose->copy_segment(out_pose->size(),*(source_pose_),1,1);

	for ( core::Size current_residue = 1; current_residue <= out_pose->size(); ++current_residue ) {
		out_pose->set_secstruct(current_residue, source_pose_->secstruct(current_residue));
	}

	if ( add_elements ) {
		out_pose->pdb_info( utility::pointer::make_shared< core::pose::PDBInfo >( *out_pose ) );
	}

	return out_pose;
}

void
PoseWithTerminalSegmentsOfKnownDSSP::store_source_pose_for_segment ( std::string const & pose_filename, bool store_mmTF /*true*/, bool  add_elements/*false*/){

	core::pose::PoseOP out_pose = create_source_pose_for_segment(add_elements);

	if ( store_mmTF ) {
		out_pose->dump_mmtf( pose_filename );
	} else {
		out_pose->dump_pdb(pose_filename );
	}
}

void
PoseWithTerminalSegmentsOfKnownDSSP::store_reference_pdb ( std::string pose_filename){
	core::pose::PoseOP out_pose = core::pose::PoseOP(new core::pose::Pose());
	core::chemical::ResidueTypeSetCOP res_type_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	core::pose::make_pose_from_sequence(*out_pose, source_pose_->sequence(),res_type_set,true,false);
	out_pose->copy_segment(out_pose->size(),*(source_pose_),1,1);
	for ( core::Size current_residue = 1; current_residue <= out_pose->size(); ++current_residue ) {
		out_pose->set_secstruct(current_residue, source_pose_->secstruct(current_residue));
	}
	//out_pose->dump_mmtf(pose_filename);
	out_pose->dump_scored_pdb(pose_filename, *(new core::scoring::ScoreFunction()));
}


core::pose::PoseCOP
PoseWithTerminalSegmentsOfKnownDSSP::get_source_pose(bool clone_if_new /*true*/){
	return get_source_pose_op(clone_if_new);
}

core::pose::PoseOP
PoseWithTerminalSegmentsOfKnownDSSP::get_source_pose_op(bool clone_if_new /*true*/){
	if ( source_pose_==nullptr ) {
		core::pose::PoseOP in_pose = core::import_pose::pose_from_file(filename_);
		for ( core::Size current_residue = 1; current_residue <= in_pose->size(); ++current_residue ) {
			in_pose->set_secstruct(current_residue,secstruct_[current_residue-1]);
		}
		core::pose::remove_variant_type_from_pose_residue(*in_pose, core::chemical::LOWER_TERMINUS_VARIANT, 1);
		core::pose::remove_variant_type_from_pose_residue(*in_pose, core::chemical::UPPER_TERMINUS_VARIANT, in_pose->size());
		if ( clone_if_new ) {
			source_pose_ = in_pose->clone();
		} else {
			source_pose_= in_pose;//Set it instead of cloning.  This removes an extra clone operation, which is expensive.
		}
		return in_pose;
	}
	return source_pose_->clone();

}

} //protocols
} //pose_sewing
} //data_storage






