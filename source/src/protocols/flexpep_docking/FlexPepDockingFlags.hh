// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file   FlexPepDockingFlags.hh
///
/// @brief flags structure for FlexPepDocking protocols
/// @date January 1, 2009
/// @author Barak Raveh

#ifndef INCLUDED_protocols_flexpep_docking_FlexPepDockingFlags_hh
#define INCLUDED_protocols_flexpep_docking_FlexPepDockingFlags_hh

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/flexpep_docking/FlexPepDockingFlags.fwd.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>
#include <map>
#include <string>



namespace protocols {
namespace flexpep_docking {


class FlexPepDockingFlags
	: public utility::VirtualBase
{
public:

	///////////////////////////////////////////////
	/// @brief
	/// initialize all flags from cmd-line options
	///////////////////////////////////////////////
	FlexPepDockingFlags();

	// TODO: documentation
	FlexPepDockingFlagsOP clone() const {
		return utility::pointer::make_shared< FlexPepDockingFlags >( *this );
	}

	// TODO: more documentation
	// validate chain info (if no info specified, choose first two chains)
	void updateChains(core::pose::Pose const& pose);

	// TODO: more documentation
	// set default chain anchors (peptide c.o.m, and nearest receptor residue)
	void setDefaultAnchors(core::pose::Pose& pose); // TODO: pose should be const, but need to fix RB_geometry 4 this

	// update chane boundaries and anchors -
	// support for old params-file format
	//
	// TODO: also take chain-ids from params file
	void updateChainsAndAnchors_fromParamsFile
	( std::string const& params_file );

	// TODO: documentation
	bool valid_chain_info() const
	{
		return (valid_receptor_chain_ || pep_fold_only) && valid_peptide_chain_ && valid_chain_bounds_;
	}

	bool valid_receptor_chain() const
	{ return valid_receptor_chain_; }

	bool valid_peptide_chain() const
	{ return valid_peptide_chain_; }

	// TODO: documentation
	bool valid_anchors() const
	{
		return
			( !peptide_anchors.empty() &&
			peptide_anchors.size() == peptide_cuts.size() + 1 &&
			(receptor_anchor_pos != 0 || pep_fold_only) );
	}


public:
	/////////////////////
	// Accessor Methods:
	/////////////////////

	bool is_ligand_present( core::pose::Pose const& pose ) const;

	utility::vector1<std::string> const & receptor_chains() const;

	std::string peptide_chain() const;

	void set_receptor_chains(utility::vector1<std::string> const & ch){
		receptor_chains_ = ch; valid_receptor_chain_ = true;
	}

	void set_peptide_chain(std::string const & ch){
		peptide_chain_ = ch; valid_peptide_chain_ = true;
	}

	void set_user_defined_receptor(bool state) {
		user_set_receptor_chain_=state;
	}

	void set_user_defined_peptide(bool state) {
		user_set_peptide_chain_=state;
	}

	bool user_defined_receptor_chain() const
	{
		return user_set_receptor_chain_;
	}

	bool user_defined_peptide_chain() const
	{
		return user_set_peptide_chain_;
	}


	core::Size receptor_first_res() const;

	core::Size receptor_last_res() const;

	core::Size receptor_nres() const;

	core::Size peptide_first_res() const;

	core::Size peptide_last_res() const;

	core::Size peptide_nres() const;

	std::string ref_start_struct() const;

	bool valid_ref_start_struct() const;


private:
	// TODO: change all ints to core::Size, and add access to basic::options::user() to check validity
	// chain info fields
	utility::vector1< std::string > receptor_chains_;
	std::string peptide_chain_;
	core::Size receptor_first_res_;
	core::Size receptor_nres_;
	core::Size peptide_first_res_;
	core::Size peptide_nres_;
	bool valid_receptor_chain_;
	bool valid_peptide_chain_;
	bool valid_chain_bounds_; // refers to fields above
	std::string ref_start_struct_;
	bool valid_ref_start_struct_;
	bool user_set_receptor_chain_;
	bool user_set_peptide_chain_;
public:
	bool pep_fold_only;
	core::Size receptor_anchor_pos; // anchor position within the receptor protein
	std::map< core::Size, core::Size > peptide_cuts; // peptide internal cuts
	std::map< core::Size, core::Size > peptide_anchors; // anchors for peptide fragments
	std::string params_file; // parameters for describing the complex, anchor residues, etc.
	bool lowres_abinitio;
	bool lowres_preoptimize;
	bool min_only;
	bool pep_refine;
	bool random_phi_psi_pert;
	core::Real random_phi_psi_pert_size;
	bool extend;
	bool place_peptide;
	int sample_pc;
	bool slideintocontact;
	bool randomRBstart;
	bool recal_foldtree;
	bool rbMCM;
	core::Real rb_trans_size;
	core::Real rb_rot_size;
	bool torsionsMCM;
	bool peptide_loop_model;
	core::Real smove_angle_range;
	bool design_peptide;
	bool backrub_opt;
	bool boost_fa_atr;
	bool ramp_fa_rep;
	bool ramp_rama;
	int rep_ramp_cycles;
	int mcm_cycles; // TODO: spearate rigid-body and torsion cycles #
	bool min_receptor_bb;
	bool score_only;
	//bool use_cen_score;
	bool ppk_only;
	bool no_prepack1;
	bool no_prepack2;
	core::Real score_filter;
	core::Size hb_filter;
	core::Size hotspot_filter;
	core::Real frag3_weight;
	core::Real frag5_weight;
	core::Real frag9_weight;
	bool pSer2Asp_centroid;
	bool pSer2Glu_centroid;
	bool dumpPDB_abinitio;
	bool dumpPDB_lowres;
	bool dumpPDB_hires;
}; // class FlexPepDockingFlags

} // namespace flexpep_docking
} // namespace protocols

#endif

