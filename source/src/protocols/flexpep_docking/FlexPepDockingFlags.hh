// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file   FlexPepDockingFlags.hh
///
/// @brief flags structure for FlexPepDocking protocols
/// @date January 1, 2009
/// @author Barak Raveh

#ifndef INCLUDED_protocols_flexpep_docking_FlexPepDockingFlags_hh
#define INCLUDED_protocols_flexpep_docking_FlexPepDockingFlags_hh

#include <core/pose/Pose.fwd.hh>
#include <protocols/flexpep_docking/FlexPepDockingFlags.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <map>
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace flexpep_docking {


class FlexPepDockingFlags
	: public utility::pointer::ReferenceCount
{
public:

	///////////////////////////////////////////////
	/// @brief
	/// initialize all flags from cmd-line options
	///////////////////////////////////////////////
	FlexPepDockingFlags();

	// TODO: documentation
	FlexPepDockingFlagsOP clone() const {
		return FlexPepDockingFlagsOP( new FlexPepDockingFlags( *this ) );
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
			(peptide_anchors.size() > 0 &&
			peptide_anchors.size() == peptide_cuts.size() + 1 &&
			(receptor_anchor_pos != -1 || pep_fold_only) );
	}


public:
	/////////////////////
	// Accessor Methods:
	/////////////////////

	bool is_ligand_present( core::pose::Pose const& pose ) const;

	std::string receptor_chain() const;

	char peptide_chain() const;

	void set_receptor_chain(std::string ch){
		receptor_chain_ = ch; valid_receptor_chain_ = true;
	}

	void set_peptide_chain(char ch){
		peptide_chain_ = ch; valid_peptide_chain_ = true;
	}

	int receptor_first_res() const;

	int receptor_last_res() const;

	int receptor_nres() const;

	int peptide_first_res() const;

	int peptide_last_res() const;

	int peptide_nres() const;

	std::string ref_start_struct() const;

	bool valid_ref_start_struct() const;


private:
	// TODO: change all ints to Size, and add access to basic::options::user() to check validity
	// chain info fields
	std::string receptor_chain_;
	char peptide_chain_;
	int receptor_first_res_;
	int receptor_nres_;
	int peptide_first_res_;
	int peptide_nres_;
	bool valid_receptor_chain_;
	bool valid_peptide_chain_;
	bool valid_chain_bounds_; // refers to fields above
	std::string ref_start_struct_;
	bool valid_ref_start_struct_;
public:
	bool pep_fold_only;
	int receptor_anchor_pos; // anchor position within the receptor protein
	std::map<int,int> peptide_cuts; // peptide internal cuts
	std::map<int,int> peptide_anchors; // anchors for peptide fragments
	std::string params_file; // parameters for describing the complex, anchor residues, etc.
	bool lowres_abinitio;
	bool lowres_preoptimize;
	bool min_only;
	bool pep_refine;
	bool random_phi_psi_pert;
	double random_phi_psi_pert_size;
	bool extend;
	bool place_peptide;
	int sample_pc;
	bool slideintocontact;
	bool randomRBstart;
	bool recal_foldtree;
	bool rbMCM;
	double rb_trans_size;
	double rb_rot_size;
	bool torsionsMCM;
	bool peptide_loop_model;
	double smove_angle_range;
	bool design_peptide;
	bool backrub_opt;
	bool boost_fa_atr;
	bool ramp_fa_rep;
	bool ramp_rama;
	int rep_ramp_cycles;
	int mcm_cycles; // TODO: spearate rigid-body and torsion cycles #
	bool min_receptor_bb;
	bool score_only;
	bool use_cen_score;
	bool ppk_only;
	bool no_prepack1;
	bool no_prepack2;
	double score_filter;
	int    hb_filter;
	int    hotspot_filter;
	double frag3_weight;
	double frag5_weight;
	double frag9_weight;
	bool pSer2Asp_centroid;
	bool pSer2Glu_centroid;
	bool dumpPDB_abinitio;
	bool dumpPDB_lowres;
	bool dumpPDB_hires;
}; // class FlexPepDockingFlags

} // namespace flexpep_docking
} // namespace protocols

#endif

