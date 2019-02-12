// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/tcr/TCRmodel.hh
/// @brief Class for the tcr modeling protocol
/// @author Ragul Gowthaman (ragul@umd.edu)

#ifndef INCLUDED_protocols_tcr_TCRmodel_hh
#define INCLUDED_protocols_tcr_TCRmodel_hh

#include <protocols/tcr/TCRmodel.fwd.hh>
#include <protocols/tcr/TCRseqInfo.fwd.hh>
#include <protocols/tcr/TCRseqInfo.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <iostream>
#include <sstream>
#include <fstream>
#include <protocols/antibody/grafting/antibody_sequence.hh>
#include <set>
///////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace tcr {

/// @brief class for modeling TCR structure with optional refinement/loop modeling
/// @details TCRmodel class uses TCRseqInfo for parsing and numbering the sequence
/// TCrmodel takes in tcr a & b chain sequence or TCRseqInfo object directly and
/// do template identification, grafting and assembling of TCR structure
/// Calls various utility functions from util files in protocols/tcr
/// Modeling protool uses either germline segments or framework segments
/// germline segments (germline + cdr3 sequence)
/// framework segments (framework + cdr1+ cdr2 + cdr3 sequence)
/// currently extended cdr2 segment is used (cdr2 extended upto hv4 loop)
/// see also: TCRmodel documentation and publication
/// see also: tcr/TCRseqInfo
class TCRmodel : public utility::pointer::ReferenceCount {

public:

	TCRmodel( std::string const & aseq, std::string const &bseq );
	TCRmodel( TCRseqInfoOP tcrinfo );
	//TCRmodel();

	/// @brief Clone the pointer
	TCRmodelOP
	clone() const;

	//Default destructor
	~TCRmodel() override;

	/// @brief tmpltinfo holds the template information for an individual tcr segment
	/// @details template pdb id(tid), pdb structure(tpdb), sequence(tseq), alignment score(tscore), pose(tpiece) and template databse used(tdb).
	struct tmpltinfo {
		std::string tid;
		std::string tpdb;
		std::string tseq;
		core::Real tscore;
		core::pose::PoseOP tpiece;
		utility::vector1<std::string> tdb;
	};

	/// @brief tcrtmplts holds the tmpltinfo values for different tcr segments
	/// @details orientation(ori), framework(fr), germline(gm), CDR1(cdr1), extended CDR2(cdr2hv4), CDR3(cdr3)
	struct tcrtmplts {
		tmpltinfo ori;
		tmpltinfo fr;
		tmpltinfo gm;
		tmpltinfo cdr1;
		tmpltinfo cdr2hv4;
		tmpltinfo cdr3;
	};

	///Setters
	/// @brief use_gma_templates uses germline segments (germline + cdr3 sequence) for modeling alpha chain
	/// @details if false, (framework + cdr1+ cdr2 + cdr3 sequence) is used for modeling
	void set_use_gma_templates(bool gma) {
		use_gma_templates_ = gma;
	}
	/// @brief use_gma_templates uses germline segments (germline + cdr3 sequence) for modeling beta chain
	/// @details if false, (framework + cdr1+ cdr2 + cdr3 sequence) is used for modeling
	void set_use_gmb_templates(bool gmb) {
		use_gmb_templates_ = gmb;
	}

	/// @brief returns the template information (tmpltinfo) for tcr alpha chain segments
	tcrtmplts atmplt() const { return atmplt_;};

	/// @brief returns the template information (tmpltinfo) for tcr alpha chain segments
	tcrtmplts btmplt() const { return btmplt_;};

	/// @brief returns true if germline alpha chain sequence is used for modeling
	bool use_gma_templates() const { return use_gma_templates_; }

	/// @brief returns true if germline beta chain sequence is used for modeling
	bool use_gmb_templates() const { return use_gmb_templates_; }

	core::pose::PoseOP tcr_model() const { return tcr_model_; }
	core::pose::PoseOP tcr_graft_model() const { return tcr_graft_model_; }
	core::pose::PoseOP tcr_loop_model() const { return tcr_loop_model_; }

private:

	/// @brief check for user provided options and assign values.
	/// @details assign values for template files, template cutoffs, model refinement options etc..
	void init_from_options();

	/// @brief set default values for database path and other members
	void set_default();

	/// @brief Find templates for tcr segments
	/// @details Checks for user provided template pdb id or structure file before searching for templates from db
	/// @details Calls individual util functions for identification and grafting of templates for TCR segments
	void setup_templates();

	/// @brief make_model parent function for assembling the model from templates
	/// @details calls other functions to build graft model, optional refinement and loop modeling
	void make_model();

	/// @brief build the crude model from the grafted templates
	void build_graft_model();

private:

	//@brief TCRseqInfo holds the parsed tcr segments and numbering of CDR's
	tcr::TCRseqInfoOP tsinfo;

	//@brief template informations for alpha chain segments
	tcrtmplts atmplt_;
	//@brief template informations for beta chain segments
	tcrtmplts btmplt_;

	/// @brief No. of overhang residues used for grafting; cterminal & nterminal
	core::Size nter_overhang_;
	core::Size cter_overhang_;

	/// @brief path of the template database seq dir, structure dir
	std::string tseq_db_;
	std::string tpdb_db_;

	/// @brief use user provide templates; no template identification from template db
	bool use_user_templates_;

	/// @brief use germline templates instead of framework templates for a/b tcr chain
	bool use_gma_templates_;
	bool use_gmb_templates_;

	/// @brief include antibody templates during template search
	bool include_ab_templates_;
	/// @brief templated db path for antibody templates
	std::string ab_db_path_;

	/// @brief list of PDB id's to blacklist as templates
	/// @details list provided by user or created by blast/identity cutoff
	std::list< std::set<std::string> > ignore_lists_;

	/// @brief cut-off values
	/// @details cutoff values used to blacklist pdb templates
	core::Real blastp_identity_cutoff_;
	core::Real template_identity_cutoff_;

	/// @brief member poseop's to store modeled poses
	core::pose::PoseOP tcr_graft_model_;
	core::pose::PoseOP tcr_loop_model_;
	core::pose::PoseOP tcr_model_;

	/// @brief dump all templates
	bool dump_templates_;

	/// @skip modeling part, used for testing
	bool skip_modeling_;

	/// @brief minimize or relax option for model structure
	bool minimize_model_;
	bool relax_model_;

	core::scoring::ScoreFunctionOP scorefxn_;

};


} //namespace tcr
} //namespace protocols


#endif //INCLUDED_protocols_loops_TCRmodel_HH
