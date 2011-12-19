// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/Splice.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_Splice_hh
#define INCLUDED_protocols_protein_interface_design_movers_Splice_hh
#include <protocols/protein_interface_design/movers/Splice.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief designs alanine residues in place of the residue identities at the interface. Retains interface glycines and prolines.
class Splice : public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;
public:
	Splice();
	void apply( Pose & pose );
	virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new Splice ); }
		void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	virtual ~Splice();
	void from_res( core::Size const f ){ from_res_ = f; }
	core::Size from_res() const { return from_res_; }
	void to_res( core::Size const t ){ to_res_ = t; }
	core::Size to_res() const { return to_res_; }
	std::string source_pdb() const { return source_pdb_; }
	void source_pdb( std::string const s ){ source_pdb_ = s; }
	void ccd( bool const c ){ ccd_ = c;}
	bool ccd() const { return ccd_; }
	void scorefxn( core::scoring::ScoreFunctionOP sf );
	core::scoring::ScoreFunctionOP scorefxn() const;
	core::Real rms_cutoff() const{ return rms_cutoff_; }
	void rms_cutoff( core::Real const r ){ rms_cutoff_ = r; }
	void res_move( core::Size const r ){ res_move_ = r; }
	core::Size res_move() const{ return res_move_; }
private:
	core::Size from_res_, to_res_;
	std::string source_pdb_;
	bool ccd_;//dflt true; do ccd?
	core::scoring::ScoreFunctionOP scorefxn_; //dflt score12 with reweighted sheet weight
	core::Real rms_cutoff_; //dflt 99999; after splicing, checks the average displacement of Ca atoms in the source and target segments. Failure leads to mover failure and no output
	core::Size res_move_; //dflt 4; how many residues to allow to move during ccd
};


} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_Splice_HH*/
