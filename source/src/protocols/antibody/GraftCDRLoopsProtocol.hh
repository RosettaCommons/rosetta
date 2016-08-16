// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/GraftCDRLoopsProtocol.hh
/// @brief Build a homology model of an antibody
/// @details
///
///
/// @author Jianqing Xu ( xubest@gmail.com )


#ifndef INCLUDED_protocols_antibody_GraftCDRLoopsProtocol_hh
#define INCLUDED_protocols_antibody_GraftCDRLoopsProtocol_hh


#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/Ab_TemplateInfo.fwd.hh>
#include <protocols/antibody/GraftCDRLoopsProtocol.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/antibody/CDRsMinPackMin.fwd.hh>


namespace protocols {
namespace antibody {

class GraftCDRLoopsProtocol: public moves::Mover {
public:
	typedef std::map < std::string, bool > GraftMap;

	// default constructor
	GraftCDRLoopsProtocol();

	// default destructor
	~GraftCDRLoopsProtocol();

	virtual protocols::moves::MoverOP clone() const;

	/// @brief Assigns default values to primitive members
	void set_default();


	virtual void apply( core::pose::Pose & pose );

	// simple inline setters
	void set_graft_l1( bool graft_l1 ) {
		graft_l1_ = graft_l1;
	}
	void set_graft_l2( bool graft_l2 ) {
		graft_l2_ = graft_l2;
	}
	void set_graft_l3( bool graft_l3 ) {
		graft_l3_ = graft_l3;
	}
	void set_graft_h1( bool graft_h1 ) {
		graft_h1_ = graft_h1;
	}
	void set_graft_h2( bool graft_h2 ) {
		graft_h2_ = graft_h2;
	}
	void set_graft_h3( bool graft_h3 ) {
		graft_h3_ = graft_h3;
	}
	void set_h3_stem_graft(bool h3_stem_graft) {
		h3_no_stem_graft_=h3_stem_graft;
	}
	void set_packonly_after_graft (bool setting) {
		packonly_after_graft_ = setting;
	}
	void set_stem_optimize (bool setting) {
		stem_optimize_ = setting;
	}
	void set_camelid( bool camelid ) {
		camelid_ = camelid;
	}
	void set_camelid_constraints( bool camelid_constraints ) {
		camelid_constraints_ = camelid_constraints;
	}
	void set_benchmark( bool benchmark ) {
		benchmark_ = benchmark;
	}
	void set_cst_weight(core::Real cst_weight) {
		cst_weight_ = cst_weight;
	}
	void set_sc_min (bool scmin) {
		sc_min_ = scmin ;
	}
	void set_rt_min (bool rtmin) {
		rt_min_ = rtmin ;
	}
	virtual std::string get_name() const;


	/// @brief Associates relevant options with the AntibodyModeler class
	static void register_options();

	void display_constraint_residues( core::pose::Pose & pose );

	void show( std::ostream & out=std::cout ) const;
	friend std::ostream & operator<<(std::ostream& out, const GraftCDRLoopsProtocol & ab_m_2 );


private:
	void finalize_setup( core::pose::Pose & pose );
	void init();
	void setup_objects();
	void init_from_options();

private:
	bool graft_l1_, graft_l2_, graft_l3_;
	bool graft_h1_, graft_h2_, graft_h3_;
	bool h3_no_stem_graft_;
	bool packonly_after_graft_;
	bool camelid_;
	bool camelid_constraints_;
	bool sc_min_;
	bool rt_min_;
	bool benchmark_;
	bool stem_optimize_;

	core::Real cst_weight_;
	core::scoring::ScoreFunctionOP scorefxn_pack_;

	AntibodyInfoOP ab_info_;
	Ab_TemplateInfoOP ab_t_info_ ;

	core::pack::task::TaskFactoryOP tf_;

	protocols::moves::SequenceMoverOP graft_sequence_ ;
	protocols::moves::SequenceMoverOP optimize_sequence_ ;
	CDRsMinPackMinOP cdrs_min_pack_min_;


}; // class GraftCDRLoopsProtocol


} // namespace antibody
} // namespace protocols

#endif
