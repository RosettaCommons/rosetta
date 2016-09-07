// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/cutoutdomain/CutOutDomain.hh
/// @author Gideon Lapidoth (glapidoth@gmail.com)

#ifndef INCLUDED_devel_cutoutdomain_CutOutDomain_hh
#define INCLUDED_devel_cutoutdomain_CutOutDomain_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector1.hh>


namespace devel {
namespace cutoutdomain {

/// @brief This mover takes a template pdb and cuts the active pose accroding to start and end position relative to the template pose
class CutOutDomain : public protocols::moves::Mover
{
public:
	CutOutDomain();
	~CutOutDomain() override;
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;
	protocols::moves::MoverOP clone() const override {
		return( protocols::moves::MoverOP( new CutOutDomain( *this ) ) );
	}
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new CutOutDomain ); }

	void start_res(core::Size s) {start_res_ = s;}
	void start_end(core::Size e) {end_res_ = e;}
	void suffix(std::string suf) {suffix_ = suf;}
	void source_pdb_name(std::string name) {source_pdb_name_ = name;}
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
private:
	core::Size start_res_;//where the domain will be cut from
	core::Size end_res_;//where the domain will be cut to
	std::string source_pdb_name_;
	std::string suffix_; //the output of this mover will be suffix+source_pdb_name
	core::Size
	find_nearest_res( core::pose::Pose const & source, core::pose::Pose const & target, core::Size const res, core::Size const chain/*=0*/ );
	core::Size delta_n_ter_;//one can cut an additioanl delta aa from designated cut site. dflt=0
	core::Size delta_c_ter_;
};


} // CutOutDomain
} // devel


#endif
