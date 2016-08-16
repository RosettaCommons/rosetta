// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/vip/VIP_Report.hh
/// @brief class for reports from vip mover implementation


#ifndef INCLUDED_protocols_vip_VIP_Report_HH
#define INCLUDED_protocols_vip_VIP_Report_HH

#include <core/pose/Pose.hh>

namespace protocols {
namespace vip {

class VIP_Report
{

public:
	VIP_Report();
	virtual ~VIP_Report();

	//  void define_report_file();
	//  void close_report_file();
	void get_GOE_repack_report( core::pose::Pose & goe_native, utility::vector1<core::Real> & goe_repack_e, utility::vector1<core::conformation::ResidueOP> & goe_repack_res, utility::vector1<core::Size> & goe_repack_pos, core::Size it, bool use_stored, core::Real stored_e );
	void get_GOE_relaxed_report( core::pose::Pose & goe_native, utility::vector1<core::Real> & goe_repack_e, utility::vector1<core::conformation::ResidueOP> & goe_repack_res, utility::vector1<core::Size> & goe_repack_pos, core::Size it, bool use_stored, core::Real stored_e );
	void get_GOE_packstat_report(   core::pose::Pose & goe_native, utility::vector1<core::pose::Pose> & goe_relax );

	//  friend class VIP_Mover;
};

}}
#endif

