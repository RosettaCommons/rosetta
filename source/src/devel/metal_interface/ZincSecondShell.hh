// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/metal_interface/ZincSecondShell.hh
/// @brief
/// @author Bryan Der

#ifndef INCLUDED_devel_metal_interface_ZincSecondShell_HH
#define INCLUDED_devel_metal_interface_ZincSecondShell_HH

#include <devel/metal_interface/ZincSecondShell.fwd.hh>

#include <devel/metal_interface/MetalSiteResidue.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
//#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <sstream>
#include <basic/MetricValue.hh>
#include <set>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/types.hh>


namespace devel {
namespace metal_interface {


/// @details
class ZincSecondShell : public utility::pointer::ReferenceCount {

public:

	typedef core::pose::Pose Pose;

  /// @brief
  ZincSecondShell( core::pose::Pose const & pose, utility::vector1< devel::metal_interface::MetalSiteResidueOP > msr );

  virtual ~ZincSecondShell();


  virtual void register_calculators();

	virtual utility::vector1< core::id::AtomID > get_first_shell_atom_ids();
	virtual utility::vector1< core::id::AtomID > get_second_shell_atom_ids();
	virtual basic::MetricValue< core::id::AtomID_Map< core::Real > > get_atom_sasa();
	virtual basic::MetricValue< core::id::AtomID_Map< core::Size > > get_atom_hbonds();

	virtual void fill_first_shell_atom_ids();
	virtual void fill_second_shell_atom_ids();

	virtual void calculate_hbonds_and_sasa( core::pose::Pose const & pose );
	virtual void report_buried_unsat( utility::vector1< core::id::AtomID > atom_ids );

	virtual core::Size satisfaction_cutoff( std::string atom_type );


private:

	core::pose::Pose pose_;
	std::string pdbname_;

	utility::vector1< devel::metal_interface::MetalSiteResidueOP > msr_;

	core::scoring::hbonds::HBondSet hbond_set_;

	basic::MetricValue< core::id::AtomID_Map< core::Real > > atom_sasa_;
	basic::MetricValue< core::id::AtomID_Map< core::Size > > atom_hbonds_;
	bool hbond_sasa_have_been_calculated_;


	utility::vector1< core::id::AtomID > first_shell_atom_ids_;
	utility::vector1< std::string > first_shell_atom_names_;

	utility::vector1< core::id::AtomID > second_shell_atom_ids_;
	utility::vector1< std::string > second_shell_atom_names_;


	utility::vector1< core::Real > second_shell_atom_hbond_energy_;

	utility::vector1< core::id::AtomID > third_shell_atom_ids_;

};//end ZincSecondShell


}//namespace metal_interface
}//namespace devel

#endif // INCLUDED_devel_metal_interface_ZincSecondShell_HH
