// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/metal_interface/LigandBurial.hh
/// @brief  Takes two protein and two match poses, grafting match onto protein, then combines the two grafted poses by overlaying the zinc atoms.
/// @author Bryan Der

#ifndef INCLUDED_devel_metal_interface_LigandBurial_HH
#define INCLUDED_devel_metal_interface_LigandBurial_HH

#include <devel/metal_interface/LigandBurial.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
//#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <sstream>
#include <basic/MetricValue.hh>
#include <set>

namespace devel {
namespace metal_interface {


/// @details
class LigandBurial : public utility::pointer::ReferenceCount {

public:

	typedef core::pose::Pose Pose;
	//typedef std::set< core::Size > SetSize;

  /// @brief
  LigandBurial( Pose const & pose, std::string ligand_3_letter_code );

  ~LigandBurial() override;

  virtual core::Size find_ligand();
  virtual void register_calculators();

	virtual basic::MetricValue< std::set<core::Size> > get_ligand_neighbors();
	virtual core::Real get_ligand_sasa();

  virtual void calculate_ligand_neighbors();
	virtual void calculate_ligand_sasa();

private:
	std::stringstream calcname_;
	core::Size ligand_resnum_;
	core::pose::Pose pose_; //does anyone ever do this?  Should I use PoseOP?
	std::string ligand_3_letter_code_;

	basic::MetricValue< std::set<core::Size> > ligand_neighbors_;
	core::Real ligand_sasa_;

};//end LigandBurial


}//namespace metal_interface
}//namespace devel

#endif // INCLUDED_devel_metal_interface_LigandBurial_HH
