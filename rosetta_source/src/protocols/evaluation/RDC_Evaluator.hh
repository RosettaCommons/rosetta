// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @detailed
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange



#ifndef INCLUDED_protocols_evaluation_RDC_Evaluator_hh
#define INCLUDED_protocols_evaluation_RDC_Evaluator_hh


// Unit Headers
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/RDC_Evaluator.fwd.hh>
// Package Headers
// AUTO-REMOVED #include <protocols/evaluation/util.hh>
// AUTO-REMOVED #include <core/scoring/ResidualDipolarCoupling.hh>

// Project Headers
#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

//// C++ headers
#include <list>

#include <core/scoring/ResidualDipolarCoupling.fwd.hh>
#include <core/scoring/rms_util.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace evaluation {

class SelectRDC_Evaluator : public SingleValuePoseEvaluator< core::Real > {
public:
  SelectRDC_Evaluator( core::scoring::ResidueSelection const& selection, std::string tag = "", std::string file ="" );
  SelectRDC_Evaluator( utility::vector1< core::Size> const& selection, std::string tag = "" , std::string file ="");

  //work it out by yourself from missing density == whacky random coords
  SelectRDC_Evaluator( core::pose::PoseCOP, std::string tag = "" );

  //work it out by yourself from missing density == whacky random coords
  SelectRDC_Evaluator( core::pose::Pose const&, std::string tag = "" );

  ///@brief evaluate pose
  virtual core::Real apply( core::pose::Pose& ) const;

private:

  void init_rdcs();

  core::scoring::ResidueSelection selection_;
  std::string tag_;
  core::scoring::ResidualDipolarCouplingOP rdc_data_;
	std::string rdc_file_;
};


}
}

#endif
