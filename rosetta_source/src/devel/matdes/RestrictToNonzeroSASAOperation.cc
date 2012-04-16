// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/matdes/RestrictToNonzeroSASAOperation.cc
/// @brief  Restrict design to only residues at inter-building block interfaces
/// @author Neil King (neilking@uw.edu)

// Unit Headers
#include <devel/matdes/RestrictToNonzeroSASAOperation.hh>
#include <devel/matdes/RestrictToNonzeroSASAOperationCreator.hh>

// Project Headers
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <devel/matdes/util.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// C++ Headers

static basic::Tracer TR("devel.matdes.RestrictToNonzeroSASAOperation" );

namespace devel {
namespace matdes {

core::pack::task::operation::TaskOperationOP
RestrictToNonzeroSASAOperationCreator::create_task_operation() const
{
	return new RestrictToNonzeroSASAOperation;
}


RestrictToNonzeroSASAOperation::RestrictToNonzeroSASAOperation( core::Real probe_radius /* = 2.2 */):
	probe_radius_(probe_radius)
{}

RestrictToNonzeroSASAOperation::~RestrictToNonzeroSASAOperation() {}

core::pack::task::operation::TaskOperationOP RestrictToNonzeroSASAOperation::clone() const
{
	return new RestrictToNonzeroSASAOperation( *this );
}

void
RestrictToNonzeroSASAOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{

  // Get the SASA for each residue in the monomeric state
	core::conformation::symmetry::SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
  core::pose::Pose mono(pose, 1, sym_info->get_nres_subunit());
  utility::vector1<Real> sc_sasa = devel::matdes::sidechain_sasa(mono, probe_radius_);

	// If the residue is totally buried, prevent_repacking
	for (core::Size ir = 1; ir <= sym_info->get_nres_subunit(); ir++) {
		if (sc_sasa[ir] <= 0.0) {
			//TR << "resi " << ir << " will not be repacked because it is buried." << std::endl;
			task.nonconst_residue_task(ir).prevent_repacking();
		}
	}

}

void
RestrictToNonzeroSASAOperation::parse_tag( TagPtr tag )
{

  probe_radius_ = tag->getOption<core::Real>("probe_radius", 2.2);
	
}

} //namespace matdes
} //namespace devel
