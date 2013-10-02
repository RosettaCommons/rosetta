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
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Residue.hh>
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
#include <ObjexxFCL/format.hh>

// C++ Headers

static basic::Tracer TR("devel.matdes.RestrictToNonzeroSASAOperation" );

namespace devel {
namespace matdes {

core::pack::task::operation::TaskOperationOP
RestrictToNonzeroSASAOperationCreator::create_task_operation() const
{
	return new RestrictToNonzeroSASAOperation;
}


RestrictToNonzeroSASAOperation::RestrictToNonzeroSASAOperation( core::Real probe_radius /* = 2.2 */, core::Size ncomp /* = 1 */ ):
	probe_radius_(probe_radius),
	ncomp_(ncomp)
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
	core::Size res_count = 0;

	core::pose::Pose pose_asu;
	core::pose::symmetry::extract_asymmetric_unit(pose, pose_asu);
	for (core::Size i = 1; i <= ncomp_; i++) {
		std::string select_buried_pos("select buried_pos, resi ");
		core::pose::PoseOP mono = pose_asu.split_by_chain(i); // Extract monomer from each component
  	utility::vector1<Real> sc_sasa = devel::matdes::sidechain_sasa(*mono, probe_radius_);
		//mono.dump_pdb("mono.pdb");
		// If the residue is totally buried, prevent_repacking
		for (core::Size ir = 1; ir <= mono->n_residue(); ir++) {
			if( !mono->residue( ir ).is_protein() ) continue;
			res_count++;
			if (sc_sasa[ir] <= 0.0) {
				select_buried_pos.append(ObjexxFCL::string_of(res_count) + "+");
				TR.Debug << "resi " << res_count << " will not be repacked because it is buried." << std::endl;
				task.nonconst_residue_task(res_count).prevent_repacking();
			}
		}
		TR << "chain " << i << ": " << select_buried_pos << std::endl;
	}
}

void
RestrictToNonzeroSASAOperation::parse_tag( TagPtr tag )
{

  probe_radius_ = tag->getOption<core::Real>("probe_radius", 2.2);
  ncomp_ = tag->getOption<core::Size>("ncomp", 1);

}

void
RestrictToNonzeroSASAOperation::parse_def( utility::lua::LuaObject const & def) {
	probe_radius_ = def["probe_radius"] ? def["probe_radius"].to<core::Real>() : 2.2;
}
} //namespace matdes
} //namespace devel
