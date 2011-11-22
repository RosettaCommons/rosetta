// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file  core/pose/symmetry/util.hh
/// @brief utility functions for handling of symmetric conformations
/// @author Ingemar Andre

// Unit headers
#include <core/pack/make_symmetric_task.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmData.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/VirtualCoordinate.hh>
// AUTO-REMOVED #include <core/scoring/symmetry/SymmetricEnergies.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <core/pack/task/PackerTask.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>

// Utility functions
// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/keys/symmetry.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/fold_and_dock.OptionKeys.gen.hh>
#include <core/id/AtomID.hh>
// AUTO-REMOVED #include <numeric/random/random.hh>
// AUTO-REMOVED #include <numeric/xyzTriple.hh>

// Package Headers
// AUTO-REMOVED #include <core/kinematics/Edge.hh>
#include <core/pose/symmetry/util.hh>

// AUTO-REMOVED #include <core/pose/PDBInfo.hh>

// AUTO-REMOVED #include <basic/Tracer.hh>

// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
// AUTO-REMOVED #include <numeric/xyzVector.io.hh>

#include <utility/vector1.hh>



namespace core {
namespace pack {

void
make_symmetric_PackerTask(
  pose::Pose const & pose,
  pack::task::PackerTaskOP task
)
{
	using namespace conformation::symmetry;
	using namespace pose::symmetry;

	assert( is_symmetric( pose ) );

	SymmetricConformation const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

  for ( Size i = 1; i <= pose.total_residue(); ++i ) {
  if ( !symm_info->chi_is_independent(i) ) {
      task->nonconst_residue_task( i ).prevent_repacking();
    }
  }
}

} // pack
} // core
