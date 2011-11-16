// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/RotamerConstraint.cc
///
/// @brief
/// @author Ian W. Davis


#include <core/pack/dunbrack/RotamerConstraint.hh>

#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <basic/Tracer.hh>


// option key includes

#include <basic/options/keys/packing.OptionKeys.gen.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/types.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
// AUTO-REMOVED #include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
// AUTO-REMOVED #include <core/id/SequenceMapping.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/io/pdb/file_data.fwd.hh>
#include <core/io/pdb/file_data.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerConstraint.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Remarks.fwd.hh>
#include <core/pose/Remarks.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/FuncFactory.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/XYZ_Func.fwd.hh>
#include <core/scoring/constraints/XYZ_Func.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/fixedsizearray1.fwd.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.fwd.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>



namespace core {
namespace pack {
namespace dunbrack {

using namespace scoring;
using namespace scoring::constraints;

static basic::Tracer TR("core.pack.dunbrack.RotamerConstraint");


void load_unboundrot(pose::Pose & pose)
{
	using namespace basic::options;
	using namespace core::pose;
	if( !option[ OptionKeys::packing::unboundrot ].active() ) return;

	static utility::vector1< PoseCOP > unboundrot_poses;
	if( unboundrot_poses.empty() ) {
		for(Size i = 1; i <= option[ OptionKeys::packing::unboundrot ]().size(); ++i) {
			std::string filename = option[ OptionKeys::packing::unboundrot ]()[i].name();
			TR << "Adding 'unbound' rotamers from " << filename << std::endl;
			PoseOP pose = new Pose();
			//core::import_pose::pose_from_pdb( *pose, filename );
			core::io::pdb::build_pose_from_pdb_as_is( *pose, filename );
			unboundrot_poses.push_back( pose );
		}
	}

	if( unboundrot_poses.empty() ) return; // guaranteed at least one pose now

	using namespace std;
	for(Size rsd_num = 1; rsd_num <= pose.total_residue(); ++rsd_num) {
		// Each constraint can contain only one ResidueType, so we have to sort them out here.
		// We should get any scoring overlap since ResidueTypes are mutually exclusive.
		map< string, RotamerConstraintOP > by_res_type;
		for(Size pose_num = 1; pose_num <= unboundrot_poses.size(); ++pose_num) {
			if( rsd_num > unboundrot_poses[pose_num]->total_residue() ) continue;
			conformation::Residue const & rsd = unboundrot_poses[pose_num]->residue(rsd_num);
			if( !rsd.is_protein() ) {
				// Can't determine rotamer number for anything but protein.
				TR << "Can't use " << rsd.type().name() << " " << rsd_num << " for residue constraint -- protein only." << std::endl;
				continue;
			}
			if( core::pack::dunbrack::RotamerLibrary::get_instance().get_rsd_library( rsd.type() )() == 0 ) continue; // no point in creating constraint
			if( by_res_type.find( rsd.type().name() ) == by_res_type.end() ) { // first one, create constraint
				TR.Debug << "Creating rotamer constraint for " << rsd.type().name() << " at " << rsd_num << std::endl;
				RotamerConstraintOP constraint = new RotamerConstraint( *unboundrot_poses[pose_num], rsd_num );
				pose.add_constraint( constraint );
				by_res_type[ rsd.type().name() ] = constraint;
			} else { // subsequent one, just add residue
				RotamerConstraintOP constraint = by_res_type[ rsd.type().name() ];
				constraint->add_residue( unboundrot_poses[pose_num]->residue(rsd_num) );
			}
		}
	}

}


RotamerConstraint::RotamerConstraint(
	pose::Pose const & pose,
	Size seqpos
):
	Constraint( core::scoring::fa_dun ), // most like a Dunbrack rotamer energy, and should share the same weight
	seqpos_( seqpos ),
	rsd_type_name_( pose.residue_type(seqpos).name() ),
	atom_ids_(),
	rotlib_( core::pack::dunbrack::RotamerLibrary::get_instance().get_rsd_library( pose.residue_type(seqpos) ) ), // may be NULL
	favored_rotamers_(),
	favored_rotamer_numbers_()
{
	// Depends on ~ all heavy atoms (all chis + phi and psi)
	// Could cause problems if this position mutates to something with fewer atoms?
	conformation::Residue const & rsd( pose.residue(seqpos_) );
	for(Size i = 1, i_end = rsd.nheavyatoms(); i <= i_end; ++i) {
		atom_ids_.push_back(AtomID( i, seqpos_ )); // atom, rsd
	}

	// Depends on neighbors' backbone, too, due to phi,psi dependence
	// This could cause problems if neighbor residues went away / mutated radically?
	for(Size conn = 1, conn_end = rsd.n_residue_connections(); conn <= conn_end; ++conn) {
		Size const nbr_seqpos( rsd.connected_residue_at_resconn(conn) );
		conformation::Residue const & nbr( pose.residue(nbr_seqpos) );
		for(Size i = 1, i_end = nbr.type().last_backbone_atom(); i <= i_end; ++i) {
			atom_ids_.push_back(AtomID( i, nbr_seqpos )); // atom, rsd
		}
	}

	add_residue( rsd );
}


RotamerConstraint::~RotamerConstraint() {}


void
RotamerConstraint::add_residue( conformation::Residue const & rsd )
{
	assert( rsd.type().name() == rsd_type_name_ );
	favored_rotamers_.push_back( rsd.chi() );
	pack::dunbrack::RotVector rot;
	pack::dunbrack::rotamer_from_chi( rsd, rot );
	favored_rotamer_numbers_.push_back( rot );
}


Size
RotamerConstraint::natoms() const
{
	return atom_ids_.size();
}


id::AtomID const &
RotamerConstraint::atom( Size const index ) const
{
	return atom_ids_[index];
}


// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
RotamerConstraint::score(
	XYZ_Func const & xyz_func,
	EnergyMap const & weights,
	EnergyMap & emap
) const
{
	if( rotlib_() == NULL ) return;
	if( weights[ this->score_type() ] == 0 ) return; // what's the point?

	conformation::Residue const & rsd( xyz_func.residue(seqpos_) );
	if( rsd.type().name() != rsd_type_name_ ) return; // residue types must match

	pack::dunbrack::RotVector rot;
	pack::dunbrack::rotamer_from_chi( rsd, rot );
	for(Size i = 1, i_end = favored_rotamer_numbers_.size(); i <= i_end; ++i) {
		if( rot == favored_rotamer_numbers_[i] ) {
			pack::dunbrack::RotamerLibraryScratchSpace scratch;
			Real const best_rotE = rotlib_->best_rotamer_energy(rsd, false /* => global min */, scratch);
			Real const this_rotE = rotlib_->best_rotamer_energy(rsd, true /* => local min */, scratch);
			assert( best_rotE <= this_rotE );
			//std::cout << "rotamer constraint active for " << seqpos_ << " thisE = " << this_rotE << " bestE = " << best_rotE << " dE = " << ( best_rotE - this_rotE ) << std::endl;
			emap[ this->score_type() ] +=  ( best_rotE - this_rotE );
			return; // quit once we find a match
		}
	}
	// no match found, don't adjust score for this rotamer
}


void
RotamerConstraint::fill_f1_f2(
	AtomID const & ,//atom,
	XYZ_Func const &,
	Vector & ,//F1,
	Vector & ,//F2,
	EnergyMap const & //weights
) const
{
	// Do nothing.
	// Derivative of this restraint is effectively zero (because it's constant
	// within each rotamer well), so we just "add zero" to F1 and F2.
}


void RotamerConstraint::show( std::ostream & out ) const
{
	out << type() << ": " << seqpos_ << " " << rsd_type_name_ << " " << favored_rotamer_numbers_.size();
	for(Size i = 1, i_end = favored_rotamer_numbers_.size(); i <= i_end; ++i) {
		out << ",";
		for(Size j = 1, j_end = favored_rotamer_numbers_[i].size(); j <= j_end; ++j) {
			out << " " << favored_rotamer_numbers_[i][j];
		}
	}
	out << std::endl;
}


} // namespace constraints
} // namespace scoring
} // namespace core
