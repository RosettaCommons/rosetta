// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/SequnceCouplingConstraint.cc
/// @brief  This is a constraint that refers to a core::sequence::SequnceCopuling? in order to influence the scoring of amino acid types based on multiple sequence alignments (i.e. for biasing amino acid choices during design). A note about the SequenceProfile::read_from_checkpoint function that is used to read in scores for amino acid types: the first line of the file will be ignored.
/// @author HetuKamisetty

#include <protocols/constraints_additional/SequenceCouplingConstraint.hh>

#include <core/conformation/Residue.hh>
#include <core/sequence/SequenceCoupling.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreType.hh>

#include <basic/Tracer.hh>

#include <utility/file/file_sys_util.hh> // file_exists, create_directory

// AUTO-REMOVED #include <basic/options/option.hh>
// option key includes
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>

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
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
// AUTO-REMOVED #include <core/id/SequenceMapping.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
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
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceCoupling.fwd.hh>
#include <core/sequence/SequenceProfile.fwd.hh>
#include <core/sequence/SequenceProfile.hh>
#include <protocols/constraints_additional/SequenceCouplingConstraint.fwd.hh>
// AUTO-REMOVED #include <utility/Bound.fwd.hh>
// AUTO-REMOVED #include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
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
// AUTO-REMOVED #include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/AnyOption.hh>
// AUTO-REMOVED #include <utility/options/AnyVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
// AUTO-REMOVED #include <utility/options/Option.hh>
// AUTO-REMOVED #include <utility/options/OptionCollection.fwd.hh>
// AUTO-REMOVED #include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/RealVectorOption.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption_T_.fwd.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/StringVectorOption.hh>
// AUTO-REMOVED #include <utility/options/VariantOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/VariantOption.hh>
// AUTO-REMOVED #include <utility/options/VectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/VectorOption.hh>
// AUTO-REMOVED #include <utility/options/VectorOption_T_.fwd.hh>
// AUTO-REMOVED #include <utility/options/VectorOption_T_.hh>
// AUTO-REMOVED #include <utility/options/mpi_stderr.hh>
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
// AUTO-REMOVED #include <ObjexxFCL/TypeTraits.hh>
// AUTO-REMOVED #include <ObjexxFCL/char.functions.hh>
// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
// AUTO-REMOVED #include <list>
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

//#include <core/id/SequenceMapping.hh>


namespace protocols {
namespace constraints_additional {

using namespace core;
using namespace chemical;
using namespace conformation;
using namespace basic::options;
using namespace scoring;
using namespace constraints;
using namespace sequence;

using basic::t_warning;
using basic::t_info;
using basic::t_debug;
using basic::t_trace;
static basic::Tracer TR("protocols.constraints_additional.SequenceCouplingConstraint");

SequenceCouplingConstraint::SequenceCouplingConstraint()
	: Constraint( res_type_constraint ),
		seqpos1_(0),
		seqpos2_(0),
		sequence_coupling_(NULL)
{}

SequenceCouplingConstraint::SequenceCouplingConstraint(
	Pose const & pose,
	Size seqpos1,
	Size seqpos2,
	SequenceCouplingOP coupling/* = NULL */
):
	Constraint( res_type_constraint ),
	seqpos1_(seqpos1),
	seqpos2_(seqpos2),
	sequence_coupling_( coupling)
{
	//Residue const & rsd( pose.residue(seqpos1_) );
	atom_ids_.push_back(AtomID(1,seqpos1_));
	//Residue const & rsd( pose.residue(seqpos2_) );
	atom_ids_.push_back(AtomID(1,seqpos2_));
	/*for( Size i(1), i_end( rsd.nheavyatoms() ); i <= i_end; ++i ) {
		atom_ids_.push_back( AtomID( i, seqpos_ ) );
	}
	*/
}

SequenceCouplingConstraint::SequenceCouplingConstraint(
	Size seqpos1,
	Size seqpos2,
	utility::vector1< AtomID > const & atoms_in,
	SequenceCouplingOP sequence_coupling /* = NULL */
):
	Constraint( res_type_constraint ),
	seqpos1_( seqpos1 ),
	seqpos2_( seqpos2 ),
	sequence_coupling_( sequence_coupling )
{
	for( utility::vector1< AtomID >::const_iterator at_it( atoms_in.begin() ), end( atoms_in.end() ); at_it != end; ++at_it ) {
		atom_ids_.push_back( *at_it );
	}
}

SequenceCouplingConstraint::~SequenceCouplingConstraint() {}

ConstraintOP
SequenceCouplingConstraint::clone() const {
	return new SequenceCouplingConstraint( *this );
}

///@details one line definition "SequenceProfile resindex profilefilename" (profilefilename can also be set to "none" in the constraints file, and specified by -in::file::pssm)
void
SequenceCouplingConstraint::read_def(
	std::istream & is,
	Pose const & pose,
	FuncFactory const &
)
{
	Size residue_index1(0);
	Size residue_index2(0);
	std::string coupling_filename;

//	note: is >> "SequenceProfile" has already occured
	is >> residue_index1 >> residue_index2 >> coupling_filename;

	TR(t_debug) << "reading: " << residue_index1 << " " << residue_index2 << " " << coupling_filename << std::endl;
	if ( residue_index1 < 1 || residue_index1 > pose.total_residue() ) {
		std::cerr << "no such residue index " << residue_index1 << " in pose!)" << std::endl;
		utility_exit();
	}
	if ( residue_index2 < 1 || residue_index2 > pose.total_residue() ) {
		std::cerr << "no such residue index " << residue_index2 << " in pose!)" << std::endl;
		utility_exit();
	}

	seqpos1_ = residue_index1;
	seqpos2_ = residue_index2;

	atom_ids_.push_back(AtomID(1,seqpos1_));
	atom_ids_.push_back(AtomID(1,seqpos2_));

	// figure out sequence profile filename
	using namespace utility::file;
	// if specified, verify file exists
	if ( coupling_filename!= "none" ) {
		if ( ! file_exists( coupling_filename ) ) {
			utility_exit_with_message( "no such file " + coupling_filename );
		}
	// if filename not specified, load from commandline option -pssm only if sequence_coupling_ is NULL
	} else {
			utility_exit_with_message("\"none\" is not a valid value for -pssm in this context!");
	}

	// if filename is not "none" by this point, read it even if sequence_coupling_ is not currently NULL
	if ( coupling_filename != "none" ) {
		sequence_coupling_= new SequenceCoupling;
		sequence_coupling_->read_from_file( FileName(coupling_filename) );
	}

	// if sequence_coupling_ is still NULL by this point, it is assumed that the user intended so

} // read_def

void
SequenceCouplingConstraint::show_def( std::ostream & os, Pose const & ) const
{
	show( os );
}

void
SequenceCouplingConstraint::show( std::ostream & os ) const {
	os << "SequenceCouplingConstraint between seqpos " << seqpos1_ << " " << seqpos2_ << ": ";
	if ( ! sequence_coupling_ ) os << "(uninitialized sequence profile)";
//	else {
//		typedef utility::vector1<Real> RealVec;
//		RealVec const & aa_scores( sequence_profile_->prof_row( seqpos_ ) );
//		runtime_assert( aa_scores.size() >= num_canonical_aas );
//		for ( Size aa(1); aa <= num_canonical_aas; ++aa ) {
//			os << aa_scores[aa] << " ";
//		}
//	}
	os << '\n';
}

void
SequenceCouplingConstraint::set_sequence_coupling( SequenceCouplingOP profile )
{
	sequence_coupling_ = profile;
}

SequenceCouplingOP
SequenceCouplingConstraint::sequence_coupling() { return sequence_coupling_; }

SequenceCouplingCOP
SequenceCouplingConstraint::sequence_coupling() const { return sequence_coupling_; }

Size
SequenceCouplingConstraint::natoms() const
{
	return atom_ids_.size();
}

id::AtomID const &
SequenceCouplingConstraint::atom( Size const index ) const
{
	return atom_ids_[index];
}

utility::vector1< id::AtomID > const &
SequenceCouplingConstraint::atom_ids() const { return atom_ids_; }

/*
 * hk: how does one fix sequencecoupling on a remap_resid call?
 *
ConstraintOP
SequenceCouplingConstraint::remap_resid(
	SequenceMapping const & seqmap
) const {
	Size newseqpos1( seqmap[ seqpos1_ ] );
	if ( newseqpos != 0 ) {
		TR(t_debug) << "Remapping resid " << seqpos1_ << " to " << newseqpos1 << std::endl;

		utility::vector1< AtomID > new_atomids;
		for ( utility::vector1< AtomID >::const_iterator at_it( atom_ids_.begin() ), end( atom_ids_.end() ); at_it != end; ++at_it ) {
			if ( seqmap[ at_it->rsd() ] != 0 ) {
				new_atomids.push_back( AtomID( at_it->atomno(), seqmap[ at_it->rsd() ] ) );
			}
		}
		return new SequenceCouplingConstraint(	newseqpos1, newseqpos2,  new_atomids, sequence_coupling_ );
	}
	else return NULL;
}
*/

// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
SequenceCouplingConstraint::score(
	XYZ_Func const & xyz_func,
	EnergyMap const & weights,
	EnergyMap & emap
) const
{
	if ( weights[ this->score_type() ] == 0 ) return; // what's the point?
	runtime_assert( sequence_coupling_ );

	chemical::AA aa1( xyz_func.residue( seqpos1_ ).type().aa() );
	chemical::AA aa2( xyz_func.residue( seqpos2_ ).type().aa() );
	if ( seqpos1_ > sequence_coupling_->npos() || seqpos2_ > sequence_coupling_->npos()) return; // safety/relevance check

	Size edgeId = sequence_coupling_->findEdgeId(seqpos1_, seqpos2_);
	Real score( 0);
	if(edgeId<=0){//direction important. if (i,j) edge exists, will not return if (j,i) asked 
		 edgeId = sequence_coupling_->findEdgeId(seqpos2_, seqpos1_);

		 //epot = sequence_coupling_->edgePotBetween(edgeId); 
		 utility::vector1< utility::vector1 < Real > > const &  epot(sequence_coupling_->edgePotBetween(edgeId)); 
		 //runtime_assert(epot);
		 /*
		 if(epot==NULL){
			std::cerr << "no such edge " << seqpos1_ << " " << seqpos2_ << " in sequence coupling !)" << std::endl;
			utility_exit();
		 }
		 */
		 score = epot[aa2][aa1];
	}else{
		utility::vector1< utility::vector1 < Real > > const &  epot(sequence_coupling_->edgePotBetween(edgeId)); 
		//runtime_assert(epot);
		score = epot[aa1][aa2];
	}

	TR(t_trace) << "seqpos1 " << seqpos1_ << " aa1 " << aa1 << " " << "seqpos2 " << seqpos2_<< " aa2 "<< aa2 << " " << score << std::endl;

	emap[ this->score_type() ] += score;//pot scores are like energies; lower better.
}

void
SequenceCouplingConstraint::fill_f1_f2(
	AtomID const & ,//atom,
	XYZ_Func const & ,//conformation,
	Vector & ,//f1,
	Vector & ,//f2,
	EnergyMap const & //weights
) const
{
	// Do nothing, as the value of this function doesn't change with respect to
	// the torsions.
}

} // namespace constraints_additional
} // namespace protocols
