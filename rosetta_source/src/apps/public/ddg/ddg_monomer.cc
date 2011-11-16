// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Liz Kellogg ekellogg@u.washington.edu

#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <protocols/scoring/Interface.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ResidueSelector.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/ddg.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <basic/database/open.hh>

#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/random/random.hh>
#include <core/pack/task/ResfileReader.hh>

#include <fstream>
#include <iostream>
#include <sstream>
// AUTO-REMOVED #include <ios>
// AUTO-REMOVED #include <utility/io/izstream.hh>
#include <ObjexxFCL/format.hh>

// C++ headers
#include <cstdlib>
// Auto-header: duplicate removed #include <fstream>
// Auto-header: duplicate removed #include <iostream>
#include <string>
// Auto-header: duplicate removed #include <sstream>
#include <protocols/moves/ddGMover.hh>

//Auto Headers
#include <platform/types.hh>
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
#include <core/chemical/ResidueSelector.fwd.hh>
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
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
// AUTO-REMOVED #include <core/conformation/signals/ConnectionEvent.fwd.hh>
// AUTO-REMOVED #include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/JumpID.fwd.hh>
#include <core/id/JumpID.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
// AUTO-REMOVED #include <core/id/SequenceMapping.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/file_data.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/types.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.fwd.hh>
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/ResfileReader.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
// AUTO-REMOVED #include <core/scoring/MinimizationData.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Constraints.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Constraints.hh>
// AUTO-REMOVED #include <core/scoring/constraints/DOF_Constraint.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/DOF_Constraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Func.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Func.hh>
// AUTO-REMOVED #include <core/scoring/constraints/FuncFactory.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/HarmonicFunc.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/HarmonicFunc.hh>
// AUTO-REMOVED #include <core/scoring/constraints/XYZ_Func.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
// AUTO-REMOVED #include <protocols/jobdist/Jobs.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/scoring/Interface.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector0.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector0_bool.hh>
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
// AUTO-REMOVED #include <utility/file/gzip_util.hh>
// AUTO-REMOVED #include <utility/io/irstream.fwd.hh>
// AUTO-REMOVED #include <utility/io/irstream.hh>
#include <utility/io/izstream.fwd.hh>
// AUTO-REMOVED #include <utility/io/ozstream.fwd.hh>
// AUTO-REMOVED #include <utility/io/zipstream.hpp>
// AUTO-REMOVED #include <utility/io/zipstream.ipp>
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
#include <utility/tag/Tag.fwd.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
// AUTO-REMOVED #include <numeric/random/random.fwd.hh>
// AUTO-REMOVED #include <numeric/random/uniform.fwd.hh>
// AUTO-REMOVED #include <numeric/random/uniform.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/Fstring.fwd.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/byte.fwd.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/ubyte.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <execinfo.h>
#include <iomanip>
#include <iosfwd>
#include <istream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/options/option.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>
// AUTO-REMOVED #include <zlib/zlib.h>
// AUTO-REMOVED #include <zlib/zutil.h>

//Auto Headers

using basic::T;
using basic::Error;
using basic::Warning;
static basic::Tracer TR("apps.public.ddg.ddg_monomer");


//numeric::random::RandomGenerator RG(54324); // <- Magic number, do not change it!!!

using namespace core;
using namespace scoring;

typedef utility::vector1< core::chemical::AA > mutations;
typedef utility::vector1< double > ddgs;

void
print_ddgs(std::string ddg_out,
	   std::string label,
	   ddgs delta_e_components,
//mjo commenting out 'mut_avg_components' because it is unused and causes a warning
	   ddgs /*mut_avg_components*/,
	   double total_ddgs,
	   protocols::moves::ddGMover& mover,
	   bool print_header,
	   bool min_cst
	   ){
  std::ofstream ddg_output(ddg_out.c_str(), std::ios_base::app);
  if(!ddg_output){
    TR << "having trouble opening output file for dumping predicted ddgs"
       << ddg_out << std::endl;
    utility::exit(EXIT_FAILURE, __FILE__, __LINE__);
  }

  utility::vector1<std::string> scorefxn_header;
  if( min_cst ){
    mover.get_scorefunction_header(mover.minimization_score_function(),scorefxn_header);
  }else{
    mover.get_scorefunction_header(mover.score_function(),scorefxn_header);
  }

  if( print_header ){
    ddg_output << "ddG: description total ";
    for(Size i =1; i <=scorefxn_header.size();i++){
      ddg_output << scorefxn_header[i] << " ";
    }
    ddg_output << "\n";
  }
  if(label.compare("") != 0){
    ddg_output << "ddG: " << label << " " << ObjexxFCL::fmt::F(9,3,total_ddgs) << " ";
    for(Size m=1;m<=delta_e_components.size();m++){
      ddg_output << ObjexxFCL::fmt::F(9,3,delta_e_components[m]) << " ";
    }
    ddg_output << "\n";
  }

  ddg_output << std::endl;
}


/// @brief The input file is a list of mutation blocks.  Usually, this will be a set of point mutations.
/// where each "block" lists a single mutation.  However, it is possible to specify multiple mutations
/// together in a single block.
///
/// The file format is:
/// "total N"
/// followed by N blocks, where each block is
/// "M"
/// specifying followed by M lines of wt/resid/mutaa triples
/// "wtaa resid mutaa"
/// N, M and resid are all supposed to be integers.
/// wtaa, and mutaa are supposed to be 1-letter amino acid codes.
void
read_in_mutations(
	utility::vector1< mutations > & res_to_mut,
	std::string filename,
	pose::Pose & pose
)
{
	std::ifstream inputstream;
	inputstream.open(filename.c_str());
	if(inputstream.is_open()) {
		int total; std::string total_keyword;
		inputstream >> total_keyword;
		assert(total_keyword.compare("total") == 0);

		inputstream >> total; //keep for cross-checking
		while (!inputstream.eof()) {
			mutations current_mutation(pose.total_residue(),core::chemical::aa_unk);
			int num_mutations;
			inputstream >> num_mutations;
			while (num_mutations > 0) {
				char wt; int resnum; char mut;
				inputstream >> wt >> resnum >> mut;
				TR << "wt is " << wt << " resnum is " << resnum << " and mut is " << mut << std::endl;
				runtime_assert(pose.residue(resnum).name1() == wt); /// APL -- never use regular asserts when it comes to user input
				runtime_assert(core::chemical::oneletter_code_specifies_aa( mut ) ); /// APL -- input should specify the 1-letter code for an amino acid.
				core::chemical::AA mutation= core::chemical::aa_from_oneletter_code(mut);
				current_mutation[resnum]=mutation;
				num_mutations--; total--;
			}
			TR << "end reading mutations for this" << std::endl;
			if (num_mutations < 0) {
				TR.Error << "number of mutations mismatch! num_mutations < 0" << std::endl;
				return;
			} else {
				res_to_mut.push_back(current_mutation);
			}
		}
		if (total < 0) {
			TR.Error << "total number of mutations mismatch! total < 0" << std::endl;
			return;
		}
	}
}

int
main( int argc, char * argv [] )
{
	using namespace pose;
	using namespace scoring;
	using namespace conformation;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack::task;
	using namespace protocols::moves;

	OPT(ddg::weight_file);
	OPT(ddg::iterations);
	OPT(ddg::debug_output);
	OPT(ddg::dump_pdbs);
	OPT(ddg::out);
	OPT(ddg::interface_ddg);
	OPT(ddg::opt_radius);
	//	OPT(score::weights);
	//	OPT(score::patch);
	OPT(in::file::s);

	// setup random numbers and options
	devel::init(argc, argv);

	bool header_printed = false;

	// read the pose
	pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, basic::options::start_file() ); // gets filename from -s option

	std::string weight_file = option[ ddg::weight_file ]();

	/// Only change the fa_max_dis parameter if it is unspecified on the command line, otherwise, prefer the
	/// command line definition of the parameter.
	if ( ! basic::options::option[ score::fa_max_dis ].user() ) {
		TR << "option score::fa_max_dis unspecified on the command ine.  Setting fa_max_dis to 9.0 A." << std::endl;
		basic::options::option[ score::fa_max_dis ](9.0); //set fa_max_dis before scorefunction is created!
	} else {
		TR << "Using command line defined score::fa_max_dis of " <<  option[ score::fa_max_dis ] << " A." << std::endl;
	}

	ScoreFunctionOP score_structure_scorefxn(ScoreFunctionFactory::create_score_function(weight_file));

	ScoreFunctionOP minimize_sfxn;
	if(basic::options::option[ddg::minimization_scorefunction].user() && basic::options::option[ddg::minimization_patch].user()){
		minimize_sfxn=ScoreFunctionFactory::create_score_function(
			basic::options::option[ddg::minimization_scorefunction](),
			basic::options::option[ddg::minimization_patch]());
	}else if(basic::options::option[ddg::minimization_scorefunction].user()){
		minimize_sfxn=ScoreFunctionFactory::create_score_function(basic::options::option[ddg::minimization_scorefunction]());
	}else{
		minimize_sfxn=ScoreFunctionFactory::create_score_function(
			basic::database::full_name("scoring/weights/standard.wts"),
			basic::database::full_name("scoring/weights/score12.wts_patch"));
	}

	int num_iterations = option[ ddg::iterations ]();
	bool opt_nbrs = false;
	double cutoff = -1;
	if(basic::options::option[ ddg::opt_radius].user()){
		opt_nbrs = true;
		cutoff = basic::options::option[ ddg::opt_radius ]();
	} else if (basic::options::option[ddg::local_opt_only]()) {
		opt_nbrs = true;
		cutoff = 8.0; //default cutoff
	}

	//initialize output options.
	//debug output?
	bool debug_output = option[ ddg::debug_output ]();
	if(debug_output){
		TR << "weights being used: " <<
		score_structure_scorefxn->weights() << "\n";
	}

	//dump repacked pdbs?
	bool dump_pdbs = option[ ddg::dump_pdbs ]();

	//output ddgs into what file?
	std::string ddg_out = option[ ddg::out ]();


	//interface mode? setting = jump number to use for interface
	Size const interface_ddg = option[ ddg::interface_ddg ]();
	runtime_assert( interface_ddg <= pose.num_jump() );

	//minimize after repacking?
	bool min_cst = option[ddg::min_cst]();

	//take mean or min energy as predicted ddg?
	bool mean = option[ddg::mean]();
	bool min = option[ddg::min]();

	ObjexxFCL::FArray2D<double> wt_scores(20,num_iterations);

	utility::vector1<core::chemical::AA> all_unk(pose.total_residue(),core::chemical::aa_unk);

	utility::vector1<double> wt_averaged_score_components;
	utility::vector1<ddgs> delta_energy_components;
	utility::vector1<double> total_ddgs;
	utility::vector1<ddgs> mutant_averaged_score_components;
	utility::vector1<std::string> delta_delta_G_label;

	ddGMover get_wt_score(score_structure_scorefxn,minimize_sfxn,all_unk);
	get_wt_score.set_min_cst(min_cst);
	get_wt_score.set_min(min);
	get_wt_score.set_mean(mean);
	if(!opt_nbrs){
		get_wt_score.restrict_to_nbrs(opt_nbrs);
		get_wt_score.neighbor_cutoff(cutoff);
		get_wt_score.num_iterations(num_iterations);
		get_wt_score.dump_pdbs(dump_pdbs);
		get_wt_score.is_interface_ddg(interface_ddg);
		get_wt_score.debug_output(debug_output);
		get_wt_score.num_iterations(num_iterations);
		get_wt_score.residues_to_mutate(all_unk);
		get_wt_score.apply(pose);
		wt_averaged_score_components=get_wt_score.get_wt_averaged_score_components();
	}

	if(option[ ddg::mut_file ].user()){//check if mutfile is specified
		TR << "reading in mutfile" << std::endl;
		std::string filename = option[ddg::mut_file]();
		utility::vector1< mutations > res_to_mut;
		read_in_mutations( res_to_mut, filename, pose);
		TR << "size of res_to_mut is: " << res_to_mut.size() << std::endl;
		//initialize wildtype scores
		for(Size i=1;  i <= res_to_mut.size(); i++){
			utility::vector1<core::chemical::AA> residues_to_mutate = res_to_mut[i];
			bool mutation_defined = false; //to check if any mutation was specified
			for(Size m =1; m<= residues_to_mutate.size(); m++){
				if(residues_to_mutate[m] != core::chemical::aa_unk){
					mutation_defined=true;
				}
			}
			if(mutation_defined){
				ddGMover point_mutation(score_structure_scorefxn,minimize_sfxn,residues_to_mutate);
				point_mutation.set_min_cst(min_cst);
				point_mutation.set_min(min);
				point_mutation.set_mean(mean);
				if(!opt_nbrs && get_wt_score.is_wt_calc_complete()){
					TR << "testing if wt calc is complete. should be complete!" << std::endl;
					point_mutation.wt_score_components(get_wt_score.wt_score_components());
				}
				point_mutation.restrict_to_nbrs(opt_nbrs);
				point_mutation.neighbor_cutoff(cutoff);
				point_mutation.dump_pdbs(dump_pdbs);
				point_mutation.debug_output(debug_output);
				point_mutation.num_iterations(num_iterations);
				point_mutation.apply(pose);
				delta_delta_G_label.push_back(point_mutation.mutation_label(pose));
				TR << "mutation label for this round is " << point_mutation.mutation_label(pose) << std::endl;
				if(point_mutation.is_wt_calc_complete() &&
					point_mutation.is_mutant_calc_complete()){
					//TR << " both calculations are complete so start storing info!" << std::endl;
					//output everything or store everything for output later
					delta_energy_components.push_back(point_mutation.get_delta_energy_components());
					mutant_averaged_score_components.push_back(point_mutation.get_mutant_averaged_score_components());
					total_ddgs.push_back(point_mutation.ddG());
					print_ddgs(ddg_out,
								  point_mutation.mutation_label(pose),
								  point_mutation.get_delta_energy_components(),
								  point_mutation.get_mutant_averaged_score_components(),
								  point_mutation.ddG(),
								  point_mutation,
								  !header_printed,
								  min_cst);
					if( ! header_printed ) {
						header_printed = true;
					}
				}
			}
		}
	}

	if(option[packing::resfile].user()){ //check is resfile is specified
		pack::task::PackerTaskOP storage_task(pack::task::TaskFactory::create_packer_task(pose));

		storage_task->initialize_from_command_line();
		pack::task::parse_resfile(pose, *storage_task);
		storage_task->or_include_current(true);

		for(Size i =1;i<=pose.total_residue();i++){
			if(storage_task->design_residue(i)){
				for(ResidueLevelTask::ResidueTypeCAPListConstIter aa_iter(storage_task->residue_task(i).allowed_residue_types_begin()),
					 aa_end(storage_task->residue_task(i).allowed_residue_types_end());
					 aa_iter != aa_end; ++aa_iter){
					utility::vector1<core::chemical::AA> residues_to_mutate(pose.total_residue(),core::chemical::aa_unk);
					residues_to_mutate[i]=((*aa_iter)->aa());
					if(residues_to_mutate[i] != core::chemical::aa_unk){
						ddGMover point_mutation(score_structure_scorefxn,minimize_sfxn,residues_to_mutate);
						point_mutation.set_min_cst(min_cst);
						point_mutation.set_min(min);
						point_mutation.set_mean(mean);
						//initialize wildtype scores
						if(!opt_nbrs && get_wt_score.is_wt_calc_complete()){
							TR << "testing if wt calc is complete. should be complete!" << std::endl;
							point_mutation.wt_score_components(get_wt_score.wt_score_components());
						}
						point_mutation.restrict_to_nbrs(opt_nbrs);
						point_mutation.neighbor_cutoff(cutoff);
						point_mutation.dump_pdbs(dump_pdbs);
						point_mutation.debug_output(debug_output);
						point_mutation.num_iterations(num_iterations);
						point_mutation.apply(pose);
						delta_delta_G_label.push_back(point_mutation.mutation_label(pose));
						if(point_mutation.is_wt_calc_complete() &&
							point_mutation.is_mutant_calc_complete()){
							//TR << " both calculations are complete so start storing info!" << std::endl;
							//output everything
							delta_energy_components.push_back(point_mutation.get_delta_energy_components());
							mutant_averaged_score_components.push_back(point_mutation.get_mutant_averaged_score_components());
							total_ddgs.push_back(point_mutation.ddG());
							//output information to file
							print_ddgs(ddg_out,
										  point_mutation.mutation_label(pose),
										  point_mutation.get_delta_energy_components(),
										  point_mutation.get_mutant_averaged_score_components(),
										  point_mutation.ddG(),
										  point_mutation,
										  !header_printed,
										  min_cst);
							if( ! header_printed ) {
								header_printed = true;
							}
						}
					}
				}
			}
		}
	}
	//INTERFACE MODE
	if(interface_ddg > 0){

		//TR << "[DEBUG]: now  in interface mode"<< std::endl;
		//detect interface residues
		using namespace core;
		using namespace core::conformation;
		using namespace core::chemical;

		utility::vector1<core::chemical::AA> residues_to_mutate;
		//set up interface object
		protocols::scoring::Interface protein_interface(interface_ddg);
		protein_interface.distance(10.0);
		protein_interface.calculate(pose);
		// protein_interface.print(pose); // unnecessary log output

		//debug statement
		for(Size i =1;i<=pose.total_residue();i++){
			if(protein_interface.is_interface(i)){
				TR.Debug << "[DEBUG]:" << i << " is in the interface " << std::endl;
			}
		}
		//debug statement end

		for(Size i =1;i<=pose.total_residue();i++){
			if(protein_interface.is_interface(i)){//is interface residue
				for(Size j =1; j <= 20 ; j++){ //iterate through all amino acids
					residues_to_mutate = all_unk; //mutate each interface residue one at a time
					core::chemical::AA curr_aa = (core::chemical::AA)j;
					if(curr_aa != pose.aa(i) && (pose.aa(i) != aa_unk)/*this hopefully will never happen?*/ ){
						residues_to_mutate[i]=curr_aa;
						ddGMover interface_mutation(score_structure_scorefxn,minimize_sfxn,residues_to_mutate);
						interface_mutation.set_min_cst(min_cst);
						interface_mutation.is_interface_ddg(interface_ddg);
						interface_mutation.set_min(min);
						interface_mutation.set_mean(mean);
						if(get_wt_score.is_wt_calc_complete()){
							TR << "testing if wt calc is complete. should be complete!" << std::endl;
							interface_mutation.wt_score_components(get_wt_score.wt_score_components());
							interface_mutation.wt_unbound_score_components(get_wt_score.wt_unbound_score_components());
						}

						if(dump_pdbs) interface_mutation.dump_pdbs(dump_pdbs);
						if(debug_output) interface_mutation.debug_output(debug_output);

						interface_mutation.num_iterations(num_iterations);
						interface_mutation.apply(pose);
						delta_delta_G_label.push_back(interface_mutation.mutation_label(pose));
						TR << "mutation label for this round is " << interface_mutation.mutation_label(pose) << std::endl;
						if(interface_mutation.is_wt_calc_complete() &&
							interface_mutation.is_mutant_calc_complete()){

							delta_energy_components.push_back(interface_mutation.get_delta_energy_components());
							mutant_averaged_score_components.push_back(interface_mutation.get_mutant_averaged_score_components());
							total_ddgs.push_back(interface_mutation.ddG());

							print_ddgs(ddg_out,
										  interface_mutation.mutation_label(pose),
										  interface_mutation.get_delta_energy_components(),
										  interface_mutation.get_mutant_averaged_score_components(),
										  interface_mutation.ddG(),
										  interface_mutation,
										  !header_printed,
										  min_cst);
							if( ! header_printed ) {
								header_printed = true;
							}
							TR << "interface mutation is complete and ddg is: " << interface_mutation.ddG() << std::endl;
						}
					}
				}//iterate through all amino acids
			}
		}
	}

	/**
		//format and output all the stored information
	 utility::vector1<std::string> scorefxn_header = get_wt_score.get_scorefunction_header(score_structure_scorefxn);


	 ddg_output <<"\n***********************************\n" <<
	 "ddG: description total ";
	 for(Size i =1; i <=scorefxn_header.size();i++){
		 ddg_output << scorefxn_header[i] << " ";
	 }
	 ddg_output << "\n***********************************\n";
	 for(Size c=1;c<=delta_delta_G_label.size();c++){
		 if(delta_delta_G_label[c].compare("") != 0){
			 ddg_output << "ddG: " << delta_delta_G_label[c] << " " << F(9,3,total_ddgs[c]) << " ";
			 ddgs ddg_score_components = delta_energy_components[c];
			 for(Size m=1;m<=ddg_score_components.size();m++){
				 ddg_output << F(9,3,ddg_score_components[m]) << " ";
			 }
			 ddg_output << "\n";
		 }
	 }
	 ddg_output << std::endl;
	 **/

}

