// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/ProteinSilentStruct.tmpl.hh
///
/// @brief Representation of rosetta++ protein silent-file structures.
/// @author James Thompson, Mike Tyka

#ifndef INCLUDED_core_io_silent_ProteinSilentStruct_tmpl_hh
#define INCLUDED_core_io_silent_ProteinSilentStruct_tmpl_hh

// Unit headers
#include <core/io/silent/ProteinSilentStruct.hh>

// C++ Headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <sstream>

// mini headers
#include <ObjexxFCL/char.functions.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>
#include <basic/prof.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedStubID.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/EnergyNames.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <numeric/model_quality/rms.hh>

// Boost headers
#include <boost/lexical_cast.hpp>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
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
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/GeneralEvent.fwd.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/signals/IdentityEvent.fwd.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.hh>
#include <core/conformation/symmetry/SymDof.fwd.hh>
#include <core/conformation/symmetry/SymSlideInfo.fwd.hh>
#include <core/conformation/symmetry/SymSlideInfo.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/types.hh>
#include <core/io/silent/SilentEnergy.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
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
#include <utility/signals/PausableSignalHub.fwd.hh>
#include <utility/signals/PausableSignalHub.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/types.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1A.fwd.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1P.fwd.hh>
#include <ObjexxFCL/FArray1P.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2P.fwd.hh>
#include <ObjexxFCL/FArray2P.hh>
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
#include <ObjexxFCL/StaticIndexRange.fwd.hh>
#include <ObjexxFCL/StaticIndexRange.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/byte.fwd.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/ubyte.fwd.hh>
#include <cassert>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <iomanip>
#include <iosfwd>
#include <istream>
#include <limits>
#include <ostream>
#include <set>
#include <time.h>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>


//using namespace ObjexxFCL;
//using namespace ObjexxFCL::format;

// option key includes
namespace core {
namespace io {
namespace silent {

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;
static basic::Tracer pss_tr("core.io.silent");

//explciit instantiation

//template class ProteinSilentStruct_Template<core::Real>;
//template class ProteinSilentStruct_Template<float>;



template <class T>
void
ProteinSilentStruct_Template<T>::fill_struct(
	core::pose::Pose const & pose,
	std::string tag
) {
	basic::ProfileThis doit( basic::SILENT_FILL_STRUCT );
	if ( !core::pose::is_ideal_pose(pose) ) {
		// Barak 07/09 - non-ideal poses are currently not supported for PSS
			// nor will they ever be!
			pss_tr.Error << "ERROR: trying to use a 'protein' type silent file for a non-ideal pose" << std::endl;
			pss_tr.Error << "consider using the '-out:file:silent_struct_type binary' flag" << std::endl;
			//			utility_exit();
	}

	SilentStruct::fill_struct( pose, tag );

	using namespace core::chemical;
	fullatom( pose.is_fullatom() );

	// conformation information
	resize( pose.total_residue() );
	static const std::string important_atom = "CA";
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		core::conformation::Residue const & resi = pose.residue(i);

		// skip VRT atoms
		if ( ( resi.aa() == core::chemical::aa_vrt ) || ( resi.aa() == core::chemical::aa_unk ) ) {
			phi      ( i, 0.0 );
			psi      ( i, 0.0 );
			omega    ( i, 0.0 );
			coords   ( i, resi.xyz(1) ); // for non-standard polymer residue, output coord of the first atom ).
			secstruct( i, 'L' );
			if ( fullatom() ) {
				//chi_[i].push_back(0.0);
				for ( Size kk = 1; kk <= max_chi(); ++kk ) {
					// give this non-standard residue 4 chi's if we're in fullatom mode.
					chi( i, kk, 0.0 );
				}
			} // if ( fullatom )
		} else {
			phi      ( i, resi.mainchain_torsion( 1 ) );
			psi      ( i, resi.mainchain_torsion( 2 ) );
			omega    ( i, resi.mainchain_torsion( 3 ) );
			coords   ( i, resi.xyz( important_atom ) );
			secstruct( i, pose.secstruct(i) );
			if ( fullatom() ) {
				//chi_[i] = resi.chi();
				for ( Size kk = 1; kk <= max_chi(); ++kk ) {
					if ( resi.nchi() >= kk ) {
						chi( i, kk, resi.chi(kk) );
					}
				}
			} // if ( fullatom )
		} // for ( unsigned int i = 1; i <= pose.total_residue(); ++i )
	} // for ( Size i = 1; i <= pose.total_residue; ++i )
	fold_tree( pose.fold_tree() );

	jumps_.clear();
	for ( Size nr = 1; nr <= fold_tree().num_jump(); nr++)  {
		add_jump( pose.jump(nr) );
	}

	fill_struct_with_residue_numbers( pose ); // grabs residue numbers from pose PDBInfo object.

	chain_endings( pose.conformation().chain_endings() );
}

template <class T>
bool ProteinSilentStruct_Template<T>::init_from_lines(
	utility::vector1< std::string > const & lines,
	SilentFileData & container
) {
	bool success( false );

	utility::vector1< std::string > energy_names;
	utility::vector1< std::string >::const_iterator iter = lines.begin();
	if ( iter->substr(0,9) != "SEQUENCE:" ) { //a full new header would change the default columns
		// get sequence and scorename data from the silent-file data object, because I don't have it!
		EnergyNamesOP enames = EnergyNamesOP(
			static_cast< EnergyNames * > ( container.get_shared_silent_data( energynames )() )
		);

		SimpleSequenceDataOP seqdata = SimpleSequenceDataOP(
			static_cast< SimpleSequenceData * > ( container.get_shared_silent_data( simplesequencedata )() )
		);

		sequence      ( seqdata->sequence()   );
		energy_names = enames ->energy_names();
	} else {
		// get sequence and scorename data from the first two lines provided, put
		// them into container for further use by other ProteinSilentStruct
		// objects.

		// first line is SEQUENCE:
		if ( !read_sequence( *iter ) ) return false;
		++iter;
		read_score_headers( *iter, energy_names, container ); ++iter;

	} // get header information

	// resize myself appropriately, according to length of sequence
	Size const total_residue = one_letter_sequence().length();
	resize( total_residue );
	utility::vector1< bool > read_flag( total_residue, false ); //to check whether all seqpos have been filled with info
	for ( utility::vector1< std::string >::const_iterator end = lines.end(); iter != end;	++iter ) {
		std::string tag;
		std::istringstream line_stream( *iter );

		if ( iter->substr(0,6) == "REMARK" ) {
// 			std::string tag;
// 			std::string comment;
// 			std::string value;
// 			runtime_assert( tag == "REMARK" );
// 			line_stream >> tag >> comment >> value;
			comment_from_line( *iter );//add_comment( comment, value );
			continue;  // don't skip comments
		}

		if ( iter->substr(0,7) == "SCORE: " ) { // SCORE: line with values from this structure.
			resize( total_residue ); // sequence_ should be defined by now.

			std::string tag;
			line_stream >> tag;
			if ( line_stream.fail() || tag != "SCORE:" ) {
				pss_tr.Error << "bad format in first score line of silent file" << std::endl;
				pss_tr.Error << "line = " << *iter << std::endl;
				pss_tr.Error << "tag = " << tag << std::endl;
			}

			parse_energies( line_stream, energy_names );
		} else { // conformation lines
			// parse fold_tree and jump lines
			if ( iter->substr(0,10) == "FOLD_TREE " ) {
				kinematics::FoldTree f;
				line_stream >> f;
				fold_tree( f ); // add fold-tree to this SilentStruct
				pss_tr.Debug << "read fold-tree " << f; //"\n" is in fold-tree output
				pss_tr.Debug << "reading " << f.num_jump() << " jumps " << std::endl;
				continue;
			} else if ( iter->substr(0,2) == "RT" ) {
				kinematics::Jump jump;
				line_stream >> jump;
				pss_tr.Debug << "read jump " << jump << std::endl;
				add_jump( jump );
				bJumps_use_IntraResStub_ = false; // modern style jumps, defined completely with the FoldTree
				continue;
			} else if ( iter->substr(0,9) == "SEQUENCE:" ) {
				pss_tr.Warning << "skipping duplicate sequence declaration " << std::endl;
				//after a SEQUENCE declaration we might find another SCORE header that should be skipped, too...
				utility::vector1< std::string >::const_iterator iter2 = ++iter;
				if (( iter2 != end ) && iter2->substr(0,7) == "SCORE: " ) {
					pss_tr.Warning << "re-reading score declaration from second line... " << std::endl;
					read_score_headers( *iter2, energy_names, container );
					++iter;
				}
				continue;
			} else if ( iter->substr(0,19) == "ANNOTATED_SEQUENCE:" ) {
				std::string annotated_seq;
				line_stream >> tag; //ANNOTATED_SEQUENCE
				line_stream >> annotated_seq;
				sequence( annotated_seq );
				continue;
			} else if ( iter->substr(0,4) == "JUMP" ) {
				// support for rosetta++ silent files
				std::string tag;
				Size nr;
				line_stream >> tag; //JUMP
				line_stream >> nr;
				if ( nr != fold_tree().num_jump() ) {
					pss_tr.Warning << "WARNING: corrupted silent file read line JUMP X -- X ";
					pss_tr.Warning << "should match number of jumps in FOLD_TREE " << std::endl;
				}
				for ( Size i = 1; i<= nr; i++ ) {
					kinematics::Jump jump;
					line_stream >> jump;
					add_jump( jump );
				}
				bJumps_use_IntraResStub_ = true; // jump is defined via N-C-CA rosetta++ style
				continue;
			} else if ( iter->substr( 0, 13 ) == "CHAIN_ENDINGS" ) {
				chain_endings_.clear();
				parse_chain_endings( line_stream );
				continue;
			} else if ( iter->substr(0,7) == "RES_NUM" ) {
				figure_out_residue_numbers_from_line( line_stream );
				continue;
			}

			// parse ss,torsions, and c-alpha coords
			Size seqpos;
			Real x, y, z, my_phi, my_psi, my_omega;
			char ss;

			line_stream >> tag;
			if ( !is_int( tag ) ) {
				pss_tr.Error << "ERROR:  !is_int( " << tag << " ) from line (" << *iter << ")" << std::endl;
				pss_tr.Error << "Are you trying to read a binary silent file ? Use -in:file:silent_struct_type binary " << std::endl;
				success = false;
				return success;
			}
			//assert( is_int( tag ) ); // this tag should represent the sequence position within the silent-file
			seqpos = int_of( tag );
			if ( seqpos <= 0 || seqpos > nres() ) {
				pss_tr.Error << "ERROR: incorrect sequence number " << seqpos << " (nres = "
					<< nres() << " ) from line (" << *iter << ")\n";
				success = false;
				return success;
			} else {
				read_flag[ seqpos ] = true;
			}
			line_stream >> ss >> my_phi >> my_psi >> my_omega >> x >> y >> z;
			numeric::xyzVector< T > temp_vec( x, y, z );

			phi      ( seqpos, my_phi   );
			psi      ( seqpos, my_psi   );
			omega    ( seqpos, my_omega );
			secstruct( seqpos, ss       );
			coords   ( seqpos, temp_vec );

			// Parse chi angles if this is a fullatom ProteinSilentStruct. Try to
			// figure out if the file format is fullatom by measuring if the tag is a
			// number.
			line_stream >> tag;
			Size chi_idx(1);
			while ( !line_stream.fail() ) {
				if ( is_float( tag ) ) {
					fullatom( true );
					//chis.push_back( float_of( tag ) );
					if ( chi_idx > max_chi() ) {
						pss_tr.Warning << "parse error (" << *iter << " ) " << tag << " has to many chi-angles" << std::endl;
						success = false;
						break;
					}
					chi( seqpos, chi_idx, float_of(tag) );
					++chi_idx;
				}
				//chi( seqpos, chis );
				line_stream >> tag;
			}

			if ( tag != decoy_tag() ) { // decoy_tag should be last tag.
				pss_tr.Warning 	<< "parse error(" << *iter << ") " << tag << " != " << decoy_tag() << std::endl;
				success = false;
				break;
			}
		} // conformation lines
	} // for ( iter ... )

	if ( fold_tree().num_jump() != jumps_.size() ) {
		pss_tr.Warning << "parse error:  found only " << jumps_.size()
			<< " RT-lines for a fold-tree with " << fold_tree().num_jump()
			<< " for decoy tag " << decoy_tag() << std::endl;
		return false;
	}

	// check if each sequence-position has been set
	utility::vector1< bool >::iterator bad_pos
		= find( read_flag.begin(), read_flag.end(), false );

	if ( bad_pos != read_flag.end() ) {
		pss_tr.Error << "ERROR: did not find coordinates for all sequence positions for "
			<< decoy_tag() << std::endl;
		core::Size idx = 1;
		for ( utility::vector1< bool >::iterator it = read_flag.begin(), end = read_flag.end();
					it != end; ++it, ++idx
		) {
			if ( !*it ) pss_tr.Error << "Couldn't read position " << idx << std::endl;
		}
		return false; //no success
	}

	// if no fold-tree available generate a standard tree
	if ( fold_tree().size() < 1 ) {
		fold_tree_.simple_tree( nres() );
		pss_tr.Debug << " generating simple fold-tree " << fold_tree();
	}

	if ( bJumps_use_IntraResStub_ ) { //for rosetta++ file-format
		//prepares of setting RT via N, CA, C
		fold_tree_.put_jump_stubs_intra_residue();
		//on could also think of making this a temporary change after read is finished
		// return to a standard fold_tree...
	}

	success = true;
	return success;
} // init_from_lines

template <class T>
void
ProteinSilentStruct_Template<T>::resize(
	Size const nres_in
) {
	nres( nres_in );
	secstruct_.resize( nres() );
	phi_      .resize( nres() );
	psi_      .resize( nres() );
	omega_    .resize( nres() );
	coords_   .resize( nres() );
	secstruct_.reserve( nres() );
	phi_      .reserve( nres() );
	psi_      .reserve( nres() );
	omega_    .reserve( nres() );
	coords_   .reserve( nres() );
	fold_tree_.simple_tree( nres() );

	// chi_ is a little bit more dangerous, so it has its own method.
	resize_chi();
} // resize

template <class T>
void
ProteinSilentStruct_Template<T>::resize_chi() {
	// make sure that there are the appropriate
	// number of chi angles in each position of
	// the chi_ vector1.
	chi_.resize( nres() );
	chi_.reserve( nres() );
	//for ( Size kk = 1; kk <= nres(); ++kk ) {
	//	if ( n_chi(kk) < max_chi() ) {
	//		chi_[ kk ].resize( max_chi(), 0.0 );
	//	}
	//} // kk
} // resize_chi

template <class T>
void ProteinSilentStruct_Template<T>::fill_pose(
	core::pose::Pose & pose
) const {
	using namespace core::chemical;
	ResidueTypeSetCAP residue_set;
	pss_tr.Debug << "fill_pose: SilentStruct is " << ( fullatom() ? "fullatom" : "centroid" ) << std::endl;
	if ( fullatom() ) {
		residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	} else {
		residue_set = ChemicalManager::get_instance()->residue_type_set( CENTROID );
	}
	fill_pose( pose, *residue_set );
} // fill_pose

template <class T>
void ProteinSilentStruct_Template<T>::fill_pose(
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const & residue_set
) const {
	basic::ProfileThis doit( basic::SILENT_FILL_POSE );
	runtime_assert( nres() != 0 );
	runtime_assert( sequence() != "" );
	//fpd  why was this called twice???
	//if ( pose.annotated_sequence() != sequence() ) {
	//	core::pose::make_pose_from_sequence( pose, sequence(), residue_set );
	//}
	//if ( pss_tr.Debug.visible() ) pss_tr.Debug << "FOLD TREE: " << fold_tree();
	core::pose::make_pose_from_sequence( pose, sequence(), residue_set );
	pss_tr.Debug << "FOLD TREE: " << fold_tree();

	// set fold_tree
	pose.fold_tree( fold_tree() );

	// set jumps
	for ( Size nr = 1; nr <= fold_tree().num_jump(); nr++)  {
		if ( !bJumps_use_IntraResStub_ ) { //default modern file-format
			pose.set_jump( nr, jump( nr ) );
		} else { //support for rosetta++ format
			Size start = fold_tree().jump_edge( nr ).start();
			Size stop = fold_tree().jump_edge( nr ).stop();
			id::StubID   up_stub( pose::named_stub_id_to_stub_id( id::NamedStubID("CA", "N", "CA", "C", start < stop ? start : stop  ), pose ));
			id::StubID down_stub( pose::named_stub_id_to_stub_id( id::NamedStubID("CA", "N", "CA", "C", start < stop ? stop  : start ), pose ));
			pose.conformation().set_stub_transform( up_stub, down_stub, jump( nr ) );
		}
	}


	if ( nres() != one_letter_sequence().length() ) {
		utility_exit_with_message( "RuntimeAssert failed: nres() == one_letter_sequence().length()" );
	}

	{ basic::ProfileThis doit( basic::SILENT_SET_POSE_COORDS );
		// angles and secondary structure
		for ( Size seqpos = 1; seqpos <= nres(); ++seqpos ) {
			if ( pose.residue_type(seqpos).is_protein() ) { //skip non-protein residues //ask for type to avoid updateing of coords.
				pose.set_phi   ( seqpos, phi  (seqpos) );
				pose.set_psi   ( seqpos, psi  (seqpos) );
				pose.set_omega ( seqpos, omega(seqpos) );
				if ( fullatom() ) {
					for ( Size jj = 1,
									end = std::min( pose.residue_type(seqpos).nchi(), n_chi(seqpos) );
								jj <= end; ++jj
					) {
						pose.set_chi( jj, seqpos, chi(seqpos,jj) );
					}
				} // if ( fullatom() )
			} // skip non-protein residues
			pose.set_secstruct( seqpos, secstruct(seqpos) );
		}

		if ( !chain_endings().empty() ) {
			pose.conformation().chain_endings( chain_endings() );
		}
	} //Profile scope

	core::pose::initialize_disulfide_bonds(pose);

	finish_pose( pose );
} // fill_pose


template <class T>
void ProteinSilentStruct_Template<T>::print_conformation( std::ostream & output ) const {
	// fold tree
	output << "REMARK PROTEIN SILENTFILE";
	if ( is_single_precision() ) output  << " SINGLE_PRECISION\n";
	else output << "\n";
	if ( fold_tree().size() > 1 ) { //assume non-trivial fold_tree only if more than one edge, i.e., EDGE 1 <nres> -1
		output << "FOLD_TREE ";
		for ( kinematics::FoldTree::const_iterator it = fold_tree().begin(), it_end = fold_tree().end();
					it != it_end; ++it ) {
			output << *it;
		}
		//		output << fold_tree(); this produces a new-line --- wrong behaviour of fold_tree but I don't want to fix 1000 u-tracer unit-tests!
		output << ' ' << decoy_tag() << "\n";
	}
	for ( Size i = 1; i <= fold_tree().num_jump(); i++ ) {
		output << jump( i ) << ' ' << decoy_tag() << "\n";
	}

	// sequence
	output << "ANNOTATED_SEQUENCE: " << sequence() << " " << decoy_tag() << "\n"; // print annotated_sequence per decoy

	// chain endings
	if ( !chain_endings().empty() ) {
		output << chain_endings_str() << ' ' << decoy_tag() << '\n';
	}

	// the actual dihedral/coordinate info
	//pss_tr.Debug << "FOLD_TREE Size: " << fold_tree().size() << " " << fold_tree() << std::endl;
	for ( Size i = 1; i <= nres(); ++i ) {
		output
			<< I( 4, i ) << ' '
			<< secstruct(i) << ' '
			<< F( 9, 3, phi(i) )
			<< F( 9, 3, psi(i) )
			<< F( 9, 3, omega(i) )
			<< F( 9, 3, coords(i).x() )
			<< F( 9, 3, coords(i).y() )
			<< F( 9, 3, coords(i).z() );

		if ( fullatom() ) {
			for ( Size chino = 1; chino <= max_chi(); ++chino ) {
				Real chi_to_print = 0.0;
				if ( chino <= n_chi(i) ) {
					chi_to_print = chi(i,chino);
				}
				output << F( 9, 3, chi_to_print );
			}
		} // if fullatom()

		output << ' ' << decoy_tag();
		output << "\n";
	} // for ( Size i = 1; i <= nres; ++i )
} // print_conformation

template <class T>
Real ProteinSilentStruct_Template<T>::get_debug_rmsd() {
	pose::Pose temp_pose;
	FArray2D< Real > rebuilt_coords (3, coords_.size() ), original_coords( 3, coords_.size() );
	static std::string atom_name = "CA";

	// build temp_pose from coordinates
	fill_pose( temp_pose );

	for ( Size i = 1; i <= temp_pose.total_residue(); ++i ) {
		for ( Size k = 1; k <= 3; ++k ) { // k = X, Y and Z
			rebuilt_coords (k,i) = temp_pose.residue(i).xyz( atom_name )[k-1];
			original_coords(k,i) = coords_[i][k-1];
		}
	}

	Real rmsd = numeric::model_quality::rms_wrapper( temp_pose.total_residue(), rebuilt_coords, original_coords );
	return rmsd;
}

template <class T>
ObjexxFCL::FArray2D< Real >
ProteinSilentStruct_Template<T>::get_CA_xyz() const {
	core::Size n_residues = nres();
	FArray2D< Real > my_coords( 3, n_residues );
	for ( Size i = 1; i <= n_residues; ++i ) { // i = n_residues
		for ( Size k = 1; k <= 3; ++k ) { // k = X, Y and Z
			my_coords(k,i) = coords_[i][k-1];
		} // k
	} // i

	return my_coords;
} // get_CA_positions

template <class T>
Real ProteinSilentStruct_Template<T>::CA_rmsd( ProteinSilentStruct_Template<T> other_pss ) {
	FArray2D< Real > my_coords    = get_CA_xyz();
	FArray2D< Real > other_coords = other_pss.get_CA_xyz();
	Real rmsd = numeric::model_quality::rms_wrapper( nres(), my_coords, other_coords );

	return rmsd;
} // ProteinSilentStruct_Template<T>::CA_rmsd

template <class T>
ProteinSilentStruct_Template<T> & ProteinSilentStruct_Template<T>::operator= (
	ProteinSilentStruct_Template<T> const & src
)
{
	resize( src.nres() );

	// per-residue conformation information
	for ( Size i = 1; i <= nres(); ++i ) {
		phi      ( i, src.phi(i)       );
		psi      ( i, src.psi(i)       );
		omega    ( i, src.omega(i)     );
		coords   ( i, src.coords(i)    );
		secstruct( i, src.secstruct(i) );
	}

	// fold-tree and jumps
	for ( Size jj = 1; jj <= src.njumps(); ++jj ) {
		add_rt( src.jump(jj) );
	}
	fold_tree( src.fold_tree() );

	chain_endings( src.chain_endings() );

	// add energy copying here

	return *this;
}

/// @brief parse the chain endings string from an input stream
template <class T>
void ProteinSilentStruct_Template<T>::parse_chain_endings( std::istream & stream ) {
	std::string s;
	stream >> s; // first column is "CHAIN_ENDINGS" tag, skip it

	utility::vector1< std::string > v;
	while ( stream.good() ) {
		stream >> s;
		v.push_back( s );
	}

	// remember to skip the last entry in the vector, which is the structure's nametag
	for ( Size i = 1, ie = v.size(); i < ie; ++i ) {
		add_chain_ending( boost::lexical_cast< Size >( v[ i ] ) );
	}
}

/// @brief return the chain endings string
template <class T>
std::string ProteinSilentStruct_Template<T>::chain_endings_str() const {
	std::ostringstream ss;
	ss << "CHAIN_ENDINGS ";

	for ( utility::vector1< Size >::const_iterator i = chain_endings().begin(), ie = chain_endings().end(); i != ie; ++i ) {
		ss << ' ' << (*i);
	}

	return ss.str();
}

template <class T>
void ProteinSilentStruct_Template<T>::chain_endings( utility::vector1< Size > const & endings ) {
	for ( utility::vector1< Size >::const_iterator i = endings.begin(), ie = endings.end(); i != ie; ++i ) {
		if ( (*i) < 1 || (*i) > nres() ) {
			pss_tr.Fatal << "ERROR: chain_endings() invalid chain ending " << (*i) << std::endl;
			utility_exit();
		}
	}

	chain_endings_ = endings;
	std::sort( chain_endings_.begin(), chain_endings_.end() ); // keep the list sorted
}

template <class T>
void ProteinSilentStruct_Template<T>::add_chain_ending( Size const seqpos ) {
	if ( seqpos < 1 || seqpos >= nres() ) {
		pss_tr.Fatal << "ERROR: add_chain_ending() invalid chain ending " << seqpos << std::endl;
		utility_exit();
	}

	chain_endings_.push_back( seqpos );
	std::sort( chain_endings_.begin(), chain_endings_.end() ); // keep the list sorted
}

template <class T>
Real ProteinSilentStruct_Template<T>::chi( Size const seqpos, Size const chi_num ) const {
	// error checking.
	if ( chi_num > max_chi() ) {
		std::string const msg(
			"Error: trying to chi " + string_of(chi_num) +
			" when max_chi is " + string_of( max_chi() ) + '\n'
		);
		utility_exit_with_message( msg );
	}

	// super-safe check to make sure that we don't return a memory access
	// violation on platforms like Windows in release mode.
	if ( chi_num > n_chi(seqpos) ) {
		pss_tr.Error << "Error: attempting to access chi that doesn't exist!"
			<< "(chi = " << chi_num << " seqpos = " << seqpos << ")"
			<< std::endl;;
		return 0.0;
	}
	return chi_[seqpos][chi_num];
}

template <class T>
void ProteinSilentStruct_Template<T>::chi( Size const seqpos, Size const chi_num, Real const chi ) {
	// error checking.
	if ( chi_num > max_chi() ) {
		std::string const msg(
			"Error: trying to set chi " + string_of(chi_num) +
			" when max_chi is " + string_of( max_chi() ) + '\n'
		);
		utility_exit_with_message( msg );
	}

	if ( chi_num > chi_[seqpos].size() ) {
		chi_[seqpos].resize( chi_num );
	}

	chi_[seqpos][chi_num] = chi;
} // chi

/// @brief returns the number of chis at this position.
template <class T>
Size ProteinSilentStruct_Template<T>::n_chi( Size const seqpos ) const {
	return chi_[seqpos].size();
}

template <class T>
void ProteinSilentStruct_Template<T>::chi( Size const seqpos, utility::vector1< T > const & chis ) {
	chi_[seqpos] = chis;
}

template <class T>
Size ProteinSilentStruct_Template<T>::max_chi() const {
	return max_chi_;
}
/// end accessors


template <class T>
core::Size ProteinSilentStruct_Template<T>::mem_footprint() const{
	typedef typename utility::vector1<  utility::vector1< T > >::const_iterator chi_iterator;
	T var;
	core::Size sum=
		sizeof(char)* secstruct_.capacity() +
		sizeof(var)* phi_.capacity() +
		sizeof(var)* psi_.capacity() +
		sizeof(var)* omega_.capacity() +
		sizeof(  numeric::xyzVector< T > )* coords_.capacity() +
		sizeof(kinematics::RT)* jumps_.capacity()+
		sizeof(Size) * chain_endings_.capacity()+
		sizeof( fold_tree_ ) +
		sizeof( symminfo_ );
	for( chi_iterator it = chi_.begin(); it < chi_.end(); ++it ) sum += (*it).capacity() * sizeof( T ) ;
	return sum;
}




} // namespace silent
} // namespace io
} // namespace core

#endif

