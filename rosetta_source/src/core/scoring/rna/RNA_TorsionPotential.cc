// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/RNA_TorsionPotential.cc
/// @brief  RNA_TorsionPotential potential class implementation
/// @author Rhiju Das

// Unit Headers
#include <core/scoring/rna/RNA_TorsionPotential.hh>


// Package Headers
#include <core/scoring/rna/RNA_Util.hh>

#include <core/chemical/VariantType.hh>
#include <core/chemical/AtomType.hh>  //Need this to prevent the compiling error: invalid use of incomplete type 'const struct core::chemical::AtomType Oct 14, 2009
#include <core/chemical/AtomTypeSet.hh>  //Need this to prevent the compiling error: invalid use of incomplete type 'const struct core::chemical::AtomType Oct 14, 2009


// Project Headers
#include <core/conformation/Conformation.hh>
#include <basic/database/open.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/CircularGeneral1D_Func.hh>
#include <core/scoring/constraints/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/FadeFunc.hh>
#include <core/scoring/constraints/SumFunc.hh>
#include <core/scoring/EnergyMap.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>
#include <numeric/deriv/dihedral_deriv.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
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
#include <core/chemical/types.hh>
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
#include <core/conformation/Residue.hh>
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
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
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
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
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
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/CircularHarmonicFunc.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/FadeFunc.fwd.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/FuncFactory.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/SumFunc.fwd.hh>
#include <core/scoring/constraints/XYZ_Func.fwd.hh>
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/rna/RNA_TorsionPotential.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
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
#include <utility/io/izstream.fwd.hh>
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
#include <numeric/NumericTraits.hh>
#include <numeric/constants.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/sphericalVector.hh>
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
#include <numeric/random/random.fwd.hh>
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
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
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
#include <cstdlib>
#include <execinfo.h>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end
///////////////////////////


static basic::Tracer tr( "core.scoring.rna.RNA_TorsionPotential" );

using namespace ObjexxFCL::fmt;
using numeric::conversions::radians;

namespace core {
namespace scoring {
namespace rna {

	RNA_TorsionPotential::RNA_TorsionPotential():
		//		path_to_torsion_files_( "scoring/rna/torsion_potentials/rd2008/" ),
		// path_to_torsion_files_( "scoring/rna/torsion_potentials/ps2009/" ),
//		path_to_torsion_files_( "scoring/rna/torsion_potentials/rd_ps_2010/" ),
		path_to_torsion_files_( basic::options::option[  basic::options::OptionKeys::score::rna_torsion_potential ]() ),
//			path_to_torsion_files_( "scoring/rna/torsion_potentials/Mar_28_revert_back_zeta/"),

		rna_tight_torsions_( true ),
		delta_fade_( 10.0 ),
		alpha_fade_( 10.0 ),
		verbose_( false )
	{
		//		init_potentials_from_gaussian_parameters();
		init_potentials_from_rna_torsion_database_files();

		init_fade_functions();

	}

	//////////////////////////////////////////////////////////////////////////
	Real const
	RNA_TorsionPotential::eval_intrares_energy(core::conformation::Residue const & rsd, pose::Pose const & pose) const
	{
		using namespace core::id;

		if (!rsd.is_RNA() ) return 0.0;

		if(verbose_){
			std::cout << std::endl;
			std::cout << "Intra_res: " << " rsd.seqpos()= " << rsd.seqpos() << std::endl;
			std::cout << std::endl;
		}
		Real score( 0.0 );

		Real const beta= numeric::principal_angle_degrees( rsd.mainchain_torsion( BETA ) );
		Real const gamma= numeric::principal_angle_degrees(  rsd.mainchain_torsion( GAMMA ) );
		Real const delta= numeric::principal_angle_degrees(  rsd.mainchain_torsion( DELTA ) );
		Real const chi= numeric::principal_angle_degrees(  rsd.chi( CHI - NUM_RNA_MAINCHAIN_TORSIONS ) );
		Real const nu2= numeric::principal_angle_degrees(  rsd.chi( NU2 - NUM_RNA_MAINCHAIN_TORSIONS ) );
		Real const nu1= numeric::principal_angle_degrees(  rsd.chi( NU1 - NUM_RNA_MAINCHAIN_TORSIONS ) );

		if(verbose_) std::cout << "Beta torsion" << std::endl;

		if(Should_score_torsion(pose, TorsionID( rsd.seqpos(), id::BB, BETA ) ) ) {
			score += beta_potential_->func( beta ); //beta
		}

		if(verbose_) std::cout << "Gamma torsion" << std::endl;

		if(Should_score_torsion(pose, TorsionID( rsd.seqpos(), id::BB, GAMMA ) ) ) {
			score += gamma_potential_->func( gamma ); //gamma
		}

		if(verbose_) std::cout << "Delta torsion" << std::endl;


		if(Should_score_torsion(pose, TorsionID( rsd.seqpos(), id::BB, DELTA ) ) ) {
			score += ( fade_delta_north_->func( delta ) * delta_north_potential_->func( delta ) +
								 fade_delta_south_->func( delta ) * delta_south_potential_->func( delta ) ); //delta
		}

		if(verbose_)  std::cout << "Chi torsion" << std::endl;

		if(Should_score_torsion(pose, TorsionID( rsd.seqpos(), id::CHI, CHI - NUM_RNA_MAINCHAIN_TORSIONS ) ) ) {
			score += ( fade_delta_north_->func( delta ) * chi_north_potential_->func( chi ) +
							 fade_delta_south_->func( delta ) * chi_south_potential_->func( chi ) ); //chi
		}

		if(verbose_)  std::cout << "nu2 torsion" << std::endl;

		if(Should_score_torsion(pose, TorsionID( rsd.seqpos(), id::CHI, NU2 - NUM_RNA_MAINCHAIN_TORSIONS ) ) ) {
			score += ( fade_delta_north_->func( delta ) * nu2_north_potential_->func( nu2 ) +
								 fade_delta_south_->func( delta ) * nu2_south_potential_->func( nu2 ) ); //nu2
		}

		if(verbose_)  std::cout << "nu1 torsion" << std::endl;

		if(Should_score_torsion(pose, TorsionID( rsd.seqpos(), id::CHI, NU1 - NUM_RNA_MAINCHAIN_TORSIONS ) ) ) {
			score += ( fade_delta_north_->func( delta ) * nu1_north_potential_->func( nu1 ) +
							 fade_delta_south_->func( delta ) * nu1_south_potential_->func( nu1 ) ); //nu1
		}

		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////



/*
		std::cout << "Intra_res: " << " rsd.seqpos()= " << rsd.seqpos() << std::endl;

		Real beta_score  = beta_potential_->func( beta );

		Real gamma_score = gamma_potential_->func( gamma );

		Real delta_score = ( fade_delta_north_->func( delta ) * delta_north_potential_->func( delta ) +
							 fade_delta_south_->func( delta ) * delta_south_potential_->func( delta ) );

		Real chi_score = ( fade_delta_north_->func( delta ) * chi_north_potential_->func( chi ) +
							 fade_delta_south_->func( delta ) * chi_south_potential_->func( chi ) );

		Real nu2_score =( fade_delta_north_->func( delta ) * nu2_north_potential_->func( nu2 ) +
							 fade_delta_south_->func( delta ) * nu2_south_potential_->func( nu2 ) );

		Real nu1_score=( fade_delta_north_->func( delta ) * nu1_north_potential_->func( nu1 ) +
							 fade_delta_south_->func( delta ) * nu1_south_potential_->func( nu1 ) );
		Real cutoff_score=0.1;

		if(beta_score>cutoff_score){
			std::cout << "beta= " << beta << " " << beta_score << std::endl;
		}
		if(gamma_score>cutoff_score){
			std::cout << "gamma= " << gamma << " " << gamma_score << std::endl;
		}
		if(delta_score>cutoff_score){
			std::cout << "delta= " << delta << " " << delta_score << std::endl;
		}
		if(chi_score>cutoff_score){
			std::cout << "chi= " << chi << " " << chi_score << std::endl;
		}
		if(nu2_score>cutoff_score){
			std::cout << "nu2= " << nu2 << " " << nu2_score << std::endl;
		}
		if(nu1_score>cutoff_score){
			std::cout << "nu1= " << nu1 << " " << nu1_score << std::endl;
		}
*/

/*
		if(rsd.seqpos()==5){

			std::ofstream outfile;
			std::string data_filename="./Delta_Chi.txt";
//			std::string data_filename="./Alpha_zeta_torsion_data.txt";
			outfile.open (data_filename.c_str(), std::ios_base::out | std::ios_base::app);




			Real delta_score=( fade_delta_north_->func( delta ) * delta_north_potential_->func( delta ) +
							 fade_delta_south_->func( delta ) * delta_south_potential_->func( delta ) );


			Real chi_score=( fade_delta_north_->func( delta ) * chi_north_potential_->func( chi ) +
							 fade_delta_south_->func( delta ) * chi_south_potential_->func( chi ) );

			Size spacing=10;

		 	outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left << delta;
		 	outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left << chi;
		 	outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left << delta_score;
		 	outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left << chi_score;
			outfile << "\n";

			outfile.flush();
			outfile.close();

		}
*/
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////


		///////////////////////////////////////////////////////////////////////////////////////


		return score;

	}

	//////////////////////////////////////////////////////////////////////////
	Real const
	RNA_TorsionPotential::residue_pair_energy(core::conformation::Residue const & rsd1, core::conformation::Residue const & rsd2, pose::Pose const & pose) const
	{

		using namespace core::id;

		// CAN WE ASSUME rsd1 < rsd2 ???

//THESE ASSERTION FAIL!! TRY RERUN ~/minirosetta/test/Mar_27_torsion_potential_fix/1zih/new_RNA_TorssionPotential/trail_2/output.txt
//		if((rsd1.seqpos() < rsd2.seqpos())==false){
//			std::cout << "rsd1.seqpos()= " << rsd1.seqpos() << " rsd2.seqpos()= " << rsd2.seqpos() << std::endl;
//			utility_exit_with_message( "(rsd1.seqpos() < rsd2.seqpos())==false" );
//		}

//		assert( rsd1.seqpos() < rsd2.seqpos() );

		if ( rsd1.seqpos() != (rsd2.seqpos() - 1) ) return 0.0;
		if (!rsd1.is_RNA() ) return 0.0;
		if (!rsd2.is_RNA() ) return 0.0;

		if(verbose_)  {
			std::cout << std::endl;
			std::cout << "Between_res= " << " rsd1.seqpos()= " << rsd1.seqpos() << " rsd2.seqpos()= " << rsd2.seqpos() << std::endl;
			std::cout << std::endl;
		}

		Real score( 0.0 );

		Real const delta= numeric::principal_angle_degrees(  rsd1.mainchain_torsion( DELTA ) );
		Real const epsilon= numeric::principal_angle_degrees(  rsd1.mainchain_torsion( EPSILON ) );
		Real const zeta= numeric::principal_angle_degrees(  rsd1.mainchain_torsion( ZETA ) );
		Real const next_alpha= numeric::principal_angle_degrees(  rsd2.mainchain_torsion( ALPHA ) );

		if(verbose_)  std::cout << "epsilon torsion" << std::endl;

		if(Should_score_torsion(pose, TorsionID( rsd1.seqpos(), id::BB, EPSILON ) ) ) {
			score += ( fade_delta_north_->func( delta ) * epsilon_north_potential_->func( epsilon ) +
								 fade_delta_south_->func( delta ) * epsilon_south_potential_->func( epsilon ) );
		}

		if(verbose_)  std::cout << "zeta torsion" << std::endl;

		if(Should_score_torsion(pose, TorsionID( rsd1.seqpos(), id::BB, ZETA ) ) ) {
			score += ( fade_alpha_sc_minus_->func( next_alpha ) * zeta_alpha_sc_minus_potential_->func( zeta ) +
							 fade_alpha_sc_plus_->func( next_alpha )  * zeta_alpha_sc_plus_potential_->func( zeta ) +
							 fade_alpha_ap_->func( next_alpha )       * zeta_alpha_ap_potential_->func( zeta ) );
		}

		if(verbose_)  std::cout << "next alpha torsion" << std::endl;

		if(Should_score_torsion(pose, TorsionID( rsd2.seqpos(), id::BB, ALPHA ) ) ) {
			score += alpha_potential_->func( next_alpha );
		}

		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////


/*

		std::cout << "Between_res= " << " rsd1.seqpos()= " << rsd1.seqpos() << " rsd2.seqpos()= " << rsd2.seqpos() << std::endl;
		Real delta_score = ( fade_delta_north_->func( delta ) * delta_north_potential_->func( delta ) +
							 fade_delta_south_->func( delta ) * delta_south_potential_->func( delta ) );

		Real epsilon_score = ( fade_delta_north_->func( delta ) * epsilon_north_potential_->func( epsilon ) +
							 fade_delta_south_->func( delta ) * epsilon_south_potential_->func( epsilon ) );

		Real zeta_score = ( fade_alpha_sc_minus_->func( next_alpha ) * zeta_alpha_sc_minus_potential_->func( zeta ) +
							 fade_alpha_sc_plus_->func( next_alpha )  * zeta_alpha_sc_plus_potential_->func( zeta ) +
							 fade_alpha_ap_->func( next_alpha )       * zeta_alpha_ap_potential_->func( zeta ) );

		Real alpha_score = alpha_potential_->func( next_alpha );

		Real cutoff_score=0.1;

		if(delta_score>cutoff_score){
			std::cout << "delta_check= " << delta << " " << delta_score << std::endl;
		}
		if(epsilon_score>cutoff_score){
			std::cout << "epsilon= " << epsilon << " " << epsilon_score << std::endl;
		}

		if(zeta_score>cutoff_score){
			std::cout << "zeta= " << zeta << " " << zeta_score << std::endl;
		}
		if(alpha_score>cutoff_score){
			std::cout << "alpha= " << next_alpha << " " << alpha_score << std::endl;
		}
*/
/*
		if(zeta_score>cutoff_score || alpha_score>cutoff_score){
			std::cout << "zeta= " << zeta << " " << zeta_score << std::endl;
			std::cout << "alpha= " << next_alpha << " " << alpha_score << std::endl;
		}

		if(rsd1.seqpos()==5){

			std::cout << "ALPHA= " << ALPHA << std::endl;
			std::cout << "BETA= " << BETA << std::endl;
			std::cout << "GAMMA= " << GAMMA << std::endl;
			std::cout << "DELTA= " << DELTA << std::endl;
			std::cout << "EPSILON= " << EPSILON << std::endl;
			std::cout << "ZETA= " << ZETA << std::endl;
			std::cout << "CHI - NUM_RNA_MAINCHAIN_TORSIONS= " << CHI - NUM_RNA_MAINCHAIN_TORSIONS << std::endl;
			std::cout << "NU2 - NUM_RNA_MAINCHAIN_TORSIONS= " << NU2 - NUM_RNA_MAINCHAIN_TORSIONS << std::endl;
			std::cout << "NU1 - NUM_RNA_MAINCHAIN_TORSIONS= " << NU1 - NUM_RNA_MAINCHAIN_TORSIONS << std::endl;
			std::cout << "NUM_RNA_MAINCHAIN_TORSIONS= " << NUM_RNA_MAINCHAIN_TORSIONS << std::endl;
			std::cout << "CHI= " << CHI << std::endl;
			std::cout << "NU2= " << NU2 << std::endl;
			std::cout << "NU1= " << NU1 << std::endl;



			std::ofstream outfile;
			std::string data_filename="./Delta_Chi.txt";
//			std::string data_filename="./Alpha_zeta_torsion_data.txt";
			outfile.open (data_filename.c_str(), std::ios_base::out | std::ios_base::app);




			Real alpha_score=alpha_potential_->func( next_alpha );


			Real zeta_score=( fade_alpha_sc_minus_->func( next_alpha ) * zeta_alpha_sc_minus_potential_->func( zeta ) +
								 fade_alpha_sc_plus_->func( next_alpha )  * zeta_alpha_sc_plus_potential_->func( zeta ) +
								 fade_alpha_ap_->func( next_alpha )       * zeta_alpha_ap_potential_->func( zeta ) );

			Size spacing=10;

		 	outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left << zeta;
		 	outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left << next_alpha;
		 	outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left << zeta_score;
		 	outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left << alpha_score;
			outfile << "\n";

			outfile.flush();
			outfile.close();

		}
*/

		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////
		/////////////////////////////////////DEBUG////////////////////////////////////////////////

		return score;
	}

	///////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_TorsionPotential::get_f1_f2( core::id::TorsionID const & torsion_id, pose::Pose const & pose, core::id::AtomID const & id, Vector & f1, Vector & f2 ) const
	{

			conformation::Conformation const & conformation( pose.conformation() );

			// Check that torsion is intraresidue.
			id::AtomID id1,id2,id3,id4;
			if  ( conformation.get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 ) ) return false;

			//Kinda hacky, but most succinct this way
			if(verbose_) {
				tr.Info << "In get_f1_f2_function: ";
				tr.Info << " atom_id: " << id << std::endl;
			}
		  if(Should_score_torsion(pose, torsion_id)==false) return false;

			Real theta( 0.0 );
			// copied from DihedralConstraint. Better work, damnit.
			if ( id == id1 ) {
				numeric::deriv::dihedral_p1_cosine_deriv( conformation.xyz( id1 ),conformation.xyz( id2 ),
																									conformation.xyz( id3 ),conformation.xyz( id4 ), theta, f1, f2 );
			} else if ( id == id2 ) {
				numeric::deriv::dihedral_p2_cosine_deriv( conformation.xyz( id1 ),conformation.xyz( id2 ),
																									conformation.xyz( id3 ),conformation.xyz( id4 ), theta, f1, f2 );
			} else if ( id == id3 ) {
				numeric::deriv::dihedral_p2_cosine_deriv( conformation.xyz( id4 ),conformation.xyz( id3 ),
																									conformation.xyz( id2 ),conformation.xyz( id1 ), theta, f1, f2 );
			} else if ( id == id4 ) {
				numeric::deriv::dihedral_p1_cosine_deriv( conformation.xyz( id4 ),conformation.xyz( id3 ),
																									conformation.xyz( id2 ),conformation.xyz( id1 ), theta, f1, f2 );
			} else {
				return false;
			}
			return true;
	}

	//////////////////////////////////////////////////////////////
	void
	RNA_TorsionPotential::eval_atom_derivative(
																						 id::AtomID const & id,
																						 pose::Pose const & pose,
																						 EnergyMap const & weights,
																						 Vector & F1,
																						 Vector & F2
																						 ) const{

		Real const radians2degrees = 1.0 / radians( 1.0 );

		Size const current_seqpos( id.rsd() );
		Size const nres( pose.total_residue() );

//		if( Should_score_atom(pose, id)==false ) return; //Early return

		// Look for torsions in current, previous, and next residues that might match up with this atom_id.
		// This isn't totally efficient, but the extra checks should not be rate limiting...
		for ( int offset = -1; offset <= +1; offset++ ) {

			Size const seqpos = current_seqpos + offset;

			//Should we check that residue at seqpos is a RNA??? Parin Mar 28, 2010

			if (seqpos < 1 ) continue;
			if (seqpos > nres ) continue;

			conformation::Residue const & rsd( pose.residue( seqpos ) );

			Real const alpha= numeric::principal_angle_degrees(   rsd.mainchain_torsion( ALPHA ) );
			Real const beta = numeric::principal_angle_degrees(   rsd.mainchain_torsion( BETA ) );
			Real const gamma= numeric::principal_angle_degrees(   rsd.mainchain_torsion( GAMMA ) );
			Real const delta= numeric::principal_angle_degrees(   rsd.mainchain_torsion( DELTA ) );
			Real const epsilon= numeric::principal_angle_degrees( rsd.mainchain_torsion( EPSILON ) );
			Real const zeta= numeric::principal_angle_degrees(    rsd.mainchain_torsion( ZETA ) );
			Real const chi  = numeric::principal_angle_degrees(   rsd.chi( CHI - NUM_RNA_MAINCHAIN_TORSIONS ) );
			Real const nu2  = numeric::principal_angle_degrees(   rsd.chi( NU2 - NUM_RNA_MAINCHAIN_TORSIONS ) );
			Real const nu1  = numeric::principal_angle_degrees(   rsd.chi( NU1 - NUM_RNA_MAINCHAIN_TORSIONS ) );
			Vector f1( 0.0 ), f2( 0.0 );

			///////////////////////////////ALPHA//////////////////////////////
			if ( seqpos > 1 && pose.residue( seqpos - 1 ).is_RNA() && get_f1_f2( id::TorsionID( seqpos, id::BB, ALPHA ),	pose, id, f1, f2 ) ){

				Real dE_dtorsion = alpha_potential_->dfunc( alpha );

				Real const previous_zeta = numeric::principal_angle_degrees( pose.residue( seqpos - 1 ).mainchain_torsion( ZETA ) );  ///NEED TO CHANGE THIS AS WELL
				dE_dtorsion += ( fade_alpha_sc_minus_->dfunc( alpha ) * zeta_alpha_sc_minus_potential_->func( previous_zeta ) +
												 fade_alpha_sc_plus_->dfunc( alpha )  * zeta_alpha_sc_plus_potential_->func( previous_zeta ) +
												 fade_alpha_ap_->dfunc( alpha )       * zeta_alpha_ap_potential_->func( previous_zeta ) );

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			/////////////////////////////////Beta/////////////////////////////////////////////////////////
			if ( get_f1_f2( id::TorsionID( seqpos, id::BB, BETA ), pose, id, f1, f2 ) ){
				Real const dE_dtorsion = beta_potential_->dfunc( beta );

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			/////////////////////////////////Gamma/////////////////////////////////////////////////////////
			if ( get_f1_f2( id::TorsionID( seqpos, id::BB, GAMMA ),	pose, id, f1, f2 ) ){
				Real const dE_dtorsion = gamma_potential_->dfunc( gamma );

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			/////////////////////////////////Delta////////////////////////////////////////////////////////
			if ( get_f1_f2( id::TorsionID( seqpos, id::BB, DELTA ),	pose, id, f1, f2 ) ){
				Real dE_dtorsion = ( fade_delta_north_->dfunc( delta ) * delta_north_potential_->func( delta ) +
														 fade_delta_north_->func( delta )  * delta_north_potential_->dfunc( delta ) +
														 fade_delta_south_->dfunc( delta ) * delta_south_potential_->func( delta ) +
														 fade_delta_south_->func( delta )  * delta_south_potential_->dfunc( delta ) );

				dE_dtorsion += ( fade_delta_north_->dfunc( delta ) * chi_north_potential_->func( chi ) +
												 fade_delta_south_->dfunc( delta ) * chi_south_potential_->func( chi ) );

				dE_dtorsion += ( fade_delta_north_->dfunc( delta ) * nu2_north_potential_->func( nu2 ) +
												 fade_delta_south_->dfunc( delta ) * nu2_south_potential_->func( nu2 ) );

				dE_dtorsion += ( fade_delta_north_->dfunc( delta ) * nu1_north_potential_->func( nu1 ) +
												 fade_delta_south_->dfunc( delta ) * nu1_south_potential_->func( nu1 ) );

				if ( seqpos < nres && pose.residue( seqpos+1 ).is_RNA() ){
					dE_dtorsion += ( fade_delta_north_->dfunc( delta ) * epsilon_north_potential_->func( epsilon ) +
																	 fade_delta_south_->dfunc( delta ) * epsilon_south_potential_->func( epsilon ) );
				}

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			////////////////////////////////////EPSILON////////////////////////////////////////////////
			if ( seqpos < nres && pose.residue( seqpos+1 ).is_RNA() && get_f1_f2( id::TorsionID( seqpos, id::BB, EPSILON ),	pose, id, f1, f2 ) ){
				Real const dE_dtorsion = ( fade_delta_north_->func( delta ) * epsilon_north_potential_->dfunc( epsilon ) +
																	 fade_delta_south_->func( delta ) * epsilon_south_potential_->dfunc( epsilon ) );

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			////////////////////////////////////ZETA////////////////////////////////////////////////
			if ( seqpos < nres && pose.residue( seqpos+1 ).is_RNA() && get_f1_f2( id::TorsionID( seqpos, id::BB, ZETA ),	pose, id, f1, f2 ) ){

				Real const next_alpha = numeric::principal_angle_degrees(   pose.residue( seqpos+1 ).mainchain_torsion( ALPHA ) );  ///NEED TO CHANGE THIS AS WELL
				Real const dE_dtorsion = ( fade_alpha_sc_minus_->func( next_alpha ) * zeta_alpha_sc_minus_potential_->dfunc( zeta ) +
																	 fade_alpha_sc_plus_->func( next_alpha )  * zeta_alpha_sc_plus_potential_->dfunc( zeta ) +
																	 fade_alpha_ap_->func( next_alpha )       * zeta_alpha_ap_potential_->dfunc( zeta ) );

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			/////////////////////////////////////CHI///////////////////////////////////////////////////
			if ( get_f1_f2( id::TorsionID( seqpos, id::CHI, CHI - NUM_RNA_MAINCHAIN_TORSIONS ),	pose, id, f1, f2 ) ){
				Real const dE_dtorsion = ( fade_delta_north_->func( delta ) * chi_north_potential_->dfunc( chi ) +
																	 fade_delta_south_->func( delta ) * chi_south_potential_->dfunc( chi ) );

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			/////////////////////////////////////NU2//////////////////////////////////////////////////////
			if ( get_f1_f2( id::TorsionID( seqpos, id::CHI, NU2 - NUM_RNA_MAINCHAIN_TORSIONS ),	pose, id, f1, f2 ) ){
				Real const dE_dtorsion = ( fade_delta_north_->func( delta ) * nu2_north_potential_->dfunc( nu2 ) +
																	 fade_delta_south_->func( delta ) * nu2_south_potential_->dfunc( nu2 ) );

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			/////////////////////////////////////NU1//////////////////////////////////////////////////////
			if ( get_f1_f2( id::TorsionID( seqpos, id::CHI, NU1 - NUM_RNA_MAINCHAIN_TORSIONS ),	pose, id, f1, f2 ) ){
				Real const dE_dtorsion = ( fade_delta_north_->func( delta ) * nu1_north_potential_->dfunc( nu1 ) +
																	 fade_delta_south_->func( delta ) * nu1_south_potential_->dfunc( nu1 ) );

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}
		}

	}


	////////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_TorsionPotential::check_intra_residue( id::TorsionID const & torsion_id, pose::Pose const & pose, Size const seqpos ) const{

		// Check that torsion is intraresidue.
		id::AtomID id1,id2,id3,id4;
		bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
		if (fail) return false;

		if ( id1.rsd() != seqpos ) return false;
		if ( id2.rsd() != seqpos ) return false;
		if ( id3.rsd() != seqpos ) return false;
		if ( id4.rsd() != seqpos ) return false;

		return true;
	}

	///////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_TorsionPotential::init_potentials_from_rna_torsion_database_files() {

		// Initialize potential functions.
		initialize_potential_from_file( alpha_potential_, "alpha_potential.txt" );
		initialize_potential_from_file( beta_potential_, "beta_potential.txt" );
		initialize_potential_from_file( gamma_potential_, "gamma_potential.txt" );
		initialize_potential_from_file( delta_north_potential_, "delta_north_potential.txt" );
		initialize_potential_from_file( delta_south_potential_, "delta_south_potential.txt" );
		initialize_potential_from_file( epsilon_north_potential_, "epsilon_north_potential.txt" );
		initialize_potential_from_file( epsilon_south_potential_, "epsilon_south_potential.txt" );
		initialize_potential_from_file( zeta_alpha_sc_minus_potential_, "zeta_alpha_sc_minus_potential.txt" );
		initialize_potential_from_file( zeta_alpha_sc_plus_potential_, "zeta_alpha_sc_plus_potential.txt" );
		initialize_potential_from_file( zeta_alpha_ap_potential_, "zeta_alpha_ap_potential.txt" );
		initialize_potential_from_file( chi_north_potential_, "chi_north_potential.txt" );
		initialize_potential_from_file( chi_south_potential_, "chi_south_potential.txt" );
		initialize_potential_from_file( nu2_north_potential_, "nu2_north_potential.txt" );
		initialize_potential_from_file( nu2_south_potential_, "nu2_south_potential.txt" );
		initialize_potential_from_file( nu1_north_potential_, "nu1_north_potential.txt" );
		initialize_potential_from_file( nu1_south_potential_, "nu1_south_potential.txt" );

	}

	///////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_TorsionPotential::initialize_potential_from_file( core::scoring::constraints::FuncOP & func,
																	std::string const & filename ) {

		std::string const full_filename = basic::database::full_name( "scoring/rna/torsion_potentials/" + path_to_torsion_files_ + "/"+filename  );
		tr << "Reading in: " << full_filename << std::endl;
		func = new core::scoring::constraints::CircularGeneral1D_Func( full_filename );

	}


	///////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_TorsionPotential::init_fade_functions()
	{

		using namespace scoring::constraints;

		RNA_FittedTorsionInfo rna_torsion_fitted_info;
		Real const DELTA_CUTOFF_ = rna_torsion_fitted_info.delta_cutoff();

		// FadeFunc initialized with min, max, fade-width, and function value.
		fade_delta_north_ = new FadeFunc(
																		 -180.0 -delta_fade_,
																		 DELTA_CUTOFF_ + 0.5*delta_fade_ ,
																		 delta_fade_ ,
																		 1.0  );
		fade_delta_south_ = new FadeFunc(
																		 DELTA_CUTOFF_ - 0.5*delta_fade_,
																		 180.0 +delta_fade_ ,
																		 delta_fade_ ,
																		 1.0  );


		// FadeFunc initialized with min, max, fade-width, and function value.
		fade_alpha_sc_minus_ = new FadeFunc(
																				-120.0 - 0.5 * alpha_fade_ ,
																				0.0 + 0.5 * alpha_fade_ ,
																				alpha_fade_ ,
																				1.0  );
		fade_alpha_sc_plus_ = new FadeFunc(
																			 0.0 - 0.5 * alpha_fade_ ,
																			 100.0 + 0.5 * alpha_fade_ ,
																			 alpha_fade_ ,
																			 1.0  );

		fade_alpha_ap_ = new SumFunc();
		fade_alpha_ap_->add_func(
														 new FadeFunc(
																					-180.0 - alpha_fade_ ,
																					-120.0 + 0.5 * alpha_fade_ ,
																					alpha_fade_ ,
																					1.0  ) );

		fade_alpha_ap_->add_func(
														 new FadeFunc(
																					 100.0 - 0.5 * alpha_fade_ ,
																					 180.0 + alpha_fade_ ,
																					 alpha_fade_ ,
																					 1.0  ) );
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Defunct -- used for checking against original
	// (non-file-based) formulation of torsion potentials...
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_TorsionPotential::init_potentials_from_gaussian_parameters() {

		RNA_FittedTorsionInfo const rna_fitted_torsion_info;

		// Initialize potential functions.
		initialize_potential( alpha_potential_, rna_fitted_torsion_info.gaussian_parameter_set_alpha() );
		initialize_potential( beta_potential_, rna_fitted_torsion_info.gaussian_parameter_set_beta() );
		initialize_potential( gamma_potential_, rna_fitted_torsion_info.gaussian_parameter_set_gamma() );
		initialize_potential( delta_north_potential_, rna_fitted_torsion_info.gaussian_parameter_set_delta_north() );
		initialize_potential( delta_south_potential_, rna_fitted_torsion_info.gaussian_parameter_set_delta_south() );
		initialize_potential( epsilon_north_potential_, rna_fitted_torsion_info.gaussian_parameter_set_epsilon_north() );
		initialize_potential( epsilon_south_potential_, rna_fitted_torsion_info.gaussian_parameter_set_epsilon_south() );
		initialize_potential( zeta_alpha_sc_minus_potential_, rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_sc_minus() );
		initialize_potential( zeta_alpha_sc_plus_potential_, rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_sc_plus() );
		initialize_potential( zeta_alpha_ap_potential_, rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_ap() );
		initialize_potential( chi_north_potential_, rna_fitted_torsion_info.gaussian_parameter_set_chi_north() );
		initialize_potential( chi_south_potential_, rna_fitted_torsion_info.gaussian_parameter_set_chi_south() );
		initialize_potential( nu2_north_potential_, rna_fitted_torsion_info.gaussian_parameter_set_nu2_north() );
		initialize_potential( nu2_south_potential_, rna_fitted_torsion_info.gaussian_parameter_set_nu2_south() );
		initialize_potential( nu1_north_potential_, rna_fitted_torsion_info.gaussian_parameter_set_nu1_north() );
		initialize_potential( nu1_south_potential_, rna_fitted_torsion_info.gaussian_parameter_set_nu1_south() );

	}

	/////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_TorsionPotential::initialize_potential( core::scoring::constraints::FuncOP & func, Gaussian_parameter_set const & gaussian_parameter_set ) {

		using namespace core::scoring::constraints;

		Real const scale_rna_torsion_tether_( 0.05 ); // THIS IS A SCALING FACTOR FOR ALL CONSTRAINTS.
		Real const scale_rna_torsion_sd_( 1.0 / std::sqrt( scale_rna_torsion_tether_ ) );

		Real const xbin( 0.1 );
		Size const num_bins = static_cast<Size>( 360.0  /xbin );
		ObjexxFCL::FArray1D< Real > potential_values( num_bins, 0.0 );
		Real const xmin = -180.0 ;

		for ( Size m = 1; m <= gaussian_parameter_set.size(); m++ ) {
			// Must be in radians...
			CircularHarmonicFunc test_func( radians( gaussian_parameter_set[ m ].center ),
																			radians( scale_rna_torsion_sd_ * gaussian_parameter_set[ m ].width ) );

			for( Size i = 1; i <= num_bins; i++ ) {
				Real const bin_value = xmin + xbin * (i - 1 );
				Real const gaussian_value = test_func.func( radians( bin_value ) );
				if ( m == 1 || gaussian_value <= potential_values( i ) ) {
					potential_values( i ) = gaussian_value;
				}
			}
		}

		func = new CircularGeneral1D_Func( potential_values, xmin, xbin );

	}

	//HACKY	/////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_TorsionPotential::Output_boolean(std::string const & tag, bool boolean) const {

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::fmt;
		std::cout << tag;

		if(boolean==true){
			std::cout << A(4,"T");
		} else {
			std::cout << A(4,"F");
		}
	}
	/////////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_TorsionPotential::Is_chain_break_close_atom(core::conformation::Residue const & rsd, id::AtomID const & id) const{
		std::string const & atom_name=rsd.type().atom_name(id.atomno());

		if(atom_name=="OVU1" || atom_name=="OVL1" || atom_name=="OVL2"){
			return true;
		}else{
			return false;
		}

	}
	/////////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_TorsionPotential::Is_chain_break_close_torsion(pose::Pose const & pose, id::TorsionID const & torsion_id) const
	{
		using namespace ObjexxFCL;

		//We have to check that torsion is not a chain_break close.

//		bool Is_chain_break_close_torsion=false;

		Size torsion_seq_num=torsion_id.rsd();
		Size lower_seq_num=0;
		Size upper_seq_num=0;


		if ( torsion_id.type() != id::BB) return false;

		if ( torsion_id.torsion() == ALPHA){ //COULD BE A UPPER RESIDUE OF A CHAIN_BREAK_CLOSE

			lower_seq_num=torsion_seq_num-1;
			upper_seq_num=torsion_seq_num;

		}else if(torsion_id.torsion() == EPSILON || torsion_id.torsion() == ZETA){
			lower_seq_num=torsion_seq_num;
			upper_seq_num=torsion_seq_num+1;
		}else{
			if( torsion_id.torsion()!=DELTA && torsion_id.torsion() != BETA && torsion_id.torsion() != GAMMA){
				utility_exit_with_message("The torsion should be DELTA(lower), BETA(upper) or GAMMA(upper) !!" );
			}
			return false;
		}

		if(upper_seq_num==1) return false;

		if(lower_seq_num==pose.total_residue()) return false;

		if(pose.residue( lower_seq_num ).has_variant_type( chemical::CUTPOINT_LOWER ) ){
			if( pose.residue( upper_seq_num ).has_variant_type( chemical::CUTPOINT_UPPER )==false ){
				utility_exit_with_message("seq_num " + string_of(lower_seq_num) + " is a CUTPOINT_LOWER but seq_num " + string_of(upper_seq_num) + " is not a cutpoint CUTPOINT_UPPER??" );
			}
			return true;
		}

		return false;
	}
/*
	bool
	RNA_TorsionPotential::Should_score_atom(pose::Pose const & pose, id::AtomID const & id) const
	{

		conformation::Residue const & rsd=pose.residue(id.rsd());

		bool Is_virtual_atom=( rsd.atom_type(id.atomno()).name()=="VIRT");

		if(Is_virtual_atom && verbose_){
			tr.Info << " In Should_score_atomfunction" << std::endl;
			tr.Info << "  atom_id: " << id << std::endl;
			tr.Info << "  name: " << rsd.type().atom_name(id.atomno()) << std::endl;
			tr.Info << "  type: " << rsd.atom_type(id.atomno()).name() << std::endl;
			tr.Info << "		atom_type_index: " << rsd.atom_type_index( id.atomno()) << std::endl;
			tr.Info << "		atomic_charge: " << rsd.atomic_charge( id.atomno()) << std::endl;
		}

		bool Is_chain_break_close_atom_SECOND_METHOD = Is_chain_break_close_atom(rsd, id);

		Size const should_score_atom= (Is_virtual_atom==true && Is_chain_break_close_atom_SECOND_METHOD == false) ? false : true;

		if(verbose_){
			tr.Info << "  atom_id: " << id << std::endl;
			Output_boolean(" should_score_atom= ", should_score_atom);
			Output_boolean(" Is_chain_break_close_torsion= ", Is_chain_break_close_atom_SECOND_METHOD);
			Output_boolean(" Is_virtual_atom= ", Is_virtual_atom);
			std::cout << std::endl;
		}

		return should_score_atom;

	}
*/

	bool
	RNA_TorsionPotential::Should_score_torsion(pose::Pose const & pose, id::TorsionID const & torsion_id) const
	{

		id::AtomID id1,id2,id3,id4;
		pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );

		if(verbose_) std::cout << "torsion_id: " << torsion_id << std::endl;

		bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
		if (fail){
			if(verbose_) std::cout << "fail to get torsion!, perhap this torsion is located at a chain_break " << std::endl;
			return false;
//		 	utility_exit_with_message( "Fail to get torsion angle!" );
		}

		/////////////////Check for virtual atom (Last updated on Mar 27, 2010 Parin)/////////////////////////////////////////////////////////////////////
		conformation::Residue const & rsd_1=pose.residue(id1.rsd());
		conformation::Residue const & rsd_2=pose.residue(id2.rsd());
		conformation::Residue const & rsd_3=pose.residue(id3.rsd());
		conformation::Residue const & rsd_4=pose.residue(id4.rsd());

		bool Is_virtual_torsion=( rsd_1.atom_type(id1.atomno()).name()=="VIRT" || rsd_2.atom_type(id2.atomno()).name()=="VIRT"  || rsd_3.atom_type(id3.atomno()).name()=="VIRT"  || rsd_4.atom_type(id4.atomno()).name()=="VIRT");

		if(Is_virtual_torsion && verbose_){
			tr.Info << " Torsion containing one or more virtual atom(s)" << std::endl;
			tr.Info << "  torsion_id: " << torsion_id;
			tr.Info << "  atom_id: " << id1 << " " << id2 << " " << id3 << " " << id4 << std::endl;
			tr.Info << "  name: " << rsd_1.type().atom_name(id1.atomno()) << " " << rsd_2.type().atom_name(id2.atomno()) << " " << rsd_3.type().atom_name(id3.atomno()) << " " << rsd_4.type().atom_name(id4.atomno()) << std::endl;
			tr.Info << "  type: " << rsd_1.atom_type(id1.atomno()).name() << " " << rsd_2.atom_type(id2.atomno()).name() << " " << rsd_3.atom_type(id3.atomno()).name() << " " << rsd_4.atom_type(id4.atomno()).name() << std::endl;
			tr.Info << "		atom_type_index: " << rsd_1.atom_type_index( id1.atomno()) << " " << rsd_2.atom_type_index( id2.atomno()) << " " << rsd_3.atom_type_index( id3.atomno())  << " " << rsd_4.atom_type_index( id4.atomno()) << std::endl;
			tr.Info << "		atomic_charge: " << rsd_1.atomic_charge( id1.atomno())	<< " " << rsd_2.atomic_charge( id2.atomno())	<< " " << rsd_3.atomic_charge( id3.atomno())	<< " " << rsd_4.atomic_charge( id4.atomno()) << std::endl;
		}


		/////////////////////Check for chain_break_close (Since these will be contain virtual atom, but want score these torsions!///////////////////////////////
//		if(verbose_) std::cout << "Check_for_chain_break_close() first method" << std::endl;
		//Method 1:
		bool Is_chain_break_close_torsion_FIRST_METHOD = Is_chain_break_close_torsion(pose, torsion_id);

//		if(verbose_) std::cout << "Check_for_chain_break_close() second method" << std::endl;
    //Method 2:
		bool Is_chain_break_close_torsion_SECOND_METHOD=false;

		if(Is_chain_break_close_atom(rsd_1, id1) || Is_chain_break_close_atom(rsd_2, id2) || Is_chain_break_close_atom(rsd_3, id3) || Is_chain_break_close_atom(rsd_4, id4) ){
			Is_chain_break_close_torsion_SECOND_METHOD=true;
		}

		if( Is_chain_break_close_torsion_FIRST_METHOD != Is_chain_break_close_torsion_SECOND_METHOD){
			Output_boolean(" Is_chain_break_close_torsion_FIRST_METHOD= ", Is_chain_break_close_torsion_FIRST_METHOD);
			Output_boolean(" Is_chain_break_close_torsion_SECOND_METHOD= ", Is_chain_break_close_torsion_SECOND_METHOD);
			Output_boolean(" Is_virtual_torsion= ", Is_virtual_torsion);
			std::cout << std::endl;
			utility_exit_with_message( "Is_chain_break_close_torsion_FIRST_METHOD != Is_chain_break_close_torsion_SECOND_METHOD !!" );
		}

		if( Is_chain_break_close_torsion_FIRST_METHOD==true && Is_virtual_torsion==false){
			utility_exit_with_message( "Is_chain_break_close_torsion_FIRST_METHOD==true && Is_virtual_torsion==false !!" );
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Size const should_score_torsion= (Is_virtual_torsion==true && Is_chain_break_close_torsion_FIRST_METHOD==false) ? false : true;

		if(verbose_){
			Output_boolean(" should_score_torsion= ", should_score_torsion);
			Output_boolean(" Is_chain_break_close_torsion= ", Is_chain_break_close_torsion_FIRST_METHOD);
			Output_boolean(" Is_virtual_torsion= ", Is_virtual_torsion);
			std::cout << std::endl;
		}

		return should_score_torsion;


	}



}
}
}
