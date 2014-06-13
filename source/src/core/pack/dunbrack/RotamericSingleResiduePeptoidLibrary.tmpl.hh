// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/dunbrack/RotamericSingleResiduePeptoidLibrary.hh
/// @brief   Declaration of rotameric peptoid rotamer libraries
/// @author  P. Douglas Renfrew (renfrew@nyu.edu)


#ifndef INCLUDED_core_pack_dunbrack_RotamericSingleResiduePeptoidLibrary_tmpl_hh
#define INCLUDED_core_pack_dunbrack_RotamericSingleResiduePeptoidLibrary_tmpl_hh

// Unit Headers
#include <core/pack/dunbrack/RotamericSingleResiduePeptoidLibrary.hh>

// Package Headers
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/ChiSet.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/graph/Graph.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <basic/basic.hh>
#include <basic/interpolate.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/LexicographicalIterator.hh>

// Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

// Boost Headers
#include <boost/cstdint.hpp>

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
//XRW_B_T1
//#include <core/coarse/Translator.fwd.hh>
//XRW_E_T1
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/PDBPoseMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/pack/dunbrack/ChiSet.fwd.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>
#include <core/pack/dunbrack/RotamericSingleResiduePeptoidLibrary.fwd.hh>
#include <core/pack/dunbrack/SingleResiduePeptoidLibrary.fwd.hh>
#include <core/pack/dunbrack/SingleResiduePeptoidLibrary.hh>
#include <basic/MetricValue.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <utility/LexicographicalIterator.fwd.hh>
#include <utility/assert.hh>
#include <utility/down_cast.hh>
#include <utility/fixedsizearray1.fwd.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/gzip_util.hh>
#include <utility/io/irstream.fwd.hh>
#include <utility/io/irstream.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/io/mpistream.hh>
#include <utility/io/mpistream.ipp>
#include <utility/io/orstream.fwd.hh>
#include <utility/io/orstream.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/io/zipstream.hpp>
#include <utility/io/zipstream.ipp>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
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
#include <numeric/NumericTraits.hh>
#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/trig.functions.hh>
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
#include <numeric/random/uniform.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray2P.fwd.hh>
#include <ObjexxFCL/FArray2P.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray3.hh>
#include <ObjexxFCL/FArray.all.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/Fmath.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/StaticIndexRange.fwd.hh>
#include <ObjexxFCL/StaticIndexRange.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <zlib/zlib.h>
#include <zlib/zutil.h>


namespace core {
namespace pack {
namespace dunbrack {

/// DOUG DOUG DOUG This use to be PHIPSI_BINRANGE
template < Size T >
Real const RotamericSingleResiduePeptoidLibrary< T >::OMG_BINRANGE = 10.0;
template < Size T >
Real const RotamericSingleResiduePeptoidLibrary< T >::PHI_BINRANGE = 10.0;
template < Size T >
Real const RotamericSingleResiduePeptoidLibrary< T >::PSI_BINRANGE = 10.0;
template < Size T >
Real const RotamericSingleResiduePeptoidLibrary< T >::CIS_OMG_LOWER_RANGE = -30.0;
template < Size T >
Real const RotamericSingleResiduePeptoidLibrary< T >::CIS_OMG_UPPER_RANGE =  30.0;
template < Size T >
Real const RotamericSingleResiduePeptoidLibrary< T >::TRANS_OMG_LOWER_RANGE =  150.0;
template < Size T >
Real const RotamericSingleResiduePeptoidLibrary< T >::TRANS_OMG_UPPER_RANGE = -150.0;
template < Size T >
Real const RotamericSingleResiduePeptoidLibrary< T >::CIS_UPPER_TRANS_LOWER_SPLIT = 90.0;
template < Size T >
Real const RotamericSingleResiduePeptoidLibrary< T >::TRANS_UPPER_CIS_LOWER_SPLIT = -90.0;

template < Size T >
RotamericSingleResiduePeptoidLibrary< T >::RotamericSingleResiduePeptoidLibrary() :
	parent( T )
	//max_rotprob_( N_PHIPSI_BINS, N_PHIPSI_BINS, 0.0 )
{}

template < Size T >
RotamericSingleResiduePeptoidLibrary< T >::~RotamericSingleResiduePeptoidLibrary()
{}


template < Size T >
Real
RotamericSingleResiduePeptoidLibrary< T >::rotamer_energy(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const
{
	return eval_rotameric_energy_deriv( rsd, scratch, false );
}


template < Size T >
Real
RotamericSingleResiduePeptoidLibrary< T >::rotamer_energy_deriv(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const
{
	/// most of the work is done in this call
	Real score = eval_rotameric_energy_deriv( rsd, scratch, true );

	if ( score != score ) { // NaN check
		score = 0;
		std::cerr << "NaN at residue rsd: " << rsd.seqpos() << " " << rsd.name() << std::endl;
	}

	if ( score > 1e16 ) { // inf check
		std::cerr << "inf at residue rsd: " << rsd.seqpos() << " " << rsd.name() <<  " " << score << std::endl;
		score = 0;
	}


	/// sum derivatives.
	Real3 & dE_dbb(  scratch.dE_dbb() );
	Real4 & dE_dchi( scratch.dE_dchi() );

	// p0 - the base probability -- not modified by the chi-dev penalty
	Real const rotprob( scratch.rotprob() );

	Size const nbb( std::min( rsd.mainchain_torsions().size(), DUNBRACK_MAX_BBTOR) );

	Real const invp( ( rotprob == Real( 0.0 ) ) ? 0.0 : -1.0 / rotprob );

	for ( Size i=1; i<= nbb; ++i ) {
		dE_dbb[ i ] = invp * scratch.drotprob_dbb()[ i ] + scratch.dchidevpen_dbb()[ i ];
	}

	for ( Size i=1; i<= T; ++i ) {
		dE_dchi[ i ] = scratch.dchidevpen_dchi()[ i ];
	}

	correct_termini_derivatives( rsd, scratch );

	return score;
}

/// DOUG DOUG DOUG this function will need to be modified
template < Size T >
Real
RotamericSingleResiduePeptoidLibrary< T >::eval_rotameric_energy_deriv(
	conformation::Residue const & /*rsd*/,
	RotamerLibraryScratchSpace & /*scratch*/,
	bool /*eval_deriv*/
) const
{

	utility_exit_with_message( "ERROR:: rotameric single residue peptoid library::eval_rotameric_energy_deriv stubbed out!");
	return 0.0;

	/*
	assert( rsd.aa() == aa() );

	// Grab data from rsd
	Size const nbb ( rsd.mainchain_torsions().size() );
	Size const nchi( rsd.nchi() );
	ChiVector const & chi( rsd.chi() );

	assert( nbb == 3 && chi.size() == nchi );

	Real4 & chimean( scratch.chimean() );
	Real4 & chisd(   scratch.chisd()   );
	Real4 & chidev(  scratch.chidev()  );
	Real4 & chidevpen( scratch.chidevpen() );
	Real3 & drotprob_dbb( scratch.drotprob_dbb() );
	Real3 & dchidevpen_dbb( scratch.dchidevpen_dbb()  );
	Real4 & dchidevpen_dchi(scratch.dchidevpen_dchi() );
	Real4 & dchimean_dphi( scratch.dchimean_dphi()  );
	Real4 & dchimean_dpsi( scratch.dchimean_dpsi()  );
	Real4 & dchisd_dphi( scratch.dchisd_dphi()  );
	Real4 & dchisd_dpsi( scratch.dchisd_dpsi()  );

	std::fill( chimean.begin(), chimean.end(), 0.0 );
	std::fill( chisd.begin(),   chisd.end(),   0.0 );
	std::fill( chidev.begin(),  chidev.end(),  0.0 );
	std::fill( chidevpen.begin(),  chidevpen.end(),  0.0 );
	std::fill( drotprob_dbb.begin(),  drotprob_dbb.end(), 0.0 );
	std::fill( dchidevpen_dbb.begin(),  dchidevpen_dbb.end(), 0.0 );
	std::fill( dchidevpen_dchi.begin(), dchidevpen_dchi.end(), 0.0 );
	std::fill( dchimean_dphi.begin(), dchimean_dphi.end(), 0.0 );
	std::fill( dchimean_dpsi.begin(), dchimean_dpsi.end(), 0.0 );
	std::fill( dchisd_dphi.begin(), dchisd_dphi.end(), 0.0 );
	std::fill( dchisd_dpsi.begin(), dchisd_dpsi.end(), 0.0 );

	scratch.fa_dun_tot() = 0;
	scratch.fa_dun_rot() = 0;
	scratch.fa_dun_semi() = 0;
	scratch.fa_dun_dev() = 0;


	// compute rotamer number from chi
	Size4 & rotwell( scratch.rotwell() );

	/// Don't use derived class's version of this function.
	//std::cout << "RSD " << rsd.seqpos() << " ";
	RotamericSingleResiduePeptoidLibrary< T >::get_rotamer_from_chi_static( rsd.chi(), rotwell );

	Size packed_rotno( rotwell_2_packed_rotno( rotwell ));
	if ( packed_rotno == 0 ) {
		// panic!  Extremely unlikely rotamer found.  Find another rotamer that has at least some probability,
		// and move this rotamer toward it -- do so in a predictable manner so that the score function is continuous
		// as it tries to move from this rotamer to another.
		packed_rotno = find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	}

	PackedDunbrackRotamer< T, Real > interpolated_rotamer;
	interpolate_rotamers( rsd, scratch, packed_rotno, interpolated_rotamer );

	if ( dun02() ) {
		for ( Size ii = 1; ii <= T; ++ii ) { chidev[ ii ] = subtract_chi_angles( chi[ ii ], chimean[ ii ], aa(), ii ); }
	} else {
		for ( Size ii = 1; ii <= T; ++ii ) { chidev[ ii ] = basic::periodic_range( chi[ ii ] - chimean[ ii ], 360 ); }
	}

	for ( Size ii = 1; ii <= T; ++ii ) {
		/// chidev penalty: define a gaussian with a height of 1 and a standard deviation of chisd[ ii ];
		/// exp( -1 * (chi_i - chi_mean_i )**2 / ( 2 * chisd_ii**2 ) )
		/// Pretend this gaussian defines a probability for having an angle at a certain deviation from the mean.
		/// Convert that probability into an energy:
		/// -ln p = -ln exp( -1* (chi_i - chi_mean_i)**2/(2 chisd_i**2) ) = (chi_i - chi_mean_i)**2/(2 chisd_i**2)
		chidevpen[ ii ] = chidev[ ii ]*chidev[ ii ] / ( 2 * chisd[ ii ] * chisd[ ii ] );

	}
	Real chidevpensum( 0.0 );
	for ( Size ii = 1; ii <= T; ++ii ) {
		chidevpensum += chidevpen[ ii ];
	}

	scratch.fa_dun_tot() = -std::log( scratch.rotprob() ) + chidevpensum;
	scratch.fa_dun_rot() = -std::log( scratch.rotprob() );
	scratch.fa_dun_dev() = chidevpensum;


	Real const score( -std::log( scratch.rotprob() ) + chidevpensum );

	if ( ! eval_deriv ) return score;


	dchidevpen_dbb[ RotamerLibraryScratchSpace::AA_PHI_INDEX ] = 0.0;
	dchidevpen_dbb[ RotamerLibraryScratchSpace::AA_PSI_INDEX ] = 0.0;
	for ( Size ii = 1; ii <= T; ++ii ) {

		/// Backbone derivatives for chi-dev penalty.
		/// Xmean_i and sd_i both depend on phi and psi.
		/// Let: f = (X_i-Xmean_i)**2
		/// Let: g = 2 sd_i**2
		/// Then, chidevpen = f/g
		/// and, dchidevpen = (f'g - fg')/(gg)
		Real const f      = chidev[ ii ]*chidev[ ii ];
		Real const fprime = -2*chidev[ ii ];
		Real const g      = 2*chisd[ ii ]*chisd[ ii ];
		Real const gprime = 4*chisd[ ii ];
		Real const invgg  = 1 / (g*g);

		dchidevpen_dbb[ RotamerLibraryScratchSpace::AA_PHI_INDEX  ] +=
			( g*fprime*dchimean_dphi[ ii ] - f*gprime*dchisd_dphi[ ii ] ) * invgg;
		dchidevpen_dbb[ RotamerLibraryScratchSpace::AA_PSI_INDEX  ] +=
			( g*fprime*dchimean_dpsi[ ii ] - f*gprime*dchisd_dpsi[ ii ] ) * invgg;

		dchidevpen_dchi[ ii ] = chidev[ ii ] / ( chisd[ ii ] * chisd[ ii ] );
	}

	return score;

	*/
}

template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::assign_random_rotamer_with_bias(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	RotamerLibraryScratchSpace & scratch,
	numeric::random::RandomGenerator & RG,
	ChiVector & new_chi_angles,
	bool perturb_from_rotamer_center
) const
{
	Size packed_rotno( 0 );
	assign_random_rotamer( rsd, pose, scratch, RG, new_chi_angles, perturb_from_rotamer_center, packed_rotno );
}

/// DOUG DOUG DOUG this function will need to be modified
template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::assign_random_rotamer(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	RotamerLibraryScratchSpace & scratch,
	numeric::random::RandomGenerator & RG,
	ChiVector & new_chi_angles,
	bool perturb_from_rotamer_center,
	Size & packed_rotno
) const
{
	Real random_prob = RG.uniform();

	/// DOUG DOUG DOUG BROKEN
	Real const omg( get_omg_from_rsd( rsd, pose ) );
	Real const phi( get_phi_from_rsd( rsd, pose ) );
	Real const psi( get_psi_from_rsd( rsd, pose ) );

	Size omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next;
	Real omg_alpha, phi_alpha, psi_alpha;
	get_omgphipsi_bins( omg, phi, psi, omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next, omg_alpha, phi_alpha, psi_alpha );

	typename ObjexxFCL::FArray4D< PackedDunbrackRotamer< T > > const & rotamers_array( rotamers( omg ) );
	ObjexxFCL::FArray4D< Size > const & packed_rotno_2_sorted_rotno_array( packed_rotno_2_sorted_rotno( omg ) );

	PackedDunbrackRotamer< T, Real > interpolated_rotamer;

	/// Go through rotamers in decreasing order by probability and stop when the
	Size count = 0;
	while ( random_prob > 0 ) {
		packed_rotno = rotamers_array( omgbin, phibin, psibin, ++count ).packed_rotno();
		interpolate_rotamers( rotamers_array, packed_rotno_2_sorted_rotno_array,
			scratch, packed_rotno,
			omgbin, phibin, psibin,
			omgbin_next, phibin_next, psibin_next,
			omg_alpha, phi_alpha, psi_alpha,
			interpolated_rotamer );
		random_prob -= interpolated_rotamer.rotamer_probability();
		PackedDunbrackRotamer< T, Real > interpolated_rotamer;
		//loop condition might end up satisfied even if we've walked through all possible rotamers
		// if the chosen random number was nearly 1
		// (and interpolation introduced a tiny bit of numerical noise).
		if ( count == rotamers_array.size3() ) break;
	}
	assign_chi_for_interpolated_rotamer( interpolated_rotamer, rsd, RG, new_chi_angles, perturb_from_rotamer_center );

}

template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::assign_chi_for_interpolated_rotamer(
	PackedDunbrackRotamer< T, Real > const & interpolated_rotamer,
	conformation::Residue const & rsd,
	numeric::random::RandomGenerator & RG,
	ChiVector & new_chi_angles,
	bool perturb_from_rotamer_center
) const
{
	new_chi_angles.resize( rsd.nchi() );
	if ( ! perturb_from_rotamer_center ) {
		for ( Size ii = 1; ii <= T; ++ii ) {
			new_chi_angles[ ii ] = interpolated_rotamer.chi_mean( ii );
		}
	} else {
		for ( Size ii = 1; ii <= T; ++ii ) {
			new_chi_angles[ ii ] = interpolated_rotamer.chi_mean(ii) + RG.gaussian() * interpolated_rotamer.chi_sd(ii);
		}
	}

	/// Set any remaining chi uniformly? ( proton chi)
	for ( Size ii = T + 1; ii <= rsd.nchi(); ++ii ) {
		new_chi_angles[ ii ] = RG.uniform()*360.0 - 180.0;
	}

}


/// @details The new rotamer library represents only 75 of 81 possible arginine rotamers,
/// and 75 of 81 lysine rotamers. In the unlikely event that a rotamer is encountered
/// that's not represented in the library, find another rotamer to represent it.
/// Ideally, this stand-in rotamer would be closest to the input rotamer in physical geometry.
/// (sidechain atom rms, e.g.)
/// The following code instead first looks through all possible rotamers with a Hamming distance
/// of one from input rotamer trying to find one that works, looking first for rotamers from
/// the furthest chi toward the closest chi.  If no such rotamer may be found, it gives up and
/// returns the most-probable rotamer for a phi-psi bin.
///
/// This function modifies the "rotwell" assigned to this rotamer so that later code that relies
/// on the consistency of the rotwell and packed_rotno information will behave correctly.
template < Size T >
Size
RotamericSingleResiduePeptoidLibrary< T >::find_another_representative_for_unlikely_rotamer(
	conformation::Residue const & rsd,
	Size4 & rotwell
) const
{
	/// Start from the furthest chi and work inwards trying to find another rotamer that could work
	for ( Size ii = T; ii >= 1; --ii ) {
		Size ii_orig_value = rotwell[ ii ];
		for ( Size jj = 1; jj <= 3; ++jj ) { /// 3 == NUMBER OF ROTAMERIC CHI BINS
			if ( jj == ii_orig_value ) continue;
			rotwell[ ii ] = jj;
			Size new_packed_rotno = rotwell_2_packed_rotno( rotwell );
			if ( new_packed_rotno != 0 ) {
				return new_packed_rotno;
			}
		}
		/// didn't find one for chi_ii, reset rotwell to original value
		rotwell[ ii ] = ii_orig_value;
	}

	// At this point, we've examined all the rotamer wells with a Hamming distance of 1 from
	// the input rotamer well and have come up empty.
	// Just go for the most likely rotamer in this well -- no guarantee that as the rotamer swings from the
	// current rotamer well to the target rotamer well that it doesn't first swing through some other

	Size omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next;
	Real omg_alpha, phi_alpha, psi_alpha;
	get_omgphipsi_bins(
		get_omg_from_rsd( rsd ), get_phi_from_rsd( rsd ), get_psi_from_rsd( rsd ),
		omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next, omg_alpha, phi_alpha, psi_alpha );

	typename ObjexxFCL::FArray4D< PackedDunbrackRotamer< T > > const & rotamers_array( rotamers( 	get_omg_from_rsd( rsd ) ) );
	ObjexxFCL::FArray4D< Size > const & packed_rotno_2_sorted_rotno_array( packed_rotno_2_sorted_rotno( get_omg_from_rsd( rsd ) ) );

	Size packed_rotno = rotamers_array( omgbin, phibin, psibin, 1 ).packed_rotno();
	packed_rotno_2_rotwell( packed_rotno, rotwell );
	return packed_rotno;

}

/// DOUG DOUG DOUG this function will need to be modified
template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::correct_termini_derivatives(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const
{
	// mt: for the termini, these derivatives are 0, because these "pseudo"
	// mt: torsions are kept fixed.
	if ( rsd.is_lower_terminus() ) {
		scratch.dE_dbb()[ RotamerLibraryScratchSpace::AA_PHI_INDEX ] = 0;
	}
	if ( rsd.is_upper_terminus() ) {
		scratch.dE_dbb()[ RotamerLibraryScratchSpace::AA_PSI_INDEX ] = 0;
	}


}

/// @brief Returns the energy of the lowest-energy rotamer accessible to the given residue
/// (based on e.g. its current omega, phi and psi values).
/// If curr_rotamer_only is true, then consider only the idealized version of the
/// residue's current rotamer (local optimum); otherwise, consider all rotamers (global optimum).
template < Size T >
Real
RotamericSingleResiduePeptoidLibrary< T >::best_rotamer_energy(
	conformation::Residue const & /*rsd*/,
	bool /*curr_rotamer_only*/,
	RotamerLibraryScratchSpace & /*scratch*/
) const
{
	utility_exit_with_message( "ERROR:: rotameric single residue peptoid library::best_rotamer_energy stubbed out!");
	return 0.0;

	/*
	Real maxprob( 0 );
	if ( curr_rotamer_only ) {
		Size4 & rotwell( scratch.rotwell() );
		RotamericSingleResiduePeptoidLibrary< T >::get_rotamer_from_chi_static( rsd.chi(), rotwell );
		Size const packed_rotno = rotwell_2_packed_rotno( rotwell );

		PackedDunbrackRotamer< T, Real > interpolated_rotamer;
		interpolate_rotamers( rsd, pose, scratch, packed_rotno, interpolated_rotamer );

		maxprob = interpolated_rotamer.rotamer_probability();
	} else {
		Real const omg( get_omg_from_rsd( rsd, pose ) );
		Real const phi( get_phi_from_rsd( rsd, pose ) );
		Real const psi( get_psi_from_rsd( rsd, pose ) );

		Size omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next;
		Real omg_alpha, phi_alpha, psi_alpha;
		get_omgphipsi_bins( omg, phi, psi, omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next, omg_alpha, phi_alpha, psi_alpha );

		Size n_packed_rotnos(8);
		utility::vector1< Size > packed_rotnos( n_packed_rotnos, 0 );

		/// check all eight bins...
		packed_rotnos[ 1 ] = rotamers_( omgbin, phibin, psibin, 1 ).packed_rotno();
		packed_rotnos[ 2 ] = rotamers_( omgbin, phibin_next, psibin, 1 ).packed_rotno();
		packed_rotnos[ 3 ] = rotamers_( omgbin, phibin, psibin_next, 1 ).packed_rotno();
		packed_rotnos[ 4 ] = rotamers_( omgbin, phibin_next, psibin_next, 1 ).packed_rotno();
		packed_rotnos[ 5 ] = rotamers_( omgbin_next, phibin, psibin, 1 ).packed_rotno();
		packed_rotnos[ 6 ] = rotamers_( omgbin_next, phibin_next, psibin, 1 ).packed_rotno();
		packed_rotnos[ 7 ] = rotamers_( omgbin_next, phibin, psibin_next, 1 ).packed_rotno();
		packed_rotnos[ 8 ] = rotamers_( omgbin_next, phibin_next, psibin_next, 1 ).packed_rotno();

		PackedDunbrackRotamer< T, Real > interpolated_rotamer;
		for ( Size ii = 1; ii <= n_packed_rotnos; ++ii ) {
			interpolate_rotamers( rsd, pose, scratch, packed_rotnos[ ii ], interpolated_rotamer );
			maxprob = ( maxprob < interpolated_rotamer.rotamer_probability() ?
				interpolated_rotamer.rotamer_probability() : maxprob );
		}

	}

	//std::cout << "packed_rotno " << curr_rotamer_only << " " << packed_rotno <<
	//	" " << packed_rotno_2_sorted_rotno_( omgbin, phibin, psibin, 1 ).packed_rotno() <<
	//	" " << packed_rotno_2_sorted_rotno_( omgbin, phibin_next, psibin, 1 ).packed_rotno() <<
	//	" " << packed_rotno_2_sorted_rotno_( omgbin, phibin, psibin_next, 1 ).packed_rotno() <<
	//	" " << packed_rotno_2_sorted_rotno_( omgbin, phibin_next, psibin_next, 1 ).packed_rotno() <<
	//	" " << packed_rotno_2_sorted_rotno_( omgbin_next, phibin, psibin, 1 ).packed_rotno() <<
	//	" " << packed_rotno_2_sorted_rotno_( omgbin_next, phibin_next, psibin, 1 ).packed_rotno() <<
	//	" " << packed_rotno_2_sorted_rotno_( omgbin_next, phibin, psibin_next, 1 ).packed_rotno() <<
	//	" " << packed_rotno_2_sorted_rotno_( omgbin_next, phibin_next, psibin_next, 1 ).packed_rotno() << std::endl;


	return -1 * std::log( maxprob );
	*/
}

/// DOUG DOUG DOUG Update this comment
/// @details The Dunbrack library's phi/psi data for a grid point
/// (x,y) collects data in the neighborhood of (x,y).  As an interpolation
/// point p moves toward a grid point (x,y), the (x,y) share in the
/// interpolated value goes to 1.  This is distinct from having interpolation wells
/// where an interpolation in the center of a well produces the maximum contribution
/// from that well.  Most of the basic::interpolation code is designed for
/// the second interpretation of interpolation, and so CTSA's Dunbrack library code did
/// funky thinks like shift by 5 degrees so that the basic::interpolation code could
/// shift it back by 5 degrees again.  The code below makes no such shift.
///
/// The alpha fraction is the distance along each axis that the interpolation point
/// has progressed from the lower grid point toward the upper grid point; it ranges
/// from 0 to 1.
template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::get_omgphipsi_bins(
	Real omg,
	Real phi,
	Real psi,
	Size & omgbin,
	Size & phibin,
	Size & psibin,
	Size & omgbin_next,
	Size & phibin_next,
	Size & psibin_next,
	Real & omg_alpha,
	Real & phi_alpha,
	Real & psi_alpha
) const
{

	// this hideous monstrocity is required because omg is not periodic, I can't think of a better way right now
	omg = basic::periodic_range( omg, 360 );

	if ( omg >= 0.00 && omg < CIS_OMG_UPPER_RANGE ) {
		bin_angle( -30.0, OMG_BINRANGE, 70.0, N_OMG_BINS, numeric::principal_angle_degrees(omg), omgbin, omgbin_next, omg_alpha );
	} else if ( omg >= CIS_OMG_UPPER_RANGE && omg < CIS_UPPER_TRANS_LOWER_SPLIT ) { // [   30,   90 ) use cis_range_upper
		bin_angle( -30.0, OMG_BINRANGE, 70.0, N_OMG_BINS, numeric::principal_angle_degrees(CIS_OMG_UPPER_RANGE), omgbin, omgbin_next, omg_alpha );
	} else if ( omg >= CIS_UPPER_TRANS_LOWER_SPLIT && omg < TRANS_OMG_LOWER_RANGE ) { // [   90,  150 ) use trans_range_lower
		bin_angle( 150.0, OMG_BINRANGE, 70.0, N_OMG_BINS, numeric::principal_angle_degrees(TRANS_OMG_LOWER_RANGE), omgbin, omgbin_next, omg_alpha );
	} else if ( omg >= TRANS_OMG_LOWER_RANGE && omg < 180.00 ) { //	[  150, 180 ) use omg
		bin_angle( 150.0, OMG_BINRANGE, 70.0, N_OMG_BINS, numeric::nonnegative_principal_angle_degrees(omg), omgbin, omgbin_next, omg_alpha );
	} else if (omg >= -180.00 && omg < TRANS_OMG_UPPER_RANGE ) { // [  180, -150 use nnpad
		bin_angle( 150.0, OMG_BINRANGE, 70.0, N_OMG_BINS, numeric::nonnegative_principal_angle_degrees(omg), omgbin, omgbin_next, omg_alpha );
	} else if ( omg >= TRANS_OMG_UPPER_RANGE && omg < TRANS_UPPER_CIS_LOWER_SPLIT ) { // [ -150,  -90 ) use trans_range_upper
		bin_angle( 150.0, OMG_BINRANGE, 70.0, N_OMG_BINS, numeric::nonnegative_principal_angle_degrees(TRANS_OMG_UPPER_RANGE), omgbin, omgbin_next, omg_alpha );
	} else if ( omg >= TRANS_UPPER_CIS_LOWER_SPLIT && omg < CIS_OMG_LOWER_RANGE) { // [  -90,  -30 ) use cis_range_lower
		bin_angle( -30.0, OMG_BINRANGE, 70.0, N_OMG_BINS, numeric::principal_angle_degrees(CIS_OMG_LOWER_RANGE), omgbin, omgbin_next, omg_alpha );
	} else if ( omg >= CIS_OMG_LOWER_RANGE && omg < CIS_OMG_UPPER_RANGE ) { // [  -30,    30 ) use omg
		bin_angle( -30.0, OMG_BINRANGE, 70.0, N_OMG_BINS, numeric::principal_angle_degrees(omg), omgbin, omgbin_next, omg_alpha );
	}

	//std::cout << "GOPPB\t" << omg << "\t" << numeric::principal_angle_degrees(omg) << "\t" << numeric::nonnegative_principal_angle_degrees(omg) << "\t" << omgbin << "\t" << omgbin_next << "\t" << omg_alpha << std::endl;

	parent::bin_angle( -180.0, PHI_BINRANGE, 360.0, N_PHI_BINS, basic::periodic_range( phi, 360 ), phibin, phibin_next, phi_alpha );
	parent::bin_angle( -180.0, PSI_BINRANGE, 360.0, N_PSI_BINS, basic::periodic_range( psi, 360 ), psibin, psibin_next, psi_alpha );

	/// DOUG DOUG DOUG double check what this does and see if it can be removed
	/// DOUG DOUG DOUG seems like it would be more appropriate as an assert anyway
	//verify_omgphipsi_bins( omg, phi, psi, omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next );
}

template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::get_omgphipsi_bins(
	Real omg,
	Real phi,
	Real psi,
	Size & omgbin,
	Size & phibin,
	Size & psibin
) const
{
	Size omgbin_next, phibin_next, psibin_next;
	Real omg_alpha, phi_alpha, psi_alpha;
	get_omgphipsi_bins( omg, phi, psi, omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next, omg_alpha, phi_alpha, psi_alpha );
}

template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::interpolate_rotamers(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	RotamerLibraryScratchSpace & scratch,
	Size const packed_rotno,
	PackedDunbrackRotamer< T, Real > & interpolated_rotamer
) const
{
	Real omg( get_omg_from_rsd( rsd, pose ) );
	Real phi( get_phi_from_rsd( rsd, pose ) );
	Real psi( get_psi_from_rsd( rsd, pose ) );

	Size omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next;
	Real omg_alpha, phi_alpha, psi_alpha;
	get_omgphipsi_bins( omg, phi, psi, omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next, omg_alpha, phi_alpha, psi_alpha );

	typename ObjexxFCL::FArray4D< PackedDunbrackRotamer< T > > const & rotamers_array( rotamers( omg ) );
	ObjexxFCL::FArray4D< Size > const & packed_rotno_2_sorted_rotno_array( packed_rotno_2_sorted_rotno( omg ) );

	interpolate_rotamers( rotamers_array, packed_rotno_2_sorted_rotno_array,
		scratch, packed_rotno, omgbin, phibin, psibin,
		omgbin_next, phibin_next, psibin_next, omg_alpha, phi_alpha, psi_alpha,
		interpolated_rotamer );
}

/// DOUG DOUG DOUG
/// DOUG DOUG DOUG The interface has been changed but the definition has not
/// DOUG DOUG DOUG
template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::interpolate_rotamers(
	typename ObjexxFCL::FArray4D< PackedDunbrackRotamer< T > > const & rotamers,
	ObjexxFCL::FArray4D< Size > const & packed_rotno_2_sorted_rotno,
	RotamerLibraryScratchSpace & scratch,
	Size const packed_rotno,
	Size const omgbin,
	Size const phibin,
	Size const psibin,
	Size const omgbin_next,
	Size const phibin_next,
	Size const psibin_next,
	Real const omg_alpha,
	Real const phi_alpha,
	Real const psi_alpha,
	PackedDunbrackRotamer< T, Real > & interpolated_rotamer
) const
{
	using namespace basic;

	interpolated_rotamer.packed_rotno() = packed_rotno;
	Size const
		sorted_rotno_000( packed_rotno_2_sorted_rotno( omgbin,      phibin     , psibin     , packed_rotno )),
		sorted_rotno_001( packed_rotno_2_sorted_rotno( omgbin,      phibin     , psibin_next, packed_rotno )),
		sorted_rotno_010( packed_rotno_2_sorted_rotno( omgbin,      phibin_next, psibin     , packed_rotno )),
		sorted_rotno_011( packed_rotno_2_sorted_rotno( omgbin,      phibin_next, psibin_next, packed_rotno )),
		sorted_rotno_100( packed_rotno_2_sorted_rotno( omgbin_next, phibin     , psibin     , packed_rotno )),
		sorted_rotno_101( packed_rotno_2_sorted_rotno( omgbin_next, phibin     , psibin_next, packed_rotno )),
		sorted_rotno_110( packed_rotno_2_sorted_rotno( omgbin_next, phibin_next, psibin     , packed_rotno )),
		sorted_rotno_111( packed_rotno_2_sorted_rotno( omgbin_next, phibin_next, psibin_next, packed_rotno ));

	PackedDunbrackRotamer< T > const
		rot000( rotamers( omgbin,      phibin     , psibin     , sorted_rotno_000 ) ),
		rot001( rotamers( omgbin,      phibin     , psibin_next, sorted_rotno_001 ) ),
		rot010( rotamers( omgbin,      phibin_next, psibin     , sorted_rotno_010 ) ),
		rot011( rotamers( omgbin,      phibin_next, psibin_next, sorted_rotno_011 ) ),
		rot100( rotamers( omgbin_next, phibin     , psibin     , sorted_rotno_100 ) ),
		rot101( rotamers( omgbin_next, phibin     , psibin_next, sorted_rotno_101 ) ),
		rot110( rotamers( omgbin_next, phibin_next, psibin     , sorted_rotno_110 ) ),
		rot111( rotamers( omgbin_next, phibin_next, psibin_next, sorted_rotno_111 ) );

		/// DOUG DOUG DOUG check scratch INDEX stuff may need to change some things
	basic::interpolate_trilinear_by_value(
		static_cast< Real >  ( rot000.rotamer_probability()),
		static_cast< Real >  ( rot100.rotamer_probability()),
		static_cast< Real >  ( rot010.rotamer_probability()),
		static_cast< Real >  ( rot110.rotamer_probability()),
		static_cast< Real >  ( rot001.rotamer_probability()),
		static_cast< Real >  ( rot101.rotamer_probability()),
		static_cast< Real >  ( rot011.rotamer_probability()),
		static_cast< Real >  ( rot111.rotamer_probability()),
		omg_alpha, phi_alpha, psi_alpha, PHI_BINRANGE, false /*treat_as_angles*/,
		scratch.rotprob(),
		scratch.drotprob_dbb()[ RotamerLibraryScratchSpace::AA_OMG_INDEX ],
		scratch.drotprob_dbb()[ RotamerLibraryScratchSpace::AA_PHI_INDEX ],
		scratch.drotprob_dbb()[ RotamerLibraryScratchSpace::AA_PSI_INDEX ]
	);

	interpolated_rotamer.rotamer_probability() = scratch.rotprob();

	for ( Size ii = 1; ii <= T; ++ii ) {
		/// DOUG DOUG DOUG check scratch INDEX stuff may need to change some things
		basic::interpolate_trilinear_by_value(
			static_cast< Real >  ( rot000.chi_mean(ii)),
			static_cast< Real >  ( rot100.chi_mean(ii)),
			static_cast< Real >  ( rot010.chi_mean(ii)),
			static_cast< Real >  ( rot110.chi_mean(ii)),
			static_cast< Real >  ( rot001.chi_mean(ii)),
			static_cast< Real >  ( rot101.chi_mean(ii)),
			static_cast< Real >  ( rot011.chi_mean(ii)),
			static_cast< Real >  ( rot111.chi_mean(ii)),
			omg_alpha, phi_alpha, psi_alpha, PHI_BINRANGE, true /*treat_as_angles*/,
			scratch.chimean()[ii], scratch.dchimean_domg()[ii], scratch.dchimean_dphi()[ii], scratch.dchimean_dpsi()[ii] );

		// // DOUG DOUG DOUG DEBUG OUTPUT
		// std::cout << std::setprecision(2) << std::fixed
		// << "INTEPOLATE_ROTAMERS: "
		// << "CHI: " << ii  << "/" << T
		// << "\t000: " << static_cast< Real >  ( rot000.chi_mean(ii))
		// << "\t100: " << static_cast< Real >  ( rot100.chi_mean(ii))
		// << "\t010: " << static_cast< Real >  ( rot010.chi_mean(ii))
		// << "\t110: " << static_cast< Real >  ( rot110.chi_mean(ii))
		// << "\t001: " << static_cast< Real >  ( rot001.chi_mean(ii))
		// << "\t101: " << static_cast< Real >  ( rot101.chi_mean(ii))
		// << "\t011: " << static_cast< Real >  ( rot011.chi_mean(ii))
		// << "\t111: " << static_cast< Real >  ( rot111.chi_mean(ii))
		// << "\tMEAN: " << scratch.chimean()[ii]
		// << "\tOMGA: " << omg_alpha << "\tPHIA: " << phi_alpha << "\tPSIA: " << psi_alpha
		// << std::endl;

		interpolated_rotamer.chi_mean( ii ) = scratch.chimean()[ii];
		/// DOUG DOUG DOUG check scratch INDEX stuff may need to change some things
		basic::interpolate_trilinear_by_value(
			static_cast< Real >  ( rot000.chi_sd(ii)),
			static_cast< Real >  ( rot100.chi_sd(ii)),
			static_cast< Real >  ( rot010.chi_sd(ii)),
			static_cast< Real >  ( rot110.chi_sd(ii)),
			static_cast< Real >  ( rot001.chi_sd(ii)),
			static_cast< Real >  ( rot101.chi_sd(ii)),
			static_cast< Real >  ( rot011.chi_sd(ii)),
			static_cast< Real >  ( rot111.chi_sd(ii)),
			omg_alpha, phi_alpha, psi_alpha, PHI_BINRANGE, false /*treat_as_angles*/,
			scratch.chisd()[ii], scratch.dchisd_domg()[ii], scratch.dchisd_dphi()[ii], scratch.dchisd_dpsi()[ii] );
		interpolated_rotamer.chi_sd( ii ) = scratch.chisd()[ii];

	}

}

/*
	get_omg/phi/psi_from_rsd() return the omega (of the preceeding residue), phi and psi dihedral angles given a residue and a pose
	there are special cases, some of which are handled properly by the conformation::Residue class and other are not (ie. cyclic
  structures)

	default nterm
	     pomg: undefined, return NEUTRAL OMG
			 phi:  undefined, return NEUTRAL PHI
			 psi:  defined, return psi from current residue
			 omg:  defined, return omg from current residue

	acetylated nterm
	     pomg: undefined, return NEUTRAL OMG (this should change as all atoms there but don't know how to handle phi needing to be bb torsion1 when ths preceeds it)
			 phi:  defined, return phi from current residue
			 psi:  defined, return psi from current residue
			 omg:  defined, return omg from current residue

	connected end term with connected cterm in same chain (cyclic chain)
	     pomg: defined, return omg from cyclic partner (non-atomtree lookup in Residue)
			 phi:  defined, return phi from current residue (non-atomtree lookup in Residue)
			 psi:  defined, return psi from current residue
			 omg:  defined, return omg from current residue

	connected cterm with connect nterm in same chain (cyclic chain)
	     pomg: defined, return omg from preceeding residue
			 phi:  defined, return phi from current residue
			 psi:  defined, return psi from current residue (non-atomtree lookup in Residue)
			 omg:  defined, return omg from current residue (non-atomtree lookup in Residue)

	methylated cterms
	     pomg: defined, return omg from preceeding residue
			 phi:  defined, return phi from current residue
			 psi:  defined, return psi from current residue
			 omg:  defined, return omg from current residue (no problems like with pomg in actylated nterm since we are appending)

	default cterm
	     pomg: defined, return omg from preceeding residue
			 phi:  defined, return phi from current residue
			 psi:  undefined, return NEUTRAL PSI
			 omg:  undefined, return NEUTRAL OMG

 */

/// @details Returns the preceeding omg.
template < Size T >
Real
RotamericSingleResiduePeptoidLibrary< T >::get_omg_from_rsd(
	conformation::Residue const & rsd,
	pose::Pose const & pose
) const
{
	assert( rsd.is_peptoid() || rsd.is_protein() );

	if ( rsd.has_variant_type( chemical::NTERM_CONNECT ) && pose.residue( pose.conformation().chain_end( rsd.chain() ) ).has_variant_type( chemical::CTERM_CONNECT ) ) {
		assert( pose.residue( pose.conformation().chain_end( rsd.chain() ) ).is_protein() || pose.residue( pose.conformation().chain_end( rsd.chain() ) ).is_peptoid() );
		return pose.residue( pose.conformation().chain_end( rsd.chain() ) ).mainchain_torsion( RSD_OMG_INDEX );

	} else if ( rsd.is_lower_terminus() ) {
		return parent::NEUTRAL_OMG;

	}	else {
		assert( pose.residue( rsd.seqpos() - 1 ).is_protein() || pose.residue( rsd.seqpos() - 1 ).is_peptoid() );
		return pose.residue( rsd.seqpos() - 1 ).mainchain_torsion( RSD_OMG_INDEX );
	}

}

/// @details Handle lower-term residues by returning a "neutral" phi value
template < Size T >
Real
RotamericSingleResiduePeptoidLibrary< T >::get_phi_from_rsd(
	conformation::Residue const & rsd,
	pose::Pose const & pose
) const
{
	assert( rsd.is_peptoid() || rsd.is_protein() );

	if ( rsd.has_variant_type( chemical::NTERM_CONNECT ) && pose.residue( pose.conformation().chain_end( rsd.chain() ) ).has_variant_type( chemical::CTERM_CONNECT ) ) {
		assert( pose.residue( pose.conformation().chain_end( rsd.chain() ) ).is_protein() || pose.residue( pose.conformation().chain_end( rsd.chain() ) ).is_peptoid() );
		return rsd.mainchain_torsion( RSD_PHI_INDEX );

	} else if ( rsd.has_variant_type( chemical::ACETYLATED_NTERMINUS ) ) {
		return rsd.mainchain_torsion( RSD_PHI_INDEX );

	} else if ( rsd.is_lower_terminus() ) {
		return parent::NEUTRAL_PHI;

	}	else {
		return rsd.mainchain_torsion( RSD_PHI_INDEX );
	}

}

/// @details Handle upper-term residues by returning a "neutral" psi value
template < Size T >
Real
RotamericSingleResiduePeptoidLibrary< T >::get_psi_from_rsd(
	conformation::Residue const & rsd,
	pose::Pose const & pose
) const
{
	assert( rsd.is_peptoid() || rsd.is_protein() );

	if ( rsd.has_variant_type( chemical::CTERM_CONNECT ) && pose.residue( pose.conformation().chain_begin( rsd.chain() ) ).has_variant_type( chemical::NTERM_CONNECT ) ) {
		assert( pose.residue( pose.conformation().chain_begin( rsd.chain() ) ).is_protein() || pose.residue( pose.conformation().chain_begin( rsd.chain() ) ).is_peptoid() );
		return rsd.mainchain_torsion( RSD_PSI_INDEX );

	} else if ( rsd.has_variant_type( chemical::METHYLATED_CTERMINUS ) ) {
		return rsd.mainchain_torsion( RSD_PSI_INDEX );

	} else if ( rsd.is_upper_terminus() ) {
		return parent::NEUTRAL_PSI;

	}	else {
		return rsd.mainchain_torsion( RSD_PSI_INDEX );
	}

}

template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::fill_rotamer_vector(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	pack::task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const & existing_residue,
	utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
	bool buried,
	RotamerVector & rotamers
) const
{
	RotamerLibraryScratchSpace scratch;

	// DEBUG DEBUG DEBUG handel for getting type info since gdb has a hard time with RT-CAPS
	core::chemical::ResidueType existing_type( *concrete_residue );

	/// Save backbone interpolation data for reuse
	Real omg( get_omg_from_rsd( existing_residue, pose ) );
	Real phi( get_phi_from_rsd( existing_residue, pose ) );
	Real psi( get_psi_from_rsd( existing_residue, pose ) );
	Size omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next;
	Real omg_alpha, phi_alpha, psi_alpha;
	get_omgphipsi_bins( omg, phi, psi, omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next, omg_alpha, phi_alpha, psi_alpha );

	typename ObjexxFCL::FArray4D< PackedDunbrackRotamer< T > > const & rotamers_array( RotamericSingleResiduePeptoidLibrary::rotamers( omg ) );
	ObjexxFCL::FArray4D< Size > const & packed_rotno_2_sorted_rotno_array( packed_rotno_2_sorted_rotno( omg ) );

	Real const requisit_probability = probability_to_accumulate_while_building_rotamers( buried ); // ( buried  ? 0.98 : 0.95 )
	Real accumulated_probability( 0.0 );

	Size const max_rots_that_can_be_built = n_packed_rots();
	Size count_rotamers_built = 0;

	// /// DOUG DOUG DOUG DEBUG
	// std::cout << "OMG: " << omg << "/"<< omgbin << "/" << omgbin_next << "/" << omg_alpha << "\t"
	//           << "PHI: " << phi << "/"<< phibin << "/" << phibin_next << "/" << phi_alpha << "\t"
	//           << "PSI: " << psi << "/"<< psibin << "/" << psibin_next << "/" << psi_alpha << "\t"
	//           << "RP: " << requisit_probability << "\tMR: " << max_rots_that_can_be_built << std::endl;

	while ( accumulated_probability < requisit_probability ) {
		// Iterate through rotamaers in decreasing order of probabilities, stopping once the requisit probility is hit.
		++count_rotamers_built;

		Size const packed_rotno00 = rotamers_array( omgbin, phibin, psibin, count_rotamers_built ).packed_rotno();

		PackedDunbrackRotamer< T, Real > interpolated_rotamer;
		interpolate_rotamers(
			rotamers_array, packed_rotno_2_sorted_rotno_array,
			scratch, packed_rotno00,
			omgbin, phibin, psibin,
			omgbin_next, phibin_next, psibin_next,
			omg_alpha, phi_alpha, psi_alpha,
			interpolated_rotamer );

		build_rotamers(
			pose, scorefxn, task, packer_neighbor_graph,
			concrete_residue, existing_residue, extra_chi_steps, buried, rotamers,
			interpolated_rotamer );

		accumulated_probability += interpolated_rotamer.rotamer_probability();
		if ( count_rotamers_built == max_rots_that_can_be_built ) break; // this shouldn't happen...
	}

	// /// DOUG DOUG DOUG
	// std::cout << "---" << std::endl;
	// for ( Size i(1); i <= rotamers.size(); ++i ) {
	// 	std::cout << "ROTAMER " << i << ":\t" << rotamers[i]->chi(1) << "\t" << rotamers[i]->chi(2) << "\t" << rotamers[i]->chi(3) << std::endl;
	// }

}

template < Size T >
utility::vector1< DunbrackRotamerSampleData >
RotamericSingleResiduePeptoidLibrary< T >::get_all_rotamer_samples(
	Real omg,
	Real phi,
	Real psi
) const
{
	RotamerLibraryScratchSpace scratch;

	Size omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next;
	Real omg_alpha, phi_alpha, psi_alpha;
	get_omgphipsi_bins( omg, phi, psi, omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next, omg_alpha, phi_alpha, psi_alpha );

	typename ObjexxFCL::FArray4D< PackedDunbrackRotamer< T > > const & rotamers_array( rotamers( omg ) );
	ObjexxFCL::FArray4D< Size > const & packed_rotno_2_sorted_rotno_array( packed_rotno_2_sorted_rotno( omg ) );

	Size const n_rots = n_packed_rots();
	utility::vector1< DunbrackRotamerSampleData > all_rots;
	all_rots.reserve( n_rots );

	for ( Size ii = 1; ii <= n_rots; ++ii ) {
		// Iterate through rotamaers in decreasing order of probabilities
		Size const packed_rotno00 = rotamers_array( omgbin, phibin, psibin, ii ).packed_rotno();

		PackedDunbrackRotamer< T, Real > interpolated_rotamer;
		interpolate_rotamers(
			rotamers_array, packed_rotno_2_sorted_rotno_array,
			scratch, packed_rotno00,
			omgbin, phibin, psibin,
			omgbin_next, phibin_next, psibin_next,
			omg_alpha, phi_alpha, psi_alpha,
			interpolated_rotamer );

		DunbrackRotamerSampleData sample( false );
		sample.set_nchi( T );
		sample.set_rotwell( parent::packed_rotno_2_rotwell( interpolated_rotamer.packed_rotno() ));
		for ( Size jj = 1; jj <= T; ++jj ) sample.set_chi_mean( jj, interpolated_rotamer.chi_mean( jj ) );
		for ( Size jj = 1; jj <= T; ++jj ) sample.set_chi_sd( jj, interpolated_rotamer.chi_sd( jj ) );
		sample.set_prob( interpolated_rotamer.rotamer_probability() );

		all_rots.push_back( sample );
	}

	return all_rots;
}

/// DOUG DOUG DOUG The interface of this function has beeen changed but the definition needs more work specifically  the interpolation
template < Size T >
Real
RotamericSingleResiduePeptoidLibrary< T >::get_probability_for_rotamer(
	Real omg,
	Real phi,
	Real psi,
	Size rot_ind
) const
{
	Size omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next;
	Real omg_alpha, phi_alpha, psi_alpha;
	get_omgphipsi_bins( omg, phi, psi, omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next, omg_alpha, phi_alpha, psi_alpha );

	typename ObjexxFCL::FArray4D< PackedDunbrackRotamer< T > > const & rotamers_array( rotamers( omg ) );
	ObjexxFCL::FArray4D< Size > const & packed_rotno_2_sorted_rotno_array( packed_rotno_2_sorted_rotno( omg ) );

	PackedDunbrackRotamer< T > const & rot000( rotamers_array( omgbin, phibin, psibin, rot_ind ) );
	Size const packed_rotno000 = rot000.packed_rotno();

	Size const
		sorted_rotno_001( packed_rotno_2_sorted_rotno_array( omgbin,      phibin     , psibin_next, packed_rotno000 )),
		sorted_rotno_010( packed_rotno_2_sorted_rotno_array( omgbin,      phibin_next, psibin     , packed_rotno000 )),
		sorted_rotno_011( packed_rotno_2_sorted_rotno_array( omgbin,      phibin_next, psibin_next, packed_rotno000 )),
		sorted_rotno_100( packed_rotno_2_sorted_rotno_array( omgbin_next, phibin     , psibin_next, packed_rotno000 )),
		sorted_rotno_101( packed_rotno_2_sorted_rotno_array( omgbin_next, phibin     , psibin_next, packed_rotno000 )),
		sorted_rotno_110( packed_rotno_2_sorted_rotno_array( omgbin_next, phibin_next, psibin     , packed_rotno000 )),
		sorted_rotno_111( packed_rotno_2_sorted_rotno_array( omgbin_next, phibin_next, psibin_next, packed_rotno000 ));

	PackedDunbrackRotamer< T > const &
		rot001( rotamers_array( omgbin,      phibin     , psibin_next, sorted_rotno_001 ) ),
		rot010( rotamers_array( omgbin,      phibin_next, psibin     , sorted_rotno_010 ) ),
		rot011( rotamers_array( omgbin,      phibin_next, psibin_next, sorted_rotno_011 ) ),
		rot100( rotamers_array( omgbin_next, phibin_next, psibin_next, sorted_rotno_100 ) ),
		rot101( rotamers_array( omgbin_next, phibin     , psibin_next, sorted_rotno_101 ) ),
		rot110( rotamers_array( omgbin_next, phibin_next, psibin     , sorted_rotno_110 ) ),
		rot111( rotamers_array( omgbin_next, phibin_next, psibin_next, sorted_rotno_111 ) );

	Real rot_prob, dummy1, dummy2, dummy3;

	basic::interpolate_trilinear_by_value(
		static_cast< Real >  ( rot000.rotamer_probability()),
		static_cast< Real >  ( rot100.rotamer_probability()),
		static_cast< Real >  ( rot010.rotamer_probability()),
		static_cast< Real >  ( rot110.rotamer_probability()),
		static_cast< Real >  ( rot001.rotamer_probability()),
		static_cast< Real >  ( rot101.rotamer_probability()),
		static_cast< Real >  ( rot011.rotamer_probability()),
		static_cast< Real >  ( rot111.rotamer_probability()),
		omg_alpha, phi_alpha, psi_alpha, PHI_BINRANGE, false /*treat_as_angles*/,
		rot_prob, dummy1, dummy2, dummy3
	);
	return rot_prob;
}

template < Size T >
DunbrackRotamerSampleData
RotamericSingleResiduePeptoidLibrary< T >::get_rotamer(
	Real omg,
	Real phi,
	Real psi,
	Size rot_ind
) const
{
	RotamerLibraryScratchSpace scratch;

	Size omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next;
	Real omg_alpha, phi_alpha, psi_alpha;
	get_omgphipsi_bins( omg, phi, psi, omgbin, phibin, psibin, omgbin_next, phibin_next, psibin_next, omg_alpha, phi_alpha, psi_alpha );

	typename ObjexxFCL::FArray4D< PackedDunbrackRotamer< T > > const & rotamers_array( rotamers( omg ) );
	ObjexxFCL::FArray4D< Size > const & packed_rotno_2_sorted_rotno_array( packed_rotno_2_sorted_rotno( omg ) );

	PackedDunbrackRotamer< T > const & rot00( rotamers_array( omgbin, phibin, psibin, rot_ind ) );
	Size const packed_rotno00 = rot00.packed_rotno();
	PackedDunbrackRotamer< T, Real > interpolated_rotamer;
	interpolate_rotamers(
		rotamers_array, packed_rotno_2_sorted_rotno_array,
		scratch, packed_rotno00,
		omgbin, phibin, psibin,
		omgbin_next, phibin_next, psibin_next,
		omg_alpha, phi_alpha, psi_alpha,
		interpolated_rotamer );

	DunbrackRotamerSampleData sample( false );
	sample.set_nchi( T );
	sample.set_rotwell( parent::packed_rotno_2_rotwell( interpolated_rotamer.packed_rotno() ));
	for ( Size jj = 1; jj <= T; ++jj ) sample.set_chi_mean( jj, interpolated_rotamer.chi_mean( jj ) );
	for ( Size jj = 1; jj <= T; ++jj ) sample.set_chi_sd( jj, interpolated_rotamer.chi_sd( jj ) );
	sample.set_prob( interpolated_rotamer.rotamer_probability() );

	return sample;
}


template < Size T >
Size
RotamericSingleResiduePeptoidLibrary< T >::nchi() const
{
	return T;
}

template < Size T >
Size
RotamericSingleResiduePeptoidLibrary< T >::n_rotamer_bins() const
{
	return parent::n_possible_rots();
}


/// @details Load interpolated rotamer data into a RotamericData object on the stack
/// so that it can be handed into the enumerate_chi_sets method, which itself relies
/// on the "template method" chisamples_for_rotamer_and_chi.  Enumerate the chi samples,
/// and build rotamers from these chi samples.
///
/// "template method" is a design pattern where a base class calls a polymorphic method
/// that can be overloaded by a derived class, usually in the middle of a function that
/// does a lot of work.  See "Design Patterns," Gamma et al.
template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::build_rotamers(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	pack::task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const& existing_residue,
	utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
	bool buried,
	RotamerVector & rotamers,
	PackedDunbrackRotamer< T, Real > const & interpolated_rotamer
) const
{
	DunbrackRotamer< T, Real > interpolated_rot( packed_rotamer_2_regular_rotamer( interpolated_rotamer ));
	RotamericData< T > rotameric_rotamer_building_data( interpolated_rot );


	/// DOUG DOUG DOUG DEBUG OUTPUT
	/*
	std::cout << "EXTRA_CHI_STEPS::build_rotamers\t" << extra_chi_steps.size() << std::endl;
	for ( Size i(1); i <= extra_chi_steps.size(); ++i ) {
			for ( Size j(1); j <= extra_chi_steps[i].size(); ++j ) {
				std::cout << i << "/" << j << ":\t" << extra_chi_steps[i][j] << "\t" << std::flush;
			}
			std::cout << std::endl;
	}
	std::cout << std::endl;
	*/

	// now build the chi sets derived from this base rotamer
	utility::vector1< ChiSetOP > chi_set_vector;
	enumerate_chi_sets(
		*concrete_residue, task, existing_residue.seqpos(), buried,
		rotameric_rotamer_building_data,
		extra_chi_steps, chi_set_vector );

	create_rotamers_from_chisets(
		pose, scorefxn, task,
		packer_neighbor_graph, concrete_residue, existing_residue,
		chi_set_vector, rotamers );
}

template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::verify_omgphipsi_bins(
	Real omg,
	Real phi,
	Real psi,
	Size const omgbin,
	Size const phibin,
	Size const psibin,
	Size const omgbin_next,
	Size const phibin_next,
	Size const psibin_next
) const
{
	if ( (omgbin < 1 || omgbin > N_OMG_BINS ) || ( phibin < 1 || phibin > N_PHI_BINS ) || ( psibin < 1 || psibin > N_PHI_BINS ) ||
		( omgbin_next < 1 || omgbin_next > N_OMG_BINS) || ( phibin_next < 1 || phibin_next > N_PHI_BINS ) || ( psibin_next < 1 || psibin_next > N_PSI_BINS )) {
		std::cerr << "ERROR: omg/phi/psi bin out of range: " << aa() << " " << omg << " " << phi << " " << psi << " " << omgbin << " " << omgbin_next << " " << phibin << " " << phibin_next << " " << psibin << " " << psibin_next <<  std::endl;
		utility_exit();
	}
}

/// @details once a list of chi samples has been enumerated, this function
/// instantiates Residue objectes and give them the correct geometry.
template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::create_rotamers_from_chisets(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	pack::task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const& existing_residue,
	utility::vector1< ChiSetOP > const & chi_set_vector,
	RotamerVector & rotamers
) const
{
	using namespace utility;
	pack::task::ResidueLevelTask const & rtask( task.residue_task( existing_residue.seqpos() ) );

	// construct real rotamers
	for ( vector1< ChiSetOP >::const_iterator chi_set( chi_set_vector.begin() );
	      chi_set != chi_set_vector.end(); ++chi_set ) {
		conformation::ResidueOP rotamer = conformation::ResidueFactory::create_residue(
			*concrete_residue, existing_residue, pose.conformation(), rtask.preserve_c_beta() );
		for ( Size jj = 1; jj <= (*chi_set)->chi.size(); ++jj ) {
			rotamer->set_chi( jj, (*chi_set)->chi[ jj ] );
		}
		// apply an operation (or a filter) to this rotamer at build time
		bool reject(false);
		for ( pack::rotamer_set::RotamerOperations::const_iterator
			    op( rtask.rotamer_operations().begin() );
		      op != rtask.rotamer_operations().end(); ++op ) {
			reject |= ! (**op)( rotamer, pose, scorefxn, rtask, packer_neighbor_graph, *chi_set );
		}
		if ( !reject ) rotamers.push_back( rotamer );
	}
}

template< Size T >
template< class P >
DunbrackRotamer< T, P >
RotamericSingleResiduePeptoidLibrary< T >::packed_rotamer_2_regular_rotamer(
	PackedDunbrackRotamer< T, P > const & packedrot
) const
{
	DunbrackRotamer< T, P > dunrot;
	for ( Size ii = 1; ii <= T; ++ii ) {
		dunrot.chi_mean( ii ) = packedrot.chi_mean( ii );
		dunrot.chi_sd(   ii ) = packedrot.chi_sd( ii );
		dunrot.rotwell(  ii ) = packed_rotno_2_rotwell( packedrot.packed_rotno() )[ ii ];
	}
	dunrot.rotamer_probability() = packedrot.rotamer_probability();
	return dunrot;
}


template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::enumerate_chi_sets(
	chemical::ResidueType const & rsd_type,
	pack::task::PackerTask const & task,
	Size const seqpos,
	bool buried,
	RotamericData< T > const & rotamer_data,
	utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
	utility::vector1< ChiSetOP > & chi_set_vector
) const
{
	Real const base_probability( rotamer_data.rotamer().rotamer_probability() );

	using namespace utility;

	Size const nchi( rsd_type.nchi() );

	vector1< vector1< Real > > total_chi( nchi );
	vector1< vector1< int  > > total_rot( nchi );
	vector1< vector1< Real > > total_ex_steps( nchi ); // remember which extra chi steps were used for each rotamer
	vector1< vector1< Real > > chisample_prob( nchi ); // let derived class define sample probability

	pack::task::ResidueLevelTask const & rtask( task.residue_task( seqpos ) );
	/// List all chi samples for each chi
	for ( Size ii = 1; ii <= nchi; ++ii ) {
		chisamples_for_rotamer_and_chi(
			rsd_type, rtask, buried, ii, rotamer_data, extra_chi_steps[ ii ],
			total_chi[ ii ], total_rot[ ii ], total_ex_steps[ ii ], chisample_prob[ ii ]
		);
	}

	// now convert from total_chi to chisets....

	Size exchi_product = 1;
	utility::vector1< Size > lex_sizes( total_chi.size() );
	for ( Size ii = 1; ii <= nchi; ++ii ) {
		lex_sizes[ ii ] = total_chi[ii].size();
		exchi_product *= total_chi[ii].size();
	}

	chi_set_vector.reserve( exchi_product );

	// previously named min_extrachi_rot_prob in rosetta++
	Real const minimum_probability_for_extra_rotamers = 0.001; // make this a virtual function call...

	bool first_rotamer( true );
	for ( utility::LexicographicalIterator lex( lex_sizes ); ! lex.at_end(); ++lex ) {
		Real chi_set_prob = base_probability;

    runtime_assert( nchi <=  lex.size() );
    runtime_assert( nchi <=  chisample_prob.size() );

		for ( Size ii = 1; ii <= nchi; ++ii ) {
			if( chisample_prob[ ii ].size() >= lex[ ii ] ){
				chi_set_prob *= chisample_prob[ ii ][ lex[ ii ] ];
			}
		}
		if ( ! first_rotamer && chi_set_prob < minimum_probability_for_extra_rotamers ) {
			continue;
		}
		first_rotamer = false;

		ChiSetOP chi_set = new ChiSet( nchi );
		for ( Size ii = 1; ii <= nchi; ++ii ) {
			chi_set->chi[ ii ] = total_chi[ ii ][ lex[ ii ] ];
			chi_set->ex_chi_steps[ ii ] = total_ex_steps[ ii ][ lex[ ii ] ];
			chi_set->rot[ ii ] = total_rot[ ii ][ lex[ ii ] ];
		}
		chi_set_vector.push_back( chi_set );
	}
}



template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::chisamples_for_rotamer_and_chi(
	chemical::ResidueType const & rsd_type,
	pack::task::ResidueLevelTask const & rtask,
	bool buried,
	Size const chi_index,
	RotamericData< T > const & rotamer_data,
	utility::vector1< Real > const & extra_steps,
	utility::vector1< Real > & total_chi,
	utility::vector1< int  > & total_rot,
	utility::vector1< Real > & total_ex_steps,
	utility::vector1< Real > & chisample_prob
) const
{
	//assert( dynamic_cast< RotamericData const & > ( rotamer_data ) );
	//RotamericData const & rotameric_data( static_cast< RotamericData const & > ( rotamer_data ) );

	// setting this value to zero disables the small-standdev filter -- is that what we want?
	DunbrackReal const min_extrachi_sd( 0.0 );


	if ( chi_index > T ) {
		// chi angle not present in dunbrack library, eg. ser,thr,tyr Hchi

		utility::vector1< std::pair< Real, Real > > const & chi_rotamers( rsd_type.chi_rotamers( chi_index ) );
		if ( ! chi_rotamers.empty() ) {
			// use the rotamer info encoded in the residue
			for ( Size j=1; j<= chi_rotamers.size(); ++j ) {
				Real const rot_chi_mean( chi_rotamers[j].first );
				Real const rot_chi_sdev( chi_rotamers[j].second );
				total_chi.push_back( rot_chi_mean );
				total_ex_steps.push_back( 0. );
				total_rot.push_back( j );
				for ( Size k=1; k<= extra_steps.size(); ++k ) {
					total_chi.push_back( rot_chi_mean + extra_steps[k] * rot_chi_sdev );
					total_ex_steps.push_back( extra_steps[k] );
					total_rot.push_back( j );
					chisample_prob.push_back( 1.0 ); // assume perfectly likely rotamer sample?
				}
			}
		} else if ( rsd_type.is_proton_chi( chi_index ) ) {
			pack::task::ExtraRotSample ex_samp_level = rtask.extrachi_sample_level( buried, chi_index, &rsd_type );

			//std::cout << "chi_index: " << chi_index << " ex_samp_level " << ex_samp_level << std::endl;
			bool const skip_extra_proton_chi_samples = ( ex_samp_level == pack::task::NO_EXTRA_CHI_SAMPLES );

			Size const proton_chi_index = rsd_type.chi_2_proton_chi( chi_index );
			utility::vector1< Real > const & samples = rsd_type.proton_chi_samples( proton_chi_index );
			utility::vector1< Real > const & extra_samples = rsd_type.proton_chi_extra_samples( proton_chi_index );
			for ( Size ii = 1; ii <= samples.size(); ++ii ) {
				total_chi.push_back( samples[ ii ] );
				total_ex_steps.push_back( 0.0 );
				total_rot.push_back( 1 );
				chisample_prob.push_back( 1.0 );

				if ( skip_extra_proton_chi_samples ) continue;

				for ( Size jj = 1; jj <= extra_samples.size(); ++jj ) {
					total_chi.push_back( samples[ ii ] + extra_samples[ jj ] );
					total_ex_steps.push_back( 0.0 );
					total_rot.push_back( 1 );
					chisample_prob.push_back( 1.0 );

					total_chi.push_back( samples[ ii ] - extra_samples[ jj ] );
					total_ex_steps.push_back( 0.0 );
					total_rot.push_back( 1 );
					chisample_prob.push_back( 1.0 );
				}
			}

		} else {
			// just use the chi of the platonic residue -- this could be made faster, but chi is not part of rsdtype intrfc
			// Is this code ever executed?  It's not going to be executed for any amino acid...
			Real const icoor_chi( numeric::dihedral(
				rsd_type.atom( rsd_type.chi_atoms( chi_index )[1] ).ideal_xyz(),
				rsd_type.atom( rsd_type.chi_atoms( chi_index )[2] ).ideal_xyz(),
				rsd_type.atom( rsd_type.chi_atoms( chi_index )[3] ).ideal_xyz(),
				rsd_type.atom( rsd_type.chi_atoms( chi_index )[4] ).ideal_xyz() ) );
			total_chi.push_back( icoor_chi );
			total_ex_steps.push_back( 0. );
			total_rot.push_back( 1 );
			chisample_prob.push_back( 1.0 );
		}
	} else {
		// use the dunbrack data

		total_chi.push_back( rotamer_data.rotamer().chi_mean( chi_index ) );
		total_ex_steps.push_back( 0. );
		total_rot.push_back( rotamer_data.rotamer().rotwell( chi_index ) );
		chisample_prob.push_back( 1.0 );

		// ctsa - eliminate extra chi angles with insignificant sd's
		if ( rotamer_data.rotamer().chi_sd( chi_index ) <= min_extrachi_sd ) return;

		for ( Size k=1; k<= extra_steps.size(); ++k ) {
			total_chi.push_back( rotamer_data.rotamer().chi_mean( chi_index )
				+ extra_steps[k] * rotamer_data.rotamer().chi_sd( chi_index ) );
			total_ex_steps.push_back( extra_steps[k] );
			total_rot.push_back( rotamer_data.rotamer().rotwell( chi_index ) );
			/// prob = exp( -( x - xmean )**2 / 2sd**2). Since we're sampling at intrvals of sd,
			/// prob = exp( -( step * sd )**2 / 2sd**2 ) = exp( - step**2/2 ).
			/// of course, a true gaussian probability would be scaled by 1/sd(2p)**0.5...
			chisample_prob.push_back( std::exp( -0.5 * (extra_steps[ k ]*extra_steps[ k ])) );

		}
	}
}

// Basically write_to_binary but to a file instead
template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::write_to_file( utility::io::ozstream & /*out*/ ) const
{
	using namespace boost;
	/*
	parent::write_to_binary( out );

	/// 1. rotamers_
	{
	Size const ntotalrot = N_PHI_BINS * N_PSI_BINS * parent::n_packed_rots();
	Size const ntotalchi = ntotalrot * T;
	///a. means
	DunbrackReal * rotamer_means = new DunbrackReal[ ntotalchi ];
	// b. standard deviations
	DunbrackReal * rotamer_stdvs = new DunbrackReal[ ntotalchi ];
	// c. rotamer probabilities
	DunbrackReal * rotamer_probs = new DunbrackReal[ ntotalrot ];
	// d. packed rotamer numbers
	DunbrackReal * packed_rotnos = new DunbrackReal[ ntotalrot ];
	Size count_chi( 0 ), count_rots( 0 );
	for ( Size ii = 1; ii <= parent::n_packed_rots(); ++ii ) {
		for ( Size jj = 1; jj <= N_PHI_BINS; ++jj ) {
			for ( Size kk = 1; kk <= N_PSI_BINS; ++kk ) {
				for ( Size ll = 1; ll <= T; ++ll ) {
					rotamer_means[ count_chi ] = rotamers_( kk, jj, ii ).chi_mean( ll );
					rotamer_stdvs[ count_chi ] = rotamers_( kk, jj, ii ).chi_sd( ll );
					++count_chi;
				}
				rotamer_probs[ count_rots ] = rotamers_( kk, jj, ii ).rotamer_probability();
				packed_rotnos[ count_rots ] = rotamers_( kk, jj, ii ).packed_rotno();
				++count_rots;
			}
		}
	}
	out.write( (char*) rotamer_means, ntotalchi * sizeof( DunbrackReal ));
	out.write( (char*) rotamer_stdvs, ntotalchi * sizeof( DunbrackReal ));
	out.write( (char*) rotamer_probs, ntotalrot * sizeof( DunbrackReal ));
	out.write( (char*) packed_rotnos, ntotalrot * sizeof( DunbrackReal ));
	delete [] rotamer_means;
	delete [] rotamer_stdvs;
	delete [] rotamer_probs;
	delete [] packed_rotnos;
	}

	/// 2. packed_rotno_to_sorted_rotno_
	{
	Size const ntotalpackedrots = N_PHI_BINS * N_PSI_BINS * parent::n_packed_rots();
	boost::int32_t * packed_rotno_2_sorted_rotno = new boost::int32_t[ ntotalpackedrots ];
	Size count( 0 );
	for ( Size ii = 1; ii <= parent::n_packed_rots(); ++ii ) {
		for ( Size jj = 1; jj <= N_PHI_BINS; ++jj ) {
			for ( Size kk = 1; kk <= N_PSI_BINS; ++kk ) {
				packed_rotno_2_sorted_rotno[ count ] = packed_rotno_2_sorted_rotno_( kk, jj, ii );
				++count;
			}
		}
	}
	out.write( (char*) packed_rotno_2_sorted_rotno, ntotalpackedrots * sizeof( boost::int32_t ) );
	delete [] packed_rotno_2_sorted_rotno;
	}

	*/
}

/// DOUG DOUG DOUG This function needs more updating to use omega
template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::write_to_binary( utility::io::ozstream & /*out*/ ) const
{
	std::cout << "Called write_to_binary() but it is not yet implimented!!! " << std::endl;
	/* DOUG DOUG DOUG Commenting out for now

	using namespace boost;

	parent::write_to_binary( out );

	/// 1. rotamers_
	{
	Size const ntotalrot = N_PHI_BINS * N_PSI_BINS * parent::n_packed_rots();
	Size const ntotalchi = ntotalrot * T;
	///a. means
	DunbrackReal * rotamer_means = new DunbrackReal[ ntotalchi ];
	// b. standard deviations
	DunbrackReal * rotamer_stdvs = new DunbrackReal[ ntotalchi ];
	// c. rotamer probabilities
	DunbrackReal * rotamer_probs = new DunbrackReal[ ntotalrot ];
	// d. packed rotamer numbers
	DunbrackReal * packed_rotnos = new DunbrackReal[ ntotalrot ];
	Size count_chi( 0 ), count_rots( 0 );
	for ( Size ii = 1; ii <= parent::n_packed_rots(); ++ii ) {
		for ( Size jj = 1; jj <= N_PHI_BINS; ++jj ) {
			for ( Size kk = 1; kk <= N_PSI_BINS; ++kk ) {
				for ( Size ll = 1; ll <= T; ++ll ) {
					rotamer_means[ count_chi ] = rotamers_( kk, jj, ii ).chi_mean( ll );
					rotamer_stdvs[ count_chi ] = rotamers_( kk, jj, ii ).chi_sd( ll );
					++count_chi;
				}
				rotamer_probs[ count_rots ] = rotamers_( kk, jj, ii ).rotamer_probability();
				packed_rotnos[ count_rots ] = rotamers_( kk, jj, ii ).packed_rotno();
				++count_rots;
			}
		}
	}
	out.write( (char*) rotamer_means, ntotalchi * sizeof( DunbrackReal ));
	out.write( (char*) rotamer_stdvs, ntotalchi * sizeof( DunbrackReal ));
	out.write( (char*) rotamer_probs, ntotalrot * sizeof( DunbrackReal ));
	out.write( (char*) packed_rotnos, ntotalrot * sizeof( DunbrackReal ));
	delete [] rotamer_means;
	delete [] rotamer_stdvs;
	delete [] rotamer_probs;
	delete [] packed_rotnos;
	}

	/// 2. packed_rotno_to_sorted_rotno_
	{
	Size const ntotalpackedrots = N_PHI_BINS * N_PSI_BINS * parent::n_packed_rots();
	boost::int32_t * packed_rotno_2_sorted_rotno = new boost::int32_t[ ntotalpackedrots ];
	Size count( 0 );
	for ( Size ii = 1; ii <= parent::n_packed_rots(); ++ii ) {
		for ( Size jj = 1; jj <= N_PHI_BINS; ++jj ) {
			for ( Size kk = 1; kk <= N_PSI_BINS; ++kk ) {
				packed_rotno_2_sorted_rotno[ count ] = packed_rotno_2_sorted_rotno_( kk, jj, ii );
				++count;
			}
		}
	}
	out.write( (char*) packed_rotno_2_sorted_rotno, ntotalpackedrots * sizeof( boost::int32_t ) );
	delete [] packed_rotno_2_sorted_rotno;
	}

	*/ 
}

/// DOUG DOUG DOUG This function needs more updating to use omega
template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::read_from_binary( utility::io::izstream & /*in*/ )
{
	std::cout << "Called write_to_binary() but it is not yet implimented!!! " << std::endl;
	/* DOUG DOUG DOUG Commenting out for now

	parent::read_from_binary( in );
	/// 1. rotamers_
	{
	Size const ntotalrot = N_PHI_BINS * N_PSI_BINS * parent::n_packed_rots();
	Size const ntotalchi = ntotalrot * T;
	rotamers_.dimension( N_PHI_BINS, N_PSI_BINS, parent::n_packed_rots() );
	///a. means
	DunbrackReal * rotamer_means = new DunbrackReal[ ntotalchi ];
	// b. standard deviations
	DunbrackReal * rotamer_stdvs = new DunbrackReal[ ntotalchi ];
	// c. rotamer probabilities
	DunbrackReal * rotamer_probs = new DunbrackReal[ ntotalrot ];
	// d. packed rotamer numbers
	DunbrackReal * packed_rotnos = new DunbrackReal[ ntotalrot ];

	in.read( (char*) rotamer_means, ntotalchi * sizeof( DunbrackReal ));
	in.read( (char*) rotamer_stdvs, ntotalchi * sizeof( DunbrackReal ));
	in.read( (char*) rotamer_probs, ntotalrot * sizeof( DunbrackReal ));
	in.read( (char*) packed_rotnos, ntotalrot * sizeof( DunbrackReal ));

	Size count_chi( 0 ), count_rots( 0 );
	for ( Size ii = 1; ii <= parent::n_packed_rots(); ++ii ) {
		for ( Size jj = 1; jj <= N_PHI_BINS; ++jj ) {
			for ( Size kk = 1; kk <= N_PSI_BINS; ++kk ) {
				for ( Size ll = 1; ll <= T; ++ll ) {
					rotamers_( kk, jj, ii ).chi_mean( ll ) = rotamer_means[ count_chi ];
					rotamers_( kk, jj, ii ).chi_sd(  ll ) = rotamer_stdvs[ count_chi ];
					++count_chi;
				}
				rotamers_( kk, jj, ii ).rotamer_probability() = rotamer_probs[ count_rots ];
				rotamers_( kk, jj, ii ).packed_rotno() = ( Size ) packed_rotnos[ count_rots ];
				++count_rots;
			}
		}
	}
	delete [] rotamer_means;
	delete [] rotamer_stdvs;
	delete [] rotamer_probs;
	delete [] packed_rotnos;
	}

	/// 2. packed_rotno_to_sorted_rotno_
	{
	Size const ntotalpackedrots = N_PHI_BINS * N_PSI_BINS * parent::n_packed_rots();
	boost::int32_t * packed_rotno_2_sorted_rotno = new boost::int32_t[ ntotalpackedrots ];
	in.read( (char*) packed_rotno_2_sorted_rotno, ntotalpackedrots * sizeof( boost::int32_t ) );
	Size count( 0 );
	packed_rotno_2_sorted_rotno_.dimension( N_PHI_BINS, N_PSI_BINS, parent::n_packed_rots() );
	for ( Size ii = 1; ii <= parent::n_packed_rots(); ++ii ) {
		for ( Size jj = 1; jj <= N_PHI_BINS; ++jj ) {
			for ( Size kk = 1; kk <= N_PSI_BINS; ++kk ) {
				packed_rotno_2_sorted_rotno_( kk, jj, ii ) = packed_rotno_2_sorted_rotno[ count ];
				++count;
			}
		}
	}
	delete [] packed_rotno_2_sorted_rotno;
	}
	*/
}

/// DOUG DOUG DOUG WORKING ON THIS FUNCTION
/// @details Returns the three letter string of the next amino acid specified in the
/// input library.
template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::read_from_file(
	utility::io::izstream & infile
)
{
	std::string const my_name( chemical::name_from_aa( aa() ) );
	std::string next_name; // empty string to start

	/// Read all the data from the first phi/psi bin and keep it temporarily.
	/// Note all the rotamer wells encountered along the way; then declare_all_rotwells_encountered,
	/// allocate the big tables, transfer the data in the temporary arrays into the big table,
	typename utility::vector1< DunbrackRotamer< T > > first_omgphipsibin_data;
	first_omgphipsibin_data.reserve( n_possible_rots() );

	DunbrackReal omg(0.0), phi(0.0), psi(0.0), probability(0.0);
	Size omgbin(0), phibin(0), psibin(0), lastomgbin(0), lastphibin(0), lastpsibin(0), count(0);
	bool very_first_rotamer( true );
	bool finished_first_omgphipsi_bin( false );
	Size count_in_this_omgphipsi_bin( 1 );
	utility::vector1< Size > rotwell( DUNBRACK_MAX_SCTOR, 0 );
	utility::vector1< DunbrackReal > chimean( DUNBRACK_MAX_SCTOR, 0.0 );
	utility::vector1< DunbrackReal > chisd( DUNBRACK_MAX_SCTOR, 0.0 );
	std::string three_letter_code;

	// read the whole file
	utility::vector1< std::string > lines;
	std::string line;
	while ( getline( infile, line ) ) {
		// if it starts with #, skip to the next line.
		if ( line[ 0 ] == '#' ) continue;
		else lines.push_back( line );
	}

	// parse the lines
	for( Size i( 1 ); i <= lines.size(); ++i ) {
		/// 1. Read the line.  Format is:
		/// a.   three-letter-code,
		/// b.   omg
		/// c.   phi,
		/// d.   psi,
		/// e.   count
		/// f..i r1..r4
		/// j.   probability
		/// k..n chimean1..chimean4
		/// o..r chisd1..chisd4

		std::istringstream l( lines[ i ] );

		l >> three_letter_code;
		l >> omg >> phi >> psi >> count;
		l >> rotwell[ 1 ] >> rotwell[ 2 ] >> rotwell[ 3 ] >> rotwell[ 4 ];
		l >> probability;
		l >> chimean[ 1 ] >> chimean[ 2 ] >> chimean[ 3 ] >> chimean[ 4 ];
		l >> chisd[ 1 ] >> chisd[ 2 ] >> chisd[ 3 ] >> chisd[ 4 ];

		if ( omg == -180 || phi == -180 || psi == -180 ) continue; // duplicated data...

		/// AVOID INF!
		for ( Size ii = 1; ii <= T; ++ii ) {
			if ( chisd[ ii ] == 0.0 ) {
				chisd[ ii ] = 5; // bogus!
			}
		}
		if ( probability == 0.0 ) {
			probability = 1e-4;
			/// APL -- On the advice of Roland Dunbrack, modifying the minimum probability to the
			/// resolution of the library.  This helps avoid overwhelmingly unfavorable energies
			/// (5 log-units difference between 1e-4 and 1e-9) for rare rotamers.
		}

		get_omgphipsi_bins( omg, phi, psi, omgbin, phibin, psibin );

		//std::cout << "RFF:\t" << omg << "\t" << phi << "\t" << psi << "\t" << omgbin << "\t" << phibin << "\t" << psibin << std::endl;


		if ( finished_first_omgphipsi_bin ) {
			Size const packed_rotno = rotwell_2_packed_rotno( rotwell );
			if ( packed_rotno < 1 || packed_rotno > parent::n_packed_rots() ) {
				std::cerr << "ERROR in converting rotwell to packed rotno: ";
				for ( Size ii = 1; ii <= T; ++ii ) std::cerr << " " << rotwell[ ii ];
				std::cerr << " " << packed_rotno << std::endl;
				utility_exit();
			}
			PackedDunbrackRotamer< T > rotamer( chimean, chisd, probability, packed_rotno );

			if ( omgbin != lastomgbin || phibin != lastphibin || psibin != lastpsibin ) {
				count_in_this_omgphipsi_bin = 1;
			}

			// decide to add to the cis or trans rotamer or packed_rotno_2_sorted_rotno arrays
			if ( is_trans_omg( omg ) ) {
				trans_rotamers_( omgbin, phibin, psibin, count_in_this_omgphipsi_bin ) = rotamer;
				trans_packed_rotno_2_sorted_rotno_( omgbin, phibin, psibin, packed_rotno ) = count_in_this_omgphipsi_bin;
			} else if ( is_cis_omg( omg ) ) {
				cis_rotamers_( omgbin, phibin, psibin, count_in_this_omgphipsi_bin ) = rotamer;
				cis_packed_rotno_2_sorted_rotno_( omgbin, phibin, psibin, packed_rotno ) = count_in_this_omgphipsi_bin;
			}

			++count_in_this_omgphipsi_bin;
			lastomgbin = omgbin; lastphibin = phibin; lastpsibin = psibin;
		} else {
			if ( !very_first_rotamer && (omgbin != lastomgbin || phibin != lastphibin || psibin != lastpsibin )) {
				// We have now read all rotamers from this phi/psi bin.
				// 1. Declare all existing rotwells encountered
				// 2. Allocate space for rotamer data
				// 3. Store data from first rotwell.
				// 4. Store the data from the rotamer we just read.

				// 1.
				declare_all_existing_rotwells_encountered();
				// 2.

				// decide if cis or trans rotamer or packed_rotno_2_sorted_rotno arrays
				//if ( is_trans_omg( omg ) ) {
					trans_rotamers_.dimension( N_OMG_BINS, N_PHI_BINS, N_PSI_BINS, n_packed_rots() );
					trans_packed_rotno_2_sorted_rotno_.dimension( N_OMG_BINS, N_PHI_BINS, N_PSI_BINS, n_packed_rots() );
					//} else if ( is_cis_omg( omg ) ) {
					cis_rotamers_.dimension( N_OMG_BINS, N_PHI_BINS, N_PSI_BINS, n_packed_rots() );
					cis_packed_rotno_2_sorted_rotno_.dimension( N_OMG_BINS, N_PHI_BINS, N_PSI_BINS, n_packed_rots() );
					//}

				// 3.
				utility::vector1< Size > first_bin_rotwell( DUNBRACK_MAX_SCTOR, 0 );
				for ( Size ii = 1; ii <= first_omgphipsibin_data.size(); ++ii ) {
					for ( Size jj = 1; jj <= T; ++jj ) {
						first_bin_rotwell[ jj ] = first_omgphipsibin_data[ ii ].rotwell( jj );
					}
					Size const packed_rotno = rotwell_2_packed_rotno( first_bin_rotwell );

					if ( packed_rotno < 1 || packed_rotno > parent::n_packed_rots() ) {
						std::cerr << "ERROR in converting rotwell to packed rotno: ";
						for ( Size ii = 1; ii <= T; ++ii ) std::cerr << " " << rotwell[ ii ];
						std::cerr << " " << packed_rotno << std::endl;
						utility_exit();
					}

					// decide if cis or trans rotamer or packed_rotno_2_sorted_rotno arrays
					if ( is_trans_omg( omg ) ) {
						trans_rotamers_( lastomgbin, lastphibin, lastpsibin, ii ) = PackedDunbrackRotamer< T >( first_omgphipsibin_data[ ii ], packed_rotno );
						trans_packed_rotno_2_sorted_rotno_( lastomgbin, lastphibin, lastpsibin, packed_rotno ) = ii;
					} else if ( is_cis_omg( omg ) ) {
						cis_rotamers_( lastomgbin, lastphibin, lastpsibin, ii ) = PackedDunbrackRotamer< T >( first_omgphipsibin_data[ ii ], packed_rotno );
						cis_packed_rotno_2_sorted_rotno_( lastomgbin, lastphibin, lastpsibin, packed_rotno ) = ii;
					}

				}

				// 4.
				assert( count_in_this_omgphipsi_bin == 1 );
				Size const packed_rotno = rotwell_2_packed_rotno( rotwell );
				PackedDunbrackRotamer< T > rotamer( chimean, chisd, probability, packed_rotno );
				if ( is_trans_omg( omg ) ) {
					trans_rotamers_( omgbin, phibin, psibin, count_in_this_omgphipsi_bin ) = rotamer;
					trans_packed_rotno_2_sorted_rotno_( omgbin, phibin, psibin, packed_rotno ) = count_in_this_omgphipsi_bin;
				} else if ( is_cis_omg( omg ) ) {
					cis_rotamers_( omgbin, phibin, psibin, count_in_this_omgphipsi_bin ) = rotamer;
					cis_packed_rotno_2_sorted_rotno_( omgbin, phibin, psibin, packed_rotno ) = count_in_this_omgphipsi_bin;
				}

				++count_in_this_omgphipsi_bin;
				lastomgbin = omgbin; lastphibin = phibin; lastpsibin = psibin;

				finished_first_omgphipsi_bin = true;
			} else {
				very_first_rotamer = false;
				mark_rotwell_exists( rotwell );
				first_omgphipsibin_data.push_back( DunbrackRotamer< T >( chimean, chisd, probability, rotwell ) );
				lastomgbin = omgbin; lastphibin = phibin; lastpsibin = psibin;
			}
		}
	}
}

/// @details called only if the library is actually an RSRDL<T> object.  Derived classes
/// should not call this function or recurse.  Accounts for the statically allocated data
/// that's part of this class.
template < Size T >
Size RotamericSingleResiduePeptoidLibrary< T >::memory_usage_static() const
{
	return sizeof( RotamericSingleResiduePeptoidLibrary< T > );
}

/// @details Measures the amount of dynamically allocated data in this class.  Must recurse to parent
/// to count parent's dynamically allocated data.
template < Size T >
Size RotamericSingleResiduePeptoidLibrary< T >::memory_usage_dynamic() const
{
	Size total = parent::memory_usage_dynamic(); // recurse to parent.
	total += trans_rotamers_.size() * sizeof( PackedDunbrackRotamer< T > );
	total += cis_rotamers_.size() * sizeof( PackedDunbrackRotamer< T > );
	total += trans_packed_rotno_2_sorted_rotno_.size() * sizeof( Size ); // could make these shorts or chars!
	total += cis_packed_rotno_2_sorted_rotno_.size() * sizeof( Size ); // could make these shorts or chars!
	//total += max_rotprob_.size() * sizeof( DunbrackReal );

	// DOUG DOUG DOUG DEBUG
	//std::cout << trans_rotamers_.size() * sizeof( PackedDunbrackRotamer< T > ) << "\t"
	//<< cis_rotamers_.size() * sizeof( PackedDunbrackRotamer< T > ) << "\t"
	//<< trans_packed_rotno_2_sorted_rotno_.size() * sizeof( Size ) << "\t"
	//<< cis_packed_rotno_2_sorted_rotno_.size() * sizeof( Size ) << std::endl;

	return total;
}

template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::get_rotamer_from_chi(
	ChiVector const & chi,
	RotVector & rot
) const
{
	Size4 rot4;
	get_rotamer_from_chi_static( chi, rot4 );
	rot.resize( chi.size() );
	for ( Size ii = 1; ii <= T; ++ii ) rot[ ii ] = rot4[ ii ];
	for ( Size ii = T+1; ii <= chi.size(); ++ii ) rot[ ii ] = 0;
}

template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::get_rotamer_from_chi_static(
	ChiVector const & chi,
	Size4 & rot
) const
{
	Real4 chi4;
	for ( Size ii = 1; ii <= T; ++ii ) chi4[ ii ] = chi[ ii ];
	get_rotamer_from_chi_static( chi4, rot );
}

template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::get_rotamer_from_chi_static(
	Real4 const & chi,
	Size4 & rot
) const
{
	if ( dun02() ) { rotamer_from_chi_02( chi, aa(), T, rot ); return; }

	assert( chi.size() >= T );

	/// compiler will unroll this loop
	for ( Size ii = 1; ii <= T; ++ii ) {
		 rot[ ii ]  = bin_rotameric_chi( basic::periodic_range( chi[ ii ], 360 ), ii );
	}
}

template < Size T >
typename ObjexxFCL::FArray4D< PackedDunbrackRotamer< T > > const &
RotamericSingleResiduePeptoidLibrary< T >::rotamers( Real omg_angle ) const {
	if ( is_trans_omg( omg_angle) ) return trans_rotamers_;
	else if ( is_cis_omg( omg_angle ) ) return cis_rotamers_;
	utility_exit_with_message("ERROR: Peptoid backbone is not in range of rotamer library");
	return trans_rotamers_; // we should never get here
}

template < Size T >
typename ObjexxFCL::FArray4D< PackedDunbrackRotamer< T > > &
RotamericSingleResiduePeptoidLibrary< T >::rotamers( Real omg_angle ) {
	if ( is_trans_omg( omg_angle) ) return trans_rotamers_;
	else if ( is_cis_omg( omg_angle ) ) return cis_rotamers_;
	utility_exit_with_message("ERROR: Peptoid backbone is not in range of rotamer library");
	return trans_rotamers_; // we should never get here
}

template < Size T >
ObjexxFCL::FArray4D< Size > const &
RotamericSingleResiduePeptoidLibrary< T >::packed_rotno_2_sorted_rotno( Real omg_angle ) const {
	if ( is_trans_omg( omg_angle) ) return trans_packed_rotno_2_sorted_rotno_;
	else if ( is_cis_omg( omg_angle ) ) return cis_packed_rotno_2_sorted_rotno_;
	utility_exit_with_message("ERROR: Peptoid backbone is not in range of rotamer library");
	return trans_packed_rotno_2_sorted_rotno_; // we should never get here
}

template < Size T >
ObjexxFCL::FArray4D< Size > &
RotamericSingleResiduePeptoidLibrary< T >::packed_rotno_2_sorted_rotno( Real omg_angle ) {
	if ( is_trans_omg( omg_angle) ) return trans_packed_rotno_2_sorted_rotno_;
	else if ( is_cis_omg( omg_angle ) ) return cis_packed_rotno_2_sorted_rotno_;
	utility_exit_with_message("ERROR: Peptoid backbone is not in range of rotamer library");
	return trans_packed_rotno_2_sorted_rotno_; // we should never get here
}

} // namespace dunbrack
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_pack_dunbrack_RotamericSingleResiduePeptoidLibrary_TMPL_HH
