// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_filter/SidechainDepthFilter.cc
/// @brief Filter for measuring minimum distance from any sidechain atom in a reside
//  to the closest water molecule
/// @author Hahnbeom Park

// Unit headers
#include <protocols/simple_filters/ResidueDepthFilter.hh>
#include <protocols/simple_filters/ResidueDepthFilterCreator.hh>

#include <protocols/rosetta_scripts/util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/database/open.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/CrystInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AA.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>

#include <core/types.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <numeric/random/random.hh>
#include <numeric/model_quality/rms.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>
#include <cmath>

namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.ResidueDepthFilter" );

bool
mycomp( const std::pair< core::Real, core::Size >& lhs,
	const std::pair< core::Real, core::Size >& rhs )
{ return lhs.first < rhs.first; }

protocols::filters::FilterOP
ResidueDepthFilterCreator::create_filter() const { return protocols::filters::FilterOP( new ResidueDepthFilter ); }

std::string
ResidueDepthFilterCreator::keyname() const { return "ResidueDepth"; }

////////////////////////////////////////////////////////////////////////

ResidueDepthFrag::ResidueDepthFrag()
{
	rdds_.resize( 9 );
}
ResidueDepthFrag::~ResidueDepthFrag()= default;

utility::vector1< Vector >
ResidueDepthFrag::get_CAcrd() const {
	utility::vector1< Vector > r_1( 9 );
	for ( core::Size k = 1; k <= 9; ++k ) r_1[k] = rdds_[k]->CAcrd;
	return r_1;
}

utility::vector1< Vector >
ResidueDepthFrag::get_CENcrd() const {
	utility::vector1< Vector > r_1( 9 );
	for ( core::Size k = 1; k <= 9; ++k ) r_1[k] = rdds_[k]->CENcrd;
	return r_1;
}

////////////////////////////////////////////////////////////////////////
ResidueDepthCalculator::ResidueDepthCalculator( core::pose::Pose const &pose )
{
	initialize( pose );
}

ResidueDepthCalculator::~ResidueDepthCalculator() = default;

void
ResidueDepthCalculator::initialize( core::pose::Pose const &pose )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	nres_ = pose.total_residue();
	niter_ = 25; // # independent runs
	//niter_ = 1; // # independent runs
	sc_depth_avrg_.resize( nres_, 0.0 );
	sc_depth_sdev_.resize( nres_, 0.0 );
	sc_depth_fvar_.resize( nres_, 0.0 );
	//waterbox_file_ = +"/protocol_data/spc216.pdb";
	waterbox_file_ = option[ in::path::database ](1).name() +"/protocol_data/SDE/spc216.pdb";
	report_crd_ = false;
	use_bb_ = true;
	use_sc_ = true;
	dcut1_ = 2.6;
	dcut2_ = 4.2;
}

// water crds are fixed throughout the protocol
utility::vector1< Vector >
ResidueDepthCalculator::read_unit_waterbox( Vector &boxwidth ) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::chemical::ResidueTypeSetCOP rsd_set =
		core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");

	TR << "water file: " << waterbox_file_ << std::endl;
	core::pose::Pose unit_waterbox_pose;
	core::import_pose::pose_from_file( unit_waterbox_pose, *rsd_set, waterbox_file_,
		core::import_pose::PDB_file );
	utility::vector1< Vector > unit_waterbox_crd;

	// Let's force it to be a cube
	core::io::CrystInfo ci = unit_waterbox_pose.pdb_info()->crystinfo();
	//runtime_assert(ci.A()*ci.B()*ci.C() != 0 && ci.A()==ci.B() && ci.A()==ci.C());
	//boxwidth[0] = ci.A();  boxwidth[1] = ci.B(); boxwidth[2] = ci.C();

	boxwidth[0] = 18.621; boxwidth[1] = 18.621; boxwidth[2] = 18.621;

	for ( core::Size iwat = 1; iwat <= unit_waterbox_pose.total_residue(); ++iwat ) {
		if ( unit_waterbox_pose.residue(iwat).is_virtual_residue() ||
				unit_waterbox_pose.residue(iwat).aa() != core::chemical::aa_h2o ) continue;
		unit_waterbox_crd.push_back( unit_waterbox_pose.residue(iwat).xyz(" O  ") );
	}
	TR << "boxwidth? " << boxwidth[0] << " " << boxwidth[1] << " " << boxwidth[2] << std::endl;

	return unit_waterbox_crd;
}

void
ResidueDepthCalculator::get_pose_crd_and_index( core::pose::Pose const &pose,
	utility::vector1< Vector > &protein_crd,
	utility::vector1< core::Size > &coarse_index,
	utility::vector1< core::Size > &res_id
) const
{
	core::Size icrd( 0 );

	for ( core::Size ires = 1; ires <= nres_; ++ires ) {
		if ( pose.residue( ires ).is_virtual_residue() ) continue;
		//TR << "ires: " << ires << " " << pose.residue(ires).last_backbone_atom()
		//  << " " << pose.residue(ires).first_sidechain_atom() << std::endl;

		// backbone
		if ( use_bb_ ) {
			for ( core::Size iatm = 1; iatm <= pose.residue(ires).last_backbone_atom(); ++iatm ) {
				icrd++;
				protein_crd.push_back( pose.residue(ires).xyz(iatm) );
				res_id.push_back( ires );
				// coarse
				//TR << "protein: icrd/ires/iatm " << protein_crd.size() << " " << ires << " " << iatm << " " << pose.residue(ires).xyz(iatm)[0] << std::endl;
				if ( pose.residue(ires).atom_name(iatm).compare(" CA ") == 0 ) {
					coarse_index.push_back( icrd );
				}
			}

		}

		// sidechain
		if ( use_sc_ ) {
			// special case when only using GLY "sidechain"
			if ( pose.residue( ires ).aa() == core::chemical::aa_gly && !use_bb_ ) {
				icrd++;
				protein_crd.push_back( pose.residue(ires).xyz(" CA ") );
				coarse_index.push_back( icrd );
				res_id.push_back( ires );
				//TR << "protein: icrd/ires/iatm " << protein_crd.size() << " " << ires << " " << "GLYCA" << " " << pose.residue(ires).xyz(" CA ")[0] << std::endl;
			} else {
				for ( core::Size iatm = pose.residue(ires).first_sidechain_atom();
						iatm <= pose.residue(ires).nheavyatoms(); ++iatm ) {
					icrd++;
					protein_crd.push_back( pose.residue(ires).xyz(iatm) );
					res_id.push_back( ires );
					//TR << "protein: icrd/ires/iatm " << protein_crd.size() << " " << ires << " " << iatm << " " << pose.residue(ires).xyz(iatm)[0] << std::endl;
					//if( pose.residue(ires).atom_name(iatm).compare(" CB ") == 0 )
				}
				//coarse
				coarse_index.push_back( pose.residue(ires).nbr_atom() );
			}

		}
	}
}

utility::vector1< Vector >
ResidueDepthCalculator::get_coarse_crd( utility::vector1< Vector > const &protein_crd,
	utility::vector1< core::Size > const &coarse_index
) const
{
	utility::vector1< Vector > coarse_crd( coarse_index.size() );
	for ( core::Size i = 1; i <= coarse_index.size(); ++i ) {
		coarse_crd[i] = protein_crd[coarse_index[i]];
	}
	return coarse_crd;
}

void
ResidueDepthCalculator::duplicate_waterbox( utility::vector1< Vector > const &unit_waterbox_crd,
	Vector const boxwidth,
	Vector const mincrds,
	Vector const maxcrds
) const
{
	// first center protein to origin
	//Vector maxcrds( 0.0 ), mincrds( 0.0 );

	utility::vector0< int > ndupl_min( 3 ), ndupl_max( 3 );
	// give 6 Ang room at each edge...
	for ( core::Size i = 0; i < 3; ++i ) {
		//ndupl_min[i] = (int)((mincrds[i]+boxwidth[i]*0.5-8.0)/boxwidth[i]);
		ndupl_min[i] = (int)((mincrds[i]-boxwidth[i]*0.5-6.0)/boxwidth[i]);
		if ( ndupl_min[i] > 0 ) ndupl_min[i] = 0;
		//ndupl_max[i] = (int)((maxcrds[i]-boxwidth[i]*0.5)/boxwidth[i]);
		ndupl_max[i] = (int)((maxcrds[i]+boxwidth[i]*0.5+6.0)/boxwidth[i]);
		if ( ndupl_max[i]< 0 ) ndupl_max[i] = 0;
	}

	TR << "x: " << mincrds[0] << " " << maxcrds[0] << std::endl;
	TR << "y: " << mincrds[1] << " " << maxcrds[1] << std::endl;
	TR << "z: " << mincrds[2] << " " << maxcrds[2] << std::endl;

	TR << "x: " << ndupl_min[0] << " " << ndupl_max[0] << std::endl;
	TR << "y: " << ndupl_min[1] << " " << ndupl_max[1] << std::endl;
	TR << "z: " << ndupl_min[2] << " " << ndupl_max[2] << std::endl;

	utility::vector1< Vector > waterbox_dupl;
	//duplicate unitbox
	for ( int i = ndupl_min[0]; i <= ndupl_max[0]; ++i ) {
		for ( int j = ndupl_min[1]; j <= ndupl_max[1]; ++j ) {
			for ( int k = ndupl_min[2]; k <= ndupl_max[2]; ++k ) {
				append_unitbox( unit_waterbox_crd, waterbox_dupl, boxwidth, i, j, k );
			}
		}
	}
	waterbox_ = waterbox_dupl;
}

Vector
ResidueDepthCalculator::bring_to_origin( utility::vector1< Vector > &protein_crd,
	Vector &maxcrds,
	Vector &mincrds
) const
{
	Vector com( 0.0 );
	for ( core::Size iatm = 1; iatm <= protein_crd.size(); ++iatm ) {
		Vector const &crd = protein_crd[iatm];
		com += protein_crd[iatm];
		for ( core::Size k = 0; k < 3; ++k ) {
			if ( crd[k] < mincrds[k] ) mincrds[k] = crd[k];
			if ( crd[k] > maxcrds[k] ) maxcrds[k] = crd[k];
		}
	}

	for ( core::Size k = 0; k < 3; ++k ) com[k] /= core::Real(protein_crd.size());

	for ( core::Size iatm = 1; iatm <= protein_crd.size(); ++iatm ) protein_crd[iatm] -= com;
	maxcrds -= com; mincrds -= com;

	return com;
}

void
ResidueDepthCalculator::append_unitbox( utility::vector1< Vector > const &unitbox,
	utility::vector1< Vector > &waterbox_dupl,
	Vector const boxwidth,
	int const i, int const j, int const k ) const
{
	//TR << "append: " << unitbox.size() << " " << i << " " << j << " " << k << std::endl;
	Vector dcom( 0.0 );
	dcom[0] = i*boxwidth[0];
	dcom[1] = j*boxwidth[1];
	dcom[2] = k*boxwidth[2];

	for ( core::Size iwat = 1; iwat <= unitbox.size(); ++iwat ) {
		waterbox_dupl.push_back( unitbox[iwat] + dcom );
		//TR << "water " << waterbox_dupl.size() << " " << unitbox[iwat][0] << " " << dcom[0] << std::endl;
	}
}

void
ResidueDepthCalculator::pert_protein( utility::vector1< Vector > &protein_crd ) const
{
	core::Real const MAX_TRANS( 2.8 );
	core::Real const MAX_ROT( 30.0 );
	core::Real const deg2rad( 57.29577951308232 );

	// temporarily bring to origin
	Vector maxcrds( 0.0 ), mincrds( 0.0 );
	Vector dcom = bring_to_origin( protein_crd, maxcrds, mincrds );

	// translation, only to xaxis (as in paper, but why?)
	core::Real dtrans = MAX_TRANS*(2.0*numeric::random::rg().uniform() - 1.0);
	dcom[0] += dtrans;

	// rotation
	// use Quaternion to pert randomly
	// angle
	core::Real pertangle = MAX_ROT*(2.0*numeric::random::rg().uniform() - 1.0 )*deg2rad;
	core::Real quatw = cos(pertangle*0.5);

	// get normalized axis vector
	Vector quat3;
	for ( core::Size k = 0; k < 3; ++k ) {
		core::Real val = (2.0*numeric::random::rg().uniform() - 1.0 );
		quat3[k] = val;
	}
	quat3 *= sin(0.5*pertangle)/std::sqrt(quat3.length_squared());

	TR << "pert: " << dtrans << " " << pertangle/deg2rad << ", quat: " << quatw << " " << quat3[0] << " " << quat3[1] << " " << quat3[2] << std::endl;

	// generate rotation matrix from quaternion
	utility::vector1< Vector > U = quat2U( quat3, quatw );

	// rotation operation
	for ( core::Size iatm = 1; iatm <=  protein_crd.size(); ++ iatm ) {
		Vector const crd( protein_crd[iatm] );
		//quaternion_operator( U, crd );
		for ( core::Size k = 0; k < 3; ++k ) {
			protein_crd[iatm][k] = crd[0]*U[k+1][0] + crd[1]*U[k+1][1] + crd[2]*U[k+1][2];
		}

		if ( report_crd_ ) {
			printf("ATOM  %5d  CA  ALA %5d    %8.3f%8.3f%8.3f 1.00  0.00 C\n",
				int(iatm), int(iatm), protein_crd[iatm][0], protein_crd[iatm][1], protein_crd[iatm][2] );
		}
	}

	if ( report_crd_ ) printf("TER\n");

	// translation
	for ( core::Size iatm = 1; iatm <= protein_crd.size(); ++ iatm ) {
		protein_crd[iatm] += dcom;
	}
}

//convert quaternion to rotation matrix
utility::vector1< Vector >
ResidueDepthCalculator::quat2U( Vector const & quat3, core::Real const &quatw ) const
{
	utility::vector1< Vector > U( 3 );

	core::Real q00, q01, q02, q03, q11, q12, q13, q22, q23, q33;
	q00 = 2.0*quatw*quatw - 1.0;
	q01 = 2.0*quatw*quat3[0];
	q02 = 2.0*quatw*quat3[1];
	q03 = 2.0*quatw*quat3[2];
	q11 = 2.0*quat3[0]*quat3[0]; q12 = 2.0*quat3[0]*quat3[1]; q13 = 2.0*quat3[0]*quat3[2];
	q22 = 2.0*quat3[1]*quat3[1]; q23 = 2.0*quat3[1]*quat3[2]; q33 = 2.0*quat3[2]*quat3[2];

	U[1][0] = q00 + q11; U[1][1] = q12 - q03; U[1][2] = q13 + q02;
	U[2][0] = q12 + q03; U[2][1] = q00 + q22; U[2][2] = q23 - q01;
	U[3][0] = q13 - q02; U[3][1] = q23 + q01; U[3][2] = q00 + q33;

	return U;
}

utility::vector1< bool >
ResidueDepthCalculator::get_exclusion_index( utility::vector1< Vector > const & protein_crd,
	utility::vector1< Vector > const & protein_coarse_crd
) const
{
	core::Real const D2_COARSE( 400.0 ); // 20 Ang
	core::Real const D2_CUT1( dcut1_*dcut1_ );
	core::Real const D2_CUT2( dcut2_*dcut2_ );
	core::Real const D2_CO( 12.0*12.0 ); // 12 Ang

	core::Size nexcl( 0 );
	utility::vector1< bool > exclude_wats( waterbox_.size(), false );

	// first filter too far away ones through coarse atoms
	for ( core::Size iwat = 1; iwat <= waterbox_.size(); ++iwat ) {
		bool found_close( false );
		for ( core::Size iatm = 1; iatm <= protein_coarse_crd.size(); ++iatm ) {
			core::Real d2( waterbox_[iwat].distance_squared( protein_coarse_crd[iatm] ) );
			//TR << iwat << " " << iatm << " " << d2 << std::endl;
			if ( d2 < D2_COARSE ) {
				found_close = true;
				break;
			}
		}

		if ( !found_close ) {
			nexcl++;
			exclude_wats[ iwat ] = true;
			//TR << "remove " << iwat << " by distance, " << nexcl << "/" << waterbox_.size() << std::endl;
		}
	}
	TR << "remaining waters after coarse filtering: " << waterbox_.size() - nexcl << std::endl;

	utility::vector1< utility::vector1< core::Size > > contact_list( protein_crd.size() );

	// then go through close enough ones
	for ( core::Size iwat = 1; iwat <= waterbox_.size(); ++iwat ) {
		if ( exclude_wats[ iwat ] ) continue;

		core::Size ncut2( 0 );
		bool exclude( false );

		for ( core::Size iatm = 1; iatm <= protein_crd.size(); ++iatm ) {
			core::Real d2 = waterbox_[iwat].distance_squared( protein_crd[iatm] );
			if ( d2 < D2_CUT1 ) { // 2.6A
				exclude = true;
				nexcl ++;
				//TR << "exclude " << iwat << " by " << iatm << ", dcut1 " << std::sqrt(d2) << ", " << nexcl << "/" << waterbox_.size() << std::endl;
				break;
			}
			if ( d2 < D2_CUT2 ) ncut2 ++; // 4.2A

			if ( ncut2 >= 2 ) {
				nexcl ++;
				exclude = true;
				//TR << "exclude " << iwat << " by dcut2, " << nexcl << "/" << waterbox_.size() << std::endl;
				break;
			}

			//
			if ( d2 < D2_CO ) contact_list[iatm].push_back( iwat );
		}

		if ( exclude ) exclude_wats[ iwat ] = true;
	}

	TR << "remaining waters after fine filtering: " << waterbox_.size() - nexcl << std::endl;

	if ( report_crd_ ) {
		// report
		for ( core::Size iwat = 1; iwat <= waterbox_.size(); ++iwat ) {
			if ( exclude_wats[ iwat ] ) continue;
			printf("HETATM%5d  O   HOH %5d    %8.3f%8.3f%8.3f 1.00  0.00 O\n",
				int(iwat), int(iwat),
				waterbox_[iwat][0], waterbox_[iwat][1], waterbox_[iwat][2]);
		}
	}

	return exclude_wats;
}

utility::vector1< core::Real >
ResidueDepthCalculator::get_scdepth( utility::vector1< bool > const &excluded_wat,
	utility::vector1< Vector > const & protein_crd,
	utility::vector1< core::Size > const & res_id
) const
{
	utility::vector1< core::Real > scdepth( nres_, 999.0 );

	for ( core::Size iatm = 1; iatm <= protein_crd.size(); ++iatm ) {
		Vector const &crd = protein_crd[iatm];
		core::Size const ires( res_id[iatm] );
		for ( core::Size iwat = 1; iwat <= excluded_wat.size(); ++iwat ) {
			if ( excluded_wat[iwat] ) continue;
			core::Real d = waterbox_[iwat].distance( crd );
			if ( d < scdepth[ires] ) scdepth[ires] = d;
		}
	}

	return scdepth;
}

bool
ResidueDepthCalculator::stack_and_getaverage(
	utility::vector1< utility::vector1< core::Real > > &sc_depth_stack,
	utility::vector1< core::Real > const &sc_depth,
	core::Size const niter
) const
{
	// validate
	core::Size nerr( 0 );
	if ( niter > 1 ) {
		for ( core::Size ires = 1; ires <= sc_depth.size(); ++ires ) {
			if ( std::abs(sc_depth[ires] - sc_depth_avrg_[ires]) > 10.0 ) nerr++;
			if ( nerr > 1 ) break;
		}
		if ( nerr > 1 ) {
			TR << "Skip iter: " << sc_depth_stack.size() << std::endl;
			return false;
		}
	}

	// then stack
	sc_depth_stack.push_back( sc_depth );

	bool converged = true;
	core::Size const nstack( sc_depth_stack.size() );

	TR << "stack! " << nstack << std::endl;

	// get averaged value & delta to current
	for ( core::Size ires = 1; ires <= sc_depth.size(); ++ires ) {
		core::Real avrg( 0.0 ), sdev( 0.0 ), fvar( 0.0 );

		for ( core::Size istack = 1; istack <= nstack; ++istack ) {
			avrg += sc_depth_stack[istack][ires];
		}
		avrg /= core::Real(nstack);

		for ( core::Size istack = 1; istack <= nstack; ++istack ) {
			sdev += (sc_depth_stack[istack][ires] - avrg)*(sc_depth_stack[istack][ires] - avrg);
		}
		sdev /= core::Real(nstack);
		sdev = std::sqrt(sdev);

		fvar = sdev/avrg;

		sc_depth_avrg_[ires] = avrg;
		sc_depth_sdev_[ires] = sdev;
		sc_depth_fvar_[ires] = fvar;

		if ( fvar > 0.25 ) converged = false;
	}

	if ( niter < 2 ) {
		return false;
	} else {
		return converged;
	}
}

// this is main
utility::vector1< core::Real >
ResidueDepthCalculator::estimate_sidechain_depth( core::pose::Pose const &pose ) const
{
	TR << "start!" << std::endl;

	utility::vector1< Vector > unit_waterbox_crd, protein_crd;
	utility::vector1< core::Size > coarse_index, res_id;

	Vector boxwidth( 0.0 );

	TR << "read water" << std::endl;
	unit_waterbox_crd = read_unit_waterbox( boxwidth );

	TR << "get pose" << std::endl;
	get_pose_crd_and_index( pose, protein_crd, coarse_index, res_id );

	TR << "duplicate box" << std::endl;
	Vector mincrds, maxcrds;
	bring_to_origin( protein_crd, maxcrds, mincrds );

	// report
	/*
	if( report_crd_ ){
	for( core::Size iatm = 1; iatm <= protein_crd.size(); ++iatm )
	printf("ATOM  %5d  CA  ALA %5d    %8.3f%8.3f%8.3f 1.00  0.00 C\n",
	int(iatm), int(iatm), protein_crd[iatm][0], protein_crd[iatm][1], protein_crd[iatm][2]);
	}
	*/

	duplicate_waterbox( unit_waterbox_crd, boxwidth, mincrds, maxcrds );

	utility::vector1< Vector > const protein_crd0( protein_crd );
	utility::vector1< utility::vector1< core::Real > > sc_depth_stack;

	for ( core::Size i = 1; i <= niter_; ++i ) {
		TR << "iter " << i << std::endl;
		utility::vector1< Vector > protein_coarse_crd
			= get_coarse_crd( protein_crd, coarse_index );
		TR << "wat/crd/coarse: " << waterbox_.size() << " " << protein_crd.size() << " " << coarse_index.size() << std::endl;

		utility::vector1< bool > exclude_wat
			= get_exclusion_index( protein_crd, protein_coarse_crd );

		utility::vector1< core::Real > sc_depth =
			get_scdepth( exclude_wat, protein_crd, res_id );

		//bool converged =
		stack_and_getaverage( sc_depth_stack, sc_depth, i );

		//for( core::Size ires = 1; ires <= sc_depth.size(); ++ires ){
		// printf("%4d %5s %8.5f\n", int(ires), pose.residue(ires).name().c_str(), sc_depth[ires] );
		//}

		//if( converged ) break;
		protein_crd = protein_crd0;
		pert_protein( protein_crd );
	}
	TR << "estimation done!" << std::endl;

	return sc_depth_avrg_;
}

////////////////////////////////////////////////////////////////////////////
///Filter

ResidueDepthFilter::ResidueDepthFilter( core::pose::Pose const &pose )
{
	initialize( pose );
}

ResidueDepthFilter::~ResidueDepthFilter()= default;

/*
ResidueDepthFilter::ResidueDepthFilter( ResidueDepthFilter const &init )
{
initialize( pose );
}
*/

bool
ResidueDepthFilter::apply( core::pose::Pose const & pose ) const
{
	using namespace ObjexxFCL::format;

	get_residue_depth( pose );
	//core::Real score = SDE.get_SDE_score( pose );

	utility::vector1< core::Real > sde = get_scdepth_avrg();
	utility::vector1< core::Real > sde_var = get_scdepth_fvar();

	/*
	long t2=clock();
	double time = ((double)t2 - t1) / CLOCKS_PER_SEC;
	std::cout << "Done in " << time << " sec." << std::endl;
	*/
	bool passed = true;

	TR << "===================================================================" << std::endl;
	TR << " Report starts " << std::endl;
	TR << "===================================================================" << std::endl;
	TR << "#Ires  Eval PDBRES Chain  AA Depth    Var" << std::endl;
	for ( core::Size ires = 1; ires <= sde.size(); ++ires ) {
		if ( !evalres_[ires] ) continue;
		if ( sde[ires] > maxdist_ || sde[ires] < mindist_ ) {
			passed = false;
			//break;
		}

		TR << I(4, ires) << "    " << evalres_[ires]
			<< "   " << I(6, pose.pdb_info()->number( ires ) )
			<< "     " << pose.pdb_info()->chain( ires )
			<< " " << pose.residue(ires).name().c_str()
			<< " " << F(8,5,sde[ires])
			<< " " << F(8,5,sde_var[ires]) << std::endl;
	}
	TR.flush();
	return passed;
}

void
ResidueDepthFilter::initialize( core::pose::Pose const &pose )
{
	RDC_ = ResidueDepthCalculator( pose );
	ncandidate1_ = 500;
	ncandidate2_ = 40;
	//similarity_mode_ = "simple";
	similarity_mode_ = "superimpose";
	//dbfile_ = "protocol_data/SDE/DEPTH2.bench";
	GUIP_matrix_file_ = "protocol_data/SDE/GUIP.matrix";
	gamma_ = 1.0;
	maxdist_ = 99.0;
	mindist_ = 0.0;

	// TODO: comparing wrt. fragment SDE data
	//read_db( dbfile_ );
	read_GUIP_matrix( GUIP_matrix_file_ );

	// by default eval no res
	evalres_ = utility::vector1< bool >( RDC_.nres(), false );
}

void
ResidueDepthFilter::read_db( std::string const infile )
{
	TR << "Read DB!" << std::endl;

	// read dbfile
	core::Size npdb( 0 );
	std::string pdbid_prv( "" );
	utility::io::izstream instream;
	std::string line;

	basic::database::open( instream, infile );

	core::Size ndat( 0 );
	while ( instream ) {
		getline( instream, line );

		core::Size ires, n8;
		std::string pdbid, aa;
		core::Real x, y, z, xcen, ycen, zcen, depth;

		std::istringstream linestream( line );
		linestream >> pdbid >> aa >> ires >> n8 >> x >> y >> z >> xcen >> ycen >> zcen >> depth;

		if ( pdbid != pdbid_prv ) {
			npdb++;
			utility::vector1< ResidueDepthDataCOP > pdb_rdd;
			db_.push_back( pdb_rdd );
			pdbid_.push_back( pdbid );

		}
		ResidueDepthDataOP rdd = ResidueDepthDataOP( new ResidueDepthData );

		rdd->CAcrd[0] = x; rdd->CAcrd[1] = y; rdd->CAcrd[2] = z;
		rdd->CENcrd[0] = xcen; rdd->CENcrd[1] = ycen; rdd->CENcrd[2] = zcen;
		rdd->depth = depth;
		rdd->aa = aa;
		rdd->n8 = n8;
		rdd->ipdb = npdb;
		rdd->ires = ires;

		ResidueDepthDataCOP rdd_c =
			utility::pointer::static_pointer_cast< ResidueDepthData const >( rdd );

		db_[npdb].push_back( rdd_c );
		pdbid_prv = pdbid;
		ndat++;
	}

	TR << "Stash DB! npdb/ndat: " << db_.size() << " " << ndat << std::endl;

	// stash into frag_db_
	std::string pdbid("");
	for ( core::Size ipdb = 1; ipdb <= db_.size(); ++ipdb ) {
		if ( db_[ipdb].size() < 9 ) continue;
		for ( core::Size ires = 5; ires <= db_[ipdb].size()-4; ++ires ) {
			core::Size const n8 = db_[ipdb][ires]->n8;

			bool add( true );
			ResidueDepthFrag frag;
			int res0 = int(db_[ipdb][ires]->ires);

			for ( int k = -4; k <= 4; ++k ) {
				ResidueDepthDataCOP rdr = db_[ipdb][core::Size(ires+k)];
				if ( int(rdr->ires) - res0 != k ) {
					add = false;
					//TR << "failed stash: " << ipdb << " " << ires << std::endl;
					break;
				}
				frag.set_rdd( rdr, core::Size(k+5) );
			}
			if ( add ) frag_db_[n8].push_back( frag );
			//TR << "stash: " << n8 << " " << frag_db_.at(n8).size()
			//  << " " << frag.ipdb() << " " << frag.ires() << std::endl;
		}
	}

	TR << "ReadDB done! " << std::endl;

}

// read GUIP matrix
void
ResidueDepthFilter::read_GUIP_matrix( std::string const infile )
{
	GUIP_matrix_.resize( 20 );

	utility::io::izstream instream;
	std::string line;

	basic::database::open( instream, infile );

	core::Size iaa( 0 );
	while ( instream ) {
		getline( instream, line );
		core::Real val;
		std::istringstream linestream( line );

		iaa++;
		if ( iaa > 20 ) break;
		for ( core::Size jaa = 1; jaa <= 20; ++jaa ) {
			linestream >> val;
			//TR << "iaa/jaa/val: " << iaa << " " << jaa << " " << val << std::endl;
			GUIP_matrix_[iaa].push_back( val );
		}
	}
}

// MAIN FUNCTION
// wrapper for residue depth calculation
utility::vector1< core::Real >
ResidueDepthFilter::get_residue_depth( core::pose::Pose const &pose ) const
{
	// calculate the input pose's depth
	utility::vector1< core::Real > sde = RDC_.estimate_sidechain_depth( pose );
	return sde;
}

// this is the main
core::Real
ResidueDepthFilter::get_SDE_score( core::pose::Pose const &pose )
{
	// calculate the input pose's depth
	utility::vector1< core::Real > sde = RDC_.estimate_sidechain_depth( pose );

	// use centroid from here
	protocols::moves::MoverOP tocen
		( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );

	core::pose::Pose pose_cen( pose );
	tocen->apply( pose_cen );

	// get context info
	utility::vector1< core::Size > n8 = get_n8( pose_cen );

	// scan over DB...
	core::Real SDEscore( 0.0 );
	for ( core::Size ires = 5; ires <= sde.size()-4; ++ires ) {

		if ( !frag_db_.count( n8[ires] ) ) {
			TR << "Run res " << ires << ", n8: " << n8[ires] << ", not in DB!" << std::endl;
			continue;
		} else {
			TR << "Run res " << ires << ", n8: " << n8[ires] << ", searching db size of " << frag_db_.at( n8[ires] ).size() << std::endl;
		}

		// make 9mer crd
		utility::vector1< Vector > frag_crd( 9 );
		for ( core::Size k = 1; k <= 9; ++k ) {
			frag_crd[k] = pose_cen.residue( ires+k-5 ).xyz( " CA " );
		}

		//utility::vector1< numeric::xyzMatrix< core::Real > > uus_top;
		//utility::vector1< utility::vector1< core::Real > > wws_top;

		// filter1: get the list of frags with similar backbones
		utility::vector1< core::Size > close_ids1 =
			search_close_frags( frag_crd, n8[ires], ncandidate1_ );

		// make context of ires
		utility::vector1< ResidueDepthDataCOP > context_rdd = make_context( pose_cen, ires );

		// get summed sidechain similarity
		core::Real SDEscore_res = get_residue_similarity( pose.residue( ires ),
			frag_crd, n8[ires], context_rdd,
			close_ids1, ncandidate2_ );
		TR << "SDEscore(ires/score): " << ires << " " << SDEscore_res << std::endl;
		SDEscore += SDEscore_res;
	}

	return SDEscore;
}

utility::vector1< core::Size >
ResidueDepthFilter::get_n8( core::pose::Pose const &pose ) const
{
	utility::vector1< core::Size > n_neigh( pose.total_residue(), 0 );

	// use centroid
	for ( core::Size ires = 1; ires <= pose.total_residue(); ++ires ) {
		if ( pose.residue(ires).is_virtual_residue() ) continue;

		utility::vector1< ResidueDepthDataCOP > context_rdds = make_context( pose, ires );
		n_neigh[ires] = context_rdds.size();
	}

	return n_neigh;
}

utility::vector1< core::Size >
ResidueDepthFilter::search_close_frags( utility::vector1< Vector > const &frag_crd,
	core::Size const n8,
	core::Size const npick
) const
{
	//utility::vector1< core::Real > ww( 27, 1.0 );
	//numeric::xyzMatrix< core::Real > uu;

	core::Real rmsd; //, ctx;
	//float rmsd_f;
	ObjexxFCL::FArray1D< core::Real > ww( 9, 1.0 );
	ObjexxFCL::FArray2D< core::Real > uu( 3, 3, 0.0 );
	ObjexxFCL::FArray2D< core::Real > r_1_ref( 3, 9 );

	for ( int j = 1; j <= 9; ++j ) {
		for ( int k = 1; k <= 3; ++k ) r_1_ref(k,j) = frag_crd[j][k-1];
	}

	utility::vector1< ResidueDepthFrag > const &fragdb_n = frag_db_.at( n8 );
	core::Size const ndb( fragdb_n.size() );

	std::vector< core::Real > rmsds( ndb );

	for ( core::Size i = 1; i <= ndb; ++i ) {
		ResidueDepthFrag const &frag_ref = fragdb_n[i];

		ObjexxFCL::FArray2D< core::Real > r_1( r_1_ref );
		ObjexxFCL::FArray2D< core::Real > r_2( 3, 9 );
		for ( int j = 1; j <= 9; ++j ) {
			Vector const & v2 = frag_ref.get_CAcrd(j);
			for ( int k = 1; k <= 3; ++k ) r_2(k,j) = v2[k-1];
		}

		// this is the time determining
		rmsd = numeric::model_quality::rms_wrapper( 9, r_1, r_2 );

		// stash
		rmsds[i-1] = rmsd;
	}

	std::vector< core::Real > sortable( rmsds );

	//sort by rmsd
	std::sort( sortable.begin(), sortable.end() );
	core::Real const rmsdcut = sortable[npick-1];

	TR << "sorted: ";
	for ( core::Size k = 1; k <= 5; ++k ) TR << " " << sortable[k];
	TR << std::endl;

	// get indices
	utility::vector1< core::Size > close_ids( npick, 0 );

	core::Size ifilt( 0 );
	for ( core::Size i = 1; i <= ndb; ++i ) {
		if ( rmsds[i] <= rmsdcut ) {
			ifilt++;
			if ( ifilt > npick ) break;

			close_ids[ifilt] = i;

			/*
			ResidueDepthFrag const &frag_ref = fragdb_n[i];
			//TR << "close? " << i << " " << frag_ref.ipdb() << " " << frag_ref.ires() <<  std::endl;

			ObjexxFCL::FArray2D< core::Real > r_1( r_1_ref );
			ObjexxFCL::FArray2D< core::Real > r_2( 3, 9 );
			for( int j = 1; j <= 9; ++j ){
			Vector const & v2 = frag_ref.get_CAcrd(j);
			for( int k = 1; k <= 3; ++k ) r_2(k,j) = v2[k-1];
			}

			numeric::model_quality::findUU( r_1, r_2, ww, 9, uu, ctx );
			numeric::model_quality::calc_rms_fast( rmsd_f, r_1, r_2, ww, 9, ctx );
			*/

		}
	}

	return close_ids;
}

core::Real
ResidueDepthFilter::get_residue_similarity( conformation::Residue const &rsd,
	utility::vector1< Vector > const &frag_crd,
	core::Size const n8,
	utility::vector1< ResidueDepthDataCOP > const context_rdd,
	utility::vector1< core::Size > const &close_ids,
	core::Size npick ) const
{
	std::vector< std::pair< core::Real, core::Size > > sortable;
	utility::vector1< ResidueDepthFrag > const & frag_db_n = frag_db_.at( n8 );

	core::Size npass( 0 );
	for ( core::Size i = 1; i <= close_ids.size(); ++i ) {
		if ( close_ids[i] == 0 || close_ids[i] > frag_db_n.size() ) continue;

		ResidueDepthFrag const & frag_ref = frag_db_n[ close_ids[i] ];
		core::Size const ipdb( frag_ref.ipdb() );
		core::Size const ires( frag_ref.ires() );

		if ( ipdb == 0 || ipdb > db_.size() ) continue;
		if ( ires == 0 || ires > db_[ipdb].size() ) continue;

		utility::vector1< Vector > frag_crd_ref = frag_ref.get_CENcrd();

		//TR << "ipdb,ires: " << ipdb << " " << ires << std::endl;

		core::Real rms( 0.0 );
		if ( similarity_mode_.compare( "superimpose") == 0 ) {

			// from db
			utility::vector1< ResidueDepthDataCOP > context_rdd_ref = make_context_ref( ipdb, ires );

			// apply transrot to context
			// returns RMS( delta_SDE )
			if ( context_rdd.size() != context_rdd_ref.size() ) {
				TR << i << ": unequal num: " << context_rdd.size() << " " << context_rdd_ref.size() << std::endl;
				rms = 99.9;
				npass++;
			} else {
				rms = compare_by_superposition( context_rdd, frag_crd,
					context_rdd_ref, frag_crd_ref );
				TR << i << ": equal num." << context_rdd.size() << ", rms = " << rms << std::endl;
			}

		} else if ( similarity_mode_.compare( "simple" ) == 0 ) {
			// TODO
		}

		sortable.push_back( std::make_pair( rms, i ) ); // trick for sorting...
	}

	TR << "skiped " << npass << " out of " << close_ids.size() << std::endl;

	/*
	std::sort( sortable.begin(), sortable.end(),
	[](std::pair< core::Real, core::Size > a,
	std::pair< core::Real, core::Size > b)
	{return a.first < b.first;}
	);
	*/
	std::sort( sortable.begin(), sortable.end(), mycomp );

	TR << "RMS(SDE) sort:";
	for ( core::Size k = 0; k < 5; ++k ) TR << " " << sortable[k].first << ":" << sortable[k].second;
	TR << std::endl;

	core::Real SDE( 0.0 );
	// sum up to n (weighted sum)
	for ( core::Size i = 0; i < npick; ++i ) { // this is rank
		core::Real weight = exp( -gamma_*core::Real(i)/core::Real(npick-1) );

		core::Size i_rdd = sortable[i].second;

		ResidueDepthFrag const & frag_ref = frag_db_n[ close_ids[i_rdd] ];
		core::Size const ipdb( frag_ref.ipdb() );

		core::Real simscore = get_simscore( rsd.aa(),
			core::chemical::aa_from_oneletter_code( *(frag_ref.aa().c_str()) ) );

		TR << "Final score from " << i << " (aa " << core::chemical::oneletter_code_from_aa( rsd.aa() )
			<< "-" << frag_db_n[close_ids[i_rdd]].aa()
			<< ", pdb/res " << pdbid_[ipdb] << "/" << frag_ref.ires()
			<< "): " << weight*simscore << " (" << weight << " * " << simscore << ") "
			<< " , RMS(SDE) = " << sortable[i].first
			<< std::endl;
		SDE += weight*simscore;
	}

	return SDE;
}

core::Real
ResidueDepthFilter::get_simscore( core::chemical::AA const aa1,
	core::chemical::AA const aa2 ) const
{
	// this shouldn't happen, but to be safe anyway...
	if ( aa1 > core::chemical::num_canonical_aas || aa1 > core::chemical::num_canonical_aas ) {
		return 0.0;
	} else {
		return GUIP_matrix_[aa1][aa2];
	}
}

core::Real
ResidueDepthFilter::compare_by_superposition( utility::vector1< ResidueDepthDataCOP > const & context_rdd,
	utility::vector1< Vector > const & frag_crd,
	utility::vector1< ResidueDepthDataCOP > const & context_rdd_ref,
	utility::vector1< Vector > const & frag_crd_ref
) const
{
	//runtime_assert( context_rdd.size() == context_rdd_ref.size() );

	//TR << "ncontext: " << context_rdd.size() << " " << context_rdd_ref.size() << " " << frag_crd.size() << " " << frag_crd_ref.size() << std::endl;

	core::Size ncontext( context_rdd.size() );
	int n = 9 + int(ncontext);
	core::Real rmsd;

	// first 9 are
	ObjexxFCL::FArray2D< core::Real > r_1( 3, n, 0.0 );
	ObjexxFCL::FArray2D< core::Real > r_2( 3, n, 0.0 );
	ObjexxFCL::FArray1D< core::Real > ww( n, 0.0 );

	// fill in r_1, r_2
	for ( int i = 1; i <= 9; ++i ) {
		ww( i ) = 1.0;
		for ( int k = 1; k <= 3; ++k ) {
			//r_1( k, i ) = frag_crd[core::Size(i)][core::Size(k-1)];
			//r_2( k, i ) = frag_crd_ref[core::Size(i)][core::Size(k-1)];
			r_1( k, i ) = frag_crd[i][k-1];
			r_2( k, i ) = frag_crd_ref[i][k-1];
		}
	}

	for ( int i = 10; i <= n; ++i ) {
		ww( i ) = 0.0;
		Vector const &crd1 = context_rdd[i-9]->CENcrd;
		Vector const &crd2 = context_rdd_ref[i-9]->CENcrd;
		for ( int k = 1; k <= 3; ++k ) {
			r_1( k, i ) = crd1[k-1];
			r_2( k, i ) = crd2[k-1];
		}
	}

	//for( core::Size i = 1; i <= ncontext; ++i )
	//TR << "depth: " << i << " " << context_rdd[i]->depth << " " << context_rdd_ref[i]->depth << std::endl;

	numeric::model_quality::rmsfitca2( n, r_1, r_2, ww, n, rmsd );

	// convert into vector1< Vector >
	utility::vector1< Vector > context_crd_aligned( ncontext ), context_crd_ref_aligned( ncontext );
	for ( core::Size i = 10; i <= core::Size(n); ++i ) {
		for ( core::Size k = 1; k <= 3; ++k ) {
			context_crd_aligned[i-9][k-1] = r_1( k, i );
			context_crd_ref_aligned[i-9][k-1] = r_2( k, i );
		}
	}

	// pick closest coordinates (from input side)
	// allow overlap
	core::Real rms( 0.0 );
	for ( core::Size i = 1; i <= ncontext; ++i ) {
		Vector const &crd1 = context_crd_aligned[i];
		core::Real const depth1( context_rdd[i]->depth );

		// first search for closest distance
		core::Real d2closest( 1e6 );
		core::Size jmin( 1 );
		for ( core::Size j = 1; j <= ncontext; ++j ) {
			Vector const &crd2 = context_crd_ref_aligned[j];

			core::Real d2( crd1.distance_squared( crd2 ) );
			if ( d2 < d2closest ) {
				d2closest = d2; jmin = j;
			}
		}

		core::Real const depth2( context_rdd_ref[jmin]->depth );
		TR << "i/j/d/depth1/depth2: " << i << " " << jmin << " " << std::sqrt(d2closest) << " " << depth1 << " " << depth2 << std::endl;

		rms += ( depth1 - depth2 )*( depth1 - depth2 );
	}

	return std::sqrt(rms/core::Real(ncontext));
}

utility::vector1< ResidueDepthDataCOP >
ResidueDepthFilter::make_context_ref( core::Size const ipdb,
	core::Size const ires ) const
{
	utility::vector1< ResidueDepthDataCOP > context_rdds;
	core::Real const D2CUT( 64.0 );

	ResidueDepthDataCOP rdd_i( db_[ipdb][ires] );
	Vector const iCENcrd( rdd_i->CENcrd );

	for ( core::Size jres = 1; jres <= db_[ipdb].size(); ++jres ) {
		if ( ires == jres ) continue;

		ResidueDepthDataCOP rdd_j( db_[ipdb][jres] );
		//core::Size jres = db_[ipdb]->ires;
		Vector const jCENcrd( rdd_j->CENcrd );

		core::Real const d2( iCENcrd.distance_squared( jCENcrd ) );
		if ( d2 < D2CUT ) {
			context_rdds.push_back( rdd_j );
		}
	}

	return context_rdds;
}

utility::vector1< ResidueDepthDataCOP >
ResidueDepthFilter::make_context( core::pose::Pose const &pose,
	core::Size const ires ) const
{
	utility::vector1< ResidueDepthDataCOP > context_rdds;
	core::Real const D2CUT( 100.0 );

	Vector const iCENcrd = pose.residue(ires).xyz( " CEN" );

	// use CEN to be consistent with other data
	for ( core::Size jres = 1; jres <= pose.total_residue(); ++jres ) {
		if ( ires == jres ) continue;
		if ( pose.residue( jres ).is_virtual_residue() ) continue;

		//TR << "context! " << ires << " " << jres << std::endl;
		Vector const &jCENcrd = pose.residue(jres).xyz( " CEN" );

		core::Real const d2( iCENcrd.distance_squared( jCENcrd ) );
		//TR << "ires/jres/d2: " << ires << " " << jres << " " << d2 << std::endl;
		if ( d2 < D2CUT ) {
			ResidueDepthDataOP rdd = ResidueDepthDataOP( new ResidueDepthData );
			// don't need others than below
			rdd->CAcrd = pose.residue(jres).xyz( " CA " );
			rdd->CENcrd = jCENcrd;
			rdd->ires = jres;
			rdd->aa = core::chemical::name_from_aa( pose.residue(jres).aa() );
			rdd->depth = RDC_.get_scdepth_avrg( jres );

			ResidueDepthDataCOP rdd_c =
				utility::pointer::static_pointer_cast< ResidueDepthData const >( rdd );
			context_rdds.push_back( rdd_c );
		}
	}

	return context_rdds;
}

void
ResidueDepthFilter::fill_neighs( ResidueDepthData &rdd,
	utility::vector1< ResidueDepthData > const & pdb_rdds ) const
{
	Vector const &cencrd( rdd.CENcrd );
	core::Real const D2CUT( 100.0 );

	for ( core::Size ires = 1; ires <= pdb_rdds.size(); ++ires ) {
		if ( ires == rdd.ires ) continue;

		core::Real const d2 = cencrd.distance_squared( pdb_rdds[ires].CENcrd );
		if ( d2 < D2CUT ) rdd.neigh_ress.push_back( ires );
	}
}

void
ResidueDepthFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose)
{
	using namespace ObjexxFCL;

	initialize( pose );

	mindist_ = tag->getOption<core::Real>( "mindist", 0.0 );
	maxdist_ = tag->getOption<core::Real>( "maxdist", 99.0 );

	if ( tag->hasOption( "dcut1" ) ) {
		RDC_.set_dcut1( tag->getOption<core::Real>( "dcut1", 2.6 ) );
	}
	if ( tag->hasOption( "dcut2" ) ) {
		RDC_.set_dcut2( tag->getOption<core::Real>( "dcut2", 4.2 ) );
	}


	evalres_ = utility::vector1< bool >( pose.total_residue(), false );

	if ( tag->hasOption("evalres") ) {
		utility::vector1< std::string> evalres_str
			= utility::string_split(tag->getOption<std::string>( "evalres" ), ',');
		for ( core::Size i = 1; i <= evalres_str.size(); ++i ) {
			Size ires = (Size)(int_of( evalres_str[i] ));
			evalres_[ires] = true;
		}
	}
}

} // namespace simple_filters
} // namespace protocols
