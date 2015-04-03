// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief protocols for folding into density
/// @details
/// @author Frank DiMaio

#include <protocols/rbsegment_relax/IdealizeHelices.hh>
#include <protocols/rbsegment_relax/IdealizeHelicesCreator.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/simple_moves/FragmentMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Remarks.hh>
#include <core/pose/util.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/symmetry/SetupNCSMover.hh>

#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>

#include <numeric/random/random.hh>

#include <protocols/hybridization/util.hh>
#include <protocols/loops/util.hh>
#include <protocols/moves/DsspMover.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/model_quality/rms.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>


using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace rbsegment_relax {

static thread_local basic::Tracer TR( "protocols.rbsegment_relax.IdealizeHelices" );

using namespace protocols;
using namespace core;


//////////////////
//////////////////
/// creator


std::string
IdealizeHelicesMoverCreator::keyname() const {
	return IdealizeHelicesMoverCreator::mover_name();
}

protocols::moves::MoverOP
IdealizeHelicesMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new IdealizeHelicesMover );
}

std::string
IdealizeHelicesMoverCreator::mover_name() {
	return "IdealizeHelices";
}


//////////////////
//////////////////
/// mover

void IdealizeHelicesMover::apply( core::pose::Pose & pose ) {
	using namespace core::scoring::constraints;
	using namespace core::optimization;
	using namespace core::id;

	// if no helices given use DSSP
	// TO DO

	// see if the pose has NCS
	simple_moves::symmetry::NCSResMappingOP ncs;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ) {
		ncs = ( utility::pointer::static_pointer_cast< simple_moves::symmetry::NCSResMapping > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ));
	}

	// pose to centroid
	bool fullatom_input = pose.is_fullatom();
	protocols::moves::MoverOP restore_sc;
	if (fullatom_input) {
	  restore_sc = protocols::moves::MoverOP( new protocols::simple_moves::ReturnSidechainMover( pose ) );
		protocols::moves::MoverOP tocen( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );
		tocen->apply( pose );
	}

	// prolines
	utility::vector1< std::pair<int,int> > corrected_helices;
	for (core::uint i = 1; i <= helices_.size(); ++i) {
		core::uint start_i = helices_[i].first, stop_i = helices_[i].second;

		for (core::uint j = start_i+2; j <= stop_i - 2; ++j) {
			if (pose.residue(j).aa() == core::chemical::aa_pro) {
				if (numeric::random::rg().uniform() <= 0.5) {
					if (start_i < j - 10) corrected_helices.push_back( std::make_pair( start_i, j - 6 ) );
					start_i = j - 2;
				}
			}
		}
		if (start_i < stop_i - 4) corrected_helices.push_back( std::make_pair( start_i, stop_i ) );

	}


	// for each res range
	for (core::uint i = 1; i <= corrected_helices.size(); ++i) {
		core::uint start_i = corrected_helices[i].first, stop_i = corrected_helices[i].second;

		// build ideal pose
		core::pose::Pose ideal_pose;
		Size len_i = stop_i - start_i + 1;
		// ideal geometry
		for (core::uint j = start_i; j <= stop_i; ++j) {
			ideal_pose.append_residue_by_bond( pose.residue( j ), true );
		}

		for (core::uint j = start_i; j <= stop_i; ++j) {
			protocols::loops::set_extended_torsions_and_idealize_loops( ideal_pose, protocols::loops::Loops() );
		}

		// helical
		for (core::uint j = start_i; j <= stop_i; ++j) {
			core::uint resid = j - start_i + 1;
			ideal_pose.set_phi( resid, -57.8 );
			ideal_pose.set_psi( resid, -47.0 );
			ideal_pose.set_omega( resid, 180.0 );
		}
		core::pose::addVirtualResAsRoot( ideal_pose );
		core::Size vrt_index = ideal_pose.total_residue();

		// superimpose
		ObjexxFCL::FArray1D< numeric::Real > ww( len_i, 1.0 );
		ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
		numeric::xyzVector< core::Real > com1(0,0,0), com2(0,0,0);

		// grab source coords
		ObjexxFCL::FArray2D< core::Real > init_coords( 3, len_i );
		for (core::uint j = start_i; j <= stop_i; ++j) {
			core::uint resid = j - start_i + 1;
			numeric::xyzVector< core::Real > x_j = ideal_pose.residue(resid).atom(" CA ").xyz();
			com1 += x_j;
			for (core::uint k = 0; k < 3; ++k) {
				init_coords(k + 1, resid) = x_j[k];
			}
		}
		com1 /= len_i;
		for (core::uint j = 0; j < len_i; ++j) {
			for ( core::uint k = 0; k < 3; ++k ) {
				init_coords(k + 1, j + 1) -= com1[k];
			}
		}

		// grab target coords
		ObjexxFCL::FArray2D< core::Real > final_coords( 3, len_i );
		for (core::uint j = start_i; j <= stop_i; ++j) {
			core::uint resid = j - start_i + 1;
			numeric::xyzVector< core::Real > x_j = pose.residue(j).atom(" CA ").xyz();

			// add CSTS
			for (Size iatom = 1; iatom <= 4; ++iatom) {
				ideal_pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID(iatom,resid), AtomID(1,vrt_index), x_j, core::scoring::func::FuncOP( new BoundFunc(0.0,cst_width_,1.0,"") ) ) ) ) );
			}

			com2 += x_j;
			for (core::uint k = 0; k < 3; ++k)
				final_coords(k + 1, resid) = x_j[k];
		}
		com2 /= len_i;
		for (core::uint j = 0; j < len_i; ++j) {
			for ( core::uint k = 0; k < 3; ++k ) {
				final_coords(k + 1, j + 1) -= com2[k];
			}
		}

		// superpose
		numeric::Real ctx; float rms;
		numeric::model_quality::findUU( init_coords, final_coords, ww, len_i, uu, ctx );
		numeric::model_quality::calc_rms_fast( rms, init_coords, final_coords, ww, len_i, ctx );
		numeric::xyzMatrix< core::Real > R;
		R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
		R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
		R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );
		for (core::uint j = start_i; j <= stop_i; ++j) {
			core::uint resid = j - start_i + 1;
			for ( Size k = 1; k <= ideal_pose.residue_type(resid).natoms(); ++k ) {
				core::id::AtomID id_src( k, resid );
				ideal_pose.set_xyz( id_src, R * ( ideal_pose.xyz(id_src) - com1) + com2 );
			}
		}

		// minimize with CSTS
		core::kinematics::MoveMap mm; mm.set_bb(true); mm.set_jump(true);
		MinimizerOptions min_options( "lbfgs_armijo", 0.0001, true, false, false );
		AtomTreeMinimizer minimizer;
		(*scorefxn_)(ideal_pose);
		minimizer.run( ideal_pose, mm, *scorefxn_, min_options );

		// replace coords
		for (core::uint j = start_i; j <= stop_i; ++j) {
			core::uint resid = j - start_i + 1;
			for ( Size k = 1; k <= ideal_pose.residue_type(resid).natoms(); ++k ) {
				core::id::AtomID id_src( k, resid );
				core::id::AtomID id_tgt( k, j );
				pose.set_xyz( id_tgt, ideal_pose.xyz(id_src) );
			}
		}

		// apply to NCS-symmetric copies
		if (ncs) {
			for (core::uint n = 1; n <= ncs->ngroups(); ++n ) {
				bool all_are_mapped = true;
				for ( Size k = start_i; k <= stop_i && all_are_mapped; ++k )
					all_are_mapped &= (ncs->get_equiv( n,k )!=0);
				if (!all_are_mapped) continue;

				core::Size remap_start = ncs->get_equiv( n, start_i );
				core::Size remap_stop = ncs->get_equiv( n, stop_i );

				if (remap_stop - remap_start != stop_i - start_i) continue;

				// superimpose
				ObjexxFCL::FArray1D< numeric::Real > ww( len_i, 1.0 );
				ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
				numeric::xyzVector< core::Real > com1(0,0,0), com2(0,0,0);

				// grab source coords
				ObjexxFCL::FArray2D< core::Real > init_coords( 3, len_i );
				for (core::uint j = start_i; j <= stop_i; ++j) {
					numeric::xyzVector< core::Real > x_j = pose.residue(j).atom(" CA ").xyz();
					com1 += x_j;
					for (core::uint k = 0; k < 3; ++k) {
						init_coords(k + 1, j - start_i + 1) = x_j[k];
					}
				}
				com1 /= len_i;
				for (core::uint j = 0; j < len_i; ++j) {
					for ( core::uint k = 0; k < 3; ++k ) {
						init_coords(k + 1, j + 1) -= com1[k];
					}
				}

				// grab target coords
				ObjexxFCL::FArray2D< core::Real > final_coords( 3, len_i );
				for (core::uint j = remap_start; j <= remap_stop; ++j) {
					numeric::xyzVector< core::Real > x_j = pose.residue(j).atom(" CA ").xyz();
					com2 += x_j;
					for (int k=0; k<3; ++k) final_coords(k+1,j-remap_start+1) = x_j[k];
				}
				com2 /= len_i;
				for (core::uint j = 0; j < len_i; ++j) {
					for ( core::uint k = 0; k < 3; ++k ) {
						final_coords(k + 1, j + 1) -= com2[k];
					}
				}

				// superpose
				numeric::Real ctx; float rms;
				numeric::model_quality::findUU( init_coords, final_coords, ww, len_i, uu, ctx );
				numeric::model_quality::calc_rms_fast( rms, init_coords, final_coords, ww, len_i, ctx );
				numeric::xyzMatrix< core::Real > R;
				R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
				R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
				R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );
				for (core::uint j = start_i; j <= stop_i; ++j) {
					core::uint offset = j - start_i;
					for ( Size k = 1; k <= pose.residue_type(j).natoms(); ++k ) {
						core::id::AtomID id_src( k, j );
						core::id::AtomID id_tgt( k, offset+remap_start );
						pose.set_xyz( id_tgt, R * ( pose.xyz(id_src) - com1) + com2 );
					}
				}
			}
		}
	}

	if (restore_sc)
		restore_sc->apply(pose);
}


void IdealizeHelicesMover::parse_my_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap & data,
			filters::Filters_map const & /*filters*/,
			moves::Movers_map const & /*movers*/,
			core::pose::Pose const & pose )
{
	if( tag->hasOption( "scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn" ) );
		scorefxn_ = (data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name ))->clone();
	}
	if( tag->hasOption( "cst_weight" ) ) {
		cst_weight_ = tag->getOption<core::Real>( "cst_weight" );
		scorefxn_->set_weight( core::scoring::coordinate_constraint , cst_weight_ );
	}

	cst_width_ = tag->getOption<core::Real>( "cst_width" , 0.0 );

	// fragments
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
	for (tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it) {
		if ( (*tag_it)->getName() == "Helix" ) {
			using namespace core::fragment;
			std::string start_i = (*tag_it)->getOption<std::string>( "start" );
			std::string end_i = (*tag_it)->getOption<std::string>( "end" );
			int startres = (int)core::pose::parse_resnum( start_i, pose );
			int stopres = (int)core::pose::parse_resnum( end_i, pose );
			helices_.push_back( std::make_pair( startres, stopres ) );
		}
	}
}

}
}
