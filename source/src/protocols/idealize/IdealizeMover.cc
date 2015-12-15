// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/idealize/IdealizeMover.cc
/// @brief
/// @author

// Unit Headers
#include <protocols/idealize/idealize.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/idealize/IdealizeMoverCreator.hh>
#include <utility/string_util.hh>

// // Rosetta Headers
#include <core/types.hh>
#include <utility/tag/Tag.hh>

#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>


#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <basic/Tracer.hh> // tracer output

// symmetry
#include <core/pose/symmetry/util.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Numeric headers

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

//Auto Headers
#include <core/chemical/AtomType.hh>
#include <core/kinematics/FoldTree.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


// // C++ Headers

namespace protocols {
namespace idealize {

using namespace core;

static THREAD_LOCAL basic::Tracer TR( "protocols.idealize.IdealizeMover" );

std::string
IdealizeMoverCreator::keyname() const
{
	return IdealizeMoverCreator::mover_name();
}

protocols::moves::MoverOP
IdealizeMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new IdealizeMover );
}

std::string
IdealizeMoverCreator::mover_name()
{
	return "Idealize";
}


/// @brief Add the idealize constraints to the pose's constraint set
void
IdealizeMover::setup_idealize_constraints( core::pose::Pose & pose ) {
	using namespace scoring::constraints;
	using namespace scoring::func;
	using namespace conformation;
	using namespace id;

	Real const heavyatom_dis2_threshold( 5.5 * 5.5 );
	Real const polarH_dis2_threshold( 2.5 * 2.5 );

	Real const atom_pair_sdev( 0.25 );
	Real const coord_sdev( 0.1 );


	Size const nres( pose.total_residue() );
	Size total_atompairs( 0 );

	//fpd symmetry
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		TR.Info << "setting up symmetric idealize " << std::endl;
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}

	if ( atom_pair_constraint_weight_ != 0.0 || symm_info ) {
		for ( Size i=1; i<= nres-1; ++i ) {
			if ( std::find( ignore_residues_in_csts().begin(), ignore_residues_in_csts().end(), i ) != ignore_residues_in_csts().end() ) continue;
			Residue const & i_rsd( pose.residue(i) );
			if ( i_rsd.aa() == core::chemical::aa_vrt ) continue;

			for ( Size j=i+1; j<= nres-1; ++j ) {
				if ( std::find( ignore_residues_in_csts().begin(), ignore_residues_in_csts().end(), j ) != ignore_residues_in_csts().end() ) continue;
				Residue const & j_rsd( pose.residue(j) );
				if ( j_rsd.aa() == core::chemical::aa_vrt ) continue;

				//fpd  for symmetry, we only need generate csts w.r.t. scoring subunit
				if ( symm_info && !symm_info->bb_is_independent( i ) && !symm_info->bb_is_independent( j ) ) continue;

				//fpd  if atom pair cst weight is 0, we _only_ generate atom pair csts across symm interface
				if ( symm_info && atom_pair_constraint_weight_ == 0.0 &&
						symm_info->bb_is_independent( i ) && symm_info->bb_is_independent( j ) ) {
					continue;
				}

				for ( Size ii = 1; ii<= i_rsd.natoms(); ++ii ) {
					chemical::AtomType const & it( i_rsd.atom_type( ii ) );

					for ( Size jj = 1; jj<= j_rsd.natoms(); ++jj ) {
						chemical::AtomType const & jt( j_rsd.atom_type( jj ) );

						Real const dis2( i_rsd.xyz( ii ).distance_squared( j_rsd.xyz( jj ) ) );

						if ( ( it.is_polar_hydrogen() && jt.is_acceptor()  && dis2 <    polarH_dis2_threshold ) ||
								( jt.is_polar_hydrogen() && it.is_acceptor()  && dis2 <    polarH_dis2_threshold ) ||
								( it.is_heavyatom()      && jt.is_heavyatom() && dis2 < heavyatom_dis2_threshold ) ) {

							pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint( AtomID(ii,i), AtomID(jj,j), core::scoring::func::FuncOP( new HarmonicFunc( std::sqrt( dis2 ), atom_pair_sdev ) ) ) ) ) );
							++total_atompairs;
						}
					} // jj
				} // ii
			} // j>=i
		} // i
	}

	TR.Info << "total atompairs: " << total_atompairs << std::endl;

	if ( coordinate_constraint_weight_ != 0.0 ) {
		// should already have setup for using coordinate constraints
		runtime_assert( pose.residue( nres ).aa() == core::chemical::aa_vrt );
		for ( Size i=1; i<= nres-1; ++i ) {
			// only put coord csts on master
			if ( symm_info &&
					( !symm_info->bb_is_independent( i ) || pose.residue(i).aa() == core::chemical::aa_vrt) ) {
				continue;
			}
			Residue const & i_rsd( pose.residue(i) );
			for ( Size ii = 1; ii<= i_rsd.natoms(); ++ii ) {
				pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID(ii,i), AtomID(1,nres), i_rsd.xyz( ii ), core::scoring::func::FuncOP( new HarmonicFunc( 0.0, coord_sdev ) ) ) ) ) );
			}
		}
	} // coordinate_constraint_weight_ != 0
	TR.flush();
} // setup_idealize_constraints

void
IdealizeMover::apply( pose::Pose & pose ) {
	using namespace scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	pose::Pose unmodified_pose = pose;
	//pose.dump_pdb("idl_initial.pdb");

	// save a copy of the pose's constraints
	scoring::constraints::ConstraintSetOP original_cst_set( pose.constraint_set()->clone() );
	if ( impose_constraints() ) {
		pose.constraint_set( NULL );
	}
	// add virtual residue at the end
	//Size const old_root( pose.fold_tree().root() );
	if ( pose.residue( pose.total_residue() ).aa() != core::chemical::aa_vrt ) {
		/// bugfix for single-residue pose: don't append residue by jump from residue 0
		Size const midpoint( pose.total_residue() == 1 ? 1 : pose.total_residue() / 2 );
		pose.append_residue_by_jump(
			*conformation::ResidueFactory::create_residue( pose.residue(1).residue_type_set()->name_map( "VRT" ) ),
			midpoint
		);

		Size const nres( pose.total_residue() ); // includes pseudo-rsd

		kinematics::FoldTree f( pose.fold_tree() );
		f.reorder( nres );
		pose.fold_tree( f );
	}

	//fpd symmetry
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}

	// setup scorefunction
	scoring::ScoreFunctionOP scorefxn;
	if ( symm_info ) {
		scorefxn = scoring::ScoreFunctionOP( new scoring::symmetry::SymmetricScoreFunction() );
	} else {
		scorefxn = scoring::ScoreFunctionOP( new scoring::ScoreFunction() );
	}
	scorefxn->set_weight( atom_pair_constraint,  atom_pair_constraint_weight_ );
	scorefxn->set_weight( coordinate_constraint, coordinate_constraint_weight_ );

	// if we're symmetric, we need to turn atom pair csts on
	if ( symm_info && atom_pair_constraint_weight_ == 0 ) {
		scorefxn->set_weight( atom_pair_constraint,  coordinate_constraint_weight_ );
	}


	if ( pose.is_fullatom() ) {
		// keep prolines closed during idealizations.
		scorefxn->set_weight( pro_close, 0.5 );
		// keep disulphides together.
		scorefxn->set_weight( dslf_ss_dst, 0.5 );//SG-SG bond length
		scorefxn->set_weight( dslf_cs_ang, 2.0 );//CB-SG-SG covalent bond angle
	}
	// setup constraints
	if ( impose_constraints() ) {
		setup_idealize_constraints( pose );
	}
	if ( constraints_only() ) {
		return;
	}

	// by default idealize everything
	//fpd  ... unless symmetric, then only idealize master
	if ( pos_list_.size() == 0 ) {
		for ( Size i = 1; i <= pose.total_residue()-1; ++i ) {
			if ( symm_info &&
					(!symm_info->bb_is_independent( i ) || pose.residue(i).aa() == core::chemical::aa_vrt ) ) {
				continue;
			}
			pos_list_.push_back( i );
		}
	}


	if ( !(option[ basic::options::OptionKeys::run::dry_run ]() ) ) {
		basic_idealize( pose, pos_list_, *scorefxn, fast_, chainbreaks_, cis_omega_ );
	}

	// remove that virtual residue now!
	pose::Pose final_pose = unmodified_pose;

	if ( symm_info ) {
		// special case for symmetry .. replace VRTs first
		for ( Size ii = unmodified_pose.total_residue(); ii>=1; --ii ) {
			if ( symm_info->bb_is_independent(ii) ) {
				final_pose.replace_residue( ii, pose.residue( ii ), false );
			}
		}
	} else {
		for ( Size ii = 1; ii <= unmodified_pose.total_residue(); ++ii ) {
			final_pose.replace_residue( ii, pose.residue( ii ), false );
		}
	}
	pose = final_pose;

	// restore the original constraint set
	pose.constraint_set( original_cst_set );

	/// Pose must be rescored after the original constraint set is restored.
	(*scorefxn)( pose );

	TR.Info << "RMS between original pose and idealised pose: ";
	if ( report_CA_rmsd_ ) TR.Info << core::scoring::CA_rmsd( unmodified_pose, pose ) << " CA RMSD, ";

	TR.Info << core::scoring::all_atom_rmsd( unmodified_pose, pose ) << " All-Atom RMSD, "
		<< std::endl;
	TR.flush();
} // apply

std::string
IdealizeMover::get_name() const {
	return "IdealizeMover";
}

void
IdealizeMover::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose ){
	atom_pair_constraint_weight( tag->getOption< core::Real >( "atom_pair_constraint_weight", 0.0 ) );
	coordinate_constraint_weight( tag->getOption< core::Real >( "coordinate_constraint_weight", 0.01 ) );
	fast( tag->getOption< bool >( "fast", false ) );
	chainbreaks( tag->getOption< bool >( "chainbreaks", false ) );
	report_CA_rmsd( tag->getOption< bool >( "report_CA_rmsd", true ) );
	if ( tag->hasOption( "ignore_residues_in_csts" ) ) {
		ignore_residues_in_csts( core::pose::get_resnum_list( tag, "ignore_residues_in_csts", pose ) );
	}
	impose_constraints( tag->getOption< bool >( "impose_constraints", 1 ) );
	constraints_only( tag->getOption< bool >( "constraints_only", 0 ) );
	if ( tag->hasOption( "pos_list" ) ) {
		pos_list_ = utility::string_split( tag->getOption< std::string >( "pos_list", "" ), ',', core::Size() );
	}

	TR<<"IdealizeMover with atom_pair_constraint_weight="<<atom_pair_constraint_weight_<<" coordinate_constraint_weight="<<coordinate_constraint_weight_<<" fast="<<fast_<<" chainbreaks="<<chainbreaks_<<" and report CA_rmsd_="<<report_CA_rmsd_<<" impose constraints "<<impose_constraints()<<std::endl;
}

void
IdealizeMover::ignore_residues_in_csts( utility::vector1< core::Size > const i ){
	ignore_residues_in_csts_ = i;
}

utility::vector1< core::Size >
IdealizeMover::ignore_residues_in_csts() const{
	return ignore_residues_in_csts_;
}

} // namespace idealize
} // namespace protocols
