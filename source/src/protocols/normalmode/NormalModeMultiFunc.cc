// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/normalmode/AtomTreeMultifunc.hh
/// @brief  Atom tree multifunction class
/// @author Phil Bradley

/// Unit headers
#include <protocols/normalmode/NormalModeMultiFunc.hh>

/// Package headers
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/atom_tree_minimize.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>
#include <core/optimization/DOF_Node.hh>
#include <core/optimization/types.hh>

/// Project headers
#include <basic/prof.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/id/types.hh>
#include <core/id/TorsionID.hh>

#include <core/scoring/ScoringManager.hh>
//#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/Energies.hh>


#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>

/// Utility headers
#include <utility/string_util.hh>

#include <utility/vector1.hh>
#include <cmath>

static basic::Tracer TR( "protocols.normalmode.NormalModeMultiFunc" );

namespace protocols {
namespace normalmode {

using namespace core;
using namespace core::optimization;

NormalModeMultifunc::~NormalModeMultifunc() = default;

NormalModeMultifunc::NormalModeMultifunc(
	pose::Pose & pose_in,
	MinimizerMap & min_map_in,
	scoring::ScoreFunction const & scorefxn_in,
	protocols::normalmode::NormalMode const & normalmode_in,
	bool const use_omega,
	bool const deriv_check_in,
	bool const deriv_check_verbose_in
) :
	pose_( pose_in ),
	min_map_( min_map_in ),
	score_function_( scorefxn_in ),
	use_omega_( use_omega ),
	pose0_( pose_in ),
	NM_( normalmode_in ),
	deriv_check_( deriv_check_in ),
	deriv_check_verbose_( deriv_check_verbose_in ),
	deriv_check_result_( /* 0 */ )
{
	// NormalMode should be TorsionalNormalMode
	debug_assert( NM_.torsion() );
	// NormalMode should be solved prior to the minimization
	debug_assert( NM_.ntor() > 0 );

	k_dampen_ = basic::options::option[
		basic::options::OptionKeys::optimization::scale_normalmode_dampen ]();

	set_default_modes( );
	TR << "Setup dofs using default modes." << std::endl;
	get_dofs_map( );
	get_dofs_for_pose0( );
}

Real
NormalModeMultifunc::operator ()( Multivec const & vars ) const {
	PROF_START( basic::FUNC );

	// change Mode variables into torsions
	Multivec dofs( vars_to_dofs( vars ) );

	/* //Just for debugging
	for( Size i = 1; i <= vars.size(); ++i ){
	std::cout << "Vars: " << i << " " << vars[i] << std::endl;
	}

	std::list< DOF_NodeOP > dof_nodes( min_map_.dof_nodes() );
	int imap = 1;
	for ( std::list< DOF_NodeOP >::const_iterator it=dof_nodes.begin(),
	it_end = dof_nodes.end(); it != it_end; ++it, ++imap ) {
	id::TorsionID const tor_id( (*it)->torsion_id() );
	}
	*/

	// min_map_ should be an atomtree min_map
	min_map_.copy_dofs_to_pose( pose_, dofs );
	Real const score( score_function_( pose_ ) );

	PROF_STOP( basic::FUNC );
	return score;
}

void
NormalModeMultifunc::dfunc( Multivec const & vars, Multivec & dE_dvars ) const
{
	PROF_START( basic::DFUNC );

	// change Mode variables into torsions
	Multivec dE_ddofs;
	Multivec dofs = vars_to_dofs( vars );

	// min_map_ should be an atomtree min_map
	// in atom_tree_minimize.cc
	atom_tree_dfunc( pose_, min_map_, score_function_, dofs, dE_ddofs );

	dE_dvars = dEddofs_to_dEdvars( dE_ddofs );

	// optional derivative checking:
	//if ( deriv_check_ ) {
	// don't do this because min_map_ dof is different from vars here
	// this will be taken care of at NormalModeMinimizer instead
	// numerical_derivative_check( min_map_, *this, vars, dE_dvars, deriv_check_result_, deriv_check_verbose_ );
	//}

	PROF_STOP( basic::DFUNC );
}

void NormalModeMultifunc::set_deriv_check_result( NumericalDerivCheckResultOP deriv_check_result )
{
	deriv_check_result_ = deriv_check_result;
}

core::pose::Pose & NormalModeMultifunc::pose() const {
	return pose_;
}

void
NormalModeMultifunc::get_dofs_for_pose0()
{

	dofs_for_pose0_.resize( min_map_.dof_nodes().size() );
	min_map_.copy_dofs_from_pose( pose0_, dofs_for_pose0_ );
	/*
	std::list< DOF_NodeOP > dof_nodes( min_map_.dof_nodes() );

	Real const rad2deg( 57.29577951308232 );

	int imap = 1;
	for ( std::list< DOF_NodeOP >::const_iterator it=dof_nodes.begin(),
	it_end = dof_nodes.end(); it != it_end; ++it, ++imap ) {
	//id::TorsionID const tor_id( (*it)->torsion_id() );
	//tors_for_pose0_.push_back( pose0_.torsion( tor_id ) );
	id::DOF_ID const dof_id( (*it)->dof_id() );
	id::DOF_Type const type( (*it)->type() );

	if( type ==  id::RB1 || type == id::RB2 || type == id::RB3 ||
	type ==  id::RB4 || type == id::RB5 || type == id::RB6  ) {
	dofs_for_pose0_.push_back( pose0_.dof( dof_id ) );

	} else {
	dofs_for_pose0_.push_back( pose0_.dof( dof_id )*rad2deg ); // convert to degee
	}
	}
	*/
}

// Converting vars into torsion values
Multivec
NormalModeMultifunc::vars_to_dofs( Multivec const & vars ) const {

	Multivec dofs( dofs_for_pose0_ );
	utility::vector1< id::TorsionID > const NM_torID = NM_.get_torID();
	std::list< DOF_NodeOP > dof_nodes( min_map_.dof_nodes() );


	for ( Size i_var = 1; i_var <= vars.size(); ++i_var ) {
		//std::cout << "i_var, type: " << i_var << " " << var_type_[i_var] << std::endl;
		// Normal Mode variables: project eigenvectors
		if ( var_type_[i_var].compare("NM") == 0 ) {

			// Get index / scale for the mode in vars
			auto it = map_var_to_modeno_.find( i_var );
			Size const modeno( it->second );
			utility::vector1< Real > const eigv = NM_.get_eigvec_tor( modeno );

			Real const scale( get_modescale( modeno ) );

			// Iter over Normal mode torsions to accummulate into dofs
			for ( Size i_tor = 1; i_tor <= NM_torID.size(); ++i_tor ) {
				if ( map_NM_to_DOF_.count( i_tor ) ) {
					auto it = map_NM_to_DOF_.find( i_tor );
					Size const i_dof( it->second );
					dofs[ i_dof ] += scale * vars[ i_var ]*eigv[ i_tor ];
				}
			}

			// Non-Normal Mode variables: just copy
		} else {
			auto it = map_var_to_DOF_.find( i_var );
			Size const i_dof( it->second );
			// Overwrite values into dofs
			dofs[ i_dof ] = vars[ i_var ];
		}
	}

	return dofs;
} // end vars_to_dofs

Multivec
NormalModeMultifunc::dofs_to_vars( Multivec const & dofs ) const
{
	// First store the difference in torsion from the starting values
	Multivec tors_NM( NM_.ntor(), 0.0 );

	for ( Size i_dof = 1; i_dof <= dofs.size(); ++i_dof ) {
		if ( map_DOF_to_NM_.count( i_dof ) ) {
			auto it = map_DOF_to_NM_.find( i_dof );
			Size const i_tor( it->second );
			tors_NM[ i_tor ] = dofs[ i_dof ] - dofs_for_pose0( i_dof );
			//std::cout << "dof! " << i_dof << " " << tors_NM[i_tor];
			//std::cout << " " << dofs[i_dof] << " " << dofs_for_pose0(i_dof) << std::endl;

		}
	}

	// Next, iter over vars to convert dofs info to vars info
	Multivec vars( nvar(), 0.0 );

	for ( Size i_var = 1; i_var <= nvar(); ++i_var ) {
		// 1. Normal Mode vars
		// Get dot product for each given mode, and add it into vars
		if ( var_type_[i_var].compare("NM") == 0 ) {

			Real dotsum( 0.0 );

			// Get index/scale for the modes
			auto it = map_var_to_modeno_.find( i_var );
			Size const modeno( it->second );
			utility::vector1< Real > const eigv = NM_.get_eigvec_tor( modeno );
			Real const scale( get_modescale( modeno ) );

			for ( Size i_tor = 1; i_tor <= tors_NM.size(); ++i_tor ) {
				dotsum += tors_NM[ i_tor ]*eigv[ i_tor ];
			}
			vars[ i_var ] = scale*dotsum;

			// 2. Non-Normal Mode vars
		} else {
			auto it = map_var_to_DOF_.find( i_var );
			Size const i_dof( it->second );
			vars[ i_var ] = dofs[ i_dof ];

		}
	}

	return vars;
} // end dof_to_vars

Multivec
NormalModeMultifunc::dEddofs_to_dEdvars( Multivec const & dEddofs ) const
{

	Multivec dEdvars( nvar(), 0.0 );

	for ( Size i_var = 1; i_var <= nvar(); ++i_var ) {
		// 1. Normal Mode vars
		// Get dot product for each given mode, and add it into vars
		if ( var_type_[i_var].compare("NM") == 0 ) {

			// Get index/scale for the modes
			auto it = map_var_to_modeno_.find( i_var );
			Size const modeno( it->second );
			utility::vector1< Real > const eigv = NM_.get_eigvec_tor( modeno );
			Real const scale( get_modescale( modeno ) );

			for ( Size i_tor = 1; i_tor <= eigv.size(); ++i_tor ) {
				if ( map_NM_to_DOF_.count( i_tor ) ) {
					auto it = map_NM_to_DOF_.find( i_tor );
					Size const i_dof( it->second );

					//std::cout << "NM map, itor/idof: " << i_tor << " " << i_dof;
					//std::cout << " " << scale * dEddofs[i_dof]*eigv[i_tor] << std::endl;

					dEdvars[ i_var ] += scale * dEddofs[i_dof]*eigv[i_tor];
				}
			}

			// 2. Non-Normal Mode vars
		} else {
			auto it = map_var_to_DOF_.find( i_var );
			Size const i_dof( it->second );
			dEdvars[ i_var ] = dEddofs[ i_dof ];
			//std::cout << "non-NM map, ivar/idof: " << i_var << " " << i_dof;
			//std::cout << " " << dEddofs[ i_dof ] << std::endl;
		}
	}

	return dEdvars;
}

void
NormalModeMultifunc::get_dofs_map()
{
	std::list< DOF_NodeOP > dof_nodes( min_map_.dof_nodes() );
	utility::vector1< id::TorsionID > const NM_torID = NM_.get_torID();

	// Iterate loop over three times independently not to mess up the order
	// 1. RB
	// 2. NormalMode
	// 3. Omega (optional)

	// Map between NormalMode torsion index and DOF
	map_NM_to_DOF_.clear();
	map_DOF_to_NM_.clear();

	// Map between var index and DOF
	map_var_to_DOF_.clear();

	// Map between var index and NormalMode number
	map_var_to_modeno_.clear();

	var_type_.resize( 0 );

	nvar_ = 0;

	// 1. Rigid-Body DOF
	int imap = 1;
	for ( auto it=dof_nodes.begin(),
			it_end = dof_nodes.end(); it != it_end; ++it, ++imap ) {
		id::DOF_Type const type( (*it)->type() );

		//std::cout << "check: " << imap << " " << dof.rsd() << " " << dof.type() << std::endl;

		if ( type ==  id::RB1 || type == id::RB2 || type == id::RB3 ||
				type ==  id::RB4 || type == id::RB5 || type == id::RB6  ) {

			nvar_++;
			map_var_to_DOF_[ nvar_ ] = imap;
			var_type_.push_back("RB");
			TR.Debug << "Map(ivar/imap/type): " << std::setw(4) << nvar_ << " " << std::setw(4) << imap;
			TR.Debug << " RB " << std::endl;
		}
	}

	// 2-1. Store indices for Normal mode numbers
	for ( Size i_mode = 1; i_mode <= modes_using_.size(); ++i_mode ) {
		Size const modeno( modes_using_[i_mode] );
		nvar_++;
		var_type_.push_back("NM");
		map_var_to_modeno_[ nvar_ ] = modeno;
		TR.Debug << "Map(ivar/imap/type): " << std::setw(4) << nvar_ << " " << std::setw(4) << imap;
		TR.Debug << " NM_" << i_mode << std::endl;
	}

	// 2-2. Phi/psis being used by NormalMode
	// In this case, iteration is independent to nvar_
	imap = 1;
	for ( std::list< DOF_NodeOP >::const_iterator it=dof_nodes.begin(),
			it_end = dof_nodes.end(); it != it_end; ++it, ++imap ) {
		id::TorsionID const tor_id( (*it)->torsion_id() );

		for ( Size itor = 1; itor <= NM_torID.size(); ++itor ) {
			id::TorsionID const &id( NM_torID[itor] );

			if ( id == tor_id ) {
				map_NM_to_DOF_[ itor ] = imap;
				map_DOF_to_NM_[ imap ] = itor;
				break;
			}
		}
	}

	// 3. (Optional) Omega angles
	if ( use_omega_ ) {
		imap = 1;
		for ( std::list< DOF_NodeOP >::const_iterator it=dof_nodes.begin(),
				it_end = dof_nodes.end(); it != it_end; ++it, ++imap ) {
			id::TorsionID const tor_id( (*it)->torsion_id() );
			if ( tor_id.torsion() == id::omega_torsion ) {
				nvar_++;
				map_var_to_DOF_[ nvar_ ] = imap;
				var_type_.push_back("omega");

				TR.Debug << "Map(ivar/imap/type): " << std::setw(4) << nvar_ << " " << std::setw(4) << imap;
				TR.Debug << " omega_" << tor_id.rsd() << std::endl;

			}
		}
	} else {
		TR << "Passing omega dofs..." << std::endl;
	}
	TR << "dof setup completed." << std::endl;

}

/// @details Useful debugging code that can be re-enabled by changing the boolean
/// variables at the top of this function.
void
NormalModeMultifunc::dump( Multivec const & vars, Multivec const & vars2 ) const {
	// AMW: cppcheck doesn't like this, but keep stuff here for convenience
	bool debug_inaccurateG = false; // disables everything below
	bool check_score_components = true;
	bool check_score_components_verbose = false;
	bool check_rama = false;
	bool check_hbonds = true;
	//bool check_nblist = true;

	if ( ! debug_inaccurateG ) return;

	static int count_dumped( 0 );
	static bool after( true ); // dump two poses, a before and an after. Note, dumping poses as pdbs is often useless.

	if ( after ) { ++count_dumped; after = false; }
	else { after = true; }

	pose::Pose pose1( pose_ );
	pose::Pose pose2( pose_ );
	min_map_.copy_dofs_to_pose( pose1, vars );
	min_map_.copy_dofs_to_pose( pose2, vars2 );

	min_map_.copy_dofs_to_pose( pose_, vars );
	Real score_vars( score_function_( pose_ ) );

	min_map_.copy_dofs_to_pose( pose_, vars2 );
	Real score_vars2( score_function_( pose_ ) );

	Real alt_score_vars = score_function_( pose1 );
	pose1.dump_pdb( "atomtree_multifunc_error_pose_before" + utility::to_string( count_dumped  ) + ".pdb" );

	Real alt_score_vars2 = score_function_( pose2 );
	pose2.dump_pdb( "atomtree_multifunc_error_pose_after" + utility::to_string( count_dumped  ) + ".pdb" );

	std::cerr << "starting pose energies: " << score_vars << " vs " << alt_score_vars << std::endl;
	pose1.energies().total_energies().show_weighted( std::cerr, score_function_.weights() );
	std::cerr << std::endl;
	std::cerr << "moved pose energies: " << score_vars2 << " vs " << alt_score_vars2 << std::endl;
	pose2.energies().total_energies().show_weighted( std::cerr, score_function_.weights() );
	std::cerr << std::endl;
	using namespace scoring;

	if ( check_score_components ) {
		// slow! Iterate through all the components and check their derivatives one by one.
		const_cast< bool & > (deriv_check_) = true;
		if ( check_score_components_verbose ) {
			const_cast< bool & > (deriv_check_verbose_) = true;
		}
		Multivec dvars( vars );
		scoring::EnergyMap orig_weights( score_function_.weights() );
		for ( Size ii = 1; ii <= scoring::n_score_types; ++ii ) {
			using namespace scoring;

			if ( score_function_.weights()[ (ScoreType ) ii ] == 0.0 ) continue;

			for ( Size jj = 1; jj <= scoring::n_score_types; ++jj ) {
				if ( jj == ii ) {
					const_cast< scoring::ScoreFunction & > (score_function_).set_weight( (scoring::ScoreType) jj, orig_weights[ (scoring::ScoreType) jj ]);
				} else if ( score_function_.weights()[ (scoring::ScoreType ) jj ] != 0.0 ) {
					const_cast< scoring::ScoreFunction & > (score_function_).set_weight( (scoring::ScoreType) jj, 1e-9 );
				}
			}
			//std::cout << "Checking score type: " << scoring::ScoreType( ii ) << std::endl;
			dfunc( vars, dvars ); // invokes numeric derivative checker.
		}
		for ( Size ii = 1; ii <= scoring::n_score_types; ++ii ) {
			if ( orig_weights[ scoring::ScoreType( ii ) ] != 0 ) {
				const_cast< scoring::ScoreFunction & > (score_function_).set_weight( (scoring::ScoreType)ii, orig_weights[ (scoring::ScoreType) ii ]);
			}
		}
		const_cast< bool & > (deriv_check_) = false;
		const_cast< bool & > (deriv_check_verbose_) = false;
	}

	if ( check_rama ) {
		// useful if rama seems to be the culprit.  This is the only piece of code
		// that invokes eval_rama_score_all, so that function may be hacked to provide
		// clearer debugging output.
		ScoringManager::get_instance()->get_Ramachandran().eval_rama_score_all( pose1, score_function_ );
	}

	if ( check_hbonds ) {
		scoring::hbonds::HBondSet hbond_set;
		fill_hbond_set( pose1, true, hbond_set );

		for ( Size ii = 1; ii <= hbond_set.nhbonds(); ++ii ) {
			scoring::hbonds::HBond const & iihbond = hbond_set.hbond( ii );
			std::cerr << "Hbond " << ii <<
				" d: " << iihbond.don_res() << " " << iihbond.don_hatm() <<
				" a: " << iihbond.acc_res() << " " << iihbond.acc_atm() <<
				" e: " << iihbond.energy() << " " << iihbond.weight() <<
				" good? " << hbond_set.allow_hbond( iihbond ) <<
				/*" f1: " << iihbond.deriv().first.x() <<
				" " << iihbond.deriv().first.y() <<
				" " << iihbond.deriv().first.z() <<
				" f2: " << iihbond.deriv().second.x() <<
				" " << iihbond.deriv().second.y() <<
				" " << iihbond.deriv().second.z() <<*/ std::endl;
		}
	}
}

// get dampening scaling factor given mode number
Real
NormalModeMultifunc::get_modescale( Size const modeno ) const
{
	Real scale = exp( -k_dampen_*( modeno - 1.0 ));
	return scale;
}

// Set default modes
void
NormalModeMultifunc::set_default_modes()
{
	modes_using_.resize( 0 );
	Size ntors = std::min( (Size)(50), NM_.ntor() );

	for ( core::Size i = 1; i <= ntors; ++i ) {
		modes_using_.push_back( i );
	}
}

void
NormalModeMultifunc::set_modes( utility::vector1< Size > modes_using_in )
{
	modes_using_ = modes_using_in;
	// Set dofs again
	TR << "Reset dofs based on user-defined modes." << std::endl;
	get_dofs_map( );
	get_dofs_for_pose0( );

}

MinimizerMap const & NormalModeMultifunc::min_map() const {
	return min_map_;
}

core::scoring::ScoreFunction const & NormalModeMultifunc::score_function() const {
	return score_function_;
}

} // namespace optimization
} // namespace core

