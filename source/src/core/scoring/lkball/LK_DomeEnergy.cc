// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/scoring/lkball/LK_DomeEnergy.cc
/// @brief  lk_dome second shell water interactions
/// @author Brian Coventry (bcov@uw.edu)

// Unit Headers
#include <core/scoring/lkball/LK_DomeEnergy.hh>
#include <core/scoring/lkball/LK_DomeEnergyCreator.hh>


// Package headers

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/MinimizationData.hh>


#include <core/conformation/residue_datacache.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/etable/Etable.hh>

// Basic Headers
#include <basic/datacache/BasicDataCache.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>

#include <numeric/conversions.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>
#include <numeric/interpolation/spline/SimpleInterpolator.hh>
#include <numeric/numeric.functions.hh>

#include <core/scoring/ResidueNeighborList.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>
#include <utility/graph/Graph.hh>


static basic::Tracer TR("core.scoring.methods.LK_DomeEnergy");

using core::Real;

namespace core {
namespace scoring {
namespace lkball {



class LK_DomeInvoker : public etable::count_pair::Invoker {
public:
	LK_DomeInvoker(
		LK_DomeEnergy const & lk_dome,
		conformation::Residue const & rsd1,
		LKD_ResidueInfo const & rsd1_info,
		conformation::Residue const & rsd2,
		LKD_ResidueInfo const & rsd2_info,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	);

	~LK_DomeInvoker() override = default;

protected:
	inline LK_DomeEnergy const & lk_dome() const;
	inline core::conformation::Residue const & rsd1() const;
	inline LKD_ResidueInfo const & rsd1_info() const;
	inline conformation::Residue const & rsd2() const;
	inline LKD_ResidueInfo const & rsd2_info() const;
	inline ScoreFunction const & sfxn() const;
	inline EnergyMap & emap() const;


private:
	LK_DomeEnergy const & lk_dome_;
	core::conformation::Residue const & rsd1_;
	LKD_ResidueInfo const & rsd1_info_;
	conformation::Residue const & rsd2_;
	LKD_ResidueInfo const & rsd2_info_;
	ScoreFunction const & sfxn_;
	EnergyMap & emap_;
};

class LK_Dome_RPE_Invoker : public LK_DomeInvoker {
public:
	LK_Dome_RPE_Invoker(
		LK_DomeEnergy const & lk_dome,
		conformation::Residue const & rsd1,
		LKD_ResidueInfo const & rsd1_info,
		conformation::Residue const & rsd2,
		LKD_ResidueInfo const & rsd2_info,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	);

	~LK_Dome_RPE_Invoker() override = default;

	void invoke( etable::count_pair::CountPairFunction const & cp ) override;
};


LK_DomeInvoker::LK_DomeInvoker(
	LK_DomeEnergy const & lk_dome,
	conformation::Residue const & rsd1,
	LKD_ResidueInfo const & rsd1_info,
	conformation::Residue const & rsd2,
	LKD_ResidueInfo const & rsd2_info,
	ScoreFunction const & sfxn,
	EnergyMap & emap_in
) :
	etable::count_pair::Invoker(),
	lk_dome_( lk_dome ),
	rsd1_( rsd1 ),
	rsd1_info_( rsd1_info ),
	rsd2_( rsd2 ),
	rsd2_info_( rsd2_info ),
	sfxn_( sfxn ),
	emap_( emap_in )
{}


LK_DomeEnergy const &
LK_DomeInvoker::lk_dome() const
{
	return lk_dome_;
}

core::conformation::Residue const &
LK_DomeInvoker::rsd1() const
{
	return rsd1_;
}

LKD_ResidueInfo const &
LK_DomeInvoker::rsd1_info() const
{
	return rsd1_info_;
}

conformation::Residue const &
LK_DomeInvoker::rsd2() const
{
	return rsd2_;
}

LKD_ResidueInfo const &
LK_DomeInvoker::rsd2_info() const
{
	return rsd2_info_;
}

ScoreFunction const &
LK_DomeInvoker::sfxn() const
{
	return sfxn_;
}

EnergyMap &
LK_DomeInvoker::emap() const
{
	return emap_;
}






Distance          // Just assume the worst that max_angle is 180
// lk_dome  Atom 2.9 LK_ball 2.65 LK_dome w_rad(1.4) fade(1.92) radius(2.2) Atom
// lk_dome = 11.07
// lk_dome_br  Atom 2.9 LK_ball 2.65 LK_dome fade(2.23) LK_ball 2.9 Atom
// lk_dome_br = 10.68
LK_DomeEnergy::atomic_interaction_cutoff() const { return packing_ ? 0 : 11.07 - 2.65 + w_dist_; }


Real                                            // only true if the above assumes max_angle is 180
LK_DomeEnergy::water_atom_interaction_cutoff() const { return 11.07-2.65 - 2.65 + w_dist_; }

Real                                            // only true if the above assumes max_angle is 180
LK_DomeEnergy::water_atom_interaction_cutoff2() const { return (11.07-2.65 - 2.65 + w_dist_)*(11.07-2.65 - 2.65 + w_dist_); }

Real                                            // only true if the above assumes max_angle is 180
LK_DomeEnergy::water_water_bridge_cutoff2() const { return (11.07-2.65*2 - 2.65 + w_dist_)*(11.07-2.65*2 - 2.65 + w_dist_); }


core::Size
LK_DomeEnergy::version() const { return 1; }


LK_BallEnergyCOP
LK_DomeEnergy::lk_ball() const { return lk_ball_; }

Real
LK_DomeEnergy::occlusion_max() const { return occlusion_max_; }

Real
LK_DomeEnergy::occlusion_min() const { return occlusion_min_; }

Real
LK_DomeEnergy::w_dist() const { return w_dist_; }

Real
LK_DomeEnergy::water_adjust() const { return water_adjust_; }






scoring::methods::EnergyMethodOP
LK_DomeEnergyCreator::create_energy_method(
	scoring::methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< LK_DomeEnergy >( options );
}

scoring::ScoreTypes
LK_DomeEnergyCreator::score_types_for_method() const {
	scoring::ScoreTypes sts;
	sts.push_back( lk_dome );
	sts.push_back( lk_dome_iso );
	sts.push_back( lk_dome_bridge );
	sts.push_back( lk_dome_bridge_uncpl );
	sts.push_back( lk_ball_bridge2 );
	sts.push_back( lk_ball_bridge_uncpl2 );
	return sts;
}

LK_DomeEnergy::LK_DomeEnergy( core::scoring::methods::EnergyMethodOptions const & options ) :
	scoring::methods::ContextDependentTwoBodyEnergy( utility::pointer::make_shared< LK_DomeEnergyCreator >() ),
	dump_dome_waters_( false ),
	dump_dome_bridge_waters_( false ),
	debug_disable_count_pair_( false ),
	packing_( false ),
	minimizing_( false )
{

	lk_ball_ = utility::pointer::make_shared< LK_BallEnergy >( options );

	Real max_angle = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_dome_max_angle ]();
	Real min_angle = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_dome_min_angle ]();

	if ( max_angle > 180 ) {
		utility_exit_with_message("lk_dome_max_angle must be less than 180.");
	}
	if ( min_angle < 0 ) {
		utility_exit_with_message("lk_dome_min_angle must be greater than 0.");
	}
	if ( min_angle > max_angle ) {
		utility_exit_with_message("lk_dome_min_angle less than or equal to lk_dome_min_angle.");
	}

	max_angle_cos_ = std::cos( numeric::conversions::radians( 180 - max_angle ) );
	// max_angle_cos_ = 5;
	max_angle_sin_ = std::sin( numeric::conversions::radians( 180 - max_angle ) );

	min_angle_cos_ = std::cos( numeric::conversions::radians( 180 - min_angle ) );
	// min_angle_cos_ = -5;
	min_angle_sin_ = std::sin( numeric::conversions::radians( 180 - min_angle ) );

	h2o_radius_ = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_dome_h2o_radius ]();    // 1.4
	ramp_width_A2_ = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_dome_ramp_width_A2 ](); // 3.9
	overlap_width_A2_ = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_dome_overlap_width_A2 ](); // 5.0
	ball_overlap_width_A2_ = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_dome_ball_overlap_width_A2 ](); // 5.0
	occlusion_max_ = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_dome_occlusion_max ]();
	occlusion_min_ = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_dome_occlusion_min ]();
	w_dist_ = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_dome_water_dist ]();
	water_adjust_ = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_dome_water_adjust ]();

	if ( occlusion_min_ >= occlusion_max_ ) {
		utility_exit_with_message("lk_dome_occlusion_max must be greater than occlusion_min.");
	}

	setup_d2_bounds();
	setup_poly_params();

}


void
LK_DomeEnergy::setup_d2_bounds()
{
	Real const h2o_radius( h2o_radius_ );
	chemical::AtomTypeSet const & atom_set( *(lk_ball_->etable().atom_set().lock() ) );
	d2_low_.resize( atom_set.n_atomtypes() );
	for ( Size i=1; i<= atom_set.n_atomtypes(); ++i ) {
		chemical::AtomType const & atype( atom_set[ i ] );
		Real const d2_high( numeric::square( h2o_radius + atype.lj_radius() ) ); // was 3.0 * 3.0
		d2_low_[ i ] = std::max( 0.0, d2_high - ramp_width_A2_ );
	}
}


scoring::methods::EnergyMethodOP
LK_DomeEnergy::clone() const {
	return utility::pointer::make_shared< LK_DomeEnergy >( *this );
}


void
LK_DomeEnergy::setup_poly_params() {
	numeric::interpolation::spline::SplineGenerator gen(
		5, lk_frac_sigmoid(5), dlk_frac_sigmoid(5),
		6, 0, 0);
	numeric::interpolation::spline::InterpolatorOP interp = gen.get_interpolator();
	numeric::interpolation::spline::SimpleInterpolatorOP sinterp =
		utility::pointer::dynamic_pointer_cast< numeric::interpolation::spline::SimpleInterpolator > ( interp );
	runtime_assert(sinterp);

	numeric::SplineParameters sparams;
	sparams.ylo  = sinterp->y()[ 1 ];
	sparams.yhi  = sinterp->y()[ 2 ];
	sparams.y2lo = sinterp->ddy()[ 1 ];
	sparams.y2hi = sinterp->ddy()[ 2 ];

	poly_params_ = numeric::cubic_polynomial_from_spline( 5, 6, sparams );
}

Real
LK_DomeEnergy::dlk_frac_sigmoid( Real distance ) const {
	return -(2 * std::exp(2 * (distance - 4.2)))/numeric::square(std::exp(2 * (distance - 4.2)) + 1);
}

Real
LK_DomeEnergy::dmy_lk_fraction( Real distance ) const {
	// if ( distance < 1000000 ) return 0;
	if ( distance > 6 ) {
		return 0;
	}
	if ( distance < 5 ) {
		return dlk_frac_sigmoid( distance );
	}
	return numeric::cubic_polynomial_deriv( distance, poly_params_ );
}

Real
LK_DomeEnergy::lk_frac_sigmoid( Real distance ) const {
	return 1.0/ ( 1.0 + std::exp(  (distance - 4.2) * 2));
}

// This should be the functional form of fa_sol basically
Real
LK_DomeEnergy::my_lk_fraction( Real distance ) const {
	// if ( distance < 10000 ) return 1;
	if ( distance > 6 ) {
		return 0;
	}
	if ( distance < 5 ) {
		return lk_frac_sigmoid( distance );
	}
	return numeric::eval_cubic_polynomial( distance, poly_params_ );
}

Real
LK_DomeEnergy::get_sol_value( chemical::AtomType const & at, Size nattached_waters ) const {
	if ( nattached_waters == 0 ) nattached_waters = 1;
	return at.lk_dgfree() * 0.01 / nattached_waters;
}


void
evaluate_lk_dome_energy_for_atom_ranges(
	LK_DomeEnergy const & lk_dome,
	conformation::Residue const & rsd1,
	LKD_ResidueInfo const & rsd1_info,
	conformation::Residue const & rsd2,
	LKD_ResidueInfo const & rsd2_info,
	ScoreFunction const & sfxn,
	etable::count_pair::CountPairFunction const & cpfxn,
	Size const res1_start_atom,
	Size const res1_end_atom,
	Size const res2_start_atom,
	Size const res2_end_atom,
	EnergyMap & emap
)
{
	(void)sfxn;

	WaterCoords const & rsd1_waters( rsd1_info.waters() );
	WaterCoords const & rsd2_waters( rsd2_info.waters() );

	// utility::vector1< AtomWeights > const & rsd1_atom_wts( rsd1_info.atom_weights() );
	// utility::vector1< AtomWeights > const & rsd2_atom_wts( rsd2_info.atom_weights() );

	utility::vector1< WaterOcclusions > const & rsd1_water_occlusions( rsd1_info.water_occlusions() );
	utility::vector1< WaterOcclusions > const & rsd2_water_occlusions( rsd2_info.water_occlusions() );

	utility::vector1< Real > const & rsd1_water_sol_values( rsd1_info.water_sol_values() );
	utility::vector1< Real > const & rsd2_water_sol_values( rsd2_info.water_sol_values() );


	// bool use_lkbr = (sfxn.get_weight( core::scoring::lk_ball_bridge )!=0);
	// bool use_lkbr_uncpl = (sfxn.get_weight( core::scoring::lk_ball_bridge_uncpl)!=0);

	// setup residue information
	for ( Size atom1 = res1_start_atom; atom1 <= res1_end_atom; ++atom1 ) {
		Size const n_atom1_waters = rsd1_info.n_attached_waters()[ atom1 ];
		// WaterCoords const & atom1_waters( rsd1_waters[ atom1 ] );
		WaterOcclusions const & atom1_occlusions( rsd1_water_occlusions[ atom1 ] );
		Real const & atom1_sol_value( rsd1_water_sol_values[ atom1 ] );
		Vector const & atom1_xyz( rsd1.xyz( atom1 ) );
		Size const atom1_type_index( rsd1.atom( atom1 ).type() );
		Size rsd1_offset = rsd1_info.water_offset_for_atom()[atom1];
		// AtomWeights const & atom1_weights( rsd1_atom_wts[atom1] );

		for ( Size atom2 = res2_start_atom; atom2 <= res2_end_atom; ++atom2 ) {
			Vector const & atom2_xyz( rsd2.xyz( atom2 ) );
			Real const d2( atom1_xyz.distance_squared( atom2_xyz ) );
			if ( d2 >= lk_dome.atomic_interaction_cutoff() * lk_dome.atomic_interaction_cutoff() ) continue;

			Size const n_atom2_waters = rsd2_info.n_attached_waters()[ atom2 ];
			if ( n_atom1_waters == 0 && n_atom2_waters == 0 ) continue;

			Real cp_weight = 1.0;
			Size pathdist;
			if ( ! lk_dome.debug_disable_count_pair_ ) {
				if ( ! cpfxn.count( atom1, atom2, cp_weight, pathdist ) ) continue;
			}

			// WaterCoords const & atom2_waters( rsd2_waters[ atom2 ] );
			WaterOcclusions const & atom2_occlusions( rsd2_water_occlusions[ atom2 ] );
			Real const & atom2_sol_value( rsd2_water_sol_values[ atom2 ] );
			Size const atom2_type_index( rsd2.atom( atom2 ).type() );
			Size rsd2_offset = rsd2_info.water_offset_for_atom()[atom2];;
			// AtomWeights const & atom2_weights( rsd2_atom_wts[atom2] );




			if ( d2 == Real(0.0) ) continue; // sanity check

			// Real lk_desolvation_of_atom1_by_atom2 = 0.0, lk_desolvation_of_atom2_by_atom1 = 0.0;

			// Replace these with the value of the lk_ball water
			Real lk_desolvation_of_atom1_by_water = cp_weight * atom1_sol_value;
			Real lk_desolvation_of_atom2_by_water = cp_weight * atom2_sol_value;

			bool save_debug = lk_dome.dump_dome_waters_;

			if (  atom2 != rsd2.nheavyatoms() ) lk_dome.dump_dome_waters_ = false;

			// if ( lk_dome.dump_dome_waters_ ) {
			//     std::cout << "DUMPING!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
			// }

			// std::cout << rsd2.seqpos() << " " << rsd2.atom_name(atom2) << " -> " << rsd1.seqpos() << " " << rsd1.atom_name(atom1) << std::endl;
			// std::cout << "K " << rsd1.seqpos() << " " << rsd2.seqpos() << " " << atom1 << " " << atom2 << std::endl;
			lk_dome.accumulate_single_atom_contributions( atom1, atom1_type_index, n_atom1_waters, rsd1_offset, rsd1_waters, atom1_occlusions,
				rsd1, atom2_type_index, atom2_xyz,
				lk_desolvation_of_atom1_by_water, emap );
			lk_dome.dump_dome_waters_ = false;
			// std::cout << rsd1.seqpos() << " " << rsd1.atom_name(atom1) << " -> " << rsd2.seqpos() << " " << rsd2.atom_name(atom2) << std::endl;

			// std::cout << "K " << rsd2.seqpos() << " " << rsd1.seqpos() << " " << atom2 << " " << atom1 << std::endl;
			lk_dome.accumulate_single_atom_contributions( atom2, atom2_type_index, n_atom2_waters, rsd2_offset, rsd2_waters, atom2_occlusions,
				rsd2, atom1_type_index, atom1_xyz,
				lk_desolvation_of_atom2_by_water, emap );

			lk_dome.dump_dome_waters_ = save_debug;

			// if ( !use_lkbr_uncpl && !use_lkbr ) continue;

			// // fpd - get lk_ball_bridge
			if ( n_atom1_waters != 0 && n_atom2_waters != 0 ) {
				core::Real bridge_frac1 = lk_dome.get_lkd_bridge_fractional_1way(
					atom1_xyz,
					n_atom1_waters, rsd1_offset, n_atom2_waters, rsd2_offset,
					rsd1_waters, rsd2_waters,
					atom1_occlusions );

				core::Real bridge_frac2 = lk_dome.get_lkd_bridge_fractional_1way(
					atom2_xyz,
					n_atom2_waters, rsd2_offset, n_atom1_waters, rsd1_offset,
					rsd2_waters, rsd1_waters,
					atom2_occlusions );

				emap[ lk_dome_bridge ] += bridge_frac1 * lk_desolvation_of_atom1_by_water + bridge_frac2 * lk_desolvation_of_atom2_by_water;
				emap[ lk_dome_bridge_uncpl ] += cp_weight * ( bridge_frac1 + bridge_frac2 );


				core::Real ball_bridge_frac1 = lk_dome.get_lkb_bridge2_fractional_1way(

					n_atom1_waters, rsd1_offset, n_atom2_waters, rsd2_offset,
					rsd1_waters, rsd2_waters,
					atom1_occlusions );

				core::Real ball_bridge_frac2 = lk_dome.get_lkb_bridge2_fractional_1way(

					n_atom2_waters, rsd2_offset, n_atom1_waters, rsd1_offset,
					rsd2_waters, rsd1_waters,
					atom2_occlusions );

				emap[ lk_ball_bridge2 ] += ball_bridge_frac1 * lk_desolvation_of_atom1_by_water + ball_bridge_frac2 * lk_desolvation_of_atom2_by_water;
				emap[ lk_ball_bridge_uncpl2 ] += cp_weight * ( ball_bridge_frac1 + ball_bridge_frac2 );
			}
		} // atom2
	} // atom1

	//PROF_STOP( basic::LK_BALL_RESIDUE_PAIR_ENERGY );

}


LK_Dome_RPE_Invoker::LK_Dome_RPE_Invoker(
	LK_DomeEnergy const & lk_dome,
	conformation::Residue const & rsd1,
	LKD_ResidueInfo const & rsd1_info,
	conformation::Residue const & rsd2,
	LKD_ResidueInfo const & rsd2_info,
	ScoreFunction const & sfxn,
	EnergyMap & emap_in
) :
	LK_DomeInvoker( lk_dome, rsd1, rsd1_info, rsd2, rsd2_info, sfxn, emap_in )
{}

void LK_Dome_RPE_Invoker::invoke( etable::count_pair::CountPairFunction const & cpfxn )
{
	evaluate_lk_dome_energy_for_atom_ranges( lk_dome(), rsd1(), rsd1_info(), rsd2(), rsd2_info(), sfxn(), cpfxn,
		1, rsd1().nheavyatoms(), 1, rsd2().nheavyatoms(), emap() );
}




void
LK_DomeEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const & sf,
	EnergyMap & emap
) const {
	if ( packing_ ) return;
	runtime_assert( ! packing_ );

	using conformation::residue_datacache::LK_DOME_INFO;
	debug_assert( dynamic_cast< LKD_ResidueInfo const * >( rsd1.data().get_raw_const_ptr( LK_DOME_INFO ) ));
	debug_assert( dynamic_cast< LKD_ResidueInfo const * >( rsd2.data().get_raw_const_ptr( LK_DOME_INFO ) ));

	residue_pair_energy( rsd1,
		*( static_cast< LKD_ResidueInfo const * >( rsd1.data().get_raw_const_ptr( LK_DOME_INFO ))),
		rsd2,
		*( static_cast< LKD_ResidueInfo const * >( rsd2.data().get_raw_const_ptr( LK_DOME_INFO ))),
		sf, emap );
}


void
LK_DomeEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	LKD_ResidueInfo const & rsd1_info,
	conformation::Residue const & rsd2,
	LKD_ResidueInfo const & rsd2_info,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	//PROF_START( basic::LK_BALL_RESIDUE_PAIR_ENERGY );

	using namespace etable::count_pair;
	//bool const verbose( false );

	LK_Dome_RPE_Invoker invoker(
		*this, rsd1, rsd1_info,
		rsd2, rsd2_info,
		sfxn, emap);

	CPCrossoverBehavior crossover = (rsd1.is_polymer_bonded(rsd2) && rsd2.is_polymer_bonded(rsd1))? CP_CROSSOVER_4 : CP_CROSSOVER_3;

	CountPairFactory::create_count_pair_function_and_invoke( rsd1, rsd2, crossover, invoker );
}







class LKD_ResPairMinData : public basic::datacache::CacheableData {
public:
	LKD_ResPairMinData();
	LKD_ResPairMinData( LKD_ResPairMinData const & src );
	~LKD_ResPairMinData() override = default;
	basic::datacache::CacheableDataOP clone() const override
	{ return utility::pointer::make_shared<LKD_ResPairMinData>( *this ) ; }

	void
	initialize(
		LKD_ResidueInfoCOP res1_data,
		LKD_ResidueInfoCOP res2_data
	);

	LKD_ResidueInfo const & res1_data() const { return *res1_data_; }
	LKD_ResidueInfo const & res2_data() const { return *res2_data_; }

	bool
	initialized() const { return initialized_; }

private:

	LKD_ResidueInfoCOP res1_data_;
	LKD_ResidueInfoCOP res2_data_;

	bool initialized_{ false };
};

using LKD_ResPairMinDataOP = utility::pointer::shared_ptr<LKD_ResPairMinData>;
using LKD_ResPairMinDataCOP = utility::pointer::shared_ptr<const LKD_ResPairMinData>;

LKD_ResPairMinData::LKD_ResPairMinData()= default;

LKD_ResPairMinData::LKD_ResPairMinData( LKD_ResPairMinData const & src ) :
	basic::datacache::CacheableData( src ),
	res1_data_( src.res1_data_ ),
	res2_data_( src.res2_data_ ),
	initialized_( src.initialized_ )
{
}


void
LKD_ResPairMinData::initialize(
	LKD_ResidueInfoCOP res1_data,
	LKD_ResidueInfoCOP res2_data
)
{
	initialized_ = true;
	res1_data_ = res1_data;
	res2_data_ = res2_data;
}







Real
LK_DomeEnergy::get_avail( Real occl ) const {
	if ( occl < occlusion_min_ ) return 1;
	if ( occl > occlusion_max_ ) return 0;

	Real span = occlusion_max_ - occlusion_min_;
	return ( occlusion_max_ - occl) / span;
}



Real
LK_DomeEnergy::eval_d_lk_fraction_dr_over_r( Real const d2_delta, Real const width ) const
{
	// if ( d2_delta < 10000 ) return 0;
	debug_assert( d2_delta >= -0.001 && d2_delta <= width + 0.001 );
	Real const inv_range( 1.0 / width );
	Real const xprime( inv_range * d2_delta );
	return -8.0 * inv_range * ( 1 - xprime * xprime ) * xprime;  //?? why 8 (and not 4) // because derivative of d2_delta brings a 2
}


// Stolen from LK_Ball
Real
LK_DomeEnergy::eval_lk_fraction( Real const d2_delta, Real const width ) const
{
	debug_assert( d2_delta >= -0.001 && d2_delta <= width + 0.001 );
	Real const inv_range( 1.0 / width );
	Real const xprime( inv_range * d2_delta );
	return ( 1 - xprime*xprime ) * ( 1 - xprime*xprime );
}

// Let d = dome_dist. l = d2_low

// delta = d**2 - l
// xprime = delta / ramp
// E = (1 - xprime**2 )**2

// dxprime = 1/ramp*ddelta
// ddelta = 2d*dd
// dxprime = 2d/ramp*dd

//
// dE = 2(1 - xprime**2)*(-2xprime)dxprime
// dE = (-4xprime + 4xprime**3)*dxprime
// dE = (-4xprime * 4xprime**3)*2d/ramp*dd
// dE/dd = -8 * ( 1 - xprime**2)*xprime * d / ramp
// dE/dd/d = -8 * ( 1 - xprime**2)*xprime / ramp

Real
LK_DomeEnergy::get_lkd_frac_dome_dist( Real dw_dist2, Real atom2_type ) const {
	// if ( dw_dist2 < 100000 ) return 1;

	Real const d2_low( d2_low_[ atom2_type ] );

	// find the closest water:
	Real const d2_delta = dw_dist2 - d2_low;
	if ( d2_delta > ramp_width_A2_ ) return 0.0;
	else if ( d2_delta < 0.0 ) return 1.0;
	else return eval_lk_fraction( d2_delta, ramp_width_A2_ );
}

Real
LK_DomeEnergy::single_water_water_fraction_1way(
	Vector const & base_atom_xyz,
	Vector const & water_xyz,
	Vector const & other_water_xyz
) const {

	Real _;
	Vector _2;
	Real dist2 = dome_water_dist2( base_atom_xyz, water_xyz, other_water_xyz, _, _2 );
	Real d2_delta = dist2;

	Real this_frac( 0.0 );
	if ( d2_delta < overlap_width_A2_ ) {
		this_frac = ( d2_delta < 0.0 ? Real( 1.0 ) : eval_lk_fraction( d2_delta, overlap_width_A2_ ) );
	}

	return this_frac;
}


// Hbonds resulting from dome1 overlapping with lk_water2
Real
LK_DomeEnergy::get_lkd_bridge_fractional_1way(
	Vector const & atom1_xyz,
	Size n_atom1_waters,
	Size atom1_start_water,
	Size n_atom2_waters,
	Size atom2_start_water,
	WaterCoords const & atom1_waters,
	WaterCoords const & atom2_waters,
	WaterOcclusions const & atom1_occlusions
) const {

	Real frac = 0;

	for ( Size i_water = 1; i_water <= n_atom1_waters; i_water++ ) {

		if ( atom1_occlusions[i_water] > occlusion_max_ ) continue;

		Real availability = get_avail( atom1_occlusions[i_water] );

		for ( Size j_water = 1; j_water <= n_atom2_waters; j_water++ ) {

			if ( atom1_waters[i_water + atom1_start_water].distance_squared( atom2_waters[j_water + atom2_start_water] ) > water_water_bridge_cutoff2() ) {
				debug_assert( single_water_water_fraction_1way( atom1_xyz, atom1_waters[i_water + atom1_start_water], atom2_waters[j_water + atom2_start_water] ) == 0 );
				continue;
			}

			frac += availability * single_water_water_fraction_1way( atom1_xyz, atom1_waters[i_water + atom1_start_water], atom2_waters[j_water + atom2_start_water] );

			// if ( dump_dome_bridge_waters_ && this_frac > 0.25 ) {
			//     dump_dome_waters_ = true;
			//     dome_water_dist( atom1_xyz, atom1_waters[i_water], atom2_waters[j_water] );
			//     dump_dome_waters_ = false;
			// }

		}
	}
	return frac;
}



Real
LK_DomeEnergy::single_water_water_bridge2_fraction_1way(
	Vector const & water_xyz,
	Vector const & other_water_xyz
) const {

	Real d2_delta = water_xyz.distance_squared(other_water_xyz);

	Real this_frac( 0.0 );
	if ( d2_delta < ball_overlap_width_A2_ ) {
		this_frac = ( d2_delta < 0.0 ? Real( 1.0 ) : eval_lk_fraction( d2_delta, ball_overlap_width_A2_ ) );
	}

	return this_frac;
}

// Hbonds from single overlapping waters
Real
LK_DomeEnergy::get_lkb_bridge2_fractional_1way(
	Size n_atom1_waters,
	Size atom1_start_water,
	Size n_atom2_waters,
	Size atom2_start_water,
	WaterCoords const & atom1_waters,
	WaterCoords const & atom2_waters,
	WaterOcclusions const & atom1_occlusions
) const {

	Real frac = 0;

	for ( Size i_water = 1; i_water <= n_atom1_waters; i_water++ ) {

		if ( atom1_occlusions[i_water] > occlusion_max_ ) continue;

		Real availability = get_avail( atom1_occlusions[i_water] );

		for ( Size j_water = 1; j_water <= n_atom2_waters; j_water++ ) {

			frac += availability * single_water_water_bridge2_fraction_1way( atom1_waters[i_water + atom1_start_water], atom2_waters[j_water + atom2_start_water] );

			// if ( dump_dome_bridge_waters_ && this_frac > 0.25 ) {
			//     dump_dome_waters_ = true;
			//     dome_water_dist( atom1_xyz, atom1_waters[i_water], atom2_waters[j_water] );
			//     dump_dome_waters_ = false;
			// }

		}
	}
	return frac;
}


Vector
LK_DomeEnergy::get_dome_water(
	Vector const & base_atom_xyz,
	Vector const & water,
	Vector const & other_atom_xyz,
	Real & water_atom_distance,
	Vector & wadj_water
) const
{

	Vector water_unit = (water - base_atom_xyz).normalized();

	wadj_water = water + water_unit * water_adjust_;

	// debug_assert( std::abs<Real>((water - base_atom_xyz).norm()) - 2.65 < 0.01 );
	// Vector water_unit = (water - base_atom_xyz) / 2.65;
	Vector water_to_other = other_atom_xyz - wadj_water;
	Real water_to_other_norm = water_to_other.norm();
	Vector water_to_other_unit = water_to_other / water_to_other_norm;

	water_atom_distance = water_to_other_norm;

	Real dot = water_unit.dot( water_to_other_unit );

	// debug_assert( dot <= 1.01 && dot >= -1.01);
	// dot = std::min<Real>(dot, 1 );
	// dot = std::max<Real>(dot, -1 );

	Vector water_pos;

	if ( dot < min_angle_cos_ ) {
		Vector out_of_plane = water_unit.cross( water_to_other_unit ); //.normalized();
		Vector sin_axis = out_of_plane.cross(water_unit).normalized();

		water_pos = wadj_water + sin_axis * min_angle_sin_ * w_dist_ + water_unit * min_angle_cos_ * w_dist_;

	} else if ( dot > max_angle_cos_ ) {
		Vector out_of_plane = water_unit.cross( water_to_other_unit ); //.normalized();
		Vector sin_axis = out_of_plane.cross(water_unit).normalized();

		water_pos = wadj_water + sin_axis * max_angle_sin_ * w_dist_ + water_unit * max_angle_cos_ * w_dist_;
	} else {
		water_pos = wadj_water + water_to_other_unit * w_dist_;
	}

	Real _;
	Vector & _2 = water_to_other;
	debug_assert( std::abs( water_pos.distance( other_atom_xyz ) - std::sqrt(dome_water_dist2( base_atom_xyz, water, other_atom_xyz, _, _2)) ) < 0.01 );

	return water_pos;
}

// Get the location of the water inside the dome which is always
// the closest point out at w_dist_Å and within the angular bounds
// Return distance to avoid unnecessary vector projections and stuff
Real
LK_DomeEnergy::dome_water_dist2(
	Vector const & base_atom_xyz,
	Vector const & water,
	Vector const & other_atom_xyz,
	Real & water_atom_distance,
	Vector & wadj_water
) const
{

	Vector water_unit = (water - base_atom_xyz).normalized();

	wadj_water = water + water_unit * water_adjust_;

	// debug_assert( std::abs<Real>((water - base_atom_xyz).norm()) - w_dist_ < 0.01 );
	// Vector water_unit = (water - base_atom_xyz) / w_dist_;
	Vector water_to_other = other_atom_xyz - wadj_water;
	Real water_to_other_norm = water_to_other.norm();
	Vector water_to_other_unit = water_to_other / water_to_other_norm;

	water_atom_distance = water_to_other_norm;

	Real dot = water_unit.dot( water_to_other_unit );

	debug_assert( dot <= 1.01 && dot >= -1.01);
	dot = std::min<Real>(dot, 1 );
	dot = std::max<Real>(dot, -1 );

	Real dist2;
	if ( dot < min_angle_cos_ ) {
		Real off_axis = water_to_other_norm * std::sqrt( 1 - dot*dot );
		Real sin_dist = std::abs( off_axis - min_angle_sin_ * w_dist_ );
		Real cos_dist = std::abs( water_to_other_norm * dot - min_angle_cos_ * w_dist_ );
		dist2 = sin_dist * sin_dist + cos_dist * cos_dist;
	} else if ( dot > max_angle_cos_ ) {
		Real off_axis = water_to_other_norm * std::sqrt( 1 - dot*dot );
		Real sin_dist = std::abs( off_axis - max_angle_sin_ * w_dist_ );
		Real cos_dist = std::abs( water_to_other_norm * dot - max_angle_cos_ * w_dist_ );
		dist2 = sin_dist * sin_dist + cos_dist * cos_dist;
	} else {
		// atoms are colinear
		Real dist = std::abs( water_to_other_norm - w_dist_ );
		dist2 = dist * dist;
	}

	// if ( dump_dome_waters_ ) {


	//     Vector out_of_plane = water_unit.cross( water_to_other_unit ).normalized();
	//     Vector sin_axis = out_of_plane.cross(water_unit).normalized();

	//     Vector water_pos;

	//     if ( dot < min_angle_cos_ ) {
	//         // std::cout << "ldot" << std::endl;

	//         water_pos = water + sin_axis * min_angle_sin_ * w_dist_ + water_unit * min_angle_cos_ * w_dist_;

	//     } else if ( dot > max_angle_cos_ ) {
	//         // std::cout << "hdot" << std::endl;

	//         water_pos = water + sin_axis * max_angle_sin_ * w_dist_ + water_unit * max_angle_cos_ * w_dist_;
	//     } else {
	//         // std::cout << "cdot" << std::endl;
	//         // atoms are colinear
	//         water_pos = water + water_to_other_unit * w_dist_;
	//     }

	//     // std::cout << "TEST: " << water_pos.distance( other_atom_xyz ) << " " << dist << std::endl;
	//     runtime_assert( std::abs( water_pos.distance( other_atom_xyz ) - dist ) < 0.01 );

	//     debug_dome_waters_.push_back( water_pos );


	// }

	return dist2;

}

void
LK_DomeEnergy::single_water_atom_fractions(
	Vector const & base_atom_xyz,
	Vector const & water_xyz,
	Vector const & other_atom_xyz,
	Size other_atom_type_index,
	Real & lk_dome_iso_frac,
	Real & lk_dome_frac
) const {

	Real ballw_dist = 0;
	Vector _2;
	Real dw_dist2 = dome_water_dist2(base_atom_xyz, water_xyz, other_atom_xyz, ballw_dist, _2 );

	lk_dome_iso_frac = my_lk_fraction( ballw_dist );
	lk_dome_frac = lk_dome_iso_frac * get_lkd_frac_dome_dist( dw_dist2, other_atom_type_index );
}


void
LK_DomeEnergy::accumulate_single_atom_contributions(
	Size const atom1,
	Size const,
	Size atom1_n_attached_waters,
	Size atom1_start_water,
	WaterCoords const & atom1_waters,
	WaterOcclusions const & atom1_occlusions,
	conformation::Residue const & rsd1,
	Size const atom2_type_index,
	Vector const & atom2_xyz,
	Real const lk_desolvation_of_atom1_by_water,
	EnergyMap & emap
) const
{
	if ( atom1_n_attached_waters == 0 ) return;

	for ( Size i_water1 = 1; i_water1 <= atom1_n_attached_waters; i_water1++ ) {

		if ( atom1_occlusions[i_water1] > occlusion_max_ ) continue;

		if ( atom1_waters[i_water1 + atom1_start_water].distance_squared( atom2_xyz) > water_atom_interaction_cutoff2() ) continue;

		Real availability = get_avail( atom1_occlusions[i_water1] );

		Real lk_dome_iso_frac = 0;
		Real lk_dome_frac = 0;

		single_water_atom_fractions( rsd1.xyz(atom1), atom1_waters[i_water1 + atom1_start_water], atom2_xyz, atom2_type_index, lk_dome_iso_frac, lk_dome_frac );

		// std::cout << i_water1 << " " << lk_desolvation_of_atom1_by_water << std::endl;

		// fraction of "fa_sol" of water
		emap[ lk_dome_iso ] += -1 * lk_dome_iso_frac * lk_desolvation_of_atom1_by_water * availability;
		emap[ lk_dome ]     += -1 * lk_dome_frac     * lk_desolvation_of_atom1_by_water * availability;
	}

}

Real
cubed_root( Real val ) {
	return std::pow( val, (3.0/2.0) );
}




DerivativeFinderWadj::DerivativeFinderWadj( Vector const & base, Vector const & water, Real w_dist, Real water_adjust )
: B( base ),
	W( water ),
	B_to_W( water - base ),
	wadj( water_adjust ),
	w_dist_( w_dist ),
	x(14)
{
	B_to_W_norm = B_to_W.norm();
	B_to_W_norm2 = B_to_W_norm*B_to_W_norm;
	B_to_W_norm3 = B_to_W_norm2*B_to_W_norm;

	x[1] = wadj/B_to_W_norm;
	x[2] = -x[1];
	x[3] = wadj/B_to_W_norm3;
	x[4] = numeric::square(B_to_W[0])*x[3];
	x[5] = B_to_W[0]*x[3];
	x[6] = B_to_W[1]*x[5];
	x[7] = B_to_W[2]*x[5];
	x[8] = numeric::square(B_to_W[1])*x[3];
	x[9] = B_to_W[1]*B_to_W[2]*x[3];
	x[10] = numeric::square(B_to_W[2])*x[3];
	x[11] = x[1] + 1;
	x[12] = -x[6];
	x[13] = -x[7];
	x[14] = -x[9];
}

numeric::xyzMatrix<Real>
DerivativeFinderWadj::dWadj_dBase() {
	numeric::xyzMatrix<Real> out;
	out(1, 1) = x[2] + x[4];
	out(2, 1) = x[6];
	out(3, 1) = x[7];
	out(1, 2) = x[6];
	out(2, 2) = x[2] + x[8];
	out(3, 2) = x[9];
	out(1, 3) = x[7];
	out(2, 3) = x[9];
	out(3, 3) = x[10] + x[2];
	return out;
}

numeric::xyzMatrix<Real>
DerivativeFinderWadj::dWadj_dWater() {
	numeric::xyzMatrix<Real> out;
	out(1, 1) = x[11] - x[4];
	out(2, 1) = x[12];
	out(3, 1) = x[13];
	out(1, 2) = x[12];
	out(2, 2) = x[11] - x[8];
	out(3, 2) = x[14];
	out(1, 3) = x[13];
	out(2, 3) = x[14];
	out(3, 3) = -x[10] + x[11];
	return out;
}




DerivativeFinder::DerivativeFinder( Vector const & base, Vector const & water, Vector const & other,
	Real min_cose, Real min_sine, Real max_cose, Real max_sine, Real w_dist, Real water_adjust ) :
	B(base),
	W(water),
	Ot(other),
	W_to_Ot(Ot - W),
	B_to_W(W - B),
	w_dist_( w_dist ),
	wadj( water_adjust ),
	x(175)
{
	W_to_Ot_norm = W_to_Ot.norm();
	W_to_Ot_norm2 = W_to_Ot_norm*W_to_Ot_norm;
	W_to_Ot_norm3 = W_to_Ot_norm*W_to_Ot_norm2;

	B_to_W_norm = B_to_W.norm();
	B_to_W_norm2 = B_to_W_norm*B_to_W_norm;
	B_to_W_norm3 = B_to_W_norm2*B_to_W_norm;

	norm_prod = W_to_Ot_norm * B_to_W_norm;

	Vector wadj_water = W + B_to_W / B_to_W_norm * wadj;
	Vector wadj_to_Ot_unit = (Ot - wadj_water).normalized();

	Real dot = wadj_to_Ot_unit.dot( B_to_W ) / B_to_W_norm;

	if ( dot < min_cose ) {
		colinear_water = false;
		cose = min_cose;
		sine = min_sine;
	} else if ( dot > max_cose ) {
		colinear_water = false;
		cose = max_cose;
		sine = max_sine;
	} else {
		colinear_water = true;
		cose = 0;
		sine = 0;
	}

	// Hopefully it's clear how this works from the following code -- bcov
	if ( colinear_water ) {
		x[1] = wadj/B_to_W_norm;
		x[2] = -x[1];
		x[3] = wadj/B_to_W_norm3;
		x[4] = numeric::square(B_to_W[0])*x[3];
		x[5] = w_dist_*x[1];
		x[6] = w_dist_*x[4];
		x[7] = B_to_W[0]*x[1];
		x[8] = W_to_Ot[0] - x[7];
		x[9] = B_to_W[1]*x[1];
		x[10] = W_to_Ot[1] - x[9];
		x[11] = B_to_W[2]*x[1];
		x[12] = W_to_Ot[2] - x[11];
		x[13] = numeric::square(x[10]) + numeric::square(x[12]) + numeric::square(x[8]);
		x[14] = 1/std::sqrt(x[13]);
		x[15] = B_to_W[0]*x[3];
		x[16] = B_to_W[1]*x[15];
		x[17] = x[10]*x[16];
		x[18] = B_to_W[2]*x[15];
		x[19] = x[12]*x[18];
		x[20] = 2*x[1];
		x[21] = 2*x[4];
		x[22] = x[8]/2;
		x[23] = w_dist_/cubed_root(x[13]);
		x[24] = x[23]*(x[17] + x[19] - x[22]*(x[20] - x[21]));
		x[25] = w_dist_*x[14];
		x[26] = x[16]*x[25];
		x[27] = x[16] - x[26];
		x[28] = x[18]*x[25];
		x[29] = x[18] - x[28];
		x[30] = x[16]*x[8];
		x[31] = B_to_W[1]*B_to_W[2]*x[3];
		x[32] = x[12]*x[31];
		x[33] = numeric::square(B_to_W[1])*x[3];
		x[34] = 2*x[33];
		x[35] = x[10]/2;
		x[36] = x[23]*(x[30] + x[32] - x[35]*(x[20] - x[34]));
		x[37] = w_dist_*x[33];
		x[38] = x[25]*x[31];
		x[39] = x[31] - x[38];
		x[40] = x[18]*x[8];
		x[41] = x[10]*x[31];
		x[42] = numeric::square(B_to_W[2])*x[3];
		x[43] = 2*x[42];
		x[44] = x[12]/2;
		x[45] = x[23]*(x[40] + x[41] - x[44]*(x[20] - x[43]));
		x[46] = w_dist_*x[42];
		x[47] = -x[5] - w_dist_;
		x[48] = -x[20] - 2;
		x[49] = x[23]*(-x[17] - x[19] - x[22]*(x[21] + x[48]));
		x[50] = x[1] + 1;
		x[51] = -x[16] + x[26];
		x[52] = -x[18] + x[28];
		x[53] = x[23]*(-x[30] - x[32] - x[35]*(x[34] + x[48]));
		x[54] = -x[31] + x[38];
		x[55] = x[23]*(-x[40] - x[41] - x[44]*(x[43] + x[48]));
		x[56] = x[23]*(-W_to_Ot[0] + x[7]);
		x[57] = x[23]*(-W_to_Ot[1] + x[9]);
		x[58] = x[23]*(-W_to_Ot[2] + x[11]);
	} else {
		x[1] = numeric::square(B_to_W[0]);
		x[2] = 1/B_to_W_norm3;
		x[3] = x[1]*x[2];
		x[4] = wadj*x[3];
		x[5] = w_dist_*cose;
		x[6] = x[3]*x[5];
		x[7] = 1/B_to_W_norm;
		x[8] = wadj*x[7];
		x[9] = B_to_W[0]*x[8];
		x[10] = W_to_Ot[0] - x[9];
		x[11] = x[10]*x[7];
		x[12] = B_to_W[0]*x[11];
		x[13] = B_to_W[1]*x[8];
		x[14] = W_to_Ot[1] - x[13];
		x[15] = x[14]*x[7];
		x[16] = B_to_W[1]*x[15];
		x[17] = B_to_W[2]*x[8];
		x[18] = W_to_Ot[2] - x[17];
		x[19] = x[18]*x[7];
		x[20] = B_to_W[2]*x[19];
		x[21] = x[12] + x[16] + x[20];
		x[22] = x[21]*x[7];
		x[23] = numeric::square(B_to_W[1]);
		x[24] = wadj/numeric::square(B_to_W_norm2);
		x[25] = B_to_W[0]*x[24];
		x[26] = x[23]*x[25];
		x[27] = numeric::square(B_to_W[2]);
		x[28] = x[25]*x[27];
		x[29] = B_to_W[0]*x[2];
		x[30] = B_to_W[1]*x[29];
		x[31] = x[14]*x[30];
		x[32] = B_to_W[2]*x[29];
		x[33] = x[18]*x[32];
		x[34] = x[10]*x[3];
		x[35] = -x[4] + x[8];
		x[36] = B_to_W[0]*x[7];
		x[37] = -x[26] - x[28] + x[31] + x[33] + x[34] + x[35]*x[36] + x[7]*(-W_to_Ot[0] + x[9]);
		x[38] = x[36]*x[37];
		x[39] = x[21]*x[3];
		x[40] = -B_to_W[0]*x[22] + x[10];
		x[41] = -B_to_W[1]*x[22] + x[14];
		x[42] = -B_to_W[2]*x[22] + x[18];
		x[43] = numeric::square(x[40]) + numeric::square(x[41]) + numeric::square(x[42]);
		x[44] = w_dist_*sine;
		x[45] = x[44]/std::sqrt(x[43]);
		x[46] = x[37]*x[7];
		x[47] = B_to_W[1]*x[46];
		x[48] = wadj*x[30];
		x[49] = 2*x[48];
		x[50] = x[21]*x[30];
		x[51] = 2*x[50];
		x[52] = -x[49] - x[51];
		x[53] = x[41]/2;
		x[54] = B_to_W[2]*x[46];
		x[55] = wadj*x[32];
		x[56] = 2*x[55];
		x[57] = x[21]*x[32];
		x[58] = 2*x[57];
		x[59] = -x[56] - x[58];
		x[60] = x[42]/2;
		x[61] = 2*x[4];
		x[62] = 2*x[39];
		x[63] = 2*x[12];
		x[64] = 2*x[16];
		x[65] = 2*x[20];
		x[66] = 2*x[8];
		x[67] = x[66] + x[7]*(x[63] + x[64] + x[65]);
		x[68] = x[40]/2;
		x[69] = x[44]/cubed_root(x[43]);
		x[70] = x[69]*(-x[53]*(-2*x[47] + x[52]) - x[60]*(-2*x[54] + x[59]) - x[68]*(-2*x[38] - x[61] - x[62] + x[67]));
		x[71] = -x[8];
		x[72] = x[5]*x[7];
		x[73] = x[71] - x[72];
		x[74] = -x[48];
		x[75] = -x[50] + x[74];
		x[76] = x[30]*x[5];
		x[77] = x[48] + x[76];
		x[78] = -x[55];
		x[79] = -x[57] + x[78];
		x[80] = x[32]*x[5];
		x[81] = x[55] + x[80];
		x[82] = B_to_W[1]*x[24];
		x[83] = x[1]*x[82];
		x[84] = x[27]*x[82];
		x[85] = x[10]*x[30];
		x[86] = B_to_W[1]*B_to_W[2];
		x[87] = x[2]*x[86];
		x[88] = x[18]*x[87];
		x[89] = x[2]*x[23];
		x[90] = x[14]*x[89];
		x[91] = wadj*x[89];
		x[92] = x[8] - x[91];
		x[93] = B_to_W[1]*x[7];
		x[94] = x[7]*(-W_to_Ot[1] + x[13]) - x[83] - x[84] + x[85] + x[88] + x[90] + x[92]*x[93];
		x[95] = x[36]*x[94];
		x[96] = B_to_W[2]*x[7];
		x[97] = x[94]*x[96];
		x[98] = wadj*x[87];
		x[99] = 2*x[98];
		x[100] = x[21]*x[87];
		x[101] = 2*x[100];
		x[102] = -x[101] - x[99];
		x[103] = x[93]*x[94];
		x[104] = 2*x[91];
		x[105] = x[21]*x[89];
		x[106] = 2*x[105];
		x[107] = x[69]*(-x[53]*(-2*x[103] - x[104] - x[106] + x[67]) - x[60]*(x[102] - 2*x[97]) - x[68]*(x[52] - 2*x[95]));
		x[108] = x[5]*x[89];
		x[109] = -x[98];
		x[110] = -x[100] + x[109];
		x[111] = x[5]*x[87];
		x[112] = x[111] + x[98];
		x[113] = B_to_W[2]*x[24];
		x[114] = x[1]*x[113];
		x[115] = x[113]*x[23];
		x[116] = x[10]*x[32];
		x[117] = x[14]*x[87];
		x[118] = x[2]*x[27];
		x[119] = x[118]*x[18];
		x[120] = wadj*x[118];
		x[121] = -x[120] + x[8];
		x[122] = -x[114] - x[115] + x[116] + x[117] + x[119] + x[121]*x[96] + x[7]*(-W_to_Ot[2] + x[17]);
		x[123] = x[122]*x[36];
		x[124] = x[122]*x[93];
		x[125] = x[122]*x[96];
		x[126] = 2*x[120];
		x[127] = x[118]*x[21];
		x[128] = 2*x[127];
		x[129] = x[69]*(-x[53]*(x[102] - 2*x[124]) - x[60]*(-2*x[125] - x[126] - x[128] + x[67]) - x[68]*(-2*x[123] + x[59]));
		x[130] = x[118]*x[5];
		x[131] = x[7]*(-x[12] - x[16] - x[20]);
		x[132] = x[71] - 1;
		x[133] = x[132] + x[4];
		x[134] = x[11] + x[133]*x[36] + x[26] + x[28] - x[31] - x[33] - x[34];
		x[135] = x[134]*x[36];
		x[136] = x[134]*x[93];
		x[137] = x[49] + x[51];
		x[138] = x[134]*x[96];
		x[139] = x[56] + x[58];
		x[140] = -x[66] + x[7]*(-x[63] - x[64] - x[65]) - 2;
		x[141] = x[69]*(-x[53]*(-2*x[136] + x[137]) - x[60]*(-2*x[138] + x[139]) - x[68]*(-2*x[135] + x[140] + x[61] + x[62]));
		x[142] = x[72] + 1;
		x[143] = x[48] + x[50];
		x[144] = x[74] - x[76];
		x[145] = x[55] + x[57];
		x[146] = x[78] - x[80];
		x[147] = x[132] + x[91];
		x[148] = x[147]*x[93] + x[15] + x[83] + x[84] - x[85] - x[88] - x[90];
		x[149] = x[148]*x[36];
		x[150] = x[148]*x[96];
		x[151] = x[101] + x[99];
		x[152] = x[148]*x[93];
		x[153] = x[69]*(-x[53]*(x[104] + x[106] + x[140] - 2*x[152]) - x[60]*(-2*x[150] + x[151]) - x[68]*(x[137] - 2*x[149]));
		x[154] = x[100] + x[98];
		x[155] = x[109] - x[111];
		x[156] = x[120] + x[132];
		x[157] = x[114] + x[115] - x[116] - x[117] - x[119] + x[156]*x[96] + x[19];
		x[158] = x[157]*x[36];
		x[159] = x[157]*x[93];
		x[160] = x[157]*x[96];
		x[161] = x[69]*(-x[53]*(x[151] - 2*x[159]) - x[60]*(x[126] + x[128] + x[140] - 2*x[160]) - x[68]*(x[139] - 2*x[158]));
		x[162] = 1/B_to_W_norm2;
		x[163] = x[1]*x[162];
		x[164] = B_to_W[0]*x[162];
		x[165] = B_to_W[1]*x[164];
		x[166] = B_to_W[2]*x[164];
		x[167] = x[69]*(x[165]*x[41] + x[166]*x[42] - x[68]*(2 - 2*x[163]));
		x[168] = -x[165]*x[45];
		x[169] = -x[166]*x[45];
		x[170] = x[162]*x[86];
		x[171] = x[162]*x[23];
		x[172] = x[69]*(x[165]*x[40] + x[170]*x[42] - x[53]*(2 - 2*x[171]));
		x[173] = -x[170]*x[45];
		x[174] = x[162]*x[27];
		x[175] = x[69]*(x[166]*x[40] + x[170]*x[41] - x[60]*(2 - 2*x[174]));
	}

}
numeric::xyzMatrix<Real>
DerivativeFinder::dDome_dBase() {
	if ( colinear_water ) {
		return colinear_dDome_dBase();
	} else {
		return fringe_dDome_dBase();
	}
}

numeric::xyzMatrix<Real>
DerivativeFinder::dDome_dWater() {
	if ( colinear_water ) {
		return colinear_dDome_dWater();
	} else {
		return fringe_dDome_dWater();
	}
}

numeric::xyzMatrix<Real>
DerivativeFinder::dDome_dOther() {
	if ( colinear_water ) {
		return colinear_dDome_dOther();
	} else {
		return fringe_dDome_dOther();
	}
}


numeric::xyzMatrix<Real>
DerivativeFinder::colinear_dDome_dBase() {
	numeric::xyzMatrix<Real> out;
	out(1, 1) = x[14]*(x[5] - x[6]) + x[2] + x[24]*x[8] + x[4];
	out(2, 1) = x[10]*x[24] + x[27];
	out(3, 1) = x[12]*x[24] + x[29];
	out(1, 2) = x[27] + x[36]*x[8];
	out(2, 2) = x[10]*x[36] + x[14]*(-x[37] + x[5]) + x[2] + x[33];
	out(3, 2) = x[12]*x[36] + x[39];
	out(1, 3) = x[29] + x[45]*x[8];
	out(2, 3) = x[10]*x[45] + x[39];
	out(3, 3) = x[12]*x[45] + x[14]*(-x[46] + x[5]) + x[2] + x[42];
	return out;
}

numeric::xyzMatrix<Real>
DerivativeFinder::colinear_dDome_dWater() {
	numeric::xyzMatrix<Real> out;
	out(1, 1) = x[14]*(x[47] + x[6]) - x[4] + x[49]*x[8] + x[50];
	out(2, 1) = x[10]*x[49] + x[51];
	out(3, 1) = x[12]*x[49] + x[52];
	out(1, 2) = x[51] + x[53]*x[8];
	out(2, 2) = x[10]*x[53] + x[14]*(x[37] + x[47]) - x[33] + x[50];
	out(3, 2) = x[12]*x[53] + x[54];
	out(1, 3) = x[52] + x[55]*x[8];
	out(2, 3) = x[10]*x[55] + x[54];
	out(3, 3) = x[12]*x[55] + x[14]*(x[46] + x[47]) - x[42] + x[50];
	return out;
}


numeric::xyzMatrix<Real>
DerivativeFinder::colinear_dDome_dOther() {
	numeric::xyzMatrix<Real> out;
	out(1, 1) = x[25] + x[56]*x[8];
	out(2, 1) = x[10]*x[56];
	out(3, 1) = x[12]*x[56];
	out(1, 2) = x[57]*x[8];
	out(2, 2) = x[10]*x[57] + x[25];
	out(3, 2) = x[12]*x[57];
	out(1, 3) = x[58]*x[8];
	out(2, 3) = x[10]*x[58];
	out(3, 3) = x[12]*x[58] + x[25];
	return out;
}

numeric::xyzMatrix<Real>
DerivativeFinder::fringe_dDome_dBase() {
	numeric::xyzMatrix<Real> out;
	out(1, 1) = x[4] + x[40]*x[70] + x[45]*(x[22] + x[35] - x[38] - x[39]) + x[6] + x[73];
	out(2, 1) = x[41]*x[70] + x[45]*(-x[47] + x[75]) + x[77];
	out(3, 1) = x[42]*x[70] + x[45]*(-x[54] + x[79]) + x[81];
	out(1, 2) = x[107]*x[40] + x[45]*(x[75] - x[95]) + x[77];
	out(2, 2) = x[107]*x[41] + x[108] + x[45]*(-x[103] - x[105] + x[22] + x[92]) + x[73] + x[91];
	out(3, 2) = x[107]*x[42] + x[112] + x[45]*(x[110] - x[97]);
	out(1, 3) = x[129]*x[40] + x[45]*(-x[123] + x[79]) + x[81];
	out(2, 3) = x[112] + x[129]*x[41] + x[45]*(x[110] - x[124]);
	out(3, 3) = x[120] + x[129]*x[42] + x[130] + x[45]*(x[121] - x[125] - x[127] + x[22]) + x[73];
	return out;
}
numeric::xyzMatrix<Real>
DerivativeFinder::fringe_dDome_dWater() {
	numeric::xyzMatrix<Real> out;
	out(1, 1) = x[141]*x[40] + x[142] + x[35] + x[45]*(x[131] + x[133] - x[135] + x[39]) - x[6];
	out(2, 1) = x[141]*x[41] + x[144] + x[45]*(-x[136] + x[143]);
	out(3, 1) = x[141]*x[42] + x[146] + x[45]*(-x[138] + x[145]);
	out(1, 2) = x[144] + x[153]*x[40] + x[45]*(x[143] - x[149]);
	out(2, 2) = -x[108] + x[142] + x[153]*x[41] + x[45]*(x[105] + x[131] + x[147] - x[152]) + x[92];
	out(3, 2) = x[153]*x[42] + x[155] + x[45]*(-x[150] + x[154]);
	out(1, 3) = x[146] + x[161]*x[40] + x[45]*(x[145] - x[158]);
	out(2, 3) = x[155] + x[161]*x[41] + x[45]*(x[154] - x[159]);
	out(3, 3) = x[121] - x[130] + x[142] + x[161]*x[42] + x[45]*(x[127] + x[131] + x[156] - x[160]);
	return out;
}
numeric::xyzMatrix<Real>
DerivativeFinder::fringe_dDome_dOther() {
	numeric::xyzMatrix<Real> out;
	out(1, 1) = x[167]*x[40] + x[45]*(1 - x[163]);
	out(2, 1) = x[167]*x[41] + x[168];
	out(3, 1) = x[167]*x[42] + x[169];
	out(1, 2) = x[168] + x[172]*x[40];
	out(2, 2) = x[172]*x[41] + x[45]*(1 - x[171]);
	out(3, 2) = x[172]*x[42] + x[173];
	out(1, 3) = x[169] + x[175]*x[40];
	out(2, 3) = x[173] + x[175]*x[41];
	out(3, 3) = x[175]*x[42] + x[45]*(1 - x[174]);
	return out;
}





/// @note  atom2 is desolvating atom1. atom1_waters is non-empty
/// @note  Pretend that atom1 is the atom whose derivs are being calculated. weight_factor may include -1 term
/// to switch the order...
///
void
LK_DomeEnergy::sum_deriv_contributions_for_heavyatom_pair_one_way(
	Size const heavyatom1,
	conformation::Residue const & rsd1,
	LKD_ResidueInfo const & rsd1_info,
	Size const heavyatom2,
	conformation::Residue const & rsd2,
	LKD_ResidueInfo const & rsd2_info,
	EnergyMap const & weights,
	Real const weight_factor,
	Real const,
	utility::vector1< DerivVectorPair > & r1_at_derivs,
	utility::vector1< DerivVectorPair > & r2_at_derivs
) const
{
	Size atom1_n_waters = rsd1_info.n_attached_waters()[ heavyatom1 ];
	if ( atom1_n_waters == 0 ) return;

	Size atom1_offset = rsd1_info.water_offset_for_atom()[heavyatom1];
	Size atom2_offset = rsd2_info.water_offset_for_atom()[heavyatom2];
	WaterCoords const & res1_waters = rsd1_info.waters();
	// AtomWeights const & atom1_wts = rsd1_info.atom_weights()[heavyatom1];
	Size atom2_n_waters = rsd2_info.n_attached_waters()[ heavyatom2 ];
	WaterCoords const & res2_waters = rsd2_info.waters();
	// (void)atom2_waters;
	// (void)atom2_n_waters;


	Real lk_dome_wt = weights[core::scoring::lk_dome];
	Real lk_dome_iso_wt = weights[core::scoring::lk_dome_iso];
	Real lk_dome_bridge_wt = weights[core::scoring::lk_dome_bridge];
	Real lk_dome_bridge_uncpl_wt = weights[core::scoring::lk_dome_bridge_uncpl];
	Real lk_ball_bridge2_wt = weights[core::scoring::lk_ball_bridge2];
	Real lk_ball_bridge_uncpl2_wt = weights[core::scoring::lk_ball_bridge_uncpl2];
	// (void)lk_dome_wt;
	// (void)lk_dome_bridge_wt;
	// (void)lk_dome_bridge_uncpl_wt;


	Vector const & heavyatom1_xyz( rsd1.xyz( heavyatom1 ) );
	Vector const & atom2_xyz( rsd2.xyz( heavyatom2 ) );
	// Vector f1( heavyatom1_xyz.cross( atom2_xyz ) ), f2( heavyatom1_xyz - atom2_xyz );
	// Real const inv_dis( 1.0 / std::sqrt( d2 ) );
	Size atom2_type_index = rsd2.atom_type_index( heavyatom2 );

	// (void)inv_dis;

	// lk_dome_iso = - avail * sol_value * iso_frac( B, O )
	// dlk_dome_iso_dO = - avail * sol_value * diso_frac_dBO_dist

	// lk_dome = - avail * sol_value * iso_frac( B, O ) * dome_frac( B, W, O )

	// dlk_dome_dR = - avail * sol_value * ( diso_frac_dR( B, O ) * dome_frac( B, W, O ) + iso_frac( B, O ) * ddome_frac_dR( B, W, O ))
	// diso_frac_dR = dlk_frac_sigmoid
	// ddome_frac_dR = 0 if (O is closer than ramp_width_A2) or (O is too far away)
	//                   else eval_d_lk_fraction_dr_over_r * r


	Vector const & base = heavyatom1_xyz;

	Real sol_value = rsd1_info.water_sol_values()[heavyatom1];

	Real water_atom_range = water_atom_interaction_cutoff();
	Real water_atom_range2 = water_atom_range * water_atom_range;


	// WaterDerivContributions dwwd2_ddi;

	// (2) waters-atom2 interaction
	// note that there's no derivative unless we're in the ramping zone:
	// we only do this if this term hasn't already been captured by one of our dependent hydrogens
	// we assume that for heavyatoms with dependent polar hydrogens, every water belongs to one of the hydrogens
	for ( Size iwat=1; iwat <= atom1_n_waters; ++iwat ) {
		Vector const & atom1_water_xyz( res1_waters[ iwat + atom1_offset ] );

		// update f1 and f2 to reflect water-atom2 as the interaction
		Vector f1w = atom1_water_xyz.cross( atom2_xyz );
		Vector f2w = atom1_water_xyz - atom2_xyz;

		if ( f2w.norm_squared() > water_atom_range2 ) continue;

		Real avail = get_avail( rsd1_info.water_occlusions()[heavyatom1][iwat] );
		if ( avail == 0 ) continue;

		Vector const & water = atom1_water_xyz;
		Vector const & other = atom2_xyz;

		Real water_other_distance = 0;
		Vector wadj_water;
		Vector dome_water = get_dome_water( base, water, other, water_other_distance, wadj_water);
		Real dw_dist2 = dome_water.distance_squared( other );
		Real dome_frac = get_lkd_frac_dome_dist( dw_dist2, atom2_type_index );
		Real iso_frac = my_lk_fraction( water_other_distance );

		Vector f1wadj = wadj_water.cross( atom2_xyz );
		Vector f2wadj = wadj_water - atom2_xyz;


		// what is the derivative of the iso_frac wrt r?
		Real diso_frac_dr =  dmy_lk_fraction( water_other_distance );



		Real d2_low = d2_low_[ atom2_type_index ];
		Real d2_delta = dw_dist2 - d2_low;

		Real ddome_frac_dr_over_r = 0;

		if ( d2_delta < ramp_width_A2_ && d2_delta > 0 ) {
			ddome_frac_dr_over_r = eval_d_lk_fraction_dr_over_r( d2_delta, ramp_width_A2_ );
		}

		// This is nightmarish. Within the partial derivatives of lk_dome, one of the terms depends on the
		// Water-other vector and not on the dome-other vector.
		// That must be summed into iso


		// nope nope nope
		// Real dlk_dome_E_dR_over_r = -avail * sol_value * ( diso_frac_dr * dome_frac / dw_dist + iso_frac * ddome_frac_dr_over_r )
		//                                     * lk_dome_wt;

		Real dlk_dome_E_dR_over_r_iso_portion = - weight_factor * avail * sol_value * diso_frac_dr * dome_frac / water_other_distance * lk_dome_wt;
		Real dlk_dome_E_dR_over_r_dome_portion = - weight_factor * avail * sol_value * iso_frac * ddome_frac_dr_over_r * lk_dome_wt;



		//////////////////////////////////////// lk_dome_iso ///////////////////////////////////////////////////

		// R = sqrt(dotsq( W - O ))
		// dR = sqrt`(dotsq( W - O ))*dotsq`(W - O)*d(W - O)
		// dR = 1/2/sqrt(dotsq( W - O ))*2*(W - O)*-dO
		// dR/dO = (O - W)/norm( W - O )
		// dR/dO_times_R = (O - W)

		// f2 = (O - W) * dE_dR_over_r




		Real dE_iso_dR_over_r = - weight_factor * avail * sol_value * lk_dome_iso_wt * diso_frac_dr / water_other_distance;

		// Derivatives for the desolvating atom
		r2_at_derivs[heavyatom2].f1() -= (dE_iso_dR_over_r + dlk_dome_E_dR_over_r_iso_portion) * f1wadj;
		r2_at_derivs[heavyatom2].f2() -= (dE_iso_dR_over_r + dlk_dome_E_dR_over_r_iso_portion) * f2wadj;

		// derivatives for the desolvated atoms
		// dR/datom1 = dR/dwater * dwater/datom1
		//core::Real dwater = f2.length();
		//Real denom = dwater*dwater*dwater;


		// R = sqrt(dotsq( W(L) - O ))
		// dR = sqrt`(dotsq( W(L) - O ))*dotsq`(W(L) - O)*d(W(L) - O)
		// dR = 1/2 /sqrt(dotsq( W(L) - O ))*2*(W(L) - O)*dwater_dl*dL
		// dR/dO = weird_dot( W - O, dwater_dL ) / norm( W - O )
		// dR/dO_times_R = weird_dot( W - O, dwater_dL )

		// f1 = ( weird_dot( W - O, dwater_dL ), L ) * dE_dR_over_r
		// f2 = weird_dot( W - O, dwater_dL ) * dE_dR_over_r



		// frank is storing
		// dE/dR / R * atom1 x ( atom1 - dRdatom1_times_R )
		// dE/dR / R * atom1 x ( atom1 - dRdatom1_times_R )
		// dE/dR / R * - atom1 x dRdatom1_times_R
		// dE/dR * - atom1 x dRdatom1

		///////////////////////////////////////////////////////
		// f1 = (dR/datom1 x atom1 ) * dE/dR
		// f2 = dR/datom1 * dE/dR
		///////////////////////////////////////////////////////

		// f1 = ((atom1 - atom2)/R x atom1) * dE/dR
		// f1 = (-atom2 x atom1 ) * dE/dR/R
		// f1 = ( atom1 x atom2 ) * dE/dR/R

		DerivativeFinderWadj wadj_finder( base, water, w_dist_, water_adjust_ );
		numeric::xyzMatrix<Real> dwadj_dbase = wadj_finder.dWadj_dBase();
		numeric::xyzMatrix<Real> dwadj_dwater = wadj_finder.dWadj_dWater();


		WaterBuilders const & rsd1_wb( rsd1_info.get_water_builder( rsd1 , heavyatom1 ) );
		{
			Size atom1 = rsd1_wb[iwat].atom1();
			Vector const & r1_atom1_xyz( rsd1.xyz( atom1 ) );

			numeric::xyzMatrix< Real >const & dwater_datom1 = rsd1_info.atom1_derivs()[iwat + atom1_offset];

			numeric::xyzMatrix< Real > dwadj_datom1 = dwadj_dwater * dwater_datom1;
			if ( atom1 == heavyatom1 ) dwadj_datom1 += dwadj_dbase;

			Vector dwadj_datom1x ( dwadj_datom1(1,1), dwadj_datom1(2,1), dwadj_datom1(3,1) );
			Vector dwadj_datom1y ( dwadj_datom1(1,2), dwadj_datom1(2,2), dwadj_datom1(3,2) );
			Vector dwadj_datom1z ( dwadj_datom1(1,3), dwadj_datom1(2,3), dwadj_datom1(3,3) );
			Vector dRdatom_times_R( f2wadj.dot( dwadj_datom1x ), f2wadj.dot( dwadj_datom1y ), f2wadj.dot( dwadj_datom1z ) );

			Vector f2t = dRdatom_times_R;
			Vector f1t = r1_atom1_xyz.cross( r1_atom1_xyz - dRdatom_times_R );

			r1_at_derivs[atom1].f1() += (dE_iso_dR_over_r + dlk_dome_E_dR_over_r_iso_portion) * f1t;
			r1_at_derivs[atom1].f2() += (dE_iso_dR_over_r + dlk_dome_E_dR_over_r_iso_portion) * f2t;
		}
		{
			Size atom2 = rsd1_wb[iwat].atom2();
			Vector const & r1_atom2_xyz( rsd1.xyz( atom2 ) );

			numeric::xyzMatrix< Real >const & dwater_datom2 = rsd1_info.atom2_derivs()[iwat + atom1_offset];

			numeric::xyzMatrix< Real > dwadj_datom2 = dwadj_dwater * dwater_datom2;
			if ( atom2 == heavyatom1 ) dwadj_datom2 += dwadj_dbase;

			Vector dwadj_datom2x ( dwadj_datom2(1,1), dwadj_datom2(2,1), dwadj_datom2(3,1) );
			Vector dwadj_datom2y ( dwadj_datom2(1,2), dwadj_datom2(2,2), dwadj_datom2(3,2) );
			Vector dwadj_datom2z ( dwadj_datom2(1,3), dwadj_datom2(2,3), dwadj_datom2(3,3) );
			Vector dRdatom_times_R( f2wadj.dot( dwadj_datom2x ), f2wadj.dot( dwadj_datom2y ), f2wadj.dot( dwadj_datom2z ) );

			Vector f2t = dRdatom_times_R;
			Vector f1t = r1_atom2_xyz.cross( r1_atom2_xyz - dRdatom_times_R );

			r1_at_derivs[atom2].f1() += (dE_iso_dR_over_r + dlk_dome_E_dR_over_r_iso_portion) * f1t;
			r1_at_derivs[atom2].f2() += (dE_iso_dR_over_r + dlk_dome_E_dR_over_r_iso_portion) * f2t;
		}
		{
			Size atom3 = rsd1_wb[iwat].atom3();
			Vector const & r1_atom3_xyz( rsd1.xyz( atom3 ) );

			numeric::xyzMatrix< Real >const & dwater_datom3 = rsd1_info.atom3_derivs()[iwat + atom1_offset];

			numeric::xyzMatrix< Real > dwadj_datom3 = dwadj_dwater * dwater_datom3;
			if ( atom3 == heavyatom1 ) dwadj_datom3 += dwadj_dbase;

			Vector dwadj_datom3x ( dwadj_datom3(1,1), dwadj_datom3(2,1), dwadj_datom3(3,1) );
			Vector dwadj_datom3y ( dwadj_datom3(1,2), dwadj_datom3(2,2), dwadj_datom3(3,2) );
			Vector dwadj_datom3z ( dwadj_datom3(1,3), dwadj_datom3(2,3), dwadj_datom3(3,3) );
			Vector dRdatom_times_R( f2wadj.dot( dwadj_datom3x ), f2wadj.dot( dwadj_datom3y ), f2wadj.dot( dwadj_datom3z ) );

			Vector f2t = dRdatom_times_R;
			Vector f1t = r1_atom3_xyz.cross( r1_atom3_xyz - dRdatom_times_R );

			r1_at_derivs[atom3].f1() += (dE_iso_dR_over_r + dlk_dome_E_dR_over_r_iso_portion) * f1t;
			r1_at_derivs[atom3].f2() += (dE_iso_dR_over_r + dlk_dome_E_dR_over_r_iso_portion) * f2t;
		}

		if ( dlk_dome_E_dR_over_r_dome_portion == 0 ) continue;

		/////////////////////////////////////////// lk_dome ////////////////////////////////////////////////////////////

		// dlk_dome_dR = - avail * sol_value * ( diso_frac_dR( B, O ) * dome_frac( B, W, O ) + iso_frac( B, O ) * ddome_frac_dR( B, W, O ))
		// diso_frac_dR = dlk_frac_sigmoid
		// ddome_frac_dR = 0 if (O is closer than ramp_width_A2) or (O is too far away)
		//                   else eval_d_lk_fraction_dr_over_r


		// ddome_frac_dO

		// check if we are inside the ramping region


		// nope nope nope, see iso
		// Real dlk_dome_E_dR_over_r = -avail * sol_value * ( diso_frac_dr * dome_frac / dw_dist + iso_frac * ddome_frac_dr_over_r )
		//                                     * lk_dome_wt;


		DerivativeFinder finder( base, water, other, min_angle_cos_, min_angle_sin_, max_angle_cos_, max_angle_sin_, w_dist_, water_adjust_ );

		// update f1 and f2 to reflect dome-atom2 as the interaction
		// Vector f1d = dome_water.cross( atom2_xyz );
		Vector f2d = dome_water - atom2_xyz;


		// Math for other

		// f2 is actually: dR/datom * dE/dR
		// not dR/dAtom * dE/dR/R

		// let R_vec = D - atom2
		// dR/datom2 = d(|D - atom2|)/datom2 = d(sqrt(dotsq(D - atom2)))/datom2
		//
		// dR/datom2 = (pR/pD + pR/patom2)
		// pR/pD = p(sqrt(dotsq(D-atom2)))/pD
		// pD/patom2 = ddome_dother
		//

		// dotsq(X) = dot( X, X ) = X[0]**2 + X[1]**2 + X[2]**2
		// ddotsq/dX = 2*X[0]*(1, 0, 0) + 2*X[1]*(0, 1, 0) + 2*X[2]*(0, 0, 1)
		// ddotsq/dX = 2*X

		// R = sqrt(dotsq( D(atom2) - atom2 ))
		// D(atom2) = dome_water( B, W, atom2 )
		// pD/patom2_B_W = ddome_dother( B, W, atom2 )

		// dR = sqrt`( dotsq(D(atom2) - atom2 )) * d( dotsq(D(atom2) - atom2 ) )
		// dR = sqrt`( dotsq(D(atom2) - atom2 )) * dotsq`(D(atom2) - atom2 ) ) * d( D(atom2) - atom2 )
		// dR = sqrt`( dotsq(D(atom2) - atom2 )) * dotsq`(D(atom2) - atom2 ) ) * ( D`(atom2)*datom2 - 1 * datom2 )
		// dR = 1/2 / sqrt( dotsq(D(atom2) - atom2 )) * 2 * (D(atom2) - atom2 ) ) * ( ddome_dother *datom2 - 1 * datom2 )
		// dR/datom2 =  weird_dot( D(atom2) - atom2, ddome_dother - 1 ) / norm( (Datom2 - atom2 ) )

		// dR/datom2_times_r = weird_dot( D(atom2) - atom2, ddome_dother - 1 )

		// f1 = ( weird_dot( D(atom2) - atom2, ddome_dother - 1 ) x atom2 ) * dlk_dome_E_dR_over_r
		// f2 = weird_dot( D(atom2) - atom2, ddome_dother - 1 ) * dlk_dome_E_dR_over_r

		// f1 = ( weird_dot( f2d, ddome_dother - 1 ) x atom2 ) * dlk_dome_E_dR_over_r
		// f2 = weird_dot( fd2, ddome_dother - 1 ) * dlk_dome_E_dR_over_r



		// dR/datom2 =  ( weird_dot( D(atom2) - atom2, ddome_dother ) + weird_dot( D(atom2) - atom2, -1 ) ) / norm( (Datom2 - atom2 ) )
		// dR/datom2 =  ( weird_dot( D(atom2) - atom2, ddome_dother ) - (D(atom2) - atom2) ) / norm( (Datom2 - atom2 ) )


		///////////////////////////////////////////////////////
		// f1 = (dR/datom1 x atom1 ) * dE/dR
		// f2 = dR/datom1 * dE/dR
		///////////////////////////////////////////////////////

		numeric::xyzMatrix<Real> ddome_dother = finder.dDome_dOther();
		// ddome_dother -= numeric::xyzMatrix<Real>::identity();

		Vector ddome_dotherx ( ddome_dother(1,1), ddome_dother(2,1), ddome_dother(3,1) );
		Vector ddome_dothery ( ddome_dother(1,2), ddome_dother(2,2), ddome_dother(3,2) );
		Vector ddome_dotherz ( ddome_dother(1,3), ddome_dother(2,3), ddome_dother(3,3) );
		Vector dRdother_times_r( f2d.dot( ddome_dotherx ), f2d.dot( ddome_dothery ), f2d.dot( ddome_dotherz ) );
		dRdother_times_r -= f2d;

		Vector f1_cross_times_r = dRdother_times_r.cross( other );

		// These are + because we've calculated dR/datom2 instead of the standard dR/datom1
		r2_at_derivs[heavyatom2].f1() += dlk_dome_E_dR_over_r_dome_portion * f1_cross_times_r;
		r2_at_derivs[heavyatom2].f2() += dlk_dome_E_dR_over_r_dome_portion * dRdother_times_r;

		// Now for Base

		// R = sqrt(dotsq( D(B, W(B)) - atom2 ))
		// D(B, W(B)) = dome_water( B, W(B), atom2 )
		// dD = ddome_dBase( B, W(B), atom2 )*dB + ddome_dwater( B, W(B), atom2)*dWdB * dB
		// dD = (ddome_dbase( B, W, atom2) + ddome_dwater( B, W, atom2 ) * dwater_dbase) dB

		// dR = sqrt`(dotsq( D(B, W(B)) - atom2 )) * dotsq`( D(B, W(B)) - atom2 ) * d(D(B, W(B)) - atom2)
		// dR = 1/2 / sqrt(dotsq( D(B, W(B)) - atom2 )) * 2( D(B, W(B)) - atom2 ) * (ddome_dbase( B, W, atom2) + ddome_dwater( B, W, atom2 ) * dwater_dbase) dB
		// dR =  ( D - atom2 ) / norm( D - atom2 ) * (ddome_dbase( B, W, atom2) + ddome_dwater( B, W, atom2 ) * dwater_dbase) dB
		// dR/dB =  1 / norm( D - atom2 ) * weird_dot(  D - atom2, ddome_dbase( B, W, atom2) + ddome_dwater( B, W, atom2 ) * dwater_dbase)

		// dR/dB =  1 / dw_dist * weird_dot( f2d, ddome_dbase + ddome_dwater * dwater_dbase)

		// dR/dB_times_R = weird_dot( f2d, ddome_dbase + ddome_dwater * dwater_dbase)

		// f1 = ( weird_dot x B ) * dlk_dome_E_dR_over_r
		// f2 = weird_dot * dlk_dome_E_dR_over_r


		// Now for non-base L

		// R = sqrt(dotsq( D(W(L)) - atom2 ))
		// D(W(L)) = dome_water( B, W(L), atom2 )
		// dD = ddome_dwater( B, W(L), atom2)*dWdL * dL
		// dD = ddome_dwater( B, W, atom2 ) * dwater_dL

		// dR = sqrt`(dotsq( D(W(L)) - atom2 )) * dotsq`( D(W(L)) - atom2 ) * d(D(W(L)) - atom2)
		// dR = 1/2 / sqrt(dotsq( D(W(L)) - atom2 )) * 2 *( D(W(L)) - atom2 ) * ddome_dwater( B, W, atom2 ) * dwater_dL * dL
		// dR/dL =  1 / norm( D - atom2 ) * weird_dot( D - atom2, ddome_dwater( B, W, atom2 ) * dwater_dL )

		// dR/dL =  1 / dw_dist * weird_dot( f2d, ddome_dwater * dwater_dbase)

		// dR/dL_times_R = weird_dot( f2d, ddome_dwater * dwater_dbase)

		// f1 = ( weird_dot x L ) * dlk_dome_E_dR_over_r
		// f2 = weird_dot * dlk_dome_E_dR_over_r

		// Only difference is that for base, you add ddome_dbase to the other matrix

		numeric::xyzMatrix<Real> ddome_dwater = finder.dDome_dWater();
		numeric::xyzMatrix<Real> ddome_dbase = finder.dDome_dBase();


		{
			Size atom1 = rsd1_wb[iwat].atom1();
			Vector const & r1_atom1_xyz( rsd1.xyz( atom1 ) );

			numeric::xyzMatrix< Real >const & dwater_datom1 = rsd1_info.atom1_derivs()[iwat+atom1_offset];
			numeric::xyzMatrix< Real > ddome_datom1 = ddome_dwater * dwater_datom1;

			if ( atom1 == heavyatom1 ) ddome_datom1 += ddome_dbase;

			Vector ddome_datom1x ( ddome_datom1(1,1), ddome_datom1(2,1), ddome_datom1(3,1) );
			Vector ddome_datom1y ( ddome_datom1(1,2), ddome_datom1(2,2), ddome_datom1(3,2) );
			Vector ddome_datom1z ( ddome_datom1(1,3), ddome_datom1(2,3), ddome_datom1(3,3) );
			Vector dRdatom_times_R( f2d.dot( ddome_datom1x ), f2d.dot( ddome_datom1y ), f2d.dot( ddome_datom1z ) );

			Vector f2t = dRdatom_times_R;
			Vector f1t = r1_atom1_xyz.cross( r1_atom1_xyz - dRdatom_times_R ); // this is just dRdatom_times_R x r1_atom1_xyz

			r1_at_derivs[atom1].f1() += dlk_dome_E_dR_over_r_dome_portion * f1t;
			r1_at_derivs[atom1].f2() += dlk_dome_E_dR_over_r_dome_portion * f2t;
		}
		{
			Size atom2 = rsd1_wb[iwat].atom2();
			Vector const & r1_atom2_xyz( rsd1.xyz( atom2 ) );

			numeric::xyzMatrix< Real >const & dwater_datom2 = rsd1_info.atom2_derivs()[iwat + atom1_offset];
			numeric::xyzMatrix< Real > ddome_datom2 = ddome_dwater * dwater_datom2;

			if ( atom2 == heavyatom1 ) ddome_datom2 += ddome_dbase;

			Vector ddome_datom2x ( ddome_datom2(1,1), ddome_datom2(2,1), ddome_datom2(3,1) );
			Vector ddome_datom2y ( ddome_datom2(1,2), ddome_datom2(2,2), ddome_datom2(3,2) );
			Vector ddome_datom2z ( ddome_datom2(1,3), ddome_datom2(2,3), ddome_datom2(3,3) );
			Vector dRdatom_times_R( f2d.dot( ddome_datom2x ), f2d.dot( ddome_datom2y ), f2d.dot( ddome_datom2z ) );

			Vector f2t = dRdatom_times_R;
			Vector f1t = r1_atom2_xyz.cross( r1_atom2_xyz - dRdatom_times_R );

			r1_at_derivs[atom2].f1() += dlk_dome_E_dR_over_r_dome_portion * f1t;
			r1_at_derivs[atom2].f2() += dlk_dome_E_dR_over_r_dome_portion * f2t;
		}
		{
			Size atom3 = rsd1_wb[iwat].atom3();
			Vector const & r1_atom3_xyz( rsd1.xyz( atom3 ) );

			numeric::xyzMatrix< Real >const & dwater_datom3 = rsd1_info.atom3_derivs()[iwat + atom1_offset];
			numeric::xyzMatrix< Real > ddome_datom3 = ddome_dwater * dwater_datom3;

			if ( atom3 == heavyatom1 ) ddome_datom3 += ddome_dbase;

			Vector ddome_datom3x ( ddome_datom3(1,1), ddome_datom3(2,1), ddome_datom3(3,1) );
			Vector ddome_datom3y ( ddome_datom3(1,2), ddome_datom3(2,2), ddome_datom3(3,2) );
			Vector ddome_datom3z ( ddome_datom3(1,3), ddome_datom3(2,3), ddome_datom3(3,3) );
			Vector dRdatom_times_R( f2d.dot( ddome_datom3x ), f2d.dot( ddome_datom3y ), f2d.dot( ddome_datom3z ) );

			Vector f2t = dRdatom_times_R;
			Vector f1t = r1_atom3_xyz.cross( r1_atom3_xyz - dRdatom_times_R );

			r1_at_derivs[atom3].f1() += dlk_dome_E_dR_over_r_dome_portion * f1t;
			r1_at_derivs[atom3].f2() += dlk_dome_E_dR_over_r_dome_portion * f2t;
		}

	}






	// // (3) waters-waters interaction
	// // only if lk_ball_bridge is turned on
	// // unlike the other terms, this is _symmetric_ -- instead of A desolv B and B desolv A, we just get A/B water overlap computed once
	// //   thus, we only compute the derivatives for the A component here
	if ( atom2_n_waters == 0 ) return;
	if ( lk_dome_bridge_wt!=0 || lk_dome_bridge_uncpl_wt!=0 ) {


		for ( Size iwat=1; iwat <= atom1_n_waters; ++iwat ) {
			Vector const & atom1_water_xyz( res1_waters[ iwat + atom1_offset] );

			if ( atom1_water_xyz.distance_squared( atom2_xyz) > water_atom_range2 ) continue;

			Real avail = get_avail( rsd1_info.water_occlusions()[heavyatom1][iwat] );
			if ( avail == 0 ) continue;

			// Vector const & base = atom1_xyz;
			Vector const & water = atom1_water_xyz;

			for ( Size jwat=1; jwat <= atom2_n_waters; ++jwat ) {

				Vector const & atom2_water_xyz( res2_waters[ jwat + atom2_offset ] );

				Vector const & other = atom2_water_xyz;



				// Real _;
				// Real dist2 = dome_water_dist2( base_atom_xyz, water_xyz, other_water_xyz, _ );
				// Real d2_delta = dist2;

				// Real this_frac( 0.0 );
				// if ( d2_delta < overlap_width_A2_ ) {
				//     this_frac = ( d2_delta < 0.0 ? Real( 1.0 ) : eval_lk_fraction( d2_delta, overlap_width_A2_ ) );
				// }

				// frac = availability * single_water_water_fraction_1way( atom1_xyz, atom1_waters[i_water], atom2_waters[j_water] );

				// emap[ lk_dome_bridge ] += bridge_frac1 * lk_desolvation_of_atom1_by_water + bridge_frac2 * lk_desolvation_of_atom2_by_water;
				// emap[ lk_dome_bridge_uncpl ] += cp_weight * ( bridge_frac1 + bridge_frac2 );


				// bridge_frac = avail * eval_lk_fraction()

				// dbridge_frac_dr = avail * deval_lk_fraction_dr()

				// bool _;
				Real water_other_distance = 0;
				Vector _2;
				Vector dome_water = get_dome_water( base, water, other, water_other_distance, _2);
				Real dw_dist2 = dome_water.distance_squared( other );


				Real d2_delta = dw_dist2;

				Real dlk_frac_dr_over_r = 0;

				if ( d2_delta < overlap_width_A2_ && d2_delta > 0 ) {
					dlk_frac_dr_over_r = eval_d_lk_fraction_dr_over_r( d2_delta, overlap_width_A2_ );
				}

				Real dbridge_frac_dr_over_r = avail * dlk_frac_dr_over_r;

				if ( dbridge_frac_dr_over_r == 0 ) continue;


				Real dbridgeE_dr_over_r = weight_factor * dbridge_frac_dr_over_r * ( lk_dome_bridge_uncpl_wt + lk_dome_bridge_wt * sol_value);


				Vector f2d = dome_water - other;

				DerivativeFinder finder( base, water, other, min_angle_cos_, min_angle_sin_, max_angle_cos_, max_angle_sin_, w_dist_, water_adjust_ );

				// Let Base2 = L
				// Let O = Water2

				// R = sqrt(normsq( D(O(L)) - O(L)))
				// dR = sqrt`(normsq( D(O(L)) - O(L))) * normsq`( D(O(L)) - O(L)) ) * d( D(O(L)) - O(L)) )
				// dR = 1/2 / norm * 2( f2d ) * ( D`(O(L))*d( O(L) ) - d( O(L)) ) )
				// dR = 1/2 / norm * 2( f2d ) * ( ddome_dother * dother_dL * dL - dother_dL * dL ) )
				// dR/dL = 1/2 / norm * 2( f2d ) * ( ddome_dother * dother_dL - dother_dL ) )
				// dR/dL_times_R = weird_dot( f2d, ddome_dother * dother_dL - dother_dL )



				numeric::xyzMatrix<Real> ddome_dother = finder.dDome_dOther();

				WaterBuilders const & rsd2_wb( rsd2_info.get_water_builder( rsd2 , heavyatom2 ) );
				{
					Size atom1 = rsd2_wb[jwat].atom1();
					Vector const & r2_atom1_xyz( rsd2.xyz( atom1 ) );

					numeric::xyzMatrix< Real >const & dwater_datom1 = rsd2_info.atom1_derivs()[jwat + atom2_offset];
					numeric::xyzMatrix< Real > ddome_datom1 = ddome_dother * dwater_datom1;
					ddome_datom1 -= dwater_datom1;

					Vector ddome_datom1x ( ddome_datom1(1,1), ddome_datom1(2,1), ddome_datom1(3,1) );
					Vector ddome_datom1y ( ddome_datom1(1,2), ddome_datom1(2,2), ddome_datom1(3,2) );
					Vector ddome_datom1z ( ddome_datom1(1,3), ddome_datom1(2,3), ddome_datom1(3,3) );
					Vector dRdatom_times_R( f2d.dot( ddome_datom1x ), f2d.dot( ddome_datom1y ), f2d.dot( ddome_datom1z ) );

					Vector f2t = dRdatom_times_R;
					Vector f1t = r2_atom1_xyz.cross( r2_atom1_xyz - dRdatom_times_R ); // this is just dRdatom_times_R x r1_atom1_xyz

					r2_at_derivs[atom1].f1() += dbridgeE_dr_over_r * f1t;
					r2_at_derivs[atom1].f2() += dbridgeE_dr_over_r * f2t;
				}
				{
					Size atom2 = rsd2_wb[jwat].atom2();
					Vector const & r2_atom2_xyz( rsd2.xyz( atom2 ) );

					numeric::xyzMatrix< Real >const & dwater_datom2 = rsd2_info.atom2_derivs()[jwat + atom2_offset];
					numeric::xyzMatrix< Real > ddome_datom2 = ddome_dother * dwater_datom2;
					ddome_datom2 -= dwater_datom2;

					Vector ddome_datom2x ( ddome_datom2(1,1), ddome_datom2(2,1), ddome_datom2(3,1) );
					Vector ddome_datom2y ( ddome_datom2(1,2), ddome_datom2(2,2), ddome_datom2(3,2) );
					Vector ddome_datom2z ( ddome_datom2(1,3), ddome_datom2(2,3), ddome_datom2(3,3) );
					Vector dRdatom_times_R( f2d.dot( ddome_datom2x ), f2d.dot( ddome_datom2y ), f2d.dot( ddome_datom2z ) );

					Vector f2t = dRdatom_times_R;
					Vector f1t = r2_atom2_xyz.cross( r2_atom2_xyz - dRdatom_times_R );

					r2_at_derivs[atom2].f1() += dbridgeE_dr_over_r * f1t;
					r2_at_derivs[atom2].f2() += dbridgeE_dr_over_r * f2t;
				}
				{
					Size atom3 = rsd2_wb[jwat].atom3();
					Vector const & r2_atom3_xyz( rsd2.xyz( atom3 ) );

					numeric::xyzMatrix< Real >const & dwater_datom3 = rsd2_info.atom3_derivs()[jwat + atom2_offset];
					numeric::xyzMatrix< Real > ddome_datom3 = ddome_dother * dwater_datom3;
					ddome_datom3 -= dwater_datom3;

					Vector ddome_datom3x ( ddome_datom3(1,1), ddome_datom3(2,1), ddome_datom3(3,1) );
					Vector ddome_datom3y ( ddome_datom3(1,2), ddome_datom3(2,2), ddome_datom3(3,2) );
					Vector ddome_datom3z ( ddome_datom3(1,3), ddome_datom3(2,3), ddome_datom3(3,3) );
					Vector dRdatom_times_R( f2d.dot( ddome_datom3x ), f2d.dot( ddome_datom3y ), f2d.dot( ddome_datom3z ) );

					Vector f2t = dRdatom_times_R;
					Vector f1t = r2_atom3_xyz.cross( r2_atom3_xyz - dRdatom_times_R );

					r2_at_derivs[atom3].f1() += dbridgeE_dr_over_r * f1t;
					r2_at_derivs[atom3].f2() += dbridgeE_dr_over_r * f2t;
				}








				numeric::xyzMatrix<Real> ddome_dwater = finder.dDome_dWater();
				numeric::xyzMatrix<Real> ddome_dbase = finder.dDome_dBase();



				WaterBuilders const & rsd1_wb( rsd1_info.get_water_builder( rsd1 , heavyatom1 ) );
				{
					Size atom1 = rsd1_wb[iwat].atom1();
					Vector const & r1_atom1_xyz( rsd1.xyz( atom1 ) );

					numeric::xyzMatrix< Real >const & dwater_datom1 = rsd1_info.atom1_derivs()[iwat + atom1_offset];
					numeric::xyzMatrix< Real > ddome_datom1 = ddome_dwater * dwater_datom1;

					if ( atom1 == heavyatom1 ) ddome_datom1 += ddome_dbase;

					Vector ddome_datom1x ( ddome_datom1(1,1), ddome_datom1(2,1), ddome_datom1(3,1) );
					Vector ddome_datom1y ( ddome_datom1(1,2), ddome_datom1(2,2), ddome_datom1(3,2) );
					Vector ddome_datom1z ( ddome_datom1(1,3), ddome_datom1(2,3), ddome_datom1(3,3) );
					Vector dRdatom_times_R( f2d.dot( ddome_datom1x ), f2d.dot( ddome_datom1y ), f2d.dot( ddome_datom1z ) );

					Vector f2t = dRdatom_times_R;
					Vector f1t = r1_atom1_xyz.cross( r1_atom1_xyz - dRdatom_times_R ); // this is just dRdatom_times_R x r1_atom1_xyz

					r1_at_derivs[atom1].f1() += dbridgeE_dr_over_r * f1t;
					r1_at_derivs[atom1].f2() += dbridgeE_dr_over_r * f2t;
				}
				{
					Size atom2 = rsd1_wb[iwat].atom2();
					Vector const & r1_atom2_xyz( rsd1.xyz( atom2 ) );

					numeric::xyzMatrix< Real >const & dwater_datom2 = rsd1_info.atom2_derivs()[iwat + atom1_offset];
					numeric::xyzMatrix< Real > ddome_datom2 = ddome_dwater * dwater_datom2;

					if ( atom2 == heavyatom1 ) ddome_datom2 += ddome_dbase;

					Vector ddome_datom2x ( ddome_datom2(1,1), ddome_datom2(2,1), ddome_datom2(3,1) );
					Vector ddome_datom2y ( ddome_datom2(1,2), ddome_datom2(2,2), ddome_datom2(3,2) );
					Vector ddome_datom2z ( ddome_datom2(1,3), ddome_datom2(2,3), ddome_datom2(3,3) );
					Vector dRdatom_times_R( f2d.dot( ddome_datom2x ), f2d.dot( ddome_datom2y ), f2d.dot( ddome_datom2z ) );

					Vector f2t = dRdatom_times_R;
					Vector f1t = r1_atom2_xyz.cross( r1_atom2_xyz - dRdatom_times_R );

					r1_at_derivs[atom2].f1() += dbridgeE_dr_over_r * f1t;
					r1_at_derivs[atom2].f2() += dbridgeE_dr_over_r * f2t;
				}
				{
					Size atom3 = rsd1_wb[iwat].atom3();
					Vector const & r1_atom3_xyz( rsd1.xyz( atom3 ) );

					numeric::xyzMatrix< Real >const & dwater_datom3 = rsd1_info.atom3_derivs()[iwat + atom1_offset];
					numeric::xyzMatrix< Real > ddome_datom3 = ddome_dwater * dwater_datom3;

					if ( atom3 == heavyatom1 ) ddome_datom3 += ddome_dbase;

					Vector ddome_datom3x ( ddome_datom3(1,1), ddome_datom3(2,1), ddome_datom3(3,1) );
					Vector ddome_datom3y ( ddome_datom3(1,2), ddome_datom3(2,2), ddome_datom3(3,2) );
					Vector ddome_datom3z ( ddome_datom3(1,3), ddome_datom3(2,3), ddome_datom3(3,3) );
					Vector dRdatom_times_R( f2d.dot( ddome_datom3x ), f2d.dot( ddome_datom3y ), f2d.dot( ddome_datom3z ) );

					Vector f2t = dRdatom_times_R;
					Vector f1t = r1_atom3_xyz.cross( r1_atom3_xyz - dRdatom_times_R );

					r1_at_derivs[atom3].f1() += dbridgeE_dr_over_r * f1t;
					r1_at_derivs[atom3].f2() += dbridgeE_dr_over_r * f2t;
				}



			}
		}
	}
	if ( lk_ball_bridge2_wt!=0 || lk_ball_bridge_uncpl2_wt!=0 ) {


		for ( Size iwat=1; iwat <= atom1_n_waters; ++iwat ) {
			Vector const & atom1_water_xyz( res1_waters[ iwat + atom1_offset] );

			if ( atom1_water_xyz.distance_squared( atom2_xyz) > water_atom_range2 ) continue;

			Real avail = get_avail( rsd1_info.water_occlusions()[heavyatom1][iwat] );

			// Vector const & base = atom1_xyz;
			Vector const & water = atom1_water_xyz;

			for ( Size jwat=1; jwat <= atom2_n_waters; ++jwat ) {

				Vector const & atom2_water_xyz( res2_waters[ jwat + atom2_offset ] );

				Real dist2 = water.distance_squared(atom2_water_xyz);

				Real d2_delta = dist2;

				Real dlk_frac_dr_over_r = 0;

				if ( d2_delta < ball_overlap_width_A2_ && d2_delta > 0 ) {
					dlk_frac_dr_over_r = eval_d_lk_fraction_dr_over_r( d2_delta, ball_overlap_width_A2_ );
				}

				Real dbridge_frac_dr_over_r = avail * dlk_frac_dr_over_r;

				if ( dbridge_frac_dr_over_r == 0 ) continue;


				Real dbridgeE_dr_over_r = weight_factor * dbridge_frac_dr_over_r * ( lk_ball_bridge_uncpl2_wt + lk_ball_bridge2_wt * sol_value);


				Vector f2w = water - atom2_water_xyz;

				// R = sqrt(normsq( W1(K) - W2(L) ) )

				// dR = sqrt`(normsq( W1(K) - W2(L) ) ) * normsq`( W1(K) - W2(L) ) * d( W1(K) - W2(L) )
				// dR = 1/2/norm * 2( W1(K) - W2(L) ) * d( W1(K) - W2(L) )
				// dR_times_R = ( W1(K) - W2(L) ) * d( W1(K) - W2(L) )

				// Fixed L
				// dR_times_R = ( W1(K) - W2(L) ) * dW1_dK *dK
				// dRdK_times_R = ( W1(K) - W2(L) ) * dW1_dK

				// Fixed K
				// dR_times_R = -( W1(K) - W2(L) ) * dW2_dL *dL
				// dRdL_times_R = -( W1(K) - W2(L) ) * dW2_dL * dL



				WaterBuilders const & rsd2_wb( rsd2_info.get_water_builder( rsd2 , heavyatom2 ) );
				{
					Size atom1 = rsd2_wb[jwat].atom1();
					Vector const & r2_atom1_xyz( rsd2.xyz( atom1 ) );

					numeric::xyzMatrix< Real >const & dwater_datom1 = rsd2_info.atom1_derivs()[jwat + atom2_offset];

					Vector dwater_datom1x ( dwater_datom1(1,1), dwater_datom1(2,1), dwater_datom1(3,1) );
					Vector dwater_datom1y ( dwater_datom1(1,2), dwater_datom1(2,2), dwater_datom1(3,2) );
					Vector dwater_datom1z ( dwater_datom1(1,3), dwater_datom1(2,3), dwater_datom1(3,3) );
					Vector dRdatom_times_R( f2w.dot( dwater_datom1x ), f2w.dot( dwater_datom1y ), f2w.dot( dwater_datom1z ) );

					Vector f2t = dRdatom_times_R;
					Vector f1t = r2_atom1_xyz.cross( r2_atom1_xyz - dRdatom_times_R ); // this is just dRdatom_times_R x r1_atom1_xyz

					r2_at_derivs[atom1].f1() -= dbridgeE_dr_over_r * f1t;
					r2_at_derivs[atom1].f2() -= dbridgeE_dr_over_r * f2t;
				}
				{
					Size atom2 = rsd2_wb[jwat].atom2();
					Vector const & r2_atom2_xyz( rsd2.xyz( atom2 ) );

					numeric::xyzMatrix< Real >const & dwater_datom2 = rsd2_info.atom2_derivs()[jwat + atom2_offset];

					Vector dwater_datom2x ( dwater_datom2(1,1), dwater_datom2(2,1), dwater_datom2(3,1) );
					Vector dwater_datom2y ( dwater_datom2(1,2), dwater_datom2(2,2), dwater_datom2(3,2) );
					Vector dwater_datom2z ( dwater_datom2(1,3), dwater_datom2(2,3), dwater_datom2(3,3) );
					Vector dRdatom_times_R( f2w.dot( dwater_datom2x ), f2w.dot( dwater_datom2y ), f2w.dot( dwater_datom2z ) );

					Vector f2t = dRdatom_times_R;
					Vector f1t = r2_atom2_xyz.cross( r2_atom2_xyz - dRdatom_times_R );

					r2_at_derivs[atom2].f1() -= dbridgeE_dr_over_r * f1t;
					r2_at_derivs[atom2].f2() -= dbridgeE_dr_over_r * f2t;
				}
				{
					Size atom3 = rsd2_wb[jwat].atom3();
					Vector const & r2_atom3_xyz( rsd2.xyz( atom3 ) );

					numeric::xyzMatrix< Real >const & dwater_datom3 = rsd2_info.atom3_derivs()[jwat + atom2_offset];

					Vector dwater_datom3x ( dwater_datom3(1,1), dwater_datom3(2,1), dwater_datom3(3,1) );
					Vector dwater_datom3y ( dwater_datom3(1,2), dwater_datom3(2,2), dwater_datom3(3,2) );
					Vector dwater_datom3z ( dwater_datom3(1,3), dwater_datom3(2,3), dwater_datom3(3,3) );
					Vector dRdatom_times_R( f2w.dot( dwater_datom3x ), f2w.dot( dwater_datom3y ), f2w.dot( dwater_datom3z ) );

					Vector f2t = dRdatom_times_R;
					Vector f1t = r2_atom3_xyz.cross( r2_atom3_xyz - dRdatom_times_R );

					r2_at_derivs[atom3].f1() -= dbridgeE_dr_over_r * f1t;
					r2_at_derivs[atom3].f2() -= dbridgeE_dr_over_r * f2t;
				}



				WaterBuilders const & rsd1_wb( rsd1_info.get_water_builder( rsd1 , heavyatom1 ) );
				{
					Size atom1 = rsd1_wb[iwat].atom1();
					Vector const & r1_atom1_xyz( rsd1.xyz( atom1 ) );

					numeric::xyzMatrix< Real >const & dwater_datom1 = rsd1_info.atom1_derivs()[iwat + atom1_offset];

					Vector dwater_datom1x ( dwater_datom1(1,1), dwater_datom1(2,1), dwater_datom1(3,1) );
					Vector dwater_datom1y ( dwater_datom1(1,2), dwater_datom1(2,2), dwater_datom1(3,2) );
					Vector dwater_datom1z ( dwater_datom1(1,3), dwater_datom1(2,3), dwater_datom1(3,3) );
					Vector dRdatom_times_R( f2w.dot( dwater_datom1x ), f2w.dot( dwater_datom1y ), f2w.dot( dwater_datom1z ) );

					Vector f2t = dRdatom_times_R;
					Vector f1t = r1_atom1_xyz.cross( r1_atom1_xyz - dRdatom_times_R ); // this is just dRdatom_times_R x r1_atom1_xyz

					r1_at_derivs[atom1].f1() += dbridgeE_dr_over_r * f1t;
					r1_at_derivs[atom1].f2() += dbridgeE_dr_over_r * f2t;
				}
				{
					Size atom2 = rsd1_wb[iwat].atom2();
					Vector const & r1_atom2_xyz( rsd1.xyz( atom2 ) );

					numeric::xyzMatrix< Real >const & dwater_datom2 = rsd1_info.atom2_derivs()[iwat + atom1_offset];

					Vector dwater_datom2x ( dwater_datom2(1,1), dwater_datom2(2,1), dwater_datom2(3,1) );
					Vector dwater_datom2y ( dwater_datom2(1,2), dwater_datom2(2,2), dwater_datom2(3,2) );
					Vector dwater_datom2z ( dwater_datom2(1,3), dwater_datom2(2,3), dwater_datom2(3,3) );
					Vector dRdatom_times_R( f2w.dot( dwater_datom2x ), f2w.dot( dwater_datom2y ), f2w.dot( dwater_datom2z ) );

					Vector f2t = dRdatom_times_R;
					Vector f1t = r1_atom2_xyz.cross( r1_atom2_xyz - dRdatom_times_R );

					r1_at_derivs[atom2].f1() += dbridgeE_dr_over_r * f1t;
					r1_at_derivs[atom2].f2() += dbridgeE_dr_over_r * f2t;
				}
				{
					Size atom3 = rsd1_wb[iwat].atom3();
					Vector const & r1_atom3_xyz( rsd1.xyz( atom3 ) );

					numeric::xyzMatrix< Real >const & dwater_datom3 = rsd1_info.atom3_derivs()[iwat + atom1_offset];

					Vector dwater_datom3x ( dwater_datom3(1,1), dwater_datom3(2,1), dwater_datom3(3,1) );
					Vector dwater_datom3y ( dwater_datom3(1,2), dwater_datom3(2,2), dwater_datom3(3,2) );
					Vector dwater_datom3z ( dwater_datom3(1,3), dwater_datom3(2,3), dwater_datom3(3,3) );
					Vector dRdatom_times_R( f2w.dot( dwater_datom3x ), f2w.dot( dwater_datom3y ), f2w.dot( dwater_datom3z ) );

					Vector f2t = dRdatom_times_R;
					Vector f1t = r1_atom3_xyz.cross( r1_atom3_xyz - dRdatom_times_R );

					r1_at_derivs[atom3].f1() += dbridgeE_dr_over_r * f1t;
					r1_at_derivs[atom3].f2() += dbridgeE_dr_over_r * f2t;
				}



			}
		}
	}

	// Real bridging_lk_factor = lk_ball_bridge_uncpl_weight + lk_ball_bridge_weight * (lk_score+other_lk_score);
	// Real weighted_d2_water_delta(0), angleterm_lkbr(0), pointterm_lkbr(0), d_angleterm_lkbr_dr(0);

	// //utility::vector1< numeric::xyzVector<core::Real> > d_weighted_d2_d_di;
	// WaterDerivVectors d_weighted_d2_d_di;

	// Real const lkbr_fraction( get_lkbr_fractional_contribution(
	//     heavyatom1_xyz, atom2_xyz, atom1_n_waters, atom2_n_waters, atom1_waters, atom2_waters,
	//     d_weighted_d2_d_di, weighted_d2_water_delta, pointterm_lkbr, angleterm_lkbr, d_angleterm_lkbr_dr
	//     ) );

	// // A: change in LK (used as a scalefactor) on _bridge but not _bridge_uncpl
	// if ( lk_ball_bridge_weight != 0 ) {
	//     Real const dE_dr_over_r( weight_factor * lk_ball_bridge_weight * lkbr_fraction * lk_deriv * inv_dis );
	//     r1_at_derivs[heavyatom1].f1() += dE_dr_over_r * f1;
	//     r1_at_derivs[heavyatom1].f2() += dE_dr_over_r * f2;
	//     r2_at_derivs[heavyatom2].f1() -= dE_dr_over_r * f1;
	//     r2_at_derivs[heavyatom2].f2() -= dE_dr_over_r * f2;
	// }

	// // A': water-angle potential (also used as a scalefactor)
	// if ( lkbridge_angle_widthscale_!=0 ) {
	//     Real const dE_dr_over_r( -0.25*weight_factor * bridging_lk_factor * pointterm_lkbr * d_angleterm_lkbr_dr * inv_dis );
	//     r1_at_derivs[heavyatom1].f1() += dE_dr_over_r * f1;
	//     r1_at_derivs[heavyatom1].f2() += dE_dr_over_r * f2;
	//     r2_at_derivs[heavyatom2].f1() -= dE_dr_over_r * f1;
	//     r2_at_derivs[heavyatom2].f2() -= dE_dr_over_r * f2;
	// }

	// // B: change in water positions
	// if ( weighted_d2_water_delta < overlap_width_A2_ && weighted_d2_water_delta > 0.0 ) {
	//     for ( Size i=1; i<=atom1_n_waters; ++i ) {
	//         // what is the derivative of the lkbr_fraction term wrt r?
	//         numeric::xyzVector<core::Real> const dE_dwi
	//             ( weight_factor * bridging_lk_factor * d_weighted_d2_d_di[i] * angleterm_lkbr * eval_d_lk_fraction_dr_over_r( weighted_d2_water_delta, overlap_width_A2_ ) );

	//         // derivatives for the desolvated atoms
	//         // dR/datom1 = dR/dwater * dwater/datom1
	//         WaterBuilders const & rsd1_wb( rsd1_info.get_water_builder( rsd1 , heavyatom1 ) );
	//         {
	//             Size atom1 = rsd1_wb[i].atom1();
	//             Vector const & r1_atom1_xyz( rsd1.xyz( atom1 ) );

	//             numeric::xyzMatrix< Real >const & dwater_datom1 = rsd1_info.atom1_derivs()[heavyatom1][i];
	//             Vector dwater_datom1x ( dwater_datom1(1,1), dwater_datom1(2,1), dwater_datom1(3,1) );
	//             Vector dwater_datom1y ( dwater_datom1(1,2), dwater_datom1(2,2), dwater_datom1(3,2) );
	//             Vector dwater_datom1z ( dwater_datom1(1,3), dwater_datom1(2,3), dwater_datom1(3,3) );
	//             Vector dRdatom( dE_dwi.dot( dwater_datom1x ), dE_dwi.dot( dwater_datom1y ), dE_dwi.dot( dwater_datom1z ) );

	//             Vector f2t = dRdatom;
	//             Vector f1t = r1_atom1_xyz.cross( r1_atom1_xyz - dRdatom );

	//             r1_at_derivs[atom1].f1() += f1t;
	//             r1_at_derivs[atom1].f2() += f2t;
	//         }
	//         {
	//             Size atom2 = rsd1_wb[i].atom2();
	//             Vector const & r1_atom2_xyz( rsd1.xyz( atom2 ) );

	//             numeric::xyzMatrix< Real >const & dwater_datom2 = rsd1_info.atom2_derivs()[heavyatom1][i];
	//             Vector dwater_datom2x ( dwater_datom2(1,1), dwater_datom2(2,1), dwater_datom2(3,1) );
	//             Vector dwater_datom2y ( dwater_datom2(1,2), dwater_datom2(2,2), dwater_datom2(3,2) );
	//             Vector dwater_datom2z ( dwater_datom2(1,3), dwater_datom2(2,3), dwater_datom2(3,3) );
	//             Vector dRdatom( dE_dwi.dot( dwater_datom2x ), dE_dwi.dot( dwater_datom2y ), dE_dwi.dot( dwater_datom2z ) );

	//             Vector f2t = dRdatom;
	//             Vector f1t = r1_atom2_xyz.cross( r1_atom2_xyz - dRdatom );

	//             r1_at_derivs[atom2].f1() += f1t;
	//             r1_at_derivs[atom2].f2() += f2t;
	//         }
	//         {
	//             Size atom3 = rsd1_wb[i].atom3();
	//             Vector const & r1_atom3_xyz( rsd1.xyz( atom3 ) );

	//             numeric::xyzMatrix< Real >const & dwater_datom3 = rsd1_info.atom3_derivs()[heavyatom1][i];
	//             Vector dwater_datom3x ( dwater_datom3(1,1), dwater_datom3(2,1), dwater_datom3(3,1) );
	//             Vector dwater_datom3y ( dwater_datom3(1,2), dwater_datom3(2,2), dwater_datom3(3,2) );
	//             Vector dwater_datom3z ( dwater_datom3(1,3), dwater_datom3(2,3), dwater_datom3(3,3) );
	//             Vector dRdatom( dE_dwi.dot( dwater_datom3x ), dE_dwi.dot( dwater_datom3y ), dE_dwi.dot( dwater_datom3z ) );

	//             Vector f2t = dRdatom;
	//             Vector f1t = r1_atom3_xyz.cross( r1_atom3_xyz - dRdatom );

	//             r1_at_derivs[atom3].f1() += f1t;
	//             r1_at_derivs[atom3].f2() += f2t;
	//         }
	//     }
	// }
}





bool
LK_DomeEnergy::defines_score_for_residue_pair(
	conformation::Residue const & ,
	conformation::Residue const & ,
	bool
) const {
	if ( packing_ ) {
		return false;
	} else {
		return true;
	}
}

inline
LKD_ResidueInfo const &
retrieve_lkd_resdata(
	conformation::Residue const & res
)
{
	using namespace core::conformation::residue_datacache;
	debug_assert( utility::pointer::dynamic_pointer_cast< LKD_ResidueInfo const > ( res.data().get_const_ptr( LK_DOME_INFO )));
	return ( static_cast< LKD_ResidueInfo const & > ( res.data().get( LK_DOME_INFO ) ) );
}



void
LK_DomeEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & , //pose,
	ScoreFunction const & , //scorefxn,
	kinematics::MinimizerMapBase const &, // min_map,
	ResSingleMinimizationData const &,// res1data,
	ResSingleMinimizationData const &,// res2data,
	ResPairMinimizationData & pair_data
) const
{
	using namespace etable::count_pair;


	CPCrossoverBehavior crossover = (rsd1.is_polymer_bonded(rsd2) && rsd2.is_polymer_bonded(rsd1))? CP_CROSSOVER_4 : CP_CROSSOVER_3;
	CountPairFunctionOP cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, crossover );

	ResiduePairNeighborListOP nblist( utility::pointer::static_pointer_cast< core::scoring::ResiduePairNeighborList > ( pair_data.get_data( lkdome_nblist ) ));
	if ( ! nblist ) nblist = utility::pointer::make_shared<ResiduePairNeighborList>();


	Real const XX2 = atomic_interaction_cutoff()*atomic_interaction_cutoff();
	nblist->initialize_from_residues( XX2, 0, 0, rsd1, rsd2, cpfxn );

	pair_data.set_data( lkdome_nblist, nblist );

	LKD_ResPairMinDataOP respair_data = utility::pointer::make_shared<LKD_ResPairMinData>();
	pair_data.set_data( lkd_respair_data, respair_data );
}


/// @note  Assumes that atom1 is the "moving" atom, ie the atom for which eval_atom_derivative was called
void
LK_DomeEnergy::sum_deriv_contributions_for_heavyatom_pair(
	Real const d2,
	Size const heavyatom1,
	conformation::Residue const & rsd1,
	LKD_ResidueInfo const & rsd1_info,
	Size const heavyatom2,
	conformation::Residue const & rsd2,
	LKD_ResidueInfo const & rsd2_info,
	pose::Pose const &,
	EnergyMap const & weights,
	Real const cp_weight,
	utility::vector1< DerivVectorPair > & r1_at_derivs,
	utility::vector1< DerivVectorPair > & r2_at_derivs
) const
{
	sum_deriv_contributions_for_heavyatom_pair_one_way(
		heavyatom1, rsd1, rsd1_info, heavyatom2, rsd2, rsd2_info, weights,
		cp_weight, d2, r1_at_derivs, r2_at_derivs );

	sum_deriv_contributions_for_heavyatom_pair_one_way(
		heavyatom2, rsd2, rsd2_info, heavyatom1, rsd1, rsd1_info, weights,
		cp_weight, d2, r2_at_derivs, r1_at_derivs );
}


void
LK_DomeEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_at_derivs,
	utility::vector1< DerivVectorPair > & r2_at_derivs
) const
{

	// std::cout << "pair derivative " << rsd1.seqpos() << " " << rsd2.seqpos() << std::endl;

	debug_assert( r1_at_derivs.size() >= rsd1.natoms() );
	debug_assert( r2_at_derivs.size() >= rsd2.natoms() );

	// retrieve some info
	LKD_ResidueInfo const & rsd1_info( retrieve_lkd_resdata( rsd1 ) );
	LKD_ResidueInfo const & rsd2_info( retrieve_lkd_resdata( rsd2 ) );

	auto const & nblist =
		static_cast< ResiduePairNeighborList const & > (min_data.get_data_ref( lkdome_nblist ));

	// bool use_lkbr_uncpl = (weights[core::scoring::lk_ball_bridge_uncpl]!=0);

	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );
	for ( Size ii = 1, iiend = neighbs.size(); ii <= iiend; ++ii ) {
		Size const heavyatom1( neighbs[ ii ].atomno1() ), heavyatom2( neighbs[ ii ].atomno2() );
		if ( rsd1.atom_is_hydrogen( heavyatom1 ) || rsd2.atom_is_hydrogen( heavyatom2 ) ) continue;
		Real const cp_weight( neighbs[ ii ].weight() );

		Real const d2( rsd1.xyz( heavyatom1 ).distance_squared( rsd2.xyz( heavyatom2 ) ) );

		if ( d2 == Real(0.0) ) continue; // sanity check
		// if ( !use_lkbr_uncpl && d2 >= fasol_max_dis2_ ) continue;
		//if ( d2 >= lkb_max_dis2_ ) continue;

		//std::cout << "deriv: " << rsd1.seqpos() << " " << rsd2.seqpos() << std::endl;

		sum_deriv_contributions_for_heavyatom_pair(
			d2, heavyatom1, rsd1, rsd1_info, heavyatom2, rsd2, rsd2_info, pose, weights, cp_weight, r1_at_derivs, r2_at_derivs );
	}
}


bool
LK_DomeEnergy::defines_intrares_energy( EnergyMap const & ) const {
	return false;
}

void
LK_DomeEnergy::eval_intrares_energy(
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const {

}


void
LK_DomeEnergy::setup_for_packing( pose::Pose &, utility::vector1< bool > const &, utility::vector1< bool > const & ) const {
	packing_ = true;
}

void
LK_DomeEnergy::setup_for_minimizing( pose::Pose & pose, ScoreFunction const & sfxn, kinematics::MinimizerMapBase const &) const {
	packing_ = false;
	minimizing_ = true;

	// std::cout << "Setup for min" << std::endl;


	// We need this before we can do the next step
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		this->my_setup_for_scoring_for_residue( pose.residue( ii ), pose, sfxn, pose.residue_data( ii ) );
	}

	fill_static_occlusions( pose );


}

void
LK_DomeEnergy::finalize_after_minimizing(
	pose::Pose &
) const {
	minimizing_ = false;
}


void
LK_DomeEnergy::update_cached_lkb_resinfo(
	conformation::Residue const & rsd,
	basic::datacache::BasicDataCache & residue_data_cache,
	bool compute_derivs,
	bool clear_occlusions
) const {
	using conformation::residue_datacache::LK_DOME_INFO;
	if ( residue_data_cache.has( LK_DOME_INFO ) ) {
		debug_assert( utility::pointer::dynamic_pointer_cast< LKD_ResidueInfo > ( residue_data_cache.get_ptr( LK_DOME_INFO )));
		auto & info( static_cast< LKD_ResidueInfo & > ( residue_data_cache.get( LK_DOME_INFO )));
		info.build_waters( rsd, compute_derivs ); // update the coordinates for the existing lkb-resinfo object
		if ( clear_occlusions ) info.clear_occlusions();
	} else {
		LKD_ResidueInfoOP info = utility::pointer::make_shared<LKD_ResidueInfo>( rsd, this );
		residue_data_cache.set( LK_DOME_INFO, info );
	}
}

bool
LK_DomeEnergy::requires_a_setup_for_scoring_for_residue_opportunity_during_regular_scoring( pose::Pose const & ) const
{
	return false;
}

void
LK_DomeEnergy::my_setup_for_scoring_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const &,
	ScoreFunction const &,
	basic::datacache::BasicDataCache & residue_data_cache
) const
{
	update_cached_lkb_resinfo( rsd, residue_data_cache, false, true );
}

bool
LK_DomeEnergy::requires_a_setup_for_derivatives_for_residue_opportunity( pose::Pose const &  ) const
{
	return true;
}

void
LK_DomeEnergy::setup_for_derivatives_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const &,
	ScoreFunction const &,
	ResSingleMinimizationData &,
	basic::datacache::BasicDataCache & residue_data_cache
) const
{
	// compute water locations and the water-positional derivatives
	update_cached_lkb_resinfo( rsd, residue_data_cache, true, false );
}

void
LK_DomeEnergy::setup_for_scoring(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn
) const {
	packing_ = false;


	// std::cout << "Setup for scoring " << (minimizing_ ? "safe" : "unsafe" ) << std::endl;

	// This can't happen by chance. (sqrt(3) sqrt(5) sqrt(6) sqrt(7)) + 0.005
	if (   std::abs( sfxn.weights()[core::scoring::lk_dome] - 1.737050808 ) < 0.000000001
			&& std::abs( sfxn.weights()[core::scoring::lk_dome_iso] - 2.241067977 ) < 0.000000001
			&& std::abs( sfxn.weights()[core::scoring::lk_dome_bridge] - 2.454489743 ) < 0.000000001
			&& std::abs( sfxn.weights()[core::scoring::lk_dome_bridge_uncpl] - 2.650751311 ) < 0.000000001
			) {
		debug_disable_count_pair_ = true;
		TR.Warning << "Countpair disabled! You didn't pick the lk_dome* weights by accident right?" << std::endl;
	}

	// We don't recalculate occlusion during minimizing so that it stays 2-body
	if ( ! minimizing_ ) {
		// We need this before we can do the next step
		for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			this->my_setup_for_scoring_for_residue( pose.residue( ii ), pose, sfxn, pose.residue_data( ii ) );
		}

		fill_static_occlusions( pose );
	} else {
		for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			update_cached_lkb_resinfo( pose.residue(ii), pose.residue_data( ii ), false, false );
		}
	}

}

void
LK_DomeEnergy::fill_static_occlusions( pose::Pose & pose ) const {

	using conformation::residue_datacache::LK_DOME_INFO;

	// This will have a neighbor distance more than we need
	EnergyGraph & energy_graph( pose.energies().energy_graph() );

	for ( Size i=1, i_end = pose.size(); i<= i_end; ++i ) {
		conformation::Residue const & resl( pose.residue( i ) );

		LKD_ResidueInfo & info_l = *( static_cast< LKD_ResidueInfo * >( pose.residue_data( i ).get_raw_ptr( LK_DOME_INFO )));

		for ( utility::graph::Graph::EdgeListIter
				iru  = energy_graph.get_node(i)->upper_edge_list_begin(),
				irue = energy_graph.get_node(i)->upper_edge_list_end();
				iru != irue; ++iru ) {
			auto & edge( static_cast< EnergyEdge & > (**iru) );

			Size const j( edge.get_second_node_ind() );
			conformation::Residue const & resu( pose.residue( j ) );

			LKD_ResidueInfo & info_u = *( static_cast< LKD_ResidueInfo * >( pose.residue_data( j ).get_raw_ptr( LK_DOME_INFO )));

			fill_occlusions_1way( resl, info_l, resu );
			fill_occlusions_1way( resu, info_u, resl );

		}

	}
}


void
LK_DomeEnergy::fill_occlusions_1way(
	conformation::Residue const & rsd1,
	LKD_ResidueInfo & rsd1_info,
	conformation::Residue const & rsd2
) const {

	if ( ! rsd1_info.has_waters() ) return;

	WaterCoords const & rsd1_waters( rsd1_info.waters() );
	utility::vector1< WaterOcclusions > & rsd1_occlusions( rsd1_info.water_occlusions() );

	for ( Size atom1 = 1; atom1 <= rsd1.nheavyatoms(); ++atom1 ) {
		Size const n_atom1_waters = rsd1_info.n_attached_waters()[ atom1 ];
		if ( n_atom1_waters == 0 ) continue;

		// WaterCoords const & atom1_waters( rsd1_waters[ atom1 ] );
		Size atom1_offset = rsd1_info.water_offset_for_atom( atom1 );

		for ( Size atom2 = 1; atom2 <= rsd2.nheavyatoms(); ++atom2 ) {
			Vector const & atom2_xyz( rsd2.xyz( atom2 ) );
			Size const atom2_type_index( rsd2.atom( atom2 ).type() );

			for ( Size i_water = 1; i_water <= n_atom1_waters; i_water++ ) {
				rsd1_occlusions[atom1][i_water] += lk_ball_->get_lk_fractional_contribution_for_single_water( atom2_xyz, atom2_type_index,
					rsd1_waters[i_water + atom1_offset ] );
			}

		}
	}

}


}
}
}
