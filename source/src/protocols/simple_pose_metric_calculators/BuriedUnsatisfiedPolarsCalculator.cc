// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// How many buried unsatisfied polars are there?
///
/// @details significantly updated in 2017 to have more generous definition of h-bonds
///  (previously, legit h-bonds were excluded because of sfxn exceptions); users can now choose
///  between legacy SASA and VSASA for burial; poses with more than 3 chains now supported; the
///  way unsats are counted and reported is now different (before, all unsats were counted as equal,
///  which is misleading); users can choose different reporting behaviours; legacy options can be
///  restored by setting legacy=true, but this is only recommended for benchmarking purposes
///
/// (not sure where these numbers came from...they don't really make sense...but leaving them in for now)
/// Buried unsatisfied polar hbonds are destabilizing for proteins. It is good to have less.
/// In a study of 2299 high resolution crystal structures of 1.5A or better, there was an average
/// 71 unsatisfied buried polar hbonds. The normalized average (normalized against aa #) was 0.30 (unpublished).
/// To get this piece of code to work, you must first load in your pdb. Then, you need the following lines:
///
/// core::pose::metrics::PoseMetricCalculatorOP sasa_calculator = new protocols::toolbox::PoseMetricCalculators::SasaCalculatorLegacy;
/// core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );
///
/// core::pose::metrics::PoseMetricCalculatorOP num_hbonds_calculator = new protocols::toolbox::PoseMetricCalculators::NumberHBondsCalculator();
/// core::pose::metrics::CalculatorFactory::Instance().register_calculator( "num_hbonds", num_hbonds_calculator );
///
/// core::pose::metrics::PoseMetricCalculatorOP unsat_calculator = new protocols::toolbox::PoseMetricCalculators::BuriedUnsatisfiedPolarsCalculator("sasa", "num_hbonds");
/// core::pose::metrics::CalculatorFactory::Instance().register_calculator( "unsat", unsat_calculator );
///
/// This segment of code sets everything up to be used in the calculator. To use this on your protein, you simply need to
/// write the following: pose.print_metric("unsat", "all_bur_unsat_polars");
///
/// @author
/// Florian Richter
/// Steven Combs - comments
/// @author Scott Boyken (sboyken@gmail.com), major updates and refactoring


// Unit headers
#include <protocols/simple_pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/Pose.hh>
#include <protocols/simple_pose_metric_calculators/NumberHBondsCalculator.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <protocols/vardist_solaccess/VarSolDRotamerDots.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/LayerSelector.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/packing/surf_vol.hh>
#include <core/pose/util.tmpl.hh>
#include <core/chemical/AtomType.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <utility/assert.hh>
#include <utility/vector1.hh>
#include <basic/MetricValue.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>
#include <basic/options/keys/bunsat_calc2.OptionKeys.gen.hh>
#include <basic/options/keys/holes.OptionKeys.gen.hh>

using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;

static basic::Tracer TR( "protocols.toolbox.PoseMetricCalculators.BuriedUnsatisfiedPolarsCalculator" );

#ifdef    SERIALIZATION
// Project serialization headers
#include <core/id/AtomID_Map.srlz.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace simple_pose_metric_calculators {

BuriedUnsatisfiedPolarsCalculator::BuriedUnsatisfiedPolarsCalculator(
	std::string const & sasa_calc,
	std::string const & hbond_calc,
	core::Real const burial_cutoff, // must be 3rd
	bool const generous, /* true */
	bool const legacy, /* false */
	bool const vsasa /* true */
) : all_bur_unsat_polars_( 0 ),
	bb_heavy_unsats_(0),
	all_heavy_unsats_(0),
	countable_nonheavy_unsats_(0),
	//special_region_bur_unsat_polars_(0),
	name_of_hbond_calc_( hbond_calc ),
	name_of_sasa_calc_( sasa_calc ),
	special_region_( /* NULL */ ),
	special_region_entire_residue_( true ),
	burial_cutoff_( burial_cutoff ),
	probe_radius_( basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius] ),
	residue_surface_cutoff_( 45.0 ),
	generous_hbonds_( generous ),
	legacy_counting_( legacy ),
	vsasa_( vsasa ),
	use_sc_neighbors_( false ),
	skip_surface_res_( false )
{
	if ( burial_cutoff < 0.0 ) {
		burial_cutoff_ = ( vsasa_ ) ? basic::options::option[basic::options::OptionKeys::bunsat_calc2::sasa_burial_cutoff] : basic::options::option[basic::options::OptionKeys::pose_metrics::atomic_burial_cutoff];
	}
	if ( vsasa_ ) residue_surface_cutoff_ = 20.0;
	atom_bur_unsat_.clear();
	residue_bur_unsat_polars_.clear();
	special_region_.clear();
	assert_calculators();
}


BuriedUnsatisfiedPolarsCalculator::BuriedUnsatisfiedPolarsCalculator(
	std::string const & sasa_calc,
	std::string const & hbond_calc,
	std::set< core::Size > const & special_region,
	core::Real const burial_cutoff, // must be 4th
	//core::Real const probe_r,
	//core::Real const residue_surface_cutoff, /* 45.0 */
	bool const generous, /* true */
	bool const legacy, /* false */
	bool const vsasa /* true */
) : all_bur_unsat_polars_(0),
	bb_heavy_unsats_(0),
	all_heavy_unsats_(0),
	countable_nonheavy_unsats_(0),
	name_of_hbond_calc_( hbond_calc ),
	name_of_sasa_calc_( sasa_calc ),
	// special_region_( special_region ),
	special_region_entire_residue_( true ),
	burial_cutoff_( burial_cutoff ),
	probe_radius_( basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius] ),
	residue_surface_cutoff_( 45.0 ),
	generous_hbonds_( generous ),
	legacy_counting_( legacy ),
	vsasa_( vsasa ),
	use_sc_neighbors_( false ),
	skip_surface_res_( false )
{
	set_special_region( special_region );
	if ( burial_cutoff < 0.0 ) {
		burial_cutoff_ = ( vsasa_ ) ? basic::options::option[basic::options::OptionKeys::bunsat_calc2::sasa_burial_cutoff] : basic::options::option[basic::options::OptionKeys::pose_metrics::atomic_burial_cutoff];
	}
	if ( vsasa_ ) residue_surface_cutoff_ = 20.0;
	atom_bur_unsat_.clear();
	residue_bur_unsat_polars_.clear();
	assert_calculators();
}

BuriedUnsatisfiedPolarsCalculator::BuriedUnsatisfiedPolarsCalculator(
	std::string const & sasa_calc,
	std::string const & hbond_calc,
	core::id::AtomID_Map< bool > const & special_region,
	core::Real const burial_cutoff,
	core::Real const probe_r,
	core::Real const residue_surface_cutoff,
	bool const generous,
	bool const legacy,
	bool const vsasa,
	bool const use_sc_neighbors,
	bool const skip_surface_res
) : all_bur_unsat_polars_(0),
	bb_heavy_unsats_(0),
	all_heavy_unsats_(0),
	countable_nonheavy_unsats_(0),
	//special_region_bur_unsat_polars_(0),
	name_of_hbond_calc_( hbond_calc ),
	name_of_sasa_calc_( sasa_calc ),
	special_region_( special_region ),
	special_region_entire_residue_( false ),
	burial_cutoff_( burial_cutoff ),
	probe_radius_( probe_r ),
	residue_surface_cutoff_( residue_surface_cutoff ),
	generous_hbonds_( generous ),
	legacy_counting_( legacy ),
	vsasa_( vsasa ),
	use_sc_neighbors_( use_sc_neighbors ),
	skip_surface_res_( skip_surface_res )
{
	// params explicitly defined so no need to check
	atom_bur_unsat_.clear();
	residue_bur_unsat_polars_.clear();
	assert_calculators();
}


void
BuriedUnsatisfiedPolarsCalculator::assert_calculators()
{
	if ( !CalculatorFactory::Instance().check_calculator_exists( name_of_hbond_calc_ ) ) {
		name_of_hbond_calc_ = ( generous_hbonds_ ) ? "bur_unsat_calc_generous_hbond_calc" : "bur_unsat_calc_legacy_hbond_calc";
		// potentially will need instances of both if want to run in same pose or XML
		// default now uses generous h-bonds (counts all h-bonds); if generous_hbonds_=false, use lagacy behavior, which only counts h-bonds that sfxn respects
		if ( !CalculatorFactory::Instance().check_calculator_exists( name_of_hbond_calc_ ) ) {
			CalculatorFactory::Instance().register_calculator( name_of_hbond_calc_, PoseMetricCalculatorOP( new NumberHBondsCalculator( generous_hbonds_ ) ) );
		}
	}
	if ( name_of_sasa_calc_ == "dalphaball" ) {
		if ( basic::options::option[ basic::options::OptionKeys::holes::dalphaball ]() == std::string("") ) {
			utility_exit_with_message( "Error! You must set -holes:dalphaball to use dalphaball_sasa!!" );
		}
	} else if ( name_of_sasa_calc_ != "nondefault" && !CalculatorFactory::Instance().check_calculator_exists( name_of_sasa_calc_ ) ) {
		name_of_sasa_calc_ = ( vsasa_ ) ? "bur_unsat_calc_vsasa_calc" : "bur_unsat_calc_legacy_sasa_calc";
		if ( !CalculatorFactory::Instance().check_calculator_exists( name_of_sasa_calc_ ) ) {
			if ( vsasa_ ) { // new default 17/09/03
				CalculatorFactory::Instance().register_calculator( name_of_sasa_calc_, protocols::vardist_solaccess::VarSolDistSasaCalculatorOP( new protocols::vardist_solaccess::VarSolDistSasaCalculator() ) );
			} else {
				name_of_sasa_calc_ = "bur_unsat_calc_legacy_sasa_calc";
				CalculatorFactory::Instance().register_calculator( name_of_sasa_calc_, PoseMetricCalculatorOP( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy( probe_radius_ ) ) );
			}
		}
	}
}


void
BuriedUnsatisfiedPolarsCalculator::lookup(
	std::string const & key,
	basic::MetricValueBase * valptr
) const
{
	if ( key == "all_bur_unsat_polars" ) { // for Legacy behavior
		basic::check_cast( valptr, &all_bur_unsat_polars_, "all_bur_unsat_polars expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( all_bur_unsat_polars_ );

	} else if ( key == "bb_heavy_unsats" ) {
		basic::check_cast( valptr, &bb_heavy_unsats_, "bb_heavy_unsats expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( bb_heavy_unsats_ );

	} else if ( key == "all_heavy_unsats" ) {
		basic::check_cast( valptr, &all_heavy_unsats_, "all_heavy_unsats expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( all_heavy_unsats_ );

	} else if ( key == "countable_nonheavy_unsats" ) {
		basic::check_cast( valptr, &countable_nonheavy_unsats_, "hpol_unsats expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( countable_nonheavy_unsats_ );

	} else if ( key == "atom_bur_unsat" ) {
		basic::check_cast( valptr, &atom_bur_unsat_, "atom_bur_unsat expects to return a id::AtomID_Map< bool >" );
		(static_cast<basic::MetricValue<core::id::AtomID_Map< bool > > *>(valptr))->set( atom_bur_unsat_ );

	} else if ( key == "residue_bur_unsat_polars" ) {
		basic::check_cast( valptr, &residue_bur_unsat_polars_, "residue_bur_unsat_polars expects to return a utility::vector1< Size >" );
		(static_cast<basic::MetricValue<utility::vector1< core::Size > > *>(valptr))->set( residue_bur_unsat_polars_ );

	} else {
		basic::Error() << "BuriedUnsatisfiedPolarsCalculator cannot compute the requested metric " << key << std::endl;
		utility_exit();
	}

} //lookup


std::string
BuriedUnsatisfiedPolarsCalculator::print( std::string const & key ) const
{
	if ( key == "all_bur_unsat_polars" ) {
		return utility::to_string( all_bur_unsat_polars_ );
	} else if ( key == "bb_heavy_unsats" ) {
		return utility::to_string( bb_heavy_unsats_ );
	} else if ( key == "all_heavy_unsats" ) {
		return utility::to_string( all_heavy_unsats_ );
	} else if ( key == "countable_nonheavy_unsats" ) {
		return utility::to_string( countable_nonheavy_unsats_ );
	} else if ( key == "atom_bur_unsat" ) {
		basic::Error() << "core::id::AtomID_Map< bool > has no output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "residue_bur_unsat_polars" ) {
		return utility::to_string( residue_bur_unsat_polars_ );
	}

	basic::Error() << "BuriedUnsatisfiedPolarsCalculator cannot compute metric " << key << std::endl;
	utility_exit();
	return "";

} //print


/// @brief this function doesn't actually recompute anything by itself, but calls the
/// @brief two member calculators and then processes the information out of the two of them
void
BuriedUnsatisfiedPolarsCalculator::recompute( Pose const & this_pose )
{
	all_bur_unsat_polars_ = 0;
	bb_heavy_unsats_ = 0;
	all_heavy_unsats_ = 0;
	countable_nonheavy_unsats_ = 0;

	if ( this_pose.size() != residue_bur_unsat_polars_.size() ) {
		residue_bur_unsat_polars_.resize( this_pose.size() );

		pose::initialize_atomid_map( atom_bur_unsat_, this_pose.conformation(), false ); // allocate for all atoms of pose
	}
	basic::MetricValue< core::id::AtomID_Map< core::Size > > atom_hbonds;
	this_pose.metric( name_of_hbond_calc_, "atom_Hbonds", atom_hbonds );

	core::id::AtomID_Map< core::Real > atom_sasa;
	utility::vector1< core::Real > residue_sasa;
	core::select::residue_selector::ResidueSubset buried_residues( this_pose.size(), false );

	// set up burial // VSASA now added by assert_calculators()
	if ( use_sc_neighbors_ ) {
		// SASA uses atom cutoff; us_sc_neighbor uses residue cutoff
		// TODO need better solution but this catches cases where forgot to update cutoff for sc_neighbor (rather than sasa) values
		if ( burial_cutoff_ < 1.0 ) burial_cutoff_ = 4.4;
		if ( residue_surface_cutoff_ > 7.0 ) burial_cutoff_ = 2.0;
		core::select::residue_selector::LayerSelectorOP core_layer( new core::select::residue_selector::LayerSelector() );
		core_layer->set_layers( true, true, false ); // we'll check for unsats in core and boundary but not surface
		core_layer->set_use_sc_neighbors( true ); // now true by default
		core_layer->set_cutoffs( burial_cutoff_ /* core */, residue_surface_cutoff_ /* surface */ );
		buried_residues = core_layer->apply( this_pose );
	} else { // use SASA
		calculate_sasa( this_pose, atom_sasa, residue_sasa );
	}

	for ( core::Size resnum = 1; resnum <= this_pose.size(); ++resnum ) {

		if ( use_sc_neighbors_ && !(buried_residues[ resnum ]) ) continue;
		if ( !use_sc_neighbors_ && skip_surface_res_ && residue_sasa[ resnum ] > residue_surface_cutoff_ ) continue;

		residue_bur_unsat_polars_[resnum] = 0;

		for ( core::Size at = 1; at <= this_pose.residue( resnum ).nheavyatoms(); ++at ) {

			if ( ! special_region_.empty() ) {
				core::Size lookup_atom = at;
				if ( special_region_entire_residue_ ) lookup_atom = 1;
				if ( ! special_region_( resnum, lookup_atom ) ) continue;
			}

			core::id::AtomID atid( at, resnum );
			bool is_buried( use_sc_neighbors_ ); // if sc_neighbors and we made it this far, it's buried; else, set false and determine SASA below
			bool this_atom_bur_unsat(false);

			if ( this_pose.residue( resnum ).atom_type( at ).is_acceptor() || this_pose.residue( resnum ).atom_type( at ).is_donor() ) {

				// have to check for backbone N on proline ( this can go away once Npro no longer incorrectly typed as DONOR
				if ( this_pose.residue( resnum ).name1() == 'P' && this_pose.residue(resnum).atom_type( at ).atom_type_name() == "Npro" ) continue;

				if ( !use_sc_neighbors_ ) {
					Real cursasa =  atom_sasa[ atid ];

					//we also have to add up the sasas for the H attached to this atom if there are any
					for ( core::Size hindex = this_pose.residue( resnum ).attached_H_begin( at ); hindex<= this_pose.residue( resnum ).attached_H_end( at ); hindex++ ) {
						cursasa = cursasa + atom_sasa[ core::id::AtomID ( hindex, resnum ) ];
					}
					is_buried = ( cursasa < burial_cutoff_ );
				}

				if ( is_buried ) {

					if ( legacy_counting_ ) {
						Size satisfac_cut = satisfaction_cutoff( this_pose.residue( resnum ).type().atom_type( at ).name() );
						Size bonded_heavyatoms = this_pose.residue( resnum ).n_bonded_neighbor_all_res( at ) - this_pose.residue( resnum ).type().number_bonded_hydrogens( at );
						if ( ( bonded_heavyatoms + atom_hbonds.value()[ atid ] ) < satisfac_cut ) {
							all_bur_unsat_polars_++;
							residue_bur_unsat_polars_[ resnum ]++;
							this_atom_bur_unsat = true;

							if ( this_pose.residue( resnum ).atom_is_backbone( at ) ) {
								bb_heavy_unsats_++;
							}
						}
					} else if ( this_pose.residue( resnum ).atom_type( at ).is_donor() && this_pose.residue( resnum ).atomic_charge( at ) != 0.0
							&& this_pose.residue( resnum ).atom_type(at).name() != "OH" ) {

						core::Size h_count(0);
						core::Size h_unsat(0);
						for ( core::Size hatm = this_pose.residue(resnum).attached_H_begin(at); hatm <= this_pose.residue(resnum).attached_H_end(at); ++hatm ) {
							h_count++;
							core::id::AtomID const hatm_id( hatm, resnum);

							if ( atom_hbonds.value().has( hatm_id ) && atom_hbonds.value()[ hatm_id ] == 0 ) {
								h_unsat++;
								countable_nonheavy_unsats_++;
								atom_bur_unsat_.set( hatm_id, true );
								// legacy behavior assumes H's will not be added, so should never get here in legacy behavior (same for BUNS2 filter)
							}
						}

						if ( h_unsat == h_count ) { // every Hpol attached to N donor is unsat, so heavy-atom N is unsat
							this_atom_bur_unsat = true;
							all_heavy_unsats_++;
							residue_bur_unsat_polars_[ resnum ]++;
							if ( this_pose.residue( resnum ).atom_is_backbone( at ) ) {
								bb_heavy_unsats_++;
							}
						}
					} else if ( this_pose.residue( resnum ).atom_type( at ).is_acceptor() && this_pose.residue( resnum ).atomic_charge( at ) != 0.0 ) {
						// != 0.0 important for ligand / small molecule case (if atom has net charge 0.0 indicates it should be ignored
						if ( atom_hbonds.value()[ atid ] == 0 ) { // if "OH" this will not be 0
							this_atom_bur_unsat = true;
							all_heavy_unsats_++;
							residue_bur_unsat_polars_[ resnum ]++;
							if ( this_pose.residue( resnum ).atom_is_backbone( at ) ) {
								bb_heavy_unsats_++;
							}
						}
					}
				}
			} // if donor or acceptor
			atom_bur_unsat_.set( atid, this_atom_bur_unsat );
		} // for all heavy atoms
	}
} //recompute

void
BuriedUnsatisfiedPolarsCalculator::calculate_sasa(
	core::pose::Pose const & pose,
	core::id::AtomID_Map< core::Real > & atom_sasa,
	utility::vector1< core::Real > & residue_sasa
) const {

	atom_sasa.clear();
	residue_sasa.clear();

	if ( name_of_sasa_calc_ == "dalphaball" ) {

		core::id::AtomID_Mask atoms;
		core::pose::initialize_atomid_map( atoms, pose, true );
		core::scoring::packing::SurfVol surf_vol = core::scoring::packing::get_surf_vol( pose, atoms, probe_radius_ );

		atom_sasa = surf_vol.surf;
		residue_sasa.resize( atom_sasa.size() );
		for ( core::Size seqpos = 1; seqpos <= atom_sasa.size(); seqpos++ ) {
			for ( core::Size atno = 1; atno < atom_sasa.n_atom( seqpos ); atno++ ) {
				residue_sasa[ seqpos ] += atom_sasa( seqpos, atno );
			}
		}

	} else if ( name_of_sasa_calc_ == "nondefault" ) { // if we need to update pore_radius, don't use default calculator with default radius
		core::scoring::calc_per_atom_sasa( pose, atom_sasa, residue_sasa, probe_radius_);
	} else {
		basic::MetricValue< core::id::AtomID_Map< core::Real > > temp_atom_sasa;
		basic::MetricValue< utility::vector1< core::Real > > temp_residue_sasa;
		pose.metric( name_of_sasa_calc_, "atom_sasa", temp_atom_sasa);
		pose.metric( name_of_sasa_calc_, "residue_sasa", temp_residue_sasa);
		atom_sasa = temp_atom_sasa.value();
		residue_sasa = temp_residue_sasa.value();
	}

}



///@brief ONLY USED FOR LEGACY BEHAVIOR (legacy=true)
core::Size
BuriedUnsatisfiedPolarsCalculator::satisfaction_cutoff( std::string atom_type )
{
	//according to jk, buried hydroxyls are often seen making only one hydrogen bond. also, ether oxygens often are bad h-bond acceptors
	if ( atom_type == "OH" ) return 2;

	//backbone oxygens also only have one h-bbond in most secondary structure elements
	else if ( atom_type == "OCbb" ) return 2;

	else if ( atom_type ==  "S" ) return 2;

	// no longer the assumption, except in legacy=true
	//everything else we expect to have 3 bonded/h-bonded neighbours to count as satisfied
	else return 3;
}

} //namespace simple_pose_metric_calculators
} //namespace protocols

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
protocols::simple_pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator::BuriedUnsatisfiedPolarsCalculator() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::simple_pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::EnergyDependentCalculator >( this ) );
	arc( CEREAL_NVP( all_bur_unsat_polars_ ) ); // core::Size
	arc( CEREAL_NVP( bb_heavy_unsats_ ) ); // core::Size
	arc( CEREAL_NVP( all_heavy_unsats_ ) ); // core::Size
	arc( CEREAL_NVP( countable_nonheavy_unsats_ ) ); // core::Size
	arc( CEREAL_NVP( atom_bur_unsat_ ) ); // core::id::AtomID_Map<_Bool>
	arc( CEREAL_NVP( residue_bur_unsat_polars_ ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( name_of_hbond_calc_ ) ); // std::string
	arc( CEREAL_NVP( name_of_sasa_calc_ ) ); // std::string
	arc( CEREAL_NVP( burial_cutoff_ ) ); // core::Real
	arc( CEREAL_NVP( probe_radius_ ) ); // core::Real
	arc( CEREAL_NVP( residue_surface_cutoff_ ) ); // core::Real
	arc( CEREAL_NVP( special_region_ ) ); // core::id::AtomID_Map<_Bool>
	arc( CEREAL_NVP( special_region_entire_residue_ ) ); // bool
	arc( CEREAL_NVP( generous_hbonds_ ) ); // bool
	arc( CEREAL_NVP( legacy_counting_ ) ); // bool
	arc( CEREAL_NVP( vsasa_ ) ); // bool
	arc( CEREAL_NVP( use_sc_neighbors_ ) ); // bool
	arc( CEREAL_NVP( skip_surface_res_ ) ); // bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::simple_pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::EnergyDependentCalculator >( this ) );
	arc( all_bur_unsat_polars_ ); // core::Size
	arc( bb_heavy_unsats_ ); // core::Size
	arc( all_heavy_unsats_ ); // core::Size
	arc( countable_nonheavy_unsats_ ); // core::Size
	arc( atom_bur_unsat_ ); // core::id::AtomID_Map<_Bool>
	arc( residue_bur_unsat_polars_ ); // utility::vector1<core::Size>
	arc( name_of_hbond_calc_ ); // std::string
	arc( name_of_sasa_calc_ ); // std::string
	arc( burial_cutoff_ ); // core::Real
	arc( probe_radius_ ); // core::Real
	arc( residue_surface_cutoff_ ); // core::Real
	arc( special_region_ ); // core::id::AtomID_Map<_Bool>
	arc( special_region_entire_residue_ ); // bool
	arc( generous_hbonds_ ); // bool
	arc( legacy_counting_ ); // bool
	arc( vsasa_ ); // bool
	arc( use_sc_neighbors_ ); // bool
	arc( skip_surface_res_ ); // bool
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::simple_pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator );
CEREAL_REGISTER_TYPE( protocols::simple_pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_simple_pose_metric_calculators_BuriedUnsatisfiedPolarsCalculator )
#endif // SERIALIZATION
