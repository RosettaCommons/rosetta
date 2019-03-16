// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/InterfaceHbondsFilter.cc
/// @brief Counts the number of hydrogen bonds across an interface defined by a jump, the salt
/// @brief bridge mode only counts the salt bridges which are defined to be formed if the distance
/// @brief between any of the oxygen atoms of acidic residues and the nitrogen atoms of basic
/// @brief residues are within the cut-off distance (default 4.0).
/// @author Longxing Cao (longxing@uw.edu)

// Unit headers
#include <protocols/simple_filters/InterfaceHbondsFilter.hh>
#include <protocols/simple_filters/InterfaceHbondsFilterCreator.hh>


// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

#include <protocols/filters/filter_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

//hbonds headers
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

// interface detection
#include <protocols/scoring/Interface.hh>


// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyzVector.hh>


namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.InterfaceHbondsFilter" );


InterfaceHbondsFilter::InterfaceHbondsFilter() :
	Filter("InterfaceHbonds"),
	hbonds_num_threshold_(0),
	hbonds_energy_cutoff_(0.0),
	jump_(1),
	scorefxn_(nullptr),
	salt_bridge_mode_(false),
	include_His_chain_terminus_(false),
	salt_bridge_distance_cutoff_(4.0)
{
	if ( ! scorefxn_ ) {
		scorefxn_ = core::scoring::get_score_function( true /* full atom score function */ );
	}
}
InterfaceHbondsFilter::InterfaceHbondsFilter(
	core::Size hbonds_num_threshold_in,
	core::Real hbonds_energy_cutoff_in,
	core::Size jump_in,
	core::scoring::ScoreFunctionOP scorefxn_in,
	bool salt_bridge_mode_in,
	bool include_His_chain_terminus_in,
	core::Real salt_bridge_distance_cutoff_in
) :
	Filter("InterfaceHbonds"),
	hbonds_num_threshold_(hbonds_num_threshold_in),
	hbonds_energy_cutoff_(hbonds_energy_cutoff_in),
	jump_(jump_in),
	scorefxn_(scorefxn_in),
	salt_bridge_mode_(salt_bridge_mode_in),
	include_His_chain_terminus_(include_His_chain_terminus_in),
	salt_bridge_distance_cutoff_(salt_bridge_distance_cutoff_in),
	charged_res_()
{
	charged_res_.max_load_factor( 0.7 );
	if ( ! scorefxn_ ) {
		scorefxn_ = core::scoring::get_score_function( true /* full atom score function */ );
	}
	if ( salt_bridge_mode_ ) {
		initialize_charged_residue_map();
	}
}
InterfaceHbondsFilter::~InterfaceHbondsFilter() = default;

core::Size InterfaceHbondsFilter::compute_hbonds( core::pose::Pose const & pose ) const {
	core::scoring::hbonds::HBondSet temp_hbond_set;
	core::scoring::hbonds::HBondOptions new_options( temp_hbond_set.hbond_options() );
	new_options.decompose_bb_hb_into_pair_energies( true );
	//new_options.use_hb_env_dep( false );
	core::scoring::hbonds::HBondSetOP full_hbond_set( new core::scoring::hbonds::HBondSet( new_options ) );

	core::scoring::hbonds::fill_hbond_set( pose, false, *full_hbond_set, false, false, false, false );

	core::Size const chain1_num( pose.chain( pose.fold_tree().upstream_jump_residue( jump_ ) ) );
	core::Size const chain2_num( pose.chain( pose.fold_tree().downstream_jump_residue( jump_ ) ));
	core::Size total_hbonds(0);
	for ( core::Size ihb=1; ihb <= full_hbond_set->nhbonds(); ++ihb ) {
		core::scoring::hbonds::HBond const & hb( full_hbond_set->hbond(ihb) );
		if ( hb.energy() > hbonds_energy_cutoff_ ) {
			continue;
		}
		core::conformation::Residue const & res1( pose.residue(hb.don_res()) );
		core::conformation::Residue const & res2( pose.residue(hb.acc_res()) );
		if ( ( chain1_num == res1.chain()  && chain2_num == res2.chain() ) ||
				( chain1_num == res2.chain() && chain2_num == res1.chain() )
				) {
			TR << "New hbond found between the atom " << res1.atom_name( hb.don_hatm() )
				<< " of residue " << res1.seqpos()
				<< " and the atom  " << res2.atom_name( hb.acc_atm() )
				<< " of residue " << res2.seqpos()
				<< std::endl;
			++total_hbonds;
		}
	}
	return total_hbonds;
}



void InterfaceHbondsFilter::initialize_charged_residue_map()
{
	// initialize the map for Arg, Lys, Asp and Glu.
	charged_res_ = std::unordered_map<std::string, std::vector<Charged_Group> > {
		{ "ARG", { {"NE",POSITIVE}, {"NH1",POSITIVE}, {"NH2",POSITIVE} } },
		{ "LYS", { {"NZ",POSITIVE} } },
		{ "ASP", { {"OD1",NEGATIVE}, {"OD2",NEGATIVE} } },
		{ "GLU", { {"OE1",NEGATIVE}, {"OE2",NEGATIVE} } }
		};
	charged_res_.max_load_factor( 0.7 );


	if ( include_His_chain_terminus_ ) {
		// HIS
		charged_res_.insert( std::make_pair<std::string, std::vector<Charged_Group> >("HIS", {
			Charged_Group("ND1", POSITIVE),
			Charged_Group("NE2", POSITIVE)
			}));
		// HIS_D
		charged_res_.insert( std::make_pair<std::string, std::vector<Charged_Group> >("HIS_D", {
			Charged_Group("ND1", POSITIVE),
			Charged_Group("NE2", POSITIVE)
			}));
		std::vector<std::string> aas = {"ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS",
			"LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL",
			"TRP", "TYR" };
		for ( std::string aa : aas ) {
			std::vector<Charged_Group > nter = {
				Charged_Group("N", POSITIVE)
				};
			std::vector<Charged_Group > cter = {
				Charged_Group("O", NEGATIVE),
				Charged_Group("OXT", NEGATIVE)
				};
			std::unordered_map<std::string, std::vector<Charged_Group > >::const_iterator got = charged_res_.find (aa);
			if ( got != charged_res_.end() ) {
				nter.insert( nter.end(), got->second.begin(), got->second.end() );
				cter.insert( cter.end(), got->second.begin(), got->second.end() );
			}
			charged_res_.insert( std::make_pair<std::string, std::vector<Charged_Group > >(aa + ":NtermProteinFull", std::move(nter) ) );
			charged_res_.insert( std::make_pair<std::string, std::vector<Charged_Group > >(aa + ":CtermProteinFull", std::move(cter) ) );
		}
	}
}

core::Size InterfaceHbondsFilter::compute_salt_bridges( core::pose::Pose const & pose ) const {

	// After checking the code, the Interface class doesn't work as I thought, as it just simply
	// split the pose into two parts by the jump. How about a multi-chain pose? Keep a note here
	// just as a waring to be careful about multi-chain pdbs.
	protocols::scoring::Interface iface(jump_);
	iface.distance( 6.5 /* maximum length of Arg?*/ + 4.5 /* maximum length of Glu */ + salt_bridge_distance_cutoff_ );
	iface.calculate( pose );

	core::Real const squared_max_dist( salt_bridge_distance_cutoff_ * salt_bridge_distance_cutoff_ );
	core::Size total_salt_bridges(0);

	for ( core::Size ires = 1; ires <= pose.size(); ++ires ) {
		for ( core::Size jres = ires + 1; jres <= pose.size(); ++jres ) {
			core::conformation::Residue const & res1( pose.residue(ires) );
			core::conformation::Residue const & res2( pose.residue(jres) );
			if ( iface.is_pair( res1, res2 ) /*just table looking up, fast enough*/ ) {
				std::string res1_name  = res1.name();
				std::string res2_name  = res2.name();
				std::unordered_map<std::string, std::vector<Charged_Group > >::const_iterator res1_iter = charged_res_.find (res1_name);
				std::unordered_map<std::string, std::vector<Charged_Group > >::const_iterator res2_iter = charged_res_.find (res2_name);
				if ( res1_iter != charged_res_.end() && res2_iter != charged_res_.end() ) {
					bool flag = false;
					for ( Charged_Group atm1 : res1_iter->second ) {
						for ( Charged_Group atm2 : res2_iter->second ) {
							if ( atm1.second == atm2.second /* same charge */ ) {
								continue;
							}
							numeric::xyzVector< core::Real> atm1_xyz = res1.xyz( atm1.first );
							numeric::xyzVector< core::Real> atm2_xyz = res2.xyz( atm2.first );
							if ( (atm1_xyz - atm2_xyz).norm_squared() <= squared_max_dist ) {
								++total_salt_bridges;
								TR << "New salt bridge found between residue " << ires << " and residue " << jres << std::endl;
								flag = true;
								break;
							}
						}
						if ( flag ) break;
					}
				}
			}
		}
	}
	return total_salt_bridges;
}

core::Size InterfaceHbondsFilter::compute( core::pose::Pose const & pose ) const {
	core::pose::Pose work_pose( pose );
	// rescore the pose using the scoring function to make sure everything is updated,
	// as this is required for both hbonds and salt birdge detection( neighbor detection ).
	scorefxn_->score( work_pose );
	runtime_assert( work_pose.energies().energies_updated() );
	if ( ! salt_bridge_mode_ ) {
		return compute_hbonds( work_pose );
	} else {
		return compute_salt_bridges( work_pose );
	}
}

void
InterfaceHbondsFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const &
) {

	if ( tag->hasOption( "threshold" ) ) {
		hbonds_num_threshold_ = tag->getOption<core::Size>( "threshold" );
	}
	if ( tag->hasOption( "hbond_energy" ) ) {
		hbonds_energy_cutoff_ = tag->getOption<core::Real>( "hbond_energy_cutoff" );
	}
	if ( tag->hasOption( "jump" ) ) {
		jump_ = tag->getOption<core::Size>( "jump" );
	}
	if ( tag->hasOption( "salt_bridge_mode" ) ) {
		salt_bridge_mode_ = tag->getOption<bool>( "salt_bridge_mode" );
	}
	if ( tag->hasOption( "include_His_chain_terminus" ) ) {
		include_His_chain_terminus_ = tag->getOption<bool>( "include_His_chain_terminus" );
	}
	if ( tag->hasOption( "salt_bridge_distance" ) ) {
		salt_bridge_distance_cutoff_ = tag->getOption<core::Real>( "salt_bridge_distance" );
	}
	if ( tag->hasOption( "scorefxn" ) ) {
		set_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, "scorefxn", data ) );
	}

	if ( salt_bridge_mode_ ) {
		initialize_charged_residue_map();
	}

}

bool
InterfaceHbondsFilter::apply( core::pose::Pose const & pose ) const {
	core::Size const interface_hbonds( compute( pose ));
	return interface_hbonds >= hbonds_num_threshold_;
}

void
InterfaceHbondsFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Size const interface_hbonds( compute( pose ));
	out<<"The number of hydrogen bonds across jump "<< jump_ << " is "<< interface_hbonds << '\n';
}

core::Real
InterfaceHbondsFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Size const interface_hbonds( compute( pose ));
	return( (core::Real) interface_hbonds );
}

std::string InterfaceHbondsFilter::name() const {
	return class_name();
}

std::string InterfaceHbondsFilter::class_name() {
	return "InterfaceHbonds";
}

void InterfaceHbondsFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {

	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute(
		"name", xs_string,
		"Name of the filter");

	rosetta_scripts::attributes_for_parse_score_function_w_description(attlist, "scorefxn",
		"Which score function to use to find the hydrogen bonds?");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"threshold", xs_integer,
		"How many hydrogen bonds should there be across the interface?",
		"0");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"hbond_energy", xsct_real,
		"score cut of the hbond energy to ignore bad hydrogen bonds",
		"0.0");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"jump", xsct_non_negative_integer,
		"jump number that defines the interface",
		"1");


	attlist + XMLSchemaAttribute::attribute_w_default(
		"salt_bridge_mode", xsct_rosetta_bool,
		"Only count the number of salt bridges?",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"salt_bridge_distance", xsct_real,
		"cutoff distance for salt bridges",
		"4.0");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"include_His_chain_terminus", xsct_rosetta_bool,
		"treat histidine and chain terminus as charge residues?",
		"false");

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(),
		"Counts the number of hydrogen bonds across the interface that defines by the jump number", attlist );
}

std::string InterfaceHbondsFilterCreator::keyname() const {
	return InterfaceHbondsFilter::class_name();
}

protocols::filters::FilterOP
InterfaceHbondsFilterCreator::create_filter() const {
	return utility::pointer::make_shared< InterfaceHbondsFilter >();
}

void InterfaceHbondsFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InterfaceHbondsFilter::provide_xml_schema( xsd );
}


}
}
