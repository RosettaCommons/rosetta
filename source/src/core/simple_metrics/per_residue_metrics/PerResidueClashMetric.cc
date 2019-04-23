// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/PerResidueClashMetric.cc
/// @brief A SimpleMetric that calculates the number of atomic clashes per residue using the LJ radius (at 0).  Can use a soft radius, which reduces it by 33%.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/per_residue_metrics/PerResidueClashMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>

// Core headers
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.simple_metrics.per_residue_metrics.PerResidueClashMetric" );


namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

using namespace core::select;
using namespace core::select::residue_selector;
using namespace core::scoring;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
PerResidueClashMetric::PerResidueClashMetric():
	core::simple_metrics::PerResidueRealMetric()
{}

PerResidueClashMetric::PerResidueClashMetric( ResidueSelectorCOP selector1, ResidueSelectorCOP selector2):
	core::simple_metrics::PerResidueRealMetric()
{
	set_residue_selector(selector1);
	set_secondary_residue_selector(selector2);
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PerResidueClashMetric::~PerResidueClashMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
PerResidueClashMetric::PerResidueClashMetric( PerResidueClashMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
PerResidueClashMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new PerResidueClashMetric( *this ) );

}

std::string
PerResidueClashMetric::name() const {
	return name_static();
}

std::string
PerResidueClashMetric::name_static() {
	return "PerResidueClashMetric";

}
std::string
PerResidueClashMetric::metric() const {
	return "atomic_clashes";
}

///@brief Should we calculate only heavy-heavy atom clashes
/// Default True
void
PerResidueClashMetric::set_use_hydrogens( bool use_hydrogens ){
	use_hydrogens_ = use_hydrogens;
}

///@brief When we calculate atom-atom distances using LJ distances,
///
///@details
///  clash if distance < (atomI_LJ + atomJ_LJ)*(1 - soft_clash)
///
void
PerResidueClashMetric::set_use_soft_clash( bool soft_clash_check ){
	use_soft_clash_ = soft_clash_check;
}

///@brief Set the dampening of the LJ to use for soft-clash.
/// Default=.33
void
PerResidueClashMetric::set_soft_dampening( core::Real dampening ){
	soft_clash_dampening_ = dampening;
}

void
PerResidueClashMetric::set_secondary_residue_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	selector2_ = selector;
}

core::Real
PerResidueClashMetric::lj_radius_to_zero_e_radius( core::Real num) const{
	return num/lj_n_;
}

void
PerResidueClashMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{

	SimpleMetric::parse_base_tag( tag );
	PerResidueRealMetric::parse_per_residue_tag( tag, datamap );


	set_use_hydrogens( tag->getOption< bool >("use_hydrogens", use_hydrogens_));
	set_use_soft_clash( tag->getOption<bool >("soft_clash", use_soft_clash_));
	set_soft_dampening( tag->getOption< core::Real >("dampening_percent", soft_clash_dampening_ ));

	if ( soft_clash_dampening_ > .99 ) {
		utility_exit_with_message("Dampening percent must be a value between 0 and .99");
	}

	if ( tag->hasOption("residue_selector2") ) {
		set_secondary_residue_selector(select::residue_selector::parse_residue_selector( tag, datamap, "residue_selector2" ));
	}
}

void
PerResidueClashMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("use_hydrogens", xsct_rosetta_bool, "Use Hydrogens OR calculate only heavy atom clashes", "false");
	attlist + XMLSchemaAttribute::attribute_w_default("soft_clash", xsct_rosetta_bool, "IF TRUE, soften the clashing by a percent (dampening_percent).  When we calculate atom-atom distances using LJ distances, clash if distance less_than (atomI_LJ + atomJ_LJ)*(1 - dampening_percent)", "true");

	attlist + XMLSchemaAttribute::attribute_w_default("dampening_percent", xsct_real, "If we use soft_clash, dampen by this percent.  Number is between 0 and .99. IE- When we calculate atom-atom distances using LJ distances, clash if distance less_then (atomI_LJ + atomJ_LJ)*(1 - dampening_percent)", ".33");

	attributes_for_parse_residue_selector( attlist, "residue_selector2", "The residue selector to use for calculating clashes TO for each residue of the primary selector.  Default is to calculate clashes of each residue to ALL OTHER residues (including those in the primary selector)" );

	std::string description =
		"Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)\n"
		"A SimpleMetric that calculates the total number of atom-atom clashes from a residue in a residue selector\n"
		" to all other residues defined in a second residue selector using the LJ radius of each atom \n"
		"Can use a soft radius, which reduces it by 33% (by default - see dampening_percent option. \n"
		"DETAILS: \n"
		"  Does NOT calculate INTRA-RESIDUE clashes!\n";


	core::simple_metrics::xsd_per_residue_real_metric_type_definition_w_attributes(xsd, name_static(),
		description, attlist);
}

std::map< core::Size, core::Real >
PerResidueClashMetric::calculate(const pose::Pose & input_pose) const {
	if ( soft_clash_dampening_ > .99 || soft_clash_dampening_ <= 0 ) {
		utility_exit_with_message("Dampening percent must be a value between 0 and .99");
	}

	// The pose must be scored for this metric to work. Here, we use a modified version of what Brian Conventry added to
	//  the Neighborhood Residue Selector.
	bool using_clone = ! input_pose.energies().energies_updated();

	//If fa_rep has a weight 0, we need to score it.
	if ( ! using_clone && input_pose.energies().weights()[fa_rep] == 0 ) {
		using_clone = true;
	}
	pose::Pose pose_clone;
	if ( using_clone ) {
		ScoreFunctionOP scorefxn = ScoreFunctionOP( new ScoreFunction());
		scorefxn->set_weight(fa_rep, 1.0 );
		pose_clone = input_pose;
		scorefxn->score(pose_clone);
		if ( TR.Warning.visible() ) {
			TR.Warning << "################ Cloning pose and Scoring! ##############################" << std::endl;
			TR.Warning << "Ensure that pose is scored " << std::endl;
			TR.Warning << "before PerResidueClashMetric for maximum performance!" << std::endl;
			TR.Warning << "##########################################################################" << std::endl;
		}
	}
	const pose::Pose & pose = using_clone ? pose_clone : input_pose;

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	std::map< Size, Real > residue_atomic_clashes;

	utility::vector1< Size > selection1 = get_residues_from_subset( get_selector()->apply(pose) );
	utility::vector1< Size > selection2;

	if ( selector2_ == nullptr ) {
		TrueResidueSelector true_selector = TrueResidueSelector();
		selection2 = get_residues_from_subset( true_selector.apply(pose));
	} else {
		selection2 = get_residues_from_subset( selector2_->apply( pose ));
	}

	for ( core::Size resA : selection1 ) {
		core::Size clashes = 0;
		if ( pose.residue_type(resA).is_virtual_residue() ) continue;
		for ( core::Size resB : selection2 ) {
			if ( resA == resB ) continue;
			if ( pose.residue_type(resB).is_virtual_residue() ) continue;

			//Check energy before going through all atoms.
			EnergyEdge const * edge = energy_graph.find_energy_edge( resA, resB );
			if ( edge == nullptr ) continue;
			core::Real fa_rep_score = edge->fill_energy_map()[fa_rep];
			if ( fa_rep_score == 0 ) continue;

			///Atoms
			for ( core::Size atomA = 1; atomA <= pose.residue( resA ).natoms(); ++atomA ) {
				if ( pose.residue( resA ).atom_type( atomA).is_virtual() ) continue;
				if ( (! use_hydrogens_) && pose.residue( resA).atom_type( atomA ).element() == "H" ) continue;

				for ( core::Size atomB = 1; atomB <= pose.residue( resB ).natoms(); ++atomB ) {
					if ( pose.residue( resB ).atom_type( atomB).is_virtual() ) continue;
					if ( (! use_hydrogens_) && pose.residue( resB).atom_type( atomB ).element() == "H" ) continue;

					//If these atoms are bonded, we skip them.
					id::AtomID atm_id1 = id::AtomID( atomA, resA);
					id::AtomID atm_id2 = id::AtomID( atomB, resB);
					if ( pose.conformation().is_bonded(atm_id1, atm_id2) ) {
						TR <<"Bonded: "<<atm_id1 << ":" << atm_id2 << std::endl;
						continue;
					}

					if ( is_clashing(pose, resA, atomA, resB, atomB) ) {
						clashes+=1;
					}
				}
			} //Atoms
		}
		residue_atomic_clashes[resA] = clashes;
	} //Residues

	return residue_atomic_clashes;
}

bool
PerResidueClashMetric::is_clashing( core::pose::Pose const & pose, core::Size resA, core::Size atomA, core::Size resB, core::Size atomB) const {

	numeric::xyzVector<core::Real> gx_xyz= pose.residue( resA ).xyz( atomA );
	numeric::xyzVector<core::Real> ox_xyz= pose.residue( resB ).xyz( atomB );

	core::Real gx_ox_dis = gx_xyz.distance_squared(ox_xyz);

	core::Real atm1_lj = pose.residue( resA ).atom_type( atomA ).lj_radius();
	core::Real atm2_lj = pose.residue( resB ).atom_type( atomB ).lj_radius();

	core::Real atm1_vdw_radii = lj_radius_to_zero_e_radius( atm1_lj  );

	core::Real atm2_vdw_radii = lj_radius_to_zero_e_radius( atm2_lj  );

	core::Real cutoff =  atm1_vdw_radii + atm2_vdw_radii;
	if ( use_soft_clash_ ) {
		cutoff  =cutoff*(1 - soft_clash_dampening_ );
	}

	//std::cout << gx_ox_dis << " " << normal_cutoff << " " << soft_cutoff << std::endl;
	if ( gx_ox_dis < pow(cutoff, 2) ) {
		TR.Debug <<std::endl << std::endl << "Clashing: " <<pose.pdb_info()->pose2pdb(resA) <<":"<<pose.pdb_info()->pose2pdb(resB) << std::endl;
		TR.Debug << "  atmA:    " << pose.residue( resA).atom_name( atomA ) ;
		TR.Debug << "  atmB:    " << pose.residue( resB).atom_name( atomB )  << std::endl;
		TR.Debug << std::endl;
		TR.Debug <<"atm1 LJ: "<<atm1_lj <<" atm2 LJ " << atm2_lj << std::endl;
		TR.Debug <<"atm1 ra: "<<atm1_vdw_radii<<" atm2 ra: "<<atm2_vdw_radii <<std::endl;
		TR.Debug <<gx_ox_dis<< "<" << cutoff << std::endl;
		return true;
	}
	return false;
}

void
PerResidueClashMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PerResidueClashMetric::provide_xml_schema( xsd );
}

std::string
PerResidueClashMetricCreator::keyname() const {
	return PerResidueClashMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
PerResidueClashMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new PerResidueClashMetric );

}

} //core
} //simple_metrics
} //per_residue_metrics


#ifdef    SERIALIZATION



template< class Archive >
void
core::simple_metrics::per_residue_metrics::PerResidueClashMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric>( this ) );
	arc( CEREAL_NVP( use_soft_clash_ ) );
	arc( CEREAL_NVP( use_hydrogens_ ) );
	arc( CEREAL_NVP( soft_clash_dampening_ ) );
	arc( CEREAL_NVP( selector2_ ));
	arc( CEREAL_NVP( lj_n_));
}

template< class Archive >
void
core::simple_metrics::per_residue_metrics::PerResidueClashMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric >( this ) );
	arc( use_soft_clash_ );
	arc( use_hydrogens_ );
	arc( soft_clash_dampening_ );

	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_selector;
	arc( local_selector ); // ResidueSelectorCOP
	selector2_ = local_selector; // copy the non-const pointer(s) into the const pointer(s)
	arc( lj_n_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::per_residue_metrics::PerResidueClashMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::per_residue_metrics::PerResidueClashMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_per_residue_metrics_PerResidueClashMetric )
#endif // SERIALIZATION




