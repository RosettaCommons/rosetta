// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/CrankshaftFlipMover.cc
/// @brief This mover performs a cyclic-aware crankshaft flip at a residue position.
/// @author P. Douglas Renfrew (doug.renfrew@gmail.com)
/// @author James Eastwood (jre318@nyu.edu)

// Unit headers
#include <protocols/cyclic_peptide/CrankshaftFlipMover.hh>
#include <protocols/cyclic_peptide/CrankshaftFlipMoverCreator.hh>

// Protocol headers
#include <protocols/rosetta_scripts/util.hh>

// Core headers
#include <core/pose/Pose.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>

#include <core/conformation/Residue.hh>

#include <core/select/residue_selector/ResidueSelector.hh>

// Basic/Utility/Numeric headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

#include <numeric/random/random_permutation.hh>
#include <numeric/angle.functions.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Citation Manager
#include <basic/citation_manager/CitationManager.hh>
#include <basic/citation_manager/CitationCollection.hh>

#include <cmath>
#include <numeric>

static basic::Tracer TR( "protocols.cyclic_peptide.CrankshaftFlipMover" );

namespace protocols {
namespace cyclic_peptide {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
CrankshaftFlipMover::CrankshaftFlipMover():
	protocols::moves::Mover( CrankshaftFlipMover::mover_name() ),
	selector_()
{
	allow_peptoid_flips_ = true;
	allow_peptide_flips_ = true;
	allow_pre_gly_flips_ = true;
	tolerance_ = 20.0;
}

/// @brief Constructor
CrankshaftFlipMover::CrankshaftFlipMover(
	core::select::residue_selector::ResidueSelectorCOP selector_in,
	bool allow_peptoid_flips_in = true,
	bool allow_peptide_flips_in = true,
	bool allow_pre_gly_flips_in = true,
	core::Real tolerance_in = 20.0
) :
	allow_peptoid_flips_( allow_peptoid_flips_in ),
	allow_peptide_flips_( allow_peptide_flips_in ),
	allow_pre_gly_flips_( allow_pre_gly_flips_in ),
	tolerance_( tolerance_in ),
	selector_() // cloned below
{
	if ( selector_in ) {
		selector_ = selector_in->clone();
	}
}


/// @brief Copy constructor
CrankshaftFlipMover::CrankshaftFlipMover( CrankshaftFlipMover const & src ):
	Mover( src ),
	allow_peptoid_flips_( src.allow_peptoid_flips_ ),
	allow_peptide_flips_( src.allow_peptide_flips_ ),
	allow_pre_gly_flips_( src.allow_pre_gly_flips_ ),
	tolerance_( src.tolerance_ ),
	selector_() // cloned below
{
	if ( src.selector_ ) {
		selector_ = src.selector_->clone();
	}
}


////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
CrankshaftFlipMover::~CrankshaftFlipMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief copy assignment operator
CrankshaftFlipMover &
CrankshaftFlipMover::operator=( CrankshaftFlipMover const & src ) {

	if ( this != &src ) {
		allow_peptoid_flips_ = src.allow_peptoid_flips_;
		allow_peptide_flips_ = src.allow_peptide_flips_;
		allow_pre_gly_flips_ = src.allow_pre_gly_flips_;
		tolerance_ = src.tolerance_;

		if ( src.selector_ ) {
			selector_ = src.selector_->clone();
		}
	}

	return *this;
}


/// @brief Apply the mover
void
CrankshaftFlipMover::apply( core::pose::Pose& pose ){

	// check selector, if null then all positions
	core::select::residue_selector::ResidueSubset selected_positions;
	if ( selector_ ) {
		selected_positions = selector_->apply( pose );
	} else {
		selected_positions = core::select::residue_selector::ResidueSubset( pose.size(), true );
	}


	// make a vector of selected residues indicies in random order
	utility::vector1< core::Size > random_positions;
	for ( core::Size i(1); i <= selected_positions.size(); ++i ) {
		if ( selected_positions[ i ] == true ) {
			random_positions.push_back(i);
		}
	}
	numeric::random::random_permutation( random_positions );

	// go through positions in random order, break afer one sucessful flip
	for ( core::Size i(1); i <= random_positions.size(); ++i ) {

		// get indicies of res1 and it's upper conectted residue, res2
		core::Size res1_index, res2_index;
		res1_index = random_positions[i];
		if ( pose.residue( res1_index ).has_upper_connect() ) {
			res2_index = pose.residue( res1_index ).connected_residue_at_resconn( pose.residue( res1_index ).type().upper_connect_id() );
			if ( res2_index == 0 ) {
				TR.Debug << "Position " << res1_index << " has an upper connect at position but it is zero, skipping" << std::endl;
				continue;
			} else {
				TR.Debug << "Position " << res1_index << " has an upper connect at position " << res2_index << std::endl;
			}
		} else {
			TR.Debug << "Position " << res1_index << " does not have an upper connect, skipping" << std::endl;
			continue;
		}

		if ( allow_peptoid_flips_ && try_peptoid_flip( pose, res1_index, res2_index ) ) {
			TR.Debug << "Performed peptoid flip at positions " << res1_index << " and " << res2_index << std::endl;
			break;
		} else if ( allow_peptide_flips_ && try_peptide_flip( pose, res1_index, res2_index ) ) {
			TR.Debug << "Performed peptide flip at positions " << res1_index << " and " << res2_index << std::endl;
			break;
		} else if ( allow_pre_gly_flips_ && try_pre_gly_flip( pose, res1_index, res2_index ) ) {
			TR.Debug << "Performed pre-gly flip at positions " << res1_index << " and " << res2_index << std::endl;
			break;
		}
	}
}

/// @brief returns true if the angle is close to the target withing some tolerance
bool
CrankshaftFlipMover::angle_degrees_in_range( core::Real angle, core::Real target ) {
	if ( std::abs( numeric::principal_angle_degrees( angle - target ) ) <= tolerance_ ) {
		return true;
	} else {
		return false;
	}
}


/// @brief returns true if res1 is Peptoid or Gly and res2 is Peptoid, Pro, or N-Methyl
bool
CrankshaftFlipMover::try_peptoid_flip( core::pose::Pose & pose, core::Size res1_index, core::Size res2_index ) {
	// residue types
	core::chemical::ResidueType const & res1_type( pose.residue( res1_index ).type() );
	core::chemical::ResidueType const & res2_type( pose.residue( res2_index ).type() );

	if ( // check residue identities
			( res1_type.is_peptoid() || res1_type.aa() == core::chemical::AA::aa_gly ) &&
			( res2_type.is_peptoid() || res2_type.aa() == core::chemical::AA::aa_pro || res2_type.is_N_substituted() )
			) {

		// flip residues
		pose.set_phi( res1_index, ( numeric::principal_angle_degrees( pose.phi( res1_index ) ) * -1.0 ) );
		pose.set_psi( res1_index, ( numeric::principal_angle_degrees( pose.psi( res1_index ) ) * -1.0 ) );
		pose.set_omega( res1_index, ( numeric::principal_angle_degrees( pose.omega( res1_index ) ) + 180.0 ) );

		return true;
	} else {
		return false;
	}
}


/// @brief returns true if res1 is Peptide and res2 is Peptide (not P) and both positions have correct dihedrals
bool
CrankshaftFlipMover::try_peptide_flip( core::pose::Pose & pose, core::Size res1_index, core::Size res2_index ) {
	// residue types
	core::chemical::ResidueType const & res1_type( pose.residue( res1_index ).type() );
	core::chemical::ResidueType const & res2_type( pose.residue( res2_index ).type() );

	// residue dihedrals
	core::Real res1_psi( numeric::principal_angle_degrees( pose.psi( res1_index ) ) );
	core::Real res2_phi( numeric::principal_angle_degrees( pose.phi( res2_index ) ) );
	core::Real res1_omg( numeric::principal_angle_degrees( pose.omega( res1_index ) ) );

	TR.Debug << "Peptide Flip res1_name / res2_name / res1_psi / res2_phi / res1_omg: " << res1_type.name() << " / " << res2_type.name() << " / " << res1_psi << " / " << res2_phi << " / " << res1_omg << std::endl;

	if (
			( // check residue identities
			( res1_type.is_protein() ) &&
			( res2_type.is_protein() && res2_type.aa() != core::chemical::AA::aa_pro )
			) &&
			( // check dihedral angles
			( ( angle_degrees_in_range( res1_psi,  120.0 ) ) && ( angle_degrees_in_range( res2_phi,  60.0 ) ) ) ||
			( ( angle_degrees_in_range( res1_psi, -120.0 ) ) && ( angle_degrees_in_range( res2_phi, -60.0 ) ) )
			)
			) {

		// flip residues
		pose.set_psi( res1_index, res1_psi * -1.0 );
		pose.set_phi( res2_index, res2_phi * -1.0 );
		pose.set_omega( res1_index, res1_omg * -1.0 );

		return true;

	} else {
		return false;
	}
}


/// @brief returns true if res1 is Peptoid or Peptide and res2 is Gly and both positions have correct dihedrals
bool
CrankshaftFlipMover::try_pre_gly_flip( core::pose::Pose & pose, core::Size res1_index, core::Size res2_index ) {
	// residue types
	core::chemical::ResidueType const & res1_type( pose.residue( res1_index ).type() );
	core::chemical::ResidueType const & res2_type( pose.residue( res2_index ).type() );

	// residue dihedrals
	core::Real res1_psi( numeric::principal_angle_degrees( pose.psi( res1_index ) ) );
	core::Real res2_phi( numeric::principal_angle_degrees( pose.phi( res2_index ) ) );

	TR.Debug << "Pre-Gly Flip res1_name / res2_name / res1_psi / res2_phi : " << res1_type.name() << " / " << res2_type.name() << " / " << res1_psi << " / " << res2_phi << std::endl;

	if ( // check residue identities
			( res1_type.is_peptoid() || res1_type.is_protein() ) &&
			( res2_type.aa() == core::chemical::AA::aa_gly )
			) {
		// check dihedral angles for different cases
		if ( // near -30,-90
				angle_degrees_in_range( res1_psi, -30.0 ) && angle_degrees_in_range( res2_phi, -90.0 )
				) {
			// flip to near 0, 180
			pose.set_psi( res1_index, numeric::principal_angle_degrees( res1_psi + 30.0 ) );
			pose.set_phi( res2_index, numeric::principal_angle_degrees( res2_phi - 90.0 ) );
			return true;

		} else if ( // near 180, 60
				angle_degrees_in_range( res1_psi, 180.0 ) && angle_degrees_in_range( res2_phi,  60.0 )
				) {
			// flip to near 0, 180
			pose.set_psi( res1_index, numeric::principal_angle_degrees( res1_psi - 180.0 ) );
			pose.set_phi( res2_index, numeric::principal_angle_degrees( res2_phi - 120.0 ) );
			return true;

		} else if ( // near 180, -60
				angle_degrees_in_range( res1_psi, 180.0 ) && angle_degrees_in_range( res2_phi, -60.0 )
				) {
			// flip to near 0, 180
			pose.set_psi( res1_index, numeric::principal_angle_degrees( res1_psi - 180.0 ) );
			pose.set_phi( res2_index, numeric::principal_angle_degrees( res2_phi + 120.0 ) );
			return true;

		} else if ( // near 0, 180
				angle_degrees_in_range( res1_psi,   0.0 ) && angle_degrees_in_range( res2_phi, 180.0 )
				) {
			// flip to near -30, -90 (case 1) or 180, +60 (case 2) or 180, -60 (case 3) randomly
			core::Size bin( numeric::random::random_range( 1, 3 ) );
			switch( bin ) {
			case 1 :
				pose.set_psi( res1_index, numeric::principal_angle_degrees( res1_psi - 30.0 ) );
				pose.set_phi( res2_index, numeric::principal_angle_degrees( res2_phi + 90.0 ) );
				break;
			case 2 :
				pose.set_psi( res1_index, numeric::principal_angle_degrees( res1_psi + 180.0 ) );
				pose.set_phi( res2_index, numeric::principal_angle_degrees( res2_phi + 120.0 ) );
				break;
			case 3 :
				pose.set_psi( res1_index, numeric::principal_angle_degrees( res1_psi + 180.0 ) );
				pose.set_phi( res2_index, numeric::principal_angle_degrees( res2_phi - 120.0 ) );
				break;
			}
			return true;
		} else { // dihedral angle check
			return false;
		}
	} else { // residue identity check
		return false;
	}
}


/// @brief Setter for selector
void
CrankshaftFlipMover::selector( core::select::residue_selector::ResidueSelectorCOP selector ){
	selector_ = selector->clone();
}


////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
CrankshaftFlipMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
CrankshaftFlipMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) {

	if ( tag->hasOption( "residue_selector" ) ) {
		core::select::residue_selector::ResidueSelectorCOP res_selector = protocols::rosetta_scripts::parse_residue_selector( tag, data );
		if ( !res_selector ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "ResidueSelector passed to CrankshaftFlipMover mover could not be found." );
		}
		selector_ = res_selector;
	}

	allow_peptoid_flips_ = tag->getOption<bool>( "allow_peptoid_flips", allow_peptoid_flips_ );
	allow_peptide_flips_ = tag->getOption<bool>( "allow_peptide_flips", allow_peptide_flips_ );
	allow_pre_gly_flips_ = tag->getOption<bool>( "allow_pre_gly_flips", allow_pre_gly_flips_ );
	tolerance_ = tag->getOption< core::Real >( "tolerance", tolerance_ );

}


void
CrankshaftFlipMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	// bools
	attlist
		+ XMLSchemaAttribute::attribute_w_default("allow_peptoid_flips", xsct_rosetta_bool, "Allow peptoid type crankshaft flips", "true" )
		+ XMLSchemaAttribute::attribute_w_default("allow_peptide_flips", xsct_rosetta_bool, "Allow peptide type crankshaft flips", "true" )
		+ XMLSchemaAttribute::attribute_w_default("allow_pre_gly_flips", xsct_rosetta_bool, "Allow pre-gly type crankshaft flips", "true" )
		+ XMLSchemaAttribute::attribute_w_default("tolerance", xsct_real, "Half-width of allowed range for peptoid and pre-gly flips", "20.0" );

	// residue selector
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector", "Previously-defined ResidueSelector, specifying the subset of residues to which the mover will be applied." );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "This mover performs a cyclic-aware crankshaft flip at a residue position.", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
CrankshaftFlipMover::fresh_instance() const
{
	return utility::pointer::make_shared< CrankshaftFlipMover >();
}


/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
CrankshaftFlipMover::clone() const
{
	return utility::pointer::make_shared< CrankshaftFlipMover >( *this );
}


std::string CrankshaftFlipMover::get_name() const {
	return mover_name();
}


std::string CrankshaftFlipMover::mover_name() {
	return "CrankshaftFlipMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
CrankshaftFlipMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< CrankshaftFlipMover >();
}

std::string
CrankshaftFlipMoverCreator::keyname() const
{
	return CrankshaftFlipMover::mover_name();
}

void CrankshaftFlipMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CrankshaftFlipMover::provide_xml_schema( xsd );
}

/// @brief This mover is unpublished.  It returns P. Douglas Renfrew as its author.
void
CrankshaftFlipMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	basic::citation_manager::CitationCollectionOP collection(
		utility::pointer::make_shared< basic::citation_manager::CitationCollection >(
		mover_name(),
		basic::citation_manager::CitedModuleType::Mover
		)
	);
	collection->add_citation( basic::citation_manager::CitationManager::get_instance()->get_citation_by_doi("10.1021/acs.jpcb.2c01669") );

	citations.add( collection );
}


////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, CrankshaftFlipMover const & mover )
{
	mover.show(os);
	return os;
}


} //cyclic_peptide
} //protocols
