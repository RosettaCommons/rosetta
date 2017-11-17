// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/bin_transitions/InitializeByBins.cc
/// @brief  This mover takes a stretch of backbone and initializes its mainchain torsions based on the probabilities of transitions from
/// one torsion bin to another.
/// @details Bin transitions are read from database files.  The algorithm is: set the first residue based on the probability of a residue
/// being in a bin.  Set subsequent residues based on the probability of a residue being in a bin given that the previous residue is in
/// a particular bin.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit headers
#include <protocols/simple_moves/bin_transitions/InitializeByBins.hh>
#include <protocols/simple_moves/bin_transitions/InitializeByBinsCreator.hh>

// Bin transition calculator headers:
#include <core/scoring/bin_transitions/BinTransitionCalculator.hh>
#include <core/scoring/bin_transitions/BinTransitionData.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
//parsing
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map
#include <protocols/rosetta_scripts/util.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>

#include <numeric/random/random.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace simple_moves {
namespace bin_transitions {

static basic::Tracer TR( "protocols.simple_moves.bin_transitions.InitializeByBins" );

// XRW TEMP std::string
// XRW TEMP InitializeByBinsCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return InitializeByBins::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP InitializeByBinsCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new InitializeByBins );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP InitializeByBins::mover_name()
// XRW TEMP {
// XRW TEMP  return "InitializeByBins";
// XRW TEMP }

/// @brief Default constructor
///
InitializeByBins::InitializeByBins() : //TODO: initialize variables here!
	protocols::moves::Mover("InitializeByBins"),
	start_res_(0),
	end_res_(0),
	binfile_("ABBA"),
	binfile_loaded_(false),
	bin_transition_calculator_( new core::scoring::bin_transitions::BinTransitionCalculator )
{}

/// @brief Copy constructor.
///
InitializeByBins::InitializeByBins( InitializeByBins const &src ) :
	protocols::moves::Mover( src ),
	start_res_(src.start_res_),
	end_res_(src.end_res_),
	binfile_(src.binfile_),
	binfile_loaded_(src.binfile_loaded_),
	bin_transition_calculator_( src.bin_transition_calculator_->clone() ) //CLONE this object.
{}

/// @brief Destructor.
///
InitializeByBins::~InitializeByBins() {}

// XRW TEMP std::string
// XRW TEMP InitializeByBins::get_name() const {
// XRW TEMP  return InitializeByBins::mover_name();
// XRW TEMP }

/// @brief Apply the mover to a pose.
///
void InitializeByBins::apply( core::pose::Pose & pose ) {
	//Check whether we've loaded bin transitions.
	if ( !bin_transition_calculator_->bin_params_loaded() ) {
		utility_exit_with_message(
			"In protocols::simple_moves::bin_transitions::InitializeByBins::apply(): Bin transition probability parameters must be loaded before calling the apply() function!");
	}

	//Number of residues in the pose
	core::Size const nres( pose.size() );
	if ( nres<1 ) utility_exit_with_message( "In protocols::simple_moves::bin_transitions::InitializeByBins::apply(): The pose has no residues!" );

	core::Size const startres( start_res_==0 ? 1 : start_res_ );
	core::Size const endres ( end_res_==0 ? nres : end_res_ );
	if ( endres < startres ) utility_exit_with_message( "In protocols::simple_moves::bin_transitions::InitializeByBins::apply(): The last residue cannot come before the first!" );
	if ( startres > nres ) utility_exit_with_message( "In protocols::simple_moves::bin_transitions::InitializeByBins::apply(): The first residue is out of range!  It must lie within [1, n_residues]." );
	if ( endres > nres ) utility_exit_with_message( "In protocols::simple_moves::bin_transitions::InitializeByBins::apply(): The last residue is out of range!  It must lie within [1, n_residues]." );

	utility::vector1 < utility::vector1 < core::Real > > mainchain_torsions; //Vector of mainchain torsion values that will be populated by the BinTransitionCalculator object.
	utility::vector1 < core::Size > res_indices; //Vector of residue indices to pass to the BinTransitionCalculator object
	res_indices.resize( endres - startres + 1, 0 );
	for ( core::Size i=1, imax=res_indices.size(); i<=imax; ++i ) res_indices[i] = startres + i - 1;

	//Actually generate the mainchain torsion values:
	bin_transition_calculator_->random_mainchain_torsions_from_bins( pose.conformation(), res_indices, mainchain_torsions) ; //pose.conformation() and loop_indices() are const inputs; mainchain_torsions is the output.

	//A thing that should be true at this point:
	debug_assert(mainchain_torsions.size() == endres-startres+1);

	{ //Scope to set the mainchain torsions:
		core::Size i=1; //mainchain_torsions index
		for ( core::Size ir=startres; ir<=endres; ++ir ) {
			debug_assert(mainchain_torsions[i].size()==pose.residue(ir).mainchain_torsions().size()); //Should be true at this point.
			for ( core::Size j=1, jmax=mainchain_torsions[i].size(); j<=jmax; ++j ) { //Set mainchain torsions
				pose.set_torsion( core::id::TorsionID( ir, core::id::BB, j ) , mainchain_torsions[i][j] );
			}

			++i; //Increment the mainchain_torsions index
		}
	}

	pose.update_residue_neighbors();

	return;
} //apply()

/// @brief Parse XML for RosettaScripts.
///
void InitializeByBins::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &//pose
)
{
	if ( tag->getName() != "InitializeByBins" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "This should be impossible");
	}

	//Get the bin params file:
	set_binfile_and_load( tag->getOption< std::string >( "bin_params_file", "ABBA" ) );

	//Get the residue ranges:
	set_residue_range( tag->getOption< core::Size >("start", 0), tag->getOption< core::Size >("end", 0) );

	if ( TR.visible() ) TR.flush();

	return;
} //parse_my_tag()

/// @brief Set the bin transition probability file.
///
void InitializeByBins::set_binfile_and_load( std::string const &name ) {
	if ( binfile_loaded_==true ) {
		utility_exit_with_message(
			"In protocols::simple_moves::bin_transitions::InitializeByBins::set_binfile_and_load(): The bin params file was already loaded!  This operation cannot be repeated.");
	}
	binfile_=name;
	if ( TR.visible() ) TR << "Set bin params file name to \"" << binfile_ << "\"." << std::endl;
	bin_transition_calculator_->load_bin_params(binfile_);
	binfile_loaded_=true;
	return;
} //set_binfile_and_load()

/// @brief Set the residue ranges.  If set to (0,0), the
/// start and end of the pose are used as the range bounds.
void InitializeByBins::set_residue_range( core::Size const start, core::Size const end )
{
	if ( (end != 0) && (end < start) ) {
		utility_exit_with_message(
			"In protocols::simple_moves::bin_transitions::InitializeByBins::set_residue_range(): The end residue cannot be before the start.");
	}
	start_res_=start;
	end_res_=end;
	if ( TR.visible() ) {
		TR << "Set bin start and end ranges to ";
		if ( start_res_!=0 ) TR << start_res_; else TR << "[start of pose]";
		TR << ", ";
		if ( end_res_!=0 ) TR << end_res_; else TR << "[end of pose]";
		TR << std::endl;
	}
	return;
} //set_residue_range

std::string InitializeByBins::get_name() const {
	return mover_name();
}

std::string InitializeByBins::mover_name() {
	return "InitializeByBins";
}

void InitializeByBins::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("bin_params_file", xs_string,
		"File identifier indicating which set of bin parameters is desired.", "ABBA" )
		+ XMLSchemaAttribute::attribute_w_default("start", xsct_non_negative_integer, "First residue of the backbone stretch of interest.", "0" )
		+ XMLSchemaAttribute::attribute_w_default("end", xsct_non_negative_integer, "Last residue of the backbone stretch of interest.", "0" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Set mainchain torsions based on the probability of occurance"
		"given the torsion angles of the preceeding amino acid.", attlist );
}

std::string InitializeByBinsCreator::keyname() const {
	return InitializeByBins::mover_name();
}

protocols::moves::MoverOP
InitializeByBinsCreator::create_mover() const {
	return protocols::moves::MoverOP( new InitializeByBins );
}

void InitializeByBinsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InitializeByBins::provide_xml_schema( xsd );
}



} // bin_transitions
} // simple_moves
} // protocols
