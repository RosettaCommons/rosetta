// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/InterfaceScoreCalculator.cc
/// @brief
/// @author Gordon Lemmon (glemmon@gmail.com)

// Unit Headers
#include <protocols/ligand_docking/InterfaceScoreCalculator.hh>
#include <protocols/ligand_docking/InterfaceScoreCalculatorCreator.hh>

#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/pose/util.hh>
#include <protocols/ligand_docking/ligand_scores.hh>
#include <protocols/qsar/scoring_grid/schema_util.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/jd2/util.hh>
#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/map_util.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

static THREAD_LOCAL basic::Tracer TR( "protocols.ligand_docking.ligand_options.InterfaceScoreCalculator" );

/// @brief
InterfaceScoreCalculator::InterfaceScoreCalculator():
	Mover("InterfaceScoreCalculator"),
	chains_(),
	native_(/* NULL */),
	score_fxn_(/* NULL */),
	normalization_function_(/* NULL */),
	compute_grid_scores_(false),
	grid_set_prototype_(/* NULL */),
	prefix_("")
{}

InterfaceScoreCalculator::InterfaceScoreCalculator(InterfaceScoreCalculator const & ) = default;

InterfaceScoreCalculator::~InterfaceScoreCalculator() = default;

protocols::moves::MoverOP InterfaceScoreCalculator::clone() const {
	return protocols::moves::MoverOP( new InterfaceScoreCalculator( *this ) );
}

protocols::moves::MoverOP InterfaceScoreCalculator::fresh_instance() const {
	return protocols::moves::MoverOP( new InterfaceScoreCalculator );
}

void InterfaceScoreCalculator::chains(std::vector<std::string> const & chains)
{
	chains_ = chains;
}


void InterfaceScoreCalculator::score_fxn(core::scoring::ScoreFunctionOP const & score_fxn)
{
	score_fxn_ = score_fxn;
}


void InterfaceScoreCalculator::grid_set_prototype(protocols::qsar::scoring_grid::GridSetCOP grid_prototype)
{
	grid_set_prototype_ = grid_prototype;
	if ( grid_set_prototype_ != nullptr ) {
		compute_grid_scores_ = true; // If we explicitly set the prototype, chances are we want the scores
	} else {
		compute_grid_scores_ = false;
	}
}
/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
InterfaceScoreCalculator::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
)
{
	if ( ! tag->hasOption("chains") ) throw utility::excn::EXCN_RosettaScriptsOption("'InterfaceScoreCalculator' requires 'chains' tag (comma separated chains to dock)");

	std::string const chains_str = tag->getOption<std::string>("chains");
	chains_= utility::string_split(chains_str, ',');

	/// Score Function ///
	if ( ! tag->hasOption("scorefxn") ) throw utility::excn::EXCN_RosettaScriptsOption("'HighResDocker' requires 'scorefxn' tag");
	std::string scorefxn_name= tag->getOption<std::string>("scorefxn");
	score_fxn_= datamap.get_ptr<core::scoring::ScoreFunction>( "scorefxns", scorefxn_name);
	debug_assert(score_fxn_);

	if ( tag->hasOption("native") ) {
		std::string const & native_str= tag->getOption<std::string>("native");
		utility::vector1<std::string> natives_strs= utility::string_split(native_str, ',');
		std::string natives_str = utility::join(natives_strs, " ");

		native_ = core::pose::PoseOP( new core::pose::Pose );
		core::import_pose::pose_from_file(*native_, natives_str, core::import_pose::PDB_file);
	} else if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		std::string const & native_str= basic::options::option[ basic::options::OptionKeys::in::file::native ]().name();
		native_ = core::pose::PoseOP( new core::pose::Pose );
		core::import_pose::pose_from_file(*native_, native_str, core::import_pose::PDB_file);
	}
	if ( tag->hasOption("normalize") ) {
		std::string const & normalization_mode = tag->getOption<std::string>("normalize");
		normalization_function_ = protocols::qsar::scoring_grid::get_score_normalization_function(normalization_mode);
	}
	compute_grid_scores_ = tag->getOption<bool>("compute_grid_scores", true);
	if ( compute_grid_scores_ ) {
		grid_set_prototype_ = protocols::qsar::scoring_grid::parse_optional_grid_set_from_tag(tag, datamap);
		if ( grid_set_prototype_ == nullptr ) {
			if ( tag->hasOption("compute_grid_scores") ) {
				// Explicitly asked for grid scores, but there are none to be had.
				utility::excn::EXCN_RosettaScriptsOption("InterfaceScoreCalculator cannot compute Grid Scores as requested, as the appropriate grid set is not present!");
			} else {
				compute_grid_scores_ = false; // Just ignore it.
			}
		}
	}

	prefix_ = tag->getOption<std::string>("prefix","");


}

void InterfaceScoreCalculator::apply(core::pose::Pose & pose) {

	// This is less than ideal, but at the least it should gracefully degrade if we're not under JD2
	std::map< std::string, std::string > string_string_pairs( protocols::jd2::get_string_string_pairs_from_current_job() );
	if ( string_string_pairs.find("native_path") != string_string_pairs.end() ) {
		std::string native_string(string_string_pairs.find("native_path")->second);
		native_ = core::pose::PoseOP( new core::pose::Pose );
		core::import_pose::pose_from_file(*native_,native_string, core::import_pose::PDB_file);
	}

	// get_scores() is really only useful if prefix_ is not empty,
	// as the output machinery will typically overwrite the un-prefixed scoreterm labels.
	StringRealMap allscores( get_scores( pose ) );

	utility::map_merge( allscores, get_ligand_docking_scores( pose ) );

	// For now, keep current behavior of appending things into the Job.
	// Ideally, these would be attached to the Pose extra scores, rather than the Job.
	for ( auto const & entry: allscores ) {
		protocols::jd2::add_string_real_pair_to_current_job( entry.first, entry.second );
	}

}

InterfaceScoreCalculator::StringRealMap
InterfaceScoreCalculator::get_scores(
	core::pose::Pose & pose
) const
{
	StringRealMap retval;

	debug_assert(score_fxn_);
	using namespace core::scoring;

	core::Real const tot_score = score_fxn_->score( pose );

	// Which score terms to use
	typedef utility::vector1<ScoreType> ScoreTypeVec;
	ScoreTypeVec score_types;
	for ( int i = 1; i <= n_score_types; ++i ) {
		ScoreType ii = ScoreType(i);
		if ( score_fxn_->has_nonzero_weight(ii) ) score_types.push_back(ii);
	}

	for ( ScoreType const & score_type : score_types ) {
		core::Real const weight = score_fxn_->get_weight(score_type);

		std::string term_label( name_from_score_type(score_type) );
		if ( prefix_ != "" ) {
			term_label = prefix_ + "_" + term_label;
		}
		retval[ term_label ] = weight * pose.energies().total_energies()[ score_type ];
	}

	std::string term_label( name_from_score_type(core::scoring::total_score) );
	if ( prefix_ != "" ) {
		term_label = prefix_ + "_" + term_label;
	}
	retval[ term_label ] = tot_score;

	return retval;
}

/// @brief For multiple ligands, append ligand docking scores for each ligand
InterfaceScoreCalculator::StringRealMap
InterfaceScoreCalculator::get_ligand_docking_scores(
	core::pose::Pose const & after
) const
{
	StringRealMap retval;

	for ( std::string const & chain : chains_ ) {
		debug_assert( chain.size() == 1 );
		TR.Debug << "Calculating ligand docking scores for chain: " << chain << std::endl;
		debug_assert( core::pose::has_chain(chain, after) );
		if ( native_ ) {
			if ( !core::pose::has_chain(chain, *native_) ) {
				utility_exit_with_message("The native pose passed to InterfaceScoreCalculator does not have chain " + utility::to_string(chain) );
			}
		}

		utility::map_merge( retval, get_interface_deltas( chain[0], after, score_fxn_, prefix_, normalization_function_ ) );
		utility::map_merge( retval, get_ligand_docking_scores( chain[0], after ) );
	}

	return retval;
}

/// @brief Scores to be output that aren't normal scorefunction terms.
InterfaceScoreCalculator::StringRealMap
InterfaceScoreCalculator::get_ligand_docking_scores(
	char chain,
	core::pose::Pose const & after
) const {

	StringRealMap retval;

	if ( ! core::pose::has_chain( chain, after ) ) {
		utility_exit_with_message("The pose does not have chain " + utility::to_string( chain ) );
	}

	if ( native_ ) {
		if ( ! core::pose::has_chain( chain, *native_ ) ) {
			utility_exit_with_message("The native pose does not have chain " + utility::to_string( chain ) );
		}

		utility::map_merge( retval, get_ligand_travel( chain, after, *native_, prefix_ ) );
		utility::map_merge( retval, get_radius_of_gyration( chain, after, prefix_ ) );
		utility::map_merge( retval, get_ligand_RMSDs(chain, after, *native_, prefix_ ) );
	}

	if ( compute_grid_scores_ ) {
		debug_assert( grid_set_prototype_ != nullptr );
		// get_ligand_grid_scores() won't double normalize if the normalization is enabled in the Grids already.
		utility::map_merge( retval, get_ligand_grid_scores( *grid_set_prototype_, chain, after, prefix_, normalization_function_ ) );
	}

	return retval;
}

std::string InterfaceScoreCalculator::get_name() const {
	return mover_name();
}

std::string InterfaceScoreCalculator::mover_name() {
	return "InterfaceScoreCalculator";
}

void InterfaceScoreCalculator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute( "chains", xs_string , "Comma separated chains to dock." )
		+ XMLSchemaAttribute::required_attribute( "scorefxn", xs_string , "Scorefxn of choice." )
		+ XMLSchemaAttribute( "native", xs_string , "This is your native pdb without interface mutations. If a native structure is specified, 4 additional score terms are calculated: ligand_centroid_travel, ligand_radious_of_gyration, ligand_rms_no_super, and ligand_rms_with_super." )
		+ XMLSchemaAttribute( "normalize", xs_string , "The normalization function you wish to use." )
		+ XMLSchemaAttribute::attribute_w_default( "compute_grid_scores", xsct_rosetta_bool , "If compute_grid_scores is true, the scores for each grid will be calculated. This may result in the regeneration of the scoring grids, which can be slow.", "false" );

	protocols::qsar::scoring_grid::attributes_for_parse_grid_set_from_tag( attlist, "The Grid Set to use when computing grid scores." );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "InterfaceScoreCalculator calculates a myriad of ligand specific scores and appends them to the output file. After scoring the complex the ligand is moved 1000 Å away from the protein. The model is then scored again. An interface score is calculated for each score term by subtracting separated energy from complex energy.", attlist );
}

std::string InterfaceScoreCalculatorCreator::keyname() const {
	return InterfaceScoreCalculator::mover_name();
}

protocols::moves::MoverOP
InterfaceScoreCalculatorCreator::create_mover() const {
	return protocols::moves::MoverOP( new InterfaceScoreCalculator );
}

void InterfaceScoreCalculatorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InterfaceScoreCalculator::provide_xml_schema( xsd );
}



} //namespace ligand_docking
} //namespace protocols
