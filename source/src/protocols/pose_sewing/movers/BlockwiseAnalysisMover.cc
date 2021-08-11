// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/movers/BlockwiseAnalysisMover.cc
/// @brief analyzes pairs of interacting elements
/// @author Frank Teets (frankdt@email.unc.edu)

// Unit headers
#include <protocols/pose_sewing/movers/BlockwiseAnalysisMover.hh>
#include <protocols/pose_sewing/movers/BlockwiseAnalysisMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/methods/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/pose_sewing/util.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <protocols/simple_filters/ShapeComplementarityFilter.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Citation Manager
#include <utility/vector1.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

static basic::Tracer TR( "protocols.pose_sewing.movers.BlockwiseAnalysisMover" );

namespace protocols {
namespace pose_sewing {
namespace movers {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
BlockwiseAnalysisMover::BlockwiseAnalysisMover():
	protocols::moves::Mover( BlockwiseAnalysisMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
BlockwiseAnalysisMover::~BlockwiseAnalysisMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
BlockwiseAnalysisMover::apply( core::pose::Pose& pose){

	core::scoring::dssp::Dssp dssp(pose);
	dssp.insert_ss_into_pose(pose);
	if ( pose.pdb_info() == nullptr ) {
		core::pose::PDBInfoOP new_info = core::pose::PDBInfoOP(new core::pose::PDBInfo(pose,true));
		pose.pdb_info(new_info);
	}

	std::map< core::Size, core::Size > element_blocks;
	calculate_helices(element_blocks,pose,5); // key is residue, value is block

	std::set<std::pair<core::Size, core::Size>> helix_pairs;

	std::pair<core::Size, core::Size> helix_pair;

	for ( core::Size upstream_res = 1; upstream_res < pose.size(); ++upstream_res ) {
		if ( pose.secstruct(upstream_res) == 'H' ) {
			for ( core::Size downstream_res = upstream_res+1; downstream_res <= pose.size(); ++downstream_res ) {
				if ( pose.secstruct(downstream_res) == 'H' ) {
					if ( element_blocks[upstream_res] != element_blocks[downstream_res] && pose.residue(upstream_res).xyz(2).distance(pose.residue(upstream_res).xyz(2)) <= crit_dist_ ) {
						helix_pair.first = element_blocks[upstream_res];
						helix_pair.second = element_blocks[downstream_res];
						helix_pairs.insert(helix_pair);
					}
				}

			}
		}
	}
	//helix_pairs now has all the interacting helix pairs we care about, so compare them here

	core::scoring::sc::ShapeComplementarityCalculator scc;
	core::Real sc;
	core::Real area;
	core::Real d_median;
	bool has_disulfide;
	std::string label;
	//core::Size label_res = 1;
	for ( auto current_pair : helix_pairs ) {
		scc.Reset(); // this may not be needed anymore, but I'm leaving it here for safety
		for ( core::Size current_res = 1; current_res < pose.size(); ++current_res ) {
			if ( element_blocks[current_res] == current_pair.first ) {
				scc.AddResidue( 0, pose.residue(current_res) );
			}
			if ( element_blocks[current_res] == current_pair.second ) {
				scc.AddResidue( 1, pose.residue(current_res) );
			}

		}
		// subsets have been set
		scc.Calc();
		core::scoring::sc::RESULTS const & r = scc.GetResults();
		sc = r.sc;
		area = r.area;
		d_median = r.distance;
		has_disulfide = false;
		core::Size last_upstream_res = 1;
		for ( core::Size upstream_res = 1; upstream_res <pose.size(); ++upstream_res ) {
			if ( pose.residue(upstream_res).type().is_disulfide_bonded() && element_blocks[upstream_res] == current_pair.first ) {
				for ( core::Size downstream_res = upstream_res+2; downstream_res <= pose.size(); ++downstream_res ) {
					if ( pose.residue(downstream_res).type().is_disulfide_bonded() && pose.residue(upstream_res).is_bonded(pose.residue(downstream_res)) && element_blocks[downstream_res] == current_pair.second ) {
						has_disulfide = true;
						last_upstream_res = upstream_res;
					}
				}
			}
		}
		if ( has_disulfide ) {
			label = "BWA_DISULFIDE ";
		} else {
			label = "BWA_NOTDISULF ";
		}
		label = label + std::to_string(sc) + " " + std::to_string(area) + " " + std::to_string(d_median);

		TR << label << std::endl;
		pose.pdb_info()->add_reslabel(last_upstream_res,label);


	}

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
BlockwiseAnalysisMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
BlockwiseAnalysisMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap
) {
	if ( tag->hasOption("scorefxn") ) {
		scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, datamap );
	}

}
void BlockwiseAnalysisMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_score_function( attlist );

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "DOCUMENTATION STRING", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BlockwiseAnalysisMover::fresh_instance() const
{
	return utility::pointer::make_shared< BlockwiseAnalysisMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BlockwiseAnalysisMover::clone() const
{
	return utility::pointer::make_shared< BlockwiseAnalysisMover >( *this );
}

std::string BlockwiseAnalysisMover::get_name() const {
	return mover_name();
}

std::string BlockwiseAnalysisMover::mover_name() {
	return "BlockwiseAnalysisMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
BlockwiseAnalysisMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< BlockwiseAnalysisMover >();
}

std::string
BlockwiseAnalysisMoverCreator::keyname() const
{
	return BlockwiseAnalysisMover::mover_name();
}

void BlockwiseAnalysisMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BlockwiseAnalysisMover::provide_xml_schema( xsd );
}

/// @brief Indicate that this mover is unpublished.
void
BlockwiseAnalysisMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"InterocitorMover", basic::citation_manager::CitedModuleType::Mover,
		"Frank Teets",
		"Institute for Protein Innovation",
		"teetsf@gmail.com",
		"Wrote the InterocitorMover."
		)
	);
}


////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, BlockwiseAnalysisMover const & mover )
{
	mover.show(os);
	return os;
}


} //movers
} //pose_sewing
} //protocols
