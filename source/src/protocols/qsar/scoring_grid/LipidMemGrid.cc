// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    src/protocols/qsar/scoring_grid/LipidMemGrid.cc
/// @brief   implementation of LipidMemGrid
/// @details Modified 12/11/17
/// @author  Brennica Marlow (brennica.marlow@vanderbilt.edu)

//Unit headers
#include <protocols/qsar/scoring_grid/LipidMemGrid.hh>
#include <protocols/qsar/scoring_grid/LipidMemGridCreator.hh>


//Package Headers
#include <protocols/qsar/scoring_grid/schema_util.hh>
#include <core/scoring/func/SplineFunc.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>
#include <protocols/membrane/OptimizeMembranePositionMover.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/membrane/AddMembraneMover.hh>
//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>

//Pose Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <utility/string_util.hh>
#include <core/conformation/Atom.hh>

//Utility headers
#include <utility/exit.hh>
#include <utility/tools/make_vector.hh>
#include <utility/json_spirit/json_spirit_value.h>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <protocols/membrane/util.hh>
#include <utility/string_util.hh>


//Numeric headers
#include <numeric/interpolation/spline/SplineGenerator.hh>
#include <numeric/interpolation/spline/Interpolator.hh>
#include <numeric/interpolation/util.hh>
#include <numeric/xyzVector.io.hh>


static basic::Tracer TR("protocols.qsar.scoring_grid.LipidMemGrid");


namespace protocols {
namespace qsar {
namespace scoring_grid {


std::string LipidMemGridCreator::keyname() const
{
	return LipidMemGrid::grid_name();
}

/// @brief Setup a grid based on RosettaScripts input
GridBaseOP LipidMemGridCreator::create_grid(utility::tag::TagCOP tag) const
{
	GridBaseOP lipid_mem_grid(new LipidMemGrid());
	lipid_mem_grid->parse_my_tag(tag);
	return lipid_mem_grid;
}

GridBaseOP LipidMemGridCreator::create_grid() const
{
	return GridBaseOP(new LipidMemGrid());
}

void LipidMemGridCreator::provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd) const
{
	LipidMemGrid::provide_xml_schema(xsd);
}

/// @brief Copy SingleGrid functions
LipidMemGrid::LipidMemGrid() : SingleGrid("LipidMemGrid"){}

/// @brief Destructor
LipidMemGrid::~LipidMemGrid() {}

/// @brief Make a copy of the grid, respecting the subclass
GridBaseOP LipidMemGrid::clone() const
{
	return GridBaseOP(new LipidMemGrid(*this));
}

/// @bried Set Grid name
std::string LipidMemGrid::grid_name()
{
	return "LipidMemGrid";
}


void LipidMemGrid::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute( "grid_name", xs_string, "The name used to insert the scoring grid into the GridSet." )
		+ XMLSchemaAttribute::required_attribute("lipid_atom", xs_string, "Lipid atom name used to represents the head group or tail group of the lipid. Example: For cholesterol the headgroup atom is O1 and the tail atom is C25")
		+ XMLSchemaAttribute::required_attribute("kbpot_file", xs_string, "Text file containing lipid z coordinate statistical potentials.")
		+ XMLSchemaAttribute::attribute_w_default("mem_weight", xsct_real, "Multiply the LipidMemGrid score with this weight.","2.0");

	xsd_type_definition_w_attributes(xsd, grid_name(), "A scoring grid that accounts for lipid directionality in the membrane space.", attributes);
}


/// @brief populate the grid with LipidMemGrid values in vector based on a passed pose
void LipidMemGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size>) {
	refresh(pose,center);
}

/// @brief populate the grid with LipidMemGrid values in vector based on a passed pose
void LipidMemGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const &) {
	refresh(pose,center);
}



/// @brief set grid and fill with lipid energies from spline
void LipidMemGrid::refresh(core::pose::Pose const & pose, core::Vector const &){
	this->fill_with_value(0.0);

	using namespace core::scoring;
	using namespace protocols::membrane;
	using namespace protocols::membrane::geometry;


	/// Copy pose
	core::pose::Pose pose_copy(pose);

	/// Check if text file with lipid energies is provided
	if ( !kbpot_spline_ ) {
		utility_exit_with_message("ERROR kbpot_file was not provided");
	} else {

		/// Get membrane potential of the pose
		/// Tried using mpframework_smooth_fa_2012.wts first, but better results with franklin2019.wts
		/// ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function("franklin2019.wts");

		/// Check if pose membrane already set
		if ( !pose_copy.membrane_info() ) {
			TR.Warning << "Membrane not set, RosettaMP will set membrane" << std::endl;

			/// Add the membrane to the pose
			AddMembraneMoverOP addmem(new AddMembraneMover());
			addmem->apply(pose_copy);
		}


		/// Compute membrane embeddingi of the pose
		core::Vector mem_center;
		core::Vector mem_normal;
		compute_structure_based_embedding( pose_copy, mem_center, mem_normal);


		/// Get dimensions of grid
		numeric::xyzVector<core::Size> dimensions = get_dimensions();


		/// Loop over grid in 3D
		for ( int x_index(0), x_stop(static_cast<int>(dimensions.x()-1)); x_index <= x_stop; ++x_index ) {
			for ( int y_index(0), y_stop(static_cast<int>(dimensions.y()-1)); y_index <= y_stop; ++y_index ) {
				for ( int z_index(0), z_stop(static_cast<int>(dimensions.z()-1)); z_index <= z_stop; ++z_index ) {

					/// Get xyz coordinates of voxels in grid
					numeric::xyzVector<int> grdPt(x_index, y_index, z_index);

					/// Get xyz coordinates at the center of each voxel
					core::Vector box_center = get_grid().coords(grdPt);


					/// Calculate membrane depth with respect to the membrane normal for each voxel
					core::Real membrane_proj = dot(box_center - mem_center, mem_normal);
					numeric::Real score, dscore;

					/// Interpolate spline with membrane depth and corresponding energy
					kbpot_spline_->interpolate(membrane_proj, score, dscore);

					/// Set the scores of each voxel in the grid with interpolated spline values
					this->set_point(box_center, score);

				}
			}
		}
	}
}

/// @brief return the current score of an the UltraLightResidue (ligand) using the LipidMemGrid
core::Real LipidMemGrid::score( core::conformation::UltraLightResidue const & residue, core::Real const max_score, qsarMapCOP /*qsar_map*/) const {

	core::Real score = 0.0;
	for ( core::Size atom_index = 1; atom_index <= residue.natoms() && score < max_score; ++atom_index ) {
		/// Set lipid coordinates
		core::Vector const & coords(residue[atom_index]);

		/// Check if lipid coordinates are in the grid
		if ( this->get_grid().is_in_grid(coords.x(), coords.y(), coords.z()) ) {
			std::string atoms(utility::stripped_whitespace(residue.residue()->atom_name(atom_index)));
			/// Find the coordinates of the input lipid atom
			if  ( atoms == lip_atom_ ) {
				/// Get the correspinding scores in the grid for the input lipid atom
				core::Real grid_value = get_point(coords.x(),coords.y(),coords.z());

				score = score + grid_value;
			}
		}
	}
	return score * mem_weight_;
}


/// @brief return the current score of an atom of the UltraLighResidue (ligand) using the LipidMemGrid
core::Real LipidMemGrid::atom_score(core::conformation::UltraLightResidue const & residue, core::Size atomno, qsarMapCOP /*qsar_map*/) const {

	core::Real score = 0.0;
	core::Vector const & coords(residue[atomno]);
	/// Check if lipid atom in grid
	if ( this->get_grid().is_in_grid(coords.x(),coords.y(),coords.z()) ) {
		std::string atoms(utility::stripped_whitespace(residue.residue()->atom_name(atomno)));
		/// Get value of the lipid atom in the grid at these coordinates
		if ( atoms == lip_atom_ ) {
			core::Real grid_value = get_point(coords.x(),coords.y(),coords.z());
			score = score + grid_value;
		}
		return score * mem_weight_;
	}
	return 0;
}


/// @brief return the current score of a Residue using the LipidMemGrid
core::Real LipidMemGrid::score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapCOP /*qsar_map*/) const {

	core::Real score = 0.0;
	for ( core::Size atom_index = 1; atom_index <= residue.natoms() && score < max_score; ++atom_index ) {
		/// Get protein residue coordinates
		core::Vector const & coords(residue.xyz(atom_index));
		// Check if protein residues in the grid
		if ( this->get_grid().is_in_grid(coords.x(), coords.y(), coords.z()) ) {
			std::string atoms(utility::stripped_whitespace(residue.atom_name(atom_index)));
			if ( atoms == lip_atom_ ) {
				core::Real grid_value = get_point(coords.x(), coords.y(), coords.z());
				score = score + grid_value;
			}
		}
	}
	return score * mem_weight_;
}


/// @brief return the current score of an atom of the Residue using LipidMemGrid
core::Real LipidMemGrid::atom_score(core::conformation::Residue const & residue, core::Size atomno, qsarMapCOP /*qsar_map*/) const {

	core::Real score = 0.0;
	/// Get protein atom coordinates
	core::Vector const & coords(residue.xyz(atomno));
	/// Check if protein atoms are in the grid
	if ( this->get_grid().is_in_grid(coords.x(),coords.y(),coords.z()) ) {
		std::string atoms(utility::stripped_whitespace(residue.atom_name(atomno)));
		if ( atoms == lip_atom_ ) {
			core::Real grid_value = get_point(coords.x(),coords.y(),coords.z());
			score = score + grid_value;
		}
		return score * mem_weight_;
	}
	return 0;
}



/// @brief Serialize the LipidMemGrid object into a json_spirit Value
utility::json_spirit::Value LipidMemGrid::serialize() const {
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair spline_data("spline",kbpot_spline_->serialize());
	Pair weight_data("mem_weight",Value(mem_weight_));
	Pair lip_data("lip_atom",Value(lip_atom_));
	Pair base_data("base_data",SingleGrid::serialize());

	return Value(utility::tools::make_vector(spline_data,base_data,weight_data,lip_data));
}

/// @brief deserialize a json spirit Value into a LipidMemGrid object
void LipidMemGrid::deserialize(utility::json_spirit::mObject data) {

	numeric::interpolation::spline::InterpolatorOP kbpot_spline(kbpot_spline_->clone());
	kbpot_spline->deserialize(data["spline"].get_obj());
	kbpot_spline_ = kbpot_spline;
	mem_weight_ = data["mem_weight"].get_real();
	lip_atom_ = data["lip_atom"].get_str();
	SingleGrid::deserialize(data["base_data"].get_obj());
}

/// @setup a LipidMemGrid based on RosettaScripts input
void LipidMemGrid::parse_my_tag(utility::tag::TagCOP const tag) {

	if ( tag->hasOption("kbpot_file") ) {

		/// Interpolate kbpot_file string into a spline for use in refresh function
		kbpot_file_ = tag->getOption<std::string>("kbpot_file");
		kbpot_file_ =  basic::database::full_name( "scoring/score_functions/lipid_orientation/"+kbpot_file_+".dat" );
		numeric::interpolation::spline::SplineGenerator kcommon_spline(numeric::interpolation::spline_from_file(kbpot_file_,1));
		kbpot_spline_ = kcommon_spline.get_interpolator();
	}
	if ( tag->hasOption("lipid_atom") ) {
		/// Set lipid atom input for scoring
		lip_atom_ = utility::stripped_whitespace(tag->getOption<std::string >("lipid_atom"));
	}

	if ( tag->hasOption("mem_weight") ) {
		/// Grid weight if used with other grids
		mem_weight_ = tag->getOption<core::Real>("mem_weight");
	}
}


//@brief Print a brief summary about this grid to the provided output stream
void LipidMemGrid::show(std::ostream & out) const {
	out << "LipidMemGrid grid set atom to: " << lip_atom_ << std::endl;
}


/// @brief Return a string representing the settings which don't change based on reinitialization
std::string LipidMemGrid::hash_fingerprint() const {
	std::stringstream ss;
	const char sep('\t');
	ss << grid_name();
	ss << sep << get_type();
	ss << sep << kbpot_file_; 
	ss << sep << kbpot_spline_;
	ss << sep << lip_atom_;
	ss << sep << mem_weight_;
	return ss.str();
}

} // namespace protocols
} // namespace qsar
} // namespace scoring_grid





