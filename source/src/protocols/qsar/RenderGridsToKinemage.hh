// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/qsar/RenderGridsToKinemage.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_RenderGridsToKinemage_HH
#define INCLUDED_protocols_qsar_RenderGridsToKinemage_HH


#include <protocols/moves/Mover.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.fwd.hh>

#include <protocols/qsar/RenderGridsToKinemage.fwd.hh>
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace qsar {

struct ColorGradient
{
	ColorGradient(numeric::xyzVector<core::Real> const & value,
		core::Real const & lower,
		core::Real const & upper,
		std::string const & name );

	numeric::xyzVector<core::Real> color_value;
	core::Real lower_bound;
	core::Real upper_bound;
	std::string color_name;
};

class RenderGridsToKinemage : public protocols::moves::Mover
{
public:

	RenderGridsToKinemage();
	RenderGridsToKinemage(RenderGridsToKinemage const & mover);
	~RenderGridsToKinemage() override;
	moves::MoverOP clone() const override;
	// XRW TEMP  std::string get_name() const override;
	void apply(core::pose::Pose & pose) override;
	void parse_my_tag(utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	void setup_colors();
	void setup_one_color_scheme();
	void setup_two_color_scheme();
	void setup_three_color_scheme();
	void write_points(utility::io::ozstream & kin_file);
	void write_colors(utility::io::ozstream & kin_file);
	void write_header(utility::io::ozstream & kin_file);

private:
	std::string filename_;
	core::Size color_mode_;
	core::Size gradient_bins_;
	core::Size stride_;
	//core::Real low_cut_;
	//core::Real high_cut_;
	std::string grid_name_;
	numeric::xyzVector<core::Real> color_;
	numeric::xyzVector<core::Real> low_color_;
	numeric::xyzVector<core::Real> zero_color_;
	numeric::xyzVector<core::Real> high_color_;
	utility::vector1<ColorGradient> color_data_;
	scoring_grid::SingleGridOP grid_;

};

}
}
#endif
