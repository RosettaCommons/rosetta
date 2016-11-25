// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/Torsion.hh
/// @brief Reports the average degree of connectivity of interface residues
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_TorsionFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_TorsionFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/protein_interface_design/filters/TorsionFilter.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <utility/vector1.hh>

// Unit headers

namespace protocols {
namespace protein_interface_design {
namespace filters {

class Torsion : public protocols::filters::Filter
{
private:
	typedef protocols::filters::Filter parent;
public:
	/// @brief default ctor
	Torsion();
	/// @brief Constructor with a single target residue
	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	protocols::filters::FilterOP clone() const override;
	protocols::filters::FilterOP fresh_instance() const override;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~Torsion();
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;
	core::Real lower() const{ return lower_; }
	void lower( core::Real const l ) { lower_ = l; }
	core::Real upper() const{ return upper_; }
	void upper( core::Real const u ) { upper_ = u; }
	core::Size resnum() const{ return resnum_; }
	void resnum( core::Size const r ) { resnum_ = r; }
	void torsion( std::string const & t ) { torsion_ = t; }
	std::string torsion() const { return torsion_; }
	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP tf );
	void task_factory_set( bool const b ) { task_factory_set_ = b; }
	bool task_factory_set() const { return task_factory_set_; }

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	core::Real lower_, upper_;
	core::Size resnum_;
	std::string torsion_;
	core::pack::task::TaskFactoryOP task_factory_;
	bool task_factory_set_;/// if you call this filter from within another routine, rather than through rosetta_scripts, and you want to use a task_factory, set this parameter to true before apply
};

} // filters
} //protein_interface_design
} // protocols

#endif //INCLUDED_protocols_Filters_Torsion_HH_
