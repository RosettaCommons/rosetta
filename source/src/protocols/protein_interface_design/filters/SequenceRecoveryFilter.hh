// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ContingentFilter.hh
/// @brief A filter that is contingent on some other mover to set its pass/fail value
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_SequenceRecoveryFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_SequenceRecoveryFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/protein_interface_design/filters/SequenceRecoveryFilter.fwd.hh>

#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/vector1.hh>

// Unit headers

namespace protocols {
namespace protein_interface_design {
namespace filters {

class SequenceRecoveryFilter : public protocols::filters::Filter
{
private:
	typedef protocols::filters::Filter parent;
public:
	/// @brief default ctor
	SequenceRecoveryFilter();
	/// @brief Constructor with a single target residue
	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	protocols::filters::FilterOP clone() const override;
	protocols::filters::FilterOP fresh_instance() const override;
	core::Real compute( core::pose::Pose const & pose, bool const & write ) const;
	void write_to_pdb(
		std::map< core::Size, std::string > const & res_names1,
		std::map< core::Size, std::string > const & res_names2 ) const;
	virtual ~SequenceRecoveryFilter();
	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	core::Real rate_threshold() const;
	void rate_threshold( core::Real const rate );
	core::Size mutation_threshold() const;
	void mutation_threshold( core::Size const mut );
	bool mutations() const;
	void mutations( bool const muts );
	bool verbose() const;
	void verbose( bool const verb );
	bool write2pdb() const;
	void write2pdb( bool const write );
	void scorefxn( core::scoring::ScoreFunctionCOP sfx );
	core::pose::PoseCOP reference_pose() const;
	void reference_pose( core::pose::PoseCOP reference_pose );
	void reference_pose( core::pose::Pose const & pose );
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::pack::task::TaskFactoryOP task_factory_;
	core::pose::PoseCOP reference_pose_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Real rate_threshold_;
	core::Size mutation_threshold_;
	bool mutations_;
	bool verbose_;
	bool write2pdb_;
};

} // filters
} //protein_interface_design
} // protocols

#endif //INCLUDED_protocols_Filters_SequenceRecoveryFilter_HH_
