// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/AddJobPairData.hh
/// @brief  Header file for the AddJobPairData Mover
/// @author Sam DeLuca <Samuel.l.deluca@vanderbilt.edu)


#ifndef INCLUDED_protocols_simple_moves_AddJobPairData_hh
#define INCLUDED_protocols_simple_moves_AddJobPairData_hh

#include <protocols/simple_moves/AddJobPairData.fwd.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace simple_moves {

/// @brief Add an arbitrary piece of data to the current Job, which will be output in the silent file, database, etc.
/// This is useful for adding metadata to keep track of data generated using multiple experimental conditions.
/// Currently ONLY RosettaScript interface supported.
/// @details The data appended to the Job consists of a key and a value. The key is a string, and the value can be either a real or string.
///
class AddJobPairData : public moves::Mover {
public:
	AddJobPairData();
	AddJobPairData(AddJobPairData const & src);
	~AddJobPairData() override;
	void apply( Pose & ) override;
	// XRW TEMP  std::string get_name() const override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap &,
		Filters_map const &,
		moves::Movers_map const &,
		Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	enum ValueType {string_value,real_value};

	std::string string_key_;
	std::string string_value_;
	core::Real real_value_;
	ValueType value_type_;
	std::string ligand_chain_;


};

}
}


#endif /* ADDJOBSTRINGPAIR_HH_ */
