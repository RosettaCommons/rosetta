// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--
/// @brief --brief--
/// @author --name-- (--email--)


#ifndef INCLUDED_--path_underscore--_--class--_hh
#define INCLUDED_--path_underscore--_--class--_hh

#include <--path--/--class--.fwd.hh>

#include <core/pack/task/operation/TaskOperation.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>


--namespace--

///@brief --brief--
class --class--: public core::pack::task::operation::TaskOperation {
public:

	--class--();

	--class--(--class-- const & src);

	virtual ~--class--();

	core::pack::task::operation::TaskOperationOP
	clone() const;

	/// @brief Configure from a RosettaScripts XML tag.
	virtual void
	parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & );


	//////////////////////

	virtual
	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const;

	/// @brief Return the name used to construct this TaskOperation from an XML file
	static std::string keyname();

	/// @brief Describe the format of XML file used to initialize this TaskOperation
	static void provide_xml_schema( XMLSchemaDefinition & xsd );

private:


};



--end_namespace--


#endif //INCLUDED_--path--_--class--_fwd_hh

