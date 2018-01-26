// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/rosetta_scripts/RosettaScriptsJobQueen.hh
/// @brief  JD3 JobQueen for RosettaScripts class declaration
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_rosetta_scripts_RosettaScriptsJobQueen_HH
#define INCLUDED_protocols_rosetta_scripts_RosettaScriptsJobQueen_HH

// protocol headers
#include <protocols/jd3/standard/StandardJobQueen.hh>
#include <protocols/jd3/deallocation/DeallocationMessage.fwd.hh>

//#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.fwd.hh>
#include <protocols/rosetta_scripts/util.hh>

// utility headers
#include <utility/vector1.hh>

// Basic headers
#include <basic/resource_manager/ResourceManager.fwd.hh>

// utility headers
#include <utility/vector1.hh>

// C++ headers
#include <list>
#include <map>
#include <set>

namespace protocols {
namespace rosetta_scripts {

class RosettaScriptsJobQueen : public jd3::standard::StandardJobQueen
{
public:
	typedef jd3::standard::StandardJobQueen parent;
public:

	RosettaScriptsJobQueen();
	~RosettaScriptsJobQueen();

	jd3::JobDigraphOP
	initial_job_dag() override;

	jd3::JobOP
	complete_larval_job_maturation(
		jd3::LarvalJobCOP larval_job,
		utility::options::OptionCollectionCOP job_options,
		utility::vector1< jd3::JobResultCOP > const &
	) override;

	/// @brief the RosettaScriptsJobQueen needs to tell remote nodes to deallocate Resources held by the ResourceManager that are no longer needed
	std::list< jd3::deallocation::DeallocationMessageOP >
	deallocation_messages() override;

	/// @brief A deallocation message first sent to the JobDistributor on this host originating from
	/// a remote JobQueen
	void
	derived_process_deallocation_message( jd3::deallocation::DeallocationMessageOP message ) override;

	void
	note_preliminary_job_node_is_complete( core::Size pjn_index ) override;

	void
	initialize_resources_for_preliminary_job_nodes();

	/// @brief Accessor to the JQs resource_manager -- only really needed for testing purposes.
	basic::resource_manager::ResourceManagerOP resource_manager();

private:
	protocols::rosetta_scripts::RosettaScriptsParserOP parser_;
	basic::resource_manager::ResourceManagerOP resource_manager_;
	std::map< std::string, std::set< core::Size > > pjns_requiring_resources_;
	utility::vector1< std::list< std::string > > resources_for_pjn_;
	std::list< jd3::deallocation::DeallocationMessageOP > deallocation_messages_to_send_;
};

}
}
#endif
