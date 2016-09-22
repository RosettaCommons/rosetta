// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/DeNovoArchitectFactory.hh
/// @brief Creates DeNovo architects
/// @author Tom Linsky (tlinsky@gmail.com)

#ifndef INCLUDED_protocols_denovo_design_architects_DeNovoArchitectFactory_hh
#define INCLUDED_protocols_denovo_design_architects_DeNovoArchitectFactory_hh

// Unit headers
#include <protocols/denovo_design/architects/DeNovoArchitectFactory.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/architects/DeNovoArchitectCreator.hh>

// Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/SingletonBase.hh>

// C++ headers
#include <map>

namespace protocols {
namespace denovo_design {
namespace architects {

///@brief Creates DeNovo architects
class DeNovoArchitectFactory : public utility::SingletonBase< DeNovoArchitectFactory > {
private:
	typedef std::map< std::string, DeNovoArchitectCreatorCOP > ArchitectCreatorMap;

public:
	DeNovoArchitectFactory();

	virtual ~DeNovoArchitectFactory();

public:
	DeNovoArchitectOP
	create_instance(
		std::string const & architect_name,
		std::string const & architect_id ) const;

	DeNovoArchitectOP
	create_from_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) const;

	void
	add_creator( DeNovoArchitectCreatorOP creator );

private:
	ArchitectCreatorMap creators_;

};

} //protocols
} //denovo_design
} //architects

#endif //INCLUDED_protocols_denovo_design_architects_DeNovoArchitectFactory_hh
