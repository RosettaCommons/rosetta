// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/mmtf/ExtraDataEnumManager.hh
/// @brief Enum string/enum functions for pose extra data we will be storing/retrieving.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_core_io_mmtf_ExtraDataEnumManager_hh
#define INCLUDED_core_io_mmtf_ExtraDataEnumManager_hh

#include <core/io/mmtf/ExtraDataEnumManager.fwd.hh>
#include <core/io/mmtf/ExtraDataEnum.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

#include <map>

namespace core {
namespace io {
namespace mmtf {

/// @brief Enum string/enum functions for pose extra data we will be storing/retrieving.
///  ALL enums must be defined in mmtf/ExtraDataEnum.hh
class ExtraDataEnumManager : public utility::VirtualBase {

public:

	ExtraDataEnumManager();
	ExtraDataEnumManager(ExtraDataEnumManager const & src);

	~ExtraDataEnumManager() override;


	ExtraDataEnum
	string_to_enum(std::string const & data_name) const;

	std::string
	enum_to_string(ExtraDataEnum const data_name) const;

	bool
	is_data_type(std::string const & data_name) const;

private:
	void setup_data_names();

private:

	// lookup map from string name to enum type
	std::map< std::string, ExtraDataEnum > string_to_enum_;
	utility::vector1< std::string >    enum_to_string_;
};


} //core
} //io
} //mmtf



#endif //INCLUDED_core_io_mmtf_ExtraDataEnumManager_hh





