// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/enzdes/EnzdesFixBBProtocol.hh
///
/// @brief
/// @author Florian Richter


#ifndef INCLUDED_protocols_enzdes_EnzdesFixBBProtocol_hh
#define INCLUDED_protocols_enzdes_EnzdesFixBBProtocol_hh


#include <protocols/enzdes/EnzdesBaseProtocol.hh>

//#include <protocols/ligand_docking/LigandBaseProtocol.hh>
//#include <core/scoring/EnergyMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


//#include <core/pack/task/PackerTask.fwd.hh>
//#include <core/chemical/ResidueTypeSet.fwd.hh>


namespace protocols {
namespace enzdes {

class EnzdesFixBBProtocol;
typedef utility::pointer::shared_ptr< EnzdesFixBBProtocol > EnzdesFixBBProtocolOP;

class EnzdesFixBBProtocol : public protocols::enzdes::EnzdesBaseProtocol
{

public:

	EnzdesFixBBProtocol();
	~EnzdesFixBBProtocol() override;

	void apply( core::pose::Pose & pose) override;
	std::string get_name() const override;

	static void register_options();

protected:

	bool start_from_random_rb_conf_;

}; //class EnzdesFixBBProtocol


} //namespace enzdes
} //namespace protocols


#endif // INCLUDED_protocols_enzdes_EnzdesFixBBProtocol_HH
