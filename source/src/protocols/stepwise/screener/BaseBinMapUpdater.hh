// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/BaseBinMapUpdater.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_BaseBinMapUpdater_HH
#define INCLUDED_protocols_stepwise_screener_BaseBinMapUpdater_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/BaseBinMapUpdater.fwd.hh>
#include <protocols/stepwise/sampling/rna/rigid_body/FloatingBaseClasses.hh>
#include <utility/vector1.fwd.hh>


using namespace protocols::stepwise::sampling::rna::rigid_body;
using namespace core;

namespace protocols {
namespace stepwise {
namespace screener {

	class BaseBinMapUpdater: public StepWiseScreener {

	public:

		//constructor
		BaseBinMapUpdater( BaseBinMap & base_bin_map );

		//destructor
		~BaseBinMapUpdater();

	public:

		virtual
		void
		get_update( rotamer_sampler::RotamerSamplerBaseOP sampler );

		virtual
		std::string
		name() const { return "BaseBinMapUpdater"; }

		virtual
		StepWiseScreenerType
		type() const { return BASE_BIN_MAP; }

	private:

		void
		update_base_bin_map( BaseBin const & base_bin );

		void
		update_base_bin_map( utility::vector1< Real > const & rigid_body_values );

	private:

		BaseBinMap & base_bin_map_;
	};

} //screener
} //stepwise
} //protocols

#endif
