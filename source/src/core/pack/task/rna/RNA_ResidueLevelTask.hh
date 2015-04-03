// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pack/task/rna/RNA_ResidueLevelTask.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pack_task_rna_RNA_ResidueLevelTask_HH
#define INCLUDED_core_pack_task_rna_RNA_ResidueLevelTask_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/pack/task/rna/RNA_ResidueLevelTask.fwd.hh>

namespace core {
namespace pack {
namespace task {
namespace rna {

	class RNA_ResidueLevelTask: public utility::pointer::ReferenceCount {

	public:

		//constructor
		RNA_ResidueLevelTask();

		//destructor
		~RNA_ResidueLevelTask();

	public:

		void set_sample_rna_chi( bool setting ){ sample_rna_chi_ = setting; }
		bool sample_rna_chi() const { return sample_rna_chi_; }

		void set_sample_five_prime_phosphate( bool setting ){ sample_five_prime_phosphate_ = setting; }
		bool sample_five_prime_phosphate() const { return sample_five_prime_phosphate_; }

		void set_sample_three_prime_phosphate( bool setting ){ sample_three_prime_phosphate_ = setting; }
		bool sample_three_prime_phosphate() const { return sample_three_prime_phosphate_; }

		void set_allow_phosphate_virtualization( bool setting ){ allow_phosphate_virtualization_ = setting; }
		bool allow_phosphate_virtualization() const { return allow_phosphate_virtualization_; }

	private:

		bool sample_rna_chi_;
		bool sample_five_prime_phosphate_;
		bool sample_three_prime_phosphate_;
		bool allow_phosphate_virtualization_;

	};

} //rna
} //task
} //pack
} //core

#endif
