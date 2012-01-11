// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file antibody2/moves/GraftMover2.hh
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)



#ifndef INCLUDED_protocols_antibody2_moves_GraftMover2_hh
#define INCLUDED_protocols_antibody2_moves_GraftMover2_hh

#include <protocols/antibody2/GraftMover2.fwd.hh>

// Rosetta headers
#include <core/pose/Pose.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/MinMover.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/simple_moves/PackRotamersMover.fwd.hh>
#include <protocols/moves/PyMolMover.fwd.hh>

#include <core/pack/task/TaskFactory.fwd.hh>

#include <protocols/antibody2/AntibodyInfo.hh>
#include <protocols/antibody2/AntibodyInfo.fwd.hh>
#include <protocols/antibody2/Ab_TemplateInfo.fwd.hh>

// ObjexxFCL Headers

// C++ Headers

// Utility Headers

namespace protocols {
namespace antibody2 {

	//////////////////////////////////////////////////////////////////////////
	/// @brief Grafts a series of CDR onto a framework
	/// @details
	class GraftMover2 : public protocols::moves::Mover {
	public:
		typedef std::map < std::string, bool > GraftMap;
		// default constructor
		GraftMover2();

		/// @brief constructor with arguments
		GraftMover2(bool l1,bool l2,bool l3,bool h1,bool h2,bool h3,bool camelid,bool benchmark);

        
        GraftMover2(AntibodyInfoCOP ab_info, Ab_TemplateInfoCOP ab_template);
		// default destructor
		~GraftMover2();

		void init(bool l1,bool l2,bool l3,bool h1,bool h2,bool h3,bool camelid,bool benchmark);

		inline void enable_graft_l1( bool setting ) { graft_l1_ = setting; }
		inline void enable_graft_l2( bool setting ) { graft_l2_ = setting; }
		inline void enable_graft_l3( bool setting ) { graft_l3_ = setting; }
		inline void enable_graft_h1( bool setting ) { graft_h1_ = setting; }
		inline void enable_graft_h2( bool setting ) { graft_h2_ = setting; }
		inline void enable_graft_h3( bool setting ) { graft_h3_ = setting; }
		inline void set_camelid( bool setting ) { camelid_ = setting; }

	        ///@brief copy ctor
        	GraftMover2( GraftMover2 const & rhs );

        	///@brief assignment operator
        	GraftMover2 & operator=( GraftMover2 const & rhs );
        
		/// @brief enable benchmark mode
		inline void enable_benchmark_mode( bool setting ) {
			benchmark_ = setting;
		}

		/// @brief relax optimized CDR grafted regions
		void relax_optimized_CDR_grafts( core::pose::Pose & pose );

		void set_default();
		virtual void apply( core::pose::Pose & pose );
		virtual std::string get_name() const;

	private:
		// Grafting options
		bool graft_l1_;
		bool graft_l2_;
		bool graft_l3_;
		bool graft_h1_;
		bool graft_h2_;
		bool graft_h3_;

		GraftMap grafts_;

		/// @brief benchmark flag
		bool benchmark_;
		bool camelid_;

		bool user_defined_;
		bool first_apply_with_current_setup_;

		// movers
		protocols::moves::SequenceMoverOP graft_sequence_, relax_sequence_;
		protocols::simple_moves::PackRotamersMoverOP packer_;
		protocols::moves::PyMolMoverOP pymol_;

		core::scoring::ScoreFunctionOP scorefxn_;

		void finalize_setup( core::pose::Pose & framework, AntibodyInfo & ab_info );
		void set_packer_default( core::pose::Pose & pose, bool include_current );
        void initForEqualOperatorAndCopyConstructor(GraftMover2 & lhs, GraftMover2 const & rhs);
	}; // class GraftMover2






} // antibody2
} // protocols


#endif
