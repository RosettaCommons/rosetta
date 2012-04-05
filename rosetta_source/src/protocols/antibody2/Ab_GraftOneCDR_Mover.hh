// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody2/Ab_GraftOneCDR_Mover.hh
/// @brief grafts a cdr onto the template of an antibody framework
/// @detailed
/// @author Jianqing Xu (xubest@gmail.com)


#ifndef INCLUDED_protocols_antibody2_Ab_GraftOneCDR_Mover_hh
#define INCLUDED_protocols_antibody2_Ab_GraftOneCDR_Mover_hh


#include <protocols/antibody2/Ab_GraftOneCDR_Mover.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>  
// Needs to be the full header so the scorefxn can default to NULL
#include <protocols/moves/Mover.hh>
#include <protocols/antibody2/AntibodyInfo.fwd.hh>
#include <protocols/antibody2/Ab_TemplateInfo.fwd.hh>

namespace protocols {
namespace antibody2 {

	/// @brief Grafts only one CDR onto a framework
	class Ab_GraftOneCDR_Mover : public protocols::moves::Mover {
	public:
		// default constructor
		Ab_GraftOneCDR_Mover();

		// constructor with arguments
		Ab_GraftOneCDR_Mover(std::string cdr_name, 
                             core::Size query_start, 
                             core::Size query_end, 
                             core::scoring::ScoreFunctionOP scorefxn );
        
		Ab_GraftOneCDR_Mover(std::string cdr_name, 
                             AntibodyInfoOP ab_info,
                             Ab_TemplateInfoOP ab_t_info,
                             core::scoring::ScoreFunctionOP scorefxn );

        
		~Ab_GraftOneCDR_Mover();

		void set_default( std::string template_name );
		virtual void apply( core::pose::Pose & pose_in );
		virtual std::string get_name() const;


        ///@brief copy ctor
    	Ab_GraftOneCDR_Mover( Ab_GraftOneCDR_Mover const & rhs );

    	///@brief assignment operator
    	Ab_GraftOneCDR_Mover & operator=( Ab_GraftOneCDR_Mover const & rhs );


		/// @brief enable benchmark mode
		inline void enable_benchmark_mode( bool setting ) {
			benchmark_ = setting;
		}

	private:
		// Limits of query loop
		core::Size query_start_;
		core::Size query_end_;
        core::Size flank_size_;
		std::string template_name_;

		// Limits of template loop
		core::Size template_start_;
		core::Size template_end_;
		core::pose::Pose template_pose_;

		/// @brief benchmark flag
		bool benchmark_;

		core::scoring::ScoreFunctionOP scorefxn_;
        void initForEqualOperatorAndCopyConstructor(Ab_GraftOneCDR_Mover & lhs, Ab_GraftOneCDR_Mover const & rhs);

	}; // class Ab_GraftOneCDR_Mover





} // antibody2
} // protocols








#endif
