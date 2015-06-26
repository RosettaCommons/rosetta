// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/movers/dumpStatsSS
/// @author TJ Brunette

#ifndef INCLUDED_protocols_simple_moves_DumpStatsSS_hh
#define INCLUDED_protocols_simple_moves_DumpStatsSS_hh

#include <devel/denovo_design/DumpStatsSS.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/datacache/DataMap.fwd.hh>


#include <protocols/denovo_design/filters/PsiPredInterface.fwd.hh>    
#include <protocols/ss_prediction/SS_predictor.fwd.hh> 
#include <protocols/jd2/parser/BluePrint.fwd.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

namespace devel {
namespace denovo_design {
	

/// @brief what you think
/// this can now be assimilated into DumpPdbMover
class DumpStatsSS : public protocols::moves::Mover
{
public:
	DumpStatsSS();
	DumpStatsSS( DumpStatsSS const &rval);
	virtual ~DumpStatsSS();
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const {
		return( protocols::moves::MoverOP( new DumpStatsSS( *this ) ) );
	}
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new DumpStatsSS ); }
	core::Real compute_boltz_sum( utility::vector1< core::Real > const & probabilities ) const;
	core::Real compute_svm_prob(std::string sequence, std::string wanted_ss); 
	core::Real compute_psipred_prob(core::pose::Pose & pose , std::string wanted_ss);
	void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn);
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

private:
	std::string fname_;
	utility::io::ozstream* output_;
	core::scoring::ScoreFunctionOP scorefxn_;
	std::string psipred_cmd_;
	protocols::denovo_design::filters::PsiPredInterfaceOP psipred_interface_;
	protocols::ss_prediction::SS_predictorOP ss_predictor_;
	protocols::jd2::parser::BluePrintOP blueprint_;
	core::Real start_time_;
};

}//devel 
}//denovo_design 


#endif /*INCLUDED_protocols_moves_DumpStatsSS_HH*/
