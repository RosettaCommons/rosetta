// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author ashworth

#ifndef INCLUDED_protocols_dna_PDBOutput_hh
#define INCLUDED_protocols_dna_PDBOutput_hh

#include <protocols/dna/PDBOutput.fwd.hh>
#include <protocols/jd2/PDBJobOutputter.hh> // base class
#include <protocols/jd2/Job.fwd.hh>

//#include <protocols/dna/typedefs.hh> // Strings

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>

#include <string>
#include <list>
#include <iosfwd>
#include <map>

namespace protocols {
namespace dna {

class PDBOutput : public jd2::PDBJobOutputter {
public:
	typedef std::list< std::string > Strings; // needs to match typedef in jd2/Job.hh
	typedef std::map< std::string, Strings > StringsMap;
	typedef jd2::JobCOP JobCOP;
	typedef jd2::JobOP JobOP;
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;
	typedef core::pose::PoseCOP PoseCOP;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;
	typedef core::conformation::Residue Residue;
	typedef utility::io::ozstream ozstream;

public:
	PDBOutput();
	~PDBOutput() override;
	/// @brief JobDistributor calls this method
	void final_pose( JobOP, Pose const &, std::string const & ) override;
	/// @brief functor for non-JobDistributor usage
	void operator() ( Pose const &, std::string const & );

	void enabled( bool value ) const { enabled_ = value; }
	bool enabled() const { return enabled_; }

	void starting_pose( Pose const & ) override;
	virtual void reference_pose( Pose const & );
	PoseCOP reference_pose() const;

	void score_function( ScoreFunction const & sf );
	ScoreFunctionCOP score_function() const;

	void add_info( std::string const &, Strings const &, bool append = true );
	bool remove_info( std::string const & );

	void designed_residue( core::Size, bool value = true );
	void note_designed_residues( PackerTaskCOP );
	bool residues_are_different( Residue const &, Residue const & ) const;

private: // methods
	void get_residue_indices_to_output();
	void output_pdb( ozstream & );
	void output_info( ozstream & );
	void output_score_info( ozstream & );
	void output_hbond_info( ozstream & );
	void output_buried_unsatisfied_hbonds( ozstream & );
	void output_design_tags( ozstream & ) const;

private: // data
	PoseOP pose_copy_; // to allow rescoring (nonconst)
	PoseCOP reference_pose_;
	ScoreFunctionOP score_function_;
	core::Real chi_diff_threshold_, mainchain_torsion_diff_threshold_;
	mutable bool enabled_;
	StringsMap info_map_; // deprecated in favor of jd2 framework (extract_data_from_Job)
	utility::vector1< core::Size > res_indices_to_output_;
	utility::vector1< bool > designed_residues_;
};

void make_subdirs( std::string const & );

} // namespace dna
} // namespace protocols

#endif
