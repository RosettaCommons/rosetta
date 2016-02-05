// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file protocols/forge/remodel/RemodelGlobalFrame.hh
/// @brief
/// @author Possu Huang ( possu@uw.edu )


#ifndef INCLUDED_protocols_forge_remodel_RemodelGlobalFrame_hh
#define INCLUDED_protocols_forge_remodel_RemodelGlobalFrame_hh

//project headers
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/forge/remodel/RemodelData.hh>
#include <protocols/forge/remodel/RemodelWorkingSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using namespace std;

namespace protocols {
namespace forge {
namespace remodel {


numeric::xyzVector< core::Real >
compute_center_of_mass( core::pose::Pose const &  pose, core::Size range_start, core::Size range_stop);

double get_RMSD(Eigen::MatrixXf &A,Eigen::MatrixXf &B);
Eigen::MatrixXf ideal_COMs(double rise, double radius, double omega, int unitn);
Eigen::Matrix3f rot_mat(Eigen::MatrixXf &A,Eigen::MatrixXf &B);


class RemodelGlobalFrame: public protocols::moves::Mover {

private: // typedefs

	typedef protocols::moves::Mover Super;

public: // typedefs

	typedef core::Real Real;
	typedef core::Size Size;

	typedef core::kinematics::MoveMap MoveMap;
	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef core::pack::task::PackerTask PackerTask;
	typedef protocols::moves::MoverOP MoverOP;
	typedef protocols::forge::remodel::RemodelData RemodelData;
	typedef protocols::forge::remodel::RemodelWorkingSet RemodelWorkingSet;
	typedef core::scoring::constraints::ConstraintSetOP ConstraintSetOP;
	typedef core::scoring::constraints::ConstraintSet ConstraintSet;


public: //constructor/destructor

	RemodelGlobalFrame();

	RemodelGlobalFrame(RemodelData const & remodel_data, RemodelWorkingSet const & working_model, ScoreFunctionOP const & sfxn);

	RemodelGlobalFrame(Size segment_size);

	virtual
	~RemodelGlobalFrame();

public: // virtual constructors

	virtual
	MoverOP clone() const;

	virtual
	MoverOP fresh_instance() const;

public: // options

public:

	virtual void apply( Pose & pose);
	virtual std::string get_name() const;

	void get_helical_params( Pose & pose );
	void align_segment( Pose & pose );
	void setup_helical_constraint( Pose & pose );
	void setup_CM_helical_constraint( Pose & pose );
	void matrix3f_to_xyzMatrix( Eigen::Matrix3f const & Re, numeric::xyzMatrix< core::Real> & R  );
	void identity_matrix( numeric::xyzMatrix< core::Real>  & R );
	void restore_original_cst( Pose & pose );
	void set_native_cst_set( ConstraintSet const & cst_set);
	void set_native_cst_set( Pose const & pose );
	void set_segment_size( Size segment_size );
	Real radius(){ return radius_; };
	Real rise(){ return rise_; };
	Real omega(){ return omega_; };


private: // data

	RemodelData remodel_data_; // design mode determined in here
	RemodelWorkingSet working_model_; // data for the remodeling pose
	ScoreFunctionOP score_fxn_;

	Real radius_;
	Real rise_;
	Real omega_;
	Size seg_size;
	int left_handed_;
	ConstraintSetOP native_cst_set;

public: // accessors

	void scorefunction(ScoreFunctionOP const & sfxn);

	// MoveMap const & movemap() const;

	// PackerTask const & packertask() const;

};

} // remodel
} // forge
} // protocols

#endif /* INCLUDED_protocols_forge_remodel_RemodelGlobalFrame_HH */
