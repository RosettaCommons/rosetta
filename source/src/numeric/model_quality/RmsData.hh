// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/model_quality/Rms_Data.hh
/// @brief  Small object to hold data for RMS calculations.
/// @author James Thompson
/// @date   Thu Jan 10 06:55:37 2008


#ifndef INCLUDED_numeric_model_quality_RmsData_HH
#define INCLUDED_numeric_model_quality_RmsData_HH

#include <numeric/types.hh>
#include <numeric/xyzVector.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2A.hh>

namespace numeric {
namespace model_quality {

/// @brief RmsData is a class intended to replace the global rms_obj namespace from rosetta++. Initial implementation
/// is with a singleton design pattern to mimic a global namespace from rosetta++.
class RmsData {

	// begin singleton stuff!
public:
	static RmsData* instance();

protected:
	RmsData();
	RmsData(const RmsData&);
	RmsData& operator= (const RmsData&);
	// end singleton stuff!

public:
	////////////////////////////////////////////////////////////////////////////////
	///
	/// @brief  computes a 3x3 matrix of cross moments between the x,y,z components of
	///  the two input vectors.
	///
	/// @details  the output is the running sum of these matricies
	///
	/// @param  i - [in/out]? -
	/// @param  xp - [in/out]? -
	/// @param  xe - [in/out]? -
	///
	/// @global_read
	///
	/// @global_write
	///
	/// @remarks
	///
	/// @references
	///
	/// @author
	///
	/////////////////////////////////////////////////////////////////////////////////
	void add_rms(
		int i,
		ObjexxFCL::FArray2A< double > xp,
		ObjexxFCL::FArray2A< double > xe
	);

	/// @brief clear the data in this RmsData
	void
	clear_rms();

	/// @brief returns the number of points in this RmsData
	int count() {
		return count_;
	}

	ObjexxFCL::FArray1D< double > xsp() {
		return xsp_;
	}

	ObjexxFCL::FArray1D< double > xse() {
		return xse_;
	}

	ObjexxFCL::FArray2D< double > xm() {
		return xm_;
	}

	double xre() {
		return xre_;
	}

	double xrp() {
		return xrp_;
	}

	//void dimension( int npoints );


	/// @brief returns the number of points
private:
	ObjexxFCL::FArray2D< double > xm_; // note zero origin
	double xre_;
	double xrp_;
	ObjexxFCL::FArray1D< double > xse_;
	ObjexxFCL::FArray1D< double > xsp_;
	int count_;

	static RmsData* pinstance_;

}; // class RmsData


} // end namespace model_quality
} // end namespace numeric

#endif
