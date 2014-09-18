// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/saxs/FormFactor.cc
/// @brief Represents an atomic scattering form factor
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <core/scoring/saxs/FormFactor.hh>

// utility headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

#include <numeric/interpolation/spline/Interpolator.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>

#include <utility/io/izstream.hh>

namespace core {
namespace scoring {
namespace saxs {

/// @details Auto-generated virtual destructor
FormFactor::~FormFactor() {}

static thread_local basic::Tracer trFormFactor( "core.scoring.saxs.FormFactor" );

FormFactor::FormFactor(std::string atom_name,std::string file_name) {

	name_ = atom_name;
	utility::vector1<Real> x;
	utility::vector1<Real> y;
	Real minX = 100000000000.0;
	Real minY = 100000000000.0;
	Real maxX = -100000000000.0;
	Real maxY = -100000000000.0;
	Real tX,tY;
	trFormFactor.Debug << "Opening: "<<file_name<<std::endl;
	utility::io::izstream input(file_name.c_str());
	std::string line;
        while( getline( input, line ) ) {
	    if ( line.substr(0,1) == "#" ) continue;
	    std::istringstream line_stream( line );
	    trFormFactor.Trace << "line " << line << std::endl;
	    line_stream >> tX >> tY;
	    if(tX<minX) minX = tX;
	    if(tY<minY) minY= tY;
	    if(tX>maxX) maxX = tX;
	    if(tY>maxY) maxY= tY;
	    x.push_back(tX);
	    y.push_back(tY);
	}
	trFormFactor.Debug << x.size() << " data points loaded for "<<atom_name<<" atomic form factor"<<std::endl;
	Real delta = 0.1;
	numeric::interpolation::spline::SplineGenerator gen( minX-delta, minY-delta, 0, maxX+delta, maxY+delta, 0 );
	for (Size i = 1; i <= x.size(); ++i)
		gen.add_known_value( x[i],y[i] );
	spline_interpolator_ = gen.get_interpolator();
}

void FormFactor::tabulate(const utility::vector1<Real> & q) {
    
    ff_values_.clear();
    for(Size i_s=1;i_s<=q.size();++i_s) {
	ff_values_.push_back(ff(q[i_s]));
    }    
}

} // saxs
} // scoring
} // core

