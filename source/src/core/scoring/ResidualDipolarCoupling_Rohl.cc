// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ResidualDipolarCoupling_Rohl.cc
/// @brief  Uses NMR RDC for scoring
/// @author Srivatsan Raman

//core headers
#include <core/scoring/ResidualDipolarCoupling_Rohl.hh>
#include <basic/options/option.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

//C++ headers
#include <iostream>
#include <fstream>
#include <string>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rdc.OptionKeys.gen.hh>

/// Utility headers
#include <utility/io/izstream.hh>

#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "core.scoring.ResidualDipolarCoupling_Rohl" );

namespace core {
namespace scoring {

//////////////////////////////////////////////////////
//@brief reads in RDC data file
//////////////////////////////////////////////////////
extern void store_RDC_ROHL_in_pose( ResidualDipolarCoupling_RohlOP rdc_info, core::pose::Pose& pose ) {
	// ////using core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA_ROHL;
	pose.data().set( core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA_ROHL, rdc_info );
}

extern ResidualDipolarCoupling_RohlCOP retrieve_RDC_ROHL_from_pose( core::pose::Pose const& pose ) {
	// ////using core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA_ROHL;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA_ROHL ) ) {
		return utility::pointer::static_pointer_cast< core::scoring::ResidualDipolarCoupling_Rohl const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA_ROHL ) );
	};
	return NULL;
}

extern ResidualDipolarCoupling_RohlOP retrieve_RDC_ROHL_from_pose( core::pose::Pose& pose ) {
	// ////using core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA_ROHL;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA_ROHL ) ) {
		return utility::pointer::static_pointer_cast< core::scoring::ResidualDipolarCoupling_Rohl > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA_ROHL ) );
	};
	return NULL;
}


void ResidualDipolarCoupling_Rohl::read_RDC_file()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	runtime_assert( option[ OptionKeys::in::file::rdc ]().size() );
	std::string filename( option[ OptionKeys::in::file::rdc ]()[ 1 ].name() );
	utility::io::izstream infile( filename.c_str() );
	std::string line;

	// std::cout << "Reading RDC file " << filename << std::endl;
	tr.Info << "Reading RDC file " << filename << std::endl;

	while ( getline( infile, line ) ) {
		std::istringstream line_stream( line );
		std::string atom1, atom2;
		Size res1, res2;
		Real Jdipolar;
		line_stream >> res1 >> atom1 >> res2 >> atom2 >> Jdipolar;
		if ( line_stream.fail() ) {
			tr.Error << "couldn't read line " << line << " in rdc-file " << filename << "\n";
			throw( utility::excn::EXCN_BadInput(" invalid line "+line+" in rdc-file "+filename));
		}

		Real weight( 1.0 );
		line_stream >> weight;
		if ( line_stream.fail() ) {
			tr.Debug << " set weight for RDC " << res1 << " to 1.0 " << std::endl;
			weight = 1.0;
		}

		if ( res1 != res2 ) {
			std::cout << "Dipolar couplings only between atoms in the same residue are acceptable !" << std::endl;
			continue;
		} else {
			Size data_type( get_RDC_data_type( atom1, atom2 ) );
			All_RDC_lines.push_back( RDC_Rohl( data_type, res1, Jdipolar, weight ) );
		}

	}

	//extra weight file?
	if ( option[ OptionKeys::rdc::weights ].user() ) {
		std::string filename( option[ OptionKeys::rdc::weights ]().name() );
		std::ifstream infile( filename.c_str() );
		while ( getline( infile, line ) ) {
			std::istringstream line_stream( line );
			Real weight; Size res1;
			line_stream >> res1 >> weight;
			if ( line_stream.fail() ) {
				tr.Error << "[Error] reading rdc-weight-file " << filename << std::endl;
				throw( utility::excn::EXCN_BadInput(" invalid line "+line+" in rdc-weight-file "+filename));
			}
			for ( RDC_lines::iterator it = All_RDC_lines.begin(); it != All_RDC_lines.end(); ++it ) {
				if ( it->res() == res1 ) it->weight( weight );
			}
		}
	}

}


//////////////////////////////////////////////////////
//@brief returns type of RDC data N-H, HA-CA etc
//////////////////////////////////////////////////////
Size ResidualDipolarCoupling_Rohl::get_RDC_data_type(
	std::string const & atom1,
	std::string const & atom2
)
{

	Size RDC_type(0);

	if ( ( atom1 == "N" && atom2 == "H" ) || ( atom1 == "H" && atom2 == "N" ) ) {
		RDC_type = 1;
	} else if ( ( atom1 == "C" && atom2 == "H" ) || ( atom1 == "H" && atom2 == "C" ) ) {
		//****************** FIX THESE LATER !! **************************
		RDC_type = 2;
	} else if ( ( atom1 == "C" && atom2 == "N" ) || ( atom1 == "N" && atom2 == "C" ) ) {
		RDC_type = 4;
	} else if ( ( atom1 == "C" && atom2 == "C" ) ) {
		RDC_type = 5;
	} else if ( atom1 == "H" && atom2 == "H" ) {
		RDC_type = 6;
	}

	// std::cout << "RDC_type " << RDC_type << std::endl;

	debug_assert(RDC_type != 0 );

	return RDC_type;

}


} //namespace Scoring
} //namespace core
