// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/ralford/color_by_lipid_type.cc
/// @brief Utility script - fills the bfactor colum with xyz-dependent fractional hydration values. Can be used to
/// visualize the hydration in PyMOL
/// @author Rebecca Alford (ralford3@jhu.edu)


// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/ImplicitLipidInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/Conformation.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/types.hh>
#include <core/pose/util.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>

// C++ Headers
#include <cstdlib>
#include <string>

// C++ Headers
#include <ostream>
#include <fstream>

static basic::Tracer TR( "apps.pilot.ralford.ColorByLipidType" );

class ColorByLipidType : public protocols::moves::Mover {

public:

	ColorByLipidType() :
		protocols::moves::Mover()
	{}

	ColorByLipidType( ColorByLipidType const & /*src*/ ) = default;

	~ColorByLipidType() override = default;

	std::string get_name() const override { return "ColorByLipidType"; }

	void
	apply( core::pose::Pose & pose ) override {

		using namespace protocols::membrane;
		using namespace core::conformation::membrane;

		// add the membrane
		AddMembraneMoverOP add_memb( new AddMembraneMover() );
		add_memb->apply( pose );

		ImplicitLipidInfoOP implicit_lipids( pose.conformation().membrane_info()->implicit_lipids() );

		// Color rsd by hyd as bfactor
		for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			for ( core::Size jj = 1; jj <= pose.residue( ii ).natoms(); ++jj ) {
				core::Real z_n( std::pow( pose.residue(ii).atom( jj ).xyz().z()/15, 10 ) );
				core::Real f_of_z( z_n / (1 + z_n ) );
				//core::Real hyd = implicit_lipids->f_hydration( pose.residue( ii ).atom( jj ).xyz() );
				pose.pdb_info()->bfactor( ii, jj, f_of_z );
			}
		}

		utility::vector1< std::string > temp( utility::string_split( pose.pdb_info()->name(), '/') );
		std::string pdb_name = temp[ temp.size() ].substr(0, temp[ temp.size() ].size()-4 );
		pose.dump_pdb( pdb_name + "_" + implicit_lipids->lipid_composition_name() + "_colored.pdb" );

	}
};

////////////////////////////////////////////////////////////////////////////////////
using ColorByLipidTypeOP = utility::pointer::shared_ptr<ColorByLipidType>;
using ColorByLipidTypeCOP = utility::pointer::shared_ptr<const ColorByLipidType>;
////////////////////////////////////////////////////////////////////////////////////


int
main( int argc, char * argv [] )
{

	using namespace protocols::moves;

	try {

		devel::init( argc, argv );
		ColorByLipidTypeOP color_by_lipid_type( new ColorByLipidType() );
		protocols::jd2::JobDistributor::get_instance()->go( color_by_lipid_type );
		return 0;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "Caught Exception " << e.msg() << std::endl;
		return -1;
	}

}
