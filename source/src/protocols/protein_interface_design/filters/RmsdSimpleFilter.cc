// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/RmsdSimpleFilter.cc
/// @brief rmsd filtering
/// @author Jacob Corn (jecorn@u.washington.edu)
#include <protocols/protein_interface_design/filters/RmsdSimpleFilter.hh>
#include <protocols/protein_interface_design/filters/RmsdSimpleFilterCreator.hh>
#include <protocols/filters/Filter.hh>


// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <core/pose/PDBInfo.hh>

#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <numeric/model_quality/rms.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/NamedAtomID.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <algorithm>
#include <list>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>

//Include Rosetta protocols
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <protocols/toolbox/superimpose.hh>

//Include Boost
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

//Include ObjexxFCL
#include <ObjexxFCL/FArray.all.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

using namespace ObjexxFCL;

RmsdSimpleFilter::RmsdSimpleFilter() :
	protocols::filters::Filter( "RmsdSimple" ),
	threshold_( 5.0 ),
	reference_pose_( /* 0 */ ),
	do_align_( 0 )
{
}

RmsdSimpleFilter::RmsdSimpleFilter(
	core::Real const threshold,
	core::pose::PoseOP reference_pose
) : protocols::filters::Filter( "RmsdSimple" )
{
	threshold_=threshold;
	reference_pose_=reference_pose;
}

RmsdSimpleFilter::~RmsdSimpleFilter() {}

protocols::filters::FilterOP
RmsdSimpleFilter::clone() const {
	return protocols::filters::FilterOP( new RmsdSimpleFilter( *this ) );
}

protocols::filters::FilterOP
RmsdSimpleFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new RmsdSimpleFilter() );
}

static thread_local basic::Tracer TR( "protocols.protein_interface_design.filters.RmsdSimpleFilter" );

core::Real
RmsdSimpleFilter::compute( core::pose::Pose const & pose ) const
{
	core::Real rmsd( 0.0 );
	core::pose::Pose this_target_pose = pose;
	core::pose::Pose this_refe_pose = *reference_pose_;

	TR << boost::format("apply to poses, Target=%s vs residues> ") % (this_target_pose.n_residue(), this_refe_pose.n_residue() )<< std::endl;
	//Apply to a particular chain
	if ( this_target_pose.n_residue() != this_refe_pose.n_residue() ) {
		utility_exit_with_message("Error, reference pose and current pose have different size ");
	}
	if ( b_target_chain_ ) {
		TR <<  "Target chain: " << target_chain_ << ", of a pose with #Chains: " << pose.conformation().num_chains() << std::endl;
		if ( target_chain_ > pose.conformation().num_chains() ) {
			utility_exit_with_message(" FragmentLookupFilter invalid chain" );
		}
		this_target_pose = *this_target_pose.split_by_chain(target_chain_);
		this_refe_pose = *this_refe_pose.split_by_chain(target_chain_);
		TR << boost::format("apply to poses, Target=%s vs residues> ") % (this_target_pose.n_residue(), this_refe_pose.n_residue() )<< std::endl;
		if ( this_target_pose.n_residue() != this_refe_pose.n_residue() ) {
			utility_exit_with_message("Error, reference pose and current pose have different size " );
		}
	}
	utility::vector1< core::Size > positions_to_alignA;
	utility::vector1< core::Size > positions_to_alignB;
	for ( core::Size i=1; i <= this_target_pose.n_residue(); ++i ) {
		positions_to_alignA.push_back(i);
		positions_to_alignB.push_back(i);
	}
	//Do simple distance
	if ( do_align_ == 0 ) {
		rmsd = dist_bb( this_target_pose,
			positions_to_alignA,
			this_refe_pose,
			positions_to_alignB );
	} else if ( do_align_ == 1 ) { //do real RMSD
		rmsd = rmsd_bb( this_target_pose,
			positions_to_alignA,
			this_refe_pose,
			positions_to_alignB );

	} else {
		utility_exit_with_message("Error, incorrect align mode" );
	}
	return rmsd;
}

/**@brief Superposition_transform wrapper (as in Alex Ford code protocols/toolbox/superimpose.[cc,hh])**/
void RmsdSimpleFilter::superposition_transform(
	core::Size natoms,
	ObjexxFCL::FArray1_double const& weights,
	ObjexxFCL::FArray2_double& ref_coords,
	ObjexxFCL::FArray2_double& coords,
	numeric::xyzMatrix< core::Real > &RotM,
	numeric::xyzVector< core::Real > &TvecA,
	numeric::xyzVector< core::Real > &TvecB) const
{
	ObjexxFCL::FArray1D_double ref_transvec( 3 );
	//Get/remove the first COM translation vector
	protocols::toolbox::reset_x( natoms, ref_coords, weights, ref_transvec );
	TvecA = numeric::xyzVector< core::Real >( ref_transvec( 1 ), ref_transvec( 2 ), ref_transvec( 3 ));

	//Get/remove the second COM translation vector
	ObjexxFCL::FArray1D_double transvec( 3 );
	protocols::toolbox::reset_x( natoms, coords, weights, transvec );
	TvecB = numeric::xyzVector< core::Real >( transvec( 1 ), transvec( 2 ), transvec( 3 ));

	//Fit centered coords, returns the rotation matrix, shifted coordinates (to the center)
	protocols::toolbox::fit_centered_coords( natoms, weights, ref_coords, coords, RotM );
}

/**@brief Function that performs alignment of the protein BB on the selected aminoacids.
**Returns the RMSD,
**Uses a reference to the rotation Matrix and Translation Vector,
**Will fail if both poses are not protein <-This can be fixed by adding a list of the atoms to align to the function, but I am not doing it now.
**Will fail if the number of residues to align is not the same in the two poses. **/
core::Real RmsdSimpleFilter::rmsd_bb(
	core::pose::Pose const & poseA,
	utility::vector1< core::Size > const & positions_to_alignA,
	core::pose::Pose const & poseB,
	utility::vector1< core::Size > const & positions_to_alignB) const
{
	core::Size sizeA=positions_to_alignA.size();
	core::Size sizeB=positions_to_alignB.size();
	//Check if the align size vectors are equivalent in size, if not die with error(-1.0)
	if ( sizeA != sizeB ) return -1.0;

	//To store the RMSD
	core::Real RMSD;
	//To store the translation and rotation matrix
	numeric::xyzMatrix< core::Real >  RotM;
	numeric::xyzVector< core::Real >  TvecA;
	numeric::xyzVector< core::Real >  TvecB;
	//To store the positions of the atoms that we want to align
	utility::vector1< numeric::xyzVector< core::Real > > matrixA;
	utility::vector1< numeric::xyzVector< core::Real > > matrixB;
	utility::vector1< bool > v_isNorC_expandedFormat;
	//Copy the BB atoms XYZ positions for pose A and B
	for ( core::Size i = 1; i <= sizeA ; ++i ) {
		matrixA.push_back(poseA.residue(positions_to_alignA[i]).xyz( "CA" ));
		matrixB.push_back(poseB.residue(positions_to_alignB[i]).xyz( "CA" ));
		matrixA.push_back(poseA.residue(positions_to_alignA[i]).xyz( "C" ));
		matrixB.push_back(poseB.residue(positions_to_alignB[i]).xyz( "C" ));
		matrixA.push_back(poseA.residue(positions_to_alignA[i]).xyz( "O" ));
		matrixB.push_back(poseB.residue(positions_to_alignB[i]).xyz( "O" ));
		matrixA.push_back(poseA.residue(positions_to_alignA[i]).xyz( "N" ));
		matrixB.push_back(poseB.residue(positions_to_alignB[i]).xyz( "N" ));
	}
	//Convert the vectors to Fortran style arrays (as in Alex Ford code)
	ObjexxFCL::FArray2D< numeric::Real > FmatrixA( 3, (sizeA*4) );
	ObjexxFCL::FArray2D< numeric::Real > FmatrixB( 3, (sizeA*4) );
	//??Alex: are this post-increments particulary Roseta-coding-standars correct?
	for ( core::Size i = 1; i <= (sizeA*4); i++ ) {
		for ( core::Size j = 1; j <= 3; j++ ) {
			FmatrixA(j,i) = matrixA[i](j);
			FmatrixB(j,i) = matrixB[i](j);
		}
	}
	//Weighted alignment (as in Alex Ford code protocols/seeded_abinitio/util.[cc,hh])
	ObjexxFCL::FArray1D< numeric::Real > weights_fa( (sizeA*4), 1);

	//The actual alignment (as in Alex Ford code protocols/seeded_abinitio/util.[cc,hh])
	superposition_transform((sizeA*4), weights_fa, FmatrixA, FmatrixB, RotM, TvecA, TvecB);

	RMSD=0.0;
	//Storage for d^2
	core::Real tmpDistSqr;
	//Calculate the RMSD
	for ( core::Size i = 1; i <= (sizeA*4); i++ ) {
		tmpDistSqr=0.0;
		for ( int j = 1; j <= 3; j++ ) {
			tmpDistSqr+=(FmatrixA(j,i)-FmatrixB(j,i))*(FmatrixA(j,i)-FmatrixB(j,i));
		}
		RMSD+=tmpDistSqr;
	}
	RMSD=std::sqrt(RMSD/(sizeA*4));
	//Return the RMSD
	return RMSD;
}

/**@brief Returns the BB distance of two poses respect to indexes**/

core::Real RmsdSimpleFilter::dist_bb(
	core::pose::Pose const & poseA,
	utility::vector1< core::Size > const & positions_to_alignA,
	core::pose::Pose const & poseB,
	utility::vector1< core::Size > const & positions_to_alignB) const
{
	core::Size sizeA=positions_to_alignA.size();
	core::Size sizeB=positions_to_alignB.size();
	//Check if the align size vectors are equivalent in size, if not die with error(-1.0)
	if ( sizeA != sizeB ) return -1.0;
	//To store the RMSD
	core::Real RMSD;
	//To store the positions of the atoms that we want to measure
	utility::vector1< numeric::xyzVector< core::Real > > matrixA;
	utility::vector1< numeric::xyzVector< core::Real > > matrixB;
	//Copy the BB atoms XYZ positions for pose A and B
	for ( core::Size i = 1; i <= sizeA ; ++i ) {
		matrixA.push_back(poseA.residue(positions_to_alignA[i]).xyz( "CA" ));
		matrixB.push_back(poseB.residue(positions_to_alignB[i]).xyz( "CA" ));
		matrixA.push_back(poseA.residue(positions_to_alignA[i]).xyz( "C" ));
		matrixB.push_back(poseB.residue(positions_to_alignB[i]).xyz( "C" ));
		matrixA.push_back(poseA.residue(positions_to_alignA[i]).xyz( "O" ));
		matrixB.push_back(poseB.residue(positions_to_alignB[i]).xyz( "O" ));
		matrixA.push_back(poseA.residue(positions_to_alignA[i]).xyz( "N" ));
		matrixB.push_back(poseB.residue(positions_to_alignB[i]).xyz( "N" ));
	}
	//Convert the vectors to Fortran style arrays (as in Alex Ford code)
	ObjexxFCL::FArray2D< numeric::Real > FmatrixA( 3, (sizeA*4) );
	ObjexxFCL::FArray2D< numeric::Real > FmatrixB( 3, (sizeA*4) );
	//??Alex: are this post-increments particulary Roseta-coding-standars correct?
	for ( core::Size i = 1; i <= (sizeA*4); i++ ) {
		for ( core::Size j = 1; j <= 3; j++ ) {
			FmatrixA(j,i) = matrixA[i](j);
			FmatrixB(j,i) = matrixB[i](j);
		}
	}

	RMSD=0.0;
	core::Real tmpDistSqr=0.0;
	//Calculate the RMSD
	for ( core::Size i = 1; i <= (sizeA*4); i++ ) {
		//Storage for d^2
		tmpDistSqr=0.0;
		for ( int j = 1; j <= 3; j++ ) {
			tmpDistSqr+=(FmatrixA(j,i)-FmatrixB(j,i))*(FmatrixA(j,i)-FmatrixB(j,i));
		}
		RMSD+=tmpDistSqr;
	}
	RMSD=std::sqrt(RMSD/(sizeA*4));
	//Return the RMSD
	return RMSD;
}


//Shared Filter Stuff
bool RmsdSimpleFilter::apply( core::pose::Pose const & pose ) const {

	TR << "Calculating BB RMSD." ;
	core::Real const rmsd( compute( pose ));
	if ( rmsd <= threshold_ ) {
		TR<<" RMSD Filter Pass: "<< rmsd <<std::endl;
		return( true );
	} else {
		TR<<" RMSD Filter Fail:"<< rmsd << std::endl;
		return( false );
	}
}

//Common Filter Stuff
void RmsdSimpleFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ));
	out << "RMSD: " << rmsd <<'\n';
	//TR<<"RMSD: " << rmsd <<'\n';
}

core::Real RmsdSimpleFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ));
	return( (core::Real) rmsd );
}

void RmsdSimpleFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data_map,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose )
{
	/// @brief
	reference_pose_ = core::pose::PoseOP( new core::pose::Pose( pose ) );
	target_chain_ = 0;
	threshold_ = 0.0;

	if ( tag->hasOption("reference_name") ) {
		reference_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag, data_map );
	} else if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		core::import_pose::pose_from_pdb( *reference_pose_, basic::options::option[ basic::options::OptionKeys::in::file::native ] );
	} else {
		utility_exit_with_message("Not reference structure defined! Use [reference_name] or [-in::file::native] ");
	}
	if ( tag->hasOption("chain") ) {
		target_chain_ = tag->getOption<core::Size>( "chain" );
		b_target_chain_ = true;
	} else {
		TR.Warning << "No chain specified, using them all for calculation" << std::endl;
	}

	if ( tag->hasOption("threshold") ) {
		threshold_ = tag->getOption<core::Real>( "threshold");
	} else {
		TR.Warning << "Threshold option not specified, using default" << threshold_ << std::endl;
	}
	if ( tag->hasOption("align") ) {
		do_align_ = tag->getOption<core::Size>( "align");
		if ( do_align_ == 0 ) {
			TR.Info << "Align mode is disabled, therefore the result will be the euclidean distance" << std::endl;
		} else if ( do_align_ == 1 ) {
			TR.Info << "Align mode is enabled, therefore the result will be RMSD" << std::endl;
		} else {
			utility_exit_with_message("Align option is boolean, allowed values are 0 and 1 (deactivated/activated)");
		}
	} else {
		TR.Warning << "Align option not specified, using default (activated)" << threshold_ << std::endl;
	}
}

protocols::filters::FilterOP RmsdSimpleFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new RmsdSimpleFilter );
}

std::string RmsdSimpleFilterCreator::keyname() const {
	return "RmsdSimple";
}


} // filters
} // protein_interface_design
} // devel
