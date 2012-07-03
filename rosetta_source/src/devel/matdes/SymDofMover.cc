// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file 
/// @brief 
/// @author Jacob Bale ( balej@uw.edu )

// Unit headers
#include <devel/matdes/SymDofMover.hh>
#include <devel/matdes/SymDofMoverCreator.hh>

// project headers
#include <protocols/moves/Mover.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/VirtualCoordinate.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <core/kinematics/Jump.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

static basic::Tracer TR("devel.matdes.SymDofMover");

namespace devel {
namespace matdes {

using namespace core;
using namespace utility;
using core::pose::Pose;

typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;

// -------------  Mover Creator -------------
std::string
SymDofMoverCreator::keyname() const
{
	return SymDofMoverCreator::mover_name();
}

protocols::moves::MoverOP
SymDofMoverCreator::create_mover() const {
	return new SymDofMover;
}

std::string
SymDofMoverCreator::mover_name()
{
	return "SymDofMover";
}
// -------------  Mover Creator -------------

SymDofMover::SymDofMover() :
	symm_file_(),
	sym_dof_names_(),
	radial_disps_(),
	angles_()
{ }

SymDofMover::SymDofMover(const SymDofMover& rval) :
	symm_file_( rval.symm_file_ ),
	sym_dof_names_( rval.sym_dof_names_ ),
	radial_disps_( rval.radial_disps_ ),
	angles_( rval.angles_ )
{ }

protocols::moves::MoverOP 
SymDofMover::clone() const {
	return new SymDofMover( *this );
}

protocols::moves::MoverOP 
SymDofMover::fresh_instance() const {
				return new SymDofMover();
}

utility::vector1<std::string>
SymDofMover::get_sym_dof_names() { 
	return sym_dof_names_;
}

utility::vector1<Real>
SymDofMover::get_angles() { 
	return angles_;
}

utility::vector1<Real>
SymDofMover::get_radial_disps() { 
	return radial_disps_;
}

// pose manipulation. Consider moving to util.cc ? //

void SymDofMover::trans_pose( Pose & pose, Vec const & trans, Size start, Size end ) {
	for(Size ir = start; ir <= end; ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, pose.xyz(aid) + (Vec)trans );
		}
	}
}
void SymDofMover::rot_pose( Pose & pose, Mat const & rot, Size start, Size end ) {
	for(Size ir = start; ir <= end; ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, rot * pose.xyz(aid) );
		}
	}
}
void SymDofMover::rot_pose( Pose & pose, Mat const & rot, Vec const & cen, Size start, Size end ) {
	trans_pose(pose,-cen,start,end);
	rot_pose(pose,rot,start,end);
	trans_pose(pose,cen,start,end);
}
void SymDofMover::rot_pose( Pose & pose, Vec const & axis, double const & ang, Size start, Size end ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang),start,end);
}
void SymDofMover::rot_pose( Pose & pose, Vec const & axis, double const & ang, Vec const & cen, Size start, Size end ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang),cen,start,end);
}
void SymDofMover::alignaxis(core::pose::Pose & pose, Vec newaxis, Vec oldaxis, Vec cen, Size start, Size end ) {
	newaxis.normalize();
	oldaxis.normalize();
	if (newaxis.normalize() != oldaxis.normalize()) {
		Vec axis = newaxis.cross(oldaxis).normalized();
		Real ang = -acos(numeric::max(-1.0,numeric::min(1.0,newaxis.dot(oldaxis))))*180/numeric::constants::d::pi;
		rot_pose(pose,axis,ang,cen,start,end);
	}
}

// Do stuff //

void
SymDofMover::apply(Pose & pose) {
	using core::pose::Pose;

// Read in user info //

	utility::vector1<Real> radial_disps = get_radial_disps();
	utility::vector1<Real> angles = get_angles();
	utility::vector1< std::string > sym_dof_names = get_sym_dof_names();
	std::string symm_file = symm_file_;

// Read in symmetry info from symmetry definition file //
	
	core::conformation::symmetry::SymmData symmdata( pose.n_residue(), pose.num_jump() );
	symmdata.read_symmetry_data_from_file(symm_file);
	std::map< std::string, core::conformation::symmetry::VirtualCoordinate > coords = symmdata.get_virtual_coordinates();
	std::map< std::string, std::pair< std::string, std::string > > virt_connects = symmdata.get_virtual_connects();

	for (int i = 1; i <= sym_dof_names.size(); i++) {

		core::Size sub_start= pose.conformation().chain_begin(i);
		core::Size sub_end= pose.conformation().chain_end(i);
	
// translate each subunit along the z-axis by user defined values //
		trans_pose(pose,Vec(0,0,radial_disps[i]),sub_start,sub_end);

// rotate each subunit along the z-axis by user defined values //
		rot_pose(pose,Vec(0,0,1),angles[i],sub_start,sub_end);

// read in the axes for each subunit //
		std::string tag (virt_connects.find( sym_dof_names[i])->second.first );
		TR << sym_dof_names[i] << tag << std::endl;
		conformation::symmetry::VirtualCoordinate virt_coord( coords.find( tag )->second );

// align the z-axis of each subunit with the appropriate axis of the symdof jump from the symmetry definition file //
		alignaxis(pose,virt_coord.get_x(),Vec(0,0,1),Vec(0,0,0),sub_start,sub_end);
	}

// symmetrize pose //

	core::pose::symmetry::make_symmetric_pose(pose, symmdata);

}

void 
SymDofMover::parse_my_tag( TagPtr const tag,
										 DataMap & data,
										 Filters_map const &,
										 Movers_map const &,
										 Pose const & ) {

	// Turn symmetry hacks on
	basic::options::option[basic::options::OptionKeys::symmetry::symmetry_definition].value( "dummy" );

	using std::string;
	symm_file_ = tag->getOption<string>( "symm_file" );
	sym_dof_names_ = utility::string_split( tag->getOption< std::string >( "sym_dof_names" ), ',' );
	utility::vector1< std::string > radial_disps_strings = utility::string_split( tag->getOption< std::string >( "radial_disps" ), ',' );
	utility::vector1<Real> radial_disps;
	Real real_disp;
	for(Size i = 1; i <= radial_disps_strings.size(); i++) {
		real_disp = std::atof( radial_disps_strings[i].c_str() );
		radial_disps.push_back( real_disp );
	}
	radial_disps_ = radial_disps;
	utility::vector1< std::string > angles_strings = utility::string_split( tag->getOption< std::string >( "angles" ), ',' );
	utility::vector1<Real> angles;
	Real real_angle;
	for(Size i = 1; i <= angles_strings.size(); i++) {
		real_angle = std::atof( angles_strings[i].c_str() );
		angles.push_back( real_angle );
	}
	angles_ = angles;
	}

} // matdes
} // devel
