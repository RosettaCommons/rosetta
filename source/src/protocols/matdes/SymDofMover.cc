// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file
/// @brief
/// @author Jacob Bale ( balej@uw.edu )

// Unit headers
#include <protocols/matdes/SymDofMover.hh>
#include <protocols/matdes/SymDofMoverCreator.hh>

// Package headers
#include <protocols/matdes/SymDofMoverSampler.hh>

// project headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/VirtualCoordinate.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <core/kinematics/Jump.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.matdes.SymDofMover" );


namespace protocols {
namespace matdes {

using namespace core;
using namespace utility;
using core::pose::Pose;

using Vec = numeric::xyzVector<Real>;
using Mat = numeric::xyzMatrix<Real>;

// -------------  Mover Creator -------------


// -------------  Mover Creator -------------

SymDofMover::SymDofMover() :
	set_sampler_(true),
	auto_range_(false),
	sampling_mode_("single_dock"),
	symm_file_(""),
	translation_axes_(),
	rotation_axes_(),
	flip_input_about_axes_(),
	align_input_axes_to_symdof_axes_(),
	sym_dof_names_(),
	radial_disps_(),
	angles_(),
	radial_offsets_(),
	radial_disps_range_min_(),
	radial_disps_range_max_(),
	angles_range_min_(),
	angles_range_max_(),
	radial_disp_steps_(),
	angle_steps_(),
	radial_disp_deltas_(),
	angle_deltas_()
{ }

protocols::moves::MoverOP
SymDofMover::clone() const {
	return utility::pointer::make_shared< SymDofMover >( *this );
}

protocols::moves::MoverOP
SymDofMover::fresh_instance() const {
	return utility::pointer::make_shared< SymDofMover >();
}

utility::vector1<std::string>
SymDofMover::get_sym_dof_names() {
	return sym_dof_names_;
}

/// @details WARNING WARNING WARNING THIS IS A THREAD-UNSAFE FUNCTION SINCE IT USES THE
/// SymDofMoverSampler, A NON-CONSTANT SINGLETON.
utility::vector1<Real>
SymDofMover::get_angles() {
	utility::vector1< std::string > sym_dof_names = get_sym_dof_names();
	utility::vector1<Real> angles;
	if ( sampling_mode_ == "grid" ) {
		utility::vector1<Real> angle_diffs = SymDofMoverSampler::get_instance()->get_angle_diffs();
		for ( core::Size i = 1; i <= sym_dof_names.size(); i++ ) {
			angles.push_back(angles_[i] + angle_diffs[i]);
		}
	} else {
		for ( core::Size i = 1; i <= sym_dof_names.size(); i++ ) {
			if ( sampling_mode_ == "uniform" ) {
				angles.push_back(angles_[i] + angles_range_min_[i] + ( angles_range_max_[i] - angles_range_min_[i]) * numeric::random::rg().uniform());
			} else if ( sampling_mode_ == "gaussian" ) {
				angles.push_back(angles_[i] + angle_deltas_[i] * numeric::random::rg().gaussian());
			} else {
				angles.push_back(angles_[i]);
			}
		}
	}
	if ( set_sampler_ ) {
		SymDofMoverSampler::get_instance()->set_angles(angles);
	}
	return angles;
}

/// @details WARNING WARNING WARNING THIS IS A THREAD-UNSAFE FUNCTION SINCE IT USES THE
/// SymDofMoverSampler, A NON-CONSTANT SINGLETON.
utility::vector1<Real>
SymDofMover::get_radial_disps() {
	utility::vector1< std::string > sym_dof_names = get_sym_dof_names();
	utility::vector1<Real> radial_disps;
	if ( sampling_mode_ == "grid" ) {
		utility::vector1<Real> radial_diffs = SymDofMoverSampler::get_instance()->get_radial_disp_diffs();
		for ( core::Size i = 1; i <= sym_dof_names.size(); i++ ) {
			radial_disps.push_back(radial_disps_[i] + radial_offsets_[i] + radial_diffs[i]);
		}
	} else {
		for ( core::Size i = 1; i <= sym_dof_names.size(); i++ ) {
			if ( sampling_mode_ == "uniform" ) {
				radial_disps.push_back(radial_disps_[i] + radial_offsets_[i] + radial_disps_range_min_[i] + ( radial_disps_range_max_[i] - radial_disps_range_min_[i]) * numeric::random::rg().uniform());
			} else if ( sampling_mode_ == "gaussian" ) {
				radial_disps.push_back(radial_disps_[i] + radial_offsets_[i] + radial_disp_deltas_[i] * numeric::random::rg().gaussian());
			} else {
				radial_disps.push_back(radial_disps_[i] + radial_offsets_[i]);
			}
		}
	}
	if ( set_sampler_ ) {
		SymDofMoverSampler::get_instance()->set_radial_disps(radial_disps);
	}
	return radial_disps;
}

// pose manipulation. Consider moving to util.cc ? //

static void trans_pose( Pose & pose, Vec const & trans, core::Size start, core::Size end ) {
	for ( core::Size ir = start; ir <= end; ++ir ) {
		for ( core::Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia ) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, pose.xyz(aid) + (Vec)trans );
		}
	}
}
static void rot_pose( Pose & pose, Mat const & rot, core::Size start, core::Size end ) {
	for ( core::Size ir = start; ir <= end; ++ir ) {
		for ( core::Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia ) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, rot * pose.xyz(aid) );
		}
	}
}
static void rot_pose( Pose & pose, Mat const & rot, Vec const & cen, core::Size start, core::Size end ) {
	trans_pose(pose,-cen,start,end);
	rot_pose(pose,rot,start,end);
	trans_pose(pose,cen,start,end);
}
static void rot_pose( Pose & pose, Vec const & axis, core::Real const & ang, core::Size start, core::Size end ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang),start,end);
}
static void rot_pose( Pose & pose, Vec const & axis, core::Real const & ang, Vec const & cen, core::Size start, core::Size end ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang),cen,start,end);
}
static void alignaxis(core::pose::Pose & pose, Vec newaxis, Vec oldaxis, Vec cen, core::Size start, core::Size end ) {
	newaxis.normalize();
	oldaxis.normalize();
	if ( newaxis.normalize() != oldaxis.normalize() ) {
		Vec axis = newaxis.cross(oldaxis).normalized();
		Real ang = -acos(numeric::max(-1.0,numeric::min(1.0,newaxis.dot(oldaxis))))*180/numeric::constants::d::pi;
		rot_pose(pose,axis,ang,cen,start,end);
	}
}

// Do stuff //

numeric::Real get_intra_contacts(Pose const & pose){
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;
	using numeric::Real;
	utility::vector1<std::string> symcomp = symmetric_components(pose);
	runtime_assert( symcomp.size() == 2 );
	std::string c1 = symcomp[1];
	std::string c2 = symcomp[2];
	core::Size beg1=get_component_lower_bound(pose,c1);
	core::Size beg2=get_component_lower_bound(pose,c2);
	core::Size end1=get_component_upper_bound(pose,c1);
	core::Size end2=get_component_upper_bound(pose,c2);
	Real ncontact = 0.0;
	for ( core::Size ir = beg1; ir <= end1; ++ir ) {
		for ( core::Size jr = beg2; jr <= end2; ++jr ) {
			ncontact += ( 49.0 >= pose.residue(ir).nbr_atom_xyz().distance_squared( pose.residue(jr).nbr_atom_xyz() ) );
		}
	}
	return ncontact;
}
void maximize_sub1_contact(Pose & pose, core::Size const & nf1, core::Size const & nf2, Vec const & ax1, Vec const & ax2 ){
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;
	using numeric::Real;
	utility::vector1<std::string> symcomp = symmetric_components(pose);
	runtime_assert( symcomp.size() == 2 );
	std::string c1 = symcomp[1];
	std::string c2 = symcomp[2];
	core::Size beg1=get_component_lower_bound(pose,c1);
	core::Size beg2=get_component_lower_bound(pose,c2);
	core::Size end1=get_component_upper_bound(pose,c1);
	core::Size end2=get_component_upper_bound(pose,c2);
	Real mx=-9e9;
	core::Size imx=0,jmx=0;
	for ( core::Size i = 0; i < nf1; ++i ) {
		for ( core::Size j = 0; j < nf2; ++j ) {
			Real ncontact = get_intra_contacts(pose);
			if ( ncontact>mx ) {
				mx = ncontact;
				imx = i;
				jmx = j;
			}
			// std::cout << i << " " << j << " " << ncontact << std::endl;
			rot_pose(pose,ax2,360.0/(Real)nf2,beg2,end2);
		}
		rot_pose(pose,ax1,360.0/(Real)nf1,beg1,end1);
	}
	rot_pose(pose,ax1,360.0/(Real)nf1 * (Real)imx,beg1,end1);
	rot_pose(pose,ax2,360.0/(Real)nf2 * (Real)jmx,beg2,end2);
}

void
SymDofMover::add_components_to_pose_if_necessary(Pose & pose){
	using namespace basic::options;
	TR << "checking for additional components" << std::endl;
	if ( option[OptionKeys::in::file::t].user() ) {
		runtime_assert_msg(option[OptionKeys::in::file::t]().size() == 1,
			"SymDofMover must have one or no inputs in -t");
		core::pose::PoseCOP b = core::import_pose::pose_from_file( option[OptionKeys::in::file::t]().front() , core::import_pose::PDB_file);
		core::Size nres1 = pose.size();
		core::pose::append_pose_to_pose( pose, *b, true );
		for ( core::uint ir =         1; ir <= nres1           ; ++ir ) {
			pose.pdb_info()->chain(static_cast<long int>(ir), "A");
		}
		for ( core::uint ir = nres1 + 1; ir <= pose.size(); ++ir ) {
			pose.pdb_info()->chain(static_cast<long int>(ir), "B");
		}
		pose.update_pose_chains_from_pdb_chains();
		// for(int ir = 1; ir <= pose.size(); ++ir) std::cout << ir << " " << pose.chain(ir) << std::endl;
		// pose.dump_pdb("combined.pdb");
		TR << "added to pose: " << option[OptionKeys::in::file::t]().front() << " num chains: " << pose.conformation().num_chains() << std::endl;
	}
}

/// @details WARNING WARNING WARNING THIS IS A THREAD-UNSAFE FUNCTION SINCE IT USES THE
/// SymDofMoverSampler, A NON-CONSTANT SINGLETON.
void
SymDofMover::apply(Pose & pose) {
	using core::pose::Pose;
	using namespace core::pose::symmetry;

	add_components_to_pose_if_necessary(pose);
	// Read in user info //

	utility::vector1< std::string > sym_dof_names = get_sym_dof_names();
	utility::vector1<Real> radial_disps = get_radial_disps();
	utility::vector1<Real> angles = get_angles();
	core::Size sym_aware_jump_id;
	Vec translation;
	Mat rotation;

	// If already symmetric pose, then apply the displacements and/or rotations to the specified symdofs //
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		for ( core::Size i = 1; i <= sym_dof_names.size(); i++ ) {
			sym_aware_jump_id = core::pose::symmetry::sym_dof_jump_num( pose, sym_dof_names[i] );
			core::kinematics::Jump j = pose.jump(sym_aware_jump_id);
			const Vec init_trans = pose.jump(sym_aware_jump_id).get_translation();
			const Mat init_rot = pose.jump(sym_aware_jump_id).get_rotation();
			/*   core::conformation::symmetry::SymDof dof;
			core::conformation::symmetry::SymmetricConformation & symm_conf( dynamic_cast<core::conformation::symmetry::SymmetricConformation & > ( split_pose.conformation()) );
			if( (radial_disps.size() > 0 && translation_axes_.size() == 0) || (angles.size() > 0 && rotation_axes_.size() == 0) ) {
			for(std::map< core::Size, core::conformation::symmetry::SymDof >::const_iterator it = symm_conf.Symmetry_Info()->get_dofs().begin(); it != symm_conf.Symmetry_Info()->get_dofs().end(); ++it){
			if( it->first == sym_aware_jump_id ){
			dof = it;
			}
			}
			if( translation_axes_.size() == 0 ) {
			for ( core::Size i=1; i<=3; ++i ) {
			if( dof.allow_dof(i) ) {
			if ( i == 1 ) translation_axes_.push_back("x");
			if ( i == 2 ) translation_axes_.push_back("y");
			if ( i == 3 ) translation_axes_.push_back("z");
			}
			}
			}
			if( rotation_axes_.size() == 0 ) {
			for ( core::Size i=4; i<=6; ++i ) {
			if( dof.allow_dof(i) ) {
			if ( i == 4 ) rotation_axes_.push_back("x");
			if ( i == 5 ) rotation_axes_.push_back("y");
			if ( i == 6 ) rotation_axes_.push_back("z");
			}
			}
			}
			*/

			if ( translation_axes_.size() > 0 ) {
				if ( translation_axes_[i] == "x" ) translation = Vec(radial_disps[i],0,0) + init_trans;
				else if ( translation_axes_[i] == "y" ) translation = Vec(0, radial_disps[i], 0) + init_trans;
				else if ( translation_axes_[i] == "z" ) translation = Vec(0,0, radial_disps[i]) + init_trans;
				else utility_exit_with_message(translation_axes_[i] + " is not a valid axis (x,y or z). Use lower case");
				TR << "displacing " << radial_disps[i] << " angstroms along symdof " << sym_dof_names[i] << std::endl;
				j.set_translation( translation );
			}
			if ( rotation_axes_.size() > 0 ) {
				if ( rotation_axes_[i] == "x" ) rotation = Mat(numeric::x_rotation_matrix_degrees( angles[i] ) * init_rot);
				else if ( rotation_axes_[i] == "y" ) rotation = Mat(numeric::y_rotation_matrix_degrees( angles[i] )* init_rot);
				else if ( rotation_axes_[i] == "z" ) rotation = Mat(numeric::z_rotation_matrix_degrees( angles[i] ) * init_rot);
				else utility_exit_with_message(translation_axes_[i] + " is not a valid axis (x,y or z). Use lower case");
				TR << "Rotating " << angles[i] << " degrees about symdof " << sym_dof_names[i] << std::endl;
				j.set_rotation( rotation );
			}
			pose.set_jump(sym_aware_jump_id,j);
		}

		// If not symmetric, then apply the displacements and/or rotations to the input subunits along the user-specified axes and then align those axes with the corresponding axes for each symdof in the symmetry definition file.  //
	} else {

		// Read in user info //

		std::string symm_file = symm_file_;
		//std::string symm_name = utility::file_basename(symm_file);

		// Read in symmetry info from symmetry definition file //

		core::conformation::symmetry::SymmData symmdata( pose.size(), pose.num_jump() );
		symmdata.read_symmetry_data_from_file(symm_file);
		std::map< std::string, core::conformation::symmetry::VirtualCoordinate > coords = symmdata.get_virtual_coordinates();
		std::map< std::string, std::pair< std::string, std::string > > virt_connects = symmdata.get_virtual_connects();

		if ( symm_file_[0]=='P' && symm_file_[2]=='_' && symm_file_[5]=='.' ) {
			TR << "Doing DOF placement for 2D lattice (based on sym file name P?_??.*" << symm_file_ << std::endl;
			for ( core::Size i = 1; i <= sym_dof_names.size(); i++ ) {

				core::Size sub_start= pose.conformation().chain_begin(i);
				core::Size sub_end= pose.conformation().chain_end(i);

				// Rotate subtunit 180 degrees about the specified axis. ie, "reverse" the component before further manipulation.
				if ( flip_input_about_axes_.size() > 0 ) {
					if ( flip_input_about_axes_[i] == "z" ) rot_pose(pose,Vec(0,0,1),180,sub_start,sub_end);
					else if ( flip_input_about_axes_[i] == "x" ) rot_pose(pose,Vec(1,0,0),180,sub_start,sub_end);
					else if ( flip_input_about_axes_[i] == "y" ) rot_pose(pose,Vec(0,1,0),180,sub_start,sub_end);
				}

				// rotate each subunit along the specified axis by user defined values //
				TR << "Rotating component " << i << " " << angles[i] << " degrees" << std::endl;
				rot_pose(pose,Vec(0,0,1),angles[i],sub_start,sub_end);

				// translate each subunit along the specified axis by user defined values //
				TR << "Translating component " << i << " " << radial_disps[i] << " angstroms" << std::endl;
				if ( 1==i ) trans_pose(pose,Vec(0,0,radial_disps[i]) ,sub_start,sub_end);
				if ( 2==i ) trans_pose(pose,Vec(-radial_disps[i],0,0),sub_start,sub_end);

			}
			// rot   1  z
			// rot   2  z
			// trans 1  z
			// trans 2 -x
			// set x jump 2 (below)

		} else {

			for ( core::Size i = 1; i <= sym_dof_names.size(); i++ ) {

				core::Size sub_start= pose.conformation().chain_begin(i);
				core::Size sub_end= pose.conformation().chain_end(i);

				// Rotate subtunit 180 degrees about the specified axis. ie, "reverse" the component before further manipulation.
				if ( flip_input_about_axes_.size() > 0 ) {
					TR << "Fliping component " << i << " about the " << flip_input_about_axes_[i] << " axis" << std::endl;
					if ( flip_input_about_axes_[i] == "z" ) rot_pose(pose,Vec(0,0,1),180,sub_start,sub_end);
					else if ( flip_input_about_axes_[i] == "x" ) rot_pose(pose,Vec(1,0,0),180,sub_start,sub_end);
					else if ( flip_input_about_axes_[i] == "y" ) rot_pose(pose,Vec(0,1,0),180,sub_start,sub_end);
				}

				// rotate each subunit along the specified axis by user defined values //
				if ( rotation_axes_.size() > 0 ) {
					TR << "Rotating component " << i << " " << angles[i] << " degrees about the " << rotation_axes_[i] << " axis" << std::endl;
					if ( rotation_axes_[i] == "z" ) rot_pose(pose,Vec(0,0,1),angles[i],sub_start,sub_end);
					else if ( rotation_axes_[i] == "x" ) rot_pose(pose,Vec(1,0,0),angles[i],sub_start,sub_end);
					else if ( rotation_axes_[i] == "y" ) rot_pose(pose,Vec(0,1,0),angles[i],sub_start,sub_end);
					else utility_exit_with_message("Specified rotation axis does not match with either x, y, or z");
				}
				// translate each subunit along the specified axis by user defined values //
				if ( translation_axes_.size() > 0 ) {
					TR << "Translating component " << i << " " << radial_disps[i] << " angstroms along the " << translation_axes_[i] << " axis" << std::endl;
					if ( translation_axes_[i] == "z" ) trans_pose(pose,Vec(0,0,radial_disps[i]),sub_start,sub_end);
					else if ( translation_axes_[i] == "x" ) trans_pose(pose,Vec(radial_disps[i],0,0),sub_start,sub_end);
					else if ( translation_axes_[i] == "y" ) trans_pose(pose,Vec(0,radial_disps[i],0),sub_start,sub_end);
					else utility_exit_with_message("Specified translation axis does not match with either x, y, or z");
				}


				// read in the axes for each subunit //

				std::string tag (virt_connects.find( sym_dof_names[i])->second.first );

				TR << sym_dof_names[i] << " -> " << tag << std::endl;

				conformation::symmetry::VirtualCoordinate virt_coord( coords.find( tag )->second );

				// align the specified axis of each subunit with the appropriate axis of the symdof jump from the symmetry definition file //
				if ( align_input_axes_to_symdof_axes_.size() > 0 ) {

					TR << "aligned_axis" << align_input_axes_to_symdof_axes_[i] << " with symdof axis " << sym_dof_names[i] << std::endl;
					if ( align_input_axes_to_symdof_axes_[i] == "z" ) alignaxis(pose,virt_coord.get_x(),Vec(0,0,1),Vec(0,0,0),sub_start,sub_end);
					else if ( align_input_axes_to_symdof_axes_[i] == "x" ) alignaxis(pose,virt_coord.get_x(),Vec(1,0,0),Vec(0,0,0),sub_start,sub_end);
					else if ( align_input_axes_to_symdof_axes_[i] == "y" ) alignaxis(pose,virt_coord.get_x(),Vec(0,1,0),Vec(0,0,0),sub_start,sub_end);
					else utility_exit_with_message("Specified input axis for alignment does not match with either x, y, or z");
				}

			}

		}


		// symmetrize pose //

		core::pose::symmetry::make_symmetric_pose(pose, symmdata);
		std::string symname = utility::file_basename(symm_file_);
		TR << symname[0] << " " << symname[2] << " " << symname[5] << std::endl;
		if ( symname[0]=='P' && symname[2]=='_' && symname[5]=='.' ) {
			TR << "try to setup P6 sym" << std::endl;
			using namespace core::conformation::symmetry;
			auto const & symm_conf ( dynamic_cast<SymmetricConformation const & > ( pose.conformation() ) );
			SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
			std::map< core::Size, SymDof > dofs ( symm_info->get_dofs() );
			std::map< core::Size, SymDof >::iterator it;
			auto it_begin = dofs.begin();
			auto it_end = dofs.end();
			for ( it=it_begin; it != it_end; ++it ) {
				core::kinematics::Jump j(pose.jump(it->first));
				Vec t = j.get_translation();
				if ( t.length() > 0.00001 ) {
					TR << "try to set P6 translation" << std::endl;
					Vec newt = Vec(0,-radial_disps[2],0);
					j.set_translation(newt);
					Vec origt = pose.xyz(core::id::AtomID(1,pose.conformation().chain_begin(2)));
					pose.set_jump(it->first,j);
					Vec delta = origt - pose.xyz(core::id::AtomID(1,pose.conformation().chain_begin(2)));
					trans_pose(pose,delta,pose.conformation().chain_begin(2),pose.conformation().chain_end(2));
					rot_pose(pose,Vec(0,0,1),180.0,Vec(0,0,0),1,symm_info->num_independent_residues()); // AAAHH!!
					// cout << "trans jump! " << j << " " << it->first << " " << dcmp1 << " " << dcmp2 << " " << newt-t << " " << delta << endl;
				}
			}
		}
	}
	if ( sampling_mode_ == "grid" ) {
		SymDofMoverSampler::get_instance()->step();
	}

}

// Parse info from xml //
void
SymDofMover::parse_my_tag( TagCOP tag,
	DataMap &
) {

	// Turn symmetry hacks on
	if ( !basic::options::option[basic::options::OptionKeys::symmetry::symmetry_definition].user() ) {
		basic::options::option[basic::options::OptionKeys::symmetry::symmetry_definition].value( "dummy" );
	}

	using std::string;
	set_sampler_ = tag->getOption<bool>( "set_sampler", true );
	auto_range_ = tag->getOption<bool>( "auto_range", false );
	symm_file_ = tag->getOption<string>( "symm_file", "");
	sym_dof_names_ = utility::string_split( tag->getOption< std::string >( "sym_dof_names","" ), ',' );
	SymDofMoverSampler::get_instance()->set_sym_dof_names(sym_dof_names_);
	sampling_mode_ = tag->getOption<std::string>("sampling_mode", "single_dock");
	if ( tag->hasOption( "flip_input_about_axes") ) {
		flip_input_about_axes_ = utility::string_split( tag->getOption<std::string>( "flip_input_about_axes"), ',' );
	}
	if ( tag->hasOption( "align_input_axes_to_symdof_axes") ) {
		align_input_axes_to_symdof_axes_ = utility::string_split(tag->getOption<std::string>( "align_input_axes_to_symdof_axes"), ',' );
	}
	if ( tag->hasOption( "radial_disps") || tag->hasOption("radial_offsets") ) {
		translation_axes_ = utility::string_split( tag->getOption<std::string>( "translation_axes"), ',' );
	}

	utility::vector1<Real> radial_disps;
	if ( tag->hasOption( "radial_disps") ) {
		utility::vector1< std::string > radial_disp_strings = utility::string_split( tag->getOption< std::string >( "radial_disps" ), ',' );
		Real real_disp;
		for ( core::Size i = 1; i <= radial_disp_strings.size(); i++ ) {
			TR << radial_disp_strings[i] << std::endl;
			real_disp = std::atof( radial_disp_strings[i].c_str() );
			radial_disps.push_back( real_disp );
		}
		if ( translation_axes_.size() != radial_disps.size() ) {
			utility_exit_with_message("The number of radial disps does not match the number of translation_axes.");
		}
	} else {
		for ( core::Size i = 1; i <= sym_dof_names_.size(); i++ ) {
			radial_disps.push_back( 0 );
		}
	}
	radial_disps_ = radial_disps;

	utility::vector1<Real> angles;
	if ( tag->hasOption( "angles") ) {
		rotation_axes_ = utility::string_split( tag->getOption<std::string>( "rotation_axes"), ',' );
		utility::vector1< std::string > angle_strings = utility::string_split( tag->getOption< std::string >( "angles" ), ',' );
		Real real_angle;
		for ( core::Size i = 1; i <= angle_strings.size(); i++ ) {
			real_angle = std::atof( angle_strings[i].c_str() );
			angles.push_back( real_angle );
		}
		if ( rotation_axes_.size() != angles.size() ) {
			utility_exit_with_message("The number of angles of rotation does not match the number of rotation_axes.");
		}
	} else {
		for ( core::Size i = 1; i <= sym_dof_names_.size(); i++ ) {
			angles.push_back( 0 );
		}
	}
	angles_ = angles;

	utility::vector1<Real> radial_offsets;
	if ( tag->hasOption("radial_offsets") ) {
		utility::vector1< std::string > radial_offset_strings = utility::string_split( tag->getOption< std::string >( "radial_offsets" ), ',' );
		if ( tag->hasOption("radial_disps") && (radial_offset_strings.size() != radial_disps_.size()) ) {
			utility_exit_with_message("The number of radial offsets does not match the number of radial disps.");
		} else {
			Real real_offset;
			for ( core::Size i = 1; i <= radial_offset_strings.size(); i++ ) {
				real_offset = std::atof( radial_offset_strings[i].c_str() );
				if ( auto_range_ && (radial_disps_[i] < 0) ) {
					radial_offsets.push_back( -real_offset );
				} else {
					radial_offsets.push_back( real_offset );
				}
			}
		}
		if ( translation_axes_.size() != radial_offsets.size() ) {
			utility_exit_with_message("The number of radial offsets does not match the number of translation_axes.");
		}
	} else {
		for ( core::Size i = 1; i <= sym_dof_names_.size(); i++ ) {
			radial_offsets.push_back( 0 );
		}
	}
	radial_offsets_ = radial_offsets;

	if ( sampling_mode_ == "grid" || sampling_mode_ == "uniform" ) {

		utility::vector1< std::string > radial_disp_range_min_strings = utility::string_split( tag->getOption< std::string >( "radial_disps_range_min" ), ',' );
		utility::vector1<Real> radial_disps_range_min;
		Real real_disps_range_min;
		for ( core::Size i = 1; i <= radial_disp_range_min_strings.size(); i++ ) {
			real_disps_range_min = std::atof( radial_disp_range_min_strings[i].c_str() );
			radial_disps_range_min.push_back( real_disps_range_min );
		}
		radial_disps_range_min_ = radial_disps_range_min;

		utility::vector1< std::string > radial_disp_range_max_strings = utility::string_split( tag->getOption< std::string >( "radial_disps_range_max" ), ',' );
		utility::vector1<Real> radial_disps_range_max;
		Real real_disps_range_max;
		for ( core::Size i = 1; i <= radial_disp_range_max_strings.size(); i++ ) {
			real_disps_range_max = std::atof( radial_disp_range_max_strings[i].c_str() );
			radial_disps_range_max.push_back( real_disps_range_max );
		}
		radial_disps_range_max_ = radial_disps_range_max;

		utility::vector1< std::string > angle_range_min_strings = utility::string_split( tag->getOption< std::string >( "angles_range_min" ), ',' );
		utility::vector1<Real> angles_range_min;
		Real real_angles_range_min;
		for ( core::Size i = 1; i <= angle_range_min_strings.size(); i++ ) {
			real_angles_range_min = std::atof( angle_range_min_strings[i].c_str() );
			angles_range_min.push_back( real_angles_range_min );
		}
		angles_range_min_ = angles_range_min;

		utility::vector1< std::string > angle_range_max_strings = utility::string_split( tag->getOption< std::string >( "angles_range_max" ), ',' );
		utility::vector1<Real> angles_range_max;
		Real real_angles_range_max;
		for ( core::Size i = 1; i <= angle_range_max_strings.size(); i++ ) {
			real_angles_range_max = std::atof( angle_range_max_strings[i].c_str() );
			angles_range_max.push_back( real_angles_range_max );
		}
		angles_range_max_ = angles_range_max;

		// If auto_range option is set to true, then the range and signs of the displacements are reversed for negative displacements.
		// This makes it so that a negative value in the range corresponds to displacement toward the origin and a positive value is away from the origin.
		if ( auto_range_ ) {
			utility::vector1<Real> new_radial_disps_range_min, new_radial_disps_range_max;
			for ( core::Size i = 1; i <= radial_disps_.size(); i++ ) {
				if ( radial_disps_[i] < 0 ) {
					new_radial_disps_range_min.push_back(-radial_disps_range_max_[i]);
					new_radial_disps_range_max.push_back(-radial_disps_range_min_[i]);
				} else {
					new_radial_disps_range_min.push_back(radial_disps_range_min_[i]);
					new_radial_disps_range_max.push_back(radial_disps_range_max_[i]);
				}
			}
			radial_disps_range_min_ = new_radial_disps_range_min;
			radial_disps_range_max_ = new_radial_disps_range_max;
		}
	}

	if ( sampling_mode_ == "grid" ) {

		utility::vector1< std::string > angle_step_strings = utility::string_split( tag->getOption< std::string >( "angle_steps" ), ',' );
		utility::vector1<Real> angle_steps;
		Real real_angle_steps;
		for ( core::Size i = 1; i <= angle_step_strings.size(); i++ ) {
			real_angle_steps = std::atof( angle_step_strings[i].c_str() );
			angle_steps.push_back( real_angle_steps );
		}
		angle_steps_ = angle_steps;

		utility::vector1< std::string > radial_disp_step_strings = utility::string_split( tag->getOption< std::string >( "radial_disp_steps" ), ',' );
		utility::vector1<Real> radial_disp_steps;
		Real real_radial_disp_steps;
		for ( core::Size i = 1; i <= radial_disp_step_strings.size(); i++ ) {
			real_radial_disp_steps = std::atof( radial_disp_step_strings[i].c_str() );
			radial_disp_steps.push_back( real_radial_disp_steps );
		}
		radial_disp_steps_ = radial_disp_steps;

		TR << "Setting the exploration grid." << std::endl;

		SymDofMoverSampler::get_instance()->set_angle_ranges(angles_range_min_, angles_range_max_, angle_steps);
		SymDofMoverSampler::get_instance()->set_radial_disp_ranges(radial_disps_range_min_, radial_disps_range_max_, radial_disp_steps);

	}

	if ( sampling_mode_ == "gaussian" ) {

		utility::vector1< std::string > angle_delta_strings = utility::string_split( tag->getOption< std::string >( "angle_deltas" ), ',' );
		utility::vector1<Real> angle_deltas;
		Real real_angle_delta;
		for ( core::Size i = 1; i <= angle_delta_strings.size(); i++ ) {
			real_angle_delta = std::atof( angle_delta_strings[i].c_str() );
			angle_deltas.push_back( real_angle_delta );
		}
		angle_deltas_ = angle_deltas;

		utility::vector1< std::string > radial_disp_delta_strings = utility::string_split( tag->getOption< std::string >( "radial_disp_deltas" ), ',' );
		utility::vector1<Real> radial_disp_deltas;
		Real real_radial_disp_delta;
		for ( core::Size i = 1; i <= radial_disp_delta_strings.size(); i++ ) {
			real_radial_disp_delta = std::atof( radial_disp_delta_strings[i].c_str() );
			radial_disp_deltas.push_back( real_radial_disp_delta );
		}
		radial_disp_deltas_ = radial_disp_deltas;
	}
}

std::string SymDofMover::get_name() const {
	return mover_name();
}

std::string SymDofMover::mover_name() {
	return "SymDofMover";
}

void SymDofMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaRestriction sampling_mode_enumeration;
	sampling_mode_enumeration.name( "sampling_mode_choices" );
	sampling_mode_enumeration.base_type( xs_string );
	sampling_mode_enumeration.add_restriction( xsr_enumeration, "grid" );
	sampling_mode_enumeration.add_restriction( xsr_enumeration, "uniform" );
	sampling_mode_enumeration.add_restriction( xsr_enumeration, "gaussian" );
	sampling_mode_enumeration.add_restriction( xsr_enumeration, "single_dock" );
	xsd.add_top_level_element( sampling_mode_enumeration );

	XMLSchemaRestriction axes_cslist;
	axes_cslist.name( "type_for_axes_cslist" );
	axes_cslist.base_type( xs_string );
	axes_cslist.add_restriction( xsr_pattern, "[xyz](,[xyz])*" );
	xsd.add_top_level_element( axes_cslist );

	XMLSchemaRestriction axes_or_zero_cslist;
	axes_or_zero_cslist.name( "type_for_axes_or_zero_cslist" );
	axes_or_zero_cslist.base_type( xs_string );
	axes_or_zero_cslist.add_restriction( xsr_pattern, "[xyz0](,[xyz0])*" );
	xsd.add_top_level_element( axes_or_zero_cslist );

	/// XRW TO DO: Might be worth adding these patterns to src/utility/tag/XMLSchemaGeneration.cc
	/// XRW TO DO: There are some detailed descriptions of when to use certain attributes based on the result of other attributes, add to documentation.
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "set_sampler" , xsct_rosetta_bool, "For use with the GetRBDOFValues filter. If set to false, then the RBDOF values will not be updated when the reported DOF values are not affected by the SymDofMoverSampler. Defaults to true." , "true" )
		+ XMLSchemaAttribute::attribute_w_default( "auto_range" , xsct_rosetta_bool, "Boolean to set the manner in which the user defined ranges for radial displacements are interpreted. If set to true, then the negative values for min or max displacement are interpreted as moving the structure closer to the origin and positive values away from the origin." , "false" )
		+ XMLSchemaAttribute::required_attribute( "symm_file" , xs_string, "Symmetry definition file." )
		+ XMLSchemaAttribute::required_attribute( "sym_dof_names" , xs_string, "Names of the sym_dofs in the symmetry definition file along which one wishes to move or rotate the input. NOTE: For multicomponent systems, the order of the displacements, angles, ranges, and steps must correspond with the the order of the sym_dof_names. Passed as a string with comma-separated list (e.g. sym_dof_names='JTP1,JDP1')" )
		+ XMLSchemaAttribute::attribute_w_default("sampling_mode", "sampling_mode_choices", "Which mode to use to sample around the initial configuration, if desired. 'grid', 'uniform', or 'gaussian'", "single_dock" )
		+ XMLSchemaAttribute( "flip_input_about_axes" , "type_for_axes_or_zero_cslist", "Rotate subunits 180 degrees about the specified axes prior to applying transtions, rotations, alignment, and symmetry. ie, 'reverse' the component before further manipulation." )
		+ XMLSchemaAttribute( "align_input_axes_to_symdof_axes", "type_for_axes_cslist", "If specified, will align the specified axis of each subunit with the corresponding axis of the symdof jump from the symmetry definition file." )
		+ XMLSchemaAttribute( "radial_disps", xs_string, "Initial displacement(s) by which to translate the input structure(s) along the user specified axis. Passed as a string with a comma-separated list (e.g. radial_disps='-65.4,109.2'). Must also specify translationa_axes to apply these along." )
		+ XMLSchemaAttribute( "radial_offsets", xs_string, "Can be used in any sampling mode. Offset value(s) are added to the corresponding radial_disps before grid, uniform, or gaussian sampling is performed. Works with auto_range. For example, if one wants to space out both symdofs a given structure by 2 angstroms, one can pass radial_offsets='2,2' and auto_range=true. Then, regardless of the sign of the radial disps, the subunits will be displaced 2 angstroms further from the origin (assuming the input subunit(s) start at the origin)." )
		+ XMLSchemaAttribute( "translation_axes", "type_for_axes_cslist", "Need to use this when using radial_disps. Axes (x, y, or z) along which to translate each subunit prior to applying symmetry." )
		+ XMLSchemaAttribute( "angles", xs_string, "Initial angle(s) by which to rotate the input structure(s) about the user specified axis. Passed as a comma-separated list (e.g. angles='-65.4,109.2'). Must also specificy rotation axes to apply these along." )
		+ XMLSchemaAttribute( "rotation_axes", "type_for_axes_cslist", "Axes (x, y, or x) along which to rotate each subunit prior to applying symmetry. Need to use this when using angles."  )
		+ XMLSchemaAttribute( "radial_disps_range_min", xs_string, "For use with grid or uniform sampling. Minimum distance(s) in Angstroms by which to modify the initial radial_disps. Passed as a string with a comma-separated list (e.g. radial_disps_range_min='-1,-1'." )
		+ XMLSchemaAttribute( "radial_disps_range_max", xs_string, "For use with grid or uniform sampling. Maximum distance(s) in Angstroms by which to modify the initial radial_disps. Passed as a string with a comma-separated list (e.g. radial_disps_range_max='1,1'." )
		+ XMLSchemaAttribute( "angles_range_min", xs_string, "For use with grid or uniform sampling. Minimum angle(s) in degrees by which to rotate the structure around the initial angle(s) provided by the user. Passed as a string with a comma-separated list (e.g. angles_range_min='-1,-1'." )
		+ XMLSchemaAttribute( "angles_range_max", xs_string, "For use with grid or uniform sampling. Maximum angle(s) in degrees by which to rotate the structure around the initial angle(s) provided by the user. Passed as a string with a comma-separated list (e.g. angles_range_max='1,1'." )
		+ XMLSchemaAttribute( "angle_steps", xs_string, "For use with grid sampling. Set the bin size(s) by which to sample within the user defined range(s). Passed as a string with a comma-separated list (e.g. angle_steps='0.5,0.5'." )
		+ XMLSchemaAttribute( "radial_disp_steps", xs_string, "For use with Gaussian sampling. The range within to sample inward and outward around the user specified initial displacement(s). Passed as a string with a comma-separated list (e.g. radial_disp_deltas='0.5,0.5'." )
		+ XMLSchemaAttribute( "angle_deltas", xs_string, "For use with Gaussian sampling. The range within to sample around the user specified initial angle(s). Passed as a string with a comma-separated list (e.g. angle_deltas='0.5,0.5'." )
		+ XMLSchemaAttribute( "radial_disp_deltas", xs_string, "For use with Gaussian sampling. The range within to sample inward and outward around the user specified initial displacement(s). Passed as a string with a comma-separated list (e.g. radial_disp_deltas='0.5,0.5'." ) ;


	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Used to setup symmetric systems in which the input structures(s) are aligned along the x, y, or z axis.", attlist );
}

std::string SymDofMoverCreator::keyname() const {
	return SymDofMover::mover_name();
}

protocols::moves::MoverOP
SymDofMoverCreator::create_mover() const {
	return utility::pointer::make_shared< SymDofMover >();
}

void SymDofMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SymDofMover::provide_xml_schema( xsd );
}


} // matdes
} // protocols
