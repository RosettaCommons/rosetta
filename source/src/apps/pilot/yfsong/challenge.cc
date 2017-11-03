// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief protocol to assign sequence
/// @details
/// @author Yifan Song

#include <map>
#include <stdio.h>

#include <devel/init.hh>
#include <core/types.hh>

#include <core/import_pose/import_pose.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Edge.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/SecondaryStructure.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/rational_mc/RationalMonteCarlo.hh>
#include <protocols/rigid/RigidBodyMotionMover.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/nonlocal/StarTreeBuilder.hh>
#include <protocols/nonlocal/util.hh>
#include <protocols/medal/MedalMover.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>
#include <utility/excn/Exceptions.hh>

#include <apps/pilot/yfsong/util.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "pilot.yfsong.challenge" );

namespace challenge {
basic::options::FileVectorOptionKey template_structure("challenge:template_structure");
basic::options::BooleanOptionKey    close_loops("challenge:close_loops");
}

class SecondaryStructureChunk
{
public:
	SecondaryStructureChunk() {
		sequence_position_start_ = 0;
		sequence_ = "";
		ss_ = 'L';
	};
	~SecondaryStructureChunk() {};

	void set_start_seqpos(core::Size seqpos) {
		sequence_position_start_ = seqpos;
	}

	void set_sequence(std::string sequence) {
		sequence_ = sequence;
	}

	void set_sequence(std::string full_sequence, core::Size start_seqpos, core::Size length) {
		sequence_ = full_sequence.substr(start_seqpos-1, length);
	}

	void set_ss(char ss) {
		ss_ = ss;
	}

	void
	append(SecondaryStructureChunk ss_chunk, std::string full_sequence) {
		core::Size seqpos_end = ss_chunk.end_seqpos();
		set_sequence(full_sequence, start_seqpos(), seqpos_end - start_seqpos() + 1 );
	}

	core::Size start_seqpos() {return sequence_position_start_;}
	core::Size end_seqpos() {return sequence_position_start_ + length() - 1;}
	core::Size length() { return sequence_.size(); }
	std::string const sequence() {return sequence_;}
private:
	core::Size sequence_position_start_;
	std::string sequence_;
	char ss_;
};

void
copy_xyz_by_chunk(core::pose::Pose & mod_pose, core::Size const mod_seqpos_start, core::pose::Pose const & ref_pose, core::Size const ref_seqpos_start, core::Size const length)
{
	utility::vector1< core::id::AtomID > ids;
	utility::vector1< numeric::xyzVector<core::Real> > positions;
	for ( core::Size ires = 1; ires <= length; ++ires ) {
		core::Size mod_seqpos = mod_seqpos_start + ires - 1;
		core::Size ref_seqpos = ref_seqpos_start + ires - 1;
		//assert(mod_pose.residue_type(mod_seqpos).name() == ref_pose.residue_type(ref_seqpos).name());

		for ( core::Size iatom = 1; iatom <= ref_pose.residue_type(ref_seqpos).natoms(); ++iatom ) {
			if ( iatom > mod_pose.residue_type(mod_seqpos).natoms() ) break;

			ids.push_back(core::id::AtomID(iatom,mod_seqpos));
			positions.push_back(ref_pose.residue(ref_seqpos).xyz(iatom));
		}
	}
	mod_pose.batch_set_xyz(ids,positions);
}

using namespace core;
using namespace ObjexxFCL;
Real
superimpose_pose(
	pose::Pose & mod_pose,
	std::list <Size> const & residue_list,
	pose::Pose const & ref_pose,
	id::AtomID_Map< id::AtomID > const & atom_map // from mod_pose to ref_pose
)
{
	using namespace numeric::model_quality;
	using namespace id;

	/// how many atoms in mod_pose?
	Size natoms(0);
	for ( Size i=1; i<= mod_pose.size(); ++i ) natoms += mod_pose.residue(i).natoms();

	// pack coords into local arrays for passing into the rms fitting routine
	FArray2D_double xx1(3,natoms);
	FArray2D_double xx2(3,natoms);
	FArray1D_double wt(natoms);


	Size nsup(0);
	{ // pack coordinates into the local arrays
		Size atomno(0);
		Vector const zero_vector(0.0);
		for ( Size i=1; i<= mod_pose.size(); ++i ) {
			for ( Size j=1; j<= mod_pose.residue(i).natoms(); ++j ) {
				++atomno;
				AtomID const & aid( atom_map[ id::AtomID( j,i) ] );
				Vector const & x1( aid.valid() ? ref_pose.xyz( aid ) : zero_vector );
				Vector const & x2( mod_pose.residue(i).xyz(j) );
				wt( atomno ) = ( aid.valid() ? 1.0 : 0.0 );
				if ( aid.valid() ) ++nsup;
				for ( Size k=1; k<= 3; ++k ) {
					xx1(k,atomno) = x1(k);
					xx2(k,atomno) = x2(k);
				}
			}
		}
		runtime_assert( atomno == natoms );
	}

	// calculate starting center of mass (COM):
	FArray1D_double COM(3);
	COMAS(xx1,wt,natoms,COM(1),COM(2),COM(3));

	// superimpose:: shifts xx1, shifts and transforms xx2;
	double rms;
	rmsfitca2(natoms,xx1,xx2,wt,nsup,rms);

	if ( true ) { // debug:
		double tmp1,tmp2,tmp3;
		COMAS(xx1,wt,natoms,tmp1,tmp2,tmp3); // store xcen,ycen,zcen vals for later
		//std::cout << "zero??: " << std::abs(tmp1) + std::abs(tmp2) + std::abs(tmp3)
		//     << std::endl;
		runtime_assert( std::abs(tmp1) + std::abs(tmp2) + std::abs(tmp3) < 1e-3 );
	}

	{ // translate xx2 by COM and fill in the new ref_pose coordinates
		Size atomno(0);
		Vector x2;
		for ( Size i=1; i<= mod_pose.size(); ++i ) {
			for ( Size j=1; j<= mod_pose.residue_type(i).natoms(); ++j ) { // use residue_type to prevent internal coord update
				++atomno;
				if ( find(residue_list.begin(), residue_list.end(), i) == residue_list.end() ) continue;
				for ( Size k=1; k<= 3; ++k ) x2(k) = xx2(k,atomno) + COM(k);
				mod_pose.set_xyz( id::AtomID( j,i), x2 );
			}
		}
		runtime_assert( atomno == natoms );
	}

	return ( static_cast < Real > ( rms ) );
}

// atom_map: from mod_pose to ref_pose
void
get_superposition_transformation(
	pose::Pose const & mod_pose,
	pose::Pose const & ref_pose,
	id::AtomID_Map< id::AtomID > const & atom_map,
	numeric::xyzMatrix< core::Real > &R, numeric::xyzVector< core::Real > &preT, numeric::xyzVector< core::Real > &postT )
{
	using namespace numeric::model_quality;
	using namespace id;

	// count number of atoms for the array
	Size total_mapped_atoms(0);
	for ( Size ires=1; ires<= mod_pose.size(); ++ires ) {
		for ( Size iatom=1; iatom<= mod_pose.residue(ires).natoms(); ++iatom ) {
			AtomID const & aid( atom_map[ id::AtomID( iatom,ires ) ] );
			if ( !aid.valid() ) continue;

			++total_mapped_atoms;
		}
	}

	preT = postT = numeric::xyzVector< core::Real >(0,0,0);
	if ( total_mapped_atoms <= 2 ) {
		R.xx() = R.yy() = R.zz() = 1;
		R.xy() = R.yx() = R.zx() = R.zy() = R.yz() = R.xz() = 0;
		return;
	}

	ObjexxFCL::FArray2D< core::Real > final_coords( 3, total_mapped_atoms );
	ObjexxFCL::FArray2D< core::Real > init_coords( 3, total_mapped_atoms );
	preT = postT = numeric::xyzVector< core::Real >(0,0,0);
	Size atomno(0);
	for ( Size ires=1; ires<= mod_pose.size(); ++ires ) {
		for ( Size iatom=1; iatom<= mod_pose.residue(ires).natoms(); ++iatom ) {
			AtomID const & aid( atom_map[ id::AtomID( iatom,ires ) ] );
			if ( !aid.valid() ) continue;
			++atomno;

			numeric::xyzVector< core::Real > x_i = mod_pose.residue(ires).atom(iatom).xyz();
			preT += x_i;
			numeric::xyzVector< core::Real > y_i = ref_pose.xyz( aid );
			postT += y_i;

			for ( int j=0; j<3; ++j ) {
				init_coords(j+1,atomno) = x_i[j];
				final_coords(j+1,atomno) = y_i[j];
			}
		}
	}

	preT /= (float) total_mapped_atoms;
	postT /= (float) total_mapped_atoms;
	for ( int i=1; i<=(int)total_mapped_atoms; ++i ) {
		for ( int j=0; j<3; ++j ) {
			init_coords(j+1,i) -= preT[j];
			final_coords(j+1,i) -= postT[j];
		}
	}

	// get optimal superposition
	// rotate >init< to >final<
	ObjexxFCL::FArray1D< numeric::Real > ww( total_mapped_atoms, 1.0 );
	ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
	numeric::Real ctx;

	numeric::model_quality::findUU( init_coords, final_coords, ww, total_mapped_atoms, uu, ctx );
	R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
	R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
	R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );
}

void
apply_transform(
	pose::Pose & mod_pose,
	std::list <Size> const & residue_list,
	numeric::xyzMatrix< core::Real > const & R, numeric::xyzVector< core::Real > const & preT, numeric::xyzVector< core::Real > const & postT
)
{ // translate xx2 by COM and fill in the new ref_pose coordinates
	utility::vector1< core::id::AtomID > ids;
	utility::vector1< numeric::xyzVector<core::Real> > positions;

	Vector x2;
	FArray2D_double xx2;
	FArray1D_double COM(3);
	for ( std::list<Size>::const_iterator it = residue_list.begin();
			it != residue_list.end();
			++it ) {
		Size ires = *it;
		for ( Size iatom=1; iatom<= mod_pose.residue_type(ires).natoms(); ++iatom ) { // use residue_type to prevent internal coord update
			ids.push_back(core::id::AtomID(iatom,ires));
			positions.push_back(postT + (R*( mod_pose.xyz(core::id::AtomID(iatom,ires)) - preT )));
		}
	}
	mod_pose.batch_set_xyz(ids,positions);
}

void
superimpose_pose_transform(
	pose::Pose & mod_pose,
	std::list <Size> const & residue_list,
	pose::Pose const & ref_pose,
	id::AtomID_Map< id::AtomID > const & atom_map // from mod_pose to ref_pose
)
{
	numeric::xyzMatrix< core::Real > R;
	numeric::xyzVector< core::Real > preT;
	numeric::xyzVector< core::Real > postT;
	get_superposition_transformation( mod_pose, ref_pose, atom_map, R, preT, postT );
	apply_transform( mod_pose, residue_list, R, preT, postT );

}


void align_backbone_by_chunk(core::pose::Pose & pose, core::Size const residue_seq_start, core::Size const residue_seq_end,
	core::pose::Pose const & ref_pose,
	std::map <core::Size, core::Size> const & seqpos_alignment,
	int const registry_shift=0,
	core::Size MAX_TRIAL = 100)
{
	using namespace ObjexxFCL::format;
	TR << "aligning the chunk: " << I(4,residue_seq_start) << I(4,residue_seq_end) << std::endl;

	std::list < Size > residue_list;
	for ( Size ires=residue_seq_start; ires<=residue_seq_end; ++ires ) {
		residue_list.push_back(ires);
	}

	//core::pose::Pose local_pose(pose, residue_seq_start, residue_seq_end);
	/*
	for ( core::Size i_res = 1; i_res <= local_pose.size(); i_res++ ) {
	if ( ! local_pose.residue_type(i_res).is_protein() ) continue;
	local_pose.set_phi(i_res,    -70);
	local_pose.set_psi(i_res,    130);
	local_pose.set_omega(i_res,  180);
	}
	*/
	core::Size counter = 0;
	while ( counter < MAX_TRIAL ) {
		++counter;
		core::Size random_res = numeric::random::rg().random_range(residue_seq_start, residue_seq_end);
		core::Size seqpos_model = random_res+registry_shift;

		if ( seqpos_alignment.find(seqpos_model) == seqpos_alignment.end() ) continue;
		int seqpos_shift = seqpos_alignment.find(seqpos_model)->second - random_res;

		core::id::AtomID_Map< core::id::AtomID > atom_map;
		core::pose::initialize_atomid_map( atom_map, pose, core::id::AtomID::BOGUS_ATOM_ID() );
		core::Size atom_map_count = 0;
		for ( Size ires=residue_seq_start; ires<=residue_seq_end; ++ires ) {
			if ( ! pose.residue_type(ires).is_protein() ) continue;
			pose.set_phi(ires, -110);
			pose.set_psi(ires,  130);
			pose.set_omega(ires, 180);

			int seqpos_ref = ires + seqpos_shift;
			if ( seqpos_ref <= 0 || seqpos_ref > ref_pose.size() ) continue;
			if ( !ref_pose.residue_type(seqpos_ref).is_protein() ) continue;
			if ( !ref_pose.secstruct(seqpos_ref) == 'L' ) continue;
			pose.set_phi(ires, ref_pose.phi(seqpos_ref));
			pose.set_psi(ires, ref_pose.psi(seqpos_ref));
			pose.set_omega(ires, ref_pose.omega(seqpos_ref));

			core::id::AtomID const id1( pose.residue_type(ires).atom_index("CA"), ires );
			core::id::AtomID const id2( ref_pose.residue_type(seqpos_ref).atom_index("CA"), seqpos_ref );
			atom_map[ id1 ] = id2;
			++atom_map_count;


			using namespace ObjexxFCL::format;
			TR << "atom mapping " << I(4, atom_map_count) << id1 << id2 << std::endl;
		}

		if ( atom_map_count < 3 ) continue;
		//pose.dump_pdb("before_superimpose.pdb");
		superimpose_pose_transform(pose, residue_list, ref_pose, atom_map);
		//pose.dump_pdb("after_superimpose.pdb");
		//using namespace ObjexxFCL::format;
		//TR << "after  superimposing model " << I(4,atom_map.size()) << I(4,atom_map_count) << F(8,3,rms) << std::endl;
		//  if (rms > 5.0) {
		//   utility_exit_with_message("Large RMSD");
		//  }
		break;
	}
}

void
extract_ss_chunks_from_seq(std::string const secstructs,
	core::pose::Pose const & pose,
	std::string const extracted_ss_types,
	utility::vector1< SecondaryStructureChunk > & ss_chunks,
	core::Size max_gap_in_continuous_chunk = 1, // if two chunks are seperated by a gap of this size (or less), consider it one big chunk
	core::Size minimum_length_of_chunk = 3,
	core::Real CA_CA_distance_cutoff = 4.0)
{
	assert(pose.size() == secstructs.size());
	ss_chunks.clear();
	core::Size last_extracted_ss_end = 0; // so that different types of secondary structure chunks are not merged together

	for ( core::Size i_ss = 0; i_ss < extracted_ss_types.size(); ++i_ss ) {
		char ss = extracted_ss_types[i_ss];

		bool ss_chunk_started = false;
		//core::Size ss_chunk_start, ss_chunk_end;

		SecondaryStructureChunk ss_chunk;

		for ( core::Size ires = 1; ires <= secstructs.size(); ++ires ) {
			if ( !ss_chunk_started ) {
				if ( secstructs[ires-1] == ss ) {

					ss_chunk_started = true;
					ss_chunk.set_start_seqpos(ires);
				}
			} else {
				if ( secstructs[ires-1] != ss ) {
					ss_chunk_started = false;
					ss_chunk.set_sequence(pose.sequence().substr(ss_chunk.start_seqpos() -1, ires-ss_chunk.start_seqpos()));
					ss_chunks.push_back(ss_chunk);
				} else if ( ! pose.residue_type(ires).is_protein() ) {
					ss_chunk_started = false;
					ss_chunk.set_sequence(pose.sequence().substr(ss_chunk.start_seqpos() -1, ires-ss_chunk.start_seqpos()));
					ss_chunks.push_back(ss_chunk);
				} else if ( pose.residue(ires).xyz("CA").distance(pose.residue(ires-1).xyz("CA")) > CA_CA_distance_cutoff ) {
					ss_chunk_started = false;
					ss_chunk.set_sequence(pose.sequence().substr(ss_chunk.start_seqpos() -1, ires-ss_chunk.start_seqpos()));
					ss_chunks.push_back(ss_chunk);
				}
			}
		}

		// if the input sequence ends with the ss to be extracted
		if ( ss_chunk_started ) {
			ss_chunk.set_sequence(pose.sequence().substr(ss_chunk.start_seqpos() -1, secstructs.size()-ss_chunk.start_seqpos()+1));
			ss_chunks.push_back(ss_chunk);
		}

		// join ss_chunks seperated by a small gap
		for ( int i_chunk=ss_chunks.size()-1; i_chunk > (int) last_extracted_ss_end; --i_chunk ) {
			if ( ss_chunks[i_chunk+1].start_seqpos() - ss_chunks[i_chunk].end_seqpos() <= max_gap_in_continuous_chunk+1 ) {
				ss_chunks[i_chunk].append(ss_chunks[i_chunk+1], pose.sequence());
				ss_chunks.erase(ss_chunks.begin() + i_chunk);
			}
		}

		//remove short ss_chunks
		for ( int i_chunk=ss_chunks.size()-1; i_chunk > (int) last_extracted_ss_end; --i_chunk ) {
			if ( ss_chunks[i_chunk].length() < minimum_length_of_chunk ) {
				ss_chunks.erase(ss_chunks.begin() + i_chunk - 1);
			}
		}

		last_extracted_ss_end = ss_chunks.size();
	}
}

void read_template_structures(utility::vector1 <core::pose::PoseOP> & template_structures)
{
	template_structures.clear();

	// read reference structures
	if ( !basic::options::option[ challenge::template_structure ].user() ) {
		utility_exit_with_message("Error! Need the -challenge::template_structure flag.");
	}

	utility::vector1 < utility::file::FileName > ref_filenames = basic::options::option[ challenge::template_structure ];
	template_structures.resize(ref_filenames.size());
	for ( core::Size i_ref=1; i_ref<= ref_filenames.size(); ++i_ref ) {
		template_structures[i_ref] = new core::pose::Pose();
		core::import_pose::pose_from_file( *(template_structures[i_ref]), ref_filenames[i_ref] , core::import_pose::PDB_file);

		core::scoring::dssp::Dssp dssp_obj( *template_structures[i_ref] );
		dssp_obj.insert_ss_into_pose( *template_structures[i_ref] );
	}

}

void
split_pdb_into_ss_chunks(core::pose::PoseOP const pose, utility::vector1< SecondaryStructureChunk > & ss_chunks) {
	core::scoring::dssp::Dssp dssp_obj( *pose );
	dssp_obj.insert_ss_into_pose( *pose );

	std::string secstruct = pose->secstruct();
	TR << "Secondary structure: " << secstruct << std::endl;
	extract_ss_chunks_from_seq(secstruct, *pose, "HE", ss_chunks);
}

using protocols::nonlocal::StarTreeBuilder;
using namespace protocols::nonlocal;
using namespace protocols::loops;
using namespace core::kinematics;

class ChallengeMover : public protocols::moves::Mover {

public:
	void
	initialize_template_structures() {
		template_ss_chunks_.clear();
		seqpos_alignments_.clear();
		read_template_structures(template_structures_);
		template_ss_chunks_.resize(template_structures_.size());

		for ( core::Size i_template=1; i_template<=template_structures_.size(); ++i_template ) {
			// utility::vector1< SecondaryStructureChunk > ss_chunks;
			// split_pdb_into_ss_chunks(template_structures_[i_template], ss_chunks);
			//using namespace ObjexxFCL::format;
			//TR << I(4,ss_chunks.size()) << std::endl;
			//for (core::Size i_chunk=1; i_chunk<=ss_chunks.size(); ++i_chunk) {
			// using namespace ObjexxFCL::format;
			// TR << I(4, i_chunk) << I(4,ss_chunks[i_chunk].start_seqpos()) << I(4,ss_chunks[i_chunk].end_seqpos()) << I(4,ss_chunks[i_chunk].length()) << " " << ss_chunks[i_chunk].sequence() << std::endl;
			//}
			//TR.flush();

			// template_ss_chunks_.push_back(ss_chunks);

			std::map <core::Size, core::Size> seqpos_alignment;
			get_alignment_from_template(template_structures_[i_template], seqpos_alignment);
			seqpos_alignments_.push_back(seqpos_alignment);

			// Build the star fold tree, identify jumps
			core::scoring::dssp::Dssp dssp_obj( *template_structures_[i_template] );
			dssp_obj.insert_ss_into_pose( *template_structures_[i_template] );

			template_ss_chunks_[i_template] = extract_secondary_structure_chunks(*template_structures_[i_template]);
			StarTreeBuilder builder;
			if ( template_ss_chunks_[i_template].num_loop() > 0 ) {
				builder.set_up(template_ss_chunks_[i_template], &(*template_structures_[i_template]));
			}
		}
	}

	utility::vector1 <core::Real>
	count_alignment(core::pose::Pose & mod_pose, utility::vector1 < std::map <core::Size, core::Size> > & seqpos_alignments) {
		utility::vector1 <core::Real> alignment_counts;
		for ( core::Size ires=1; ires<=mod_pose.size(); ++ires ) {
			core::Size num_aligned = 0;
			for ( core::Size i_template=1; i_template<=seqpos_alignments.size(); ++i_template ) {
				if ( seqpos_alignments[i_template].find(ires) == seqpos_alignments[i_template].end() ) continue;
				++num_aligned;
			}
			alignment_counts.push_back((core::Real)num_aligned/(core::Real)seqpos_alignments.size());
			using namespace ObjexxFCL::format;
			TR << I(4,ires) << F(8,3, (core::Real)num_aligned/(core::Real)seqpos_alignments.size()) << std::endl;
			TR.flush();
		}
		return alignment_counts;
	}

	protocols::loops::Loops
	extract_chunks_from_alignment_counts(utility::vector1 <core::Real> alignment_counts)
	{
		core::Real average=0.;
		for ( core::Size ires=1; ires<=alignment_counts.size(); ++ires ) {
			average+=alignment_counts[ires];
		}
		if ( alignment_counts.size() > 0 ) {
			average /= (core::Real) alignment_counts.size();
		}

		Loops chunks;

		bool chunk_started = false;
		Size chunk_start_seqpos(0);
		Size chunk_end_seqpos(0);

		for ( core::Size ires = 1; ires <= alignment_counts.size(); ++ires ) {
			if ( !chunk_started ) {
				if ( alignment_counts[ires] >= average ) {

					chunk_started = true;
					chunk_start_seqpos = ires;
				}
			} else {
				if ( alignment_counts[ires] < average ) {
					chunk_started = false;
					chunk_end_seqpos = ires - 1;
					chunks.add_loop( chunk_start_seqpos, chunk_end_seqpos);
				}
			}
		}

		// if the chunk is not ended, add the last chunk
		if ( chunk_started ) {
			chunk_end_seqpos = alignment_counts.size();
			chunks.add_loop( chunk_start_seqpos, chunk_end_seqpos);
		}
		return chunks;
	}

	void
	extract_ss_chunks_from_alignment_counts(utility::vector1 <core::Real> alignment_counts, std::string const full_sequence, utility::vector1< SecondaryStructureChunk > & ss_chunks)
	{
		ss_chunks.clear();

		core::Real average=0.;
		for ( core::Size ires=1; ires<=alignment_counts.size(); ++ires ) {
			average+=alignment_counts[ires];
		}
		if ( alignment_counts.size() > 0 ) {
			average /= (core::Real) alignment_counts.size();
		}

		bool ss_chunk_started = false;
		SecondaryStructureChunk ss_chunk;
		for ( core::Size ires=1; ires<=alignment_counts.size(); ++ires ) {
			if ( !ss_chunk_started ) {
				if ( alignment_counts[ires] >= average ) {

					ss_chunk_started = true;
					ss_chunk.set_start_seqpos(ires);
				}
			} else {
				if ( alignment_counts[ires] < average ) {
					ss_chunk_started = false;
					ss_chunk.set_sequence(full_sequence.substr(ss_chunk.start_seqpos() -1, ires-ss_chunk.start_seqpos()));
					ss_chunks.push_back(ss_chunk);
				}
			}
		}

		// if the chunk is not ended, add the last chunk
		if ( ss_chunk_started ) {
			ss_chunk.set_sequence(full_sequence.substr(ss_chunk.start_seqpos() -1, alignment_counts.size()-ss_chunk.start_seqpos()+1));
			ss_chunks.push_back(ss_chunk);
		}

		// join ss_chunks seperated by less than 3 residues (1 or 2)
		for ( int i_chunk=ss_chunks.size()-1; i_chunk >= 1; --i_chunk ) {
			if ( ss_chunks[i_chunk+1].start_seqpos() - ss_chunks[i_chunk].end_seqpos() <= 0 ) {
				ss_chunks[i_chunk].append(ss_chunks[i_chunk+1], full_sequence);
				ss_chunks.erase(ss_chunks.begin() + i_chunk);
			}
		}
	}

	core::Size
	pick_random_template() {
		core::Size i_template = numeric::random::rg().random_range(1, template_structures_.size());
		if ( template_ss_chunks_[i_template].size() != 0 ) return i_template;
		else return 0;
	}

	void
	get_alignment_from_template(core::pose::PoseOP const template_pose, std::map <core::Size, core::Size> & seqpos_alignment) {
		// specific to this case, alignment comes from residue number
		for ( core::Size ires=1; ires<=template_pose->size(); ++ires ) {
			seqpos_alignment[template_pose->pdb_info()->number(ires)] = ires;
		}
	}

	ChallengeMover()
	{
		initialize_template_structures();
	}

	void initialize_pose(core::pose::Pose & pose)
	{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		utility::vector1 <core::Real> alignment_stats = count_alignment(pose, seqpos_alignments_);
		//extract_ss_chunks_from_alignment_counts(alignment_stats, pose.sequence(), ss_chunks_target_);

		if ( option[ OptionKeys::in::file::psipred_ss2 ].user() ) {
			bool check_psipred = set_secstruct_from_psipred_ss2(pose);
			assert (check_psipred);

			// Build the star fold tree, identify jumps
			ss_chunks_target_ = extract_secondary_structure_chunks( pose );
			loops_target_ = ss_chunks_target_.invert(pose.size());

			TR << ss_chunks_target_;

			StarTreeBuilder builder;
			if ( ss_chunks_target_.num_loop() > 0 ) {
				builder.set_up(ss_chunks_target_, &pose);
			}
			protocols::nonlocal::add_cutpoint_variants(&pose);
		} else {
			utility_exit_with_message("Error in reading psipred_ss2 file, is the -in:file:psipred_ss2 flag set correctly?");
		}

		using namespace ObjexxFCL::format;
		TR << I(4,ss_chunks_target_.size()) << std::endl;
		for ( core::Size i_chunk=1; i_chunk<=ss_chunks_target_.size(); ++i_chunk ) {
			using namespace ObjexxFCL::format;
			TR << I(4, i_chunk) << I(4,ss_chunks_target_[i_chunk].start()) << I(4,ss_chunks_target_[i_chunk].stop()) << I(4,ss_chunks_target_[i_chunk].length()) << std::endl;
		}
		TR.flush();


		//protocols::nonlocal::add_cutpoint_variants(&pose);

		core::Size i_template(0);
		while ( !i_template ) {
			i_template = pick_random_template();
		}

		using namespace ObjexxFCL::format;
		TR << pose.fold_tree() << std::endl;
		TR << "template " << I(4, i_template) << std::endl;

		using core::kinematics::Jump;
		for ( core::Size jump_number=1; jump_number<=pose.num_jump(); ++jump_number ) {

			std::list < Size > downstream_residues = downstream_residues_from_jump(pose, jump_number);
			Size seqpos_start=downstream_residues.front();
			Size seqpos_end=downstream_residues.back();
			align_backbone_by_chunk(pose, seqpos_start, seqpos_end, *(template_structures_[i_template]), seqpos_alignments_[i_template]);
		}
	}

	void realign(core::pose::Pose & pose) {
		core::Size i_template(0);
		while ( !i_template ) {
			i_template = pick_random_template();
		}
		using namespace ObjexxFCL::format;
		TR << "template " << I(4, i_template) << std::endl;

		core::Size jump_number = numeric::random::rg().random_range(1, pose.num_jump());

		std::list < Size > downstream_residues = downstream_residues_from_jump(pose, jump_number);
		Size seqpos_start=downstream_residues.front();
		Size seqpos_end=downstream_residues.back();

		int max_registry_shift = 1;
		int registry_shift = numeric::random::rg().random_range(-max_registry_shift, max_registry_shift);
		align_backbone_by_chunk(pose, seqpos_start, seqpos_end, *(template_structures_[i_template]), seqpos_alignments_[i_template], registry_shift);
	}

	std::list < Size >
	downstream_residues_from_jump(core::pose::Pose & pose, Size jump_number) {
		std::list < Size > residue_list;
		utility::vector1< Edge > edges = pose.fold_tree().get_outgoing_edges(pose.fold_tree().jump_edge(jump_number).stop());

		for ( Size i_edge = 1; i_edge <= edges.size(); ++i_edge ) {
			if ( !edges[i_edge].is_polymer() ) continue;
			Size start = edges[i_edge].start() <= edges[i_edge].stop() ? edges[i_edge].start() : edges[i_edge].stop();
			Size stop  = edges[i_edge].start() <= edges[i_edge].stop() ? edges[i_edge].stop()  : edges[i_edge].start();
			for ( Size ires = start; ires <= stop; ++ires ) {
				residue_list.push_back(ires);
			}

		}
		residue_list.sort();
		residue_list.unique();
		return residue_list;
	}


	void
	apply ( core::pose::Pose & pose )
	{
		initialize_pose(pose);

		protocols::simple_moves::ConstraintSetMoverOP const_set = new protocols::simple_moves::ConstraintSetMover();
		const_set->apply(pose);

		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		(*scorefxn)(pose);

		sample(pose, scorefxn);

		// loop closure
		if ( basic::options::option[ challenge::close_loops ] ) {
			scorefxn->set_weight(core::scoring::vdw, 1);
			scorefxn->set_weight(core::scoring::chainbreak, 1);

			read_loop_fragments(frag_libs_);
			protocols::comparative_modelling::LoopRelaxMoverOP lr_mover( new  protocols::comparative_modeling::LoopRelaxMover );
			lr_mover->frag_libs(frag_libs_);
			lr_mover->loops(loops_target_);
			lr_mover->cen_scorefxn(scorefxn);
			//lr_mover->n_rebuild_tries( 10 );
			lr_mover->apply( pose );
		}
	}

	core::Size num_trials() {
		return 1000;
	}

	typedef boost::unordered_map<int, core::kinematics::Jump> Jumps;
	void jumps_from_pose(const core::pose::Pose& pose, Jumps* jumps) const {
		using core::kinematics::Jump;
		assert(jumps);

		for ( core::Size i = 1; i <= pose.num_jump(); ++i ) {
			const Jump& jump = pose.jump(i);
			(*jumps)[i] = jump;
			TR.Debug << "Added jump_num " << i << ": " << jump << std::endl;
		}
	}

	// use cmiles's RB mover
	void do_rigid_body_moves(const core::scoring::ScoreFunctionOP& score,
		core::pose::Pose* pose) const {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using protocols::moves::MoverOP;
		using protocols::simple_moves::rational_mc::RationalMonteCarlo;
		using protocols::rigid::RigidBodyMotionMover;
		assert(pose);

		Jumps jumps;
		jumps_from_pose(*pose, &jumps);

		// Rigid body motion
		MoverOP rigid_mover = new RationalMonteCarlo(
			new rigid::RigidMotionMover(jumps),
			score,
			option[OptionKeys::rigid::rigid_body_cycles](),
			option[OptionKeys::rigid::temperature](),
			true);

		TR <<"Beginning rigid body perturbation phase..." << std::endl;
		rigid_mover->apply(*pose);
	}

	void sample( core::pose::Pose & pose, core::scoring::ScoreFunctionOP & scorefxn) {
		protocols::moves::MonteCarloOP mc_ = new protocols::moves::MonteCarlo(*scorefxn, 2.0);
		mc_->reset(pose);
		mc_->reset_counters();

		for ( Size i = 1; i <= num_trials(); ++i ) {
			// retain a copy of the pose in the event that the move is rejected
			core::pose::Pose copy(pose);

			realign(pose);

			//do_rigid_body_moves(scorefxn, &pose);

			if ( !mc_->boltzmann(pose) ) {
				pose = copy;
			}

			mc_->score_function().show(TR, pose);
			TR.flush();
		}

		// optionally recover the low-scoring pose
		mc_->recover_low(pose);

		// show simulation statistics
		mc_->show_counters();
		mc_->score_function().show(TR, pose);
		TR.flush();
	}

	std::string
	get_name() const {
		return "ChallengeMover";
	}

private:
	utility::vector1 < core::pose::PoseOP > template_structures_;
	utility::vector1 < Loops > template_ss_chunks_;

	// utility::vector1 < utility::vector1 < SecondaryStructureChunk > > template_ss_chunks_;
	utility::vector1 < std::map <core::Size, core::Size> > seqpos_alignments_;

	Loops ss_chunks_target_;
	Loops loops_target_;
	utility::vector1< core::fragment::FragSetOP > frag_libs_;
};

/*
void *
challenge_main( void* ) {
protocols::jd2::JobDistributor::get_instance()->go(new ChallengeMover());
return 0;
}
*/

void *
my_main( void* ) {
	using namespace protocols::moves;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( !basic::options::option[ challenge::template_structure ].user() ) {
		utility_exit_with_message("Error! Need the -challenge::template_structure flag.");
	}
	utility::vector1 < utility::file::FileName > template_filenames = basic::options::option[ challenge::template_structure ];

	SequenceMoverOP whole_sequence( new SequenceMover() );
	if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_file ].user() ) {
		whole_sequence->add_mover(new protocols::simple_moves::ConstraintSetMover());
	}
	apps::pilot::SampleSecondaryStructureAlignmentMover * sample_ss ( new apps::pilot::SampleSecondaryStructureAlignmentMover(numeric::random::rg(), template_filenames) );
	whole_sequence->add_mover(sample_ss);

	if ( basic::options::option[ challenge::close_loops ] ) {
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		scorefxn->set_weight(core::scoring::rama, 0.4);
		scorefxn->set_weight(core::scoring::vdw, 0.8);

		//Loops ss_chunks_target = extract_secondary_structure_chunks( pose );
		//Loops loops_target = ss_chunks_target.invert(pose.size());

		utility::vector1< core::fragment::FragSetOP > frag_libs;
		read_loop_fragments(frag_libs);
		protocols::comparative_modeling::LoopRelaxMoverOP lr_mover( new protocols::comparative_modeling::LoopRelaxMover );
		lr_mover->frag_libs(frag_libs);
		//lr_mover->loops(loops_target);
		lr_mover->cen_scorefxn(scorefxn);
		whole_sequence->add_mover(lr_mover);
	}

	try{
		protocols::jd2::JobDistributor::get_instance()->go( whole_sequence );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return 0;
}

int
main( int argc, char * argv [] )
{
	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		option.add( challenge::template_structure, "template pdb structures" );
		option.add( challenge::close_loops, "run loop relax to close loops" );
		option.add( challenge::ss, "secondary structure" );
		option.add( challenge::aligned, "alignment from template resSeq" );
		option.add( challenge::max_registry_shift, "maximum registry shift allowed" ).def(0);
		option.add( challenge::virtual_loops, "set loops to virtual residues" ).def(false);
		option.add( challenge::revert_real_loops, "reset loops to real residues" ).def(false);
		option.add( challenge::chunk_mapping, "a vector of chain numbers mapping onto secondary structure chunks" );

		devel::init( argc, argv );

		//protocols::viewer::viewer_main( challenge_main );
		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
