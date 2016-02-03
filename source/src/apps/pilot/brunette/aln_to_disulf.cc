// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file aln_to_disulf.cc
/// @brief
/// @author TJ Brunette


#include <core/types.hh>
#include <devel/init.hh>

#include <basic/Tracer.hh>

#include <core/chemical/util.hh>

#include <core/id/SequenceMapping.hh>
#include <basic/options/option.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <utility/io/ozstream.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>


#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/annotated_sequence.hh>


#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <map>
#include <iostream>
#include <string>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>

// Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>

using core::Size;
using core::Real;
namespace aln_to_disulf {
  basic::options::FileVectorOptionKey coordCstFiles("minimalCstHomology:coordCstFiles");
	basic::options::BooleanOptionKey only_res_out("minimalCstHomology:only_res_out");
}
utility::vector1< std::pair< Size, Size> > get_disulfides_from_aln( core::pose::Pose templatePose, core::sequence::SequenceAlignment aln){
	using std::string;
	using core::pose::Pose;
	using utility::vector1;
	using core::id::SequenceMapping;
	vector1<std::pair <Size,Size> > disulfides;
	core::id::SequenceMapping aln_map( aln.sequence_mapping(1,2) );;
	vector1<Size> potential_disulf;
	Size seqPos = aln.sequence(1)->start()-1;
	for ( core::Size ii = 1; ii <= aln.sequence(1)->length(); ++ii ) {
		if(aln.sequence(1)->at(ii) != '-'){
			seqPos++;
		}
		if((aln.sequence(1)->at(ii) == 'C') && (aln.sequence(2)->at(ii) == 'C')){
			potential_disulf.push_back(seqPos);
		}
	}
	typedef vector1< Size >::const_iterator iter;
	bool fullatom(templatePose.is_fullatom());
	Real const typical_disulfide_distance = fullatom? 2.02 : 3.72;
	Real const tolerance = fullatom? 0.25 : 1.0;
	std::string const distance_atom = fullatom? "SG" : "CB";
	for(iter it1 = potential_disulf.begin(), end = potential_disulf.end(); it1 != end; ++it1){
		for(iter it2 = it1+1; it2 != end; ++it2){
			Real dist = templatePose.residue(aln_map[*it1]).xyz(distance_atom).distance(templatePose.residue(aln_map[*it2]).xyz(distance_atom));
			if(dist < typical_disulfide_distance+tolerance){
				disulfides.push_back(std::pair<Size,Size>(*it1,*it2));
			}
		}
	}
	return(disulfides);
}

std::map< std::string, core::pose::Pose >
poses_from_cmd_line(
	utility::vector1< std::string > const & fn_list
) {
	using std::map;
	using std::string;
	using core::pose::Pose;
	using utility::file::file_exists;
	using utility::vector1;
	using core::import_pose::pose_from_file;
	using namespace core::chemical;

	ResidueTypeSetCOP rsd_set( rsd_set_from_cmd_line() );
	map< string, Pose > poses;

	typedef vector1< string >::const_iterator iter;
	for ( iter it = fn_list.begin(), end = fn_list.end(); it != end; ++it ) {
		if ( file_exists(*it) ) {
			Pose pose;
			core::import_pose::pose_from_file( pose, *rsd_set, *it , core::import_pose::PDB_file);
			string name = utility::file_basename( *it );
			name = name.substr( 0, 5 );
			poses[name] = pose;
		}
	}

	return poses;
}

int
main( int argc, char* argv [] ) {
	try {
	// options, random initialization
	devel::init( argc, argv );

	using std::map;
	using std::multimap;
	using std::string;
	using core::Real;
	using core::Size;
	using core::pose::Pose;
	using core::pose::PoseOP;
	using utility::vector1;
	using core::import_pose::pose_from_file;
	using core::pose::make_pose_from_sequence;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::sequence;
	basic::Tracer tr( "aln_to_disulf" );


	vector1< string > align_fns = option[ in::file::alignment ]();

	map< string, Pose > poses = poses_from_cmd_line(
			option[ in::file::template_pdb ]()
	);

	//-----Loop through alignments and get potential disulfides
	typedef vector1< string >::const_iterator aln_iter;
	vector1<std::pair <Size,Size> > disulfides;
	for ( aln_iter aln_fn = align_fns.begin(), aln_end = align_fns.end();
				aln_fn != aln_end; ++aln_fn
	) {
		vector1< SequenceAlignment > alns = core::sequence::read_aln(
			option[ cm::aln_format ](), *aln_fn
		);

		for ( vector1< SequenceAlignment >::iterator it = alns.begin(),
				end = alns.end();
				it != end; ++it
		) {
			string const template_id( it->sequence(2)->id().substr(0,5) );
			tr << *it << std::endl;
			tr << "id " << it->sequence(2)->id() << " => " << template_id
				<< std::endl;
			map< string, Pose >::iterator pose_it = poses.find( template_id );
			if ( pose_it == poses.end() ) {
				string msg( "Error: can't find pose (id = "
					+ template_id + ")"
				);
				utility_exit_with_message(msg);
				//tr.Error << msg << std::endl;
			} else {
				Pose template_pose;
				template_pose = pose_it->second;
				vector1<std::pair <Size,Size> > tmp_disulfides = get_disulfides_from_aln(template_pose,*it);
				disulfides.insert(disulfides.begin(),tmp_disulfides.begin(),tmp_disulfides.end());
			}
		}
	}
	typedef vector1< std::pair < Size, Size > >::const_iterator disulf_iter;
	typedef map< std::pair < Size, Size >, Size>::iterator disulf_map_iter;
	typedef multimap< Size , std::pair < Size, Size > >::reverse_iterator ct_disulf_map_iter;
	map< std::pair<Size ,Size> , Size> disulfides_ct_map;
	//Store in map to get a count for each observed contact
	for ( disulf_iter disulf_it = disulfides.begin(), disulf_end = disulfides.end(); disulf_it != disulf_end; ++disulf_it){
		disulf_map_iter disulfide_map_loc = disulfides_ct_map.find(*disulf_it);
		if( disulfide_map_loc == disulfides_ct_map.end())
			disulfides_ct_map[*disulf_it] = 1;
		else{
			disulfide_map_loc->second++;
		}
	}
	//convert to multimap so the cts are ordered. Then go through the list and get rid of lower ct contacts.
	multimap< Size , std::pair < Size, Size > > ct_disulf_map;
	for ( disulf_map_iter disulf_map_it = disulfides_ct_map.begin(); disulf_map_it != disulfides_ct_map.end(); ++disulf_map_it)
		ct_disulf_map.insert(std::pair<Size , std::pair<Size , Size > >(disulf_map_it->second, std::pair<Size,Size>(disulf_map_it->first.first,disulf_map_it->first.second)));
	vector1<std::pair < Size,Size> > final_disulf_pairs;

	for (ct_disulf_map_iter ct_disulf_map_it = 	ct_disulf_map.rbegin(); ct_disulf_map_it  !=  ct_disulf_map.rend(); ++ct_disulf_map_it){
		bool found = false;
		if(ct_disulf_map_it->first <= 1)
			found = true;
		for ( disulf_iter disulf_it = final_disulf_pairs.begin(), disulf_end = final_disulf_pairs.end(); disulf_it != disulf_end; ++disulf_it){
			if ( ct_disulf_map_it->second.first == disulf_it->first || ct_disulf_map_it->second.second == disulf_it->first || ct_disulf_map_it->second.first == disulf_it->second || ct_disulf_map_it->second.second == disulf_it->second)
				found = true;
		}
		if(found == false)
			final_disulf_pairs.push_back(ct_disulf_map_it->second);
	}
	if(final_disulf_pairs.size() > 0){
		string out_nametag = "disulf.txt";
		if ( basic::options::option[ out::file::o ].user() ) {
			out_nametag = option[out::file::o ]();
		}

		utility::io::ozstream output(out_nametag);

		for ( disulf_iter disulf_it = final_disulf_pairs.begin(), disulf_end = final_disulf_pairs.end(); disulf_it != disulf_end; ++disulf_it){

			output << disulf_it->first << " " <<  disulf_it->second << std::endl;
		}
	}
	tr << "disulfide detection completed successfully" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
