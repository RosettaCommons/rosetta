// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Oliver Lange

// Unit Headers
#include <protocols/abinitio/Templates.hh>

// Package Headers
#include <protocols/abinitio/Template.hh>
#include <protocols/abinitio/TemplateJumpSetup.hh>
#include <protocols/abinitio/PairingStatistics.hh>

// Project Headers
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID.hh>


// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragData.hh> //to get secondary structure
#include <core/fragment/SecstructSRFD.hh> //to get secondary structure
#include <core/fragment/FragID_Iterator.hh>

#include <core/scoring/constraints/AtomPairConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

// AUTO-REMOVED #include <core/sequence/util.hh>
#include <core/id/SequenceMapping.hh>

#include <core/fragment/SecondaryStructure.hh>
#include <protocols/jumping/JumpSample.hh>


//numeric headers
#include <numeric/random/random.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <cstdlib>
#include <string>
#include <vector>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/constraints/NamedAtomPairConstraint.hh>
#include <core/util/SwitchResidueTypeSet.hh>



static thread_local basic::Tracer tr( "protocols.abinitio.Templates" );
using namespace core;
using namespace basic;
using namespace basic::options;



void protocols::abinitio::Templates::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	Template::register_options();
	PairingStatistics::register_options();
	option.add_relevant( templates::no_culling );
	option.add_relevant( templates::helix_pairings );
	option.add_relevant( templates::prefix );

	option.add_relevant( templates::cst::topN );
	option.add_relevant( templates::cst::wTopol );
	option.add_relevant( templates::cst::wExtern );

	option.add_relevant( templates::fragsteal::topN );
	option.add_relevant( templates::fragsteal::wTopol );
	option.add_relevant( templates::fragsteal::wExtern );

}


namespace protocols {
namespace abinitio {


using namespace core;
using namespace jumping;
using namespace basic::options;
using namespace basic::options::OptionKeys;
class StructureStore : public utility::pointer::ReferenceCount {
public:
	bool has( std::string const& name ) const {
		StructureMap::const_iterator iter = poses_.find( name );
		return ( iter != poses_.end() );
	};

	pose::PoseCOP operator[] ( std::string const& name ) const{
		StructureMap::const_iterator iter = poses_.find( name );
		if ( iter != poses_.end() ) {
			return (iter->second);
		} else return NULL;
	};

	void add( std::string const& file_name );

private:
	typedef std::map<std::string, pose::PoseCOP > StructureMap;
	StructureMap poses_;
};

void StructureStore::add( std::string const& file_name ) {
	if ( !has( file_name) ) {
		pose::PoseOP pose = new pose::Pose;

		//read structure
		core::import_pose::pose_from_pdb( *pose, file_name );

		// switch to centroid --- such that constraints get correct atom numbers assigned
		core::util::switch_to_residue_type_set( *pose, chemical::CENTROID );

		// set ss structure -- good for fragpicking
		pose::set_ss_from_phipsi( *pose );

		poses_[ file_name ]=pose;
	}
}

Templates::~Templates() {}

using namespace core;
using namespace fragment;
Templates::Templates( std::string const& config_file, pose::PoseCOP native ) :
	native_( native )
{
	utility::io::izstream in( config_file);
	tr.Info << "read homolog-template information from " << config_file << std::endl;
	if ( !in ) {
		utility_exit_with_message("ERROR:: Unable to open template-config file: "+config_file);
	}
	StructureStore pose_store;
	std::string line;
	good_ = true;
	while ( getline( in, line) ) {
		std::istringstream line_stream( line );
		std::string name, pdb, align;
		int offset; Real score;
		line_stream >> name >> pdb >> align >> offset >> score;
		TemplateOP theTemplate;
		if ( !line_stream.fail() && ( pdb.size() >= 1 ) ) {
			if ( option[ templates::prefix ].user() ) {
				std::string const prefix( option[ templates::prefix ]() );
				pdb = prefix+"/"+pdb;
				align = prefix+"/"+align;
			}
			tr.Info << "add Template " << name << " ..." << std::endl;
			if ( !pose_store.has( pdb ) ) pose_store.add( pdb );
			tr.Info << "template  " << name << " read template structure " << pdb << " with offset " << offset << std::endl;

			theTemplate = new Template( name, pose_store[pdb], align, offset, score );
			if( !theTemplate->is_good() ){
				good_ = false;
				continue;
			}
			if ( target_sequence_.size() ) {
				if ( target_sequence_ != theTemplate->query_sequence() ) {
					tr.Warning << "[WARNING] the query sequence " << theTemplate->query_sequence() << " is different than previous " << std::endl;
				}
			} else target_sequence_ = theTemplate->query_sequence();
			templates_.insert( TemplateMap::value_type( name, theTemplate ) );
		}

		std::string tag;
		while ( line_stream >> tag ) {
			if ( tag == "C" || tag == "SC" ) cull_list_.push_back( theTemplate );
			if ( tag == "F" ) fragpick_list_.push_back( theTemplate );
			if ( tag == "H" ) helixjump_pick_list_.push_back( theTemplate );
			if ( tag == "CST" || tag == "SCST" ) {
				std::string cst_file;
				line_stream >> cst_file;
				if ( option[ templates::prefix ].user() ) {
					std::string const prefix( option[ templates::prefix ]() );
					cst_file = prefix + "/" + cst_file;
				}
				if ( !line_stream.fail() ) {
					tr.Info << "read constraints " << cst_file << "..." << std::endl;
					theTemplate->read_constraints( cst_file );
				}
			}
		}
	}
	PairingStatisticsOP strand_stats = new PairingStatistics( *this );
	if ( native_ ) strand_stats->set_native_topology( core::scoring::dssp::StrandPairingSet( *native_ ) );
	strand_stats_ = strand_stats;

	tr.Info << "statistics of pairings: \n\n\n" << *strand_stats_ << std::endl;

	//copy topology score into templates:
	for ( Size i = 1; i <= strand_stats->nr_models(); i++ ) {
		templates_[ strand_stats->ranked_model( i ) ]->topology_score( strand_stats->weight( i ) );
	}

	if ( option[ templates::fragsteal::topN ].user() ) {
		scored_fragpick_list( fragpick_list_ );
		helixjump_pick_list_ = fragpick_list_;
		tr.Info << "helixjump_list \n";
		for ( TemplateList::const_iterator it = helixjump_pick_list_.begin(), eit = helixjump_pick_list_.end();
					it != eit; ++it ) {
			tr.Info << (*it)->name() << "\n";
		}
		tr.Info << std::endl;
	}
}

void Templates::set_native( core::pose::PoseCOP native ) {
	native_ = native;
}


void
Templates::get_cst_list( TemplateList& cst_list, TemplateList& cull_list ) const {
	cst_list.clear();

	bool bScoreFilter = false;
	Real wTopol( 0.0 ), wExtern( 0.0 );
	if ( option[ templates::cst::topN ].user() ) {
		bScoreFilter = true;
		wTopol = option[ templates::cst::wTopol ];
		wExtern = option[ templates::cst::wExtern ];
	}

	// first get list of Templates with constraints
	for ( TemplateMap::const_iterator it=templates_.begin(),
					eit = templates_.end(); it!=eit; ++it ) {
		// if template has constraints
		TemplateCOP aTemplate( it->second );
		if ( aTemplate->has_constraints() ) {
			cst_list.push_back( aTemplate );
		}
	}
	if ( !bScoreFilter ) return;

	_get_scored_list( cst_list, option[ templates::cst::topN ], wTopol, wExtern );
	cull_list = cst_list;
}


void
Templates::scored_fragpick_list( TemplateList& frag_list ) const {
	frag_list.clear();

	//bool bScoreFilter = false;
	Real wTopol( 0.0 ), wExtern( 0.0 );
	if ( option[ templates::fragsteal::topN ].user() ) {
		//bScoreFilter = true;  // set but never used ~Labonte
		wTopol = option[ templates::fragsteal::wTopol ];
		wExtern = option[ templates::fragsteal::wExtern ];
	}

	for ( TemplateMap::const_iterator it=templates_.begin(),
					eit = templates_.end(); it!=eit; ++it ) {
		frag_list.push_back( it->second );
	}
	_get_scored_list( frag_list, option[ templates::fragsteal::topN ], wTopol, wExtern );
}



void Templates::_get_scored_list( TemplateList& cst_list, Size topN, Real wTopol, Real wExtern) const {

	Real sum_extern = 0;
	Real sum_extern2 = 0;
	Real sum_topol = 0;
	Real sum_topol2 = 0;
	Size n = 0;

	bool bScoreFilter = true;

	// first get list of Templates with constraints
	for ( TemplateMap::const_iterator it=templates_.begin(),
					eit = templates_.end(); it!=eit; ++it ) {
		TemplateCOP aTemplate( it->second );
		if ( bScoreFilter ) {
			sum_extern += aTemplate->external_score();
			sum_extern2 += aTemplate->external_score()*aTemplate->external_score();
			sum_topol += aTemplate->topology_score();
			sum_topol2 += aTemplate->topology_score()*aTemplate->topology_score();
			n++;
		}
	}
	if ( n == 0 ) return;
	if ( !bScoreFilter ) return;
	// compute std-dev and mean:
	Real mean_topol = sum_topol / n;
	Real mean_extern = sum_extern / n;
	Real std_topol = std::sqrt( sum_topol2 / n - mean_topol*mean_topol );
	Real std_extern = std::sqrt( sum_extern2 / n - mean_extern*mean_extern );

	std::list< std::pair< core::Real, TemplateCOP > > weight_list;
	for ( TemplateList::const_iterator it = cst_list.begin(), eit = cst_list.end();
				it != eit; ++it ) {
		Real score = ( (*it)->topology_score() - mean_topol ) / std_topol * wTopol
			+ ( (*it)->external_score() - mean_extern ) / std_extern * wExtern;
		weight_list.push_back( std::make_pair( score, *it ) );
	}
	weight_list.sort();
	weight_list.reverse();
	cst_list.clear();
	std::list< std::pair< core::Real, TemplateCOP > >::const_iterator iter = weight_list.begin();
	for ( Size i = 1; i <= topN && iter != weight_list.end(); i++, iter++ ) {
		cst_list.push_back( iter->second );
	}
}


Size
Templates::pick_frags( FragSet& frag_set, core::fragment::FragDataCOP frag_type, Size ncopies /* default 1 */) const {

	Size nframes = target_total_residue() - frag_type->size() + 1;
	FrameList frames;
	Size total( 0 );
	tr.Info << "pick frames for target position 1 -> " << nframes << " from templates"  << std::endl;
	for ( Size pos =1; pos<=nframes; pos ++ ) {
		FrameOP frame = new Frame( pos, frag_type );
		frames.push_back( frame );
	}

	for ( TemplateList::const_iterator it=fragpick_list_.begin(),
					eit = fragpick_list_.end(); it!=eit; ++it ) {
		tr.Info << "pick from template " << (*it)->name() << std::endl;
		Size nr_frags = (*it)->steal_frags( frames, frag_set, ncopies );
		tr.Info << "found " << nr_frags << " new fragments " << std::endl;
		total+=nr_frags;
	}

	return total;
}


FragSetOP
Templates::pick_frags( FragSetOP frag_set, core::fragment::FragDataCOP frag_type, Size min_nr_frags, Size ncopies /* default 1 */ ) const {

	ConstantLengthFragSet template_frags;
	Size total( 0 );
	total = pick_frags( template_frags, frag_type, ncopies );

	// template_frags contains only homolog fragments, .. needs filling up in gappy regions
	FragSetOP merged_frags = frag_set->empty_clone();

	Size total_fill( 0 );
	//merge fragments:
	for ( Size pos = 1; pos<=target_total_residue(); pos++ ) {
		FrameList template_frames;
		template_frags.frames( pos, template_frames );
		merged_frags->add( template_frames );
		Size nr_frags( template_frames.flat_size() );
		tr.Info << nr_frags << " fragments at pos " << pos << ". required: " << min_nr_frags << std::endl;
		if ( nr_frags < min_nr_frags ) {
			Size nr_fill ( min_nr_frags - nr_frags );
			FrameList standard_frames;
			frag_set->frames( pos, standard_frames );
			if ( standard_frames.size() ) {
				tr.Info << "attempt to fill up with " << nr_fill << " frags at position " << pos << " ... ";
				for ( FragID_Iterator it = standard_frames.begin(), eit = standard_frames.end();
							it != eit && nr_fill; ++it, --nr_fill ) {
					merged_frags->add( *it );
					++total_fill;
				}
				if ( nr_fill ) {
					tr.Info << nr_fill << " fragments short " << std::endl;
				} else {
					tr.Info << "succeeded! " << std::endl;
				}
			} // standard frags present
		} //fill up
	}
	tr.Info << "found " << total << " fragments from homologs. supplemented by " << total_fill << " frags from standard library " << std::endl;
	return merged_frags;
}
//   Size total( 0 );
//   for ( TemplateMap::const_iterator it=templates_.begin(),
// 	  eit = templates_.end(); it!=eit; ++it ) {
//     total += it->second->pick_frags( frag_set, frag_type );
//   }
//  return total;
//}

Size Templates::pick_large_frags(
			core::fragment::FragSet& frag_set,
			core::fragment::SingleResidueFragDataOP frag_type,
			core::Size ncopies /*default = 1*/
) const {
	Size total( 0 );
	for ( TemplateList::const_iterator it=fragpick_list_.begin(),
					eit = fragpick_list_.end(); it!=eit; ++it ) {
		tr.Info << "pick large frag from template " << (*it)->name() << std::endl;
		Size nr_frags = (*it)->pick_large_frags( frag_set, frag_type, ncopies );
		tr.Info << "found " << nr_frags << " new fragments " << std::endl;
		total+=nr_frags;
	}
	return total;
}


void
Templates::read_pairings( std::string const& filename, core::scoring::dssp::PairingsList& pairings ) const {
	pairings.clear();

	utility::io::izstream in( filename );
	if ( !in ) {
		utility_exit_with_message("ERROR Unable to open pairings file "+filename);
	}

	std::string pdb;
	in >> pdb;

	core::scoring::dssp::PairingsList raw_pairings;
	read_pairing_list( in, raw_pairings );
	tr.Debug << " read pairings for template " << pdb << "\n" << raw_pairings << std::endl;
	const_iterator iter = templates_.find( pdb );
	if ( iter == templates_.end() ) {
		utility_exit_with_message("unrecognized template name "+pdb+" --- this name has to be in template:config file");
	}
	iter->second->map_pairings2target( raw_pairings, pairings );
	tr.Debug << " mapped pairings for target\n " << pairings_ << std::endl;
}

TemplateJumpSetupOP Templates::create_jump_def( core::fragment::SecondaryStructureCOP ss_def ) const {
	if ( !ss_def ) {
		using namespace fragment;
		tr.Info << "TemplateJumpSetup will be initialized with secondary structure from homologs " << std::endl;
		ConstantLengthFragSet fragset;
		pick_frags( fragset, new FragData( new SecstructSRFD, 1 ) ); //for ss-structure 1mers are enough
		ss_def = new core::fragment::SecondaryStructure( fragset, target_total_residue() );
	}
// 	utility::io::ozstream dump("ss_def_for_jumps");
// 	for ( Size i = 1; i<=ss_def->total_residue(); i++ ) {
// 		dump << i << " " << ss_def->loop_fraction()(i) << std::endl;
// 	}
	core::scoring::dssp::PairingsList helix_pairings;
	if ( option[ templates::helix_pairings ].user() ) read_pairings( option[ templates::helix_pairings ], helix_pairings );
	return new TemplateJumpSetup( this, ss_def, strand_stats_, helix_pairings );
}



void
Templates::add_target_constraints( scoring::constraints::ConstraintSetOP cstset, pose::Pose const& pose ) const {
	using namespace scoring::constraints;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	typedef Template::NamedAtomPairConstraintList NamedAtomPairConstraintList;
	typedef Template::AtomPairConstraintList AtomPairConstraintList;
	AtomPairConstraintList full_list;
	// take constraint sets from each template
	TemplateList cst_list;
	TemplateList cull_list = cull_list_;

	get_cst_list( cst_list, cull_list ); // evaluates options to see if score-ranking is used...
	tr.Info << "pick constraints from " << cst_list.size() << " models " << std::endl;

	for ( TemplateList::const_iterator it=cst_list.begin(),
					eit = cst_list.end(); it!=eit; ++it ) {
		tr.Info << "pick constraints from template " << (*it)->name() << std::endl;

		// if template has constraints
		Template const& aTemplate( **it );
		// map them to target sequence
		NamedAtomPairConstraintList new_constraints;
		aTemplate.map2target( aTemplate.constraints(), new_constraints );

		tr.Info << "have " << new_constraints.size() << " constraints; start culling... " << std::endl;

		if ( !option[ templates::no_culling ] ) {
			//throw out all constraints that violate any template structure
			for ( TemplateList::const_iterator it=cull_list.begin(),
							eit = cull_list.end(); it!=eit; ++it ) {
				NamedAtomPairConstraintList culled_constraints;
				tr.Info << "cull with template " << (*it)->name() << std::endl;
				(*it)->cull_violators( new_constraints, culled_constraints );
				new_constraints = culled_constraints;
				tr.Info << (*it)->name() << " leaves " << new_constraints.size() << " unviolated " << std::endl;
			}
		}

		// add them to full_list if they are not out of bounds
		for ( NamedAtomPairConstraintList::const_iterator it = new_constraints.begin(),
						eit = new_constraints.end(); it!=eit; ++it ) {
				AtomPairConstraintOP cst = (*it)->mapto( pose );
				if ( cst ) full_list.push_back( cst );
			}
	}
	cstset->add_constraints( full_list );
	// run them over all other templates and throw out violated constraints
	// Todo: consolidate: redundant constraints --> single constraint with higher weight
}


} //abinitio
} //protocols
