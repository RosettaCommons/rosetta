// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   rosetta/benchmark/InstantiateModules.bench.cc
///
/// @brief  Benchmark the instantiation of movers, filters, etc.  These should be fast to
/// instantiate, and should load any large data from disk lazily and once.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).

#ifndef INCLUDED_apps_benchmark_InstantiateModulesBenchmark_bench_hh
#define INCLUDED_apps_benchmark_InstantiateModulesBenchmark_bench_hh


#include <apps/benchmark/performance/performance_benchmark.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pack/task/operation/TaskOperationCreator.hh>
#include <core/pack/task/operation/ResLvlTaskOperationFactory.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/ResidueSelectorCreator.hh>
#include <core/simple_metrics/SimpleMetricFactory.hh>
#include <core/simple_metrics/SimpleMetricCreator.hh>
#include <core/pack/palette/PackerPaletteFactory.hh>
#include <core/pack/palette/PackerPaletteCreator.hh>

#include <protocols/moves/MoverFactory.hh>
#include <protocols/filters/FilterFactory.hh>
#include <protocols/moves/Mover.hh>

#include <chrono>
#include <map>

enum class ModuleType {
	MOVER = 1, //Keep first
	FILTER,
	TASK_OPERATION,
	RES_LEVEL_TASK_OPERATION,
	RESIDUE_SELECTOR,
	SIMPLE_METRIC,
	PACKER_PALETTE, //Keep second to last
	N_TYPES = PACKER_PALETTE //Keep last
};

class InstantiateModulesBenchmark : public PerformanceBenchmark
{

public:

	InstantiateModulesBenchmark(
		std::string const & name,
		ModuleType const module_type
	) :
		PerformanceBenchmark(name),
		module_type_(module_type)
	{}

public:

	void
	instantiate_movers(
		core::Real const scalefactor
	) const {
		using namespace protocols::moves;
		MoverFactory::MoverMap const & movermap( MoverFactory::get_instance()->mover_creator_map() );

		utility::vector1< core::Real > timings( movermap.size() );
		utility::vector1< std::string > names( movermap.size() );
		core::Size counter(0);

		core::Size repeats( static_cast<core::Size>( scalefactor ) );
		if ( repeats == 0 ) repeats = 1;

		for ( MoverFactory::MoverMap::const_iterator it( movermap.begin() ); it!=movermap.end(); ++it ) {
			++counter;
			std::chrono::time_point<std::chrono::high_resolution_clock> const starttime( std::chrono::high_resolution_clock::now() );
			for ( core::Size j(1); j<=repeats; ++j ) {
				protocols::moves::MoverOP mover( MoverFactory::get_instance()->newMover( it->first ) );
			}
			std::chrono::time_point<std::chrono::high_resolution_clock> const endtime( std::chrono::high_resolution_clock::now() );
			std::chrono::duration< double > diff( endtime - starttime );
			timings[counter] = static_cast< core::Real >( diff.count() );
			names[counter] = it->first;
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << "\nMOVER\tTIME_(s)\n";
			for ( core::Size i(1), imax(timings.size()); i<=imax; ++i ) {
				TR.Debug << names[i] << "\t" << timings[i] << "\n";
			}
			TR.Debug.flush();
		}
	}

	void
	instantiate_filters(
		core::Real const scalefactor
	) const {
		using namespace protocols::filters;
		FilterFactory::FilterMap const & filtermap( FilterFactory::get_instance()->filter_creator_map() );

		utility::vector1< core::Real > timings( filtermap.size() );
		utility::vector1< std::string > names( filtermap.size() );
		core::Size counter(0);

		core::Size repeats( static_cast<core::Size>( scalefactor ) );
		if ( repeats == 0 ) repeats = 1;

		for ( FilterFactory::FilterMap::const_iterator it( filtermap.begin() ); it!=filtermap.end(); ++it ) {
			++counter;
			std::chrono::time_point<std::chrono::high_resolution_clock> const starttime( std::chrono::high_resolution_clock::now() );
			for ( core::Size j(1); j<=repeats; ++j ) {
				FilterOP filter( FilterFactory::get_instance()->newFilter( it->first ) );
			}
			std::chrono::time_point<std::chrono::high_resolution_clock> const endtime( std::chrono::high_resolution_clock::now() );
			std::chrono::duration< double > diff( endtime - starttime );
			timings[counter] = static_cast< core::Real >( diff.count() );
			names[counter] = it->first;
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << "\nFILTER\tTIME_(s)\n";
			for ( core::Size i(1), imax(timings.size()); i<=imax; ++i ) {
				TR.Debug << names[i] << "\t" << timings[i] << "\n";
			}
			TR.Debug.flush();
		}
	}

	void
	instantiate_taskops(
		core::Real const scalefactor
	) const {
		using namespace core::pack::task::operation;
		TaskOperationFactory::TaskOperationCreatorMap const & taskopmap( TaskOperationFactory::get_instance()->creator_map() );

		utility::vector1< core::Real > timings( taskopmap.size() );
		utility::vector1< std::string > names( taskopmap.size() );
		core::Size counter(0);

		core::Size repeats( static_cast<core::Size>( scalefactor ) );
		if ( repeats == 0 ) repeats = 1;

		for ( TaskOperationFactory::TaskOperationCreatorMap::const_iterator it( taskopmap.begin() ); it!=taskopmap.end(); ++it ) {
			++counter;
			std::chrono::time_point<std::chrono::high_resolution_clock> const starttime( std::chrono::high_resolution_clock::now() );
			for ( core::Size j(1); j<=repeats; ++j ) {
				TaskOperationOP taskop( it->second->create_task_operation() );
			}
			std::chrono::time_point<std::chrono::high_resolution_clock> const endtime( std::chrono::high_resolution_clock::now() );
			std::chrono::duration< double > diff( endtime - starttime );
			timings[counter] = static_cast< core::Real >( diff.count() );
			names[counter] = it->first;
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << "\nTASKOP\tTIME_(s)\n";
			for ( core::Size i(1), imax(timings.size()); i<=imax; ++i ) {
				TR.Debug << names[i] << "\t" << timings[i] << "\n";
			}
			TR.Debug.flush();
		}
	}

	void
	instantiate_reslvl_taskops(
		core::Real const scalefactor
	) const {
		using namespace core::pack::task::operation;
		std::map< std::string, ResLvlTaskOperationCreatorOP > const & reslvl_taskopmap( ResLvlTaskOperationFactory::get_instance()->creator_map() );

		utility::vector1< core::Real > timings( reslvl_taskopmap.size() );
		utility::vector1< std::string > names( reslvl_taskopmap.size() );
		core::Size counter(0);

		core::Size repeats( static_cast<core::Size>( scalefactor ) );
		if ( repeats == 0 ) repeats = 1;

		for ( std::map< std::string, ResLvlTaskOperationCreatorOP >::const_iterator it( reslvl_taskopmap.begin() ); it!=reslvl_taskopmap.end(); ++it ) {
			++counter;
			std::chrono::time_point<std::chrono::high_resolution_clock> const starttime( std::chrono::high_resolution_clock::now() );
			for ( core::Size j(1); j<=repeats; ++j ) {
				ResLvlTaskOperationOP reslvl_taskop( ResLvlTaskOperationFactory::get_instance()->newRLTO( it->first ) );
			}
			std::chrono::time_point<std::chrono::high_resolution_clock> const endtime( std::chrono::high_resolution_clock::now() );
			std::chrono::duration< double > diff( endtime - starttime );
			timings[counter] = static_cast< core::Real >( diff.count() );
			names[counter] = it->first;
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << "\nRESLVL_TASKOP\tTIME_(s)\n";
			for ( core::Size i(1), imax(timings.size()); i<=imax; ++i ) {
				TR.Debug << names[i] << "\t" << timings[i] << "\n";
			}
			TR.Debug.flush();
		}
	}

	void
	instantiate_res_selectors(
		core::Real const scalefactor
	) const {
		using namespace core::select::residue_selector;
		std::map< std::string, ResidueSelectorCreatorOP > const & selectormap( ResidueSelectorFactory::get_instance()->creator_map() );

		utility::vector1< core::Real > timings( selectormap.size() );
		utility::vector1< std::string > names( selectormap.size() );
		core::Size counter(0);

		core::Size repeats( static_cast<core::Size>( scalefactor ) );
		if ( repeats == 0 ) repeats = 1;

		for ( std::map< std::string, ResidueSelectorCreatorOP >::const_iterator it( selectormap.begin() ); it!=selectormap.end(); ++it ) {
			++counter;
			std::chrono::time_point<std::chrono::high_resolution_clock> const starttime( std::chrono::high_resolution_clock::now() );
			for ( core::Size j(1); j<=repeats; ++j ) {
				ResidueSelectorOP selector( it->second->create_residue_selector() );
			}
			std::chrono::time_point<std::chrono::high_resolution_clock> const endtime( std::chrono::high_resolution_clock::now() );
			std::chrono::duration< double > diff( endtime - starttime );
			timings[counter] = static_cast< core::Real >( diff.count() );
			names[counter] = it->first;
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << "\nRES_SELECTOR\tTIME_(s)\n";
			for ( core::Size i(1), imax(timings.size()); i<=imax; ++i ) {
				TR.Debug << names[i] << "\t" << timings[i] << "\n";
			}
			TR.Debug.flush();
		}
	}

	void
	instantiate_simple_metrics(
		core::Real const scalefactor
	) const {
		using namespace core::simple_metrics;
		SimpleMetricFactory::CreatorMap const & simple_metric_map( SimpleMetricFactory::get_instance()->creator_map() );

		utility::vector1< core::Real > timings( simple_metric_map.size() );
		utility::vector1< std::string > names( simple_metric_map.size() );
		core::Size counter(0);

		core::Size repeats( static_cast<core::Size>( scalefactor ) );
		if ( repeats == 0 ) repeats = 1;

		for ( SimpleMetricFactory::CreatorMap::const_iterator it( simple_metric_map.begin() ); it!=simple_metric_map.end(); ++it ) {
			++counter;
			std::chrono::time_point<std::chrono::high_resolution_clock> const starttime( std::chrono::high_resolution_clock::now() );
			for ( core::Size j(1); j<=repeats; ++j ) {
				SimpleMetricOP metric( it->second->create_simple_metric() );
			}
			std::chrono::time_point<std::chrono::high_resolution_clock> const endtime( std::chrono::high_resolution_clock::now() );
			std::chrono::duration< double > diff( endtime - starttime );
			timings[counter] = static_cast< core::Real >( diff.count() );
			names[counter] = it->first;
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << "\nSIMPLE_METRIC\tTIME_(s)\n";
			for ( core::Size i(1), imax(timings.size()); i<=imax; ++i ) {
				TR.Debug << names[i] << "\t" << timings[i] << "\n";
			}
			TR.Debug.flush();
		}
	}

	void
	instantiate_packer_palettes(
		core::Real const scalefactor
	) const {
		using namespace core::pack::palette;
		PackerPaletteFactory::PackerPaletteCreatorMap const & packer_palette_map( PackerPaletteFactory::get_instance()->packer_palette_creator_map() );

		utility::vector1< core::Real > timings( packer_palette_map.size() );
		utility::vector1< std::string > names( packer_palette_map.size() );
		core::Size counter(0);

		core::Size repeats( static_cast<core::Size>( scalefactor ) );
		if ( repeats == 0 ) repeats = 1;

		for ( PackerPaletteFactory::PackerPaletteCreatorMap::const_iterator it( packer_palette_map.begin() ); it!=packer_palette_map.end(); ++it ) {
			++counter;
			std::chrono::time_point<std::chrono::high_resolution_clock> const starttime( std::chrono::high_resolution_clock::now() );
			for ( core::Size j(1); j<=repeats; ++j ) {
				PackerPaletteOP palette( it->second->create_packer_palette() );
			}
			std::chrono::time_point<std::chrono::high_resolution_clock> const endtime( std::chrono::high_resolution_clock::now() );
			std::chrono::duration< double > diff( endtime - starttime );
			timings[counter] = static_cast< core::Real >( diff.count() );
			names[counter] = it->first;
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << "\nPACKER_PALETTE\tTIME_(s)\n";
			for ( core::Size i(1), imax(timings.size()); i<=imax; ++i ) {
				TR.Debug << names[i] << "\t" << timings[i] << "\n";
			}
			TR.Debug.flush();
		}
	}

public:

	void setUp() override {}

	void run(core::Real scaleFactor) override {
		switch( module_type_ ) {
		case ModuleType::MOVER :
			instantiate_movers(scaleFactor);
			break;
		case ModuleType::FILTER :
			instantiate_filters(scaleFactor);
			break;
		case ModuleType::TASK_OPERATION :
			instantiate_taskops(scaleFactor);
			break;
		case ModuleType::RES_LEVEL_TASK_OPERATION :
			instantiate_reslvl_taskops(scaleFactor);
			break;
		case ModuleType::RESIDUE_SELECTOR :
			instantiate_res_selectors(scaleFactor);
			break;
		case ModuleType::SIMPLE_METRIC :
			instantiate_simple_metrics(scaleFactor);
			break;
		case ModuleType::PACKER_PALETTE :
			instantiate_packer_palettes(scaleFactor);
			break;
		}
	}

	void tearDown() override {}



private:

	ModuleType module_type_ = ModuleType::MOVER;

};

#endif // include guard
