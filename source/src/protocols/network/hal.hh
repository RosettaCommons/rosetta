// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/network/hal.hh
/// @brief: HAL implementaion
///
/// @author Sergey Lyskov

#ifdef ZEROMQ

#pragma once

#include <core/pose/Pose.hh>
#include <protocols/network/util.hh>
#include <utility/json_utilities.hh>
#include <utility/exit.hh>

#include <json.hpp>

#include <functional>
#include <string>

namespace protocols {
namespace network {

struct CommandLineArguments
{
	int argc;
	char * * argv;
};


using SpecificationCallBack = std::function< json(void) >;
using ExecutionerCallBack   = std::function< json(json const &) >;

/// start HAL listener
/// use main function `argc` and `argv` if you want auto-restart/abort functionality to work properly
void hal(SpecificationCallBack const &, ExecutionerCallBack const &, CommandLineArguments const &args);


struct HA;
using HAs = std::initializer_list< std::pair<const std::string, HA> >;

void hal(HAs a, CommandLineArguments const &args);


/// auxiliary structure to specify HAL action
struct HA
{
	template <typename R, typename... As>
	HA( std::function< R(core::pose::Pose &, As...) >, json::array_t const &);

	template <typename R, typename... As>
	HA( R (*)(core::pose::Pose &, As...), json::array_t const &);


	template <typename... As>
	HA( std::function< core::pose::Pose (As...) >, json::array_t const &);

	template <typename... As>
	HA( core::pose::Pose (*)(As...), json::array_t const &);

	json specification;
	std::function< json(json const &) > execute;
};


namespace {


/// HAL_Type info, map from C++ type to various auxiliary values
template <typename T, typename enabled = void> struct HAL_Type;

template <>
struct HAL_Type<bool>
{
	using value_type = bool;
	static auto constexpr name = _t_boolean_;
};

template <typename Integer>
struct HAL_Type<Integer, typename std::enable_if<std::is_integral<Integer>::value>::type>
{
	using value_type = int;
	static auto constexpr name = _t_integer_;
};

template <typename Float>
struct HAL_Type<Float, typename std::enable_if<std::is_floating_point<Float>::value>::type>
{
	using value_type = double;
	static auto constexpr name = _t_float_;
};

template <>
struct HAL_Type<std::string>
{
	using value_type = std::string;
	static auto constexpr name = _t_string_;
};



template <typename...>
struct TypeGroup;


template <int, typename R, typename ...Groups>
struct ApplyTuple;


template <
	typename R,
	template <typename...> class G,
	typename...Fs,
	typename...Ts,
	typename...As
>

struct ApplyTuple< 0, R, G<Fs...>, G<Ts...>, G<As...> >
{
	static R apply( std::function< R(Fs...) > fun, std::tuple<Ts...> const &, As... as) {
		return fun(as...);
	}
};


template <
	int n,
	typename R,
	template <typename...> class G,
	typename...Fs,
	typename...Ts,
	typename...As
>

struct ApplyTuple<n, R, G<Fs...>, G<Ts...>, G<As...> >
{
	static R apply( std::function< R(Fs...) > fun, std::tuple<Ts...> const &tpl, As... args) {
		return ApplyTuple<
			n-1,
			R,
			TypeGroup<Fs...>,
			TypeGroup<Ts...>,
			TypeGroup<As..., typename std::tuple_element<sizeof...(Ts) - n, std::tuple<Ts...> >::type>
		>::apply(fun, tpl, args..., std::get< sizeof...(Ts) - n >(tpl) );
	}
};



template <int N>
using ArgumentNames = std::array<std::string, N>;


/// extract_arguments from json and pack then into tuple
template <int N, int I, typename... Ts>
struct ExtractArguments {};


template <>
struct ExtractArguments<0, 0>
{
	static std::tuple<> arguments(json const &, ArgumentNames<0> const &) {
		return std::tuple<>();
	}
};

template <int N, int I, typename T>
struct ExtractArguments<N, I, T>
{
	static std::tuple<T> arguments(json const &command, ArgumentNames<N> const &args) {
		static_assert( I < N  and  I >= 0, "extract_arguments: index is out of range!");

		//std::cout << "N=" << N << " I=" << I << std::endl;
		//return std::tuple<T>( SpecificationHelper<T>::extract_argument(command, args[I]) );

		typename HAL_Type<T>::value_type v{};
		utility::extract_value_if_present(command, args[I], v);
		//std::cout << "SpecificationHelper::extract_argument key:" << args[I] << " value:" << v << std::endl;
		return std::tuple<T>(v);
	}
};

template <int N, int I, typename T, typename... Ts>
struct ExtractArguments<N, I, T, Ts...>
{
	static std::tuple<T, Ts...> arguments(json const &j, ArgumentNames<N> const &args)
	{
		return std::tuple_cat(ExtractArguments<N, I, T>::arguments(j, args), ExtractArguments<N, I + 1, Ts...>::arguments(j, args) );
	}
};


template <typename Result, typename... As>
std::tuple<json, ArgumentNames<sizeof...(As)> > generate_specification( std::function< Result (core::pose::Pose &, As...) >, json::array_t const &args)//std::array<json, N> const &args)
{
	//static_assert( N <= sizeof...(As), "Too many arguments specified!");
	if( args.size() > sizeof...(As) ) {
		utility_exit_with_message("generate_specification: too many arguments provided!");
		return json();
	}

	ArgumentNames<sizeof...(As)> type_names = { {HAL_Type< typename std::decay<As>::type >::name...} };

	auto s = json::object();
	s[_f_pose_] = { {_f_type_, _t_pose_}, };

	std::array<std::string, sizeof...(As) > names;

	for(uint i=0; i < sizeof...(As); ++i) {
		std::string name = std::string("arg-") + std::to_string(i) + "(" + type_names.begin()[i] + ")";
		std::string type = type_names.begin()[i];
		json spec = { {_f_type_, type}, };

		if( i < args.size() ) {
			json const & src = args[i];

			if( src.is_string() ) name = args[i];
			else if( src.is_object() ) {
				utility::extract_value_if_present(src, _f_name_, name);
				utility::extract_value_if_present(src, _f_type_, type);
				spec = src;
				spec[_f_type_] = type;
			}
		}

		s[name] = spec;

		names[i] = name;
	}
	return make_tuple(s, names);
}


template <typename... As>
std::tuple<json, ArgumentNames<sizeof...(As)> > generate_specification( std::function< core::pose::Pose(As...) >, json::array_t const &args)//std::array<json, N> const &args)
{
	//static_assert( N <= sizeof...(As), "Too many arguments specified!");
	if( args.size() > sizeof...(As) ) {
		utility_exit_with_message("generate_specification: too many arguments provided!");
		return json();
	}

	ArgumentNames<sizeof...(As)> type_names = { {HAL_Type< typename std::remove_const< typename std::remove_reference<As>::type >::type >::name...} };

	auto s = json::object();

	std::array<std::string, sizeof...(As) > names;

	for(uint i=0; i < sizeof...(As); ++i) {
		std::string name = std::string("arg-") + std::to_string(i) + "(" + type_names.begin()[i] + ")";
		std::string type = type_names.begin()[i];
		json spec = { {_f_type_, type}, };

		if( i < args.size() ) {
			json const & src = args[i];

			if( src.is_string() ) name = args[i];
			else if( src.is_object() ) {
				utility::extract_value_if_present(src, _f_name_, name);
				utility::extract_value_if_present(src, _f_type_, type);
				spec = src;
				spec[_f_type_] = type;
			}
		}

		s[name] = spec;

		names[i] = name;
	}
	return make_tuple(s, names);
}

} // namespace


template <typename R, typename... As>
HA::HA( std::function< R (core::pose::Pose &, As...) > f, json::array_t const &args_specification)
{
	std::array<std::string, sizeof...(As)> names;

	std::tie(specification, names) = generate_specification(f, args_specification);

	//std::cout << "Specification: " << specification << std::endl;

	execute = [names, f](json const &args) {
				  std::string pose_data;
				  utility::extract_value_if_present(args, _f_pose_, pose_data);
				  core::pose::PoseOP pose = protocols::network::bytes_to_pose(pose_data);
				  if(pose) {

					  //{ std::cout << "names:"; for(auto const & n : names) std::cout << n << ' '; std::cout << std::endl; }

				  	  auto tpl = ExtractArguments<sizeof...(As), 0, typename std::decay<As>::type...>::arguments(args, names);
				  	  //apply_function<R, TypeGroup<As...>, TypeGroup<typename std::decay<As>::type...> >(f, std::tuple_cat( std::make_tuple(std::ref(*pose)), tpl) );
					  ApplyTuple<
						  sizeof...(As) + 1,
						  void,
						  TypeGroup<core::pose::Pose &, As...>,
						  TypeGroup<core::pose::Pose &, typename std::decay<As>::type...>,
						  TypeGroup<>
					  >::apply(f, std::tuple_cat( std::make_tuple(std::ref(*pose)), tpl) );

					  auto pose_binary = protocols::network::pose_to_bytes(*pose);
					  json result;
					  result[_f_pose_] = pose_binary;
					  return result;
				  }
				  else {
				  	  std::cout << "hal::HA: could not find Pose in given command, ignoring execute request..." << std::endl;
					  return json();
				  }
			  };
}

template <typename R, typename... As>
HA::HA( R (*f)(core::pose::Pose &, As...), json::array_t const &a) : HA( std::function< R(core::pose::Pose &, As...) >(f), a) {}


template <typename... As>
HA::HA( std::function< core::pose::Pose(As...)> f, json::array_t const &args_specification)
{
	std::array<std::string, sizeof...(As)> names;

	std::tie(specification, names) = generate_specification(f, args_specification);
	//std::cout << "Specification: " << specification << std::endl;

	//execute = [](json const &args) { return json(); };
	execute = [names, f](json const &args) {

				  //{ std::cout << "names:"; for(auto const & n : names) std::cout << n << ' '; std::cout << std::endl; }

				  auto tpl = ExtractArguments<sizeof...(As), 0, typename std::decay<As>::type...>::arguments(args, names);
				  //apply_function<R, TypeGroup<As...>, TypeGroup<typename std::decay<As>::type...> >(f, std::tuple_cat( std::make_tuple(std::ref(*pose)), tpl) );
				  core::pose::Pose pose = ApplyTuple<
					  sizeof...(As),
					  core::pose::Pose,
					  TypeGroup<As...>,
					  TypeGroup<typename std::decay<As>::type...>,
					  TypeGroup<>
				  >::apply(f, tpl);

				  auto pose_binary = protocols::network::pose_to_bytes(pose);
				  json result;
				  result[_f_pose_] = pose_binary;
				  return result;
			  };
}

template <typename... As>
HA::HA( core::pose::Pose (*f)(As...), json::array_t const &a) : HA( std::function< core::pose::Pose(As...) >(f), a) {}



} // namespace network
} // namespace protocols

#endif // ZEROMQ
