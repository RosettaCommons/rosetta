#ifndef FEM_MAIN_HPP
#define FEM_MAIN_HPP

#include <fem/intrinsics_extra.hpp>
#include <fem/stop.hpp>
#include <fem/utils/int_types.hpp>
#include <fem/utils/string.hpp>
#include <vector>
#include <cstdio>

#if defined(_MSC_VER)
# include <fcntl.h>
# include <io.h>
#endif

namespace fem {

  inline
  bool
  check_fem_utils_int_types()
  {
    bool result = true;
#define FABLE_LOC(bits, expected_sz) \
    { \
      int sz = static_cast<int>(sizeof(fem::utils::int##bits##_t)); \
      if (sz != expected_sz) { \
        std::fprintf(stderr, \
          "FATAL: sizeof(fem::utils::int%d_t) is %d but should be %d.\n", \
            bits, sz, expected_sz); \
        result = false; \
      } \
    }
    FABLE_LOC( 8, 1)
    FABLE_LOC(16, 2)
    FABLE_LOC(32, 4)
    FABLE_LOC(64, 8)
#undef FABLE_LOC
    if (!result) {
      std::fprintf(stderr,
        "NOTE: fem/utils/int_types.hpp"
        " needs to be adjusted for this platform.\n");
    }
    return result;
  }

  inline
  int
  main_with_catch(
    int argc,
    char const* argv[],
    void (*callable)(int argc, char const* argv[]))
  {
    user_plus_system_time();
    if (!check_fem_utils_int_types()) {
      return 255;
    }
#if defined(_MSC_VER)
    if (_setmode(_fileno(stdin), _O_BINARY) == -1) {
      std::cerr
        << "FATAL: error switching stdin to binary mode ("
        << __FILE__ << ", line " << __LINE__ << ")" << std::endl;
      return 255;
    }
#endif
    try {
      callable(argc, argv);
    }
    catch (fem::stop_info const& info) {
      std::fflush(stdout);
      std::fprintf(stderr, "%s\n", info.what());
      std::fflush(stderr);
    }
    catch (std::exception const& e) {
      std::fflush(stdout);
      char const* what = e.what();
      if (what == 0) what = "null";
      std::fprintf(stderr, "std::exception what(): %s\n", what);
      std::fflush(stderr);
      return 1;
    }
    catch (...) {
      std::fflush(stdout);
      std::fprintf(stderr, "Terminated by unknown C++ exception.\n");
      std::fflush(stderr);
      return 2;
    }
    return 0;
  }

  static const std::string
    dynamic_parameters_option(
      "--fem-dynamic-parameters");

  struct command_line_arguments
  {
    std::vector<std::string> buffer;
    std::vector<std::string> dynamic_parameters_fields;

    command_line_arguments() {}

    command_line_arguments(
      int argc,
      char const* argv[])
    {
      for(int i=0;i<argc;i++) {
        char const* arg = argv[i];
        if (utils::starts_with(
              arg,
              /*start*/ 0,
              /*stop*/ dynamic_parameters_option.size(),
              dynamic_parameters_option.c_str())) {
          size_t j = dynamic_parameters_option.size();
          if (arg[j] == '=') j++;
          utils::split_comma_separated(dynamic_parameters_fields, arg+j);
        }
        else {
          buffer.push_back(std::string(arg));
        }
      }
    }
  };

  struct dynamic_parameters_from
  {
    command_line_arguments const& command_line_args;
    int i_fld;

    dynamic_parameters_from(
      command_line_arguments const& command_line_args_,
      int max_number_of_flds)
    :
      command_line_args(command_line_args_),
      i_fld(0)
    {
      int n_flds = static_cast<int>(
        command_line_args.dynamic_parameters_fields.size());
      if (n_flds > max_number_of_flds) {
        std::ostringstream o;
        o << "Too many " << dynamic_parameters_option << " fields (given: "
          << n_flds
          << ", max. expected: "
          << max_number_of_flds
          << ")";
        throw std::runtime_error(o.str());
      }
    }

    template <typename T>
    dynamic_parameters_from&
    reset_if_given(
      T& value)
    {
      int n_flds = static_cast<int>(
        command_line_args.dynamic_parameters_fields.size());
      if (i_fld < n_flds) {
        std::string const&
          fld = command_line_args.dynamic_parameters_fields[i_fld];
        std::istringstream i(fld);
        i >> value;
        if (i.fail()) {
          std::ostringstream o;
          o << "Invalid " << dynamic_parameters_option << " field (field "
            << (i_fld+1)
            << "): \""
            << fld
            << "\"";
          throw std::runtime_error(o.str());
        }
        i_fld++;
      }
      return *this;
    }
  };

  template <typename D>
  struct dynamic_parameters_capsule
  {
    D dynamic_params;

    dynamic_parameters_capsule(
      command_line_arguments const& command_line_args)
    :
      dynamic_params(command_line_args)
    {}
  };

} // namespace fem

#endif // GUARD
