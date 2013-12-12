#ifndef FEM_INTRINSICS_EXTRA_HPP
#define FEM_INTRINSICS_EXTRA_HPP

#include <fem/str.hpp>
#include <cstdio>
#include <cstdlib>
#include <ctime>

namespace fem {

  inline
  void
  getenv(
    str_cref key,
    str_ref result)
  {
    std::string k = utils::strip_leading_and_trailing_blank_padding(key);
    char* v = std::getenv(k.c_str());
    result = (v == 0 ? " " : v);
  }

  inline
  void
  date(
    str_ref result)
  {
    static const char* months[] = {
      "Jan", "Feb", "Mar", "Apr", "May", "Jun",
      "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
    std::time_t now = std::time(0);
    std::tm const* tm(std::localtime(&now));
    str<10> buf;
    std::sprintf(buf.elems, "%02d-%s-%02d",
      tm->tm_mday,
      months[tm->tm_mon],
      tm->tm_year % 100);
    result = buf;
  }

  inline
  void
  time(
    str_ref result)
  {
    std::time_t now = std::time(0);
    std::tm const* tm(std::localtime(&now));
    str<10> buf;
    std::sprintf(buf.elems, "%02d:%02d:%02d",
      tm->tm_hour,
      tm->tm_min,
      tm->tm_sec);
    result = buf;
  }

  inline
  double
  user_plus_system_time()
  {
    static std::clock_t t_start = std::clock();
    return static_cast<double>(std::clock() - t_start)
         / static_cast<double>(CLOCKS_PER_SEC);
  }

  inline
  void
  cpu_time(
    float& result)
  {
    result = static_cast<float>(user_plus_system_time());
  }

  inline
  void
  cpu_time(
    double& result)
  {
    result = user_plus_system_time();
  }

  inline
  int
  system(
    str_cref command)
  {
    return std::system(std::string(command).c_str());
  }

} // namespace fem

#endif // GUARD
