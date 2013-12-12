#ifndef FEM_UTILS_PATH_HPP
#define FEM_UTILS_PATH_HPP

#include <fem/utils/char.hpp>

#include <string>
#include <cstdio>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#if defined(_MSC_VER)

#include <direct.h> // for _getcwd
#include <io.h> // for _chsize

#define FEM_UTILS_PATH_STRUCT_STAT struct _stat
#define FEM_UTILS_PATH_STAT _stat
#define FEM_UTILS_PATH_GETCWD _getcwd
#define FEM_UTILS_PATH_FTRUNCATE _chsize

#else

#include <unistd.h> // for ftruncate

#define FEM_UTILS_PATH_STRUCT_STAT struct stat
#define FEM_UTILS_PATH_STAT stat
#define FEM_UTILS_PATH_GETCWD getcwd
#define FEM_UTILS_PATH_FTRUNCATE ftruncate

#endif

namespace fem { namespace utils { namespace path {

  inline
  bool
  exists(
    char const* path)
  {
    FEM_UTILS_PATH_STRUCT_STAT buf;
    int stat_result = FEM_UTILS_PATH_STAT(path, &buf);
    return (stat_result == 0 || errno != ENOENT);
  }

  inline
  char const*
  split_drive(
    char const* path)
  {
#if defined(_MSC_VER)
    if (is_a_through_z(path[0]) && path[1] == ':') return path+2;
#endif
    return path;
  }

  inline
  bool
  is_absolute(
    char const* path,
    bool drive_split_already=false)
  {
    if (!drive_split_already) path = split_drive(path);
#if defined(_MSC_VER)
    return (path[0] == '/' || path[0] == '\\'); // emulates Python
#else
    return (path[0] == '/');
#endif
  }

  inline
  const char*
  separator()
  {
#if defined(_MSC_VER)
    return "\\";
#else
    return "/";
#endif
  }

  inline
  std::string
  absolute(
    char const* path)
  {
    char const* path_without_drive = split_drive(path);
    if (is_absolute(path_without_drive, /*drive_split_already*/ true)) {
      return std::string(path);
    }
    std::string result;
    static const size_t buf_size = 10000; // ad-hoc
    char buf[buf_size];
    char* getcwd_result = FEM_UTILS_PATH_GETCWD(buf, buf_size);
    if (getcwd_result == 0) {
      int en = errno;
      std::string msg = "fem::utils::path::absolute(): ";
      if (en != 0) {
        msg += std::strerror(en);
      }
      else {
        msg += "unknown error";
      }
      throw std::runtime_error(msg);
    }
    else {
      if (path != path_without_drive) {
        result += std::string(path, path_without_drive);
      }
      result += getcwd_result;
      result += separator();
      result += path_without_drive;
    }
    return result;
  }

  inline
  bool
  truncate_file_at_current_position(
    std::FILE* fp)
  {
    long curr_pos = std::ftell(fp);
    if (curr_pos < 0) return false;
    if (FEM_UTILS_PATH_FTRUNCATE(fileno(fp), curr_pos) != 0) return false;
    std::fflush(fp);
    if (std::fseek(fp, 0L, SEEK_END) != 0) return false;
    return true;
  }

}}} // namespace fem::utils::path

#endif // GUARD
