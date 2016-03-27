
#list(REMOVE_ITEM CMAKE_CXX_FLAGS "-fno-exceptions")
#set(CMAKE_CXX_FLAGS "-fexceptions ${CMAKE_CXX_FLAGS}")

set(LLVM_LINK_COMPONENTS support)
set(LLVM_USED_LIBS clangTooling clangBasic clangAST)

add_clang_executable(binder
  binder.cpp

  context.hpp
  context.cpp

  class.hpp
  class.cpp

  config.hpp
  config.cpp

  enum.hpp
  enum.cpp

  function.hpp
  function.cpp

  type.hpp
  type.cpp

  util.hpp
  util.cpp

  fmt/format.cc
  fmt/posix.cc
  )

target_link_libraries(binder
  clangTooling
  clangBasic
  clangASTMatchers
  )
