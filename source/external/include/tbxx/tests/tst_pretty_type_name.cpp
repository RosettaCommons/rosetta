#include <iostream>
#include <tbxx/pretty_type_name.hpp>
#include <tbxx/error_utils.hpp>

template <class T, class S>
class foo
{};

int main() {
  TBXX_ASSERT(tbxx::pretty_type_name<int>() == "int");
  TBXX_ASSERT((tbxx::pretty_type_name< foo<int, double> >() == "foo<int, double>"));
  std::cout << "OK\n";
  return 0;
}
