//
//  base64 encoding and decoding with C++.
//  Version: 1.01.00
//  from: https://github.com/ReneNyffenegger/cpp-base64 see base64.cc for licensing details

#ifndef BASE64_H_C0CE2A47_D10E_42C9_A27C_C883944E704A
#define BASE64_H_C0CE2A47_D10E_42C9_A27C_C883944E704A

#include <string>

namespace utility {
namespace io {

std::string base64_encode(unsigned char const* , unsigned int len);
std::string base64_decode(std::string const& s);

}  // namespace io
}  // namespace utility

#endif /* BASE64_H_C0CE2A47_D10E_42C9_A27C_C883944E704A */
