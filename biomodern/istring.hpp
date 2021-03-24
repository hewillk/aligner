#pragma once

#include <algorithm>
#include <array>
#include <istream>
#include <ostream>
#include <string>

namespace biomodern::utility {

using namespace std::string_literals;
using istring = std::basic_string<std::int8_t>;
using istring_view = std::basic_string_view<std::int8_t>;

inline auto operator""_s(const char* s) {
  auto is = istring{};
  for (const auto c : std::string_view{s}) is.push_back(c - '0');
  return is;
}

struct Codec {
  constexpr static auto ints = [] {
    auto ints = std::array<std::int8_t, 128>{};
    ints.fill(4);
    ints['a'] = ints['A'] = 0;
    ints['c'] = ints['C'] = 1;
    ints['g'] = ints['G'] = 2;
    ints['t'] = ints['T'] = 3;
    return ints;
  }();
  constexpr static auto to_int(char c) noexcept { return ints[c]; }
  constexpr static auto is_valid(char c) noexcept { return to_int(c) != 4; }
  constexpr static auto chars = std::array{'A', 'C', 'G', 'T', 'N'};
  constexpr static auto to_char(int i) noexcept { return chars[i]; }

  constexpr static auto hash(istring_view seq) noexcept {
    auto key = 0ull;
    for (auto i = 0; i < seq.size(); i++) key |= (seq[seq.size() - i - 1] & 3ull) << 2 * i;
    return key;
  };

  static auto rhash(std::size_t key, std::size_t size) {
    auto seq = istring{};
    for (auto i = 0; i < size; i++) {
      const auto shift = (size - i - 1) * 2;
      seq += (key & 3ull << shift) >> shift;
    }
    return seq;
  }

  static auto rev_comp(istring_view seq) {
    auto res = istring{};
    res.reserve(seq.size());
    std::transform(seq.rbegin(), seq.rend(), std::back_inserter(res), [](const auto c) {
      return c == 4 ? 4 : 3 - c;
    });
    return res;
  }

  static auto to_istring(std::string_view seq) {
    auto res = istring{};
    res.reserve(seq.size());
    std::transform(seq.begin(), seq.end(), std::back_inserter(res), to_int);
    return res;
  }

  constexpr static auto comp = [](char c) {
    switch (c) {
      case 'a':
      case 'A': return 'T';
      case 'c':
      case 'C': return 'G';
      case 'g':
      case 'G': return 'C';
      case 't':
      case 'T': return 'A';
      default: return 'N';
    }
  };

  static auto rev_comp(std::string_view seq) {
    auto res = std::string{};
    res.reserve(seq.size());
    std::transform(seq.rbegin(), seq.rend(), std::back_inserter(res), comp);
    return res;
  }
};

}  // namespace biomodern::utility

inline auto& operator<<(std::ostream& out, biomodern::utility::istring_view s) {
  for (const auto c : s) out << +c;
  return out;
}

inline auto& operator>>(std::istream& in, biomodern::utility::istring& is) {
  auto s = std::string{};
  in >> s;
  is.clear();
  for (const auto c : s) is.push_back(c - '0');
  return in;
}
