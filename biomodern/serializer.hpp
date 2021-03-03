#pragma once

#include <cassert>
#include <fstream>
#include <ranges>
#include <vector>

namespace biomodern::utility {

struct Serializer {
 private:
  template <std::ranges::random_access_range R>
  constexpr static auto get_bytes(const R& r) noexcept {
    if constexpr (requires { r.num_blocks(); })
      return r.num_blocks() * sizeof(typename R::block_type);
    else if constexpr (std::same_as<R, std::vector<bool>>)
      return ((r.size() - 1) / 64 + 1) * sizeof(std::uint64_t);
    else
      return r.size() * sizeof(typename R::value_type);
  }

  template <std::ranges::random_access_range R>
  static auto get_data(R& r) noexcept {
    if constexpr (std::same_as<std::remove_const_t<R>, std::vector<bool>>)
      return *reinterpret_cast<std::uint64_t* const*>(&r);
    else
      return r.data();
  }

 public:
  template <std::ranges::random_access_range R, auto BUF_SIZE = 8192>
  requires std::is_trivially_copyable_v<std::ranges::range_value_t<R>> static auto save(
      std::ofstream& fout, const R& r) {
    const auto size = r.size();
    if (size == 0) return;
    // write size of range
    fout.write(reinterpret_cast<const char*>(&size), sizeof(size));
    const auto bytes = get_bytes(r);
    // write contents
    auto begin = reinterpret_cast<const char*>(get_data(r));
    for (auto batch = 0u; batch < bytes / BUF_SIZE; batch++, begin += BUF_SIZE)
      fout.write(begin, BUF_SIZE);
    // write remains
    if (const auto remains = bytes % BUF_SIZE; remains) fout.write(begin, remains);
  }

  template <std::ranges::random_access_range R, auto BUF_SIZE = 8192>
  requires std::is_trivially_copyable_v<std::ranges::range_value_t<R>> static auto load(
      std::ifstream& fin, R& r) {
    const auto size = [&] {
      auto size = r.size();
      // read size of results
      fin.read(reinterpret_cast<char*>(&size), sizeof(size));
      return size;
    }();
    if (size == 0) return;
    r.resize(size);

    const auto bytes = get_bytes(r);
    // read contents
    auto begin = reinterpret_cast<char*>(get_data(r));
    for (auto batch = 0u; batch < bytes / BUF_SIZE; batch++, begin += BUF_SIZE) {
      fin.read(begin, BUF_SIZE);
      assert(fin.gcount() == BUF_SIZE);
    }
    // read remains
    if (const auto remains = bytes % BUF_SIZE; remains) {
      fin.read(begin, remains);
      assert(fin.gcount() == remains);
    }
  }
};

}  // namespace biomodern::utility
