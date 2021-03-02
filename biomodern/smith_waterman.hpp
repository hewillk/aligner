#pragma once

#include <array>
#include <string>

#include "ssw.h"

namespace biomodern {

struct SWAligner {
  constexpr static auto w_match = 1;
  constexpr static auto w_mismatch = 4;
  constexpr static auto w_open = 6;
  constexpr static auto w_extend = 1;
  constexpr static auto w_ambig = 1;

  constexpr static auto mat = [] {
    auto mat = std::array<std::int8_t, 25>{};
    auto k = 0;
    for (auto i = 0; i < 4; i++) {
      for (auto j = 0; j < 4; j++) mat[k++] = i == j ? w_match : -w_mismatch;
      mat[k++] = -w_ambig;
    }
    for (auto i = 0; i < 5; i++) mat[k++] = -w_ambig;
    return mat;
  }();

  struct SWResult {
    int score{};
    int score2{};
    int ref_beg{};
    int ref_end{};
    int read_beg{};
    int read_end{};
    int ref_end2{};
    std::string cigar;
  };

 private:
  static auto extract_result(s_align* res, int read_size) {
    const auto
        [score1, score2, ref_beg1, ref_end1, read_beg1, read_end1, ref_end2, icigar, cigar_size] =
            *res;
    auto result = SWResult{score1, score2, ref_beg1, ref_end1, read_beg1, read_end1, ref_end2, {}};
    if (auto& cigar = result.cigar; cigar_size > 0) {
      if (read_beg1 > 0) cigar += std::to_string(read_beg1) + 'S';

      for (auto i = 0; i < cigar_size; i++)
        cigar += std::to_string(cigar_int_to_len(icigar[i])) + cigar_int_to_op(icigar[i]);

      if (const auto end = read_size - read_end1 - 1; end > 0) cigar += std::to_string(end) + 'S';
    }
    align_destroy(res);
    return result;
  }

 public:
  template <std::ranges::contiguous_range R>
  requires std::same_as<std::ranges::range_value_t<R>, std::int8_t> static auto get_profile(
      const R& read) {
    return ssw_init(read.data(), read.size(), mat.data(), 5, 0);
  }

  template <std::ranges::contiguous_range R>
  requires std::same_as<std::ranges::range_value_t<R>, std::int8_t> static auto align(
      const s_profile& profile, const R& ref, bool report_beg = true, bool report_cigar = true,
      int min_score = 0) {
    auto flag = 0;
    if (report_beg) flag |= 0x08;
    if (report_cigar) flag |= 0x0F;
    auto res = ssw_align(
        &profile, ref.data(), ref.size(), w_open, w_extend, flag, min_score, 32767,
        profile.readLen / 2);
    return extract_result(res, profile.readLen);
  }
};

}  // namespace biomodern