#pragma once

#include "biomodern/cigar.hpp"
#include "biomodern/fm_index.hpp"
#include "biomodern/istring.hpp"
#include "biomodern/mapq.hpp"
#include "biomodern/smith_waterman.hpp"

namespace biomodern::inline algorithm {

using utility::istring_view;

struct AlignerBase {
  const int INSERT_MEAN = 550;
  const int INSERT_VAR = 150;

  constexpr static auto MAX_HIT_CNT = 512;
  constexpr static auto MAX_EM_CNT = 128;
  constexpr static auto MAX_SW_CNT = 32;
  constexpr static auto MAX_RESCUE_CNT = 128;
  constexpr static auto MAX_SEED_CNT = 4;

  constexpr static auto SEED_LEN = 19;
  constexpr static auto SEED_OVERLAP = 4;
  constexpr static auto EXTEND = 100;
  constexpr static auto SW_THRESHOLD = 30;
  constexpr static auto KMER_SIZE = 8;
  constexpr static auto MIN_FIND_CNT = 4;
  constexpr static auto MAX_FIND_CNT_DIFF = 4;
  constexpr static auto MAX_SW_DIFF = 30;
  constexpr static auto PEN_UNPAIRED = 19;

  auto print() {
    std::cout << "================== constexpr argument ==================\n";
    std::cout << "MAX_HIT_CNT: " << MAX_HIT_CNT << "\n";
    std::cout << "MAX_EM_CNT: " << MAX_EM_CNT << "\n";
    std::cout << "MAX_SW_CNT: " << MAX_SW_CNT << "\n";
    std::cout << "MAX_RESCUE_CNT: " << MAX_RESCUE_CNT << "\n";
    std::cout << "MAX_SEED_CNT: " << MAX_SEED_CNT << "\n";

    std::cout << "SEED_LEN: " << SEED_LEN << "\n";
    std::cout << "SEED_OVERLAP: " << SEED_OVERLAP << "\n";
    std::cout << "EXTEND: " << EXTEND << "\n";
    std::cout << "SW_THRESHOLD: " << SW_THRESHOLD << "\n";
    std::cout << "KMER_SIZE: " << KMER_SIZE << "\n";
    std::cout << "MIN_FIND_CNT: " << MIN_FIND_CNT << "\n";
    std::cout << "MAX_FIND_CNT_DIFF: " << MAX_FIND_CNT_DIFF << "\n";
    std::cout << "MAX_SW_DIFF: " << MAX_SW_DIFF << "\n";
    std::cout << "PEN_UNPAIRED: " << PEN_UNPAIRED << "\n";

    std::cout << "================== const argument ==================\n";
    std::cout << "INSERT_MEAN: " << INSERT_MEAN << "\n";
    std::cout << "INSERT_VAR: " << INSERT_VAR << "\n";
    std::cout << "PAIR_DIST: " << INSERT_MEAN + 4 * INSERT_VAR + 50 << "\n";
  }

  constexpr static auto DIFF = [](const auto a, const auto b) {
    return (a > b) ? (a - b) : (b - a);
  };
  istring_view ref_;
  FMIndex fm_idx_;

  struct Aln {
    std::uint32_t pos{};
    std::uint8_t score{};
    std::uint8_t score2{};
    bool forward{1};
    std::uint8_t read_end{};
    std::uint32_t ref_end{};
    std::uint8_t find_cnt{};
    std::uint8_t align_len{};
    std::uint8_t mapq{};
    std::uint8_t sub_score{};
    bool rescued{};
    std::string cigar;

    friend auto& operator<<(std::ostream& os, Aln aln) {
      return os << "(pos: " << aln.pos << (aln.forward ? "(->)" : "(<-)")
                << ", score: " << +aln.score << ", score2: "
                << +aln.score2
                // << ", align_len: " << +aln.align_len << ", find_cnt: " << +aln.find_cnt
                << ", rescued:" << aln.rescued << ", mapq:" << +aln.mapq << ")";
    }

    bool operator==(const Aln& other) const { return pos == other.pos && forward == other.forward; }
    auto operator<=>(const Aln& other) const {
      if (auto cmp = pos <=> other.pos; cmp != 0) return cmp;
      if (auto cmp = forward <=> other.forward; cmp != 0) return cmp;
      if (auto cmp = other.score <=> score; cmp != 0) return cmp;
      return other.cigar.size() <=> cigar.size();
    }
  };

  struct Anchor {
    std::uint32_t ref_pos{};
    std::uint8_t seed_pos{};
    std::uint8_t seed_size{};
    bool forward{};
    bool repeat{};
    auto operator<=>(const Anchor& other) const noexcept = default;
    friend auto& operator<<(std::ostream& os, const Anchor& an) {
      return os << "(" << +an.seed_pos << ": " << an.ref_pos << ", " << +an.seed_size << ")";
    }
  };

  struct AlnPair {
    Aln aln1;
    Aln aln2;
    auto dist() const noexcept { return DIFF(aln1.pos, aln2.pos); }
    auto score() const noexcept { return aln1.score + aln2.score; }
  };

  struct SeedSpan {
    istring_view seed;
    std::span<const std::uint32_t> span;

    auto operator<=>(const SeedSpan& other) const noexcept {
      if (auto cmp = span.size() <=> other.span.size(); cmp != nullptr) return cmp;
      return other.seed.size() <=> seed.size();
    }
    bool operator==(const SeedSpan& other) const noexcept {
      return seed.size() == other.seed.size() && span.size() == other.span.size();
    }
  };

  static auto filter_alns(std::vector<Aln>& alns) {
    const auto best_score = alns.front().score;
    const auto last = std::ranges::find_if(
        alns, [best_score](const auto score) { return score < best_score - MAX_SW_DIFF; },
        &Aln::score);
    const auto final_size = std::ranges::distance(alns.begin(), last);
    alns.resize(final_size);

    for (auto i = 0; const auto& aln : alns | std::views::drop(final_size)) {
#ifdef DEBUG
      std::cout << "\n************ filtered results ************\n";
      std::cout << "[" << i++ << "] " << aln << '\n';
#endif
    }
  }

  static auto finalize_alns(std::vector<Aln>& alns) {
    if (alns.size() <= 1) return;
    std::ranges::sort(alns);
    const auto [begin, end] = std::ranges::unique(alns);
    alns.erase(begin, end);
    std::ranges::sort(alns, std::ranges::greater{}, &Aln::score);
    filter_alns(alns);
  }

  static auto print_alns(const std::vector<Aln>& alns) {
    for (auto i = 0; const auto& aln : alns | std::views::take(32))
      std::cout << "[" << i++ << "] " << aln << '\n';
  }

  static auto split_read(istring_view read) {
    constexpr auto npos = istring_view::npos;
    auto frags = std::vector<istring_view>{};
    for (auto start = npos, i = 0ul; i <= read.size(); i++) {
      if (i == read.size() || read[i] == 4) {
        if (start != npos && i - start >= SEED_LEN) frags.push_back(read.substr(start, i - start));
        start = npos;
      } else if (start == npos)
        start = i;
    }
    return frags;
  }

  /*
  static auto need_revert(istring_view read, istring_view ref) {
    const auto max_mis_cnt = (read.size() + 4) / 5;
    auto mis_cnt = 0;
    for (auto i = 0; mis_cnt <= max_mis_cnt && i < read.size(); i++)
      if (read[i] != ref[i] && read[i] != 4 && ref[i] != 4) mis_cnt++;
    return mis_cnt <= max_mis_cnt;
  }

  auto try_revert_front(
      Aln& aln, istring_view read, std::uint8_t read_beg, std::uint32_t ref_beg) const {
#ifdef DEBUG
    std::cout << "============== revert front back ==============\n";
#endif

    if (auto cigar = static_cast<format::Cigar>(aln.cigar); cigar.front().op == 'S') {
      const auto sub_read = read.substr(0, read_beg);
      const auto sub_ref = ref_.substr(ref_beg - sub_read.size(), sub_read.size());
      if (need_revert(sub_read, sub_ref)) {
        cigar.front().op = 'M';
        aln.pos -= cigar.front().size;
        cigar.compact();
        aln.cigar = cigar;
      }
    }
  }

  auto try_revert_back(
      Aln& aln, istring_view read, std::uint8_t read_end, std::uint32_t ref_end) const {
#ifdef DEBUG
    std::cout << "============== revert cigar back ==============\n";
#endif

    if (auto cigar = static_cast<format::Cigar>(aln.cigar); cigar.back().op == 'S') {
      const auto sub_read = read.substr(read_end + 1);
      const auto sub_ref = ref_.substr(ref_end + 1, sub_read.size());
      if (need_revert(sub_read, sub_ref)) {
        cigar.back().op = 'M';
        cigar.compact();
        aln.cigar = cigar;
      }
    }
  }
  */

  auto set_cigar(Aln& aln, istring_view read, auto& profile) const {
#ifdef DEBUG
    std::cout << "============== compute cigar ==============\n";
#endif

    if (!aln.cigar.empty()) return;
    if (aln.score == read.size()) {
      aln.cigar = std::to_string(read.size()) + 'M';
#ifdef DEBUG
      std::cout << "full score: " << aln.cigar << "\n";
#endif
      aln.align_len = read.size();
      return;
    }

    const auto sw_pos = aln.ref_end - read.size() - EXTEND;
    const auto subref = ref_.substr(sw_pos, read.size() + EXTEND + 1);

    if (profile.read == nullptr) profile = SWAligner::get_profile(read);
    auto sw = SWAligner::align(profile, subref, true, true, SW_THRESHOLD);

    const auto ref_beg = sw_pos + sw.ref_beg;
    const auto ref_end = sw_pos + sw.ref_end;
    aln.pos = ref_beg;
    aln.score = sw.score;
    aln.cigar = std::move(sw.cigar);
    aln.align_len = std::max(sw.ref_end - sw.ref_beg + 1, sw.read_end - sw.read_beg + 1);

#ifdef DEBUG
    std::cout << "pos: " << aln.pos << "\n";
    std::cout << "raw cigar: " << aln.cigar << "\n";
#endif

    /*
    if (aln.forward) try_revert_back(aln, read, sw.read_end, ref_end);
    else
      try_revert_front(aln, read, sw.read_beg, ref_beg);
    */
  }

  auto display_chains(const auto& chains) const {
    std::cout << "************ seed chain ************\n";
    for (const auto& anchors : chains | std::views::take(8)) {
      const auto& front = anchors.front();
      std::cout << "(" << anchors.size() << ", " << (front.forward ? "->): " : "<-): ");
      std::cout << front;
      for (const auto anchor : anchors | std::views::drop(1)) std::cout << "->" << anchor;
      std::cout << "\n";
    }
  }

  static auto get_score(istring_view read, istring_view ref, bool forward)
      -> std::pair<std::uint16_t, format::Cigar> {
    const auto full_score = read.size();

    const auto read_beg = read.substr(0, 5);
    const auto ref_beg = ref.substr(0, 5);
    const auto read_end = read.substr(read.size() - 5);
    const auto ref_end = ref.substr(ref.size() - 5);

    read.remove_prefix(5);
    ref.remove_prefix(5);
    read.remove_suffix(5);
    ref.remove_suffix(5);

    constexpr auto only_one_mismatch = [](istring_view read, istring_view ref) {
      const auto [in1, in2] = std::ranges::mismatch(read, ref);
      return read.substr(in1 - read.begin() + 1) == ref.substr(in2 - ref.begin() + 1);
    };

    if (read != ref) {
      if (read_beg != ref_beg || read_end != ref_end) return {};
      if (!only_one_mismatch(read, ref)) return {};

      return {full_score - 5, std::to_string(full_score) + 'M'};
    }

    if (read_beg == ref_beg && read_end == ref_end) {
      return {full_score, std::to_string(full_score) + 'M'};
    }

    constexpr auto get_beg_score = [](auto read_beg, auto ref_beg,
                                      bool forward) -> std::pair<std::uint16_t, std::string> {
      constexpr static auto cigars = std::array{"5S", "4S1M", "3S2M", "2S3M", "1S4M"};
      for (auto clip_idx = 4; clip_idx >= 0; clip_idx--) {
        if (read_beg[clip_idx] != ref_beg[clip_idx]) {
          auto cigar = cigars[4 - clip_idx];
          // if (!forward && read_beg.substr(0, clip_idx) == ref_beg.substr(0, clip_idx)) cigar =
          // "5M";
          return {4 - clip_idx, cigar};
        }
      }
      return {5, "5M"};
    };

    constexpr auto get_end_score = [](auto read_end, auto ref_end,
                                      bool forward) -> std::pair<std::uint16_t, std::string> {
      constexpr static auto cigars = std::array{"5S", "1M4S", "2M3S", "3M2S", "4M1S"};
      for (auto clip_idx = 0; clip_idx <= 4; clip_idx++) {
        if (read_end[clip_idx] != ref_end[clip_idx]) {
          auto cigar = cigars[clip_idx];
          /*
          if (forward && read_end.substr(clip_idx + 1) == ref_end.substr(clip_idx + 1))
            cigar = "5M";
          */
          return {clip_idx, cigar};
        }
      }
      return {5, "5M"};
    };

    auto score = read.size();
    auto mid_cigar = std::to_string(score) + 'M';
    const auto [end_score, end_cigar] = get_end_score(read_end, ref_end, forward);
    const auto [beg_score, beg_cigar] = get_beg_score(read_beg, ref_beg, forward);
    score += (beg_score + end_score);
    if (score < full_score - 5) return {};
    auto cigar = static_cast<format::Cigar>(beg_cigar + mid_cigar + end_cigar);
    cigar.compact();
    return {score, std::move(cigar)};
  }

  static auto compute_se_mapq(const std::vector<Aln>& alns, float frac_rep) {
    const auto [opt_score, sub_score, sub_cnt] =
        utility::get_opt_subopt_count(alns | std::views::transform(&Aln::score));
    const auto& front = alns.front();
    const auto mapq = utility::mem_approx_mapq_se(
        {opt_score, front.score2, sub_score, front.align_len, sub_cnt, frac_rep});
    return std::pair{mapq, sub_score};
  }

  auto get_best_one(
      const std::vector<Aln>& alns, istring_view read, istring_view rread, auto& profile,
      auto& rprofile, float frac_rep) const {
    auto aln = alns.front();
    if (aln.forward) set_cigar(aln, read, profile);
    else
      set_cigar(aln, rread, rprofile);
    const auto [mapq, sub_score] = compute_se_mapq(alns, frac_rep);
    aln.mapq = mapq;
    aln.sub_score = sub_score;
    return aln;
  }

  // O(NlogN) pairing
  static auto pairing2(std::vector<Aln> alns1, std::vector<Aln> alns2, int PAIR_DIST) {
    std::ranges::sort(alns1);
    std::ranges::sort(alns2);
    auto aln_pairs = std::vector<AlnPair>{};
    auto begin = alns2.begin();
    auto end = begin;
    for (const auto& aln1 : alns1) {
      while (begin != alns2.end() && begin->pos < aln1.pos - PAIR_DIST) ++begin;
      while (end != alns2.end() && end->pos < aln1.pos + PAIR_DIST) ++end;
      for (auto it = begin; it != end; ++it) {
        const auto& aln2 = *it;
        if (aln1.forward != aln2.forward) aln_pairs.emplace_back(aln1, aln2);
      }
    }
    std::ranges::sort(aln_pairs, std::ranges::greater{}, &AlnPair::score);
    return aln_pairs;
  }

  static auto print_paires(const std::vector<AlnPair>& aln_pairs) {
    for (auto i = 0; const auto& aln_pair : aln_pairs) {
      std::cout << "[" << i++ << "] "
                << "{ (" << aln_pair.aln1 << ") <-- " << aln_pair.dist() << " --> ("
                << aln_pair.aln2 << ") } (score: " << aln_pair.score() << ")\n";
    }
  }
};

}  // namespace biomodern::inline algorithm
