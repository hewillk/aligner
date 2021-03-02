#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <ranges>

namespace biomodern::utility {

constexpr auto get_opt_subopt_count(auto scores, int diff = 5) noexcept
    -> std::tuple<int, int, int> {
  if (scores.empty()) return {};
  if (scores.size() == 1) return {scores.front(), 0, 0};
  const auto opt_score = scores.front();
  const auto sub_score = scores[1];
  const auto min_score = sub_score - diff;
  const auto remains = scores | std::views::drop(1);
  const auto it =
      std::ranges::find_if(remains, [min_score](const auto score) { return score < min_score; });
  const auto sub_cnt = std::ranges::distance(remains.begin(), it);
  return {opt_score, sub_score, sub_cnt};
}

struct MemAln {
  int score{};      // local best SW score
  int score2{};     // local 2nd SW score
  int sub_score{};  // global 2nd SW score
  int align_len{};  // alignment length
  int sub_n{};      // number of global 2nd SW score hits
  float frac_rep{};
};

static auto mem_approx_mapq_se(MemAln aln) {
  const auto [score, csub, sub_score, l, sub_n, frac_rep] = aln;
  auto sub = sub_score ? sub_score : 20;
  sub = csub > sub ? csub : sub;
  if (sub >= score) return 0;
  const auto identity = 1. - (l - score) / 5. / l;
  auto tmp = l < 50 ? 1. : 3.912 / std::log(l);
  tmp *= identity * identity;
  auto mapq = static_cast<int>(6.02 * (score - sub) * tmp * tmp + .499);
  if (sub_n > 0) mapq -= static_cast<int>(4.343 * std::log(sub_n + 1) + .499);
  if (mapq > 60) mapq = 60;
  if (mapq < 0) mapq = 0;
  mapq = static_cast<int>(mapq * (1. - frac_rep) + .499);
  return mapq;
}

static auto mem_mapq_pe(MemAln p0, MemAln p1, int score_un, int o, int subo, int n_sub) {
  constexpr auto raw_mapq = [](int diff) { return static_cast<int>(6.02 * diff + .499); };

  subo = subo > score_un ? subo : score_un;
  auto q_pe = raw_mapq(o - subo);
  if (n_sub > 0) q_pe -= static_cast<int>(4.343 * std::log(n_sub + 1) + .499);
  if (q_pe < 0) q_pe = 0;
  if (q_pe > 60) q_pe = 60;
  q_pe = static_cast<int>(q_pe * (1. - .5 * (p0.frac_rep + p1.frac_rep)) + .499);
  auto q_se0 = mem_approx_mapq_se(p0);
  auto q_se1 = mem_approx_mapq_se(p1);
  q_se0 = q_se0 > q_pe ? q_se0 : q_pe < q_se0 + 40 ? q_pe : q_se0 + 40;
  q_se1 = q_se1 > q_pe ? q_se1 : q_pe < q_se1 + 40 ? q_pe : q_se1 + 40;
  q_se0 = q_se0 < raw_mapq(p0.score - p0.score2) ? q_se0 : raw_mapq(p0.score - p0.score2);
  q_se1 = q_se1 < raw_mapq(p1.score - p1.score2) ? q_se1 : raw_mapq(p1.score - p1.score2);
  if (q_se0 > 60) q_se0 = 60;
  if (q_se1 > 60) q_se1 = 60;
  return std::pair{q_se0, q_se1};
}

}  // namespace biomodern::utility