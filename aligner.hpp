#pragma once

#include <map>

#include "aligner_base.hpp"

namespace biomodern::inline algorithm {

struct Aligner : AlignerBase {
  auto get_equal_span(
      std::span<const std::uint32_t> span, std::uint32_t offset, istring_view equal_seed) const {
    const auto [begin, end] = std::ranges::equal_range(
        span, equal_seed, {}, [this, offset, equal_size = equal_seed.size()](const auto pos) {
          return ref_.substr(pos + offset, equal_size);
        });
    return span.subspan(begin - span.begin(), end - begin);
  }

  auto get_spans(istring_view read) const {
    auto seed_spans = std::vector<SeedSpan>{};
    auto repeat_size = 0;
    auto origin_read = read;
    while (read.size() >= SEED_LEN) {
      const auto seed = read.substr(read.size() - SEED_LEN);
      const auto [begin, end, offset] = fm_idx_.get_range(seed, 0);
      const auto span = fm_idx_.get_offsets(begin, end);
      if (span.size() <= MAX_HIT_CNT) {
#ifdef DEBUG
        const auto cur_seed = seed.substr(offset);
        std::cout << "seed: " << cur_seed << "(" << cur_seed.size() << ") -> (" << span.size()
                  << ")\n";
#endif
        if (!span.empty()) seed_spans.emplace_back(seed.substr(offset), span);
        read.remove_suffix(SEED_LEN - offset - SEED_OVERLAP);
      } else {
        assert(offset == 0);
        const auto seed2 = read.substr(0, read.size() - SEED_LEN);
        const auto [begin2, end2, offset2] = fm_idx_.get_range(seed2, begin, end, MAX_HIT_CNT);
        const auto span2 = fm_idx_.get_offsets(begin2, end2);
        if (span2.size() <= MAX_HIT_CNT) {
#ifdef DEBUG
          const auto cur_seed = read.substr(offset2);
          std::cout << "seed: " << cur_seed << "(" << cur_seed.size() << ") -> (" << span2.size()
                    << ")\n";
#endif
          if (!span2.empty()) seed_spans.emplace_back(read.substr(offset2), span2);
          read = read.substr(0, offset2 + SEED_OVERLAP);
        } else {
          assert(offset2 == 0);
#ifdef DEBUG
          std::cout << "seed: " << read << "(" << read.size() << ") -> (" << span2.size() << ")\n";
          std::cout << "\n************ backtrace extend ************\n";
#endif

          auto remain_seed = origin_read.substr(read.size());
          auto equal_offset = read.size();
          auto equal_span = span2;
          while (!remain_seed.empty()) {
            const auto equal_size = std::min(
                static_cast<int>(remain_seed.size()),
                static_cast<int>(std::log2(equal_span.size()) / 2));
            const auto equal_seed = remain_seed.substr(0, equal_size);
            equal_span = get_equal_span(equal_span, equal_offset, equal_seed);

#ifdef DEBUG
            std::cout << "equal_seed: " << equal_seed << "(" << equal_seed.size() << ") -> ("
                      << equal_span.size() << ")\n";
#endif
            remain_seed.remove_prefix(equal_size);
            equal_offset += equal_size;
            if (equal_span.size() <= MAX_HIT_CNT) {
              const auto extend_seed = origin_read.substr(0, equal_offset);
              seed_spans.emplace_back(extend_seed, equal_span);

#ifdef DEBUG
              std::cout << "seed: " << extend_seed << "(" << extend_seed.size() << ") -> ("
                        << equal_span.size() << ")\n";

#endif
              break;
            }
          }

          repeat_size = read.size();
          read = {};
          break;
        }
      }
    }
    // recover
    if (read.size() >= SEED_LEN - SEED_OVERLAP) {
      const auto [begin, end, offset] = fm_idx_.get_range(read, 0);
      const auto span = fm_idx_.get_offsets(begin, end);
      if (span.size() <= MAX_HIT_CNT) {
#ifdef DEBUG
        std::cout << "seed: " << read.substr(offset) << " -> (" << span.size() << ")\n";
#endif
        if (!span.empty()) seed_spans.emplace_back(read.substr(offset), span);
      }
    }
    return std::pair{std::move(seed_spans), repeat_size};
  }

  auto seeding_impl(istring_view read, bool forward) const {
    const auto frags = split_read(read);
#ifdef DEBUG
    std::cout << "\n************ read seeds ************\n";
#endif

    auto seed_spans = std::vector<SeedSpan>{};
    auto repeats = 0;
    for (auto i = 0; const auto frag : frags) {
      const auto [spans, repeat_size] = get_spans(frag);
      repeats += repeat_size;
      for (const auto span : spans) seed_spans.push_back(span);
    }
    std::ranges::sort(seed_spans);

    auto anchors = std::vector<Anchor>{};
    for (const auto [seed, span] : seed_spans | std::views::take(MAX_SEED_CNT)) {
      const auto seed_pos = seed.data() - read.data();
      for (const auto ref_pos : span) anchors.emplace_back(ref_pos, seed_pos, seed.size(), forward);
    }

    auto chain_map = std::map<std::uint32_t, std::vector<Anchor>>{};
    for (const auto anchor : anchors) {
      const auto read_pos = anchor.ref_pos - anchor.seed_pos;
      const auto lower = chain_map.lower_bound(read_pos - SEED_LEN);
      const auto upper = chain_map.upper_bound(read_pos + SEED_LEN);
      if (lower != upper) {
        for (auto it = lower; it != upper; ++it) it->second.push_back(anchor);
        continue;
      }
      chain_map[read_pos].push_back(anchor);
    }

    auto chains = std::vector<std::vector<Anchor>>{};
    for (auto& chain : chain_map | std::views::values) {
      std::ranges::sort(chain);
      chains.push_back(std::move(chain));
    }
    return std::pair{std::move(chains), repeats};
  }

  auto seeding(istring_view read, istring_view rread) const {
    auto [chains, repeats] = seeding_impl(read, true);

#ifdef DEBUG
    std::cout << "\n************ reverse ************\n";
    std::cout << "rread: ";
    for (auto i = 0; i < rread.size(); i++) std::cout << (rread[i] != 'N' ? rread[i] : '|');
    std::cout << "\n";
#endif

    const auto [rchains, rrepeats] = seeding_impl(rread, false);
    std::ranges::copy(rchains, std::back_inserter(chains));
    std::ranges::sort(chains, std::ranges::greater{}, &std::vector<Anchor>::size);
    return std::tuple{std::move(chains), repeats, rrepeats};
  }

  auto exact_match(auto& chains, istring_view read, istring_view rread, int find_cnt) const {
    auto alns = std::vector<Aln>{};
    auto sw_chains = std::vector<std::vector<Anchor>>{};
    for (auto& chain : chains) {
      const auto [ref_pos, seed_pos, seed_size, forward, repeat] = chain.front();
      auto read_pos = ref_pos - seed_pos;
      const auto sw_read = forward ? read : rread;
      if (const auto [score, cigar] =
              get_score(sw_read, ref_.substr(read_pos, read.size()), forward);
          score) {
        if (const auto [size, op] = cigar.front(); op == 'S') { read_pos += size; }
        alns.emplace_back(read_pos, score, 0, forward).cigar = cigar;
        alns.back().find_cnt = find_cnt;
        alns.back().align_len = cigar.ref_size();
      } else {
        sw_chains.push_back(std::move(chain));
      }
    }
    if (sw_chains.size() > MAX_EM_CNT) sw_chains.resize(MAX_EM_CNT);
    return std::pair{std::move(alns), std::move(sw_chains)};
  }

  static auto get_kmers(istring_view read) {
    auto kmers = std::vector<std::uint32_t>{};
    // overlap one base
    for (auto i = 0; i < read.size() / (KMER_SIZE - 1); i++) {
      auto needle = read.substr(i * (KMER_SIZE - 1), KMER_SIZE);
      kmers.push_back(Codec::hash(needle));
    }
    return kmers;
  }

  static auto find_kmers(
      const std::vector<std::uint32_t>& kmers, istring_view ref, std::vector<bool>& table) {
    table.clear();
    table.resize(1 << KMER_SIZE * 2);

    for (auto i = 0; i < ref.size() - KMER_SIZE + 1; i++)
      table[Codec::hash(ref.substr(i, KMER_SIZE))] = true;

    auto find_cnt = 0;
    for (const auto kmer : kmers) find_cnt += table[kmer];
    return find_cnt;
  }

  auto get_sw_alns(
      const auto& chains, int read_size, const std::vector<std::uint32_t>& kmers,
      const std::vector<std::uint32_t>& rkmers, std::vector<bool>& table, int min_find_cnt) const {
    auto alns = std::vector<Aln>{};
    for (const auto& chain : chains) {
      const auto [ref_pos, seed_pos, seed_size, forward, repeat] = chain.front();
      const auto read_pos = ref_pos - seed_pos;
      const auto front_pad = seed_pos <= EXTEND / 2 ? seed_pos * 2 : EXTEND;
      const auto sw_pos = read_pos - front_pad;
      const auto subref = ref_.substr(sw_pos, read_size + 2 * EXTEND);

      const auto find_cnt = find_kmers(forward ? kmers : rkmers, subref, table);
      if (find_cnt < min_find_cnt) continue;
      min_find_cnt = std::max(find_cnt - MAX_FIND_CNT_DIFF, min_find_cnt);

      alns.emplace_back(sw_pos, 0, 0, forward, 0, 0, find_cnt);
    }
    std::ranges::sort(alns, std::ranges::greater{}, &Aln::find_cnt);
    return std::pair{std::move(alns), min_find_cnt};
  }

  auto get_sw_candidates(
      bool alns_empty, const auto& chains, int read_size, const std::vector<std::uint32_t>& kmers,
      const std::vector<std::uint32_t>& rkmers, std::vector<bool>& table) const {
    if (alns_empty) {
      return get_sw_alns(chains, read_size, kmers, rkmers, table, MIN_FIND_CNT);
    } else {
      const auto indel_chains = chains | std::views::take_while([](const auto& chain) {
                                  return chain.size() >= MAX_SEED_CNT / 2;
                                });
      return get_sw_alns(
          indel_chains, read_size, kmers, rkmers, table, kmers.size() - MAX_FIND_CNT_DIFF);
    }
  }

  auto extending(
      std::vector<Aln>& alns, const std::vector<Aln>& sw_alns, istring_view read,
      istring_view rread) const {
    auto profile = s_profile{}, rprofile = s_profile{};
    for (auto min_score = SW_THRESHOLD; const auto& sw_aln : sw_alns) {
      const auto subref = ref_.substr(sw_aln.pos, read.size() + 2 * EXTEND);
      auto sw = SWAligner::SWResult{};
      if (sw_aln.forward) {
        if (profile.read == nullptr) profile = SWAligner::get_profile(read);
        sw = SWAligner::align(profile, subref, false, false, min_score);
      } else {
        if (rprofile.read == nullptr) rprofile = SWAligner::get_profile(rread);
        sw = SWAligner::align(rprofile, subref, false, false, min_score);
      }

      const auto score = static_cast<int>(sw.score);
      if (score < min_score) continue;
      const auto ref_end = sw_aln.pos + sw.ref_end;
      alns.emplace_back(
          ref_end - sw.read_end, score, sw.score2, sw_aln.forward, sw.read_end, ref_end,
          sw_aln.find_cnt);

      min_score = std::max(min_score, score - MAX_SW_DIFF);
    }
    return std::pair{std::move(profile), std::move(rprofile)};
  }

  static auto release_memory(auto& v) {
    v.clear();
    v.shrink_to_fit();
  }

  auto rescue(
      const std::vector<Aln>& alns1, const std::vector<Aln>& alns2, istring_view read2,
      istring_view rread2, auto& profile2, auto& rprofile2,
      const std::vector<std::uint32_t>& kmers2, const std::vector<std::uint32_t>& rkmers2,
      std::vector<bool>& table, int min_find_cnt, int PAIR_DIST) const {
    const auto [opt_score, sub_score, sub_cnt] =
        utility::get_opt_subopt_count(alns1 | std::views::transform(&Aln::score));
    const auto rescue_cnt = std::min(sub_cnt + 1, MAX_RESCUE_CNT);

#ifdef DEBUG
    std::cout << "============== rescue count: " << rescue_cnt << " ==============\n";
#endif

    auto rescues = std::vector<Aln>{};
    for (auto min_score = SW_THRESHOLD; const auto& aln1 : alns1 | std::views::take(rescue_cnt)) {
      const auto pos1 = aln1.pos;
      if (std::ranges::any_of(
              alns2, [pos1, PAIR_DIST](const auto pos2) { return DIFF(pos1, pos2) <= PAIR_DIST; },
              &Aln::pos)) {
#ifdef DEBUG
        std::cout << "pos: " << pos1 << " already seen.\n";
#endif
        continue;
      }

      const auto forward1 = aln1.forward;
      const auto sw_pos = forward1 ? pos1 - EXTEND : pos1 - PAIR_DIST;
      const auto subref = ref_.substr(sw_pos, EXTEND + read2.size() + PAIR_DIST);

      const auto find_cnt = find_kmers(forward1 ? rkmers2 : kmers2, subref, table);
      if (find_cnt < min_find_cnt) continue;
      min_find_cnt = std::max(find_cnt - MAX_FIND_CNT_DIFF, min_find_cnt);

      auto sw = SWAligner::SWResult{};
      if (forward1) {
        if (rprofile2.read == nullptr) rprofile2 = SWAligner::get_profile(rread2);
        sw = SWAligner::align(rprofile2, subref, false, false, min_score);
      } else {
        if (profile2.read == nullptr) profile2 = SWAligner::get_profile(read2);
        sw = SWAligner::align(profile2, subref, false, false, min_score);
      }

      const auto score = static_cast<int>(sw.score);
      if (score < min_score) continue;
      const auto ref_end = sw_pos + sw.ref_end;
      const auto ref_pos = ref_end - sw.read_end;
#ifdef DEBUG
      std::cout << "{ ref pos: " << ref_pos << "(" << static_cast<std::int32_t>(ref_pos - sw_pos)
                << "), score: " << score << " }\n";
#endif

      rescues.emplace_back(ref_pos, score, sw.score2, !forward1, sw.read_end, ref_end, find_cnt)
          .rescued = true;
      min_score = std::max(min_score, score - MAX_SW_DIFF);
    }

    return rescues;
  }

  static auto shrink_sw_size(int em_size1, auto& sw_alns1, int em_size2, auto& sw_alns2) {
    if (sw_alns1.size() > MAX_SW_CNT) sw_alns1.resize(MAX_SW_CNT);
    if (em_size1 > MAX_EM_CNT) sw_alns1.clear();
    if (sw_alns2.size() > MAX_SW_CNT) sw_alns2.resize(MAX_SW_CNT);
    if (em_size2 > MAX_EM_CNT) sw_alns2.clear();

    if (em_size1 == 0 || em_size2 == 0) return;

    const auto sw_size1 = static_cast<int>(sw_alns1.size());
    const auto sw_size2 = static_cast<int>(sw_alns2.size());
    const auto total_size1 = em_size1 + sw_size1;
    const auto total_size2 = em_size2 + sw_size2;

    const auto shrink_size = std::min(total_size1, total_size2);
    if (shrink_size < total_size1) {
      const auto final_sw_size = std::max(0, shrink_size - em_size1);
      if (sw_alns1.size() > final_sw_size) sw_alns1.resize(final_sw_size);
    }
    if (shrink_size < total_size2) {
      const auto final_sw_size = std::max(0, shrink_size - em_size2);
      if (sw_alns2.size() > final_sw_size) sw_alns2.resize(final_sw_size);
    }
  }

  auto insert_penalty(int dist) const {
    const auto ns = (dist - INSERT_MEAN) / static_cast<double>(INSERT_VAR);
    return static_cast<int>(std::pow(ns, 2));
    // return static_cast<int>(.721 * std::log(2. * std::erfc(std::fabs(ns) * M_SQRT1_2)) + .499);
  }

  auto get_best_pair(
      const std::vector<Aln>& alns1, const std::vector<Aln>& alns2,
      const std::vector<AlnPair>& aln_pairs, istring_view read1, istring_view rread1,
      auto& profile1, auto& rprofile1, istring_view read2, istring_view rread2, auto& profile2,
      auto& rprofile2, float frac_rep1, float frac_rep2) const {
#ifdef DEBUG
    std::cout << "************ pairing results ************\n";
    print_paires(aln_pairs);
#endif

    const auto [opt_score1, sub_score1, sub_cnt1] =
        utility::get_opt_subopt_count(alns1 | std::views::transform(&Aln::score));
    const auto [opt_score2, sub_score2, sub_cnt2] =
        utility::get_opt_subopt_count(alns2 | std::views::transform(&Aln::score));
    const auto [opt_score, sub_score, sub_cnt] =
        utility::get_opt_subopt_count(aln_pairs | std::views::transform(&AlnPair::score));

    auto [aln1, aln2] = aln_pairs.front();
    const auto score_un = opt_score1 + opt_score2 - PEN_UNPAIRED;
    const auto success = (opt_score > score_un);
    if (!success) {
      aln1 = alns1.front();
      aln2 = alns2.front();
    }
    if (aln1.forward) set_cigar(aln1, read1, profile1);
    else
      set_cigar(aln1, rread1, rprofile1);
    if (aln2.forward) set_cigar(aln2, read2, profile2);
    else
      set_cigar(aln2, rread2, rprofile2);
    if (!success) {
      aln1.mapq = utility::mem_approx_mapq_se(
          {aln1.score, aln1.score2, sub_score1, aln1.align_len, sub_cnt1, frac_rep1});
      aln2.mapq = utility::mem_approx_mapq_se(
          {aln2.score, aln2.score2, sub_score2, aln2.align_len, sub_cnt2, frac_rep2});
    } else {
      const auto [mapq1, mapq2] = utility::mem_mapq_pe(
          {aln1.score, aln1.score2, sub_score1, aln1.align_len, sub_cnt1, frac_rep1},
          {aln2.score, aln2.score2, sub_score2, aln2.align_len, sub_cnt2, frac_rep2}, score_un,
          opt_score, sub_score, sub_cnt);
#ifdef DEBUG
      std::cout << "(raw mapq1:" << mapq1 << ", raw mapq2: " << mapq2 << ")\n";
#endif
      const auto pen_paired = insert_penalty(aln_pairs.front().dist());
      aln1.mapq = std::max(mapq1 - pen_paired, 0);
      aln2.mapq = std::max(mapq2 - pen_paired, 0);
    }
    aln1.sub_score = (aln1.score == opt_score1 ? sub_score1 : opt_score1);
    aln2.sub_score = (aln2.score == opt_score2 ? sub_score2 : opt_score2);
    return AlnPair{std::move(aln1), std::move(aln2)};
  }

  auto map(istring_view read1, istring_view rread1, istring_view read2, istring_view rread2) const
      -> AlnPair {
    const auto PAIR_DIST = INSERT_MEAN + 4 * INSERT_VAR + 50;

#ifdef DEBUG
    std::cout << "--------------- seeding read1 ---------------\n";
#endif

    const auto [chains1, repeats1, rrepeats1] = seeding(read1, rread1);

#ifdef DEBUG
    std::cout << "\n--------------- seeding read2 ---------------\n";
#endif

    const auto [chains2, repeats2, rrepeats2] = seeding(read2, rread2);

    const auto frac_rep1 = (repeats1 + rrepeats1) / (read1.size() * 2.f);
    const auto frac_rep2 = (repeats2 + rrepeats2) / (read2.size() * 2.f);

    auto table = std::vector<bool>(1 << KMER_SIZE * 2);
    const auto kmers1 = get_kmers(read1);
    const auto rkmers1 = get_kmers(rread1);
    const auto kmers2 = get_kmers(read2);
    const auto rkmers2 = get_kmers(rread2);

#ifdef DEBUG
    std::cout << "\n--------------- exact match read1 ---------------\n";
#endif

    auto [alns1, sw_chains1] = exact_match(chains1, read1, rread1, kmers1.size());

#ifdef DEBUG
    std::cout << "\n--------------- exact match read2 ---------------\n";
#endif

    auto [alns2, sw_chains2] = exact_match(chains2, read2, rread2, kmers2.size());

#ifdef DEBUG
    std::cout << "\n--------------- sw read1 ---------------\n";
#endif

    auto [sw_alns1, min_find_cnt1] =
        get_sw_candidates(alns1.empty(), sw_chains1, read1.size(), kmers1, rkmers1, table);

#ifdef DEBUG
    std::cout << "\n--------------- sw read2 ---------------\n";
#endif

    auto [sw_alns2, min_find_cnt2] =
        get_sw_candidates(alns2.empty(), sw_chains2, read2.size(), kmers2, rkmers2, table);

    release_memory(sw_chains1);
    release_memory(sw_chains2);

    shrink_sw_size(alns1.size(), sw_alns1, alns2.size(), sw_alns2);

    auto [profile1, rprofile1] = extending(alns1, sw_alns1, read1, rread1);
    auto [profile2, rprofile2] = extending(alns2, sw_alns2, read2, rread2);

    release_memory(sw_alns1);
    release_memory(sw_alns2);

    if (alns1.empty() && alns2.empty()) [[unlikely]] return {};

    finalize_alns(alns1);
    finalize_alns(alns2);

#ifdef DEBUG
    std::cout << "\n************ force rescue read1 ************\n";
#endif
    auto rescues1 = rescue(
        alns2, alns1, read1, rread1, profile1, rprofile1, kmers1, rkmers1, table, min_find_cnt1,
        PAIR_DIST);
#ifdef DEBUG
    std::cout << "\n************ force rescue read2 ************\n";
#endif
    auto rescues2 = rescue(
        alns1, alns2, read2, rread2, profile2, rprofile2, kmers2, rkmers2, table, min_find_cnt2,
        PAIR_DIST);

    if (!rescues1.empty()) {
      std::ranges::copy(rescues1, std::back_inserter(alns1));
      finalize_alns(alns1);
    }
#ifdef DEBUG
    std::cout << "\n************ read1 final result (" << alns1.size() << ") ************\n";
    print_alns(alns1);
#endif

    if (!rescues2.empty()) {
      std::ranges::copy(rescues2, std::back_inserter(alns2));
      finalize_alns(alns2);
    }
#ifdef DEBUG
    std::cout << "\n************ read2 final result (" << alns2.size() << ") ************\n";
    print_alns(alns2);
#endif

    release_memory(table);
    release_memory(rescues1);
    release_memory(rescues2);

    if (alns2.empty())
      return {get_best_one(alns1, read1, rread1, profile1, rprofile1, frac_rep1), {}};
    if (alns1.empty())
      return {{}, get_best_one(alns2, read2, rread2, profile2, rprofile2, frac_rep2)};

#ifdef DEBUG
    std::cout << "\n--------------- pairing ---------------\n";
#endif

    auto aln_pairs = pairing2(alns1, alns2, PAIR_DIST);
    if (aln_pairs.empty()) {
#ifdef DEBUG
      std::cout << "\n--------------- failed ---------------\n";
#endif
      return {
          get_best_one(alns1, read1, rread1, profile1, rprofile1, frac_rep1),
          get_best_one(alns2, read2, rread2, profile2, rprofile2, frac_rep2)};
    }
    return get_best_pair(
        alns1, alns2, aln_pairs, read1, rread1, profile1, rprofile1, read2, rread2, profile2,
        rprofile2, frac_rep1, frac_rep2);
  }

  auto map(std::string_view read1, std::string_view read2) const {
    const auto iread1 = Codec::to_istring(read1);
    const auto riread1 = Codec::rev_comp(iread1);
    const auto iread2 = Codec::to_istring(read2);
    const auto riread2 = Codec::rev_comp(iread2);
    return map(iread1, riread1, iread2, riread2);
  }
};

}  // namespace biomodern::inline algorithm
