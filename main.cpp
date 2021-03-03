#include <mutex>

#include "aligner.hpp"
#include "biomodern/hs37d5.hpp"
#include "biomodern/sam.hpp"

using namespace std::string_literals;
using namespace biomodern;
using namespace biomodern::format;
using namespace biomodern::utility;
using namespace std::chrono;

auto total = 0u;
auto mutex = std::mutex{};
auto mutex2 = std::mutex{};

auto encode_fa(std::string fa_path) {
  auto ref = ""_is;
  ref.reserve(3137454505);
  auto fin = std::ifstream{fa_path};
  assert(fin);

  std::cout << "read ref...\n";
  auto start = high_resolution_clock::now();
  for (auto line = ""s; std::getline(fin, line);) {
    if (line.front() == '>') {
      std::cout << line << "\n";
      continue;
    }
    ref += Codec::to_istring(line);
  }
  assert(ref.size() == 3137454505);
  auto end = high_resolution_clock::now();
  std::cout << "elapsed time: " << duration_cast<seconds>(end - start).count() << " s.\n";

  const auto ref_path = fa_path.substr(0, fa_path.size() - 2) + "iref";
  std::cout << "encoded ref path: " << ref_path << '\n';

  auto fout = std::ofstream{ref_path, std::ios::binary};
  assert(fout);
  std::cout << "save encoded ref...\n";
  start = high_resolution_clock::now();
  Serializer::save(fout, ref);
  end = high_resolution_clock::now();
  std::cout << "elapsed time: " << duration_cast<seconds>(end - start).count() << " s.\n";
  return ref;
}

auto build_fmi(istring_view ref, std::string fa_path) {
  const auto fmi = biomodern::FMIndex{ref};
  const auto fmi_path = fa_path.substr(0, fa_path.size() - 2) + "fmi";
  std::cout << "fm index path: " << fmi_path << '\n';
  auto fout = std::ofstream{fmi_path, std::ios::binary};
  assert(fout);
  fmi.save(fout);
}

auto index(std::string fa_path) {
  assert(fa_path.ends_with("hs37d5.fa"));
  std::cout << "fa path: " << fa_path << "\n";
  auto ref = encode_fa(fa_path);

  // fm-index only support four characters so we need change 'N' to 'A'
  std::ranges::replace(ref, 4, 0);
  build_fmi(ref, fa_path);
}

struct Read {
  std::string name;
  std::string seq1, qual1;
  std::string seq2, qual2;
};

const auto chrs =
    std::vector<std::string>{"1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10", "11", "12",
                             "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X",  "Y"};

void do_work(
    std::ifstream& fq1, std::ifstream& fq2, const Aligner& aligner, const int i,
    std::map<std::string, std::ofstream>& fouts, std::string RG_string) {
  constexpr auto batch = 3000;
  for (auto j = std::uint8_t{};;) {
    auto reads = std::vector<Read>{};
    auto chr_records = std::map<std::string, std::vector<SamRecord>>{};
    for (const auto chr : chrs) chr_records[chr];

    auto count = 0u;
    auto name = ""s, seq = ""s, split = ""s, qual = ""s;
    {
      const auto lock = std::lock_guard{mutex};
      while (std::getline(fq1, name)) {
        std::getline(fq1, seq);
        std::getline(fq1, split);
        std::getline(fq1, qual);
        reads.emplace_back();
        reads.back().name = std::move(name);
        reads.back().seq1 = std::move(seq);
        reads.back().qual1 = std::move(qual);

        if (++count == batch) break;
      }
      count = 0u;
      while (std::getline(fq2, name)) {
        std::getline(fq2, seq);
        std::getline(fq2, split);
        std::getline(fq2, qual);
        assert(reads[count].name == name);
        reads[count].seq2 = std::move(seq);
        reads[count].qual2 = std::move(qual);

        if (++count == batch) break;
      }
      total += reads.size();
    }

    if (reads.empty()) return;
    for (const auto& [name, read1, qual1, read2, qual2] : reads) {
      const auto [aln1, aln2] = aligner.map(read1, read2);
      auto
          [gpos1, score1, score21, forward1, read_end1, ref_end1, find_cnt1, align_len1, mapq1,
           sub_score1, rescued1, cigar1] = aln1;
      auto
          [gpos2, score2, score22, forward2, read_end2, ref_end2, find_cnt2, align_len2, mapq2,
           sub_score2, rescued2, cigar2] = aln2;
      if (cigar1.empty()) cigar1 = "*";
      if (cigar2.empty()) cigar2 = "*";
      // flag1
      auto flag1 = SamRecord::READ_PAIRED + SamRecord::FIRST_OF_PAIR,
           flag2 = SamRecord::READ_PAIRED + SamRecord::SECOND_OF_PAIR;
      if (!forward1) {
        flag1 += SamRecord::READ_REVERSE_STRAND;
        flag2 += SamRecord::MATE_REVERSE_STRAND;
      }
      if (!forward2) {
        flag1 += SamRecord::MATE_REVERSE_STRAND;
        flag2 += SamRecord::READ_REVERSE_STRAND;
      }
      auto [chr1, pos1] = Hs37d5::get_chr_pos(gpos1);
      auto [chr2, pos2] = Hs37d5::get_chr_pos(gpos2);
      auto rname1 = std::string{chr1};
      auto rname2 = std::string{chr2};
      if (score1 == 0) {
        flag1 += SamRecord::READ_UNMAPPED;
        flag2 += SamRecord::MATE_UNMAPPED;
        rname1 = "*";
      }
      if (score2 == 0) {
        flag1 += SamRecord::MATE_UNMAPPED;
        flag2 += SamRecord::READ_UNMAPPED;
        rname2 = "*";
      }
      auto rnext1 = rname2, rnext2 = rname1;
      auto pnext1 = pos2, pnext2 = pos1;
      auto tlen1 = 0, tlen2 = 0;
      auto proper_pair = false;
      if (score1 != 0 && score2 != 0 && rname1 == rname2) {
        rnext1 = "=";
        rnext2 = "=";
        tlen1 = SamUtil::compute_tlen(pos1, cigar1, forward1, pos2, cigar2, forward2);
        tlen2 = -tlen1;
        if (forward1 != forward2 && std::abs(tlen1) <= 1300) {
          flag1 += SamRecord::PROPER_PAIR;
          flag2 += SamRecord::PROPER_PAIR;
          proper_pair = true;
        }
      }

      const auto qname = name.substr(0, name.find_first_of(" \t"));
      auto record1 = SamRecord{
          // remove @
          qname.substr(1),
          static_cast<std::uint16_t>(flag1),
          rname1,
          pos1 + 1,
          mapq1,
          cigar1,
          rnext1,
          pos2 + 1,
          tlen1,
          forward1 ? std::move(read1) : biomodern::utility::Codec::rev_comp(read1),
          forward1 ? qual1 : std::string{qual1.rbegin(), qual1.rend()},
          std::string{"AS:i:"} + std::to_string(score1) + "\tXS:i:" + std::to_string(sub_score1) +
              RG_string + (rescued1 ? "\trs:i:1" : "")};

      if (chr_records.contains(rname1)) {
        chr_records[rname1].push_back(std::move(record1));
      } else {
        // chr_records["decoy"].push_back(std::move(record1));
      }

      auto record2 = SamRecord{
          qname.substr(1),
          static_cast<std::uint16_t>(flag2),
          rname2,
          pos2 + 1,
          mapq2,
          cigar2,
          rnext2,
          pos1 + 1,
          tlen2,
          forward2 ? std::move(read2) : biomodern::utility::Codec::rev_comp(read2),
          forward2 ? qual2 : std::string{qual2.rbegin(), qual2.rend()},
          std::string{"AS:i:"} + std::to_string(score2) + "\tXS:i:" + std::to_string(sub_score2) +
              RG_string + (rescued2 ? "\trs:i:1" : "")};

      if (chr_records.contains(rname2)) {
        chr_records[rname2].push_back(std::move(record2));
      } else {
        // chr_records["decoy"].push_back(std::move(record2));
      }
    }

    {
      const auto lock = std::lock_guard{mutex2};
      for (const auto& [chr, records] : chr_records)
        for (const auto& record : records) fouts[chr] << record << "\n";
    }

    if (i == 0 && (j++ % 16) == 0)
      std::cout << "current read name: " << name << " (" << total << ")\n";
  }
}

auto align(
    std::string fa_path, std::string fq1_path, std::string fq2_path, std::string sam_prefix,
    std::string SM, std::string RGID, int thread_num) {
  std::cout << "fa path: " << fa_path << "\n";
  std::cout << "fq1 path: " << fq1_path << "\n";
  std::cout << "fq2 path: " << fq2_path << "\n";

  std::cout << "assume read length: ~148bp\n";
  std::cout << "assume read insert size mean: ~550bp\n";
  std::cout << "assume read insert size var: ~150bp\n";

  assert(fa_path.ends_with("hs37d5.fa"));
  const auto ref_path = fa_path.substr(0, fa_path.size() - 2) + "iref";
  auto ref = ""_is;
  {
    auto fin = std::ifstream{ref_path, std::ios::binary};
    assert(fin);
    std::cout << "load ref...\n";
    Serializer::load(fin, ref);
  }
  std::cout << "ref size: " << ref.size() << '\n';

  // add padding base to prevent out of bounds in a rare case
  for (auto i = 0; i < 1300; i++) ref.push_back(0);
  auto aligner = Aligner{};
  aligner.print();
  aligner.ref_ = ref;

  const auto fmi_path = fa_path.substr(0, fa_path.size() - 2) + "fmi";
  {
    auto fin = std::ifstream{fmi_path, std::ios::binary};
    assert(fin);
    aligner.fm_idx_.load(fin);
  }

  auto fouts = std::map<std::string, std::ofstream>{};
  for (const auto& chr : chrs) {
    fouts.emplace(chr, sam_prefix + "." + chr + ".sam");
    fouts[chr] << "@HD\tVN:1.0\tSO:unsorted\n";
    for (const auto [chr_name, begin, size] : Hs37d5::chr_begin_sizes)
      fouts[chr] << "@SQ\tSN:" << chr_name << "\tLN:" << size << "\n";
    fouts[chr] << "@RG\tID:" << RGID << "\tSM:" << SM << "\n";
  }

  auto fq1 = std::ifstream{fq1_path};
  auto fq2 = std::ifstream{fq2_path};
  assert(fq1);
  assert(fq2);

  std::cout << "using " << thread_num << " threads.\n";
  std::cout << "start mapping...\n";
  auto workers = std::vector<std::thread>{};
  const auto start = high_resolution_clock::now();
  for (auto i = 0; i < thread_num; i++)
    workers.emplace_back(
        do_work, std::ref(fq1), std::ref(fq2), std::cref(aligner), i, std::ref(fouts),
        "\tRG:Z:" + RGID);
  for (auto& worker : workers) worker.join();
  const auto end = high_resolution_clock::now();
  const auto dur = duration_cast<seconds>(end - start);
  std::cout << "elapsed time: " << dur.count() << " s.\n";
  std::cout << total << '\n';
}

auto main(int argc, char* argv[]) -> int {
  if (argc == 1) {
    std::cout << R"(
Contact: Hewill Kang <hewillk@gmail.com>
Command:
         ./hewill index <hs37d5.fa>
         ./hewill align <hs37d5.fa> <in1.fq> <in2.fq> <sam_prefix> <sample_name> <read_group_id> <thread_num>

)";
    return 0;
  }
  const auto cmd = argv[1];
  if (cmd == "index"s) index(argv[2]);  // fa
  else if (cmd == "align"s)
    // fa fq1 fq2 sam_prefix SM RGID thread_num
    align(argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], std::stoi(argv[8]));
}