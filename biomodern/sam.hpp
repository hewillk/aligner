#pragma once

#include <istream>
#include <ostream>
#include <string>

namespace biomodern::format {

struct SamRecord {
  std::string qname;
  std::uint16_t flag{};
  std::string rname;
  std::uint32_t pos{};
  std::uint16_t mapq{};
  std::string cigar;
  std::string rnext;
  std::uint32_t pnext{};
  std::int32_t tlen{};
  std::string seq;
  std::string qual;
  std::string optionals;

  enum flag {
    READ_PAIRED = 0x1,
    PROPER_PAIR = 0x2,
    READ_UNMAPPED = 0x4,
    MATE_UNMAPPED = 0x8,
    READ_REVERSE_STRAND = 0x10,
    MATE_REVERSE_STRAND = 0x20,
    FIRST_OF_PAIR = 0x40,
    SECOND_OF_PAIR = 0x80,
    SECONDARY_ALIGNMENT = 0x100,
    READ_FAILS_QUALITY_CHECK = 0x200,
    DUPLICATE_READ = 0x400,
    SUPPLEMENTARY_ALIGNMENT = 0x800
  };

  auto read_paired() const noexcept { return !!(flag & READ_PAIRED); }
  auto proper_pair() const noexcept { return !!(flag & PROPER_PAIR); }
  auto read_unmapped() const noexcept { return !!(flag & READ_UNMAPPED); }
  auto mate_unmapped() const noexcept { return !!(flag & MATE_UNMAPPED); }
  auto read_reverse_strand() const noexcept { return !!(flag & READ_REVERSE_STRAND); }
  auto mate_reverse_strand() const noexcept { return !!(flag & MATE_REVERSE_STRAND); }
  auto first_of_pair() const noexcept { return !!(flag & FIRST_OF_PAIR); }
  auto second_of_pair() const noexcept { return !!(flag & SECOND_OF_PAIR); }
  auto secondary_alignment() const noexcept { return !!(flag & SECONDARY_ALIGNMENT); }
  auto read_fails_quality_check() const noexcept { return !!(flag & READ_FAILS_QUALITY_CHECK); }
  auto duplicate_read() const noexcept { return !!(flag & DUPLICATE_READ); }
  auto supplementary_alignment() const noexcept { return !!(flag & SUPPLEMENTARY_ALIGNMENT); }
};

auto& operator<<(std::ostream& os, const SamRecord& r) {
  os << r.qname << '\t' << r.flag << '\t' << r.rname << '\t' << r.pos << '\t' << r.mapq << '\t'
     << r.cigar << '\t' << r.rnext << '\t' << r.pnext << '\t' << r.tlen << '\t' << r.seq << '\t'
     << r.qual << '\t' << r.optionals;
  return os;
}

struct SamUtil {
  enum Orientation { FR, FF, RR, RF };
  static auto compute_ori(bool read_forward, bool mate_forward) noexcept {
    if (read_forward != mate_forward) return read_forward ? FR : RF;
    return read_forward ? FF : RR;
  }
  static auto compute_tlen(
      std::int32_t read_pos, const Cigar& read_cigar, bool read_forward, std::int32_t mate_pos,
      const Cigar& mate_cigar, bool mate_forward) noexcept -> std::int32_t {
    if (read_pos > mate_pos)
      return -compute_tlen(mate_pos, mate_cigar, mate_forward, read_pos, read_cigar, read_forward);

    switch (const auto read_ori = compute_ori(read_forward, mate_forward); read_ori) {
      case FR: return mate_pos + mate_cigar.ref_size() - read_pos;
      case FF:
        if (const auto tlen =
                mate_pos + mate_cigar.read_size() - (read_pos + read_cigar.read_size());
            tlen != 0)
          return tlen + (tlen > 0 ? 1 : -1);
        else
          return 0;
      case RR:
        if (const auto tlen = mate_pos + mate_cigar.ref_size() - (read_pos + read_cigar.ref_size());
            tlen != 0)
          return tlen + (tlen > 0 ? 1 : -1);
        else
          return 0;
      case RF:
        if (const auto tlen = mate_pos - (read_pos + read_cigar.ref_size()) + 1; tlen != 0)
          return tlen + (tlen > 0 ? 1 : -1);
        else
          return 0;
      default: return 0;
    }
  }
};

}  // namespace biomodern::format