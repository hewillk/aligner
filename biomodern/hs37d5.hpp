#pragma once

#include <algorithm>
#include <string_view>
#include <tuple>
#include <vector>

namespace biomodern::utility {

using namespace std::string_view_literals;

struct Hs37d5 {
  constexpr static auto whole_size = 3137454505u;

  static auto get_chr_pos(std::uint32_t pos) {
    auto upper = std::ranges::upper_bound(
        chr_begin_sizes, pos, {}, [](const auto t) { return std::get<1>(t); });
    const auto [chr, begin, size] = *(--upper);
    return std::pair{chr, pos - begin};
  }

  static inline const auto chr_begin_sizes =
      std::vector{std::tuple{"1"sv, 0u, 249250621u},  {"2", 249250621, 243199373},
                  {"3", 492449994, 198022430},        {"4", 690472424, 191154276},
                  {"5", 881626700, 180915260},        {"6", 1062541960, 171115067},
                  {"7", 1233657027, 159138663},       {"8", 1392795690, 146364022},
                  {"9", 1539159712, 141213431},       {"10", 1680373143, 135534747},
                  {"11", 1815907890, 135006516},      {"12", 1950914406, 133851895},
                  {"13", 2084766301, 115169878},      {"14", 2199936179, 107349540},
                  {"15", 2307285719, 102531392},      {"16", 2409817111, 90354753},
                  {"17", 2500171864, 81195210},       {"18", 2581367074, 78077248},
                  {"19", 2659444322, 59128983},       {"20", 2718573305, 63025520},
                  {"21", 2781598825, 48129895},       {"22", 2829728720, 51304566},
                  {"X", 2881033286, 155270560},       {"Y", 3036303846, 59373566},
                  {"MT", 3095677412, 16569},          {"GL000207.1", 3095693981, 4262},
                  {"GL000226.1", 3095698243, 15008},  {"GL000229.1", 3095713251, 19913},
                  {"GL000231.1", 3095733164, 27386},  {"GL000210.1", 3095760550, 27682},
                  {"GL000239.1", 3095788232, 33824},  {"GL000235.1", 3095822056, 34474},
                  {"GL000201.1", 3095856530, 36148},  {"GL000247.1", 3095892678, 36422},
                  {"GL000245.1", 3095929100, 36651},  {"GL000197.1", 3095965751, 37175},
                  {"GL000203.1", 3096002926, 37498},  {"GL000246.1", 3096040424, 38154},
                  {"GL000249.1", 3096078578, 38502},  {"GL000196.1", 3096117080, 38914},
                  {"GL000248.1", 3096155994, 39786},  {"GL000244.1", 3096195780, 39929},
                  {"GL000238.1", 3096235709, 39939},  {"GL000202.1", 3096275648, 40103},
                  {"GL000234.1", 3096315751, 40531},  {"GL000232.1", 3096356282, 40652},
                  {"GL000206.1", 3096396934, 41001},  {"GL000240.1", 3096437935, 41933},
                  {"GL000236.1", 3096479868, 41934},  {"GL000241.1", 3096521802, 42152},
                  {"GL000243.1", 3096563954, 43341},  {"GL000242.1", 3096607295, 43523},
                  {"GL000230.1", 3096650818, 43691},  {"GL000237.1", 3096694509, 45867},
                  {"GL000233.1", 3096740376, 45941},  {"GL000204.1", 3096786317, 81310},
                  {"GL000198.1", 3096867627, 90085},  {"GL000208.1", 3096957712, 92689},
                  {"GL000191.1", 3097050401, 106433}, {"GL000227.1", 3097156834, 128374},
                  {"GL000228.1", 3097285208, 129120}, {"GL000214.1", 3097414328, 137718},
                  {"GL000221.1", 3097552046, 155397}, {"GL000209.1", 3097707443, 159169},
                  {"GL000218.1", 3097866612, 161147}, {"GL000220.1", 3098027759, 161802},
                  {"GL000213.1", 3098189561, 164239}, {"GL000211.1", 3098353800, 166566},
                  {"GL000199.1", 3098520366, 169874}, {"GL000217.1", 3098690240, 172149},
                  {"GL000216.1", 3098862389, 172294}, {"GL000215.1", 3099034683, 172545},
                  {"GL000205.1", 3099207228, 174588}, {"GL000219.1", 3099381816, 179198},
                  {"GL000224.1", 3099561014, 179693}, {"GL000223.1", 3099740707, 180455},
                  {"GL000195.1", 3099921162, 182896}, {"GL000212.1", 3100104058, 186858},
                  {"GL000222.1", 3100290916, 186861}, {"GL000200.1", 3100477777, 187035},
                  {"GL000193.1", 3100664812, 189789}, {"GL000194.1", 3100854601, 191469},
                  {"GL000225.1", 3101046070, 211173}, {"GL000192.1", 3101257243, 547496},
                  {"NC_007605", 3101804739, 171823},  {"hs37d5", 3101976562, 35477943}};
};

}  // namespace biomodern::utility