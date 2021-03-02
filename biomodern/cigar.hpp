#pragma once

#include <algorithm>
#include <concepts>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

namespace biomodern::format {

struct Cigar {
  struct Element {
    std::uint32_t size{};
    char op{};
    operator std::string() const { return std::to_string(size) + op; }
  };

 private:
  std::vector<Element> elements;
  static auto to_elements(std::string_view cigar_string) {
    auto elements = std::vector<Element>{};
    for (auto i = 0u; i < cigar_string.size(); i++) {
      auto size = cigar_string[i] - '0';
      for (i++; std::isdigit(cigar_string[i]); i++) size = size * 10 + cigar_string[i] - '0';
      elements.emplace_back(size, cigar_string[i]);
    }
    return elements;
  }

 public:
  Cigar() = default;
  Cigar(std::convertible_to<std::string_view> auto const& cigar_string)
      : elements(to_elements(cigar_string)) {}
  auto& operator=(std::convertible_to<std::string_view> auto const& cigar_string) {
    elements = to_elements(cigar_string);
    return *this;
  }

  auto compact() {
    if (elements.size() <= 1) return;
    auto res = std::vector<Element>{};
    auto [last_size, last_op] = elements.front();
    for (const auto [cur_size, cur_op] : elements | std::views::drop(1)) {
      if (cur_op == last_op) last_size += cur_size;
      else {
        res.emplace_back(last_size, last_op);
        last_size = cur_size;
        last_op = cur_op;
      }
    }
    res.emplace_back(last_size, last_op);
    elements.swap(res);
  }

  auto emplace_back(std::uint32_t size, char op) { elements.emplace_back(size, op); }

  auto push_back(Element element) { elements.push_back(element); }

  auto append(const Cigar& other) {
    for (const auto element : other.elements) elements.push_back(element);
  }

  auto swap(Cigar& other) { elements.swap(other.elements); }

  auto ref_size() const noexcept {
    auto ref_size = 0;
    for (auto [size, op] : elements) {
      switch (op) {
        case 'M':
        case 'D':
        case 'N':
        case '=':
        case 'X': ref_size += size; break;
        default: break;
      }
    }
    return ref_size;
  }

  auto read_size() const noexcept {
    auto read_size = 0;
    for (auto [size, op] : elements) {
      switch (op) {
        case 'M':
        case 'I':
        case 'S':
        case '=':
        case 'X': read_size += size; break;
        default: break;
      }
    }
    return read_size;
  }

  auto clip_size() const noexcept {
    auto clip_size = 0;
    for (auto [size, op] : elements) {
      switch (op) {
        case 'S':
        case 'H': clip_size += size; break;
        default: break;
      }
    }
    return clip_size;
  }

  auto begin() { return elements.begin(); }

  auto begin() const { return elements.begin(); }

  auto end() { return elements.end(); }

  auto end() const { return elements.end(); }

  auto& front() { return elements.front(); }

  auto& front() const { return elements.front(); }

  auto& back() { return elements.back(); }

  auto& back() const { return elements.back(); }

  auto& operator[](std::uint32_t i) { return elements[i]; }

  auto& operator[](std::uint32_t i) const { return elements[i]; }

  operator std::string() const {
    auto cigar_string = std::string{};
    for (const auto element : elements) cigar_string += element;
    return cigar_string;
  }

  auto pop_front() { elements.erase(elements.begin()); }

  auto pop_back() { elements.pop_back(); }

  auto reverse() { std::ranges::reverse(elements); }

  auto contains(char key) const noexcept {
    for (auto [size, op] : elements)
      if (key == op) return true;
    return false;
  }

  auto contains(std::string_view keys) const noexcept {
    for (auto [size, op] : elements)
      for (auto key : keys)
        if (key == op) return true;
    return false;
  }

  auto size() const noexcept { return elements.size(); }

  auto operator==(std::string_view cigar_string) const noexcept {
    return static_cast<std::string>(*this) == cigar_string;
  }

  friend auto& operator<<(std::ostream& os, const Cigar& cigar) {
    for (const auto [size, op] : cigar) os << size << op;
    return os;
  }

  friend auto& operator>>(std::istream& is, Cigar& cigar) {
    auto cigar_string = std::string{};
    is >> cigar_string;
    cigar = cigar_string;
    return is;
  }
};

}  // namespace biomodern::format