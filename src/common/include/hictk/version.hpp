// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// This file was generated automatically by CMake.

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <string_view>

#include "hictk/git.hpp"

namespace hictk::config::version {

// clang-format off
[[nodiscard]] constexpr std::uint_fast8_t major() noexcept { return 0; }
[[nodiscard]] constexpr std::uint_fast8_t minor() noexcept { return 0; }
[[nodiscard]] constexpr std::uint_fast8_t patch() noexcept { return 1; }
// clang-format on

namespace internal {

[[nodiscard]] inline std::string generate_suffix() {
  namespace git = hictk::config::git;
  constexpr std::string_view raw_suffix{""};
  if (!git::state_available()) {
    return std::string{raw_suffix};
  }

  const std::string is_dirty_suffix = git::is_dirty() ? "-dirty" : "";
  const auto short_hash_suffix = std::string{"-"} + std::string{git::describe()};
  const auto short_version = fmt::format(FMT_STRING("{}.{}.{}"), major(), minor(), patch());

  std::string buff{};
  const auto expected_release_tag =
      "v" + short_version + (raw_suffix.empty() ? "" : "-" + std::string{raw_suffix});
  if (git::tag() == expected_release_tag) {
    buff = std::string{raw_suffix};
  } else {
    buff = std::string{raw_suffix} + short_hash_suffix;
  }
  buff.append(is_dirty_suffix);

  if (buff.front() == '-') {
    buff.erase(0, 1);
  }

  return buff;
}

}  // namespace internal

[[nodiscard]] inline std::string_view suffix() noexcept {
  static const std::string buff{internal::generate_suffix()};
  return buff;
}

[[nodiscard]] inline std::string_view str() noexcept {
  static const std::string buff =
      suffix().empty()
          ? fmt::format(FMT_STRING("{}.{}.{}"), major(), minor(), patch())
          : fmt::format(FMT_STRING("{}.{}.{}-{}"), major(), minor(), patch(), suffix());
  return buff;
}

[[nodiscard]] inline std::string_view str_long(std::string_view prefix = "hictk") noexcept {
  assert(!prefix.empty());
  static const std::string buff{fmt::format(FMT_STRING("{}-v{}"), prefix, str())};
  return buff;
}

}  // namespace hictk::config::version
