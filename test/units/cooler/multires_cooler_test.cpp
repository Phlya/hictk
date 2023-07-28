// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler/multires_cooler.hpp"

#include <fmt/format.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "tmpdir.hpp"

namespace hictk::cooler::test::multires_cooler_file {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MultiResCooler: open read-only", "[cooler][short]") {
  const auto path = datadir / "multires_cooler_test_file.mcool";

  auto mclr = MultiResFile(path.string());

  CHECK(mclr.resolutions().size() == 7);
  CHECK(mclr.attributes().format == MCOOL_MAGIC);
  CHECK(mclr.attributes().format_version == 2);
  CHECK(!mclr.attributes().bin_type);

  CHECK(utils::is_cooler(mclr.open(1600000).uri()));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MultiResCooler: init files", "[cooler][short]") {
  const auto base_path = datadir / "cooler_test_file.cool";
  const auto base_resolution = File(base_path.string()).bin_size();

  const auto path = testdir() / "test_init.mcool";
  const std::array<std::uint32_t, 4> resolutions{
      // clang-format off
      base_resolution,
      base_resolution * 2,
      base_resolution * 4,
      base_resolution * 8
      // clang-format on
  };

  SECTION("coarsen on construction") {
    SECTION("valid resolutions") {
      std::ignore = MultiResFile::create(path.string(), File(base_path.string()),
                                         resolutions.begin(), resolutions.end(), true);

      CHECK(utils::is_multires_file(path.string()));
    }

    SECTION("invalid resolutions") {
      std::vector<std::uint32_t> resolutions_{base_resolution / 2};
      CHECK_THROWS(MultiResFile::create(path.string(), File(base_path.string()),
                                        resolutions_.begin(), resolutions_.end(), true));
      resolutions_ = {base_resolution + 1};
      CHECK_THROWS(MultiResFile::create(path.string(), File(base_path.string()),
                                        resolutions_.begin(), resolutions_.end(), true));
    }
  }

  SECTION("construct then initialize") {
    const Reference chroms{Chromosome{0, "chr1", 10000}, Chromosome{1, "chr2", 5000}};
    auto mclr = MultiResFile::create(path.string(), chroms, true);

    for (const auto res : resolutions) {
      std::ignore = File::create(fmt::format(FMT_STRING("{}::/resolutions/{}"), path.string(), res),
                                 chroms, res);
    }

    CHECK(utils::is_multires_file(path.string()));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MultiResCooler: create resolutions", "[cooler][short]") {
  const auto base_path = datadir / "cooler_test_file.cool";
  File base_clr(base_path.string());
  const auto base_resolution = base_clr.bin_size();

  const auto path = testdir() / "test_create_resolutions.mcool";
  const std::array<std::uint32_t, 3> resolutions{
      // clang-format off
      base_resolution * 2,
      base_resolution * 4,
      base_resolution * 8
      // clang-format on
  };

  auto mclr = MultiResFile::create(path.string(), base_clr.chromosomes(), true);
  mclr.copy_resolution(base_clr);

  SECTION("valid resolutions") {
    for (const auto res : resolutions) {
      std::ignore = mclr.create_resolution(res);
    }

    CHECK(mclr.resolutions().size() == resolutions.size() + 1);
  }
  SECTION("invalid resolutions") {
    CHECK_THROWS(mclr.create_resolution(base_resolution / 2));
    CHECK_THROWS(mclr.create_resolution(base_resolution + 1));
  }
}

}  // namespace hictk::cooler::test::multires_cooler_file
