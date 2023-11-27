// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/transformers.hpp"

#include <fmt/format.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "hictk/cooler/cooler.hpp"
#include "hictk/hic.hpp"

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::test::transformers {

using namespace hictk::transformers;

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Transformers (cooler)", "[transformers][short]") {
  SECTION("join genomic coords") {
    const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());
    auto sel = clr.fetch("chr1", 5'000'000, 10'000'000);
    const auto jsel =
        JoinGenomicCoords(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), clr.bins_ptr());
    constexpr std::array<std::uint32_t, 3> expected{5'000'000, 5'000'000, 7'500'000};
    const auto pixels = jsel.read_all();
    REQUIRE(pixels.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(pixels[i].coords.bin1.start() == expected[i]);
    }
  }

  SECTION("coarsen") {
    const auto path1 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/100000";
    const auto path2 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/200000";
    const cooler::File clr1(path1.string());
    const cooler::File clr2(path2.string());

    auto sel = clr1.fetch("1");
    auto sel1 =
        CoarsenPixels(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), clr1.bins_ptr(), 2);
    auto sel2 = clr2.fetch("1");

    const auto v1 = sel1.read_all();
    const auto v2 = sel2.read_all<std::int32_t>();
    REQUIRE(v1.size() == v2.size());

    for (std::size_t i = 0; i < v1.size(); ++i) {
      CHECK(v1[i] == v2[i].to_thin());
    }
  }
  SECTION("coarsen recursive") {
    const auto path1 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/100000";
    const auto path2 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/400000";
    const cooler::File clr1(path1.string());
    const cooler::File clr2(path2.string());

    auto sel = clr1.fetch("1");
    auto sel1 =
        CoarsenPixels(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), clr1.bins_ptr(), 2);
    auto sel2 = CoarsenPixels(sel1.begin(), sel1.end(), sel1.dest_bins_ptr(), 2);
    auto sel3 = clr2.fetch("1");

    const auto v1 = sel2.read_all();
    const auto v2 = sel3.read_all<std::int32_t>();
    REQUIRE(v1.size() == v2.size());

    for (std::size_t i = 0; i < v1.size(); ++i) {
      CHECK(v1[i] == v2[i].to_thin());
    }
  }

  SECTION("coarsen gw") {
    const auto path1 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/100000";
    const auto path2 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/200000";
    const cooler::File clr1(path1.string());
    const cooler::File clr2(path2.string());

    auto sel = clr1.fetch();
    auto sel1 =
        CoarsenPixels(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), clr1.bins_ptr(), 2);
    auto sel2 = clr2.fetch();

    const auto v1 = sel1.read_all();
    const auto v2 = sel2.read_all<std::int32_t>();
    REQUIRE(v1.size() == v2.size());

    for (std::size_t i = 0; i < v1.size(); ++i) {
      CHECK(v1[i] == v2[i].to_thin());
    }
  }

  SECTION("stats") {
    const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());
    auto sel = clr.fetch("chr1");
    auto first = sel.begin<std::int32_t>();
    auto last = sel.end<std::int32_t>();

    CHECK_THAT(avg(first, last), Catch::Matchers::WithinRel(25231.981858902574));
    CHECK(nnz(first, last) == 4'465);
    CHECK(max(first, last) == 1'357'124);
    CHECK(sum(first, last) == 112'660'799);
  }

  SECTION("random sample") {
    const auto path = datadir / "cooler/cooler_test_file.cool";
    const std::uint64_t seed = 1234567890;
    const double frac = 0.75;

    const cooler::File clr(path.string());
    const auto sel = clr.fetch("1");

    std::mt19937_64 rand_eng{seed};
    const auto pixels = sel.read_all<std::int32_t>();
    std::vector<ThinPixel<std::int32_t>> subsampled_pixels_expected{};

    std::for_each(sel.template begin<std::int32_t>(), sel.template end<std::int32_t>(),
                  [&](ThinPixel<std::int32_t> p) {
                    p.count = std::binomial_distribution<std::int32_t>{p.count, frac}(rand_eng);
                    if (p.count != 0) {
                      subsampled_pixels_expected.push_back(p);
                    }
                  });

    const transformers::PixelRandomSampler sampler(sel.template begin<std::int32_t>(),
                                                   sel.template end<std::int32_t>(), frac, seed);

    const auto subsampled_pixels = sampler.read_all();
    REQUIRE(subsampled_pixels.size() == subsampled_pixels_expected.size());
    for (std::size_t i = 0; i < subsampled_pixels.size(); ++i) {
      CHECK(subsampled_pixels_expected[i] == subsampled_pixels[i]);
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Transformers (hic)", "[transformers][short]") {
  auto path = datadir / "hic/4DNFIZ1ZVXC8.hic8";

  SECTION("join genomic coords") {
    const hic::File hf(path.string(), 2'500'000);
    auto sel = hf.fetch("chr2L", 5'000'000, 10'000'000);
    const auto jsel =
        JoinGenomicCoords(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), hf.bins_ptr());
    constexpr std::array<std::uint32_t, 3> expected{5'000'000, 5'000'000, 7'500'000};
    const auto pixels = jsel.read_all();
    REQUIRE(pixels.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(pixels[i].coords.bin1.start() == expected[i]);
    }
  }

  SECTION("coarsen") {
    const hic::File hf1(path.string(), 500'000);
    const hic::File hf2(path.string(), 2'500'000);

    auto sel = hf1.fetch("chr2R");
    auto sel1 =
        CoarsenPixels(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), hf1.bins_ptr(), 5);
    auto sel2 = hf2.fetch("chr2R");

    const auto v1 = sel1.read_all();
    const auto v2 = sel2.read_all<std::int32_t>();
    REQUIRE(v1.size() == v2.size());

    for (std::size_t i = 0; i < v1.size(); ++i) {
      CHECK(v1[i] == v2[i].to_thin());
    }
  }

  SECTION("random sample") {
    const std::uint64_t seed = 1234567890;
    const double frac = 0.75;

    const hic::File hf(path.string(), 500'000);
    const auto sel = hf.fetch("chr2L");

    std::mt19937_64 rand_eng{seed};
    const auto pixels = sel.read_all<std::int32_t>();
    std::vector<ThinPixel<std::int32_t>> subsampled_pixels_expected{};

    std::for_each(sel.template begin<std::int32_t>(), sel.template end<std::int32_t>(),
                  [&](ThinPixel<std::int32_t> p) {
                    p.count = std::binomial_distribution<std::int32_t>{p.count, frac}(rand_eng);
                    if (p.count != 0) {
                      subsampled_pixels_expected.push_back(p);
                    }
                  });

    const transformers::PixelRandomSampler sampler(sel.template begin<std::int32_t>(),
                                                   sel.template end<std::int32_t>(), frac, seed);

    const auto subsampled_pixels = sampler.read_all();
    REQUIRE(subsampled_pixels.size() == subsampled_pixels_expected.size());
    for (std::size_t i = 0; i < subsampled_pixels.size(); ++i) {
      CHECK(subsampled_pixels_expected[i] == subsampled_pixels[i]);
    }
  }
}

}  // namespace hictk::test::transformers
