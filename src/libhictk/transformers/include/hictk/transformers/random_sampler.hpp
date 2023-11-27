// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cstdint>
#include <memory>
#include <random>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/pixel.hpp"

namespace hictk::transformers {

template <typename PixelIt, typename RandomEngine = std::mt19937_64>
class PixelRandomSampler {
  using PixelT = typename std::iterator_traits<PixelIt>::value_type;
  using N = decltype(std::declval<PixelT>().count);
  static_assert(std::is_integral_v<N>);

  PixelIt _first{};
  PixelIt _last{};
  double _sampling_fraction{};
  RandomEngine _rand_eng{};

 public:
  class iterator;

  PixelRandomSampler(PixelIt first_pixel, PixelIt last_pixel, double fraction,
                     std::uint64_t seed = 682096614238392213);
  PixelRandomSampler(PixelIt first_pixel, PixelIt last_pixel, double fraction,
                     RandomEngine rand_eng);

  auto begin() const -> iterator;
  auto end() const -> iterator;

  auto cbegin() const -> iterator;
  auto cend() const -> iterator;

  [[nodiscard]] auto read_all() const -> std::vector<ThinPixel<N>>;

  class iterator {
    PixelIt _pixel_it{};
    PixelIt _pixel_last{};
    ThinPixel<N> _buff{};
    std::shared_ptr<const RandomEngine> _rand_eng_initial{};
    std::shared_ptr<RandomEngine> _rand_eng{};

    double _sampling_fraction{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = ThinPixel<N>;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;
    iterator(PixelIt first, PixelIt last, std::shared_ptr<const std::mt19937_64> rand_eng,
             double fraction);
    static auto at_end(PixelIt last, const RandomEngine &rand_eng, double fraction) -> iterator;

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const -> const_reference;
    [[nodiscard]] auto operator->() const -> const_pointer;

    auto operator++() -> iterator &;
    // Current implementation does not allow for an efficient implementation of it++
    auto operator++(int) -> iterator;

   private:
    [[nodiscard]] auto subsample_pixel() -> ThinPixel<N>;
  };
};

}  // namespace hictk::transformers

#include "./impl/random_sampler_impl.hpp"
