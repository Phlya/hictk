// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include "hictk/pixel.hpp"

namespace hictk::transformers {

template <typename PixelIt, typename RandomEngine>
inline PixelRandomSampler<PixelIt, RandomEngine>::PixelRandomSampler(PixelIt first_pixel,
                                                                     PixelIt last_pixel,
                                                                     double fraction,
                                                                     std::uint64_t seed)
    : PixelRandomSampler(first_pixel, last_pixel, fraction, RandomEngine{seed}) {}

template <typename PixelIt, typename RandomEngine>
inline PixelRandomSampler<PixelIt, RandomEngine>::PixelRandomSampler(PixelIt first_pixel,
                                                                     PixelIt last_pixel,
                                                                     double fraction,
                                                                     RandomEngine rand_eng)
    : _first(std::move(first_pixel)),
      _last(std::move(last_pixel)),
      _sampling_fraction(fraction),
      _rand_eng(std::move(rand_eng)) {}

template <typename PixelIt, typename RandomEngine>
inline auto PixelRandomSampler<PixelIt, RandomEngine>::begin() const -> iterator {
  return cbegin();
}

template <typename PixelIt, typename RandomEngine>
inline auto PixelRandomSampler<PixelIt, RandomEngine>::cbegin() const -> iterator {
  return {_first, _last, std::make_shared<RandomEngine>(_rand_eng), _sampling_fraction};
}

template <typename PixelIt, typename RandomEngine>
inline auto PixelRandomSampler<PixelIt, RandomEngine>::end() const -> iterator {
  return cend();
}

template <typename PixelIt, typename RandomEngine>
inline auto PixelRandomSampler<PixelIt, RandomEngine>::cend() const -> iterator {
  return iterator::at_end(_last, _rand_eng, _sampling_fraction);
}

template <typename PixelIt, typename RandomEngine>
inline auto PixelRandomSampler<PixelIt, RandomEngine>::read_all() const
    -> std::vector<ThinPixel<N>> {
  std::vector<ThinPixel<N>> buffer{};
  std::copy(begin(), end(), std::back_inserter(buffer));
  return buffer;
}

template <typename PixelIt, typename RandomEngine>
inline PixelRandomSampler<PixelIt, RandomEngine>::iterator::iterator(
    PixelIt first, PixelIt last, std::shared_ptr<const std::mt19937_64> rand_eng, double fraction)
    : _pixel_it(first),
      _pixel_last(std::move(last)),
      _rand_eng_initial(std::move(rand_eng)),
      _rand_eng(std::make_shared<RandomEngine>(*_rand_eng_initial)),
      _sampling_fraction(fraction) {
  while (_pixel_it != _pixel_last) {
    _buff = subsample_pixel();
    if (_buff.count == 0) {
      ++_pixel_it;
    } else {
      break;
    }
  }
}

template <typename PixelIt, typename RandomEngine>
inline auto PixelRandomSampler<PixelIt, RandomEngine>::iterator::at_end(
    PixelIt last, const RandomEngine &rand_eng, double fraction) -> iterator {
  return {last, last, std::make_shared<const RandomEngine>(rand_eng), fraction};
}

template <typename PixelIt, typename RandomEngine>
inline bool PixelRandomSampler<PixelIt, RandomEngine>::iterator::operator==(
    const iterator &other) const noexcept {
  return _pixel_it == other._pixel_it && _sampling_fraction == other._sampling_fraction &&
         *_rand_eng_initial == *other._rand_eng_initial;
}

template <typename PixelIt, typename RandomEngine>
inline bool PixelRandomSampler<PixelIt, RandomEngine>::iterator::operator!=(
    const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename PixelIt, typename RandomEngine>
inline auto PixelRandomSampler<PixelIt, RandomEngine>::iterator::operator*() const
    -> const_reference {
  return _buff;
}

template <typename PixelIt, typename RandomEngine>
inline auto PixelRandomSampler<PixelIt, RandomEngine>::iterator::operator->() const
    -> const_pointer {
  return &(**this);
}

template <typename PixelIt, typename RandomEngine>
inline auto PixelRandomSampler<PixelIt, RandomEngine>::iterator::operator++() -> iterator & {
  while (++_pixel_it != _pixel_last) {
    _buff = subsample_pixel();
    if (_buff.count != 0) {
      break;
    }
  }
  return *this;
}

template <typename PixelIt, typename RandomEngine>
inline auto PixelRandomSampler<PixelIt, RandomEngine>::iterator::operator++(int) -> iterator {
  auto it = *this;
  it._rand_eng = std::make_shared<RandomEngine>(*_rand_eng);
  std::ignore = ++(*this);
  return it;
}

template <typename PixelIt, typename RandomEngine>
inline auto PixelRandomSampler<PixelIt, RandomEngine>::iterator::subsample_pixel() -> ThinPixel<N> {
  assert(_rand_eng);
  if (_pixel_it == _pixel_last) {
    return {};
  }

  auto pxl = *_pixel_it;
  pxl.count = std::binomial_distribution<N>{pxl.count, _sampling_fraction}(*_rand_eng);
  return pxl;
}

}  // namespace hictk::transformers
