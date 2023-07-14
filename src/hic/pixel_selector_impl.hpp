// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
#include <random>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/common.hpp"
#include "hictk/hic/block_cache.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/file_reader.hpp"
#include "hictk/hic/footer.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic {

inline PixelSelector::PixelSelector(std::shared_ptr<internal::HiCFileReader> hfs_,
                                    std::shared_ptr<const internal::HiCFooter> footer_,
                                    std::shared_ptr<internal::BlockCache> cache_,
                                    std::shared_ptr<const BinTable> bins_,
                                    PixelCoordinates coords) noexcept
    : PixelSelector(std::move(hfs_), std::move(footer_), std::move(cache_), std::move(bins_),
                    coords, std::move(coords)) {}

inline PixelSelector::PixelSelector(std::shared_ptr<internal::HiCFileReader> hfs_,
                                    std::shared_ptr<const internal::HiCFooter> footer_,
                                    std::shared_ptr<internal::BlockCache> cache_,
                                    std::shared_ptr<const BinTable> bins_, PixelCoordinates coord1_,
                                    PixelCoordinates coord2_) noexcept
    : _reader(std::move(hfs_), footer_->index(), std::move(bins_), std::move(cache_)),
      _footer(std::move(footer_)),
      _coord1(std::move(coord1_)),
      _coord2(std::move(coord2_)),
      _block_idx(std::make_shared<const internal::Index::Overlap>(
          _reader.index().find_overlaps(coord1(), coord2()))) {
  for (const auto &bi : *_block_idx) {
    fmt::print(FMT_STRING("{} ({}:{})\n"), bi.id(), bi.coords().row, bi.coords().col);
  }

  //if (_footer->metadata().url.back() == '9') {
  //  assert(false);
  //}
}

inline bool PixelSelector::operator==(const PixelSelector &other) const noexcept {
  return _reader.index().chrom1() == _reader.index().chrom2() && _coord1 == other._coord1 &&
         _coord2 == other._coord2;
}

inline bool PixelSelector::operator!=(const PixelSelector &other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline auto PixelSelector::cbegin() const -> iterator<N> {
  return iterator<N>(*this);
}

template <typename N>
inline auto PixelSelector::cend() const -> iterator<N> {
  return iterator<N>::at_end(*this);
}

template <typename N>
inline auto PixelSelector::begin() const -> iterator<N> {
  return this->cbegin<N>();
}

template <typename N>
inline auto PixelSelector::end() const -> iterator<N> {
  return this->cend<N>();
}

inline SerializedPixel PixelSelector::transform_pixel(SerializedPixel pixel) const {
  const auto &c1Norm = _footer->c1Norm();
  const auto &c2Norm = _footer->c2Norm();
  const auto &expected = _footer->expectedValues();

  const auto bin1 = static_cast<std::size_t>(pixel.bin1_id);
  const auto bin2 = static_cast<std::size_t>(pixel.bin2_id);

  assert(is_inter() || bin1 <= bin2);

  const auto skipNormalization =
      normalization() == NormalizationMethod::NONE || matrix_type() == MatrixType::expected;

  if (!skipNormalization) {
    assert(bin1 < c1Norm.size());
    assert(bin2 < c2Norm.size());
    pixel.count /= static_cast<float>(c1Norm[bin1] * c2Norm[bin2]);
  }

  if (matrix_type() == MatrixType::observed) {
    return pixel;
  }

  const auto expectedCount = [&]() {
    if (is_inter()) {
      return float(_reader.avg());
    }

    const auto i = (bin2 - bin1);
    assert(i < expected.size());
    return float(expected[i]);
  }();

  if (matrix_type() == MatrixType::expected) {
    pixel.count = expectedCount;
    return pixel;
  }

  assert(matrix_type() == MatrixType::oe);
  pixel.count /= expectedCount;

  return pixel;
}

template <typename N>
inline std::vector<Pixel<N>> PixelSelector::read_all() const {
  // We push_back into buff to avoid traversing pixels twice (once to figure out the vector size,
  // and a second time to copy the actual data)
  std::vector<Pixel<N>> buff{};
  std::transform(begin<N>(), end<N>(), std::back_inserter(buff), [&](const ThinPixel<N> &p) {
    return Pixel<N>{{bins().at_hint(p.bin1_id, coord1().bin1.chrom()),
                     bins().at_hint(p.bin2_id, coord2().bin1.chrom())},
                    p.count};
  });
  return buff;
}

inline const PixelCoordinates &PixelSelector::coord1() const noexcept { return _coord1; }
inline const PixelCoordinates &PixelSelector::coord2() const noexcept { return _coord2; }
inline MatrixType PixelSelector::matrix_type() const noexcept { return metadata().matrix_type; }
inline NormalizationMethod PixelSelector::normalization() const noexcept {
  return metadata().normalization;
}
inline MatrixUnit PixelSelector::unit() const noexcept { return _reader.index().unit(); }
inline std::uint32_t PixelSelector::resolution() const noexcept {
  return _reader.index().resolution();
}

inline const Chromosome &PixelSelector::chrom1() const noexcept { return _coord1.bin1.chrom(); }
inline const Chromosome &PixelSelector::chrom2() const noexcept { return _coord2.bin1.chrom(); }

inline const std::vector<double> &PixelSelector::chrom1_norm() const noexcept {
  return _footer->c1Norm();
}
inline const std::vector<double> &PixelSelector::chrom2_norm() const noexcept {
  return _footer->c2Norm();
}

inline const BinTable &PixelSelector::bins() const noexcept { return _reader.bins(); }

inline const internal::HiCFooterMetadata &PixelSelector::metadata() const noexcept {
  assert(!!this->_footer);
  return this->_footer->metadata();
}

inline bool PixelSelector::is_intra() const noexcept { return chrom1() == chrom2(); }

inline bool PixelSelector::is_inter() const noexcept { return !is_intra(); }

template <typename N>
inline N PixelSelector::sum() const noexcept {
  return _reader.sum();
}
inline double PixelSelector::avg() const noexcept { return _reader.avg(); }

inline std::size_t PixelSelector::estimate_optimal_cache_size(
    [[maybe_unused]] std::size_t num_samples) const {
  return 10'000'000;
  /*
  if (_reader.index().empty()) {
    return 0;  // should we throw instead?
  }

  std::seed_seq sseq({_reader.index().size()});
  std::mt19937_64 rand_eng(sseq);

  // Try to guess the average block size
  std::size_t max_block_size = 0;

  std::vector<std::size_t> block_sizes{};
  std::vector<internal::BlockIndex> blocks(std::min(_reader.index().size(), num_samples));
  std::sample(_reader.index().begin(), _reader.index().end(), blocks.begin(), blocks.size(),
              rand_eng);
  for (const auto &blki : blocks) {
    auto blk = _reader.read(chrom1(), chrom2(), blki);
    if (blk) {
      max_block_size = (std::max)(blk->size(), max_block_size);
      _reader.evict(*blk);
    }
  }

  // Try to guess how many blocks overlap a single row of pixels
  std::size_t max_blocks_per_row = 0;
  const auto bin_size = bins().bin_size();

  const std::size_t first_bin_id = 0;
  const std::size_t last_bin_id =
      bins().at(coord1().bin1.chrom(), coord1().bin1.chrom().size()).rel_id() - 1;
  const auto samples = (std::min)(num_samples, bins().subset(coord1().bin1.chrom()).size());
  for (std::size_t i = 0; i < samples; ++i) {
    const auto bin_id =
        std::uniform_int_distribution<std::size_t>{first_bin_id, last_bin_id}(rand_eng);

    const auto pos1 = static_cast<std::uint32_t>(bin_id * bin_size);
    const auto bin1 = bins().at(coord1().bin1.chrom(), pos1);

    std::vector<internal::BlockIndex> buffer{};
    auto overlap = _reader.index().find_overlaps(bin1, coord2());
    const auto num_blocks = static_cast<std::size_t>(std::distance(overlap.first, overlap.last));
    max_blocks_per_row = (std::max)(max_blocks_per_row, num_blocks);
  }

  return max_blocks_per_row * max_block_size * sizeof(SerializedPixel);
   */
}

template <typename N>
inline PixelSelector::iterator<N>::iterator(const PixelSelector &sel)
    : _sel(&sel), _block_it(_sel->_block_idx->begin()), _buffer(std::make_shared<BufferT>()) {
  if (_sel->_reader.index().empty()) {
    *this = at_end(sel);
    return;
  }

  while (!!_buffer && _buffer->empty()) {
    read_next_chunk();
  }
}

template <typename N>
inline auto PixelSelector::iterator<N>::at_end(const PixelSelector &sel) -> iterator<N> {
  iterator it{};

  it._sel = &sel;
  it._buffer = nullptr;  // end of queue

  return it;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator==(const iterator &other) const noexcept {
  if (_sel != other._sel) {
    return false;
  }

  const auto bin11 = bin1_id();
  const auto bin12 = bin2_id();
  const auto bin21 = other.bin1_id();
  const auto bin22 = other.bin2_id();

  return bin11 == bin21 && bin12 == bin22;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator<(const iterator &other) const noexcept {
  const auto bin11 = bin1_id();
  const auto bin21 = other.bin1_id();

  if (bin11 != bin21) {
    return bin11 < bin21;
  }

  const auto bin12 = bin2_id();
  const auto bin22 = other.bin2_id();

  return bin12 < bin22;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator>(const iterator &other) const noexcept {
  const auto bin11 = bin1_id();
  const auto bin21 = other.bin1_id();

  if (bin11 != bin21) {
    return bin11 > bin21;
  }

  const auto bin12 = bin2_id();
  const auto bin22 = other.bin2_id();

  return bin12 > bin22;
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator*() const -> const_reference {
  assert(!!_buffer);
  assert(_buffer_i < _buffer->size());
  return (*_buffer)[_buffer_i];
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator->() const -> const_pointer {
  assert(!!_buffer);
  assert(_buffer_i < _buffer->size());
  return &(*_buffer)[_buffer_i];
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator++() -> iterator & {
  assert(!!_buffer);

  ++_buffer_i;
  while (!is_at_end() && _buffer_i >= size()) {
    read_next_chunk();
  }

  return *this;
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

template <typename N>
inline bool PixelSelector::iterator<N>::is_at_end() const noexcept {
  return _buffer == nullptr;
}

template <typename N>
inline const BinTable &PixelSelector::iterator<N>::bins() const noexcept {
  assert(!!_sel);
  return _sel->bins();
}
template <typename N>
inline const PixelCoordinates &PixelSelector::iterator<N>::coord1() const noexcept {
  assert(!!_sel);
  return _sel->coord1();
}
template <typename N>
inline const PixelCoordinates &PixelSelector::iterator<N>::coord2() const noexcept {
  assert(!!_sel);
  return _sel->coord2();
}

template <typename N>
inline std::size_t PixelSelector::iterator<N>::size() const noexcept {
  return !_buffer ? 0 : _buffer->size();
}

template <typename N>
inline std::uint64_t PixelSelector::iterator<N>::bin1_id() const noexcept {
  return !is_at_end() ? (*this)->bin1_id : std::numeric_limits<std::size_t>::max();
}
template <typename N>
inline std::uint64_t PixelSelector::iterator<N>::bin2_id() const noexcept {
  return !is_at_end() ? (*this)->bin2_id : std::numeric_limits<std::size_t>::max();
}

template <typename N>
inline void PixelSelector::iterator<N>::read_next_chunk() {
  assert(!!_sel);

  if (_block_it == _sel->_block_idx->end()) {
    *this = at_end(*_sel);
    return;
  }

  if (_buffer.use_count() != 1) {
    _buffer = std::make_shared<BufferT>(_buffer->capacity());
  }
  _buffer->clear();
  _buffer_i = 0;

  auto first_blki = *_block_it;
  while (_block_it->coords().row == first_blki.coords().row) {
    // fmt::print(FMT_STRING("Fetching {} ({}:{})\n"), _block_it->id(), _block_it->coords().row,
    //            _block_it->coords().col);
    for (auto p : *_sel->_reader.read(coord1().bin1.chrom(), coord2().bin1.chrom(), *_block_it)) {
      if (static_cast<std::size_t>(p.bin1_id) < coord1().bin1.rel_id() ||
          static_cast<std::size_t>(p.bin1_id) > coord1().bin2.rel_id() ||
          static_cast<std::size_t>(p.bin2_id) < coord2().bin1.rel_id() ||
          static_cast<std::size_t>(p.bin2_id) > coord2().bin2.rel_id()) {
        continue;
      }

      p = _sel->transform_pixel(p);
      const auto bin1_id =
          static_cast<std::size_t>(p.bin1_id) + bins().at(coord1().bin1.chrom()).id();
      const auto bin2_id =
          static_cast<std::size_t>(p.bin2_id) + bins().at(coord2().bin1.chrom()).id();
      if constexpr (std::is_integral_v<N>) {
        _buffer->emplace_back(
            ThinPixel<N>{bin1_id, bin2_id, conditional_static_cast<N>(std::round(p.count))});
      } else {
        _buffer->emplace_back(ThinPixel<N>{bin1_id, bin2_id, conditional_static_cast<N>(p.count)});
      }
    }
    if (++_block_it == _sel->_block_idx->end()) {
      break;
    }
  }

  // fmt::print(FMT_STRING("Sorting...\n"));
  std::sort(_buffer->begin(), _buffer->end());
}

inline PixelSelectorAll::PixelSelectorAll(std::vector<PixelSelector> selectors_) noexcept
    : _selectors(std::move(selectors_)) {}

template <typename N>
inline auto PixelSelectorAll::begin() const -> iterator<N> {
  return cbegin<N>();
}
template <typename N>
inline auto PixelSelectorAll::cbegin() const -> iterator<N> {
  return iterator<N>(*this);
}

template <typename N>
inline auto PixelSelectorAll::end() const -> iterator<N> {
  return cend<N>();
}
template <typename N>
inline auto PixelSelectorAll::cend() const -> iterator<N> {
  return iterator<N>{};
}

template <typename N>
inline std::vector<Pixel<N>> PixelSelectorAll::read_all() const {
  // We push_back into buff to avoid traversing pixels twice (once to figure out the vector size,
  // and a second time to copy the actual data)
  std::vector<Pixel<N>> buff{};
  std::transform(begin<N>(), end<N>(), std::back_inserter(buff), [&](const ThinPixel<N> &p) {
    return Pixel<N>{{bins().at(p.bin1_id), bins().at(p.bin2_id)}, p.count};
  });

  return buff;
}

inline MatrixType PixelSelectorAll::matrix_type() const noexcept {
  return _selectors.front().matrix_type();
}
inline NormalizationMethod PixelSelectorAll::normalization() const noexcept {
  return _selectors.front().normalization();
}
inline MatrixUnit PixelSelectorAll::unit() const noexcept { return _selectors.front().unit(); }
inline std::uint32_t PixelSelectorAll::resolution() const noexcept {
  return _selectors.front().resolution();
}

inline const BinTable &PixelSelectorAll::bins() const noexcept { return _selectors.front().bins(); }

template <typename N>
inline bool PixelSelectorAll::iterator<N>::Pair::operator<(const Pair &other) const noexcept {
  return first < other.first;
}

template <typename N>
inline bool PixelSelectorAll::iterator<N>::Pair::operator>(const Pair &other) const noexcept {
  return first > other.first;
}

template <typename N>
inline PixelSelectorAll::iterator<N>::iterator(const PixelSelectorAll &selector)
    : _selectors(std::make_shared<SelectorQueue>()),
      _its(std::make_shared<ItPQueue>()),
      _buff(std::make_shared<std::vector<ThinPixel<N>>>()) {
  std::for_each(selector._selectors.begin(), selector._selectors.end(),
                [&](const PixelSelector &sel) { _selectors->push(&sel); });

  _chrom1_id = _selectors->front()->chrom1().id();
  init_iterators();
  read_next_chunk();
}

template <typename N>
inline bool PixelSelectorAll::iterator<N>::operator==(const iterator<N> &other) const noexcept {
  if (!_buff || !other._buff) {
    return _buff == other._buff;
  }

  assert(_i < _buff->size());
  assert(other._i < other._buff->size());
  return (*_buff)[_i] == (*other._buff)[_i];
}

template <typename N>
inline bool PixelSelectorAll::iterator<N>::operator!=(const iterator<N> &other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline auto PixelSelectorAll::iterator<N>::operator*() const -> const_reference {
  assert(_buff);
  assert(_i < _buff->size());
  return (*_buff)[_i];
}

template <typename N>
inline auto PixelSelectorAll::iterator<N>::operator->() const -> const_pointer {
  return &*(*this);
}

template <typename N>
inline auto PixelSelectorAll::iterator<N>::operator++() -> iterator & {
  assert(_buff);
  if (++_i == _buff->size()) {
    read_next_chunk();
  }
  return *this;
}

template <typename N>
inline auto PixelSelectorAll::iterator<N>::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

template <typename N>
inline void PixelSelectorAll::iterator<N>::init_iterators() {
  assert(_its->empty());
  if (_selectors->empty()) {
    _buff = nullptr;
    return;
  }

  if (_its.use_count() != 1) {
    _its = std::make_shared<ItPQueue>();
  }

  if (_selectors.use_count() != 1) {
    _selectors = std::make_shared<SelectorQueue>(*_selectors);
  }

  while (!_selectors->empty() && _selectors->front()->chrom1().id() == _chrom1_id) {
    auto *sel = _selectors->front();
    _selectors->pop();
    _its->emplace(Pair{sel->begin<N>(), sel->end<N>()});
  }
}

template <typename N>
inline void PixelSelectorAll::iterator<N>::read_next_chunk() {
  if (_selectors->empty() && _its->empty()) {
    _buff = nullptr;  // signal end
    return;
  }

  if (_its->empty()) {
    ++_chrom1_id;
    init_iterators();
    return read_next_chunk();
  }

  auto [first, last] = _its->top();
  _its->pop();

  if (first == last) {
    return read_next_chunk();
  }

  if (_buff.use_count() != 1) {
    _buff = std::make_shared<std::vector<ThinPixel<N>>>();
  }
  _buff->clear();
  _i = 0;

  const auto bin1_id = first->bin1_id;
  while (first != last && first->bin1_id == bin1_id) {
    _buff->push_back(*first++);
  }
  _its->emplace(Pair{first, last});
}

}  // namespace hictk::hic
