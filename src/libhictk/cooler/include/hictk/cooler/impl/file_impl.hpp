// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <filesystem>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5Utility.hpp>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <variant>

#include "hictk/bin_table.hpp"
#include "hictk/cooler/attribute.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/uri.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/numeric_utils.hpp"
#include "hictk/suppress_warnings.hpp"
#include "hictk/type_pretty_printer.hpp"
#include "hictk/variant_buff.hpp"

namespace hictk::cooler {

inline File::File(RootGroup entrypoint, unsigned int mode, std::size_t cache_size_bytes, double w0,
                  bool validate)
    : _mode(mode),
      _root_group(std::move(entrypoint)),
      _groups(open_groups(_root_group)),
      _datasets(open_datasets(_root_group, cache_size_bytes, w0)),
      _attrs(read_standard_attributes(_root_group)),
      _pixel_variant(detect_pixel_type(_root_group)),
      _bins(std::make_shared<BinTable>(
          import_chroms(_datasets.at("chroms/name"), _datasets.at("chroms/length"), false),
          bin_size())),
      _index(std::make_shared<Index>(init_index(_datasets.at("indexes/chrom_offset"),
                                                _datasets.at("indexes/bin1_offset"), _bins,
                                                _datasets.at("pixels/count").size(), false))) {
  assert(mode == HighFive::File::ReadOnly || mode == HighFive::File::ReadWrite);
  if (validate) {
    validate_bins();
  }
}

template <typename PixelT>
inline File::File(RootGroup entrypoint, Reference chroms, [[maybe_unused]] PixelT pixel,
                  Attributes attributes, std::size_t cache_size_bytes, double w0)
    : _mode(HighFive::File::ReadWrite),
      _root_group(std::move(entrypoint)),
      _groups(create_groups(_root_group)),
      _datasets(create_datasets<PixelT>(_root_group, chroms, cache_size_bytes, w0)),
      _attrs(std::move(attributes)),
      _pixel_variant(PixelT(0)),
      _bins(std::make_shared<const BinTable>(std::move(chroms), bin_size())),
      _index(std::make_shared<Index>(_bins)),
      _finalize(true) {
  assert(bin_size() != 0);
  assert(!_bins->empty());
  assert(!chromosomes().empty());
  assert(!_index->empty());
  assert(std::holds_alternative<PixelT>(_pixel_variant));
  write_chromosomes();
  write_bin_table();

  write_sentinel_attr();
}

template <typename PixelT>
inline File::File(RootGroup entrypoint, [[maybe_unused]] PixelT pixel, Attributes attributes,
                  std::size_t cache_size_bytes, double w0)
    : _mode(HighFive::File::ReadWrite),
      _root_group(std::move(entrypoint)),
      _attrs(std::move(attributes)),
      _pixel_variant(PixelT(0)),
      _finalize(true) {
  _groups = open_groups(_root_group);
  _datasets = open_datasets(_root_group, cache_size_bytes, w0);

  _bins = std::make_shared<BinTable>(
      import_chroms(_datasets.at("chroms/name"), _datasets.at("chroms/length"), false), bin_size());
  _index = std::make_shared<Index>(_bins);

  assert(std::holds_alternative<PixelT>(_pixel_variant));
  assert(bin_size() != 0);
  assert(!_bins->empty());
  assert(!chromosomes().empty());
  assert(!_index->empty());

  write_sentinel_attr();
}

inline File::File(std::string_view uri, std::size_t cache_size_bytes, bool validate)
    : File(open_or_create_root_group(open_file(uri, HighFive::File::ReadOnly, validate), uri),
           HighFive::File::ReadOnly, cache_size_bytes, DEFAULT_HDF5_CACHE_W0, validate) {}

inline File::File(RootGroup entrypoint, std::size_t cache_size_bytes, bool validate)
    : File(entrypoint, HighFive::File::ReadOnly, cache_size_bytes, DEFAULT_HDF5_CACHE_W0,
           validate) {}

inline File File::open_random_access(std::string_view uri, std::size_t cache_size_bytes,
                                     bool validate) {
  return File(uri, cache_size_bytes, validate);
}

inline File File::open_read_once(std::string_view uri, std::size_t cache_size_bytes,
                                 bool validate) {
  return File(open_or_create_root_group(open_file(uri, HighFive::File::ReadOnly, validate), uri),
              HighFive::File::ReadOnly, cache_size_bytes, 1.0, validate);
}

template <typename PixelT>
inline File File::create(std::string_view uri, const Reference &chroms, std::uint32_t bin_size,
                         bool overwrite_if_exists, Attributes attributes,
                         std::size_t cache_size_bytes) {
  try {
    const auto [file_path, root_path] = parse_cooler_uri(uri);
    const auto uri_is_file_path = root_path.empty() || root_path == "/";

    // URI is like myfile.mcool::/resolutions/100, but myfile.mcool does not exist
    if (!uri_is_file_path && !std::filesystem::exists(file_path)) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("parent file \"{}\" does not exist.\n"
                     "Did you forget to create the parent file with e.g. init_mcool()?"),
          uri, file_path));
    }

    // URI points to an existing file, but overwrite_if_exists=false
    if (!overwrite_if_exists && uri_is_file_path && std::filesystem::exists(file_path)) {
      throw std::runtime_error("URI points to an existing file");
    }

    auto mode = overwrite_if_exists ? HighFive::File::Overwrite : HighFive::File::Create;

    // File exists but cooler may not
    if (std::filesystem::exists(file_path) && !uri_is_file_path) {
      mode = HighFive::File::ReadWrite;
    }

    {
      auto fp = open_file(uri, mode, false);
      auto root_group = open_or_create_root_group(fp, uri);
      if (!uri_is_file_path && utils::is_cooler(root_group())) {
        if (overwrite_if_exists) {
          throw std::runtime_error(fmt::format(
              FMT_STRING("overwriting cooler nested inside .mcool or .scool is not yet supported.\n"
                         "Path to parent file: \"{}\"\""
                         "Path to nested cooler: \"{}\""),
              file_path, root_path));
        }
      }
      assert(!utils::is_cooler(root_group()));
    }

    return create<PixelT>(
        open_or_create_root_group(open_file(uri, HighFive::File::ReadWrite, false), uri), chroms,
        bin_size, attributes, cache_size_bytes);
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Cannot create cooler at the following URI: \"{}\". Reason: {}"),
                    uri, e.what()));
  }
}

inline File File::open_random_access(RootGroup entrypoint, std::size_t cache_size_bytes,
                                     bool validate) {
  return File(entrypoint, cache_size_bytes, validate);
}

inline File File::open_read_once(RootGroup entrypoint, std::size_t cache_size_bytes,
                                 bool validate) {
  return File(entrypoint, HighFive::File::ReadOnly, cache_size_bytes, 1.0, validate);
}

template <typename PixelT>
inline File File::create(RootGroup entrypoint, const Reference &chroms, std::uint32_t bin_size,
                         Attributes attributes, std::size_t cache_size_bytes) {
  static_assert(std::is_arithmetic_v<PixelT>);
  if (bin_size == 0) {
    throw std::logic_error("bin_size cannot be zero.");
  }
  attributes.bin_size = bin_size;
  try {
    if (utils::is_cooler(entrypoint())) {
      throw std::runtime_error("URI points to an already existing cooler.");
    }
    // At this point the parent file is guaranteed to exist, so we can always open it in ReadWrite
    // mode
    return File(entrypoint, chroms, PixelT(0), attributes, cache_size_bytes, true);

  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Cannot create cooler at the following URI: \"{}\". Reason: {}"),
                    entrypoint.uri(), e.what()));
  }
}

inline File::~File() noexcept {
  try {
    finalize();
  } catch (const std::exception &e) {
    fmt::print(stderr, FMT_STRING("{}\n"), e.what());
  } catch (...) {
    fmt::print(stderr,
               FMT_STRING("An unknown error occurred while closing file {}. File is likely "
                          "corrupted or incomplete."),
               path());
  }
}

inline File::operator bool() const noexcept { return !!_bins; }

inline void File::close() {
  finalize();
  *this = File{};
}

inline void File::finalize() {
  if (!_bins || !_finalize) {
    assert(!_index == !_bins);
    return;
  }

  assert(_bins);
  assert(_index);
  try {
    assert(_attrs.nnz.has_value());
    _index->set_nnz(static_cast<std::uint64_t>(*_attrs.nnz));
    write_indexes();
    write_attributes();

  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error occurred while closing file {}: {}\n"
                               "File is likely corrupted or incomplete"),
                    path(), e.what()));
  }
}

inline HighFive::File File::open_file(std::string_view uri, unsigned int mode, bool validate) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  const auto [file_path, root_grp] = parse_cooler_uri(uri);

  const auto new_file = !std::filesystem::exists(file_path);
  HighFive::File f(file_path, mode);
  if (!validate || new_file) {
    return f;
  }

  const auto status = utils::is_cooler(f, root_grp);
  if (!status) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("\"{}\" does not look like a valid Cooler file:\n"
                               "Validation report:\n{}"),
                    uri, status));
  }

  return f;
}

inline auto File::open_or_create_root_group(HighFive::File f, std::string_view uri) -> RootGroup {
  if (f.exist(parse_cooler_uri(uri).group_path)) {
    return open_root_group(f, uri);
  }
  return create_root_group(f, uri);
}

DISABLE_WARNING_PUSH
DISABLE_WARNING_UNREACHABLE_CODE
namespace internal {
template <typename Variant, std::size_t i = 0>
[[nodiscard]] inline Variant read_pixel_variant(const HighFive::DataSet &dset) {
  if constexpr (i < std::variant_size_v<Variant>) {
    using T = std::variant_alternative_t<i, Variant>;
    if (dset.getDataType() != HighFive::create_datatype<T>()) {
      return read_pixel_variant<Variant, i + 1>(dset);
    }
    return T{};
  }

  constexpr bool variant_has_monostate =
      std::is_same_v<std::monostate, std::variant_alternative_t<0, Variant>>;
  if constexpr (variant_has_monostate) {
    return std::monostate();
  } else {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unsupported type for dataset \"{}\""), dset.getPath()));
  }
}
DISABLE_WARNING_POP
}  // namespace internal

inline hictk::internal::NumericVariant File::detect_pixel_type(const RootGroup &root_grp,
                                                               std::string_view path) {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  auto dset = root_grp().getDataSet(std::string{path});
  return internal::read_pixel_variant<hictk::internal::NumericVariant>(dset);
}

template <typename N, bool cis>
inline void File::update_pixel_sum(N partial_sum) {
  static_assert(std::is_arithmetic_v<N>);

  auto &buff = cis ? _attrs.cis : _attrs.sum;
  assert(buff.has_value());
  if constexpr (std::is_floating_point_v<N>) {
    std::get<double>(*buff) += conditional_static_cast<double>(partial_sum);
  } else {
    assert(std::is_integral_v<N>);
    std::get<std::int64_t>(*buff) += conditional_static_cast<std::int64_t>(partial_sum);
  }
}

}  // namespace hictk::cooler