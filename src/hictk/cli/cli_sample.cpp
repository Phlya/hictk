// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <fmt/std.h>

#include <CLI/CLI.hpp>
#include <cassert>
#include <cstdint>
#include <string>

#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

void Cli::make_sample_subcommand() {
  auto& sc = *_cli.add_subcommand("sample", "Perform random sampling on cooler files.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(_config.index() == 0);
                    _config = SampleConfig{};
                  });

  _config = SampleConfig{};
  auto& c = std::get<SampleConfig>(_config);

  // clang-format off
  sc.add_option(
      "input-cooler",
      c.uri,
      // TODO support .scool
      "Path to the Cooler file to sample (URI syntax is supported).")
      ->check(IsValidCoolerFile)
      ->required();

  sc.add_option(
      "output-cooler",
       c.output_uri,
       // TODO support .scool
       "Path where to store the Cooler resulting from random sampling.")
       ->required();

  sc.add_flag(
      "-f,--force",
      c.force,
      "Force overwrite output cooler.")
      ->capture_default_str();

  sc.add_option(
      "--fraction",
      c.fraction,
      "Fraction used for random sampling.")
      ->check(CLI::PositiveNumber);

  sc.add_option(
      "--count",
      c.count,
      "Count used for random sampling.")
      ->check(CLI::PositiveNumber);

  sc.add_option(
      "--seed",
      c.seed,
      "Seed used to initialize the PRNG used for sampling.")
      ->capture_default_str();

  sc.add_option(
      "-v,--verbosity",
      c.verbosity,
      "Set verbosity of output to the console.")
      ->check(CLI::Range(1, 4))
      ->capture_default_str();

  // clang-format on

  sc.get_option("--fraction")->excludes(sc.get_option("--count"));

  // Hidden options to spawn child processes
  auto* grp = sc.add_option_group("");
  grp->add_flag("--spawn-reader-process", c.spawn_reader_process);
  grp->add_flag("--spawn-writer-process", c.spawn_writer_process);
  grp->add_option("--queue-identifier", c.queue_identifier);
  grp->get_option("--spawn-reader-process")->needs(grp->get_option("--queue-identifier"));
  grp->get_option("--spawn-writer-process")->needs(grp->get_option("--queue-identifier"));

  _config = std::monostate{};
}

void Cli::validate_sample_subcommand() const {
  assert(_cli.get_subcommand("sample")->parsed());

  std::vector<std::string> errors;
  const auto& c = std::get<SampleConfig>(_config);

  if (c.fraction == 0 && c.count == 0) {
    errors.emplace_back("Please specify either --fraction or --count.");
  }

  if (!c.force && std::filesystem::exists(c.output_uri)) {
    errors.emplace_back(fmt::format(
        FMT_STRING("Refusing to overwrite file {}. Pass --force to overwrite."), c.output_uri));
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}"),
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::transform_args_sample_subcommand() {
  auto& c = std::get<SampleConfig>(_config);
  c.argv0 = _exec_name;

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
}

}  // namespace hictk::tools
