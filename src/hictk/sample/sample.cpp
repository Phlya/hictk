// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <cista.h>
#include <fmt/format.h>

#include <boost/interprocess/ipc/message_queue.hpp>
#include <boost/process/v2.hpp>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/transformers/random_sampler.hpp"

namespace hictk::tools {

template <typename BufferT>
static void sample_writer(const SampleConfig& c) {
  boost::interprocess::message_queue queue(boost::interprocess::open_only,
                                           c.queue_identifier.c_str());
  BufferT vect{};

  constexpr auto mode = cista::mode::NONE | cista::mode::CAST;
  constexpr auto buffer_capacity = cista::serialized_size<BufferT>();
  std::vector<unsigned char> buffer(buffer_capacity);

  const cooler::File clr(c.uri);
  const transformers::PixelRandomSampler sampler(clr.begin<std::int32_t>(), clr.end<std::int32_t>(),
                                                 c.fraction, c.seed);
  std::size_t j = 0;
  for (const auto& p : sampler) {
    if (j == vect.size()) {
      buffer = cista::serialize<mode>(vect);
      queue.send(buffer.data(), buffer.size(), 0);
      j = 0;
    }
    vect[j++] = p;
  }

  std::fill(vect.begin() + static_cast<std::ptrdiff_t>(j), vect.end(), ThinPixel<std::int32_t>{});
  buffer = cista::serialize<mode>(vect);
  queue.send(buffer.data(), buffer.size(), 0);

  if (j == vect.size()) {
    vect[0] = {};
    buffer = cista::serialize<mode>(vect);
    queue.send(buffer.data(), buffer.size(), 0);
  }

  SPDLOG_DEBUG("Writer process is returning...");
}

template <typename BufferT>
static void sample_reader(const SampleConfig& c) {
  boost::interprocess::message_queue queue(boost::interprocess::open_only,
                                           c.queue_identifier.c_str());
  const auto [chroms, bin_size] = [&]() {
    const cooler::File clr(c.uri);
    return std::make_pair(clr.chromosomes(), clr.bin_size());
  }();

  auto clr = cooler::File::create(c.output_uri, chroms, bin_size, c.force);

  constexpr auto mode = cista::mode::NONE | cista::mode::CAST | cista::mode::UNCHECKED;
  constexpr auto buffer_capacity = cista::serialized_size<BufferT>();
  std::vector<unsigned char> buffer(buffer_capacity);
  std::size_t received_size{};
  unsigned int priority{};

  while (true) {
    queue.receive(buffer.data(), buffer.size(), received_size, priority);
    assert(received_size == buffer.size());
    const auto pixels = cista::offset::deserialize<BufferT, mode>(buffer);

    if (!!pixels->back()) {
      clr.append_pixels(pixels->begin(), pixels->end(), false);
    } else {
      auto first = pixels->begin();
      auto last = std::find_if(pixels->begin(), pixels->end(), [](const auto& p) { return !p; });
      clr.append_pixels(first, last);
      break;
    }
  }

  SPDLOG_DEBUG("Reader process is returning...");
}

std::string generate_queue_name() {
  std::random_device rd{};

  return fmt::format(FMT_STRING("hictk_sample_{}.queue"), rd());
}

int sample_subcmd(const SampleConfig& c) {
  using BufferT = cista::offset::array<ThinPixel<std::int32_t>, 64 * 1024>;

  if (c.spawn_reader_process) {
    assert(!c.queue_identifier.empty());
    SPDLOG_DEBUG(FMT_STRING("Spawning reader process..."));
    sample_reader<BufferT>(c);
    return 0;
  }

  if (c.spawn_writer_process) {
    assert(!c.queue_identifier.empty());
    SPDLOG_DEBUG(FMT_STRING("Spawning writer process..."));
    sample_writer<BufferT>(c);
    return 0;
  }

  const auto queue_name = generate_queue_name();

  try {
    constexpr auto buffer_capacity = cista::serialized_size<BufferT>();
    boost::interprocess::message_queue queue(boost::interprocess::create_only, queue_name.c_str(),
                                             8, buffer_capacity);

    assert(!c.argv0.empty());

    // TODO handle count
    boost::asio::io_context ctx{};
    boost::process::v2::process writer(
        ctx, c.argv0,
        {"sample", c.uri, c.output_uri, "--fraction", fmt::to_string(c.fraction), "--seed",
         fmt::to_string(c.seed), c.force ? "--force" : "", "--verbosity", "3",
         "--spawn-writer-process", "--queue-identifier", queue_name});
    boost::process::v2::process reader(
        ctx, c.argv0,
        {"sample", c.uri, c.output_uri, "--fraction", fmt::to_string(c.fraction), "--seed",
         fmt::to_string(c.seed), c.force ? "--force" : "", "--verbosity", "3",
         "--spawn-reader-process", "--queue-identifier", queue_name});

    // TODO improve error handling
    writer.wait();
    reader.wait();

  } catch (const std::exception& e) {
    boost::interprocess::message_queue::remove(queue_name.c_str());
    throw;
  }
  boost::interprocess::message_queue::remove(queue_name.c_str());
  return 0;
}

}  // namespace hictk::tools
