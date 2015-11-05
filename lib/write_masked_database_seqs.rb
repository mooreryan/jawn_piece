#!/usr/bin/env ruby

Signal.trap("PIPE", "EXIT")

methods = File.join(File.expand_path("~"), "lib", "ruby", "ryan.rb")
require_relative methods
require_relative "const"
require_relative "pintail_methods"

Ryan.req *%w[parse_fasta fail_fast parallel]

include FailFast::Assertions


opts = Trollop.options do
  banner <<-EOS

  Write the database as one sequence per line, for easy C parsing. The
  first line contains the number of sequences then a space, then total
  number of posistion per sequence.

  Mask posns file -> first line, number of posns in mask. The rest of
  the lines are the 0 based posns.

  Options:
  EOS
end

mask = get_mask_info

# File.open(Const::MASK_NAMES, "w") do |f|
#   count = mask[:names].count
#   max_len = mask[:names].map(&:length).max
#   f.printf "%s %s\n", count, max_len

#   mask[:names].each do |name|
#     f.puts name
#   end
# end

# File.open(Const::MASK_POSNS, "w") do |f|
#   f.puts [mask[:posns].count, 0].join " "
#   mask[:posns].each do |posn|
#     f.puts posn
#   end
# end

File.open(Const::MASKED_DATABASE, "w") do |f|
  num_seqs = mask[:seqs].count

  first_len = mask[:seqs].first.length

  unless mask[:seqs].all? { |seq| seq.length == first_len }
    abort "ERROR: Not all the sequences are the same length"
  end

  f.printf "%s\t%s\n", num_seqs, first_len

  mask[:seqs].each do |seq|
    f.puts seq
  end
end
