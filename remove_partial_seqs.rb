#!/usr/bin/env ruby

Signal.trap("PIPE", "EXIT")

methods = File.join(File.expand_path("~"), "lib", "ruby", "ryan.rb")
require_relative methods
require_relative "const"
require_relative "pintail_methods"

Ryan.req *%w[parse_fasta fail_fast parallel set]

include FailFast::Assertions


opts = Trollop.options do
  banner <<-EOS

  Breaks the input file into partial and full length sequences.

  Options:
  EOS

  opt(:input, "input file", type: :string)
  opt(:outdir, "Output directory", type: :string,
      default: "output")
end

input = Ryan.check_file(opts[:input], :input)
Ryan.try_mkdir(opts[:outdir])
part_f = File.join opts[:outdir], "#{input[:base]}.partial_length.fa"
full_f = File.join opts[:outdir], "#{input[:base]}.full_length.fa"

num_full = 0
num_partial = 0
total = 0

mask_posns = []

FastaFile.open(Const::MASK).each_record do |head, seq|
  seq.each_char.with_index do |char, idx|
    if char == '*'
      mask_posns << idx
    end
  end
end

first_five = mask_posns.take 5
last_five = mask_posns.reverse.take(5).sort

File.open(part_f, "w") do |pf|
  File.open(full_f, "w") do |ff|
    FastaFile.open(opts[:input]).each_record do |head, seq|
      missing_start = first_five.all? do |posn|
        seq[posn] == "." || seq[posn] == "-"
      end

      missing_end = last_five.all? do |posn|
        seq[posn] == "." || seq[posn] == "-"
      end

      if !missing_start && !missing_end
        num_full += 1
        ff.printf ">%s\n%s\n", head, seq
      else
        pf.printf ">%s\n%s\n", head, seq
        num_partial += 1
      end
    end
  end
end

total = num_full + num_partial.to_f
s = "full: %.3f, partial: %.3f, total: %d\n" % [num_full / total,
                                                num_partial / total,
                                                total]
Ryan.log s
