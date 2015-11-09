#!/usr/bin/env ruby

Signal.trap("PIPE", "EXIT")

methods = File.join(File.expand_path("~"), "lib", "ruby", "ryan.rb")
require_relative methods
require_relative "const"
require_relative "pintail_methods"

Ryan.req *%w[parse_fasta]

opts = Trollop.options do
  banner <<-EOS

  Mask the seqs.

  Options:
  EOS

  opt(:input, "input file", type: :string)
  opt(:outdir, "Output directory", type: :string,
      default: "output")
end

input = Ryan.check_file(opts[:input], :input)
Ryan.try_mkdir(opts[:outdir])

mask_posns = []
n = 0
FastaFile.open(Const::MASK).each_record do |head, seq|
  if n.zero? && head == "Mask"
    seq.each_char.with_index do |char, idx|
      mask_posns << idx if char == "*"
    end
  end
end

File.open(File.join(opts[:outdir], input[:base] + ".masked.fa"), "w") do |f|
  FastaFile.open(opts[:input]).each_record do |head, seq|
    this_seq = mask_posns.map { |posn| seq[posn] }.join

    f.printf ">%s\n%s\n", head, this_seq

    # printf "%s\n", this_seq.gsub(".", "").length
  end
end
