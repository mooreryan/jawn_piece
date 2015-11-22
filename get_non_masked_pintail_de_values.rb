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

  Get the DE distribution of the database.

  Options:
  EOS

  opt(:threads, "Number of threads", type: :int, default: 3)
  opt(:alignment, "Database", type: :string)
  opt(:outf, "Outfile", type: :string)
  opt(:force, "Force overwrite of the outf",
      default: true)
end

outf = Ryan.parse_fname opts[:outf]
Ryan.try_mkdir(outf[:dir])

if opts[:threads] < 1
  abort "ERROR: --threads must be >= 1"
end

if File.exists?(opts[:outf]) && !opts[:force]
  abort "ERROR: File: #{opts[:outf]} already exists!"
elsif File.exists?(opts[:outf])
  warn "WARNING: Overwriting #{opts[:outf]}!"
end

Ryan.log "Using alignment: #{opts[:alignment]}"

seqs = FastaFile.open(opts[:alignment]).to_hash.values

de_values = get_non_masked_DE_dist seqs, opts[:threads]

File.open(opts[:outf], "w") do |f|
  # f.puts %w[query subj dist de].join "\t"
  f.puts de_values.count

  de_values.sort_by { |dist, de| dist }.each do |arr|
    f.puts arr.join "\t"
  end
end

Ryan.log "Wrote #{opts[:outf]}"
