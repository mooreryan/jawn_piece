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

  Pintail! Jawns

  Options:
  EOS

  opt(:infile, "input file", type: :string,
      default: File.join("test_files", "pintail.test.fa"))
  opt(:database, "Input file", type: :string,
      default: File.join("test_files", "database.fa"))
  opt(:outdir, "Output directory", type: :string, default: "pintail")
end

infile = Ryan.check_file(opts[:infile], :infile)
database = Ryan.check_file(opts[:database], :database)
Ryan.try_mkdir(opts[:outdir])

RSCRIPT  = File.join opts[:outdir], "pintail_plot.r"
RPDF     = File.join opts[:outdir], "pintail_plot.pdf"
DATABASE = opts[:database]

mask_posns = []
masked_seqs = [] # from the database
seq_num = 0
query = nil
subj = nil
n = 0
masked_seq_names = []

# get subject and query
FastaFile.open(opts[:infile]).each_record do |head, seq|
  n += 1
  if n == 1
    query = seq
  elsif n == 2
    subj = seq
  end
end

unless n == 2
  warn "WARNING: There were more than two sequences in " +
       "#{opts[:infile]}, the first was set as query, the second as " +
       "subject, and the rest ignored."
end

# do the stuff with the expected percentage differences
# get just the masked bases
FastaFile.open(DATABASE).each_record do |head, seq|
  check_len seq
  check_first_head head, seq_num

  if head == "Mask"
    mask_posns = get_mask_posns seq
  elsif seq.match(/[-actgACTG]/) # ignore seqs with ambiguous bases
    masked_seqs << mask_seq(mask_posns, u_to_t(seq))
    masked_seq_names << head
  end

  seq_num += 1
end

## do the stuff with the observed

query = mask_seq(mask_posns, u_to_t(query))
subj  = mask_seq(mask_posns, u_to_t(subj))

puts "query"
puts query
puts "subj"
puts subj

obs_perc_diffs = windowed_str_mismatch(query, subj)
obs_evol_dist = mean_obs_perc_dist obs_perc_diffs

# special jawn
# all_obs_perc_diffs = []
# xs = nil
# masked_seqs.each do |subj|
#   opd = windowed_str_mismatch(query, subj)
#   xs ||= opd.map(&:first)
#   diffs = opd.map(&:last)
#   all_obs_perc_diffs << diffs
# end

# avg_obs_perc_diffs =
#   all_obs_perc_diffs.
#   transpose.
#   map.
#   with_index { |arr,idx| [xs[idx], mean(arr)] }


#### for the expected
positional_variability = get_positional_variability masked_seqs

# START algorithm 2
windowed_avg_probs = get_windowed_avg_probs positional_variability

database_overall_evol_dist = mean(windowed_avg_probs.map(&:last))
fitting_coefficient = obs_evol_dist / database_overall_evol_dist

exp_perc_diffs =
  windowed_avg_probs.map { |start, prob| [start,
                                          fitting_coefficient * prob] }

plot_diffs exp_perc_diffs, obs_perc_diffs

Ryan.log "\nDist: %.2f" % obs_evol_dist
de = get_de(obs_perc_diffs, exp_perc_diffs)
Ryan.log "DE: %.2f" % de

lower = upper = nil
File.open(Const::DE_DIST_PERCENTILES).each_line do |line|
  unless line.start_with? "dist"
    less_than, _, _, lower, upper = line.chomp.split

    if obs_evol_dist < less_than.to_f
      break
    end
  end
end

Ryan.log "lower DE bound: %.3f, upper DE bound: %.3f\n" % [lower, upper]

flag = !(lower.to_f <= de && de <= upper.to_f)

Ryan.log "Flagged? %s" % flag
