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

  Pintail! Take all the input sequences and run them one by one
  against each of the 85 sequences in the database.

  Options:
  EOS

  opt(:queries, "input file", type: :string,
      default: File.join("test_files", "pintail.test.fa"))
  opt(:outdir, "Output directory", type: :string,
      default: "output/pintail")
end

queries = Ryan.check_file(opts[:queries], :queries)
Ryan.try_mkdir(opts[:outdir])

RSCRIPT  = File.join opts[:outdir], "pintail_plot.r"
RPDF     = File.join opts[:outdir], "pintail_plot.pdf"

mask_posns = []
masked_db_seqs = [] # from the database
seq_num = 0
query = nil
subj = nil
n = 0
masked_db_seq_names = []
queries = {}

# get just the masked bases
FastaFile.open(Const::DATABASE).each_record do |head, seq|
  check_len seq
  check_first_head head, seq_num

  if head == "Mask"
    mask_posns = get_mask_posns seq
  elsif !seq.match(/[^-actgACTG]/) # ignore seqs with ambiguous bases
    masked_db_seqs << mask_seq(mask_posns, u_to_t(seq))
    masked_db_seq_names << head
  end

  seq_num += 1
end

# the infile contains all the queries, assumes they are all full
# length
FastaFile.open(opts[:queries]).each_record do |head, seq|
  if queries.has_key? head
    abort "#{head} is repeated in #{opts[:queires]}"
  else
    queries[head] = mask_seq(mask_posns, u_to_t(seq))
  end
end

#### for the expected
positional_variability = get_positional_variability masked_db_seqs
windowed_avg_probs = get_windowed_avg_probs positional_variability
database_overall_evol_dist = mean(windowed_avg_probs.map(&:last))

too_low_f = File.open(File.join(opts[:outdir], "pintail.too_low.txt"), "w")
too_high_f = File.open(File.join(opts[:outdir], "pintail.too_high.txt"), "w")
just_right_f = File.open(File.join(opts[:outdir], "pintail.just_right.txt"), "w")

[too_low_f, too_high_f, just_right_f].each do |f|
  f.puts %w[query subject dist de de.lower.bound de.upper.bound flag].
          join "\t"
end

## do the stuff with the observed
total_comparisons = queries.count * masked_db_seqs.count.to_f
current_comparison = 0
queries.each_with_index do |(query_name, query), qi|
  masked_db_seqs.each_with_index do |subj, si|
    subj_name = masked_db_seq_names[si]

    obs_perc_diffs = windowed_str_mismatch(query, subj)
    obs_evol_dist = mean_obs_perc_dist obs_perc_diffs

    fitting_coefficient = obs_evol_dist / database_overall_evol_dist
    exp_perc_diffs =
      windowed_avg_probs.map do |start, prob|
      [start, fitting_coefficient * prob]
    end

    de = get_de(obs_perc_diffs, exp_perc_diffs)

    lower = upper = nil
    File.open(Const::DE_DIST_PERCENTILES).each_line do |line|
      unless line.start_with? "dist"
        less_than, _, _, lower, upper = line.chomp.split

        if obs_evol_dist < less_than.to_f
          break
        end
      end
    end

    if de < lower.to_f
      flag = "too_low"
      too_low_f.puts [query_name,
                      subj_name,
                      "%.3f" % obs_evol_dist,
                      "%.3f" % de,
                      lower,
                      upper,
                      flag].join "\t"
    elsif de > upper.to_f
      flag = "too_high"
      too_high_f.puts [query_name,
                       subj_name,
                       "%.3f" % obs_evol_dist,
                       "%.3f" % de,
                       lower,
                       upper,
                       flag].join "\t"
    else
      flag = "just_right"
      just_right_f.puts [query_name,
                         subj_name,
                         "%.3f" % obs_evol_dist,
                         "%.3f" % de,
                         lower,
                         upper,
                         flag].join "\t"
    end


    current_comparison += 1
    $stderr.printf "Progress: %.2f%%\r",
                   current_comparison / total_comparisons * 100
  end
end

# TODO ensure these close
too_low_f.close
too_high_f.close
just_right_f.close
