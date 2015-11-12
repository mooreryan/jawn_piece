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

  Pintail! Take all the input sequences and run them one by one
  against each of the non-outgroup sequences in the database.

  For partial sequences, ignores any posistion with a ".".

  Length cutoff: If the number of non-"." mask positions is less than
  Const::WINDOW_SIZE * 2, that sequence will not be checked. TODO,
  this shouldn't use the mask length as a basis?

  Any seq that is two short by this is not included in anything
  downstream of Zintail.

  TODO: only check the sequence with highest similarity?

  TODO the mask or not mask option works, however, you need to have a
  percentile file that works with not masked sequences to actually get
  a meaningful result.

  Options:
  EOS

  opt(:queries, "input file", type: :string,
      default: File.join("test_files", "pintail.test.fa"))
  opt(:percentiles, "Percentiles file", type: :string,
      default: Const::DE_DIST_PERCENTILES)
  opt(:outdir, "Output directory", type: :string,
      default: "output")
  opt(:graph, "Output pintail graphs", type: :boolean,
      default: false)
  opt(:mask, "Mask the seqs?", type: :boolean)
end

queries = Ryan.check_file(opts[:queries], :queries)
percentiles = Ryan.check_file(opts[:percentiles], :percentiles)
Ryan.try_mkdir(opts[:outdir])

OUTPUT_FOLDER = opts[:outdir]
PINTAIL_FOLDER = File.join OUTPUT_FOLDER, "pintail"
PINTAIL_GRAPHS_FOLDER = File.join PINTAIL_FOLDER, "graphs"
# PINTAIL_GRAPHS_LOW_FOLDER = File.join PINTAIL_GRAPHS_FOLDER, "too_low"
PINTAIL_GRAPHS_HIGH_FOLDER = File.join PINTAIL_GRAPHS_FOLDER, "too_high"
NOT_CHECKED = File.join PINTAIL_FOLDER, "seqs_not_checked.txt"

# Ryan.try_mkdir PINTAIL_GRAPHS_LOW_FOLDER
Ryan.try_mkdir PINTAIL_GRAPHS_HIGH_FOLDER

LEN_CUTOFF = Const::WINDOW_SIZE * 2

mask_posns = []
db_seqs = [] # from the database
seq_num = 0
query = nil
subj = nil
n = 0
masked_db_seq_names = []
queries = {}
flagged_queries = Set.new

Ryan.log "Using database: #{Const::DATABASE}"
Ryan.log "Using percentiles: #{opts[:percentiles]}"

# get just the masked bases
FastaFile.open(Const::DATABASE).each_record do |head, seq|
  check_len seq
  check_first_head head, seq_num

  if head == "Mask"
    mask_posns = get_mask_posns seq
  # skip outgroups and ignore seqs with ambiguous bases
  elsif !head.match("outgroup") && !seq.match(/[^-actgACTG]/)
    if opts[:mask]
      db_seqs << mask_seq(mask_posns, u_to_t(seq))
    else
      db_seqs << u_to_t(seq)
    end

    masked_db_seq_names << head
  end

  seq_num += 1
end

File.open(NOT_CHECKED, "w") do |f|
  # the infile contains all the queries, assumes they are all full
  # length
  FastaFile.open(opts[:queries]).each_record do |head, seq|
    if queries.has_key? head
      abort "#{head} is repeated in #{opts[:queires]}"
    end

    masked_seq = mask_seq(mask_posns, u_to_t(seq))
    len_masked_seq_no_dots = masked_seq.gsub(".", "").length

    if len_masked_seq_no_dots >= LEN_CUTOFF
      if opts[:mask]
        queries[head] = masked_seq
      else
        queries[head] = u_to_t(seq)
      end
    else
      flagged_queries << head
      f.puts [head, len_masked_seq_no_dots].join "\t"
    end
  end
end

#### for the expected
positional_variability = get_positional_variability db_seqs
windowed_avg_probs = get_windowed_avg_probs positional_variability
database_overall_evol_dist = mean(windowed_avg_probs.map(&:last))

too_high_f =
  File.open(File.join(PINTAIL_FOLDER, "pintail.too_high.txt"), "w")
just_right_f =
  File.open(File.join(PINTAIL_FOLDER, "pintail.just_right.txt"), "w")

[too_high_f, just_right_f].each do |f|
  f.puts %w[query.name query.seq subject.name subject.seq dist de de.99th flag].
          join "\t"
end

## do the stuff with the observed
total_comparisons = queries.count * db_seqs.count.to_f
current_comparison = 0
queries.each_with_index do |(query_name, query), qi|
  db_seqs.each_with_index do |subj, si|
    subj_name = masked_db_seq_names[si]

    obs_perc_diffs = windowed_str_mismatch(query, subj)
    obs_evol_dist = mean_obs_perc_dist obs_perc_diffs

    fitting_coefficient = obs_evol_dist / database_overall_evol_dist
    exp_perc_diffs =
      windowed_avg_probs.map do |start, prob|
      [start, fitting_coefficient * prob]
    end

    de = get_de(obs_perc_diffs, exp_perc_diffs)

    # TODO read this into memory once
    pgroup = count = dist = p5 = nil
    File.open(opts[:percentiles]).each_line do |line|
      unless line.start_with? "dist"
        pgroup, count, dist, p1, p2, p3, p4, p5 = line.chomp.split "\t"

        if obs_evol_dist < dist.to_f
          break
        end
      end
    end

    if de > p5.to_f
      flag = "too_high"
      flagged_queries << query_name
      too_high_f.puts [query_name,
                       query,
                       subj_name,
                       subj,
                       "%.3f" % obs_evol_dist,
                       "%.3f" % de,
                       p5,
                       flag].join "\t"
      if opts[:graph]
        plot_diffs exp_perc_diffs, obs_perc_diffs, :too_high, query_name, subj_name
      end

      warn "WARNING: #{query_name} - #{subj_name} pair flagged too high"
    else
      flag = "just_right"
      just_right_f.puts [query_name,
                         query,
                         subj_name,
                         subj,
                         "%.3f" % obs_evol_dist,
                         "%.3f" % de,
                         p5,
                         flag].join "\t"
    end


    current_comparison += 1
    $stderr.printf "Progress: %.2f%%\r",
                   current_comparison / total_comparisons * 100
  end
end
$stderr.puts

# TODO ensure these close
too_high_f.close
just_right_f.close


File.open(File.join(Const::ASSETS_FOLDER,
                    "queries.pintail_flagged.txt"),
          "w") do |f|
  flagged_queries.sort.each do |query|
    f.puts query
  end
end
