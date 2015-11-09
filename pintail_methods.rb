require_relative "const"
require_relative "c_methods"

def get_de obs_perc_diffs, exp_perc_diffs
  opd_vals = obs_perc_diffs.map(&:last)
  epd_vals = exp_perc_diffs.map(&:last)
  assert opd_vals.count == epd_vals.count

  sum = opd_vals.count.times.map do |idx|
    (opd_vals[idx] - epd_vals[idx])**2
  end.reduce(:+)

  Math.sqrt(sum / opd_vals.count.to_f)
end

def windows str
  range = (0..str.length - Const::WINDOW_SIZE)

  range.step(Const::WINDOW_STEP).map do |start|
    xposn = ((start+1) + (start+1 + Const::WINDOW_SIZE-1)) / 2.0
    [xposn, str[start, Const::WINDOW_SIZE]]
  end
end

# N doesn't match anything except N
def iupac_match? char1, char2
  char1 = char1.sub(/U/, "T").sub(/u/, "t").upcase
  char2 = char2.sub(/U/, "T").sub(/u/, "t").upcase

  if char1 == char2
    true
  else
    case char1
    when "A"
      %w[R W M D H V N].include? char2
    when "C"
      %w[Y S M B H V N].include? char2
    when "G"
      %w[R S K B D V N].include? char2
    when "T"
      %w[Y W K B D H N].include? char2
    when "R"
      %w[N A G].include? char2
    when "Y"
      %w[N C T].include? char2
    when "S"
      %w[N G C].include? char2
    when "W"
      %w[N A T].include? char2
    when "K"
      %w[N G T].include? char2
    when "M"
      %w[N A C].include? char2
    when "B"
      %w[N C G T].include? char2
    when "D"
      %w[N A G T].include? char2
    when "H"
      %w[N A C T].include? char2
    when "V"
      %w[N A C G].include? char2
    when "N"
      %w[A C T G R Y S W K M B D H V].include? char2
    else
      false
    end
  end
end

def get_bases_by_posn str1, str2
  [str1.chars, str2.chars].transpose
end

def perc_mismatch str1, str2
  assert str1.length == str2.length

  bases_by_posn = get_bases_by_posn str1, str2

  num_mismatches = bases_by_posn.
                   reject { |b1, b2| iupac_match?(b1,b2) }.count

  num_mismatches / str1.length.to_f * 100
end

# This is algorithm 1 from
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1317345/
# the index returned is 1 based
def windowed_str_mismatch str1, str2
  win          = windows(str1)
  str1_starts  = win.map(&:first)
  str1_windows = win.map(&:last)

  win          = windows(str2)
  str2_starts  = win.map(&:first)
  str2_windows = win.map(&:last)


  # TODO wat?
  [str1_windows, str2_windows].
    transpose.map.with_index do |(win1, win2), idx|
    [str1_starts[idx]+1, CMethods.perc_mismatch(win1, win2)]
  end
end

def mean arr
  arr.reduce(:+) / arr.count.to_f
end

def mean_obs_perc_dist windowed_str_mismatch
  diffs = windowed_str_mismatch.map(&:last)

  mean(diffs)
end

def rpdf which, qname, sname
  if which == :too_low
    File.join PINTAIL_GRAPHS_LOW_FOLDER, "pintail_plot.query_#{qname}.subject_#{sname}.pdf"
  elsif which == :too_high
    File.join PINTAIL_GRAPHS_HIGH_FOLDER, "pintail_plot.query_#{qname}.subject_#{sname}.pdf"
  else
    abort "ARGUMENT ERROR: must be :too_low or :too_high, got #{which.inspect}"
  end
end

def write_rscript exp_perc_diffs, obs_perc_diffs, which, qname, sname
  e_xs = exp_perc_diffs.map(&:first).join(",")
  e_ys = exp_perc_diffs.map(&:last).join(",")

  o_xs = obs_perc_diffs.map(&:first).join(",")
  o_ys = obs_perc_diffs.map(&:last).join(",")

  rscript = %Q{
pdf("#{rpdf(which, qname, sname)}", height = 8, width = 12)
e.xs <- c(#{e_xs})
e.ys <- c(#{e_ys})

o.xs <- c(#{o_xs})
o.ys <- c(#{o_ys})

max.y <- ifelse(max(o.ys) >= max(e.ys), max(o.ys), max(e.ys))

plot(ylim = c(0, 100), #c(0, max.y),
     x = e.xs,
     y = e.ys,
     col = "gray48",
     type = "l",
     xlab = "Position",
     ylab = "Percent difference",
     main = "#{qname} vs. #{sname}",
     lwd = 2)

points(o.xs, o.ys, col = "red", type = "l", lwd = 3)
points(e.xs, e.ys + 5, col = "gray70", type = "l", lwd = 1, lty = 2)
points(e.xs, e.ys - 5, col = "gray70", type = "l", lwd = 1, lty = 2)
invisible(dev.off())
}

  File.open(Const::RSCRIPT, "w") do |f|
    f.puts rscript
  end
end

def run_rscript
  Ryan.run_it "Rscript #{Const::RSCRIPT}"
end

def clean_rscript
  Ryan.run_it "rm #{Const::RSCRIPT}"
end

def plot_diffs exp_perc_diffs, obs_perc_diffs, which, qname, sname
  write_rscript exp_perc_diffs, obs_perc_diffs, which, qname, sname
  run_rscript
  clean_rscript
end

# zero based
def get_mask_posns seq
  seq.chars.map.with_index do |char, idx|
    if char == "*"
      idx
    else
      nil
    end
  end.compact
end

# apply the mask to the seq
def mask_seq mask, seq
  mask.map { |idx| seq[idx] }.join
end

def check_len seq
  if seq.length != Const::SILVA_LEN
    abort "ERROR: #{head} is not #{Const::SILVA_LEN}" +
          "Are your line endings correct?"
  end
end

def check_first_head head, seq_num
  if seq_num.zero? && head != "Mask"
    abort "ERROR: The first sequence of the database should be " +
          "'Mask', was '#{head}'"
  end
end

def u_to_t seq
  seq.gsub(/u/, "T").gsub(/U/, "T")
end

# START get the positional variability for the masked database
def get_char_freqs arr
  arr_len = arr.length.to_f

  arr.group_by(&:itself).values.map do |bases|
    bases.count / arr_len
  end
end

def most_common_residue_prob arr
  char_freqs = get_char_freqs arr

  (char_freqs.max - 0.25) / 0.75
end

def variability prob
  (1 - prob) * 100 # make it percent
end

def get_single_posn_variability arr
  prob = most_common_residue_prob arr

  variability prob
end

def get_positional_variability masked_seqs
  masked_seqs.map(&:chars).transpose.map do |bases_at_posn|
    get_single_posn_variability bases_at_posn
  end
end
# END get the positional variability for the masked database

def get_windowed_avg_probs positional_variability
  windows(positional_variability).map do |start, window|
    [start, mean(window)]
  end
end

# get the mask database seqs
def get_mask_info db
  seq_num          = 0
  mask_posns       = []
  masked_seqs      = []
  masked_seq_names = []

  FastaFile.open(db).each_record do |head, seq|
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

  { posns: mask_posns, seqs: masked_seqs, names: masked_seq_names }
end

def get_database_pintail_info masked_seqs
  positional_variability = get_positional_variability masked_seqs

  windowed_avg_probs = get_windowed_avg_probs positional_variability

  database_overall_evol_dist = mean(windowed_avg_probs.map(&:last))

  [windowed_avg_probs, database_overall_evol_dist]
end

def get_DE_dist mask, threads
  # get a distribution of DE values
  names_combos = mask[:names].combination(2).to_a
  total = names_combos.count.to_f

  windowed_avg_probs, database_overall_evol_dist = get_database_pintail_info mask[:seqs]

  de_values =
    Parallel.map_with_index(mask[:seqs].combination(2),
                            in_processes: threads) do |(query, subj), idx|
    # mask[:seqs].combination(2).map.with_index do |(query, subj), idx|
    progress = (idx+1) / total * 100
    # $stderr.printf "Progress: %.2f%%\r", progress

    obs_perc_diffs =
      CMethods.windowed_str_mismatch_wrapper(query, subj)
    obs_evol_dist = mean_obs_perc_dist obs_perc_diffs

    fitting_coefficient = obs_evol_dist / database_overall_evol_dist
    exp_perc_diffs = windowed_avg_probs.map do |start, prob|
      [start, fitting_coefficient * prob]
    end

    query_name = names_combos[idx].first
    subj_name  = names_combos[idx].last

    [#query_name,
     #subj_name,
     obs_evol_dist,
     CMethods.get_de_wrapper(obs_perc_diffs, exp_perc_diffs)]
  end
  # $stderr.puts

  de_values
end
