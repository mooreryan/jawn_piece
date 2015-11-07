# TODO, for percentiles, when the dist is the same, make them one rank?

methods = File.join(File.expand_path("~"), "lib", "ruby", "ryan.rb")
require_relative methods
require_relative "const"
require "fail_fast"
require "set"

include FailFast::Assertions

class PintailPair
  attr_accessor :dist,
                :de,
                :dist_percentile,
                :dist_percentile_group,
                :perentile_group_de_range,
                :min_de,
                :max_de

  # edit these
  STEP = 5
  CONFIDENCE = 95

  # dont edit these
  CUTOFFS = (STEP..100).step(STEP)
  LOW_DE_PERCENTILE = (100 - CONFIDENCE) / 2.0
  HIGH_DE_PERCENTILE = CONFIDENCE + ((100 - CONFIDENCE) / 2.0)

  def initialize dist, de
    @dist  = dist.to_f
    @de    = de.to_f
  end

  # sorts on the dist value
  def self.sort pintail_pairs
    pintail_pairs.sort_by! { |pair| pair.dist }
  end

  # takes an array of all pintail pairs
  def self.set_percentiles pintail_pairs
    self.sort pintail_pairs

    dists = pintail_pairs.map(&:dist)
    percentile_info = self.percentile_map dists

    pintail_pairs.each do |pair|
      assert percentile_info.has_key?(pair.dist)

      pair.dist_percentile = percentile_info[pair.dist]
    end
  end

  def self.set_percentile_group pintail_pairs
    self.set_percentiles pintail_pairs

    pintail_pairs.each do |pair|
      CUTOFFS.each do |cutoff|
        if pair.dist_percentile < cutoff
          pair.dist_percentile_group = cutoff - STEP
          break
        end
      end
    end
  end

  # def self.set_percentile_group pintail_pairs
  #   self.set_percentiles pintail_pairs

  #   pintail_pairs.each do |pair|
  #     if pair.dist_percentile < 10
  #       group = 0
  #     else
  #       group = (pair.dist_percentile.to_s[0] + "0").to_i
  #     end

  #     pair.dist_percentile_group = group
  #   end
  # end

  def self.check_percentiles percentiles
    gt = percentiles.any? do |_, percentile|
      percentile >= HIGH_DE_PERCENTILE
    end

    lt = percentiles.any? do |_, percentile|
      percentile <= LOW_DE_PERCENTILE
    end

    if !gt || !lt
      abort "ERROR: Not enough datapoints to get " +
            "%sth and %sth percentiles" % [LOW_DE_PERCENTILE,
                                           HIGH_DE_PERCENTILE]
    else
      :good
    end
  end

  def self.percentile_map arr
    new_arr = arr.uniq.sort
    count = new_arr.count.to_f

    percentile_info = {}

    new_arr.each_with_index do |item, idx|
      percentile_info[item] = idx / count * 100
    end

    percentile_info
  end

  def self.percentiles arr
    percentile_info = self.percentile_map arr

    percentiles = arr.map do |item|
      assert percentile_info.has_key?(item)

      [item, percentile_info[item]]
    end

    assert self.check_percentiles(percentiles) == :good

    percentiles
  end

  def self.set_percentile_group_de_range pintail_pairs
    min = nil
    max = nil
    # something wrong here
    pintail_pairs.group_by { |pair| pair.dist_percentile_group }.
      each do |group, pairs|
      de_vals_w_perc =
        self.percentiles(pairs.map { |pair| pair.de })

      de_vals_w_perc.each do |de, perc|
        if perc <= LOW_DE_PERCENTILE
          min = de
        end

        if perc >= HIGH_DE_PERCENTILE
          max = de
          break
        end
      end

      assert !min.nil?
      assert !max.nil?

      pairs.each do |pair|
        pair.min_de = min
        pair.max_de = max
      end
    end
  end
end

keys = Set.new
pintail_pairs = []
lineno = 0
File.open(Const::DE_DIST).each_line do |line|
  unless lineno.zero?
    dist, de = line.chomp.split "\t"

    # check for repeats
    # key = [query, subj].join " "
    # if keys.include? key
    #   abort "ERROR #{query} #{subj} pair exists more than once " +
    #         "in #{Const::DE_DIST}"
    # end

    # pintail_pair = PintailPair.new(query, subj, dist, de)
    pintail_pair = PintailPair.new(dist, de)
    # TODO does this actuall work?
    # if pintail_pairs.include? pintail_pair
    #   abort "ERROR #{query} #{subj} pair exists more than once " +
    #         "in #{Const::DE_DIST}"
    # end
    pintail_pairs << pintail_pair
  end
  lineno += 1
end

PintailPair.set_percentile_group pintail_pairs
PintailPair.set_percentile_group_de_range pintail_pairs

# pintail_pairs.each do |pair|
#   puts [pair.dist,
#         pair.dist_percentile,
#         pair.dist_percentile_group,
#         pair.de,
#         pair.min_de,
#         pair.max_de].join "\t"
# end

File.open(Const::DE_DIST_PERCENTILES, "w") do |f|
  f.puts ["dist.cutoff",
          "percentile.group",
          "pairs.in.percentile",
          "low.de",
          "high.de"].join "\t"
  groups = pintail_pairs.group_by { |pair| pair.dist_percentile_group }
  groups.each do |group, pairs|
    count = pairs.count
    first = pairs.sort_by { |pair| pair.dist }.reverse.first
    f.puts [first.dist, group, count, first.min_de, first.max_de].join "\t"
  end
end
