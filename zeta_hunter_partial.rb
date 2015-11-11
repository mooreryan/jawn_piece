#!/usr/bin/env ruby

def parse_dist_f dist_f
  lineno = 1
  num_seqs = 0
  name = nil
  dists = []
  first_seq = true
  dist_info = {}
  seq_names = {}
  seqno = 0

  File.open(dist_f).each_line do |line|
    if lineno == 1
      num_seqs = line.strip.to_i
    else
      if line.match(/^[a-z]+/) && lineno == 2 # first seq start
        arr = line.chomp.split
        name = arr.shift
        seq_names[seqno] = name
        seqno += 1
        dists << arr.map { |s| s.to_f }
      elsif line.match(/^[a-z]+/)
        # add the previous info to the hash
        dists.flatten!
        assert dists.count == num_seqs

        dist_info[name] = dists.flatten

        # get this info
        arr = line.chomp.split
        name = arr.shift
        seq_names[seqno] = name
        seqno += 1
        dists = arr.map { |s| s.to_f }
      elsif line.start_with? " "
        arr = line.chomp.split
        dists << arr.map { |s| s.to_f }
      else
        abort "ERROR: Something looks off in #{ARGV.first}"
      end
    end

    lineno += 1
  end

  return dist_info, seq_names
end

def mean arr
  arr.reduce(:+) / arr.count
end

Signal.trap("PIPE", "EXIT")

methods = File.join(File.expand_path("~"), "lib", "ruby", "ryan.rb")
require_relative methods
require_relative "const"

Ryan.req *%w[parse_fasta set fail_fast]

include FailFast::Assertions

opts = Trollop.options do
  banner <<-EOS

  ZetaHunter

  Breaks with SEANIE line endings.

  Options:
  EOS

  opt(:alignment,
      "Input alignment from SILVA",
      type: :string,
      default: "test_files/zetas_aligned.fasta")
  opt(:database, "MSA database",
      type: :string,
      default: File.join(File.dirname(__FILE__),
                         "assets",
                         "database.fa"))
  opt(:db_otu_info,
      "Accession to OTU info for the database",
      type: :string,
      default: File.join(File.dirname(__FILE__),
                         "assets",
                         "otus.txt"))
  opt(:outdir,
      "Output directory",
      type: :string,
      default: "output")
end

alignment = Ryan.check_file(opts[:alignment], :alignment)
database = Ryan.check_file(opts[:database], :database)
Ryan.try_mkdir(opts[:outdir])

# use this in place of opts[:alignment]
queries_no_chimeras = File.join opts[:outdir], "queries_no_chimeras.fa"
queries_no_chimeras_full = File.join opts[:outdir], "queries_no_chimeras.full_length.fa"
queries_no_chimeras_partial = File.join opts[:outdir], "queries_no_chimeras.partial_length.fa"

masked_seqs = File.join opts[:outdir], "masked_alignment.fa"
degapped_seqs = File.join opts[:outdir], "degapped_alignment.fa"
degapped_mask = File.join opts[:outdir], "degapped_mask.fa"
masked_for_phylip = File.join opts[:outdir], "masked_for_phylip.fa"
for_phylip_map = File.join opts[:outdir], "for_phylip_map.txt"
phylip_infile = File.join opts[:outdir], "DNAdist_params.txt"
this_dir = File.dirname(__FILE__)
phylip_outfile = File.join this_dir, "outfile"
masked_dist = File.join opts[:outdir], "masked_aln.dist"
dnadist = File.join this_dir, "bin", "osx", "dnadist"
dotur = File.join this_dir, "bin", "osx", "dotur"
dotur_outdir = File.join opts[:outdir], "dotur"
otu_calls = File.join dotur_outdir, "de_novo.fn.list"
closed_ref_otu_info = File.join opts[:outdir], "closed_ref_otu_info.txt"
closed_ref_otu_calls = File.join opts[:outdir], "closed_ref_otu_calls.txt"
de_novo_dist = File.join opts[:outdir], "de_novo.dist"
de_novo_otu_calls = File.join opts[:outdir], "de_novo_otu_calls.txt"
final_otu_calls = File.join opts[:outdir], "final_otu_calls.txt"
tree = File.join opts[:outdir], "fasttree.tre"
tree_nex = File.join opts[:outdir], "fasttree.nex"
possible_non_zetas = File.join opts[:outdir], "possible_non_zetas.txt"
color_tree_name_map = File.join opts[:outdir], "color_tree_name_map.txt"
color_tree_colors = File.join opts[:outdir], "color_tree_colors.txt"

# remove phylip outfile if it exists
if File.exists? phylip_outfile
  warn "Removing #{phylip_outfile}"
  Ryan.run_it "rm #{phylip_outfile}"
end

# external scripts
pintail = File.join this_dir, "pintail_2.rb"
generate_report = File.join this_dir, "generate_html_report.rb"
separate_partial_seqs = File.join this_dir, "remove_partial_seqs.rb"
fasttree = File.join this_dir, "bin", "osx", "FastTree"
flag_non_zetas = File.join this_dir, "flag_non_zetas.rb"
color_tree = File.join this_dir, "bin", "color_tree"

Ryan.try_mkdir dotur_outdir

SILVA_ALN_LEN = 50000
OTU_LEVEL = 0.03

mask = ""
mask_posns = []
partial_mask_posns = []
n = 0
database = {}
gap_posns = []
db_otu_info = {}
new_to_orig_name = {} # used to prevent names from being too long for phylip
otu_calls_info = {}
user_provided_headers = Set.new
outgroups = Set.new
which_cutoff = -1
which_dotur_otu_group = {}
otu_group_actual_otus = {}
seqs_to_cluster_de_novo = nil
dist_info = nil
idx_to_seq_name = nil
closed_otu_info = nil
flagged_seqs = Set.new
entropy = []
total_entropy = 0
query_entropy = {}

Ryan.time_it("Read outgroups") do
  File.open(Const::OUTGROUPS).each_line do |line|
    outgroups << line.chomp
  end
end

Ryan.time_it("Read entropy info") do
  # TODO if file doesn't exist, create it
  n = 0
  File.open(Const::ENTROPY).each_line do |line|
    posn, ent = line.chomp.split "\t"

    unless posn.to_i == n
      abort "ERROR: #{Const::ENTROPY} is malformed"
    end

    entropy[posn.to_i] = ent.to_f

    n += 1
  end

  assert n == 1282
end

total_entropy = entropy.reduce(:+)

Ryan.time_it("Check for chimeras") do
  # TODO run this every time?
  # if !File.exists? Const::FLAGGED_SEQS
  Ryan.time_it("Run Zintail") do
    cmd = "ruby #{pintail} " +
          "--queries #{opts[:alignment]} " +
          "--outdir #{opts[:outdir]}"
    Ryan.run_it cmd
    # end
  end

  File.open(Const::FLAGGED_SEQS).each_line do |line|
    flagged_seqs << line.chomp
  end
end

def dot_to_dash seq
  seq.gsub ".", "-"
end

Ryan.time_it("Remove chimeras") do
  File.open(queries_no_chimeras, "w") do |f|
    FastaFile.open(opts[:alignment]).each_record do |head, seq|
      unless flagged_seqs.include? head
        f.printf ">%s\n%s\n", head, dot_to_dash(seq)
      end
    end
  end
end

Ryan.time_it("Separate full and partial seqs") do
  cmd = "ruby #{separate_partial_seqs} " +
        "--input #{queries_no_chimeras} " +
        "--outdir #{opts[:outdir]}"
  Ryan.run_it cmd
end

Ryan.time_it("Get database OTU metadata info") do
  lineno = 0
  File.open(opts[:db_otu_info]).each_line do |line|
    unless lineno.zero?
      acc, otu, clone, num = line.chomp.split "\t"

      if db_otu_info.has_key? acc
        abort "ERROR: #{acc} is repeated in #{opts[:db_otu_info]}"
      end

      db_otu_info[acc] = { otu: otu, clone: clone, num: num.to_i }
    end

    lineno += 1
  end
end

Ryan.time_it("Read the database for gaps and mask") do
  n = 0
  FastaFile.open(opts[:database]).each_record do |head, seq|
    # ensure each sequence in the database is the proper length
    unless seq.length == SILVA_ALN_LEN
      abort "ERROR: #{head} in #{opts[:database]} is #{seq.length} bases, should be #{SILVA_ALN_LEN}"
    end

    if n.zero? # in the mask
      mask = seq

      # read in the positions to keep into mask_posns array
      mask.each_char.with_index do |c, i|
        mask_posns << i if c == '*'
      end
    else
      # ensure the sequence has OTU info in the metadata file
      unless db_otu_info.has_key? head
        abort "ERROR: #{head} in #{opts[:database]} is not present in #{opts[:db_otu_info]}"
      end

      # ensure sequence headers are unique in the database
      if database.has_key? head
        abort "ERROR: #{head} is duplicated in #{opts[:database]}"
      end

      database[head] = seq.gsub(/U/, "T").gsub(/u/, "t")

      # note position of gaps
      these_gap_posns = Set.new
      seq.each_char.with_index do |c, i|
        these_gap_posns << i if c == '-'
      end

      # gaps array contians a Set for each sequence containing the
      # posn of all gaps for that sequence
      gap_posns << these_gap_posns
    end

    n += 1
  end
end

# TODO we trim the mask even for full sequences since the partial seqs
# script counts everything that has less than five bases missing from
# each end
Ryan.time_it("Update gap posns from user seqs") do
  # first run through the sequences to adjust the mask
  FastaFile.open(queries_no_chimeras).each_record do |head, seq|
    unless seq.length == SILVA_ALN_LEN
      warn "ERROR: #{head} in #{queries_no_chimeras} is #{seq.length} bases, should be #{SILVA_ALN_LEN}"
    end

    if user_provided_headers.include? head
      abort "ERROR: #{head} duplicated in #{queries_no_chimeras}"
    end

    user_provided_headers << head

    these_gap_posns = Set.new
    seq.each_char.with_index do |c, i|
      these_gap_posns << i if c == '-'
    end

    gap_posns << these_gap_posns

    first_seq_posn = seq.index /[^-\.]/
    last_seq_posn = seq.length - seq.reverse.index(/[^-\.]/) - 1

    # while first_seq_posn > mask_posns.first
    #   # the sequence starts after the first mask posn
    #   # trim the mask
    #   mask_posns.shift
    # end

    # while last_seq_posn < mask_posns.last
    #   mask_posns.pop
    # end
  end
end

Ryan.time_it("Apply mask to seqs & write") do
  # apply the mask, figure out gaps
  File.open(masked_seqs, "w") do |f|
    FastaFile.open(queries_no_chimeras).each_record do |head, seq|
      masked_seq = mask_posns.map { |i| seq[i] }.join("").gsub(/U/, "T").gsub(/u/, "t")

      f.printf ">%s\n%s\n", head, masked_seq
    end

    database.each do |head, seq|
      masked_seq = mask_posns.map { |i| seq[i] }.join("")

      f.printf ">%s\n%s\n", head, masked_seq
    end
  end
end

def get_seq_entropy entropy, seq
  non_zero_posns = Set.new
  seq_entropy = seq.each_char.map.with_index do |char, idx|
    if char.match(/[^-\.]/)
      non_zero_posns << idx
      entropy[idx]
    else
      0
    end
  end

  assert seq_entropy.count == 1282

  total_entropy = entropy.reduce(:+).to_f

  { entropy: (seq_entropy.reduce(:+) / total_entropy * 100).round(3),
    non_zero_posns: non_zero_posns }
end


Ryan.time_it("Get entropy for masked seqs") do
  FastaFile.open(masked_seqs).each_record do |head, seq|
    query_entropy[head] = get_seq_entropy entropy, seq
  end
end

Ryan.time_it("Write degapped alignment & mask") do

  shared_gap_posns = gap_posns.reduce(:&)

  File.open(degapped_mask, "w") do |f|
    FastaFile.open(Const::MASK).each_record do |head, seq|
      f.printf ">%s\n", head
      seq.each_char.with_index do |c, i|
        f.print c unless shared_gap_posns.include? i
      end
      f.print "\n"
    end
  end

  File.open(degapped_seqs, "w") do |f|
    FastaFile.open(queries_no_chimeras).each_record do |head, seq|
      seq = seq.gsub(/U/, "T").gsub(/u/, "t")

      f.printf ">%s\n", head

      seq.each_char.with_index do |c, i|
        f.print c unless shared_gap_posns.include? i
      end
      f.print "\n"
    end

    database.each do |head, seq|
      f.printf ">%s\n", head

      seq.each_char.with_index do |c, i|
        f.print c unless shared_gap_posns.include? i
      end
      f.print "\n"
    end
  end
end

Ryan.time_it("Make masked_for_phylip.fa") do

  aln_len = nil
  num_seqs = 0
  FastaFile.open(masked_seqs).each_record do |head, seq|

    num_seqs += 1

    # all alignments must have same length
    # TODO move this to step 3
    if aln_len && seq.length != aln_len
      abort("Error: alignments in #{masked_seqs} are not the same length")
    end

    aln_len ||= seq.length
  end

  # write the phylip file and the map

  File.open(masked_for_phylip, "w") do |phy_aln|
    phy_aln.printf "%s %s\n", num_seqs, aln_len

    File.open(for_phylip_map, "w") do |phy_map|

      new_name = "a"
      FastaFile.open(masked_seqs).each_record do |head, seq|
        if user_provided_headers.include? head
          type = "query"
          col = "blue"
        elsif outgroups.include? head
          type = "outgroup"
          col = "red"
        elsif (database.keys - outgroups.to_a).include? head
          type = "database"
          col = "green"
        else
          abort "ERROR: #{head} can't be grouped"
        end

        phy_map.puts [new_name, head, type, col].join "\t"
        new_to_orig_name[new_name] = head

        phy_aln.printf "%-10.10s%s\n", new_name, seq

        new_name = new_name.next
      end
    end
  end
end

Ryan.time_it("Build tree") do
  cmd = "#{fasttree} -nt #{masked_for_phylip} > #{tree}"
  Ryan.run_it cmd
end

Ryan.time_it("Color tree") do
  Ryan.time_it("Make color_tree name map") do
    File.open(color_tree_colors, "w") do |cf|
      File.open(color_tree_name_map, "w") do |nmf|
        File.open(for_phylip_map).each_line do |line|
          new_name, old_name, type, col = line.chomp.split "\t"

          nmf.puts [new_name, old_name].join "\t"

          cf.puts [old_name, col].join "\t"
        end
      end
    end
  end
  cmd = "#{color_tree} -bte -n #{color_tree_name_map} " +
        "-p #{color_tree_colors} #{tree} > #{tree_nex}"
  Ryan.run_it cmd
end

Ryan.time_it("Check for potentially bad sequences") do
  cmd = "ruby #{flag_non_zetas} #{tree} #{for_phylip_map} " +
        "> #{possible_non_zetas}"
  Ryan.run_it cmd
end

Ryan.time_it("Write the DNAdist params file") do
  # write the phylip input file

  File.open(phylip_infile, 'w') do |f|
    f.puts masked_for_phylip
    f.puts "I"
    f.puts "Y"
  end
end

Ryan.time_it("DNAdist") do
  #TODO this will die if infile/outfile is names is messed up TODO
  #phylip errors are weird and Ryan.run_it doesn't realize its an
  #error
  Ryan.run_it "#{dnadist} < #{phylip_infile}"
  Ryan.run_it "mv #{phylip_outfile} #{masked_dist}"
end

Ryan.time_it("Closed reference OTU assignment") do
  closed_otu_info = {}

  dist_info, idx_to_seq_name = parse_dist_f masked_dist

  # only check user provided sequences
  dist_info = dist_info.select do |name, dists|
    user_provided_headers.include? new_to_orig_name[name]
  end

  # go through and check which of the distsances are less than 0.03
  # and then check if those are in the database are not to assign the
  # closed reference OTUs
  dist_info.each_with_index do |(name, dists), idx1|
    dists.each_with_index do |dist, idx2|
      query_seq_mapped_name = name
      query_seq_name = new_to_orig_name[query_seq_mapped_name]

      hit_seq_mapped_name = idx_to_seq_name[idx2]
      hit_seq_name = new_to_orig_name[hit_seq_mapped_name]

      if dist.to_f <= Const::OTU_CUTOFF
        # if the hit_seq_name is in the original database
        if db_otu_info.has_key? hit_seq_name
          the_otu = db_otu_info[hit_seq_name][:otu]

          if closed_otu_info.has_key? query_seq_name
            if closed_otu_info[query_seq_name].has_key? the_otu
              closed_otu_info[query_seq_name][the_otu] << dist
            else
              closed_otu_info[query_seq_name][the_otu] = [dist]
            end
          else
            closed_otu_info[query_seq_name] = { the_otu => [dist] }
          end
        else # the hit_seq_name is NOT in the original database
          the_otu = "USR"
        end
      end
    end
  end

  File.open(closed_ref_otu_info, "w") do |f|
    f.puts %w[query otu otu.size min.dist max.dist mean.dist percent.entropy].join "\t"
    closed_otu_info.each do |query_seq_name, otu_info|
      otu_info.each do |otu, dists|
        assert query_entropy[query_seq_name]

        f.puts [query_seq_name,
                otu,
                dists.count,
                dists.min,
                dists.max,
                mean(dists),
                query_entropy[query_seq_name][:entropy]].join "\t"

      end
    end
  end

  File.open(closed_ref_otu_calls, "w") do |f|
    f.puts %w[query otu perc.entropy].join "\t"
    closed_otu_info.each do |query_seq_name, otu_info|
      if otu_info.count > 1
        # TODO this will break if there is a tie in the distance
        # between two OTUs -> consider using count to break
        # ties. Currently it will just take the last one
        min_dists = otu_info.map { |_, dists| dists.min }

        # pick min dist
        min = 2
        min_idx = -1
        min_dists.each_with_index do |dist, idx|
          if dist < min
            min = dist
            min_idx = idx
          end
        end

        assert min != 2
        assert min_idx != -1

        otu = otu_info.to_a[min_idx].first
      else
        otu = otu_info.to_a.first.first
      end

      assert query_entropy[query_seq_name]
      f.puts [query_seq_name,
              otu,
              query_entropy[query_seq_name][:entropy]].join "\t"
    end
  end
end

Ryan.time_it("Write the dist file for de novo clustering") do
  # note which user seqs need to be clustered de novo
  seqs_to_remove = Set.new closed_otu_info.keys
  seqs_to_cluster_de_novo = user_provided_headers - seqs_to_remove
  puts

  # puts "To be clustered de novo"
  # seqs_to_cluster_de_novo.each do |name|
  #   puts name
  # end

  # get the index of the de novo seqs as they appear in the dist array
  # as returned by parse_dist_f
  assert(idx_to_seq_name.values.count ==
         idx_to_seq_name.values.uniq.count)
  seq_name_to_idx = idx_to_seq_name.invert

  lines = []
  seqs_to_cluster_de_novo.each do |name|
    # attach the seq name to the dist
    dists = dist_info[new_to_orig_name.invert[name]].map.with_index do |dist, idx|
      [dist, new_to_orig_name[idx_to_seq_name[idx]]]
    end

    # keep only seqs that didn't hit the database
    dists.select! do |dist, seq_name|
      seqs_to_cluster_de_novo.include? seq_name
    end

    assert dists.count == seqs_to_cluster_de_novo.count

    lines << dists
  end

  # puts
  # puts "Dist matrix for de novo clustering"

  File.open(de_novo_dist, "w") do |f|
    f.printf "   %d\n", lines.count
    lines.each_with_index do |dists, idx|
      name = dists[idx].last
      these_dists = dists.map { |dist, name| dist }.join(" ")
      f.printf "%-10.10s %s\n", new_to_orig_name.invert[name], these_dists
    end
  end
end

Ryan.time_it("DOTUR on the de novo seqs") do
  Ryan.run_it "#{dotur} #{de_novo_dist}"
  Ryan.run_it "mv #{de_novo_dist.sub(/dist$/, "fn")}.* #{dotur_outdir}"
end

Ryan.time_it("Read DOTUR OTU calls") do
  File.open(otu_calls).each_line do |line|
    unless line.start_with? "unique"
      cutoff, num, *rest = line.chomp.split "\t"

      otus = []

      rest.each do |otu_group|
        otus << otu_group.split(",").map do |header|
          assert new_to_orig_name[header]

          new_to_orig_name[header]
        end
      end

      otu_calls_info[cutoff.to_f] = otus
    end
  end
end


Ryan.time_it("Pick cutoff closest to #{OTU_LEVEL}") do
  # which cutoff to use? get closest to OTU_LEVEL without going over
  otu_calls_info.each do |cutoff, otu_group|
    # TODO optimize this with a break

    if cutoff <= OTU_LEVEL && cutoff >= which_cutoff
      which_cutoff = cutoff
    end
  end

  # ensure cutoff actually was set
  assert which_cutoff != -1
end

Ryan.time_it("Report de novo OTU calls") do
  File.open(de_novo_otu_calls, "w") do |f|
    f.puts %w[query otu seq.perc.entropy otu.group.shared.entropy].join "\t"
    otu_calls_info[which_cutoff].each_with_index do |otu_group, otu_num|

      intersection = otu_group.map do |seq_name|
        assert query_entropy[seq_name]
        query_entropy[seq_name][:non_zero_posns] # a set of non zero posns
      end.reduce(&:intersection)

      if intersection.empty?
        group_shared_entropy = 0.0
      else
        group_shared_entropy = intersection.map do |shared_non_zero_posn|
          assert entropy[shared_non_zero_posn]
          entropy[shared_non_zero_posn] # the entropy at that posn
        end.reduce(:+) / total_entropy.to_f * 100
      end

      otu_group.each do |name|
        assert query_entropy[name]
        f.printf "%s\tnew.otu.%d\t%s\t%s\n",
                 name,
                 otu_num + 1,
                 query_entropy[name][:entropy],
                 group_shared_entropy.round(3)
      end
    end
  end
end

Ryan.time_it("Write final OTU calls") do
  File.open(final_otu_calls, "w") do |f|
    queries = {}
    File.open(closed_ref_otu_calls).each_line do |line|
      unless line.start_with? "query"
        query, otu, pent = line.chomp.split "\t"

        if queries.has_key? query
          abort "ERROR: #{query} is repeated"
        else
          queries[query] = [otu, pent, "NA"]
        end
      end
    end

    File.open(de_novo_otu_calls).each_line do |line|
      unless line.start_with? "query"
        query, otu, pent, group_ent = line.chomp.split "\t"

        if queries.has_key? query
          abort "ERROR: #{query} is repeated"
        else
          queries[query] = [otu, pent, group_ent]
        end
      end
    end

    f.puts %w[query otu perc.entropy group.ent].join "\t"
    queries.each do |query, otu|
      f.puts [query, otu].flatten.join "\t"
    end
  end
end

Ryan.time_it("Generate HTML report") do
  cmd = "ruby #{generate_report} --directory #{opts[:outdir]}"
  Ryan.run_it cmd
end

# Ryan.time_it("10 Assign OTU numbers to OTU groups from the " +
#              "#{opts[:db_otu_info]}") do
#   assert otu_calls_info[which_cutoff]

#   # otu_group is an array of headers contained in that OTU
#   #
#   # the idx is the number of that OTU from DOTUR, could be different
#   # each run
#   otu_calls_info[which_cutoff].each_with_index do |otu_group, idx|
#     otu_group_actual_otus[idx] = []

#     otu_group.each do |header|
#       if which_dotur_otu_group.has_key? header
#         abort "ERROR: #{header} repeated in the OTU groups"
#       else
#         which_dotur_otu_group[header] = idx
#         # if the sequence is in the database
#         if db_otu_info[header]
#           otu_group_actual_otus[idx] << db_otu_info[header][:otu]
#         end
#       end
#     end
#   end
# end

# otu_group_actual_otus maps the new DOTUR numbers to an array of
# acutal OTU numbers. Things go in this array only if the sequence has
# a "true" OTU call from the database. It could happen that one of the
# DOTUR OTUs contains sequences from mulitplie true db OTUs, if this
# happens, we can give a confidence to it using otu_group_actual_otus

# Ryan.time_it("11 Calculate OTU assignment confidences") do
#   # TODO what will happen with sequences that have no match in the DB
#   otu_group_actual_otus = otu_group_actual_otus.map { |otu, group|
#     total = group.count
#     [otu,
#      group.each_with_object(Hash.new 0) { |elem, counts| counts[elem] += 1 }.
#        map { |otu, count| [otu, count / total] }]
#   }
# end

# Ryan.time_it("12 Final OTU calls written to #{final_otus}") do
#   File.open(final_otus, "w") do |f|
#     f.puts %w[header otu confidence].join "\t"
#     user_provided_headers.each do |header|
#       assert which_dotur_otu_group[header]
#       assert otu_group_actual_otus[which_dotur_otu_group[header]]

#       this_header_info =
#         otu_group_actual_otus[which_dotur_otu_group[header]]

#       this_header_info.last.each do |otu, percent|
#         f.puts [header, otu, percent].join "\t"
#       end
#     end
#   end
# end
