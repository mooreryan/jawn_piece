#!/usr/bin/env ruby

Signal.trap("PIPE", "EXIT")

methods = File.join(File.expand_path("~"), "lib", "ruby", "ryan.rb")
require_relative methods

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
      default: "test_files/database.fa")
  opt(:db_otu_info,
      "Accession to OTU info for the database",
      type: :string,
      default: "test_files/otus.txt")
  opt(:outdir,
      "Output directory",
      type: :string,
      default: "output")
end

alignment = Ryan.check_file(opts[:alignment], :alignment)
database = Ryan.check_file(opts[:database], :database)
Ryan.try_mkdir(opts[:outdir])

masked_seqs = File.join opts[:outdir], "masked_alignment.fa"
degapped_seqs = File.join opts[:outdir], "degapped_alignment.fa"
masked_for_phylip = File.join opts[:outdir], "masked_for_phylip.fa"
for_phylip_map = File.join opts[:outdir], "for_phylip_map.txt"
phylip_infile = File.join opts[:outdir], "DNAdist_params.txt"
this_dir = File.dirname(__FILE__)
phylip_outfile = File.join this_dir, "outfile"
masked_dist = File.join opts[:outdir], "masked_aln.dist"
dnadist = File.join this_dir, "bin", "dnadist"
dotur = File.join this_dir, "bin", "dotur"
dotur_outdir = File.join opts[:outdir], "dotur"
otu_calls = File.join dotur_outdir, "masked_aln.fn.list"
final_otus = File.join opts[:outdir], "final_otu_calls.txt"

Ryan.try_mkdir dotur_outdir

SILVA_ALN_LEN = 50000
OTU_LEVEL = 0.03


mask = ""
mask_posns = []
n = 0
database = {}
gap_posns = []
db_otu_info = {}
head_map = {} # used to prevent names from being too long for phylip
otu_calls_info = {}
user_provided_headers = Set.new
which_cutoff = -1
which_dotur_otu_group = {}
otu_group_actual_otus = {}


Ryan.time_it("0 Get database OTU metadata info") do
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

Ryan.time_it("1 Read the database for gaps and mask") do
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

Ryan.time_it("2 Update gap posns from user seqs & trim mask") do
  # first run through the sequences to adjust the mask
  FastaFile.open(opts[:alignment]).each_record do |head, seq|
    unless seq.length == SILVA_ALN_LEN
      warn "ERROR: #{head} in #{opts[:alignment]} is #{seq.length} bases, should be #{SILVA_ALN_LEN}"
    end

    if user_provided_headers.include? head
      abort "ERROR: #{head} duplicated in #{opts[:alignment]}"
    end

    user_provided_headers << head

    these_gap_posns = Set.new
    seq.each_char.with_index do |c, i|
      these_gap_posns << i if c == '-'
    end

    gap_posns << these_gap_posns

    first_seq_posn = seq.index /[^-]/
    last_seq_posn = seq.length - seq.reverse.index(/[^-]/) - 1

    while first_seq_posn > mask_posns.first
      # the sequence starts after the first mask posn
      # trim the mask
      mask_posns.shift
    end

    while last_seq_posn < mask_posns.last
      mask_posns.pop
    end
  end
end

Ryan.time_it("3 Apply mask to seqs & write") do
  # apply the mask, figure out gaps
  File.open(masked_seqs, "w") do |f|
    FastaFile.open(opts[:alignment]).each_record do |head, seq|
      masked_seq = mask_posns.map { |i| seq[i] }.join("").gsub(/U/, "T").gsub(/u/, "t")

      f.printf ">%s\n%s\n", head, masked_seq
    end

    database.each do |head, seq|
      masked_seq = mask_posns.map { |i| seq[i] }.join("")

      f.printf ">%s\n%s\n", head, masked_seq
    end
  end
end

Ryan.time_it("3.5 Write degapped alignment") do

  shared_gap_posns = gap_posns.reduce(:&)

  File.open(degapped_seqs, "w") do |f|
    FastaFile.open(opts[:alignment]).each_record do |head, seq|
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

Ryan.time_it("4 Make masked_for_phylip.fa") do

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
        phy_map.puts [new_name, head].join "\t"
        head_map[new_name] = head

        phy_aln.printf "%-10.10s%s\n", new_name, seq

        new_name = new_name.next
      end
    end
  end
end

Ryan.time_it("5 Write the DNAdist params file") do
  # write the phylip input file

  File.open(phylip_infile, 'w') do |f|
    f.puts masked_for_phylip
    f.puts "I"
    f.puts "Y"
  end
end

Ryan.time_it("6 DNAdist") do
  Ryan.run_it "#{dnadist} < #{phylip_infile}"
  Ryan.run_it "mv #{phylip_outfile} #{masked_dist}"
end

Ryan.time_it("7 DOTUR") do
  Ryan.run_it "#{dotur} #{masked_dist}"
  Ryan.run_it "mv #{masked_dist.sub(/dist$/, "fn")}.* #{dotur_outdir}"
end

Ryan.time_it("8 Read DOTUR OTU calls") do
  File.open(otu_calls).each_line do |line|
    unless line.start_with? "unique"
      cutoff, num, *rest = line.chomp.split "\t"

      otus = []

      rest.each do |otu_group|
        otus << otu_group.split(",").map do |header|
          assert head_map[header]

          head_map[header]
        end
      end

      otu_calls_info[cutoff.to_f] = otus
    end
  end
end


Ryan.time_it("9 Pick cutoff closest to #{OTU_LEVEL}") do
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


Ryan.time_it("10 Assign OTU numbers to OTU groups from the " +
             "#{opts[:db_otu_info]}") do
  assert otu_calls_info[which_cutoff]

  # otu_group is an array of headers contained in that OTU
  #
  # the idx is the number of that OTU from DOTUR, could be different
  # each run
  otu_calls_info[which_cutoff].each_with_index do |otu_group, idx|
    otu_group_actual_otus[idx] = []

    otu_group.each do |header|
      if which_dotur_otu_group.has_key? header
        abort "ERROR: #{header} repeated in the OTU groups"
      else
        which_dotur_otu_group[header] = idx
        # if the sequence is in the database
        if db_otu_info[header]
          otu_group_actual_otus[idx] << db_otu_info[header][:otu]
        end
      end
    end
  end
end

# otu_group_actual_otus maps the new DOTUR numbers to an array of
# acutal OTU numbers. Things go in this array only if the sequence has
# a "true" OTU call from the database. It could happen that one of the
# DOTUR OTUs contains sequences from mulitplie true db OTUs, if this
# happens, we can give a confidence to it using otu_group_actual_otus

Ryan.time_it("11 Calculate OTU assignment confidences") do
  # TODO what will happen with sequences that have no match in the DB
  otu_group_actual_otus = otu_group_actual_otus.map { |otu, group|
    total = group.count
    [otu,
     group.each_with_object(Hash.new 0) { |elem, counts| counts[elem] += 1 }.
       map { |otu, count| [otu, count / total] }]
  }
end

Ryan.time_it("12 Final OTU calls written to #{final_otus}") do
  File.open(final_otus, "w") do |f|
    f.puts %w[header otu confidence].join "\t"
    user_provided_headers.each do |header|
      assert which_dotur_otu_group[header]
      assert otu_group_actual_otus[which_dotur_otu_group[header]]

      this_header_info =
        otu_group_actual_otus[which_dotur_otu_group[header]]

      this_header_info.last.each do |otu, percent|
        f.puts [header, otu, percent].join "\t"
      end
    end
  end
end
