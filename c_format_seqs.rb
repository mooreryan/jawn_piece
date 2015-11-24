require "parse_fasta"

abort "USAGE: ruby c_format_seqs seqs.fa\n" if ARGV.empty?

seqs = []
FastaFile.open(ARGV.first).each_record do |head, seq|
  # TODO i think parse_fasta removes spaces...
  seqs << seq.upcase.tr("T", "U").gsub(/ +/, "")
end
total_seqs = seqs.count

seq_len = seqs.first.length

unless seqs.all? { |seq| seq.length == seq_len }
  abort "ERROR: Not all sequences are the same length!"
end

printf "%d %d\n", total_seqs, seq_len

seqs.each { |seq| puts seq }
