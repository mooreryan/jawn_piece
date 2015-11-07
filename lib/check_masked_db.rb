methods = File.join(File.expand_path("~"), "lib", "ruby", "ryan.rb")
require_relative methods
require_relative "../const"

lineno = 0
actual_num_seqs =0
num_seqs = 0
File.open(Const::MASKED_DATABASE).each_line do |line|
  if lineno.zero?
    num_seqs, seq_len = line.chomp.split("\t").map(&:to_f)
  else
    seq = line.chomp
    actual_num_seqs += 1

    unless seq.length == seq_len
      abort "ERROR: #{num_seqs} is #{seq.length} not #{seq_len}"
    end
  end
end

unless actual_num_seqs == num_seqs
  abort "ERROR: got #{actual_num_seqs} should have been #{num_seqs}"
end

puts "No issues in #{Const::MASKED_DATABASE}"
