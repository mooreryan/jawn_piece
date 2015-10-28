#!/usr/bin/env ruby

def get_char_freqs str
  str_len = str.length.to_f

  str.chars.group_by(&:itself).values.map do |arr|
    arr.count / str_len
  end
end

# Shannon entropy (see http://rosettacode.org/wiki/Entropy#Ruby)
def entropy str
  char_freqs = get_char_freqs str

  char_freqs.reduce(0) do |entropy, char_freq|
    entropy - char_freq * Math.log2(char_freq)
  end
end

def mean arr
  arr.reduce(:+) / arr.count.to_f
end

def windows arr, size, step
  # I checked the arr.length - size for off by one errors
  (0..arr.length-size).step(step).map do |start|
    [(start*2 + size)/2.0, mean(arr[start, size])]
  end
end

Signal.trap("PIPE", "EXIT")

methods = File.join(File.expand_path("~"), "lib", "ruby", "ryan.rb")
require_relative methods

Ryan.req *%w[parse_fasta]

opts = Trollop.options do
  banner <<-EOS

  Calculate positional entropy in the ZetaMask.

  The first sequence of the --database file must be Mask.

  Options:
  EOS

  opt(:database, "Input database (SILVA aligned)", type: :string,
      default: File.join("test_files", "database.fa"))
  opt(:outdir, "Output directory", type: :string,
      default: File.join("output", "entropy"))
end

# database = Ryan.check_file(opts[:database], :database)
Ryan.try_mkdir(opts[:outdir])

mask_posns = [] # 0-based posns
masked_sequences = []
posn_entropy = []
windowed_entropy = []

entropy_f = File.join opts[:outdir], "entropy.txt"
entropy_pdf = File.join opts[:outdir], "entropy.pdf"
windowed_entropy_f = File.join opts[:outdir], "windowed_entropy.txt"
rscript = ""
rscript_f = File.join opts[:outdir], "entropy.r"

WINDOW = 25
STEP = 10
SILVA_LEN = 50000

Ryan.time_it("Create the masked sequences") do
  n = 0
  FastaFile.open(opts[:database]).each_record do |head, seq|
    if seq.length != SILVA_LEN
      abort "ERROR: #{head} is not #{SILVA_LEN}" +
            "Are your line endings correct?"
    end

    if n.zero? && head != "Mask"
      abort "ERROR: The first sequence should be Mask, was #{head}"
    end

    if head == "Mask"
      seq.each_char.with_index do |char, idx|
        if char == "*"
          mask_posns << idx
        end
      end
    else
      masked_sequence = mask_posns.map { |idx| seq[idx] }.join

      masked_sequences << masked_sequence
    end

    n += 1
  end
end

desc = "Calculate positional entropy of masked sequences"
Ryan.time_it(desc) do
  masked_sequences.map(&:chars).
    transpose.
    each_with_index do |chars, idx|
    pos = idx + 1
    posn_entropy << entropy(chars.join)
  end
end

desc = "Get entropy windows"
Ryan.time_it(desc) do
  #posn_entropy.each_cons(20)

  windowed_entropy = windows posn_entropy, WINDOW, STEP
end

Ryan.time_it("Plot entropy") do
  File.open(entropy_f, "w") do |f|
    posn_entropy.each_with_index do |ent, idx|
      f.puts [idx, ent].join "\t"
    end
  end

  File.open(windowed_entropy_f, "w") do |f|
    windowed_entropy.each do |x, y|
      f.puts [x, y].join "\t"
    end
  end

  File.open(rscript_f, "w") do |f|
    rscript = %Q{
pdf("#{entropy_pdf}", height = 8, width = 12)
dat <- read.table("#{entropy_f}", header = F, sep = "\t")
plot(dat$V1, dat$V2, col = "gray48", pch = 20, xlab = "Position",
     ylab = "Entropy", type = "l")
dat <- read.table("#{windowed_entropy_f}", header = F, sep = "\t")
points(dat$V1, dat$V2, col = "red", type = "l", lwd = 3)
invisible(dev.off())
}
    f.puts rscript
  end

  Ryan.run_it "Rscript #{rscript_f}"
  Ryan.run_it "open #{entropy_pdf}"

  zero_posns = posn_entropy.select { |x| x.zero? }.count
  total_posns = posn_entropy.count.to_f
  percent_zero = zero_posns / total_posns * 100
  $stderr.printf "\n\nZero entropy positions in the mask: %.2f%\n",
                 percent_zero
end
