#!/usr/bin/env ruby

MAX_SEQS_FOR_COLOR = 50

def round_down num
  rnum = num.round(-1)
  if rnum > num
    rnum -= 10
  end

  rnum
end

def steps last
  (10..last).step(10)
end

def make_step_str total_len
  last = round_down total_len

  the_steps = steps last

  str = "~" * total_len

  start_str = "Sequence"
  start_str.each_char.with_index do |char, idx|
    str[idx] = char
  end

  the_steps.each do |pos| # one based
    str_pos = pos.to_s
    len = str_pos.length

    (pos-1 .. pos-1+len).each_with_index do |n, idx|
      if idx.zero?
        str[n] = "|"
      else
        str[n] = str_pos[idx-1]
      end
    end
  end

  str
end

def round_if_needed str
  if str.match(/^[0-9]+\.[0-9]+$/)
    str.to_f.round(5)
  else
    str
  end
end

def check_f fname
  unless File.exists? fname
    abort "ERROR: File does not exist -- #{fname}"
  end
end

def file_to_table fname
  html_string = "<table>"
  File.open(fname).each_line do |line|
    if line.start_with? "query"
      tag = "th"
    else
      tag = "td"
    end

    str = "<tr>"
    line.chomp.split("\t").each do |elem|
      elem = round_if_needed elem

      str << "<%s>%s</%s>" % [tag, elem, tag]
    end
    str << "</tr>"

    html_string << str
  end
  html_string << "</table>"

  html_string
end

Signal.trap("PIPE", "EXIT")

methods = File.join(File.expand_path("~"), "lib", "ruby", "ryan.rb")
require_relative methods
require_relative "const"

Ryan.req *%w[parse_fasta]

opts = Trollop.options do
  banner <<-EOS

  Pass in the outdir from ZetaHunter and get a nice html report.

  Options:
  EOS

  opt(:directory, "The ZetaHunter outdir", type: :string)
end

directory = Ryan.check_file(opts[:directory], :directory)

degapped_mask_f =
  File.join opts[:directory], "degapped_mask.fa"
degapped_alignment_f =
  File.join opts[:directory], "degapped_alignment.fa"
otu_calls_f =
  File.join opts[:directory], "final_otu_calls.txt"
closed_ref_otu_info_f =
  File.join opts[:directory], "closed_ref_otu_info.txt"
final_calls_f =
  File.join opts[:directory], "final_otu_calls.txt"
possible_non_zetas_f =
  File.join opts[:directory], "possible_non_zetas.txt"
seanie_f =
  File.join File.dirname(__FILE__), "assets", "seanie.html"

check_f degapped_mask_f
check_f degapped_alignment_f
check_f otu_calls_f
check_f closed_ref_otu_info_f
check_f final_calls_f
check_f possible_non_zetas_f
check_f seanie_f

degapped_seqs = FastaFile.open(degapped_alignment_f).to_hash
otu_info = {}

html_f = File.join opts[:directory], "run_info.html"

possible_non_zetas = []
File.open(possible_non_zetas_f).each_line do |line|
  possible_non_zetas << line.chomp
end

File.open(otu_calls_f).each_line do |line|
  unless line.start_with? "query"
    query, otu = line.chomp.split "\t"

    if otu_info.has_key? query
      abort "ERROR: #{query} is repeated in #{otu_calls_f}"
    end

    otu_info[query] = otu
  end
end

# get mask posns
mask = ""
FastaFile.open(degapped_mask_f).each_record do |head, seq|
  mask = seq.gsub("-", " ")
end

# # add database seqs to otu_info
# FastaFile.open(Const::DATABASE).each_record do |head, seq|
#   unless head == "Mask"
#     otu_info[head] = "zdb"
#   end
# end

# # get colors for the degapped thing
all_seqs = []
otu_info.each_with_index do |(seq_name, otu), idx|
  unless degapped_seqs.has_key? seq_name
    abort "ERROR: #{degapped_alignment_f} doesn't contain #{seq_name}, but #{otu_calls_f} or #{Const::DATABASE} does"
  end

  seq = degapped_seqs[seq_name]
  all_seqs << seq.chars
end

num_seqs = all_seqs.count.to_f
if num_seqs <= MAX_SEQS_FOR_COLOR
  positional_color = all_seqs.transpose.map.with_index do |arr, which_pos|
    # count number of each sequence
    base_counts = arr.each_with_object(Hash.new(0)) do |base, hash|
      hash[base] += 1
    end

    base_freqs = base_counts.map do |base, count|
      [base, count / num_seqs]
    end

    max_freq = base_freqs.sort_by { |_, freq| freq }.reverse.first.last

    (1 - max_freq).round(2)
  end
end

# positional_color.each_with_index do |n, i|
#   puts [i+1, n].join "\t"
# end


# make degapped string for html file
degapped_seqs_html_string = ""
if num_seqs <= MAX_SEQS_FOR_COLOR
  color = "#000000"
  bg_color = "#ffffff"
else
  color = "#ffffff"
  bg_color = "#11BA0D"
end
# TODO change to giant table and then use <colgroup> tag.. might be
# faster for big alignments
otu_info.sort_by { |name, otu| otu }.each_with_index do |(seq_name, otu), idx|
  unless degapped_seqs.has_key? seq_name
    abort "ERROR: #{degapped_alignment_f} doesn't contain #{seq_name}, but #{otu_calls_f} does"
  end

  if (idx % 20).zero?
    step_str = make_step_str degapped_seqs[seq_name].length
    inner_str = "%-20.20s | %-20.20s | %s\n"
    str = "<span style=\"color: #{color}; background-color: " +
          "#{bg_color};\">%-20.20s   %-20.20s   %s\n#{inner_str}</span>"
    step_html_str = str % ["", "", mask,
                           "OTU",
                           "Sequence name",
                           step_str]
  else
    step_html_str = ""
  end

  seq = degapped_seqs[seq_name]

  seq_html_str = ""
  if num_seqs <= MAX_SEQS_FOR_COLOR
    seq.each_char.with_index do |char, idx|
      color = positional_color[idx]
      seq_html_str << %Q|<span style="background-color: rgba(17,186,138,#{color}">#{char}</span>|
    end
  else
    seq_html_str << seq
  end

  this_string = "%-20.20s | %-20.20s | %s\n" % [otu, seq_name, seq_html_str]
  degapped_seqs_html_string << (step_html_str + this_string)
end

# all_alignments_html =
#   degapped_seqs_html_string.each_slice(10).map do |strs|
#   first_part = %Q{<div class="container" id="degapped-alignment">\n}
#   last_part  = "\n</div>\n"

#   first_part + strs.join + last_part
# end.join "\n"

# make closed ref info html file
closed_ref = file_to_table closed_ref_otu_info_f

# get the final OTU calls
final_calls = file_to_table final_calls_f

seanie = File.open("assets/seanie.html").read

html = %Q{<!DOCTYPE html>
<html>

  <head>
    <style type="text/css">

      body {
        font-family: sans-serif;
      }

      h2 {
        font-family: sans-serif;
      }

      .container {
        overflow: scroll;
        width: 1600px;
        white-space: pre;
        font-family: monospace;
      }

      table, th, td {
        border: 2px solid black;
        border-collapse: collapse;
      }

      th, td {
        padding: 10px;
      }
    </style>
  </head>

  <body>
    <h1>ZetaHunter v0.0.1</h1>
    <h2>Degapped alignment</h2>
    <div class="container" id="degapped-alignment">
#{degapped_seqs_html_string}
    </div>

    <h2>Possible non-Zetas</h2>
#{possible_non_zetas.map { |n| "<p>#{n}" }.join("\n")}

    <h2>Closed reference OTU info</h2>
    <div id="closed-ref-otu-info">
#{closed_ref}
    </div>

    <h2>Final OTU calls</h2>
    <div id="final-calls">
#{final_calls}
    </div>

    <h2>Have a nice day!</h2>
    <div id="thank-you">
#{seanie}
    </div>
  </body>
}

File.open(html_f, "w") do |f|
  f.puts html
end

`open #{html_f}`
