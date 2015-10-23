#!/usr/bin/env ruby

def check_f fname
  unless File.exists? fname
    abort "ERROR: File does not exist -- #{f}"
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

Ryan.req *%w[parse_fasta]

opts = Trollop.options do
  banner <<-EOS

  Info

  Options:
  EOS

  opt(:directory, "The ZetaHunter outdir", type: :string)
end

directory = Ryan.check_file(opts[:directory], :directory)

degapped_alignment_f =
  File.join opts[:directory], "degapped_alignment.fa"
otu_calls_f =
  File.join opts[:directory], "final_otu_calls.txt"
closed_ref_otu_info_f =
  File.join opts[:directory], "closed_ref_otu_info.txt"
final_calls_f =
  File.join opts[:directory], "final_otu_calls.txt"
seanie_f =
  File.join File.dirname(__FILE__), "assets", "seanie.html"

check_f degapped_alignment_f
check_f otu_calls_f
check_f closed_ref_otu_info_f
check_f final_calls_f
check_f seanie_f

degapped_seqs = FastaFile.open(degapped_alignment_f).to_hash
otu_info = {}

html_f = File.join opts[:directory], "run_info.html"

File.open(otu_calls_f).each_line do |line|
  unless line.start_with? "query"
    query, otu = line.chomp.split "\t"

    if otu_info.has_key? query
      abort "ERROR: #{query} is repeated in #{otu_calls_f}"
    end

    otu_info[query] = otu
  end
end

# make degapped string for html file
degapped_seqs_html_string = "%-20.20s | %-20.20s | %s\n" % ["OTU", "Sequence name", "Sequence"]
otu_info.each do |seq_name, otu|
  unless degapped_seqs.has_key? seq_name
    abort "ERROR: #{degapped_alignment_f} doesn't contain #{seq_name}, but #{otu_calls_f} does"
  end

  seq = degapped_seqs[seq_name]
  this_string = "%-20.20s | %-20.20s | %s\n" % [otu, seq_name, seq]
  degapped_seqs_html_string << this_string
end

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
    <h2>Degapped alignment</h2>
    <div class="container" id="degapped-alignment">
#{degapped_seqs_html_string}
    </div>

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
