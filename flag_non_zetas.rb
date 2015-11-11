require "bio"
require "set"
require "fail_fast"

include FailFast::Assertions

infile = ARGV[0]
phymap = ARGV[1]

treeio = Bio::FlatFile.open(Bio::Newick, infile)

newick = treeio.next_entry
TREE = newick.tree

bad = "f"
m1 = "o"
m2 = "j"
m3 = "d"

name_map = {}
queries = []
database = []
outgroups = []
all = []
File.open(phymap).each_line do |line|
  new_name, old_name, type = line.chomp.split "\t"

  name_map[new_name] = old_name

  if type == "query"
    queries << new_name
  elsif type == "outgroup"
    outgroups << new_name
  elsif type == "database"
    database << new_name
  else
    abort "ERROR: malformed type, was #{type}"
  end

  all << new_name
end

ALL = all

QUERIES   = queries
DATABASE  = database
OUTGROUPS = outgroups

# TODO if node name is nil <- does that work?
def leaf? node
  TREE.children(node).empty?
end

# if the parent of the node has all descenants that are out groups,
# will flag as bad
def bad_node_by_siblings? node
  parent = TREE.parent node

  TREE.descendents(parent).select do |nd|
    leaf? nd
  end.map do |nd|
    if nd.name == node.name
      nil
    elsif OUTGROUPS.include?(nd.name)
      true
    else
      false
    end
  end.compact.any?
end

# get the top ten closest nodes to each node and see if any of the
# outgroups are in there, if so, then flag it
all_dists = []
total = (ALL.count * 2).to_f / 2
ALL.each_with_index do |name, idx|
  $stderr.printf "%.2f%%\r", idx / total * 100
  node = TREE.get_node_by_name name

  (ALL - [name]).each do |n2|
    node2 = TREE.get_node_by_name n2
    dist = 0
    TREE.each_edge_in_path(node, node2) do |s, t, e|
      dist += e.distance
    end

    all_dists << [name, n2, dist]
  end
end
$stderr.print "\n"

outgroups_in_closest_ten = Set.new
all_dists.group_by { |name, name2, dist| name }.each do |key, arr|
  arr.sort_by { |n1, n2, d| d }.take(10).each do |n1, n2, d|
    # puts [n1, OUTGROUPS.include?(n2), d].join "\t"
    outgroups_in_closest_ten << n1 if OUTGROUPS.include?(n2)
  end
end

outgroups_in_closest_ten.each do |name|
  unless OUTGROUPS.include? name
    assert name_map[name]
    puts name_map[name]
  end
end
