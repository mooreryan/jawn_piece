require "bio"
require "set"
require "fail_fast"

include FailFast::Assertions

infile = "/Users/moorer/projects/ZetaHunter/output/small/trees/masked_for_phylip.tre"
phymap = "/Users/moorer/projects/ZetaHunter/output/small/for_phylip_map.txt"

treeio = Bio::FlatFile.open(Bio::Newick, infile)

newick = treeio.next_entry
TREE = newick.tree

bad = "f"
m1 = "o"
m2 = "j"
m3 = "d"

ALL = ("a" .. "dw").to_a
db_start = ALL.index("u")
outgroups_start = ALL.index("da")

QUERIES   = ALL[0 .. db_start-1]
DATABASE  = ALL[db_start .. outgroups_start - 1]
OUTGROUPS = ALL[outgroups_start .. ALL.length-1]

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
  puts name
end
