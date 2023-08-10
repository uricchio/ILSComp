from dendropy import Tree

tree = "[&R] ((53402-Pv:0.030704,(53428-Pp:0.008057,53495-Pm:0.009725):0.010835):0.002275,(53477-Ps:0.015629,53468-Pc:0.013995):0.004497,((53473-Pc:0.357422,53409-Pb:0.149353):0.076321,53503-Pc:0.052807):0.030323);"
myTree = Tree.get(data=tree,schema="Newick")

outgroup_node = myTree.find_node_with_taxon_label("53473-Pc")
myTree.reroot_at_midpoint()  #to_outgroup_position(outgroup_node, update_bipartitions=False)

print(myTree.as_string(schema="Newick"))

