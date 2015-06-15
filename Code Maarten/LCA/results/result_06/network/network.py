#!/usr/bin/python
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

data= np.loadtxt("PairsPerShell.12.1.1", usecols=[0,1,2,4,5,6,9] )

G=nx.Graph()
for list in data:
  shell1=100*int(list[0])+10*int(list[1])+int(list[2])
  shell2=100*int(list[3])+10*int(list[4])+int(list[5])
  G.add_edge( shell1, shell2, weight=list[6] )
  

print G.nodes()
print G.edges(data=True)
pos= nx.spring_layout(G)

nx.draw(G,pos)
edge_labels=dict([((u,v,),d['weight']) for u,v,d in G.edges(data=True)])
nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels)



plt.show()


