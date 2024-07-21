import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


# Load the CSV file
data = pd.read_csv('A_Finite_Element_project_ClassDependencies.csv')

# Generate a directed graph from the CSV file
# the starting point are the 'From Class' and the end point are the 'To Class'
# the 'References' column is used to set the width of the edges
# the edges are to be inserted in the exact order as they appear in the CSV file
G = nx.from_pandas_edgelist(data, 'From Class', 'To Class', ['References'])

#Print the edges
print(G.edges())

# Create a layout for the graph
pos = nx.fruchterman_reingold_layout(G, k = 1, iterations = 1000)

# Set the figure background color to black
fig = plt.figure(figsize=(20, 10), facecolor='black')

# Draw the nodes with white labels and rectangular shape
nx.draw_networkx_nodes(G, pos, node_color='black' , node_shape='s', node_size = 500, edgecolors = (255/255, 140/255, 0/255))
nx.draw_networkx_labels(G, pos, font_color='white', font_size=12, font_weight='bold', font_family='serif')

# Draw the edges with gray color and widths based on the number of references
# Access the edges trough the graph and the data=True parameter and set the width of the edges
# to be stored in a vector based on the number of references for the current edge
edge_widths = [0.6 * G.edges[edge]['References'] for edge in G.edges]
print(edge_widths)
# Make the edges directed and use arrows
nx.draw_networkx_edges(G, pos, width=edge_widths, edge_color='gray', alpha=0.7, arrows=True)

# Remove the axis labels
plt.axis('off')

plt.savefig('class_dependencies.jpg', dpi = 200, bbox_inches='tight')
plt.show()