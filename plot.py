import os
import matplotlib.pyplot as plt
import imageio

def line_intersection(edge, x=None, y=None):
    """Calculate the intersection point of a line and a vertical or horizontal line."""
    (x1, y1), (x2, y2) = edge
    if x is not None:
        y = y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    elif y is not None:
        x = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
    return (x, y)

# Read edges
edges = []
with open('edges.txt', 'r') as file:
    for line in file:
        start, end = line.strip().split()
        start = tuple(map(float, start.split(',')))
        end = tuple(map(float, end.split(',')))
        edges.append((start, end))

# Read points
points = []
with open('input.txt', 'r') as file:
    lines = file.readlines()
    for i in range(1, len(lines), 2):
        x = float(lines[i].strip())
        y = float(lines[i+1].strip())
        points.append((x, y))

# Define the bounding box
min_x = min(point[0] for point in points) - 10
max_x = max(point[0] for point in points) + 10
min_y = min(point[1] for point in points)  - 10
max_y = max(point[1] for point in points) + 10

# Truncate the edges according to the bounding box
truncated_edges = []
for edge in edges:
    start, end = edge
    if start[0] < min_x:
        start = line_intersection(edge, x=min_x)
    elif start[0] > max_x:
        start = line_intersection(edge, x=max_x)
    if start[1] < min_y:
        start = line_intersection(edge, y=min_y)
    elif start[1] > max_y:
        start = line_intersection(edge, y=max_y)
    if end[0] < min_x:
        end = line_intersection(edge, x=min_x)
    elif end[0] > max_x:
        end = line_intersection(edge, x=max_x)
    if end[1] < min_y:
        end = line_intersection(edge, y=min_y)
    elif end[1] > max_y:
        end = line_intersection(edge, y=max_y)
    truncated_edges.append((start, end))

# Create a list to store the filenames of the images
filenames = []

# Plot edges one by one and save each plot as an image
for i, edge in enumerate(truncated_edges):
    for e in truncated_edges[:i+1]:
        plt.plot(*zip(*e), 'k-')
    for point in points:
        plt.plot(*point, 'ro', markersize=3)
    plt.plot([min_x, min_x, max_x, max_x, min_x], [min_y, max_y, max_y, min_y, min_y], 'b--')  # Draw bounding box
    plt.xlim(min_x, max_x)
    plt.ylim(min_y, max_y)
    plt.axis('equal')
    filename = f'plot_{i}.png'
    plt.savefig(filename)
    filenames.append(filename)
    plt.clf()

# Plot all edges and save the final plot
for edge in truncated_edges:
    plt.plot(*zip(*edge), 'k-')
for point in points:
    plt.plot(*point, 'ro', markersize=3)
plt.plot([min_x, min_x, max_x, max_x, min_x], [min_y, max_y, max_y, min_y, min_y], 'b--')  # Draw bounding box
plt.xlim(min_x, max_x)
plt.ylim(min_y, max_y)
plt.axis('equal')
plt.savefig('final_plot.png')

# Create a gif from the images
with imageio.get_writer('plot.gif', mode='I') as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)

# Remove the image files
for filename in filenames:
    os.remove(filename)