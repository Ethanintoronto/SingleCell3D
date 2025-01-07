import pyvista as pv
import numpy as np
import os

def write_polyhedron_to_vtk(solid_name, file_path):
    """
    Write the vertices, edges, and polygons of a Platonic solid to a VTK file.

    Parameters:
        solid_name (str): Name of the Platonic solid (e.g., 'tetrahedron', 'cube').
        file_path (str): Path to the output VTK file.
    """
    # Create the Platonic solid
    try:
        solid = pv.PlatonicSolid(solid_name)
    except ValueError:
        print(f"Error: Invalid solid name '{solid_name}'. Choose from 'tetrahedron', 'cube', 'octahedron', 'dodecahedron', 'icosahedron'.")
        return
    
    # scale to volume of 1:
    solid.points = solid.points*(1/solid.volume)**(1/3)



    # Extract the vertices and faces
    vertices = solid.points
    FV = solid.faces[0]
    faces = solid.faces.reshape(-1, FV+1)[:, 1:]  # Remove face sizes for VTK compatibility
    print(solid.volume)
    print(solid.area)

    # Write to a VTK file
    with open(file_path, 'w') as vtk_file:
        # Write header
        vtk_file.write("# vtk DataFile Version 3.0\n")
        vtk_file.write(f"{solid_name.capitalize()} Solid\n")
        vtk_file.write("ASCII\n")
        vtk_file.write("DATASET POLYDATA\n")

        # Write vertices
        vtk_file.write(f"POINTS {len(vertices)} float\n")
        for vertex in vertices:
            vtk_file.write(f"{' '.join(map(str, vertex))}\n")

        # Write polygons (faces)
        vtk_file.write(f"POLYGONS {len(faces)} {len(faces) * FV+1}\n")
        for face in faces:
            vtk_file.write(f"{FV} {' '.join(map(str, face))}\n")

# Usage
# Create an output directory if it doesn't exist
output_dir = "vtk_in"
os.makedirs(output_dir, exist_ok=True)

# List of Platonic solids
solids = ["tetrahedron", "cube", "octahedron", "dodecahedron", "icosahedron"]

# Write each solid to a VTK file
for solid_name in solids:
    output_path = os.path.join(output_dir, f"{solid_name}.vtk")
    write_polyhedron_to_vtk(solid_name, output_path)
