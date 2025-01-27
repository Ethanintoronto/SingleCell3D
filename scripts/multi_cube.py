def generate_vtk_for_stacked_cubes(filename):
    """
    Generate a VTK file for 4 stacked cubes where adjacent cubes share vertices and faces.

    Parameters:
        filename (str): The output VTK file name.
    """
    # Hardcoded vertices for 4 stacked cubes with shared vertices
    vertices = [
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],  # Bottom face of the first cube
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],  # Top face of the first cube, shared with bottom face of the second cube
        [0, 0, 2], [1, 0, 2], [1, 1, 2], [0, 1, 2],  # Top face of the second cube, shared with bottom face of the third cube
        [0, 0, 3], [1, 0, 3], [1, 1, 3], [0, 1, 3],  # Top face of the third cube, shared with bottom face of the fourth cube
        [0, 0, 4], [1, 0, 4], [1, 1, 4], [0, 1, 4]   # Top face of the fourth cube
    ]

    # Hardcoded faces for all cubes
    faces = [
        [0, 1, 5, 4],  # Front face of the first cube
        [1, 2, 6, 5],  # Right face of the first cube
        [2, 3, 7, 6],  # Back face of the first cube
        [3, 0, 4, 7],  # Left face of the first cube
        [0, 1, 2, 3],  # Bottom face of the first cube
        [4, 5, 6, 7],  # Top face of the first cube (shared with bottom of second cube)
        [4, 5, 9, 8],  # Front face of the second cube
        [5, 6, 10, 9],  # Right face of the second cube
        [6, 7, 11, 10],  # Back face of the second cube
        [7, 4, 8, 11],  # Left face of the second cube
        [8, 9, 10, 11],  # Top face of the second cube (shared with bottom of third cube)
        [8, 9, 13, 12],  # Front face of the third cube
        [9, 10, 14, 13],  # Right face of the third cube
        [10, 11, 15, 14],  # Back face of the third cube
        [11, 8, 12, 15],  # Left face of the third cube
        [12, 13, 14, 15],  # Top face of the third cube (shared with bottom of fourth cube)
        [12, 13, 17, 16],  # Front face of the fourth cube
        [13, 14, 18, 17],  # Right face of the fourth cube
        [14, 15, 19, 18],  # Back face of the fourth cube
        [15, 12, 16, 19],  # Left face of the fourth cube
        [16, 17, 18, 19]   # Top face of the fourth cube
    ]

    # Define cubes by referencing face indices
    cubes = [
        [0, 1, 2, 3, 4, 5],  # First cube
        [5, 6, 7, 8, 9, 10],  # Second cube
        [10, 11, 12, 13, 14, 15],  # Third cube
        [15, 16, 17, 18, 19, 20]   # Fourth cube
    ]

    # Labels for each cube
    cube_labels = [1, 2, 3, 4]

    with open(filename, 'w') as vtk_file:
        # Header for VTK file
        vtk_file.write("# vtk DataFile Version 3.0\n")
        vtk_file.write("Stacked Cubes\n")
        vtk_file.write("ASCII\n")
        vtk_file.write("DATASET POLYDATA\n")

        # Write points (vertices)
        vtk_file.write(f"POINTS {len(vertices)} float\n")
        for vertex in vertices:
            vtk_file.write(f"{vertex[0]} {vertex[1]} {vertex[2]}\n")

        # Write polygons (faces)
        total_faces = len(faces)
        total_nums = sum(len(face) + 1 for face in faces)  # Add 1 for the vertex count in each face
        vtk_file.write(f"POLYGONS {total_faces} {total_nums}\n")
        for face in faces:
            vtk_file.write(f"{len(face)} {' '.join(map(str, face))}\n")

    with open(filename[:-4]+"_helper.txt",'w') as vtk_helper:
        # Write cell data for grouping faces into cubes
        total_cubes = len(cubes)
        total_cube_nums = sum(len(cube) + 1 for cube in cubes)  # Add 1 for the vertex count in each face
        vtk_helper.write(f"CELLS {total_cubes} {total_cube_nums}\n")
        for cube in cubes:
            vtk_helper.write(f"{len(cube)} {' '.join(map(str, cube))}\n")

# Usage
output_filename = ".\\vtk_in\\stacked_cubes.vtk"
generate_vtk_for_stacked_cubes(output_filename)
print(f"VTK file '{output_filename}' has been generated.")
