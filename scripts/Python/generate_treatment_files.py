#
# Author: Stefan Pszczolkowski P. <mszspp@nottingham.ac.uk>
#
import transforms3d as tf
import nibabel as nib
import numpy as np
import argparse
import os
import subprocess


# Use Freesurfer's mri_info to obtain one of the image's transformation matrices (defined in transformation type)
def get_image_transformation_matrix(file, transformation_type):
    output = subprocess.check_output(['mri_info', file, '--{}'.format(transformation_type)])
    num_output = [[float(x) for x in line.split()] for line in output.decode('utf8').strip().split('\n') if len(line)>0]

    return np.matrix(num_output)


# Estimate the rotation an translation that brings the point set p1 into point set p2 using singular value decomposition
def rigid_transform_3D(p1, p2):
    assert len(p1) == len(p2)

    num_points = p1.shape[0]

    # Compute centroids
    centroid_p1 = np.mean(p1, axis=0)
    centroid_p2 = np.mean(p2, axis=0)

    # Compute centered data
    p1_centered = p1 - np.tile(centroid_p1, (num_points, 1))
    p2_centered = p2 - np.tile(centroid_p2, (num_points, 1))

    # Run SVD
    U, S, Vt = np.linalg.svd(np.transpose(p1_centered) * p2_centered)

    # Get rotation matrix
    R = Vt.T * U.T

    # Take care of special reflection case
    if np.linalg.det(R) < 0:
        R[:, 2] *= -1

    # Compute translation
    t = -R * centroid_p1.T + centroid_p2.T

    return R, t


# Compute intersection point where a line L1 defined by line_point1 -> line_point2 crosses a line L2 that passes through
# ref_point
def compute_intersection_point(line_point1, line_point2, ref_point):
    dir_vec = (line_point1 - line_point2) / np.linalg.norm(line_point1 - line_point2)
    aux_vec = ref_point - line_point2

    d = np.dot(dir_vec, np.transpose(aux_vec))

    inter_point = line_point2 + d * dir_vec

    return inter_point


# Transform vertex points using the provided rotation and translation
def transform_point_set(vertex_points, R, t):
    num_vertices = vertex_points.shape[0]

    transf_vertex_points = np.zeros((num_vertices, 3))

    for vertex_id in np.arange(num_vertices):
        vertex_coord = np.matrix(vertex_points[vertex_id], dtype=np.float)
        transf_vertex_points[vertex_id, :] = np.transpose(R * np.transpose(vertex_coord) + t)

    return transf_vertex_points


# Compute the id and coordinates of the vertex closest to point_native
def compute_closest_vertex(vertex_points, point_native):
    best_vertex_id = None
    best_vertex_coord = None
    best_distance = np.Inf

    num_vertices = vertex_points.shape[0]

    for vertex_id in np.arange(num_vertices):
        vertex_coord = np.matrix(vertex_points[vertex_id, 0:3], dtype=np.float)
        dist = np.linalg.norm(vertex_coord - point_native)

        if dist < best_distance:
            best_distance = dist
            best_vertex_id = vertex_id
            best_vertex_coord = vertex_coord

    return best_vertex_id, best_vertex_coord


# Compute the normal vector to the surface in vertex_id
def compute_positive_vertex_normal(vertex_points, vertex_id, triangles):
    normal_vector = np.zeros((1, 3))
    triangle_count = 0

    num_triangles = triangles.shape[0]

    vertex_point = vertex_points[vertex_id]

    for triangle_id in np.arange(num_triangles):
        if triangles[triangle_id, 0] == vertex_id:
            normal_vector += np.cross(vertex_points[triangles[triangle_id, 2]] - vertex_point, vertex_points[triangles[triangle_id, 1]] - vertex_point)
            triangle_count += 1
        elif triangles[triangle_id, 1] == vertex_id:
            normal_vector += np.cross(vertex_points[triangles[triangle_id, 0]] - vertex_point, vertex_points[triangles[triangle_id, 2]] - vertex_point)
            triangle_count += 1
        elif triangles[triangle_id, 2] == vertex_id:
            normal_vector += np.cross(vertex_points[triangles[triangle_id, 1]] - vertex_point, vertex_points[triangles[triangle_id, 0]] - vertex_point)
            triangle_count += 1

    normal_vector = normal_vector / triangle_count

    return np.abs(normal_vector)


def compute_mrt_transform():
    mid_coord_mrt = np.matrix([256, 182, 256], dtype=np.float) / 2

    # Nasion-ear-points in MRT coordinates
    nasion_coord_mrt = np.matrix([220, 92, 109], dtype=np.float) - mid_coord_mrt
    left_coord_mrt = np.matrix([124, 171, 88], dtype=np.float) - mid_coord_mrt
    right_coord_mrt = np.matrix([121, 12, 86], dtype=np.float) - mid_coord_mrt

    # Create vectors for NLR system
    yvec = (left_coord_mrt - right_coord_mrt) / np.linalg.norm(left_coord_mrt - right_coord_mrt)
    origin = right_coord_mrt + np.dot(yvec, np.transpose(nasion_coord_mrt - right_coord_mrt)) * yvec
    xvec = (nasion_coord_mrt - origin) / np.linalg.norm(nasion_coord_mrt - origin)
    zvec = np.cross(xvec, yvec)

    R = np.matrix(np.column_stack((np.transpose(xvec), np.transpose(yvec), np.transpose(zvec))), dtype=np.float)
    t = np.transpose(origin)

    return R, t


# Write StimGuide files
def write_stimguide_files(file_path, subject_id, nasion_ear_coord, mrt_coord, mrt_quat):
    target_name = 'Treatment_{}'.format(subject_id)

    file_name_xml = os.path.join(file_path, subject_id + '.xml')
    with open(file_name_xml, 'w') as sg_xml_file:
        sg_xml_file.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        sg_xml_file.write('<Targets>\n')
        sg_xml_file.write('\t<Target>\n')
        sg_xml_file.write('\t\t<Name>{}</Name>\n'.format(target_name))
        sg_xml_file.write('\t</Target>\n')
        sg_xml_file.write('</Targets>\n')

    file_name_asc = os.path.join(file_path, subject_id + '.asc')
    with open(file_name_asc, 'w') as sg_asc_file:
        sg_asc_file.write('NumberPositions=\t1\n')
        sg_asc_file.write('NumberPositionOrientations_TMSCoilTarget=\t1\n')
        sg_asc_file.write('UnitPosition\tmm\n')
        sg_asc_file.write('HSPTransformed\tfalse\n')
        sg_asc_file.write('Positions\n')
        sg_asc_file.write('{}:\t{:3.1f} {:3.1f} {:3.1f}\n'.format(target_name, nasion_ear_coord.item(0), nasion_ear_coord.item(1), nasion_ear_coord.item(2)))
        sg_asc_file.write('PositionOrientations_TMSCoilTarget\n')
        sg_asc_file.write('{}:\t{:3.8f} {:3.8f} {:3.8f} {:3.8f} {:3.8f} {:3.8f} {:3.8f}\n'.format(target_name, mrt_coord.item(0), mrt_coord.item(1), mrt_coord.item(2), mrt_quat[1], mrt_quat[2], mrt_quat[3], mrt_quat[0]))
        sg_asc_file.write('Labels\n')
        sg_asc_file.write('{}\n'.format(target_name))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create StimGuide coordinate file for the specified subject', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('subject_id', help='Subject ID')
    parser.add_argument('input_coord_file', help='Input text file with coordinates in native T1 image space')
    parser.add_argument('output_file_dir', help='Directory where the output StimGuide file with coordinate information should be written. If it does not exist, it will be created')
    parser.add_argument('--subjects-dir', dest='subjects_dir', default=None, help='Specify alternative FreeSurfer subjects directory')

    args = parser.parse_args()

    subject_id = args.subject_id
    input_coord_file = args.input_coord_file
    output_file_dir = args.output_file_dir
    subjects_dir = args.subjects_dir

    # Default SUBJECTS_DIR is the environ variable SUBJECTS_DIR or, if this is unavailable, a directory called 'subjects' in the user's home directory
    if subjects_dir is None:
        subjects_dir = os.environ['SUBJECTS_DIR']
        if subjects_dir is None:
            subjects_dir = os.path.join(os.environ['HOME'], 'subjects')

    mgz_t1_file = os.path.join(subjects_dir, subject_id, 'mri', 'T1.mgz')
    gifti_scalp_file = os.path.join(subjects_dir, subject_id, 'surf', 'lh.seghead.surf.gii')

    if not os.path.exists(output_file_dir):
        os.makedirs(output_file_dir)

    # Load landmarks text file
    native_coords = np.loadtxt(input_coord_file)

    # Point set in image native coordinates
    target_coord_native = np.matrix(native_coords[0], dtype=np.float)
    nasion_coord_native = np.matrix(native_coords[1], dtype=np.float)
    left_coord_native = np.matrix(native_coords[2], dtype=np.float)
    right_coord_native = np.matrix(native_coords[3], dtype=np.float)

    gifti_surface = nib.load(gifti_scalp_file)

    tkr_to_native_full = get_image_transformation_matrix(mgz_t1_file, 'tkr2scanner')
    R_tkr_to_native = tkr_to_native_full[0:3, 0:3]
    t_tkr_to_native = tkr_to_native_full[0:3, 3]

    vertex_coords_native = transform_point_set(gifti_surface.darrays[0].data, R_tkr_to_native, t_tkr_to_native)

    # Compute vertices id and scalp coordinate
    target_scalp_vertex_id, target_scalp_coord_native = compute_closest_vertex(vertex_coords_native, target_coord_native)
    nasion_scalp_vertex_id, nasion_scalp_coord_native = compute_closest_vertex(vertex_coords_native, nasion_coord_native)
    left_scalp_vertex_id, left_scalp_coord_native = compute_closest_vertex(vertex_coords_native, left_coord_native)
    right_scalp_vertex_id, right_scalp_coord_native = compute_closest_vertex(vertex_coords_native, right_coord_native)

    print("Target scalp vertex ID is " + str(target_scalp_vertex_id))
    print("Nasion scalp vertex ID is " + str(nasion_scalp_vertex_id))
    print("Left scalp vertex ID is " + str(left_scalp_vertex_id))
    print("Right scalp vertex ID is " + str(right_scalp_vertex_id))
    print("================================================================================================")
    print("Target scalp coordinate in native T1 space is " + str(target_scalp_coord_native))
    print("Nasion scalp coordinate in native T1 space is " + str(nasion_scalp_coord_native))
    print("Left scalp coordinate in native T1 space is " + str(left_scalp_coord_native))
    print("Right scalp coordinate in native T1 space is " + str(right_scalp_coord_native))

    origin_coord_native = compute_intersection_point(left_coord_native, right_coord_native, nasion_coord_native)

    # Distances from origin to landmarks in mm
    no_dist = np.linalg.norm(nasion_coord_native - origin_coord_native)
    lo_dist = np.linalg.norm(left_coord_native - origin_coord_native)
    ro_dist = np.linalg.norm(right_coord_native - origin_coord_native)

    # Point set in NLR coordinates
    nasion_coord_nasion_ear = np.matrix([no_dist, 0, 0], dtype=np.float)
    left_coord_nasion_ear = np.matrix([0, lo_dist, 0], dtype=np.float)
    right_coord_nasion_ear = np.matrix([0, -ro_dist, 0], dtype=np.float)
    origin_coord_nasion_ear = np.matrix([0, 0, 0], dtype=np.float)

    # Stack point sets
    coord_set_native = np.vstack((nasion_coord_native, left_coord_native, right_coord_native, origin_coord_native))
    coord_set_nasion_ear = np.vstack((nasion_coord_nasion_ear, left_coord_nasion_ear, right_coord_nasion_ear, origin_coord_nasion_ear))

    # Compute transformation between point set in native coordinates to point set in NLR coordinates
    R_native_to_nasion_ear, t_native_to_nasion_ear = rigid_transform_3D(coord_set_native, coord_set_nasion_ear)

    # Apply transformation to vertices
    vertex_coords_nasion_ear = transform_point_set(vertex_coords_native, R_native_to_nasion_ear, t_native_to_nasion_ear)

    # Apply transformation to target point
    target_scalp_coord_nasion_ear = np.transpose(R_native_to_nasion_ear * np.transpose(target_scalp_coord_native) + t_native_to_nasion_ear)

    # Compute surface normal at the target scalp coordinate
    target_surface_normal = compute_positive_vertex_normal(vertex_coords_nasion_ear, target_scalp_vertex_id, gifti_surface.darrays[1].data)

    # If (u,v,w) is the surface normal direction at the target in NLR coordinates and (a,b,c) or (-a,-b,-c) is an
    # orthogonal vector to (u,v,w) defining the direction towards the front of the coil (y coil axis), we have the
    # following three conditions:
    #
    # 1.- u*a + v*b + w*c = 0  (orthogonality condition)
    #
    # 2.- dot((a,b,0), (1,0,0)) = cos(45) = sqrt(2)/2  ->  a/sqrt(a^2+b^2) = sqrt(2)/2  ->  a^2/(a^2+b^2) = 1/2
    #
    # 3.- dot((a,b,0), (0,-1,0)) = sin(45) = sqrt(2)/2  ->  -b/sqrt(a^2+b^2) = sqrt(2)/2  ->  b^2/(a^2+b^2) = 1/2
    #
    # Conditions 2 and 3 mean that the projection of (a,b,c) into the XY plane in NLR coordinates (i.e. (a,b,0)) forms
    # a 45 degree angle with the NLR X axis.
    #
    # The solutions to this system of equations are:
    #
    # a = -w / (u - v)
    #
    # b = w / (u - v)
    #
    # c = 1
    #
    # Hence, aux_nasion_ear = b = -a
    aux_nasion_ear = target_surface_normal.item(2) / (target_surface_normal.item(0) - target_surface_normal.item(1))

    if aux_nasion_ear < 0:
        # Y coil axis is point (a,b,c) in NLR coordinates
        y_coil_axis_nasion_ear = np.matrix([-aux_nasion_ear, aux_nasion_ear, 1], dtype=np.float)
    else:
        # Y coil axis is point (-a,-b,-c) in NLR coordinates
        y_coil_axis_nasion_ear = np.matrix([aux_nasion_ear, -aux_nasion_ear, -1], dtype=np.float)

    # Z coil axis is surface normal in NLR coordinates
    z_coil_axis_nasion_ear = target_surface_normal

    # X coil axis is cross product between Y and Z
    x_coil_axis_nasion_ear = np.cross(y_coil_axis_nasion_ear, z_coil_axis_nasion_ear)

    # Corresponding X, Y and Z axes in coil coordinates
    x_axis_coil = np.matrix([np.linalg.norm(x_coil_axis_nasion_ear), 0, 0], dtype=np.float)
    y_axis_coil = np.matrix([0, np.linalg.norm(y_coil_axis_nasion_ear), 0], dtype=np.float)
    z_axis_coil = np.matrix([0, 0, np.linalg.norm(z_coil_axis_nasion_ear)], dtype=np.float)

    # Origin in coil coordinates is (0,0,0)
    origin_coord_coil = np.matrix([0, 0, 0], dtype=np.float)

    # Stack point sets
    coord_set_coil_axes_nasion_ear = np.vstack((x_coil_axis_nasion_ear, y_coil_axis_nasion_ear, z_coil_axis_nasion_ear, origin_coord_nasion_ear))
    coord_set_coil_axes_coil = np.vstack((x_axis_coil, y_axis_coil, z_axis_coil, origin_coord_coil))

    # Compute transformation between point set in coil coordinates to point set in NLR coordinates
    # t_coil_to_nasion_ear should always be very close to (0,0,0)
    R_coil_to_nasion_ear, t_coil_to_nasion_ear = rigid_transform_3D(coord_set_coil_axes_coil, coord_set_coil_axes_nasion_ear)

    # Compute quaternions of the coil to NLR rotation matrix
    q_coil_to_nasion_ear = tf.quaternions.mat2quat(R_coil_to_nasion_ear)

    # Compute transformation between nasion ear coordinates and MRT coordinates
    R_nasion_ear_to_mrt, t_nasion_ear_to_mrt = compute_mrt_transform()

    # Apply transformation to target point
    scalp_coord_mrt = np.transpose(R_nasion_ear_to_mrt * np.transpose(target_scalp_coord_nasion_ear) + t_nasion_ear_to_mrt)

    # Compute quaternions of the NLR to MRT (StimGuide's coordinate system) rotation matrix
    q_nasion_ear_to_mrt = tf.quaternions.mat2quat(R_nasion_ear_to_mrt)

    # Multiply quaternions
    q_coil_to_mrt = tf.quaternions.qmult(q_nasion_ear_to_mrt, q_coil_to_nasion_ear)

    # Write files
    write_stimguide_files(output_file_dir, subject_id, target_scalp_coord_nasion_ear, scalp_coord_mrt, q_coil_to_mrt)
