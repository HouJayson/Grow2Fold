# =============================================================================
# Copyright (c) 2026 Jixin Hou et al.
# All rights reserved.
#
# This code is provided as part of the research software accompanying:
# [Grow2Fold: mapping heterogeneous developmental growth to human brain folding]
#
# Use, modification, and redistribution are permitted under the terms of the
# license provided in the LICENSE file of this repository.
#
# Repository: [https://github.com/BioDMX-UGA/Grow2Fold]
# =============================================================================

import sys
sys.path.append(r"path-to-abaPython\Lib\site-packages")

from odbAccess import openOdb
import numpy as np
import os
from scipy.io import savemat


# ============================================================
# User settings
# ============================================================
odb_path = r"path-to-odb file"
file_path = r"path-to-excel"

instance_name = "BRAIN-1"

user_id = "*"

# Extract every N frames:
frame_interval = 1

step1_name = "Step-1"
step2_name = "Step-2"

# Step-1 has frames 0--8, so n_frames_step1 = 9
# Step-2 has frames 0--17, so n_frames_step2 = 18
n_frames_step1 = 9
n_frames_step2 = 18

# Ignore Step-2 frame 0 because it duplicates/end-starts from Step-1
skip_step2_frame0 = True


# ============================================================
# Helper functions
# ============================================================

def build_frame_plan(
    n_frames_step1,
    n_frames_step2,
    frame_interval=1,
    skip_step2_frame0=True
):

    if frame_interval < 1:
        raise ValueError("frame_interval must be >= 1")

    frame_plan = []

    start_step2 = 1 if skip_step2_frame0 else 0

    # Total number of exported candidate frames after skipping duplicated Step-2 frame 0
    total_available_frames = n_frames_step1 + (n_frames_step2 - start_step2)

    for global_frame_id in range(0, total_available_frames, frame_interval):

        if global_frame_id < n_frames_step1:
            step_name = step1_name
            local_frame_id = global_frame_id
        else:
            step_name = step2_name
            local_frame_id = start_step2 + (global_frame_id - n_frames_step1)

        frame_plan.append((step_name, local_frame_id, global_frame_id))

    return frame_plan


def extract_point_set_label_from_odb(odb, instance_name, point_set_name):
    assembly = odb.rootAssembly
    instance = assembly.instances[instance_name]
    point_set = instance.nodeSets[point_set_name]

    label_sequence = []

    for i in range(len(point_set.nodes)):
        label_sequence.append(point_set.nodes[i].label)

    return sorted(label_sequence)


def extract_point_set_coordinates_global_from_frame(
    odb,
    instance_name,
    point_set_name,
    step_name,
    frame_id
):
    assembly = odb.rootAssembly
    instance = assembly.instances[instance_name]
    point_set = instance.nodeSets[point_set_name]

    frame = odb.steps[step_name].frames[frame_id]
    coords = frame.fieldOutputs["COORD"]

    mySetCoord = coords.getSubset(region=point_set)

    coordinates = {}

    for v in mySetCoord.values:
        node_label = v.nodeLabel
        node_coords_vec = [v.data[0], v.data[1], v.data[2]]
        coordinates[node_label] = node_coords_vec

    return coordinates


def extract_point_set_displacement_global_from_frame(
    odb,
    instance_name,
    point_set_name,
    step_name,
    frame_id
):
    assembly = odb.rootAssembly
    instance = assembly.instances[instance_name]
    point_set = instance.nodeSets[point_set_name]

    frame = odb.steps[step_name].frames[frame_id]
    disp = frame.fieldOutputs["U"]

    mySetDisp = disp.getSubset(region=point_set)

    displacement = {}

    for v in mySetDisp.values:
        node_label = v.nodeLabel
        node_disp_vec = [v.data[0], v.data[1], v.data[2]]
        displacement[node_label] = np.linalg.norm(node_disp_vec)

    return displacement


# ============================================================
# Main
# ============================================================

def main():

    frame_plan = build_frame_plan(
        n_frames_step1=n_frames_step1,
        n_frames_step2=n_frames_step2,
        frame_interval=frame_interval,
        skip_step2_frame0=skip_step2_frame0
    )

    num_frames = len(frame_plan)

    print("User ID:", user_id)
    print("Frame interval:", frame_interval)
    print("Total output frames:", num_frames)

    odb = openOdb(odb_path)

    try:
        for region in ["GRAY", "WHITE"]:

            point_set_name = "SET-NODE-{}-OUTERMOST".format(region)

            mat_name = "{}.{}.par_Clustering.Hex_thk14_skull_2steps_stress.mat".format(
                user_id,
                region.lower()
            )

            print("\nProcessing region:", region)
            print("Node set:", point_set_name)

            # Extract node labels once
            extraction_sequence = extract_point_set_label_from_odb(
                odb,
                instance_name,
                point_set_name
            )

            num_nodes = len(extraction_sequence)

            print("Number of nodes:", num_nodes)

            # Preallocate: num_nodes × 4 × num_frames
            # 4 = x, y, z, displacement magnitude
            data_all = np.zeros((num_nodes, 4, num_frames), dtype=np.float32)

            # Save step/frame information for checking later in MATLAB
            output_frame_ids = np.array(
                [item[2] for item in frame_plan],
                dtype=np.int32
            )

            step_ids = np.zeros(num_frames, dtype=np.int32)
            local_frame_ids = np.zeros(num_frames, dtype=np.int32)

            for out_frame_id, item in enumerate(frame_plan):

                step_name, local_frame_id, global_frame_id = item

                if step_name == step1_name:
                    step_ids[out_frame_id] = 1
                elif step_name == step2_name:
                    step_ids[out_frame_id] = 2

                local_frame_ids[out_frame_id] = local_frame_id

                print("{} output frame {}: total frame {} -> {} frame {}".format(
                    region,
                    out_frame_id,
                    global_frame_id,
                    step_name,
                    local_frame_id
                ))

                coords = extract_point_set_coordinates_global_from_frame(
                    odb=odb,
                    instance_name=instance_name,
                    point_set_name=point_set_name,
                    step_name=step_name,
                    frame_id=local_frame_id
                )

                disp = extract_point_set_displacement_global_from_frame(
                    odb=odb,
                    instance_name=instance_name,
                    point_set_name=point_set_name,
                    step_name=step_name,
                    frame_id=local_frame_id
                )

                for i, node_label in enumerate(extraction_sequence):

                    x, y, z = coords[node_label]
                    d = disp[node_label]

                    # Handle NaN/Inf safely
                    if not np.isfinite(x):
                        x = 0.0
                    if not np.isfinite(y):
                        y = 0.0
                    if not np.isfinite(z):
                        z = 0.0
                    if not np.isfinite(d):
                        d = 0.0

                    data_all[i, :, out_frame_id] = [x, y, z, d]

            save_path = os.path.join(file_path, mat_name)

            savemat(save_path, {
                "data": data_all,
                "node_labels": np.array(extraction_sequence, dtype=np.int32),
                "output_frame_ids": output_frame_ids,   # now stores 0, 5, 10, 15, 20, 25
                "step_ids": step_ids,
                "local_frame_ids": local_frame_ids,
                "user_id": user_id,
                "frame_interval": frame_interval,
            })

            print("Saved:", save_path)

    finally:
        odb.close()

    print("\nMAT files saved.")


if __name__ == "__main__":
    main()