import io
import argparse
import numpy as np
import pyvista as pv
import pyproj
from meshfem3d_constants import R_EARTH


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert GMT multisegment file (e.g. pscoast -M) to VTK",
        epilog="Generate shoreline file using: gmt pscoast -W1 -Di -R20/50/30/45 -M > shoreline.txt"
        " and convert to VTK using: python gmt_multiseg2vtk.py shoreline.txt shoreline.vtk",
    )
    parser.add_argument("input_file", help="filename of GMT multisegment file")
    parser.add_argument(
        "output_vtk",
        help="filename of the output VTK file",
    )
    parser.add_argument(
        "--sep",
        default=">",
        help="separator at the beginning of the comment line of each segment",
    )
    return parser.parse_args()


def process(input_file, output_file, delimiter=">"):

    gps2ecef = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:4978")

    # Read the entire file content as a single string
    with open(input_file, "r", encoding="utf-8") as f:
        file_content = f.read()

    # Split the content by the special character/delimiter
    # Note: The delimiter itself is not included in the resulting list
    sections = file_content.split(delimiter)

    # Filter out any empty strings that might result from the split
    # sections = [section.strip() for section in sections if section.strip()]

    # Write each section to a new file
    vtk_lines = []
    for section in sections:
        lines = section.split("\n")[1:] # remove the first line with the separator
        if not lines:
            continue

        data_string = "\n".join(lines)  
        f = io.StringIO(data_string)
        data = np.loadtxt(f)
        if data.size == 0:
            continue

        data = np.reshape(data, (-1, 2)) # in case only one line is given
        lons = data[:, 0]
        lats = data[:, 1]

        xx, yy, zz = gps2ecef.transform(lats, lons, np.zeros_like(lats))
        points = np.stack([xx, yy, zz], axis=1) / R_EARTH
        line = pv.lines_from_points(points)
        vtk_lines.append(line)

    mesh = pv.merge(vtk_lines)
    mesh.save(output_file)


def main():
    """Main execution function."""
    args = parse_arguments()
    print(args)

    process(args.input_file, args.output_vtk, delimiter=args.sep)


if __name__ == "__main__":
    main()
