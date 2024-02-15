import os
import sys
import string
import random
import pydot
import xml.etree.ElementTree as ET
import cairosvg
import argparse
import shutil

from eugene.solution import BaseSolution


ALLELE_0 = "blue"
ALLELE_1 = "yellow"
FILE_SEED = "".join(random.choices(string.hexdigits, k=10))


def visualise_cs(d: BaseSolution, output_file: str):
    graph = pydot.Dot("G")

    nodes = {}
    for i in range(d.n_plants):
        if d.tree_type[i] != "Null":
            node = _make_node(i, d.tree_data[i])
            nodes[i] = node
            graph.add_node(node)

    for i in range(d.n_plants):
        if d.tree_type[i] == "Node":
            graph.add_edge(pydot.Edge(nodes[d.tree_left[i] - 1], nodes[i]))
            graph.add_edge(pydot.Edge(nodes[d.tree_right[i] - 1], nodes[i]))

    graph.write_png(output_file)
    # graph.write_svg("/tmp/output.svg")


def _make_node(i, array):
    generate_svg_from_array(i, array)
    node = pydot.Node(
        i, shape="none", label="", image=_node_file(i), width="1", height="1"
    )
    return node


def _node_file(i):
    return f"/tmp/crossing_schedule_{FILE_SEED}/node_{i}.png"


def generate_svg_from_array(node_id, arr):
    # Define the size of each square and spacing
    square_size = 20
    spacing = 0

    # Calculate the width and height of the SVG based on the array dimensions
    width = len(arr[0]) * (square_size + spacing)
    height = len(arr) * (square_size + spacing)

    # Create the root SVG element
    svg = ET.Element(
        "svg", width=str(width), height=str(height), xmlns="http://www.w3.org/2000/svg",
    )

    for i, row in enumerate(arr):
        for j, value in enumerate(row):
            x = j * (square_size + spacing)
            y = i * (square_size + spacing)

            color = ALLELE_1 if value == 1 else ALLELE_0
            ET.SubElement(
                svg,
                "rect",
                x=str(x),
                y=str(y),
                width=str(square_size),
                height=str(square_size),
                fill=color,
                stroke="black",
            )

    # Convert the SVG element to a string and return
    image_bytes = ET.tostring(svg, encoding="unicode").encode("utf-8")
    res = cairosvg.svg2png(image_bytes, write_to=_node_file(node_id))
    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser("csvis")
    parser.add_argument(
        "--output", "-o", default="/tmp/output.png", help="file to write to"
    )
    args = parser.parse_args()
    output_file = args.output
    os.mkdir(f"/tmp/crossing_schedule_{FILE_SEED}")
    d = eval(sys.stdin.read())
    print(d)
    visualise_cs(d, output_file)
    shutil.rmtree(f"/tmp/crossing_schedule_{FILE_SEED}")
