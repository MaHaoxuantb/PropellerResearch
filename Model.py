import numpy as np
import trimesh
from shapely.geometry import Polygon


def c(theta, h=1.0, a=2.0):
    return h * np.cos(a * theta)


def make_profile(h=1.0, a=2.0, n=500):
    theta = np.linspace(
        -np.pi / (2 * a),
         np.pi / (2 * a),
        n
    )

    r = c(theta, h, a)

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    points = np.column_stack([x, y])

    # Close the region back to the origin
    points = np.vstack([
        [0.0, 0.0],
        points,
        [0.0, 0.0]
    ])

    return Polygon(points)


profile = make_profile(h=1, a=2)

mesh = trimesh.creation.extrude_polygon(
    profile,
    height=0.2
)

print(mesh.is_watertight)
print(mesh.volume)

mesh.export("shape.stl")
