#!/bin/env -S uv run --script
# /// script
# requires-python = ">=3.14"
# dependencies = ["classy_blocks", "typing_extensions"]
# ///

# Generator script for the system/blockMeshDict

import classy_blocks as cb
import pathlib

x0 = 0
x1 = 2
x2 = 3
x3 = 5
x4 = 6

y0 = 0
y1 = 1
y2 = 2

# Blocks


def blockFor(l, r, b, t):
    return cb.Extrude(cb.Face([
        [l, b, 0],
        [r, b, 0],
        [r, t, 0],
        [l, t, 0]]), 1)


lb = blockFor(x0, x1, y0, y1)
lt = blockFor(x0, x1, y1, y2)

mt = blockFor(x1, x2, y1, y2)

rt1 = blockFor(x2, x3, y1, y2)
rb1 = blockFor(x2, x3, y0, y1)
rt2 = blockFor(x3, x4, y1, y2)
rb2 = blockFor(x3, x4, y0, y1)

# Patches


def patch(side: str, name: str, *items):
    for i in items:
        i.set_patch(side, name)


patch("left", "inlet", lb, lt)
patch("right", "outlet", rb2, rt2)

patch("back", "upperWall", lt, mt, rt1, rt2)
patch("front", "lowerWall", lb, rb1, rb2)

patch("right", "obstacle", lb)
patch("front", "obstacle", mt)
patch("left", "obstacle", rb1)

for b in [lb, lt, mt, rt1, rt2, rb1, rb2]:
    b.set_patch(["top", "bottom"], "frontAndBack")

# Cells


def chop(x, *items):
    for i in items:
        i.chop(0, count=x)
        i.chop(1, count=16)
        i.chop(2, count=1)


chop(20, lb, lt, rb1, rt1)
chop(10, mt)

rb2.chop(0, count=20, start_size=0.1)
rb2.chop(1, count=16)
rb2.chop(2, count=1)

rt2.chop(0, count=20, start_size=0.1)
rt2.chop(1, count=16)
rt2.chop(2, count=1)

# Mesh

m = cb.Mesh()
m.add(lb)
m.add(lt)
m.add(mt)
m.add(rb1)
m.add(rb2)
m.add(rt1)
m.add(rt2)

m.modify_patch("obstacle", "wall")
m.modify_patch("upperWall", "wall")
m.modify_patch("lowerWall", "wall")
m.modify_patch("frontAndBack", "empty")

m.write(pathlib.Path(__file__).parent / "system" / "blockMeshDict")
