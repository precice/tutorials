#!/bin/python3

import precice
import numpy as np
import math
import sys
from mpi4py import MPI

import argparse


def split(num):
    for a in range(math.isqrt(num), 0, -1):
        if num % a == 0:
            return a, num // a


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("participant", choices=["A", "B"])
    parser.add_argument("--config", "-c", default="../precice-config.xml")
    parser.add_argument("--no-remesh", dest="remesh", action="store_false")
    parser.add_argument("--nx", type=int, default=512)
    parser.add_argument("--ny", type=int, default=512)
    args = parser.parse_args()

    participant_name = args.participant
    remote_name = "A" if participant_name == "B" else "B"

    # x is partitioned per rank and doesn't change
    nx = args.nx
    ny = args.ny
    x = 0.0, 1.0
    y = 0.0, 1.0

    # y grows over time
    newNodesPerEvent = 2
    eventFrequency = 3  # time windows
    dz = 0.1

    # Handle partitioning
    world = MPI.COMM_WORLD
    size: int = world.size
    rank: int = world.rank
    xr, yr = split(size)

    assert nx % xr == 0, f"Cannot split {nx=} by {xr=} ranks"
    assert ny % yr == 0, f"Cannot split {ny=} by {yr=} ranks"

    pnx = nx // xr
    pny = ny // yr

    px = rank % xr
    py = rank // xr
    print(f"Rank {rank}/{size} has partition ({px}, {py})/({xr}, {yr})")
    if rank == 0:
        print(
            f"Each of {size} partitions has node size {pnx}x{pny} = {
                pnx *
                pny} for a total of {
                nx *
                ny} nodes on the base"
        )

    def getMesh(nz):
        basex = np.linspace(0, 1, nx)[px * pnx: (px + 1) * pnx]
        basey = np.linspace(0, 1, ny)[py * pny: (py + 1) * pny]
        z = np.array(range(nz)) * dz
        return np.stack(np.meshgrid(basex, basey, z, indexing="ij"), axis=-1).reshape(
            -1, 3
        )

    def requiresEvent(tw):
        return tw % eventFrequency == 0

    assert not requiresEvent(eventFrequency - 1)
    assert requiresEvent(eventFrequency)
    assert not requiresEvent(eventFrequency + 1)

    def eventsAt(tw):
        # First event block at tw=0, second at eventFrequency
        return 1 + math.floor(tw / eventFrequency)

    assert eventsAt(0) == 1
    assert eventsAt(eventFrequency - 1) == 1
    assert eventsAt(eventFrequency) == 2
    assert eventsAt(eventFrequency + 1) == 2

    def getMeshAtTimeWindow(tw):
        znodes = eventsAt(tw) * newNodesPerEvent
        return getMesh(znodes)

    def dataAtTimeWindow(tw):
        znodes = eventsAt(tw) * newNodesPerEvent
        total = pnx * pny * znodes
        return np.full(total, tw)

    participant = precice.Participant(participant_name, args.config, rank, size)

    mesh_name = participant_name + "-Mesh"
    read_data_name = "Data-" + remote_name
    write_data_name = "Data-" + participant_name

    coords = getMeshAtTimeWindow(0)
    vertex_ids = participant.set_mesh_vertices(mesh_name, coords)
    participant.initialize()

    tw = 1
    while participant.is_coupling_ongoing():
        dt = participant.get_max_time_step_size()

        data = participant.read_data(mesh_name, read_data_name, vertex_ids, dt)
        if rank == 0:
            print(f"{participant_name}: data: {data[0]} * {len(data)}")

        if args.remesh and requiresEvent(tw):
            oldCount = len(coords)
            coords = getMeshAtTimeWindow(tw)
            if rank == 0:
                print(
                    f"{participant_name}: Event grows local mesh from {oldCount} to {
                        len(coords)} and global mesh from {
                        oldCount *
                        size} to {
                        len(coords) *
                        size}"
                )
            participant.reset_mesh(mesh_name)
            vertex_ids = participant.set_mesh_vertices(mesh_name, coords)

        participant.write_data(
            mesh_name, write_data_name, vertex_ids, dataAtTimeWindow(tw)
        )

        participant.advance(dt)
        tw += 1


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e)
        sys.exit(1)
