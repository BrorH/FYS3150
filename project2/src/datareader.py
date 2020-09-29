import numpy as np


class SolvedSystem:
    """HOW DO YOU USE THIS PLEASE WRITE"""

    def __init__(self, datablock):
        self.datablock = datablock
        self.parse()

    def parse(self):
        self.name = self.datablock[0].strip()
        n, self.pmax, self.eps, self.transformations = tuple(
            [float(a) for a in self.datablock[1].split(",")]
        )
        self.n = int(n)
        self.diag = [float(a) for a in self.datablock[2][:-2].split(",")]
        self.eigvals = np.zeros(self.n)
        self.eigvecs = np.zeros((self.n, self.n))
        for i, line in enumerate(self.datablock[3:]):
            splitted = line.split(",")

            self.eigvals[i] = float(splitted[0])
            self.eigvecs[:, i] = np.asarray([float(a) for a in splitted[1:-1]])


def read_data(filename="data.dat"):
    with open(f"data/{filename}") as file:
        lines = file.readlines()
    instances = {}
    blocks = []
    last = 0
    for idx, line in enumerate(lines):
        if "*" in line:
            blocks.append(lines[last:idx])
            last = idx + 1
    for block in blocks:
        instance = SolvedSystem(block)
        instances[instance.name] = instance
    return instances
