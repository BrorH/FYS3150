
class SolvedSystem:
    """HOW DO YOU USE THIS PLEASE WRITE"""
    def __init__(self, datablock):
        self.datablock = datablock
        self.parse()
    def parse(self):
        self.n , self.pmax, self.eps, self.counts = tuple([ float(a) for a in self.datablock[0].split(",")])
        self.diag = [float(a) for a in self.datablock[1][:-2].split(",")]
        self.eigvals = []
        self.eigvecs = []
        for line in self.datablock[2:]:
            splitted = line.split(",")
            self.eigvals.append(float(splitted[0]))
            self.eigvecs.append([float(a) for a in splitted[1:-1]])

        
def read_data(filename = "data.dat"):
    with open(filename) as file:
        lines = file.readlines()
    instances = []
    blocks = []
    last = 0
    for idx,line in enumerate(lines):
        if "*" in line:
            blocks.append(lines[last:idx])
            last = idx+1
    for block in blocks:
        instances.append(SolvedSystem(block))
    return instances
