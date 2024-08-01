class Bulge():
    def __init__(self, unpaired, opposite):
        self.unpaired = unpaired
        self.opposite = opposite

    def len_unpaired(self):
        return self.unpaired[1] - self.unpaired[0] + 1
