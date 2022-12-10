class Alignment:
    def __init__(self, seqA, seqB, score):
        self.seqA = seqA
        self.seqB = seqB
        self.score = score
        self.start = 0
        self.end = max(len(self.seqA), len(self.seqB))