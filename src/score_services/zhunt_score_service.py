import numpy as np
import tempfile
import re
import os
from subprocess import check_call, DEVNULL

from .score_service import ScoreService


class ZHuntScoreService(ScoreService):
    def __init__(self):
        self.with_reversed = True

    def calculate_scores(self, seq):
        fd, tmp = tempfile.mkstemp()
        with open(tmp, 'w') as file:
            file.write(seq)

        check_call(['/home/pipeline/zrna_pipeline/score_services/zhunt3', '8', '6', '8', tmp], stdout=DEVNULL, stderr=DEVNULL)
        seq_zhunt = self._get_zhscore_list(tmp + '.Z-SCORE')

        assert len(seq_zhunt) == len(seq), f"For sequence {seq} \n Z-score has length {len(seq_zhunt)} vs sequence length {len(seq)}"
        os.close(fd)
        os.unlink(tmp)
        os.unlink(tmp + ".Z-SCORE")

        return seq_zhunt

    def _get_zhscore_list(self, filename):
        with open(filename) as f:
            lines = f.readlines()
        lines = lines[1:]
        return np.array([float(re.split('\s+', line)[5]) for line in lines])
