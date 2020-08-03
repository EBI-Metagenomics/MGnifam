from src.pipeline.batch import Batch


class Build(Batch):
    """ Make alignments against MGnifam

    This pipeline works similarly to Batch, from which inherits all methods
    (except for constructor and run which are overwritten), but searches
    against MGnifam dataset instead of UniPRot, then retrieves significant
    matches and uses them to provide wider alignments with respect to
    SEED ones, in terms of sequences number.
    """

    # Constructor
    def __init__(self, *args, **kwargs):
        raise NotImplementedError

    # Run
    def run(self, *args, **kwargs):
        raise NotImplementedError
