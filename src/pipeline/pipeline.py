class Pipeline(object):
    """ Abstract pipeline class

    Define methods for instancing and running a generic pipeline
    """

    # Constructor
    def __init__(self, *args, **kwargs):
        raise NotImplementedError

    # Run method
    def run(self, *args, **kwargs):
        raise NotImplementedError

    # Call method (wrap run method)
    def __call__(self, *args, **kwargs):
        # Call run method
        return self.run(*args, **kwargs)
