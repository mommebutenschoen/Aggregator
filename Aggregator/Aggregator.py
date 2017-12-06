from numpy import where,ones
from numpy.ma import getmaskarray, masked_where, median
from matplotlib.path import Path
from csv import writer,reader
from gzip import open as opengz
try:
    from itertools import izip as zip
except:
    pass
try:
    from itertools import irange as irange
except:
    pass

class AggregatorMapping:
    """"Base class for aggregation of high resolution data on coarser grids.
    Based on list of lists contaning the indices in a 1D representation of the
    high resolution data for each point of the low resolution data. The actual
    definition of mappings is achieved via the subclass "Aggregator", while this
    class may be used to read already defined mappings.

    Attributes:
        indices (list of lists of integers): Mappings of high resolution grid points to coarse grid points.
        size (integer): Number of coarse grid points.
    """

    def __init__(self,indices):

        """Defines mappings of points to polygon paths.

        Args:
            indices (sequence of list of integers): Mappings of high resolution grid points to coarse grid points.
        """

        #convert mapping sequence to list
        self.indices = [idx for idx in indices] #convert mapping sequence to list
        self.size=len(self.indices)

    def __call__(self,data,method=median,fv=1.e36,progress=False):

        """Aggregates high resolution data (in 1D) on coarse resolution 1D structure.

        Args:
            data(numpy.array or numpy.ma.array in 1D): High resolution data.
            method(function): Method used for aggregating reducing a sequence of values passed as input arguments to a single value.
            fv(same array dtype as data): Fill value to be used for coarse data where no high resolution data is available.
            progress(integer,): Interval in which to report mapping progress (message printing each "progress" polygons)

        Returns:
            coarse resolution data as 1D structure of same type as input data.
        """

        coarseData=fv*ones(self.size)
        for n,idn in enumerate(self.indices):
            if progress:
                if n%progress==0: print("Retrieved mappings for {} of {} polygons".format(n,self.size))
            if idn:
                aggregates=method(data[idn])
                coarseData[n]=where(getmaskarray(aggregates),fv,aggregates)
        if progress: print("Aggregation completed.")
        return coarseData

    def save_csv_gz(self,filename):
        """Save mapping to gzipped csv file."""
        with opengz(filename,'wt') as fid:
            csv=writer(fid)
            csv.writerows(self.indices)


class Aggregator(AggregatorMapping):

    """"Class for aggragation of high resolution data on coarser grids.
    The coarse grid needs to be defined by matplotlib.path.Path objects and
    high resolution data by an array of points. Initialisation of an Aggregator
    defines the mapping between the two grid and an object call aggregates a
    data field. Both, high and low resolution data need to be passed as 1D
    structures that are reshaped afterwards externally.

    Attributes:
        indices (list of lists of integers): Mappings of high resolution grid points to coarse grid polygons.
        size (integer): Number of coarse grid polygons.
        paths(sequence of matplotlib.path.Path objects): Polygons of coarse resolution grids.
        points(sequence of coordinate pairs): Cell centre points of high resolution grid pixels.
    """

    def __init__(self,paths,points,progress=False):

        """Defines mappings of points to polygon paths.

        Args:
            paths(sequence of matplotlib.path.Path objects): Polygons of coarse resolution grids.
            points(sequence of coordinate pairs): Cell centre points of high resolution grid pixels.
            progress(integer,): Interval in which to report mapping progress (message printing each "progress" polygons)
        """

        self.size=len(paths)
        idx=[ [] for n in range(self.size) ]
        if progress: n=0
        for path,idn in zip(paths,idx):
            if progress:
                if n%progress==0: print("Retrieved mappings for {} of {} polygons".format(n,self.size))
                n+=1
            if path.get_extents().size.any():
                inside=list(where(path.contains_points(points))[0])
                if len(inside)>0:
                    idn.extend(inside)
        self.indices = idx
        self.points = points
        self.paths = paths
        if progress: print("Mapping completed.")


def load_csv_gz(filename):
    """Load mapping from gzipped csv file.
    Args:
        filename(str): name of csv.gz file to read

    Returns:
        AggregatorMapping object with mapping defined in csv file."""

    with opengz(filename,'rt') as fid:
        csv=reader(fid)
        idx=[[int(n) for n in line] for line in csv]
    agg=AggregatorMapping(idx)

    return agg
