# wavelet Earth Mover Distance
from .waveletEMB import WaveletEMD

# kmeans
from pyclustering.cluster.kmeans import kmeans# kmeans_observer, kmeans_visualizer
from pyclustering.cluster.center_initializer import kmeans_plusplus_initializer

# utils 
from pyclustering.utils.metric import distance_metric, type_metric
from pyclustering.utils import timedcall

# from pyclustering.cluster import cluster_visualizer_multidim

# set waveletEMD
wemd = distance_metric(type_metric.USER_DEFINED, func=WaveletEMD())

class KmeansOT:

    def __init__(self, n_clusters: int, metric = wemd, random_state: int = 2,):
        self.n_clusters = n_clusters
        self.metric = metric
        self.random_state = random_state # reproducibility
        self.kmeans = None # Kmeans instance
        self.start_centers = []

    def fit(self, data, tolerance=0.25, ccore=False):
        # initialize centers
        self.start_centers = kmeans_plusplus_initializer(
                                data, 
                                amount_centers = self.n_clusters, 
                                random_state=self.random_state
                                ).initialize()
        
        # kmeans
        self.kmeans = kmeans(data, self.start_centers, tolerance, ccore, metric=self.metric)
        self.kmeans.process() #train 
        # (ticks, _) = timedcall(self.kmeans.process) # train saving time execution

        # results
        self.clusters = self.kmeans.get_clusters()
        self.centers = self.kmeans.get_centers()

    def predict(self, data):
        return self.kmeans.predict(data)

