
import numpy as np
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
import geopandas as gpd
from shapely.geometry import Point


# CLUSTERING POINTS THAT ARE LESS THAN 15 MIN (1300 M) APART

# code from https://stackoverflow.com/questions/53075481/how-do-i-cluster-a-list-of-geographic-points-by-distance


def distance_euclid(origin, destination):
    """
    Calculates the Euclidean distance between two points in 2D space.

    Parameters:
    ----------
    origin : tuple or list
        Coordinates of the origin point as (x, y).
    destination : tuple or list
        Coordinates of the destination point as (x, y).

    Returns:
    -------
    float
        The Euclidean distance between the origin and destination points.
    """
    d = np.sqrt((destination[0]-origin[0])**2 + (destination[1]-origin[1])**2)
    return(d)
    
def create_clusters(number_of_clusters,points):
    """
    Creates clusters from a set of points using the K-Means algorithm.

    Parameters:
    ----------
    number_of_clusters : int
        The desired number of clusters to create.
    points : np.ndarray
        A 2D numpy array of shape (n_points, 2) representing the (x, y) coordinates of the points.

    Returns:
    -------
    tuple (np.ndarray, np.ndarray)
        clusters : np.ndarray
            A 2D array where each row contains the (x, y) coordinates of a point and its cluster label.
        kmeans.cluster_centers_ : np.ndarray
            A 2D array of shape (n_clusters, 2) containing the coordinates of the cluster centers.
    """
    kmeans = KMeans(n_clusters=number_of_clusters, n_init = 'auto', random_state=0).fit(points)
    l_array = np.array([[label] for label in kmeans.labels_])
    clusters = np.append(points,l_array, axis=1)
    return(clusters, kmeans.cluster_centers_)

def validate_solution(max_dist, clusters, centers):
    """
    Validates whether a clustering solution satisfies a maximum distance constraint.

    Parameters:
    ----------
    max_dist : float
        The maximum allowable distance between any two points within a cluster.
    clusters : np.ndarray
        A 2D array where each row contains the (x, y) coordinates of a point and its cluster label.
    centers : np.ndarray
        A 2D array of shape (n_clusters, 2) representing the coordinates of the cluster centers.

    Returns:
    -------
    tuple (bool, np.ndarray, np.ndarray)
        A tuple containing:
        - bool: True if the solution is valid, False otherwise.
        - centers: The input array of cluster centers.
        - clusters: The input array of clusters.
    """
    _, __, n_clust = clusters.max(axis=0)
    n_clust = int(n_clust)
    for i in range(n_clust):
        two_d_cluster=clusters[clusters[:,2] == i][:,np.array([True, True, False])]
        if not validate_cluster(max_dist,two_d_cluster):
            return(False, centers, clusters)
        else:
            continue
    return(True, centers, clusters)

def validate_cluster(max_dist,cluster):
    """
    Validates whether a single cluster satisfies the maximum distance constraint.

    Parameters:
    ----------
    max_dist : float
        The maximum allowable distance between any two points within the cluster.
    cluster : np.ndarray
        A 2D array of shape (n_points, 2) representing the (x, y) coordinates of the points in the cluster.

    Returns:
    -------
    bool
        True if all distances between points in the cluster are within `max_dist`, False otherwise.
    """
    distances = cdist(cluster,cluster, lambda ori,des: int(round(distance_euclid(ori,des))))
    if np.any(np.array(distances) > max_dist):
        return False
    return True

def create_centroids_gdf(facility_gdf, centers, clusters):
    """
    Creates a GeoDataFrame of cluster centroids and associates facilities to their nearest clusters.

    This function calculates centroids for each cluster, associates facilities to their nearest cluster, 
    and counts the number of facilities within each cluster. The centroids and facility data are returned 
    as GeoDataFrames in the original CRS.

    Parameters:
    ----------
    facility_gdf : gpd.GeoDataFrame
        A GeoDataFrame containing the facilities and their locations in geographic coordinates.
    centers : np.ndarray
        A 2D array of shape (n_clusters, 2) representing the coordinates of the cluster centers.
    clusters : np.ndarray
        A 2D array where each row contains the (x, y) coordinates of a point and its cluster label.

    Returns:
    -------
    tuple (gpd.GeoDataFrame, gpd.GeoDataFrame)
        centroids_gdf : gpd.GeoDataFrame
            A GeoDataFrame containing the cluster centroids, their associated supply counts, and cluster indices.
        facility_gdf : gpd.GeoDataFrame
            The input GeoDataFrame with additional cluster assignments and distances to the nearest cluster.
    """
    orig_crs = facility_gdf.crs
    utm_crs = facility_gdf.estimate_utm_crs()
    points_clust = gpd.GeoDataFrame({'cluster': np.array(map(int, clusters[:,2])), 'geometry': [Point(xy) for xy in list(zip(clusters[:,0], clusters[:,1]))]})
    points_clust = points_clust.set_crs(utm_crs)

    facility_gdf_cop = facility_gdf.to_crs(utm_crs).sjoin_nearest(points_clust, distance_col = 'dist')
    facility_gdf_cop = facility_gdf_cop.drop_duplicates(subset = 'ID')
    facility_gdf_cop = facility_gdf_cop.to_crs(orig_crs)

    # Counts how many points are inside each cluster
    supplies = dict(facility_gdf_cop.cluster.value_counts())

    centroids_gdf = gpd.GeoDataFrame({'geometry':[Point(xy) for xy in centers], 
                                    'supply':[supplies[key] for key in sorted(supplies.keys())],
                                    'clust_ind': sorted(list(map(int, supplies.keys())))}, geometry = 'geometry')

    centroids_gdf['ID'] = np.arange(len(centroids_gdf))
    
    centroids_gdf = centroids_gdf.set_crs(utm_crs)
    centroids_gdf = centroids_gdf.to_crs(orig_crs)
    
    return(centroids_gdf, facility_gdf_cop)
