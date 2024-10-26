"""
This file contains functions and constants needed for the geometry of 
a sphere
"""
import numpy as np


__author__ = "Paolo Scaccia, Paolo Antonelli and Stefano Piani"
__copyright__ = "Copyright 2017, AdaptiveMeteo S.r.l."
__credits__ = ["Paolo Scaccia", "Paolo Antonelli"]
__license__ = "--"
__version__ = "0.0.1"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"
__status__ = "Development"

EARTH_RADIUS = 6372.795
f = 1/298.257223563
a = 6378137.0    # in m
b = (1-f)*a      # in m
e_squared = 1 - (b/a)**2

def define_spatial_domain(lats,lons):

    lon_min,lon_max = lons.min(),lons.max()
    lat_min,lat_max = lats.min(),lats.max()
    
    # define domain
    domain = [ (lon_min , lats[lons==lon_min][0] ),   # bottom left
		(lons[lats==lat_max][0],lat_max),   # upper left
		(lon_max,lats[lons==lon_max][0]),   # upper right
		(lons[lats==lat_min][0],lat_min), # bottom right
		(lon_min , lats[lons==lon_min][0])]   # starting point

    # from np.float to float (needed by osgeo module)
    domain = [ (float(x[0]),float(x[1])) for x in domain ]
    
    return domain



def point_in_polygon(point, polygon):
    """
    point : [longitude, latitude]

    polygon : [(lon1, lat1), (lon2, lat2), ..., (lonn, latn)]

    """
    from osgeo import ogr,osr

    # Create spatialReference
    spatialReference = osr.SpatialReference()
    spatialReference.SetWellKnownGeogCS("WGS84")

    # Create ring
    ring = ogr.Geometry(ogr.wkbLinearRing)

    # Add points
    for lon, lat in polygon:
        ring.AddPoint(lon, lat)

    # Create polygon
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AssignSpatialReference(spatialReference)
    poly.AddGeometry(ring)

    # Create point
    pt = ogr.Geometry(ogr.wkbPoint)
    pt.AssignSpatialReference(spatialReference)
    pt.SetPoint(0, point[0], point[1])

    return pt.Within(poly)



def haversine(lon1, lat1, lon2, lat2):    
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    EARTH_RADIUS = 6372.795

    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    km = EARTH_RADIUS * c
    
    return km


def _wp84_to_cartesian(lon, lat):
    """
    Starting from longitude and latitude of a point, this function returns
    the cartesian cardinate of the point as if it was on a sphere with radius=1.
    The point (0,0,1) is the north pole, the point (1,0,0) has both longitude
    and latitude equal to zero.

    Args:
        - lon: the longitude of the point
        - lat: the latitude of the point

    Returns:
        - a tuple (x,y,z) with the cartesian coordinates of the point
    """
    azimuth_term = np.cos(lat * np.pi / 180.)
    x = np.cos(lon * np.pi / 180.) * azimuth_term
    y = np.sin(lon * np.pi / 180.) * azimuth_term
    z = np.sin(lat * np.pi / 180.)

    return x, y, z


def _cartesian_to_wp84(x, y, z):
    """
    Starting from the cartesian coordinates of a point, this function returns
    the longitude and the latitude of the point.
    The point (0,0,1) is the north pole, the point (1,0,0) has both longitude
    and latitude equal to zero.

    Args:
        - x: the x coordinate of the point
        - y: the y coordinate of the point
        - z: the z coordinate of the point

    Returns:
        - a tuple (lon, lat) with the longitude and the latitude of the point
    """
    r = np.sqrt(x*x + y*y + z*z)
    lat = np.arcsin(z/r)
    lon = np.arctan2(y, x)

    return np.rad2deg(lon), np.rad2deg(lat)


class GreatCircle(object):
    """
    Given two points that are not antipodal on a sphere, it returns the great
    circle that connect the two points (the circle which has the center of the
    sphere as a center)
    """
    def __init__(self, lon1, lat1, lon2, lat2):
        x1, y1, z1 = _wp84_to_cartesian(lon1, lat1)
        x2, y2, z2 = _wp84_to_cartesian(lon2, lat2)

        # Compute the cross product of the two vectors
        self._cross = np.cross(np.array([x1, y1, z1]), np.array([x2, y2, z2]))

        # Check if the points are antipodal
        if np.linalg.norm(self._cross) < 1e-5:
            raise ValueError('Two antipodals point can not be used to define a '
                             'great circle (({}, {}) and ({}, {}))'
                             .format(lat1, lon1, lat2, lon2))

    def __eq__(self, other):
        other_p_lon, other_p_lat = other.poles()[0]
        other_cross = np.array(_wp84_to_cartesian(other_p_lon, other_p_lat))
        if np.linalg.norm(np.cross(self._cross, other_cross)) < 1e-5:
            return True
        return False


    def poles(self):
        """
        Return two points (a point is a tuple with longitude and latitude) that
        are the center of the two semispheres defined by the GreatCircle (for
        example, if this GreatCircle is the equator, it returns the North Pole
        and the South Pole)
        """
        p1 = _cartesian_to_wp84(*self._cross)
        p2 = _cartesian_to_wp84(*-self._cross)
        return p1, p2

    def perpendicular(self, lon, lat):
        """
        Given one point, return a GreatCircle that contains the point and the GreatCircle
	perpendicular to this one
        """
        # We have to find a plane that is perpendicular to the plane where this
        # GreatCircle lies and that pass trough the point
        x, y, z = _wp84_to_cartesian(lon, lat)

        # The plane we are looking for contains the (x, y, z) vector because it
        # contains the point and contains the self._cross vector because it is
        # perpendicular to the current GreatCircle
        v1 = np.array([x, y, z])
        v2 = self._cross

        # If v1 == v2 then the point is on the center of one of the two
        # semisphere defined by the current great circle (like the north pole
        # with the equator) and, therefore, every great circle that pass trough
        # the point is perpendicular to the actual one
        v3 = np.cross(v1, v2)
        v3_norm = np.linalg.norm(v3)

        # If v3 is really small v1 and v2 are on the same line! Therefore, we
        # can change v2 with any random vector such that v3 is not too small
        while v3_norm < 1e-8:
            v2 = np.random.rand(3)
            while np.linalg.norm(v2) < 1e-4:
                v2 = np.random.rand(3)
            v2 /= np.linalg.norm(v2)
            v3 = np.cross(v1, v2)
            v3_norm = np.linalg.norm(v3)

        v3 /= v3_norm

        # If we found a point that is perpendicuar to v3 and lies to the
        # sphere, then it is a point of the great circle we are looking for
        v4 = np.cross(v1, v3)

        v4_lon, v4_lat = _cartesian_to_wp84(v4[0], v4[1], v4[2])

        return GreatCircle(lon, lat, v4_lon, v4_lat)

    def intersect(self, other):
        """
        Return the two points that two GreatCircles share

        Args:
            - other: another GreatCircle

        Returns:
            Two tuples. Each one of them contains the longitude and the latitude
            of one of the two points that belongs to both the GreatCircles
        """
        if self == other:
            raise ValueError('Intersect between two object that represent the'
                             'same great circle is not allowed')
        other_p_lon, other_p_lat = other.poles()[0]
        other_cross = np.array(_wp84_to_cartesian(other_p_lon, other_p_lat))

        p1 = np.cross(self._cross, other_cross)
        p2 = -p1

        p1_wp84 = _cartesian_to_wp84(*p1)
        p2_wp84 = _cartesian_to_wp84(*p2)

        return p1_wp84, p2_wp84

    def distance_from(self, lon, lat):
        """
        Return the distance between this GreatCircle and the point as they were
        on the Earth

        Args:
            - lon: the longitude of the point
            - lat: the latitude of the point
        """
        gc = self.perpendicular(lon, lat)
        p1, p2 = self.intersect(gc)

        d1 = haversine(lon, lat, p1[0], p1[1])
        d2 = haversine(lon, lat, p2[0], p2[1])

        return min(d1, d2)

class SphericalSquare(object):

	def __init__(self, points):

	# First check on the input: it has to be a tuple of four
	# tuples containing lon and lat of the four angles of the
	# spherical square

		if type(points) is not tuple:
			raise ValueError('Input has to be a tuple of points'
					 'each one containing lat and lon')
		for point in points:
			if type(point) is not tuple:
				raise ValueError('Input has to be a tuple of points '
	  					 'each one containing lat and lon')

		first     = points[0]
		second    = points[1]
		third     = points[2]
		fourth    = points[3]


		self._great_circle1 = GreatCircle(first[0],
        	                                  first[1],
                	                          second[0],
                        	                  second[1])

		self._great_circle2 = GreatCircle(second[0],
                        	                  second[1],
						  third[0],
						  third[1])

		self._great_circle3 = GreatCircle(third[0],
						  third[1],
						  fourth[0],
						  fourth[1])

		self._great_circle4 = GreatCircle(fourth[0],
						  fourth[1],
						  first[0],
						  first[1])

		self.d1 = haversine(first[0],first[1],second[0],second[1])
		self.d2 = haversine(second[0],second[1],third[0],third[1])
		self.d3 = haversine(third[0],third[1],fourth[0],fourth[1])
		self.d4 = haversine(first[0],first[1],fourth[0],fourth[1])

	def contains(self,lon,lat):

		flag = 0


		dist = self._great_circle1.distance_from(lon,lat)
		if dist > self.d4 or dist > self.d2 :
			flag = 1

		dist = self._great_circle2.distance_from(lon,lat)
		if dist > self.d1 or dist > self.d3 :
			flag = 1

		dist = self._great_circle3.distance_from(lon,lat)
		if dist > self.d4 or dist > self.d2 :
			flag = 1

		dist = self._great_circle4.distance_from(lon,lat)
		if dist > self.d1 or dist > self.d3 :
			flag = 1

		return not(flag)

def compute_CRadius(lats):
        """
        Compute Prime Vertical Radius of Curvature for the elliptical model of Earth.
        Ref: WGS84
        """        
        lats = np.array(lats)
        sin_lat = np.sin(lats*np.pi/180)
        
        return a/np.sqrt( 1 - (e_squared*sin_lat*sin_lat) )
    
def compute_FOV_position_ECEF(dataset):
    """
    Compute G VECTOR using a VIIRS dataset
    """
    from geometry.utilities.array_reshapers import array_1d

    r_squared = (b*b)/(a*a)    
    datashape = dataset.longs.shape

#    N = compute_CRadius(dataset.lats)
    
    sin_lat = np.sin(dataset.lats  * np.pi /180)
    sin_lon = np.sin(dataset.longs * np.pi /180)
    cos_lat = np.cos(dataset.lats  * np.pi /180)
    cos_lon = np.cos(dataset.longs * np.pi /180)
       
    N = a/np.sqrt( 1 - (e_squared*sin_lat*sin_lat) )
    
    h =  dataset.height
    GX = (N+h)*cos_lat*cos_lon         # (N(lat) + h) * cos(lat) * cos(lon)
    GY = (N+h)*cos_lat*sin_lon         # (N(lat) + h) * cos(lat) * sin(lon)
    GZ =  ((r_squared*N) + h)*sin_lat  # (r_squared * N(lat) + h) * sin_lat
    
    return np.array([ vect for vect in zip(array_1d(GX),
                                           array_1d(GY),
                                           array_1d(GZ)) ]).reshape( datashape + (3,)  ) # return it with dataset shape 