import numpy as np


def point_in_polygon(point, polygon):
    """
    point : [longitude, latitude]

    polygon : [(lon1, lat1), (lon2, lat2), ..., (lonn, latn)]

    """
    from osgeo import ogr,osr

    # Create spatialReference
    spatialReference = osgeo.osr.SpatialReference()
    spatialReference.SetWellKnownGeogCS("WGS84")

    # Create ring
    ring = osgeo.ogr.Geometry(osgeo.ogr.wkbLinearRing)

    # Add points
    for lon, lat in polygon:
        ring.AddPoint(lon, lat)

    # Create polygon
    poly = osgeo.ogr.Geometry(osgeo.ogr.wkbPolygon)
    poly.AssignSpatialReference(spatialReference)
    poly.AddGeometry(ring)

    # Create point
    pt = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
    pt.AssignSpatialReference(spatialReference)
    pt.SetPoint(0, point[0], point[1])

    return pt.Within(poly)

def plot_points(lats,lons,r_points):

	min_lon = np.min(lons)
	min_lat = np.min(lats)
	max_lon = np.max(lons)
	max_lat = np.max(lats)


	DOMAIN = np.array([
		[min_lon-5, max_lon+5],
		[min_lat-5, max_lat+5]])

	m = Basemap(resolution='l',
		projection='merc',
		llcrnrlon=DOMAIN[0][0],     # lower-left corner longitude
		llcrnrlat=DOMAIN[1][0],     # lower-left corner latitude
		urcrnrlon=DOMAIN[0][1],     # upper-right corner longitude
		urcrnrlat=DOMAIN[1][1],     # upper-right corner latitude
		area_thresh=1000.0,
		suppress_ticks=True)
	m.drawcoastlines()

	# Fill the globe with a blue color
	m.drawmapboundary(fill_color='aqua')
	# Fill the continents with the land color
	m.fillcontinents(color='coral', lake_color='aqua')

	delta_lon = 5
	delta_lat = 2
	#meridians = np.arange(DOMAIN[0][0], DOMAIN[0][1] + delta_lon, delta_lon)
	#m.drawmeridians(meridians, labels=[1, 0, 0, 1])
	#parallels = np.arange(DOMAIN[1][0], DOMAIN[1][1] + delta_lat, delta_lat)
	#m.drawparallels(parallels, labels=[1, 0, 0, 1])
	x, y = m(r_points[:,0], r_points[:,1])
	m.plot(x, y, 'o', color='Indigo', markersize=4)

	x_corner, y_corner = m(lons,lats)
	m.plot(x_corner,y_corner, 'p', color='black', markersize=6)



def generate_points(N,lons,lats):

	min_lon = np.min(lons)
	min_lat = np.min(lats)
	max_lon = np.max(lons)
	max_lat = np.max(lats)

	random_points = np.random.rand(N,2)

	random_points[:,0] = random_points[:,0]*(max_lon-min_lon)+min_lon
	random_points[:,1] = random_points[:,1]*(max_lat-min_lat)+min_lat

#	plot_points(lats,lons,r_points)

	return random_points

br = (25.767368,-80.18930)
bl = (34.088808,-118.40612)
ul = (45.727093,-120.97864)
ur = (40.727093,-73.97864)

polygon = [br,ur,ul,bl,br]

lats = [x[0] for x in polygon]
lons = [x[1] for x in polygon]

r_points = generate_points(100,lats,lons)

booleans = [ point_in_polygon(tuple(x),polygon) for x in r_points ]

for x in range(0,len(r_points)):
	print(*r_points[x],sep=',')


print(np.sum(booleans))

for x in range(0,len(r_points)):
	if booleans[x]:
		print(*r_points[x],sep=',')


