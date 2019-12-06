import rasterio
from rasterio.windows import Window
from rasterio import features
from rasterio.io import MemoryFile
from rasterio.crs import CRS
from rasterio import Affine
from osgeo import osr
from shapely.geometry import Polygon
import geopandas
import numpy

def raster_subset(fil, bounds, resolution):

    with rasterio.open(fil) as src:
        prof = src.profile
    
    col_offset = (bounds[0] - prof['transform'][2])/resolution[1]
    row_offset = (prof['transform'][5] - bounds[3])/resolution[1]

    w = int((bounds[2]-bounds[0])/resolution[1])
    h = int((bounds[3]-bounds[1])/resolution[1])

    sub = Window(col_offset, row_offset, w, h)

    with rasterio.open(fil) as src:
        subset = src.read(1, window=sub)
    
    profile = prof
    
    profile['transform'] = Affine(resolution[1], 0.0, bounds[0], 0.0, resolution[0], bounds[3])
    profile['width'] = w
    profile['height'] = h
    
    return subset, profile

def polygon_till_raster(fil, bounds, resolution):
    mask = geopandas.read_file(fil)
    xmin, xmax, ymin, ymax = bounds[0], bounds[2], bounds[1], bounds[3] 
    boundsxy = Polygon( [(xmin,ymin), (xmin, ymax), (xmax, ymax), (xmax,ymin)] )

    mask['geometry'] = mask['geometry'].intersection(boundsxy)

    w = int((bounds[2]-bounds[0])/resolution[1])
    h = int((bounds[3]-bounds[1])/resolution[1])

    memfile = MemoryFile()

    profile = {'driver': 'GTiff', 'dtype': 'uint8', 'nodata': None, 'width': w, 'height': h, 'count': 1, 
               'crs': CRS.from_epsg(3006), 'transform': Affine(resolution[1], 0.0, bounds[0],
                0.0, resolution[0], bounds[3]), 'tiled': False, 'interleave': 'band'}

    with rasterio.open(memfile, 'w+', **profile) as out:
        out_arr = out.read(1)
        shapes = ((geom, value) for geom, value in zip(mask.geometry, mask.Id))
        try:
            burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
            out.write_band(1, burned)
        except ValueError:
            pass
    
    mask_fil = memfile.open(**profile).read(1)
    
    return mask_fil, profile

def omvandla(original_system, nytt_system, x, y):
    
    source = osr.SpatialReference()
    target = osr.SpatialReference()
    source.ImportFromEPSG(original_system)
    target.ImportFromEPSG(nytt_system)
    transform = osr.CoordinateTransformation(source, target)    
    ny_ref = transform.TransformPoint(x, y)

    return ny_ref

def jamforelse(tidsserie, vattenveg, vm_lm, vm_ndwi, ndvi_threshold, ndwi_threshold):
    
    maximum = tidsserie.max(dim='time').values
    maximum[vm_lm == 0] = 0
    maximum[vm_ndwi >= ndwi_threshold] = 0
    maximum[maximum <= ndvi_threshold] = 0
    maximum[maximum > 0] = 1
    
    vveg = vattenveg.astype('float')
    maximum = maximum.astype('float')
    
    temp = vveg+maximum
    temp[temp == 2] = 3
    temp[temp < 3] = 0
    
    mask1 = numpy.copy(maximum)
    mask2 = numpy.copy(vveg)
    mask = mask1+mask2
    mask[mask > 0] = 1
    
    resultat = temp+(numpy.where(vveg > 0, vveg, vveg+1)+maximum)
    resultat[resultat > 3] = 3
    resultat = resultat*mask
    
    return resultat