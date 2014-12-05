# -*- coding: utf-8 -*-
import argparse
import os
from lib import ozi_reader
from lib import map_record


def find_image_file(ozi_image_filename, base_dir):
    ozi_image_filename = ozi_image_filename.split('\\')[-1].lower()
    path = None
    for fn in os.listdir(base_dir):
        if os.path.isfile(os.path.join(base_dir, fn)) and fn.lower() == ozi_image_filename:
            if path is not None:
                raise Exception('Ambigios file name "%s"' % fn)
            path = os.path.abspath(os.path.join(base_dir, fn))
    if path is None:
        raise Exception('Image "%s" not found' % ozi_image_filename)
    return path
    
def get_srs_as_proj4(ozi_datum_string, ozi_projection=None):
    datums = {
        'wgs84': '+datum=WGS84',
        'sk42': '+ellps=krass +towgs84=23.9,-141.3,-80.9,0,-0.37,-0.85,-0.12 +no_defs'
        #'sk42': '+ellps=krass +no_defs'
        }
    datum_map = {
        # http://www.hemanavigator.com.au/Products/TopographicalGPS/OziExplorer/tabid/69/Default.aspx
        'WGS 84': datums['wgs84'],
        'Pulkovo 1942 (1)': datums['sk42'],
        'Pulkovo 1942 (2)': datums['sk42'],
        'Pulkovo 1942': datums['sk42']
        }
    proj_map = {
        'Latitude/Longitude':
            ('+proj=latlong',),
        'Transverse Mercator':
            ('+proj=tmerc +units=m',
            ('+lat_0=', 'lat_origin'),
            ('+lon_0=', 'lon_origin'),
            ('+k=', 'k_factor'),
            ('+x_0=', 'false_easting'),
            ('+y_0=', 'false_northing'))}
    try:
        srs = [datum_map[ozi_datum_string]]
    except KeyError as e:
        raise Exception("Unsupported datum: '%s'" % (e.message))
    if ozi_projection is None:
        proj_params = ('+proj=latlong',)
    else:
        try:
            proj_params = proj_map[ozi_projection.name]
        except KeyError as e:
            raise Exception("Unsupported projection: '%s'" % (e.message))
    srs.append(proj_params[0])
    for proj_param, ozi_param in proj_params[1:]:
        try:
            srs.append('%s%s' % (proj_param, ozi_projection.get(ozi_param)))
        except KeyError:
            raise Exception('Projection parameters not fully defined')
    return ' '.join(srs)

def convert_gcps(ozi_gcps):
    gcps = []
    for ozi_gcp in ozi_gcps:
        gcp = {
            'pixel': {'x': ozi_gcp['pixel']['x'], 'y': ozi_gcp['pixel']['y']},
            'ground': {'x': ozi_gcp['ref']['x'], 'y': ozi_gcp['ref']['y']},
            'is_projected': ozi_gcp.type == 'proj'
        }
        gcps.append(gcp)
    return gcps
        
def convert_cutline(ozi_cutline):
    return [{'x': x, 'y': y} for x, y in ozi_cutline]    
    
def get_maprecord_from_ozi_data(ozi_map_file):
    maprecord = {}
    ozi_map = ozi_reader.read_map(open(ozi_map_file).readlines())
    maprecord['image_path'] = find_image_file(ozi_map.file_name, os.path.dirname(ozi_map_file) or '.')
    maprecord['srs'] = get_srs_as_proj4(ozi_map.datum, ozi_map.projection)
    maprecord['gcps'] = convert_gcps(ozi_map.gcps)
    maprecord['cutline'] = {
        'srs': get_srs_as_proj4(ozi_map.datum),
        'points': convert_cutline(ozi_map.cutline)
        }
    return maprecord

def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument('in_file', metavar='ozi_file.map')
    parser.add_argument('out_file', metavar='output.maprec')
    parser.add_argument('--abs-path', action='store_true', help='Write absolute path to image file')
    return parser.parse_args()

def main():
    args = parse_command_line()
#    print args
    maprecord = get_maprecord_from_ozi_data(args.in_file)
#    print maprecord
    maprecord = map_record.Maprecord(args.in_file, maprecord)
    maprecord.write(args.out_file, image_path_relative=not args.abs_path)
#def write(filename, image_path, srs, gcps, cutline_srs=None, cutline=None, image_path_relative=True):    

    
if __name__ == '__main__':
    main()