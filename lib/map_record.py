# -*- coding: utf-8 -*-
import os
import yaml
import pyproj
import math
import thinplatespline
from copy import deepcopy
import hashlib
import json
        
def densify_linestring(points):
    x, y = zip(*points)
    max_range = min(max(x) - min(x), max(y) - min(y)) / 20
    result = []
    x0, y0 = points[0]
    for x, y in points:
        dx = x - x0
        dy = y - y0
        r = math.sqrt(dx * dx + dy * dy)
        result.append((x0, y0))
        if r > max_range:
            ax = dx / r
            ay = dy / r
            t = max_range
            while t < r:
                xi = x0 + ax * t
                yi = y0 + ay * t
                result.append((xi, yi))
                t += max_range
        x0, y0 = x, y
    result.append((x0, y0))
    return result

class Maprecord(object):
    def __init__(self, filename, data=None):
        self._filename = os.path.abspath(filename)
        self._data = data
        if data is not None:
            self._check_data()
        self._gcps = None
        self._projected_cutline = None
        self._image_path = None
        self._gcp_transformer = None
        self._inv_gcp_transformer = None
        self._fingerprint = None
    
    @property
    def data(self):
        if self._data is None:
            self._data = yaml.load(open(self._filename, 'rb'), Loader=yaml.CLoader)
            self._check_data()
        return self._data
        
        
    @property
    def image_path(self):
        if self._image_path is None:
            base_dir = os.path.dirname(self._filename)
            self._image_path = os.path.abspath(os.path.join(base_dir, self.data['image_path']))
        return self._image_path

    @property
    def mask_path(self):
        mask_path = self.data.get('mask_path')
        if mask_path:
            base_dir = os.path.dirname(self._filename)
            mask_path = os.path.abspath(os.path.join(base_dir, mask_path))
        return mask_path or None


    @property    
    def gcps(self):
        if self._gcps is None:
            self._gcps = []
            for gcp in self.data['gcps']:
                pixel = gcp['pixel']['x'], gcp['pixel']['y']
                ground = gcp['ground']['x'], gcp['ground']['y']
                if not gcp['is_projected']:
                    ground = self.proj(*ground)
                self._gcps.append({'ground': ground, 'pixel': pixel})
        return self._gcps
        
    @property
    def srs(self):
        return self.data['srs']
        
    @property
    def proj(self):
        return pyproj.Proj(self.srs)        
        
    @property
    def projected_cutline(self):
        if self._projected_cutline is None:
            cutline_srs = self.data['cutline']['srs']
            cutline = [(p['x'], p['y']) for p in self.data['cutline']['points']]
            if cutline[0] != cutline[-1]:
                cutline.append(cutline[0])
            if self.srs != cutline_srs:
                cutline = densify_linestring(cutline)
                if cutline_srs == 'RAW':
                    cutline = [self.gcp_transformer.get_point(x, y) for x, y in cutline]
                else:
                    proj_src = pyproj.Proj(cutline_srs)
                    proj_dst = pyproj.Proj(self.srs)
                    cutline = zip(*pyproj.transform(proj_src, proj_dst, *zip(*cutline)))
            self._projected_cutline  = cutline
        return self._projected_cutline
        
    @property
    def gcp_transformer(self):
        if self._gcp_transformer is None:
            points = [[gcp['pixel'][0], gcp['pixel'][1], gcp['ground']] for gcp in self.gcps]
            self._gcp_transformer = thinplatespline.TPS(points, 2)
        return self._gcp_transformer
        
    @property
    def inv_gcp_transformer(self):
        if self._inv_gcp_transformer is None:
#            try:
            points = [[gcp['ground'][0], gcp['ground'][1], gcp['pixel']] for gcp in self.gcps]
#            except KeyError:
#                 import ipdb; ipdb.set_trace()
            self._inv_gcp_transformer = thinplatespline.TPS(points, 2)
        return self._inv_gcp_transformer

    def _check_data(self):
        gcps = self.data['gcps']
        cutline = self.data.get('cutline', None)
        if len(gcps) < 3:
            raise ValueError('Too few gcps')
        for gcp in gcps:
            if sorted(gcp.keys()) != ['ground', 'is_projected', 'pixel']:
                raise ValueError('Wrong gcp format')
            if sorted(gcp['ground'].keys()) != ['x', 'y']:
                raise ValueError('Wrong gcp format (ground field)')
            if sorted(gcp['pixel'].keys()) != ['x', 'y']:
                raise ValueError('Wrong gcp format (pixel field)')
        if cutline is not None:
            cutline_srs = cutline.get('srs')
            if not cutline_srs:
                raise Exception('Cutline must have srs or RAW for pixel coordinates')
            cutline_points = cutline['points']
            if len(cutline_points) < 3:
                raise ValueError('Too few points in cutline')
            for point in cutline_points:
                if sorted(point.keys()) != ['x', 'y']:
                    raise ValueError('Wrong cutline format')
        

    def write(self, filename, image_path_relative=True):
        data = deepcopy(self.data)
        image_path = self.image_path
        if image_path_relative:
            image_path = os.path.relpath(image_path, os.path.dirname(filename))
        data['image_path'] = image_path
        s = yaml.dump(data, Dumper=yaml.CSafeDumper, indent=4, width=999)
        with open(filename, 'wb') as f:
            f.write(s)
            
    @property
    def fingerprint(self):
        if self._fingerprint is None:
            data = deepcopy(self.data)
            del data['image_path']
            data['image_size'] = os.path.getsize(self.image_path)
            data['image_mtime'] = os.path.getmtime(self.image_path)
            data = json.dumps(data)
            self._fingerprint = hashlib.sha1(data).hexdigest()
        return self._fingerprint
        
#def write(filename, image_path, srs, gcps, cutline=None, image_path_relative=True):
#    if not os.path.isabs(image_path):
#        raise ValueError('Image path must be absolute: %s' % image_path)
#    maprec_dir = os.path.dirname(filename)
#    if image_path_relative:
#        image_path = os.path.relpath(image_path, maprec_dir)
#    if len(gcps) < 3:
#        raise ValueError('Too few gcps')
#    for gcp in gcps:
#        if sorted(gcp.keys()) != ['ground', 'is_projected', 'pixel']:
#            raise ValueError('Wrong gcp format')
#    cutline = cutline or None
#    if cutline:
#        cutline_srs = cutline.get('srs')
#        if not cutline_srs:
#            raise Exception('Cutline must have srs or RAW for pixel coordinates')
#        cutline_points = cutline['points']
#        if len(cutline_points) < 3:
#            raise ValueError('Too few points in cutline')
#        for point in cutline_points:
#            if sorted(point.keys()) != ['x', 'y']:
#                raise ValueError('Wrong cutline format')
#    data = {
#            'image_path': image_path,
#            'srs': srs,
#            'gcps': gcps,
#            'cutline': {'points': cutline, 'srs': cutline_srs},
#            }
#    s = yaml.dump(data, Dumper=yaml.CSafeDumper, indent=4, width=999)
#    with open(filename, 'wb') as f:
#        f.write(s)

#if __name__ == '__main__':
#    filename = '../001m--a49.maprec'
#    mr = MaprecordReader(filename)
#    print mr.projected_cutline
