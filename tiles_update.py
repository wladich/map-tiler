#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse
import time
from lib import map_record
from lib.mpimap import mpstarimap
import ozi_to_maprec
import pyproj
import collections
import image_store
import Image
import ImageDraw
import ImageChops
import warnings
from itertools import chain

warnings.simplefilter('ignore', Image.DecompressionBombWarning)

DEBUG = False
proj_gmerc = pyproj.Proj('+init=epsg:3785')
proj_gmerc_180 = pyproj.Proj(
    '+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=180.0 +x_0=20037508.342789244 +y_0=0 +k=1.0 '
    '+units=m +nadgrids=@null +wktext  +no_defs')
max_gmerc_coord = 20037508.342789244
METATILE_DELTA = 3

config = None
tile_store = None


class UserInputError(Exception):
    pass


def tile_size_in_gmerc_meters(level):
    return max_gmerc_coord * 2 / 2 ** level


def tile_from_gmerc_meters(x, y, level):
    tile_size = tile_size_in_gmerc_meters(level)
    tx = (x + max_gmerc_coord) / tile_size
    ty = (-y + max_gmerc_coord) / tile_size
    return int(tx), int(ty)


def tile_nw_corner(tile_x, tile_y, level):
    tile_size = tile_size_in_gmerc_meters(level)
    x = tile_x * tile_size - max_gmerc_coord
    y = -tile_y * tile_size + max_gmerc_coord
    return x, y


def open_map_reference(filename):
    """file can be  in formats:
       1) maprecord yaml (or json)
       2) ozi .map
    """
    if filename.endswith('.map'):
        data = ozi_to_maprec.get_maprecord_from_ozi_file(filename)
        maprecord = Maprecord(filename, data)
    else:
        maprecord = Maprecord(filename)
    return maprecord


def reproject_cutline_gmerc(src_proj, points):
    if points:
        test_point = points[0]
        test_point = pyproj.transform(src_proj, proj_gmerc, *test_point)
        if abs(test_point[0]) > max_gmerc_coord / 2:
            dest_proj = proj_gmerc_180
        else:
            dest_proj = proj_gmerc
        return zip(*pyproj.transform(src_proj, dest_proj, *zip(*points)))
    return []


class JobManager(object):
    prev_state_filename = 'tiles_sources'

    def __init__(self, maps_references, tile_zoom):
        self.tile_zoom = tile_zoom
        self._tiles_needing_update = None
        self._removed_tiles = None

        self.old_maps_fingerprints = set()
        self.new_maps_fingerprints = set()

        self.old_maps_fingerprints_for_tiles = {}
        for x, y, z, fingerprints in self._load_state_file():
            self.old_maps_fingerprints_for_tiles[(x, y, z)] = fingerprints
            self.old_maps_fingerprints.update(fingerprints)

        self.new_maps_fingerprints_for_tiles = collections.defaultdict(list)
        self.maps_references_for_tiles = collections.defaultdict(list)
        self.maps_references_for_tiles = collections.defaultdict(list)
        for map_reference in maps_references:
            maprecord = open_map_reference(map_reference)
            fingerprint = maprecord.fingerprint
            tiles = self._get_tiles_for_maprecord(maprecord)
            for (x, y, z) in tiles:
                self.new_maps_fingerprints_for_tiles[(x, y, z)].append(fingerprint)
                self.maps_references_for_tiles[(x, y, z)].append(map_reference)
                self.new_maps_fingerprints.add(fingerprint)

    def _get_tiles_for_maprecord(self, maprecord):
        cutline = maprecord.projected_cutline
        cutline = reproject_cutline_gmerc(maprecord.proj, cutline)
        x, y = zip(*cutline)
        tx1, ty2 = tile_from_gmerc_meters(min(x), min(y), self.tile_zoom)
        tx2, ty1 = tile_from_gmerc_meters(max(x), max(y), self.tile_zoom)
        return [(x_, y_, self.tile_zoom) for x_ in xrange(tx1, tx2+1) for y_ in xrange(ty1, ty2+1)]

    def get_tiles_needing_update(self):
        if self._tiles_needing_update is None:
            self._tiles_needing_update = []
            for k, v in self.new_maps_fingerprints_for_tiles.items():
                if self.old_maps_fingerprints_for_tiles.get(k) != v:
                    self._tiles_needing_update.append((k, self.maps_references_for_tiles[k]))
        return self._tiles_needing_update

    def get_added_maps_n(self):
        return len(self.new_maps_fingerprints - self.old_maps_fingerprints)

    def get_removed_maps_n(self):
        return len(self.old_maps_fingerprints - self.new_maps_fingerprints)

    def get_removed_tiles(self):
        if self._removed_tiles is None:
            self._removed_tiles = []
            for k in self.old_maps_fingerprints_for_tiles:
                if k not in self.new_maps_fingerprints_for_tiles:
                    self._removed_tiles.append(k)
        return self._removed_tiles

    def apply_update(self):
        s = []
        for tile, fingerprints in self.new_maps_fingerprints_for_tiles.items():
            s.append('%s,%s,%s,%s\n' % (tile + ('-'.join(fingerprints),)))
        s = ''.join(s)
        tile_store.write_metadata(self.prev_state_filename, s)

    def _load_state_file(self):
        state = tile_store.read_metadata(self.prev_state_filename)
        if state:
            for line in state.splitlines():
                x, y, z, fingerprints = line.strip().split(',')
                x = int(x)
                y = int(y)
                z = int(z)
                fingerprints = fingerprints.split('-')
                yield x, y, z, fingerprints


def get_reprojected_image(tile_x, tile_y, level, map_reference):
    metatile_delta = config.max_level - level
    tile_size = 256 * (2 ** metatile_delta)
    tile_origin = tile_nw_corner(tile_x, tile_y, level)
    dest_pixel_size = tile_size_in_gmerc_meters(level) / tile_size
    maprecord = open_map_reference(map_reference)

    def transform_dest_to_src_pixel((x, y)):
        x = x * dest_pixel_size + tile_origin[0]
        y = tile_origin[1] - y * dest_pixel_size
        x, y = pyproj.transform(proj_gmerc, maprecord.proj, x, y)
        x, y = maprecord.inv_gcp_transformer.transform(x, y)
        return x, y
    im_src = Image.open(maprecord.image_path)
    src_has_alpha = im_src.mode.endswith('A')
    if src_has_alpha:
        im_src = im_src.convert('RGBA')
    else:
        im_src = im_src.convert('RGB')
    if maprecord.mask_path is not None:
        im_mask = Image.open(maprecord.mask_path)
        im_src.putalpha(im_mask)
        src_has_alpha = True
    cell_size = 64
    mesh = []
    for cell_x in xrange(tile_size / cell_size):
        for cell_y in xrange(tile_size / cell_size):
            x1 = cell_x * cell_size
            y1 = cell_y * cell_size
            x2 = x1 + cell_size
            y2 = y1 + cell_size
            quad = (x1, y1), (x1, y2), (x2, y2), (x2, y1)
            quad = map(transform_dest_to_src_pixel, quad)
            quad = sum(quad, tuple())
            mesh.append(((x1, y1, x2, y2), quad))
    im = im_src.transform((tile_size, tile_size), Image.MESH, mesh, Image.BICUBIC)
    mask = Image.new('L', (tile_size, tile_size))
    cutline = maprecord.projected_cutline
    cutline = reproject_cutline_gmerc(maprecord.proj, cutline)
    cutline = [((x - tile_origin[0]) / dest_pixel_size, (tile_origin[1] - y) / dest_pixel_size) for x, y in cutline]
    draw = ImageDraw.Draw(mask)
    draw.polygon(cutline, fill=255, outline=255)
    if not src_has_alpha:
        im.putalpha(mask)
    im.paste(0, (0, 0), ImageChops.invert(mask))
    return im


def slice_metatile(im, metatile_x, metatile_y, dest_level):
    meta_delta = dest_level - config.metatile_level
    meta_q = 2**meta_delta
    assert im.size == (256 * (meta_q),) * 2
    tile_x0 = metatile_x * meta_q
    tile_y0 = metatile_y * meta_q
    max_tile = 2 ** dest_level
    for d_tile_y in xrange(meta_q):
        y0 = d_tile_y * 256
        for d_tile_x in xrange(meta_q):
            x0 = d_tile_x * 256
            im2 = im.crop([x0, y0, x0+256, y0+256])
            tile_store.write(im2, (tile_x0 + d_tile_x) % max_tile, tile_y0 + d_tile_y, dest_level)


def process_metatile((tile_x, tile_y, metatile_level), map_references):
    # 1. render
    # 2. cut tiles
    # 3. shrink and cut overviews
    im = None
    for map_reference in map_references:
        im2 = get_reprojected_image(tile_x, tile_y, metatile_level, map_reference)
        if im is None:
            im = im2
        else:
            im.paste(im2, (0, 0), im2)
    for level in xrange(config.max_level, metatile_level, -1):
        slice_metatile(im, tile_x, tile_y, level)
        im = im.resize((im.size[0] / 2, )*2, Image.ANTIALIAS)
    max_tile = 2 ** metatile_level
    tile_store.write(im, tile_x % max_tile, tile_y, metatile_level)


def remove_tiles(metatiles):
    tiles_overviews = set()
    for tile in metatiles:
        tile_store.remove(*tile)
        tx, ty, level = tile
        while level > 0:
            level -= 1
            tx /= 2
            ty /= 2
            tiles_overviews.add((tx, ty, level))
        tx, ty, level = tile
        w = 1
        for level in xrange(level + 1, config.max_level + 1):
            tx *= 2
            ty *= 2
            w *= 2
            for x in xrange(tx, tx + w):
                for y in xrange(ty, ty + w):
                    tile_store.remove(x, y, level)

    for tile in tiles_overviews:
        tile_store.remove(*tile)


def make_tiles_from_metalevel_to_maxlevel(tiles):
    n = 0
    for _ in mpstarimap(process_metatile, tiles, _nomp=DEBUG):
        n += 1
        print ('\r%.1f%%' % (n * 100. / len(tiles))),
        sys.stdout.flush()
    print


def build_overview(x, y, z):
    im = Image.new('RGBA', (512, 512))
    im2 = tile_store.open(x * 2, y * 2, z + 1)
    if im2 is not None:
        im.paste(im2.convert('RGBA'), (0, 0))
    im2 = tile_store.open(x * 2 + 1, y * 2, z + 1)
    if im2 is not None:
        im.paste(im2.convert('RGBA'), (256, 0))
    im2 = tile_store.open(x * 2, y * 2 + 1, z + 1)
    if im2 is not None:
        im.paste(im2.convert('RGBA'), (0, 256))
    im2 = tile_store.open(x * 2 + 1, y * 2 + 1, z + 1)
    if im2 is not None:
        im.paste(im2.convert('RGBA'), (256, 256))
    im = im.resize((256, 256), Image.ANTIALIAS)
    tile_store.write(im, x, y, z)


def build_overviews(altered_tiles):
    need_update = set()
    for x, y, z in altered_tiles:
        if z > 0:
            need_update.add((x / 2, y / 2, z - 1))
    n = 0
    for _ in mpstarimap(build_overview, need_update, _nomp=DEBUG):
        n += 1
        print '\r%.1f%%' % (n * 100. / len(need_update)),
        sys.stdout.flush()
    print
    if need_update:
        build_overviews(need_update)


def filename_arg_type(s):
    return s.decode(sys.getfilesystemencoding())


def parse_image_format(s):
    fields = s.split(',')
    fmt = fields[0].upper()
    params = fields[1:]
    params = (p.split('=', 1) for p in params)
    params = dict((k, int(v)) for k, v in params)
    if fields[0] == 'PNG8':
        if not (2 <= params.get('colors', None) <= 256):
            raise ValueError
        if not (1 <= params.get('speed', None) <= 10):
            raise ValueError
        if not (1 <= params.get('compression', 6) <= 9):
            raise ValueError
    elif fields[0] == 'PNG32':
        if not (1 <= params.get('compression', 6) <= 9):
            raise ValueError
    elif fields[0] == 'JPEG':
        if not (1 <= params.get('quality', None) <= 100):
            raise ValueError
    else:
        raise ValueError('Unsupported image format "%s"' % fields[0])
    return fmt, params


class MyArgumentParser(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_line):
        arg_line = arg_line.strip()
        if arg_line and arg_line[0] == '#':
            return []
        if '=' in arg_line and not arg_line.startswith('-'):
            arg_line = '--' + arg_line
        if arg_line:
            return [arg_line]
        else:
            return []


def parse_command_line():
    # TODO: сделать возможны включение файлов при вызове не из текущей директории
    # TODO: по умолчанию включать файл tiles.cfg из текущего каталога
    # TODO: добавить описание форматов изображений
    parser = MyArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('--storage-format', choices=['files', 'mbtiles'], required=True)
    parser.add_argument('--max-level', type=int, required=True)
    parser.add_argument('maps', metavar='FILE', type=filename_arg_type, nargs='+')
    parser.add_argument('--image-format', type=parse_image_format, required=True)
    parser.add_argument('--image-border-format', type=parse_image_format)
    parser.add_argument('--out', metavar='PATH', dest='out_path', required=True,
                        help='Filename of mbtiles container or tiles dir')
    parser.add_argument('-f', dest='do_update', action='store_true', help='actualy do update')
    return parser.parse_args()


def print_job_info(job):
    print 'Maps added: %s' % job.get_added_maps_n()
    print 'Maps removed: %s' % job.get_removed_maps_n()
    print 'Metatiles to rebuild: %s' % len(job.get_tiles_needing_update())


def configure_output_storage():
    global tile_store
    image_encoder = image_store.get_image_encoder(config.image_format, config.image_border_format)
    tile_store_class = {'files': image_store.FilesWriter, 'mbtiles': image_store.MBTilesWriter}[config.storage_format]
    tile_store = tile_store_class(config.out_path, image_encoder)


def main():
    global config
    config = parse_command_line()
    config.metatile_level = max(config.max_level - METATILE_DELTA, 0)
    configure_output_storage()
    job = JobManager(config.maps, config.metatile_level)
    print_job_info(job)
    if config.do_update:
        t = time.time()
        remove_tiles(job.get_removed_tiles())
        make_tiles_from_metalevel_to_maxlevel(job.get_tiles_needing_update())
        print 'Building overviews'
        overview_tiles_needing_update = (tile for (tile, _) in job.get_tiles_needing_update())
        overview_tiles_needing_update = chain(overview_tiles_needing_update, job.get_removed_tiles())
        overview_tiles_needing_update = ((x % (2**z), y, z) for (x, y, z) in overview_tiles_needing_update)
        build_overviews(overview_tiles_needing_update)
        job.apply_update()
        print 'Done in %.1f seconds' % (time.time() - t)
    else:
        print 'To actualy update mbtiles db, specify "-f" option'

if __name__ == '__main__':

    try:
        main()
    except UserInputError as e:
        print e
        exit(1)
