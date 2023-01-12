#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import collections
import hashlib
import json
import math
import os
import sys
import time
import warnings
from itertools import chain
from multiprocessing import Pool, cpu_count

import pyproj
from PIL import Image, ImageChops, ImageDraw, ImageFile, ImageFilter
from maprec import Maprecord, densify_linestring
from ozi_map import ozi_to_maprec

from . import image_store
from . import attribution

PREV_STATE_FILENAME = "tiles_sources"
ImageFile.LOAD_TRUNCATED_IMAGES = True
Image.MAX_IMAGE_PIXELS = None


DEBUG = False
crs_gmerc = pyproj.CRS("+init=epsg:3857")

crs_gmerc_180_dict = crs_gmerc.to_json_dict()
lon_patched = False
easting_patched = False
for rec in crs_gmerc_180_dict["conversion"]["parameters"]:
    if rec["name"] == "Longitude of natural origin":
        rec["value"] = 180
        lon_patched = True

    if rec["name"] == "False easting":
        rec["value"] = 20037508.342789244
        easting_patched = True
assert lon_patched and easting_patched


crs_gmerc_180 = pyproj.CRS.from_user_input(crs_gmerc_180_dict)


max_gmerc_coord = 20037508.342789244
METATILE_DELTA = 3

highlight_color = 0xDB, 0x5A, 0x00

config = None
tile_store = None


class UserInputError(Exception):
    pass


def tile_size_in_gmerc_meters(level):
    return max_gmerc_coord * 2 / 2**level


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
    if filename.endswith(".map"):
        data = ozi_to_maprec.get_maprecord_from_ozi_file(filename)
        maprecord = Maprecord(filename, data)
    else:
        maprecord = Maprecord(filename)
    return maprecord


def calc_area(points):
    if not points:
        return 0
    area = 0
    points = points + [points[0]]
    p1 = points[1]
    for p2 in points[1:]:
        x1, y1 = p1
        x2, y2 = p2
        area += (x2 - x1) * (y1 + y2) / 2
        p1 = p2
    return abs(area)


def reproject_cutline_gmerc(src_crs, points):
    if points:
        points = densify_linestring(points)
        points1 = list(
            zip(
                *pyproj.transform(
                    src_crs, crs_gmerc, *list(zip(*points)), always_xy=True
                )
            )
        )
        points2 = list(
            zip(
                *pyproj.transform(
                    src_crs, crs_gmerc_180, *list(zip(*points)), always_xy=True
                )
            )
        )
        return points1 if calc_area(points1) <= calc_area(points2) else points2
    return []


class JobManager(object):
    def __init__(self, maps_references, tile_zoom, prev_state=None):
        self.tile_zoom = tile_zoom
        self._tiles_needing_update = None
        self._removed_tiles = None

        self.old_maps_fingerprints = set()
        self.new_maps_fingerprints = set()

        self.old_maps_fingerprints_for_tiles = {}
        if prev_state:
            for x, y, z, fingerprints in self._load_state(prev_state):
                self.old_maps_fingerprints_for_tiles[(x, y, z)] = fingerprints
                self.old_maps_fingerprints.update(fingerprints)

        self.new_maps_fingerprints_for_tiles = collections.defaultdict(list)
        self.maps_references_for_tiles = collections.defaultdict(list)
        self.maps_references_for_tiles = collections.defaultdict(list)
        for map_reference in maps_references:
            maprecord = open_map_reference(map_reference)
            fingerprint = self.get_map_fingerprint(maprecord)
            tiles = self._get_tiles_for_maprecord(maprecord)
            for (x, y, z) in tiles:
                self.new_maps_fingerprints_for_tiles[(x, y, z)].append(fingerprint)
                self.maps_references_for_tiles[(x, y, z)].append(map_reference)
                self.new_maps_fingerprints.add(fingerprint)

    def get_map_fingerprint(self, maprecord):
        fingerprint = maprecord.fingerprint
        attrib_filename = attribution.get_attrib_path(maprecord.image_path)
        if os.path.exists(attrib_filename):
            fingerprint = hashlib.sha1(fingerprint.encode())
            fingerprint.update(b":~:" + open(attrib_filename, "rb").read())
            info_filename = attribution.get_info_path(maprecord.image_path)
            fingerprint.update(b":~:" + open(info_filename, "rb").read())
            fingerprint = fingerprint.hexdigest()
        return fingerprint

    def _get_tiles_for_maprecord(self, maprecord):
        cutline = maprecord.projected_cutline
        cutline = reproject_cutline_gmerc(maprecord.crs, cutline)
        x, y = list(zip(*cutline))
        tx1, ty2 = tile_from_gmerc_meters(min(x), min(y), self.tile_zoom)
        tx2, ty1 = tile_from_gmerc_meters(max(x), max(y), self.tile_zoom)
        return [
            (x_, y_, self.tile_zoom)
            for x_ in range(tx1, tx2 + 1)
            for y_ in range(ty1, ty2 + 1)
        ]

    def get_tiles_needing_update(self):
        if self._tiles_needing_update is None:
            self._tiles_needing_update = []
            for k, v in list(self.new_maps_fingerprints_for_tiles.items()):
                if self.old_maps_fingerprints_for_tiles.get(k) != v:
                    self._tiles_needing_update.append(
                        (k, self.maps_references_for_tiles[k])
                    )
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
        for tile, fingerprints in list(self.new_maps_fingerprints_for_tiles.items()):
            s.append("%s,%s,%s,%s\n" % (tile + ("-".join(fingerprints),)))
        s = "".join(s)
        tile_store.write_metadata(self.prev_state_filename, s)

    def _load_state(self, state_data):
        for line in state_data.splitlines():
            x, y, z, fingerprints = line.strip().split(",")
            x = int(x)
            y = int(y)
            z = int(z)
            fingerprints = fingerprints.split("-")
            yield x, y, z, fingerprints


def apply_attribution(im, maprecord, src_to_dest_transformer, dest_meters_in_pixel):
    attrib_filename = attribution.get_attrib_path(maprecord.image_path)
    if not os.path.exists(attrib_filename):
        return im
    attrib_data = json.load(open(attrib_filename))
    image_info = json.load(open(attribution.get_info_path(maprecord.image_path)))
    attr_im_info = attribution.make_attribution_image(
        attrib_data, image_info, src_to_dest_transformer, dest_meters_in_pixel
    )
    x, y = src_to_dest_transformer(
        attrib_data["anchor"]["x"], attrib_data["anchor"]["y"]
    )
    x += attr_im_info["x_offset"]
    y += attr_im_info["y_offset"]
    attr_im = attr_im_info["image"]
    im.paste(attr_im, (int(round(x)), int(round((y)))), attr_im)
    return im


def get_reprojected_image(tile_x, tile_y, level, map_reference, tile_size):
    tile_origin = tile_nw_corner(tile_x, tile_y, level)
    dest_pixel_size = tile_size_in_gmerc_meters(level) / tile_size
    maprecord = open_map_reference(map_reference)
    proj_transformer = pyproj.Transformer.from_crs(
        crs_gmerc, maprecord.crs, always_xy=True
    )

    def transform_dest_to_src_pixel(xxx_todo_changeme):
        (x, y) = xxx_todo_changeme
        x = x * dest_pixel_size + tile_origin[0]
        y = tile_origin[1] - y * dest_pixel_size
        x, y = proj_transformer.transform(x, y)
        x, y = maprecord.inv_gcp_transformer.transform(x, y)
        return x, y

    def transform_src_to_dest_pixel(x, y):
        x, y = maprecord.gcp_transformer.transform(x, y)
        x, y = proj_transformer.transform(
            x, y, direction=pyproj.transformer.TransformDirection.INVERSE
        )
        x = (x - tile_origin[0]) / dest_pixel_size
        y = (tile_origin[1] - y) / dest_pixel_size
        return x, y

    im_src = Image.open(maprecord.image_path)
    src_has_alpha = im_src.mode.endswith("A")
    im_src = im_src.convert("RGBA")
    if maprecord.mask_path is not None:
        im_mask = Image.open(maprecord.mask_path)
        im_src.putalpha(im_mask)
        src_has_alpha = True
    cell_size = 64
    mesh = []
    for cell_x in range(tile_size // cell_size):
        for cell_y in range(tile_size // cell_size):
            x1 = cell_x * cell_size
            y1 = cell_y * cell_size
            x2 = x1 + cell_size
            y2 = y1 + cell_size
            quad = (x1, y1), (x1, y2), (x2, y2), (x2, y1)
            quad = list(map(transform_dest_to_src_pixel, quad))
            quad = sum(quad, tuple())
            mesh.append(((x1, y1, x2, y2), quad))
    im = im_src.transform((tile_size, tile_size), Image.MESH, mesh, Image.BICUBIC)
    cutline_mask = Image.new("L", (tile_size, tile_size))
    cutline = maprecord.projected_cutline
    cutline = reproject_cutline_gmerc(maprecord.crs, cutline)
    cutline = [
        ((x - tile_origin[0]) / dest_pixel_size, (tile_origin[1] - y) / dest_pixel_size)
        for x, y in cutline
    ]
    draw = ImageDraw.Draw(cutline_mask)
    draw.polygon(cutline, fill=255, outline=255)
    mask = Image.new("L", (tile_size, tile_size))
    mask.paste(cutline_mask, (0, 0), im)
    if not src_has_alpha:
        im.putalpha(mask)
    im.paste(0, (0, 0), ImageChops.invert(mask))

    mid_point_y = tile_origin[1] - tile_size_in_gmerc_meters(level) / 2
    mid_point_lat = pyproj.transform(
        crs_gmerc, crs_gmerc.geodetic_crs, 0, mid_point_y, always_xy=True
    )[1]
    dest_meters_per_pixel = dest_pixel_size * math.cos(math.radians(mid_point_lat))
    im = apply_attribution(
        im, maprecord, transform_src_to_dest_pixel, dest_meters_per_pixel
    )
    return im


def slice_metatile(im, metatile_x, metatile_y, dest_level):
    meta_delta = dest_level - config.metatile_level
    meta_q = 2**meta_delta
    assert im.size == (256 * (meta_q),) * 2
    tile_x0 = metatile_x * meta_q
    tile_y0 = metatile_y * meta_q
    max_tile = 2**dest_level
    for d_tile_y in range(meta_q):
        y0 = d_tile_y * 256
        for d_tile_x in range(meta_q):
            x0 = d_tile_x * 256
            im2 = im.crop([x0, y0, x0 + 256, y0 + 256])
            tile_store.write(
                im2, (tile_x0 + d_tile_x) % max_tile, tile_y0 + d_tile_y, dest_level
            )


def render_tile(tile_x, tile_y, tile_level, map_references, tile_size):
    im = None
    for map_reference in map_references:
        im2 = get_reprojected_image(
            tile_x, tile_y, tile_level, map_reference, tile_size
        )
        if im is None:
            im = im2
        else:
            im.paste(im2, (0, 0), im2)
    return im


def process_metatile(xxx_todo_changeme1):
    # 1. render
    # 2. cut tiles
    # 3. shrink and cut overviews
    ((tile_x, tile_y, metatile_level), map_references) = xxx_todo_changeme1
    metatile_delta = config.max_level - metatile_level
    tile_size = 256 * (2**metatile_delta)
    im = render_tile(tile_x, tile_y, metatile_level, map_references, tile_size)
    for level in range(config.max_level, metatile_level, -1):
        slice_metatile(im, tile_x, tile_y, level)
        im = im.resize((im.size[0] // 2,) * 2, Image.ANTIALIAS)
    max_tile = 2**metatile_level
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
        for level in range(level + 1, config.max_level + 1):
            tx *= 2
            ty *= 2
            w *= 2
            for x in range(tx, tx + w):
                for y in range(ty, ty + w):
                    tile_store.remove(x, y, level)

    for tile in tiles_overviews:
        tile_store.remove(*tile)


def make_tiles_from_metalevel_to_maxlevel(tiles):
    n = 0
    if DEBUG:
        imap_func = map
    else:
        pool = Pool(cpu_count() + 1)
        imap_func = pool.imap_unordered
    for _ in imap_func(process_metatile, tiles):
        n += 1
        print(("\r%.1f%%" % (n * 100.0 / len(tiles))), end=" ")
        sys.stdout.flush()
    print()


def build_overview(xxx_todo_changeme2):
    (x, y, z) = xxx_todo_changeme2
    im = Image.new("RGBA", (512, 512))
    im2 = tile_store.open(x * 2, y * 2, z + 1)
    if im2 is not None:
        im.paste(im2.convert("RGBA"), (0, 0))
    im2 = tile_store.open(x * 2 + 1, y * 2, z + 1)
    if im2 is not None:
        im.paste(im2.convert("RGBA"), (256, 0))
    im2 = tile_store.open(x * 2, y * 2 + 1, z + 1)
    if im2 is not None:
        im.paste(im2.convert("RGBA"), (0, 256))
    im2 = tile_store.open(x * 2 + 1, y * 2 + 1, z + 1)
    if im2 is not None:
        im.paste(im2.convert("RGBA"), (256, 256))
    im = im.resize((256, 256), Image.ANTIALIAS)
    # FIXME: добавить закрашивание для уровней, которые делаются из metatile
    if config.highlight_level is not None and z <= config.highlight_level:
        mask = im.split()[-1]
        mask = mask.filter(ImageFilter.MaxFilter(5))
        im.paste(highlight_color, (0, 0), mask)
    tile_store.write(im, x, y, z)


def build_overviews(altered_tiles):
    need_update = set()
    for x, y, z in altered_tiles:
        if z > 0:
            need_update.add((x // 2, y // 2, z - 1))
    if DEBUG:
        imap_func = map
    else:
        pool = Pool(cpu_count() + 1)
        imap_func = pool.imap_unordered
    n = 0
    for _ in imap_func(build_overview, need_update):
        n += 1
        print("\r%.1f%%" % (n * 100.0 / len(need_update)), end=" ")
        sys.stdout.flush()
    print()
    if need_update:
        build_overviews(need_update)


def parse_image_format(s):
    fields = s.split(",")
    fmt = fields[0].upper()
    params = fields[1:]
    params = (p.split("=", 1) for p in params)
    params = dict((k, int(v)) for k, v in params)
    if fields[0] == "PNG8":
        if not (2 <= params.get("colors", None) <= 256):
            raise ValueError
        if not (1 <= params.get("speed", None) <= 10):
            raise ValueError
        if not (1 <= params.get("compression", 6) <= 9):
            raise ValueError
    elif fields[0] == "PNG32":
        if not (1 <= params.get("compression", 6) <= 9):
            raise ValueError
    elif fields[0] == "JPEG":
        if not (1 <= params.get("quality", None) <= 100):
            raise ValueError
    else:
        raise ValueError('Unsupported image format "%s"' % fields[0])
    return fmt, params


class MyArgumentParser(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_line):
        arg_line = arg_line.strip()
        if arg_line and arg_line[0] == "#":
            return []
        if "=" in arg_line and not arg_line.startswith("-"):
            arg_line = "--" + arg_line
        if arg_line:
            return [arg_line]
        else:
            return []


def parse_command_line():
    # TODO: сделать возможны включение файлов при вызове не из текущей директории
    # TODO: по умолчанию включать файл tiles.cfg из текущего каталога
    # TODO: добавить описание форматов изображений
    parser = MyArgumentParser(fromfile_prefix_chars="@")
    parser.add_argument("--storage-format", choices=["files", "mbtiles"], required=True)
    parser.add_argument("--max-level", type=int, required=True)
    parser.add_argument("--highlight-level", type=int)
    parser.add_argument("maps", metavar="FILE", type=str, nargs="+")
    parser.add_argument("--image-format", type=parse_image_format, required=True)
    parser.add_argument("--image-border-format", type=parse_image_format)
    parser.add_argument(
        "--out",
        metavar="PATH",
        dest="out_path",
        required=True,
        help="Filename of mbtiles container or tiles dir",
    )
    parser.add_argument(
        "-f", dest="do_update", action="store_true", help="actualy do update"
    )
    return parser.parse_args()


def print_job_info(job):
    print("Maps added: %s" % job.get_added_maps_n())
    print("Maps removed: %s" % job.get_removed_maps_n())
    print("Metatiles to rebuild: %s" % len(job.get_tiles_needing_update()))


def configure_output_storage():
    global tile_store
    image_encoder = image_store.get_image_encoder(
        config.image_format, config.image_border_format
    )
    tile_store_class = {
        "files": image_store.FilesWriter,
        "mbtiles": image_store.MBTilesWriter,
    }[config.storage_format]
    tile_store = tile_store_class(config.out_path, image_encoder)


def main():
    warnings.filterwarnings(action="ignore", category=FutureWarning)
    global config
    config = parse_command_line()
    config.metatile_level = max(config.max_level - METATILE_DELTA, 0)
    configure_output_storage()
    prev_state = tile_store.read_metadata(PREV_STATE_FILENAME)
    job = JobManager(config.maps, config.metatile_level, prev_state)
    print_job_info(job)
    if config.do_update:
        t = time.time()
        remove_tiles(job.get_removed_tiles())
        make_tiles_from_metalevel_to_maxlevel(job.get_tiles_needing_update())
        print("Building overviews")
        overview_tiles_needing_update = (
            tile for (tile, _) in job.get_tiles_needing_update()
        )
        overview_tiles_needing_update = chain(
            overview_tiles_needing_update, job.get_removed_tiles()
        )
        overview_tiles_needing_update = (
            (x % (2**z), y, z) for (x, y, z) in overview_tiles_needing_update
        )
        build_overviews(overview_tiles_needing_update)
        job.apply_update()
        print("Done in %.1f seconds" % (time.time() - t))
    else:
        print('To actualy update mbtiles db, specify "-f" option')


if __name__ == "__main__":
    try:
        main()
    except UserInputError as e:
        print(e)
        exit(1)
