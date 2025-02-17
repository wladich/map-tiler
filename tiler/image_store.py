# -*- coding: utf-8 -*-
import multiprocessing
import os
import sqlite3 as sqlite
from array import array
from functools import partial
from io import BytesIO

import imagequant
import png
from PIL import Image

db_lock = multiprocessing.Lock()


def open_image(s):
    if s.startswith(b"\x89PNG"):
        r = png.Reader(bytes=s)
        w, h, pixels, meta = r.asRGBA()
        ar = array("B")
        for row in list(pixels):
            ar.extend(row)
        data = ar.tobytes()
        im = Image.frombytes("RGBA", (w, h), data)
        return im
    else:
        return Image.open(BytesIO(s))


def validate_tile_index(tile_x, tile_y, level):
    assert isinstance(tile_x, int)
    assert isinstance(tile_y, int)
    assert isinstance(level, int)
    if level < 0:
        raise ValueError(f"Negative value of level: {tile_x=} {tile_y=} {level=}")
    err = []
    max_coord = 2**level - 1
    if tile_x > max_coord or tile_x < 0:
        err = [f"Invalid value of tile_x: {tile_x=} {tile_y=} {level=}"]
    if tile_y > max_coord or tile_y < 0:
        err.append(f"Invalid value of tile_y: {tile_x=} {tile_y=} {level=}")
    if err:
        raise ValueError(";".join(err))


class MBTilesWriter(object):
    SCHEME = """
        CREATE TABLE tiles(
            zoom_level integer, tile_column integer, tile_row integer, tile_data blob,
            UNIQUE(zoom_level, tile_column, tile_row) ON CONFLICT REPLACE);
       
        CREATE TABLE metadata (name text, value text, UNIQUE(name) ON CONFLICT REPLACE);
    """

    PRAGMAS = """
        PRAGMA journal_mode = WAL;
        PRAGMA synchronous = 0;
        PRAGMA busy_timeout = 10000;
    """

    def __init__(self, path, image_encoder):
        need_init = not os.path.exists(path)
        self.encoder = image_encoder
        self.path = path
        if need_init:
            with db_lock, self.conn() as conn:
                conn.executescript(self.SCHEME)

    def conn(self):
        conn = sqlite.connect(self.path)
        conn.executescript(self.PRAGMAS)
        return conn

    def write(self, im, tile_x, tile_y, level):
        validate_tile_index(tile_x, tile_y, level)
        encoder = self.encoder(im)
        image_not_empty = next(encoder)
        if image_not_empty:
            tile_y = 2**level - tile_y - 1
            s = BytesIO()
            encoder.send(s)
            with db_lock, self.conn() as conn:
                conn.execute(
                    """
                    INSERT INTO tiles (zoom_level, tile_column, tile_row, tile_data) VALUES (?,?,?,?)""",
                    (level, tile_x, tile_y, s.getvalue()),
                )
        else:
            self.remove(tile_x, tile_y, level)

    def remove(self, tile_x, tile_y, level):
        validate_tile_index(tile_x, tile_y, level)
        tile_y = 2**level - tile_y - 1
        with db_lock, self.conn() as conn:
            conn.execute(
                "DELETE FROM tiles WHERE zoom_level=? AND tile_column=? AND tile_row=?",
                (level, tile_x, tile_y),
            )

    def open(self, tile_x, tile_y, level):
        validate_tile_index(tile_x, tile_y, level)
        tile_y = 2**level - tile_y - 1
        with db_lock, self.conn() as conn:
            row = conn.execute(
                "SELECT tile_data FROM tiles WHERE zoom_level=? AND tile_column=? AND tile_row=?",
                (level, tile_x, tile_y),
            ).fetchone()
        if row:
            return open_image(row[0])

    def write_metadata(self, key, value):
        with db_lock, self.conn() as conn:
            conn.execute(
                "INSERT INTO metadata (name, value) VALUES (?, ?)", (key, value)
            )

    def read_metadata(self, key):
        with db_lock, self.conn() as conn:
            row = conn.execute(
                "SELECT value FROM metadata WHERE name=?", (key,)
            ).fetchone()
        if row:
            return row[0]

    def close(self):
        #        self.conn.commit()
        with db_lock, self.conn() as conn:
            conn.execute("PRAGMA journal_mode = off")

    def __del__(self):
        self.close()


class FilesWriter(object):
    def __init__(self, path, image_encoder):
        if not os.path.isdir(path):
            os.makedirs(path)
        self.path = path
        self.encoder = image_encoder

    def _get_tile_file_name(self, tile_x, tile_y, level):
        filename = "%s_%s_%s" % (level, tile_y, tile_x)
        filename = os.path.join(self.path, filename)
        return filename

    def write(self, im, tile_x, tile_y, level):
        validate_tile_index(tile_x, tile_y, level)
        encoder = self.encoder(im)
        image_not_empty = next(encoder)
        if image_not_empty:
            filename = self._get_tile_file_name(tile_x, tile_y, level)
            with open(filename, "wb") as f:
                encoder.send(f)
        else:
            self.remove(tile_x, tile_y, level)

    def remove(self, tile_x, tile_y, level):
        validate_tile_index(tile_x, tile_y, level)
        filename = self._get_tile_file_name(tile_x, tile_y, level)
        if os.path.exists(filename):
            os.remove(filename)

    def open(self, tile_x, tile_y, level):
        validate_tile_index(tile_x, tile_y, level)
        filename = self._get_tile_file_name(tile_x, tile_y, level)
        if os.path.exists(filename):
            with open(filename, 'rb') as f:
                return open_image(f.read())

    def _get_metadata_file_name(self, key):
        filename = "_meta_%s" % key
        filename = os.path.join(self.path, filename)
        return filename

    def write_metadata(self, key, value):
        with open(self._get_metadata_file_name(key), "w") as f:
            f.write(value)

    def read_metadata(self, key):
        filename = self._get_metadata_file_name(key)
        if os.path.exists(filename):
            with open(filename) as f:
                return f.read()

    def close(self):
        pass


def save_png_rgba(im, fd, compression=None):
    has_alpha = im.mode[-1] == "A"
    ar = array("B", im.tostring())
    pngw = png.Writer(size=im.size, alpha=has_alpha, compression=compression)
    pngw.write_array(fd, ar)


def save_png_with_palette(im, fd, colors, speed, compression=None):
    quantized = imagequant.quantize_image(im, colors=colors, speed=speed)
    palette = list(quantized["palette"])
    palette = list(zip(palette[::4], palette[1::4], palette[2::4], palette[3::4]))
    if all(c[3] == 255 for c in palette):
        palette = [c[:3] for c in palette]
    pngw = png.Writer(size=im.size, palette=palette, compression=compression)
    pngw.write_array(fd, quantized["image"])


def save_jpeg(im, fd, quality):
    if im.mode[-1] == "A":
        im = im.convert("RGB")
    im.save(fd, "JPEG", quality=quality)


def get_image_type(im):
    if im.mode[-1] == "A":
        alpha_min, alpha_max = im.split()[-1].getextrema()
        if alpha_min == 255:
            return "full"
        if alpha_max == 0:
            return "empty"
        return "border"
    else:
        "full"


def get_image_encoder(image_format, image_border_format):
    encoder_functions = {
        "PNG8": save_png_with_palette,
        "PNG32": save_png_rgba,
        "JPEG": save_jpeg,
    }
    image_encoder = partial(encoder_functions[image_format[0]], **image_format[1])
    if image_border_format is None:
        image_border_encoder = image_encoder
    else:
        image_border_encoder = partial(
            encoder_functions[image_border_format[0]], **image_border_format[1]
        )

    def _image_encoder(im):
        image_type = get_image_type(im)
        fd = yield image_type != "empty"
        if image_type == "full":
            image_encoder(im, fd)
        elif image_type == "border":
            image_border_encoder(im, fd)
        else:
            pass
        yield

    return _image_encoder
