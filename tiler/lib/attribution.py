# coding: utf-8
import math

import cairo
import gi
from PIL import Image

gi.require_version("PangoCairo", "1.0")
gi.require_version("Pango", "1.0")
from gi.repository import Pango, PangoCairo


# FONTNAME = 'Droid Sans'
FONTNAME = "Noto Sans"
# FONTNAME = 'Sans'


_font_installed = None


def check_font_installed():
    global _font_installed
    if _font_installed is None:
        for font in PangoCairo.font_map_get_default().list_families():
            if font.get_name() == FONTNAME:
                _font_installed = True
                break
        else:
            _font_installed = False
    if not _font_installed:
        raise Exception('Font "%s" not installed' % FONTNAME)


def create_path(ctx, text, font_size_px, rotate):
    check_font_installed()
    ctx.set_line_width(font_size_px / 5.0)
    ctx.set_line_join(cairo.LINE_JOIN_ROUND)
    ctx.save()

    ctx.rotate(rotate)
    ctx.set_antialias(cairo.ANTIALIAS_SUBPIXEL)
    layout = PangoCairo.create_layout(ctx)
    font = Pango.FontDescription(FONTNAME)
    font.set_absolute_size(font_size_px * Pango.SCALE)
    layout.set_alignment(Pango.Alignment.LEFT)
    layout.set_font_description(font)
    layout.set_line_spacing(0.8)
    layout.set_text(text)

    PangoCairo.layout_path(ctx, layout)
    ctx.restore()


def draw(text, font_size_px, rotate):
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, 500, 500)
    ctx = cairo.Context(surface)
    create_path(ctx, text, font_size_px, rotate)
    x1, y1, x2, y2 = ctx.stroke_extents()
    x_offset, y_offset = ctx.fill_extents()[:2]

    width = int(math.ceil(x2 - x1))
    height = int(math.ceil(y2 - y1))
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
    ctx = cairo.Context(surface)
    ctx.move_to(-x1, -y1)
    create_path(ctx, text, font_size_px, rotate)

    ctx.set_source_rgb(0, 0, 0)
    ctx.stroke_preserve()
    ctx.set_source_rgb(255, 255, 255)
    ctx.fill_preserve()
    im = Image.frombuffer(
        "RGBA", (width, height), surface.get_data().tobytes(), "raw", "BGRA", 0, 1
    )
    return {"image": im, "x_offset": x_offset, "y_offset": y_offset}


def make_attribution_text(image_info):
    text = []
    date = image_info.get("date")
    if date:
        text.append(image_info["date"].split("-")[0])
    scale = image_info.get("scale")
    if scale:
        s = "1:{:,}".format(int(image_info["scale"])).replace(",", " ")
        h = image_info.get("h")
        if h:
            s += "  h=%s м" % image_info["h"]
        text.append(s)
    event = image_info.get("event")
    if event:
        text.append(event)
    text = "\n".join(text)
    return text


src_font_size_mm = 3


def make_attribution_image(
    attrib_data, image_info, src_pixel_to_dest_pixel_transformer, dest_meters_in_pixel
):
    text = make_attribution_text(image_info)
    font_size_ground_meters = (
        float(src_font_size_mm) * attrib_data["scale_denom"] / 1000
    )
    dest_font_size = font_size_ground_meters / dest_meters_in_pixel

    xy_src_1 = attrib_data["anchor"]
    xy_src_2 = {
        "x": xy_src_1["x"] + math.cos(attrib_data["angle"]) * 100,
        "y": xy_src_1["y"] + math.sin(attrib_data["angle"]) * 100,
    }

    xy_dest_1 = src_pixel_to_dest_pixel_transformer(xy_src_1["x"], xy_src_1["y"])
    xy_dest_2 = src_pixel_to_dest_pixel_transformer(xy_src_2["x"], xy_src_2["y"])

    dest_angle = math.atan2(xy_dest_2[1] - xy_dest_1[1], xy_dest_2[0] - xy_dest_1[0])
    # print 'image angle', math.degrees(attrib_data['angle']), 'proj anfle', math.degrees(dest_angle)
    return draw(text, dest_font_size, dest_angle)


def get_text_size(attrib_data, image_info, dest_meters_in_pixel):
    text = make_attribution_text(image_info)
    font_size_ground_meters = (
        float(src_font_size_mm) * attrib_data["scale_denom"] / 1000
    )
    font_size_px = font_size_ground_meters / dest_meters_in_pixel
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, 500, 500)
    ctx = cairo.Context(surface)
    create_path(ctx, text, font_size_px, rotate=0)
    x1, y1, x2, y2 = ctx.stroke_extents()

    width = int(math.ceil(x2 - x1))
    height = int(math.ceil(y2 - y1))
    return width, height


def get_attrib_path(image_path):
    return image_path + ".attrib.json"


def get_info_path(image_path):
    return image_path + ".info.json"


if __name__ == "__main__":
    # im = draw('1:25 000 Жаркий июль', 1, 0 )
    # print im.size
    # im.show()
    info = {"scale": "25000", "h": 10, "date": "2014-1-2"}
    transformer = lambda x, y: (x, y)
    data = {
        "anchor": {"x": 542, "y": 333},
        "scale_denom": 7500,
        "angle": 0.0,
    }

    make_attribution_image(data, info, transformer, 1).show()
