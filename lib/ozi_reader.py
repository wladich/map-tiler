# -*- coding: utf-8 -*-

from validate  import validate_number, validate_float, validate_notempty,\
    validate_value, validate_values, ValidationError, validate_string_start

from copy import deepcopy
class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        #self.__dict__ = self


    def __getattr__(self, name):
        return self[name]

    def __setattr__(self, name, value):
        self[name] = value

    def __deepcopy__(self, memo):
        return AttrDict(deepcopy(dict(self), memo))
        
class OziFormatError(Exception):
    def __init__(self, message):
        self.messages = [message]

    def __str__(self):
        return ': '.join(['Error in Ozi map file'] + self.messages)

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type == ValidationError:
            self.messages.append(str(exc_value))
            raise self
        elif exc_type == OziFormatError:
            self.messages.extend(exc_value.messages)
            raise self

def fields(line, count):
    fields = [s.strip() for s in line.split(',')]
    if len(fields) < count:
        fields += [''] * (count - len(fields))
    return fields

def read_map(lines):
    lines = [l.strip(' \n\r\x09') for l in lines]
    if not lines:
        raise OziFormatError('Document empty.')
    ozi_map = AttrDict()
    try:
        with OziFormatError('line 1'):
            validate_string_start(lines[0], 'OziExplorer Map Data File Version 2.')
        ozi_map.title = lines[1].decode('cp1251')
        ozi_map.file_name = lines[2].decode('cp1251').split('\\')[-1]
        with OziFormatError('line 5, datum name'):
            ozi_map.datum = validate_notempty(fields(lines[4], 1)[0])
        proj_params = fields(lines[8], 2)
        with OziFormatError('line 9'):
            validate_value(proj_params[0], 'Map Projection')
        ozi_map.projection = AttrDict()
        with OziFormatError('line 9, projection name'):
            ozi_map.projection.name = validate_notempty(proj_params[1])
        #proj_params = [v if v else '0' for v in fields(lines[39], 10)]
        proj_params =  fields(lines[39], 10)
        with OziFormatError('line 40'):
            validate_value(proj_params[0], 'Projection Setup')
            if proj_params[1] != '': ozi_map.projection.lat_origin = validate_float(proj_params[1])
            if proj_params[2] != '': ozi_map.projection.lon_origin = validate_float(proj_params[2])
            if proj_params[3] != '': ozi_map.projection.k_factor = validate_float(proj_params[3])
            if proj_params[4] != '': ozi_map.projection.false_easting = validate_float(proj_params[4])
            if proj_params[5] != '': ozi_map.projection.false_northing = validate_float(proj_params[5])
            if proj_params[6] != '': ozi_map.projection.lat1 = validate_float(proj_params[6])
            if proj_params[7] != '': ozi_map.projection.lat2 = validate_float(proj_params[7])
            if proj_params[8] != '': ozi_map.projection.height = validate_float(proj_params[8])

        ozi_map.gcps = []
        for i in xrange(1, 31):
            ozi_gcp = fields(lines[i + 8], 17)
            with OziFormatError('line %d' % (i + 9)):
                validate_value(ozi_gcp[0], 'Point%02d' % i)
                validate_value(ozi_gcp[1], 'xy')
                validate_value(ozi_gcp[5], 'deg')
                validate_value(ozi_gcp[12], 'grid')
                validate_values(ozi_gcp[4], ['in', 'ex'])
                if ozi_gcp[4] == 'in' and ozi_gcp[2] and ozi_gcp[3]:
                    gcp = AttrDict(pixel = AttrDict(), ref = AttrDict())
                    gcp.pixel.x = validate_number(ozi_gcp[2])
                    gcp.pixel.y = validate_number(ozi_gcp[3])
                    if ozi_gcp[6] and ozi_gcp [7] and ozi_gcp[9] and ozi_gcp[10]:
                        validate_values(ozi_gcp[8], ['S', 'N'])
                        validate_values(ozi_gcp[11], ['W', 'E'])
                        gcp.type = 'latlon'
                        gcp.ref.x = validate_number(ozi_gcp[9]) + validate_float(ozi_gcp[10]) / 60
                        gcp.ref.y = validate_number(ozi_gcp[6]) + validate_float(ozi_gcp[7]) / 60
                        if ozi_gcp[11] == 'W':
                            gcp.ref.x *= -1
                        if ozi_gcp[8] == 'S':
                            gcp.ref.y *= -1
                    elif ozi_gcp[14] and ozi_gcp[15]:
                        validate_values(ozi_gcp[16], ['N', 'S'])
                        gcp.type = 'proj'
                        gcp.ref.x = validate_float(ozi_gcp[14])
                        gcp.ref.y = validate_float(ozi_gcp[15])
                        if ozi_gcp[16] == 'S':
                            gcp.ref.y *= -1
                        if ozi_gcp[13]:
                            gcp.zone = validate_number(ozi_gcp[13])
                    else:
                        raise OziFormatError('incomplete gcp definition')
                    ozi_map.gcps.append(gcp)

        ozi_map.cutline = []
        ozi_map.cutline_pixels = []
        for line in lines[40:]:
            point = fields(line, 4)
            if point[0] == 'MMPLL':
                lat = validate_float(point[3])
                lon = validate_float(point[2])
                ozi_map.cutline.append((lon, lat))
            elif point[0] == 'MMPXY':
                x = validate_number(point[2])
                y = validate_number(point[3])
                ozi_map.cutline_pixels.append((x, y))
    except IndexError:
        raise OziFormatError('Document too short.')
    return ozi_map

if __name__ == '__main__':
    import sys
    import pprint
    print sys.argv[1]
    pprint.pprint(read_map(open(sys.argv[1]).readlines()))
    