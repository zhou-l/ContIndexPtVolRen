"""NRRD / NHDR reading and writing.

Python port of ``nhdr_nrrd_read.m``, ``nrrd_getMatlabDataType.m``,
``nrrd_getSpaceDimensions.m`` and ``nrrdWriter.m``.

Arrays are kept in MATLAB (column-major / Fortran) axis convention: the first
array axis is the fastest-varying axis in the file, exactly as the MATLAB
reader returns them.
"""
import gzip
import os
import re
import sys

import numpy as np

_NRRD_TO_NUMPY = {
    ('signed char', 'int8', 'int8_t'): np.int8,
    ('uchar', 'unsigned char', 'uint8', 'uint8_t'): np.uint8,
    ('short', 'short int', 'signed short', 'signed short int', 'int16', 'int16_t'): np.int16,
    ('ushort', 'unsigned short', 'unsigned short int', 'uint16', 'uint16_t'): np.uint16,
    ('int', 'signed int', 'int32', 'int32_t'): np.int32,
    ('uint', 'unsigned int', 'uint32', 'uint32_t'): np.uint32,
    ('longlong', 'long long', 'long long int', 'signed long long',
     'signed long long int', 'int64', 'int64_t'): np.int64,
    ('ulonglong', 'unsigned long long', 'unsigned long long int', 'uint64', 'uint64_t'): np.uint64,
    ('float',): np.float32,
    ('double',): np.float64,
}

_NUMPY_TO_NRRD = {
    np.dtype(np.int8): 'int8', np.dtype(np.uint8): 'uint8',
    np.dtype(np.int16): 'int16', np.dtype(np.uint16): 'uint16',
    np.dtype(np.int32): 'int32', np.dtype(np.uint32): 'uint32',
    np.dtype(np.int64): 'int64', np.dtype(np.uint64): 'uint64',
    np.dtype(np.float32): 'float', np.dtype(np.float64): 'double',
}


def nrrd_get_numpy_dtype(meta_type):
    """Map an NRRD type string to a numpy dtype (``nrrd_getMatlabDataType``)."""
    if meta_type == 'block':
        raise ValueError("Data type 'block' is currently not supported.")
    for names, dt in _NRRD_TO_NUMPY.items():
        if meta_type in names:
            return np.dtype(dt)
    raise ValueError('Unknown datatype: %s' % meta_type)


def nrrd_get_space_dimensions(space_descriptor):
    """Number of world-space dimensions for a ``space`` descriptor (-1 if unknown)."""
    d = space_descriptor.lower()
    three = ['right-anterior-superior', 'ras', 'left-anterior-superior', 'las',
             'left-posterior-superior', 'lps', 'scanner-xyz', '3d-right-handed', '3d-left-handed']
    four = ['right-anterior-superior-time', 'rast', 'left-anterior-superior-time', 'last',
            'left-posterior-superior-time', 'lpst', 'scanner-xyz-time',
            '3d-right-handed-time', '3d-left-handed-time']
    if d in three:
        return 3
    if d in four:
        return 4
    return -1


def _extract_numbers(s, tokens=('(', ')', ',', '"', "'")):
    for t in tokens:
        s = s.replace(t, ' ')
    return [float(v) for v in s.split()]


def _parse_datafiles(field_value):
    """Return (filelist, list_mode, subdim) for a 'data file:' field."""
    field_value = field_value.strip()
    if field_value[:4].upper() == 'LIST':
        rest = field_value[4:].split()
        return [], True, (int(rest[0]) if rest else None)
    parts = field_value.split()
    if len(parts) == 1:
        return [parts[0]], False, None
    if len(parts) in (4, 5):
        fmt = parts[0]
        id_min, id_max, step = int(parts[1]), int(parts[2]), int(parts[3])
        if step == 0:
            raise ValueError('detached data files: step should not be zero.')
        stop = id_max + (1 if step > 0 else -1)
        files = [fmt % i for i in range(id_min, stop, step)]
        subdim = int(parts[4]) if len(parts) == 5 else None
        return files, False, subdim
    raise ValueError("Invalid 'data file' field: '%s'" % field_value)


def _read_payload(raw_bytes, n_elems, encoding, dtype, endian):
    enc = encoding.lower()
    if enc in ('gzip', 'gz'):
        raw_bytes = gzip.decompress(raw_bytes)
        enc = 'raw'
    if enc == 'raw':
        n_bytes = n_elems * dtype.itemsize
        # data is the trailing part of the file (byte skip -1 semantics)
        raw_bytes = raw_bytes[len(raw_bytes) - n_bytes:]
        bo = '<' if endian.lower() == 'little' else '>'
        data = np.frombuffer(raw_bytes, dtype=dtype.newbyteorder(bo)).astype(dtype)
    elif enc in ('txt', 'text', 'ascii'):
        data = np.array(raw_bytes.decode('ascii').split(), dtype=np.float64).astype(dtype)
    else:
        raise ValueError('Unsupported encoding: %s' % encoding)
    if data.size != n_elems:
        raise ValueError('Detected %d elements instead of %d announced in header.'
                         % (data.size, n_elems))
    return data


def nhdr_nrrd_read(filename, read_data=True):
    """Read an NRRD file or detached NHDR header.

    Returns a dict of header fields; when ``read_data`` is True the array is
    stored under ``'data'`` with shape ``sizes`` in Fortran order (same as
    MATLAB).
    """
    main_dir = os.path.dirname(filename)
    base, ext = os.path.splitext(os.path.basename(filename))
    if ext.lower() not in ('.nhdr', '.nrrd'):
        print('WARNING: %s looks like a %s file, not a nhdr or nrrd file.' % (filename, ext),
              file=sys.stderr)
    info = {'content': base, 'data': None}

    with open(filename, 'rb') as f:
        content = f.read()

    # split header / data at the first empty line
    pos = 0
    lines = []
    while pos < len(content):
        nl = content.find(b'\n', pos)
        if nl == -1:
            nl = len(content)
        line = content[pos:nl].decode('latin-1').rstrip('\r')
        pos = nl + 1
        if not lines:
            if not line.startswith('NRRD') or len(line) < 8:
                raise ValueError('Bad signature. First line should be magic line "NRRD000X".')
            if int(line[4:]) > 5:
                raise ValueError('Only NRRD versions up to 5 are supported.')
            lines.append(line)
            continue
        if line == '':
            break
        lines.append(line)
    payload = content[pos:]

    list_mode = False
    subdim = None
    i = 1
    while i < len(lines):
        cs = lines[i]
        i += 1
        if cs.lstrip().startswith('#'):
            continue
        if ':=' in cs:
            key, val = cs.split(':=', 1)
            ku = key.strip().upper()
            if ku == 'MODALITY':
                info['modality'] = val.strip()
            elif ku == 'DWMRI_B-VALUE':
                info['bvalue'] = float(val)
            elif ku.startswith('DWMRI_GRADIENT_'):
                info.setdefault('gradients', {})[int(ku[len('DWMRI_GRADIENT_'):])] = \
                    [float(v) for v in val.split()]
            else:
                info.setdefault('keyvalues', {})[key.strip()] = val
            continue
        if ':' not in cs:
            print("WARNING nhdr_nrrd_read: Could not parse input line: '%s'" % cs, file=sys.stderr)
            continue
        key, val = cs.split(':', 1)
        key = key.strip().lower()
        val = val.strip()
        if key == 'content':
            info['content'] = val
        elif key == 'type':
            info['type'] = val
        elif key == 'endian':
            info['endian'] = val
        elif key == 'encoding':
            info['encoding'] = val
        elif key == 'dimension':
            info['dimension'] = int(val)
        elif key == 'sizes':
            info['sizes'] = [int(v) for v in val.split()]
        elif key == 'kinds':
            info['kinds'] = val.split()
        elif key == 'space':
            info['space'] = val
            if nrrd_get_space_dimensions(val) <= 0:
                raise ValueError("unrecognized 'space' descriptor '%s'" % val)
        elif key == 'space dimension':
            info['spacedimension'] = int(val)
        elif key == 'space directions':
            info['spacedirections'] = val.split()
            nums = _extract_numbers(val.replace('none', ''))
            sd = int(round(np.sqrt(len(nums))))
            if sd * sd == len(nums) and sd > 0:
                info['spacedirections_matrix'] = np.array(nums).reshape((sd, sd)).T
        elif key == 'space units':
            info['spaceunits'] = [u.strip() for u in val.replace('none', '').split('"')
                                  if u.strip()]
        elif key == 'space origin':
            info['spaceorigin'] = np.array(_extract_numbers(val))
        elif key == 'measurement frame':
            nums = _extract_numbers(val.replace('none', ''))
            sd = int(round(np.sqrt(len(nums))))
            info['measurementframe'] = np.array(nums).reshape((sd, sd)).T
        elif key == 'thicknesses':
            info['thicknesses'] = [float(v) for v in val.split()]
        elif key == 'centerings':
            info['centerings'] = val.split()
        elif key in ('line skip', 'lineskip'):
            info['lineskip'] = int(val)
        elif key in ('byte skip', 'byteskip'):
            info['byteskip'] = int(val)
        elif key in ('data file', 'datafile'):
            files, list_mode, subdim = _parse_datafiles(val)
            if list_mode:
                files = [ln for ln in lines[i:] if ln]
                i = len(lines)
            info['datafiles'] = files
        else:
            print("WARNING nhdr_nrrd_read: Could not parse input line: '%s'" % cs, file=sys.stderr)

    for req in ('sizes', 'dimension', 'encoding', 'type'):
        if req not in info:
            raise ValueError("Missing required '%s' field in header" % req)
    dtype = nrrd_get_numpy_dtype(info['type'])
    info.setdefault('endian', 'little')

    if not read_data:
        return info

    sizes = info['sizes']
    n_tot = int(np.prod(sizes))
    if 'datafiles' in info:
        if subdim is not None and subdim != info['dimension']:
            raise ValueError('Reading slices along an axis other than the last is not supported.')
        files = info['datafiles']
        if n_tot % len(files) != 0:
            raise ValueError('Number of data files does not divide the number of values.')
        n_per = n_tot // len(files)
        flat = np.zeros(n_tot, dtype=dtype)
        for k, df in enumerate(files):
            path = df if os.path.isabs(df) else os.path.join(main_dir, df)
            if path.lower().endswith('.nhdr'):
                raise ValueError('nhdr file should not be used as detached data file.')
            if path.lower().endswith('.nrrd'):
                sub = nhdr_nrrd_read(path, True)
                flat[k * n_per:(k + 1) * n_per] = sub['data'].ravel(order='F')
            else:
                with open(path, 'rb') as fd:
                    flat[k * n_per:(k + 1) * n_per] = _read_payload(
                        fd.read(), n_per, info['encoding'], dtype, info['endian'])
        data = flat
    else:
        data = _read_payload(payload, n_tot, info['encoding'], dtype, info['endian'])
    info['data'] = data.reshape(sizes, order='F')
    return info


def _fmt_num(v):
    """Rough equivalent of MATLAB num2str for scalars."""
    v = float(v)
    return str(int(v)) if v.is_integer() else '%.4g' % v


def nrrd_writer(filename, matrix, pixelspacing=(1, 1, 1), origin=(0, 0, 0), encoding='raw'):
    """Write a 2-D/3-D array to .nrrd or .nhdr (port of ``nrrdWriter.m``).

    Like the MATLAB version, the first two axes are swapped before writing.
    """
    matrix = np.asarray(matrix)
    base, ext = os.path.splitext(filename)
    fmt = ext[1:].lower()
    encoding = encoding.lower()
    assert encoding in ('ascii', 'raw', 'gzip'), 'Unsupported encoding'
    assert fmt in ('nhdr', 'nrrd'), 'Unexpected format'

    if matrix.ndim >= 2:
        axes = list(range(matrix.ndim))
        axes[0], axes[1] = 1, 0
        matrix = np.transpose(matrix, axes)  # undo permute of index in nrrd reader
    dims = matrix.shape
    nd = len(dims)
    if matrix.dtype not in _NUMPY_TO_NRRD:
        raise ValueError('Unknown datatype')
    outtype = _NUMPY_TO_NRRD[matrix.dtype]

    hdr = ['NRRD0004', 'type: ' + outtype, 'dimension: %d' % nd]
    if nd == 2:
        hdr.append('space: left-posterior')
    elif nd == 3:
        hdr.append('space: left-posterior-superior')
    hdr.append('sizes: ' + ' '.join(str(d) for d in dims))
    ps = [_fmt_num(p) for p in pixelspacing]
    if nd == 2:
        hdr.append('space directions: (%s,0) (0,%s)' % (ps[0], ps[1]))
        hdr.append('kinds: domain domain')
    elif nd == 3:
        hdr.append('space directions: (%s,0,0) (0,%s,0) (0,0,%s)' % (ps[0], ps[1], ps[2]))
        hdr.append('kinds: domain domain domain')
    hdr.append('encoding: ' + encoding)
    hdr.append('endian: ' + ('big' if sys.byteorder == 'big' else 'little'))
    org = [_fmt_num(o) for o in origin]
    if nd == 2:
        hdr.append('space origin: (%s,%s)' % (org[0], org[1]))
    elif nd == 3:
        hdr.append('space origin: (%s,%s,%s)' % (org[0], org[1], org[2]))

    flat = matrix.ravel(order='F')
    if encoding == 'raw':
        body = flat.tobytes()
    elif encoding == 'gzip':
        body = gzip.compress(flat.tobytes())
    else:
        body = (' '.join('%u' % v if np.issubdtype(flat.dtype, np.integer) else '%g' % v
                         for v in flat) + ' ').encode('ascii')

    if fmt == 'nhdr':
        data_name = os.path.basename(base) + '.' + encoding
        hdr.append('data file: ' + data_name)
        with open(filename, 'wb') as f:
            f.write(('\n'.join(hdr) + '\n').encode('ascii'))
        with open(os.path.join(os.path.dirname(filename), data_name), 'wb') as f:
            f.write(body)
    else:
        with open(filename, 'wb') as f:
            f.write(('\n'.join(hdr) + '\n\n').encode('ascii'))
            f.write(body)
    return flat.size
