"""Repair unreadable ASCII subnormals in copies of native solution fields."""
import hashlib
from pathlib import Path
import re
import sys

SCIENTIFIC = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)[eE][+-]?\d+")


def normalize_ascii_fields(directory):
    changes = []
    for path in sorted(Path(directory).iterdir()):
        if not path.is_file():
            continue
        data = path.read_bytes()
        header = data[:4096].decode("ascii", errors="ignore")
        if not re.search(r"class\s+(?:vol|surface)\w+Field\s*;", header):
            continue
        if re.search(r"format\s+binary\s*;", header):
            continue
        text = data.decode("ascii")
        magnitudes = []
        def replace(match):
            token = match.group()
            value = float(token)
            mantissa = token.lower().split('e')[0]
            if abs(value) < sys.float_info.min and any(c in '123456789' for c in mantissa):
                magnitudes.append(abs(value))
                return "0"
            return token
        normalized = SCIENTIFIC.sub(replace, text)
        if magnitudes:
            path.write_text(normalized, encoding="ascii")
            changes.append({"field":path.name,"count":len(magnitudes),
                            "maximum_absolute_replacement":max(magnitudes),
                            "original_sha256":hashlib.sha256(data).hexdigest(),
                            "normalized_sha256":hashlib.sha256(path.read_bytes()).hexdigest()})
    return {"rule":"ASCII subnormal values below IEEE-754 minimum normal are zeroed in a copy",
            "absolute_threshold":sys.float_info.min,"changes":changes}
