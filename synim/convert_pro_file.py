"""
Stand-alone tool to convert IDL configuration files (.pro)
to the standard YAML format used by SPECULA and SynIM.
"""

import re
import yaml
import argparse
from pathlib import Path

def parse_pro_file(pro_file_path):
    """Parse a .pro file and extract its structure into a Python dictionary."""
    data = {}
    current_section = None

    with open(pro_file_path, 'r') as file:
        for line_num, line in enumerate(file, 1):
            # Remove comments (everything after ';')
            comment_pos = line.find(';')
            if comment_pos != -1:
                line = line[:comment_pos]

            line = line.strip()
            if not line:
                continue

            try:
                # Detect the beginning of a new section (e.g., {main, {dm1, etc.)
                section_match = re.match(r'^\{\s*(\w+)\s*,?', line)
                if section_match:
                    current_section = section_match.group(1).lower()
                    data[current_section] = {}
                    continue

                # Detect the end of a section
                if line == '}':
                    current_section = None
                    continue

                # Process key-value pairs if we are inside a section
                if current_section:
                    key_value_match = re.match(r'(\w+)\s*[:=]\s*(.+)', line)
                    if key_value_match:
                        key = key_value_match.group(1).strip()
                        value = key_value_match.group(2).strip()

                        if value.endswith(','):
                            value = value[:-1].strip()

                        data[current_section][key] = _parse_pro_value(value)

            except Exception as e:
                print(f"Warning: Error parsing line {line_num}: '{line}' - {e}")
                continue

    return data

def _parse_pro_value(value):
    """Parse a single value handling IDL-specific syntax."""
    value = value.strip()

    if value == '!VALUES.F_INFINITY':
        return float('inf')

    if value.lower() in ['0b', 'false']:
        return False
    elif value.lower() in ['1b', 'true']:
        return True

    if (value.startswith("'") and value.endswith("'")) or \
       (value.startswith('"') and value.endswith('"')):
        return value[1:-1]

    if value.startswith('[') and value.endswith(']'):
        return _parse_pro_array(value)

    replicate_match = re.match(r'replicate\(([^,]+),\s*(\d+)\)', value, re.IGNORECASE)
    if replicate_match:
        val = _parse_pro_value(replicate_match.group(1))
        num = int(replicate_match.group(2))
        return [val] * num

    if re.match(r'^-?\d+[lL]?$', value):
        return int(value.rstrip('lL'))

    if re.match(r'^-?\d*\.?\d*([eE][+-]?\d+)?[dD]?$', value) and any(c.isdigit() for c in value):
        try:
            return float(value.rstrip('dD'))
        except ValueError:
            pass

    if re.match(r'^-?\d+[eE][+-]?\d+$', value):
        return float(value)

    if '/' in value and re.match(r'^[\d\.\+\-\*/\(\)\s]+$', value):
        try:
            return eval(value)
        except:
            pass

    if value.lower() in ['auto']:
        return value.lower()

    return value

def _parse_pro_array(array_str):
    """Parse a PRO array string such as [val1, val2, val3]."""
    content = array_str[1:-1].strip()
    if not content:
        return []

    elements = []
    parts = content.split(',')

    for part in parts:
        part = part.strip()
        if part:
            if 'replicate(' in part.lower():
                replicate_match = re.match(r'replicate\(([^,]+),\s*(\d+)\)', part, re.IGNORECASE)
                if replicate_match:
                    val = _parse_pro_value(replicate_match.group(1))
                    num = int(replicate_match.group(2))
                    elements.extend([val] * num)
                    continue
            elements.append(_parse_pro_value(part))

    return elements

def convert_file(pro_path, overwrite=False):
    """Convert a file from .pro to .yaml."""
    yaml_path = pro_path.with_suffix('.yaml')

    if yaml_path.exists() and not overwrite:
        print(f"Skipped: {yaml_path.name} already exists. Use --overwrite to overwrite.")
        return

    print(f"Converting: {pro_path.name} -> {yaml_path.name} ...")
    data = parse_pro_file(pro_path)

    with open(yaml_path, 'w') as f:
        # Use default_flow_style=False for a nicely readable multi-line YAML
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)

def main():
    parser = argparse.ArgumentParser(
        description="Convert IDL configuration files (.pro) to YAML."
        )
    parser.add_argument(
        "target", help="A single .pro file or a folder containing .pro files to convert"
        )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing .yaml files"
        )

    args = parser.parse_args()
    target_path = Path(args.target)

    if target_path.is_file() and target_path.suffix == '.pro':
        convert_file(target_path, args.overwrite)
    elif target_path.is_dir():
        pro_files = list(target_path.rglob('*.pro'))
        if not pro_files:
            print(f"No .pro files found in directory {target_path}")
            return

        print(f"Found {len(pro_files)} .pro file(s). Starting conversion...")
        for pro_file in pro_files:
            convert_file(pro_file, args.overwrite)
        print("Conversion completed!")
    else:
        print("Error: the specified target is neither a valid .pro file nor a directory.")

if __name__ == "__main__":
    main()
