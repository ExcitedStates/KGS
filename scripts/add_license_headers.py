#!/usr/bin/python
# coding: utf-8
"""
Application to add license to header of source code files.
The source files will be analyzed to see if they contain a /*..*/ comment in the head. If they do the contents
are replaced with a license, otherwise a license will be added.
"""

import argparse as ap
import re
import os

parser = ap.ArgumentParser(description=__doc__)

parser.add_argument('--license',
                    dest='licensefile',
                    required=True,
                    metavar='<licensefile>',
                    help='License file to include')

parser.add_argument('--no-backup',
                    dest='backup',
                    action='store_false',
                    help='Disable backing up to .bak files')
parser.set_defaults(backup=True)

parser.add_argument('--sources',
                    dest='sourcefiles',
                    nargs='+',
                    metavar='<sourcefile>',
                    help='Source code files in which to place license'
)

args = parser.parse_args()


# Read license file
with open(args.licensefile, 'r') as f:
    license_contents = f.read()

comment_start = "/*\n"
comment_end = "*/\n"

# Define header pattern and option re.S = DOTALL which means '.*' should match across line-breaks
header_pat = re.compile(r"^\s*/\*.*?\*/", re.S)

# Count how many source files were modified
source_count = 0

# For each source file check for existing license, remove it and add new license
for sourcefname in args.sourcefiles:
    try:
        # Read from source file
        with open(sourcefname, 'r') as f:
            source_contents = f.read()

        # Make backup
        if args.backup:
            with open(sourcefname+".bak", "w") as f:
                f.write(source_contents)

        # Strip old header
        source_contents = header_pat.sub("", source_contents)

        # Add new header
        source_contents = comment_start + license_contents + comment_end + source_contents

        # Write to source file
        with open(sourcefname, "w") as f:
            f.write(source_contents)

        source_count += 1
    except IOError as ioe:
        print("Could not read or write "+ioe.filename)

if args.backup and source_count > 0:
    print("Made backups to .bak files")

print("Done, added licenses to "+str(source_count)+" source files")

