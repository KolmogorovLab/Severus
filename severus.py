#!/usr/bin/env python3

#(c) 2023 by Authors
#This file is a part of Severus program.
#Released under the BSD license (see LICENSE file)

"""
This script sets up environment paths
and invokes Severus without installation.
"""

import os
import sys

#BIN_DIR = "bin"

def main():
    #Setting executable paths
    severus_root = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, severus_root)
    #bin_absolute = os.path.join(severus_root, BIN_DIR)
    #os.environ["PATH"] = bin_absolute + os.pathsep + os.environ["PATH"]

    #Severus entry point
    from severus.main import main
    sys.exit(main())


if __name__ == "__main__":
    main()
