#!/usr/bin/env python3
from moke import * #@UnusedWildImport
import re
import os 

@task
def fix_tags(i=stdin, o=stdout):
    for line in i:
        if "tag" in line:
            res = re.findall("tag \"(\S*)\";", line)
            o.write(line.replace("\n", " tags \"%s\";\n" % ",".join(res)))
        else:
            o.write(line)

if __name__ == "__main__":
    task()
