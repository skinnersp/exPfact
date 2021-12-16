"""
Copyright (C) 2019-2020 Emanuele Paci, Simon P. Skinner, Michele Stofella

This program is free software: you can redistribute it and/or modify
it under the terms of version 2 of the GNU General Public License as published
by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import os

list1 = [x.strip() for x in open('none.dat', 'r').readlines()]

for f in os.listdir():
    if f.endswith('.pfact'):
        lines = open(f, 'r').readlines()
        fout = open(f.split('.p')[0]+'.cpfact', 'w')
        for line in lines:
            if not line.strip().split()[0] in list1:
                fout.write(line)
        fout.close()
