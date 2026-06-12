"""Copyright (C) 2026 UKATC

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>."""

from typing import Final, Dict, Tuple
from starbug2.constants import NIRCAM, SHORT, LONG, NULL, STAR_BUG_MIRI


class FilterStruct: #(struct) containing JWST filter info
    # noinspection SpellCheckingInspection, PyPep8Naming
    def __init__(self, wavelength, aFWHM, pFWHM, instr, length):
        self.wavelength = wavelength
        self.aFWHM = aFWHM
        self.pFWHM = pFWHM
        self.instr = instr
        self.length = length

# as of 08/06/2023
STAR_BUG_FILTERS: Final[Dict[str, FilterStruct]] = {
    "F070W": FilterStruct(0.704, 0.023, 0.742, NIRCAM, SHORT),
    "F090W": FilterStruct(0.901, 0.030, 0.968, NIRCAM, SHORT),
    "F115W": FilterStruct(1.154, 0.037, 1.194, NIRCAM, SHORT),
    "F140M": FilterStruct(1.404, 0.046, 1.484, NIRCAM, SHORT),
    "F150W": FilterStruct(1.501, 0.049, 1.581, NIRCAM, SHORT),
    "F162M": FilterStruct(1.626, 0.053, 1.710, NIRCAM, SHORT),
    "F164N": FilterStruct(1.644, 0.054, 1.742, NIRCAM, SHORT),
    "F150W2": FilterStruct(1.671, 0.045, 1.452, NIRCAM, SHORT),
    "F182M": FilterStruct(1.845, 0.060, 1.935, NIRCAM, SHORT),
    "F187N": FilterStruct(1.874, 0.061, 1.968, NIRCAM, SHORT),
    "F200W": FilterStruct(1.990, 0.064, 2.065, NIRCAM, SHORT),
    "F210M": FilterStruct(2.093, 0.068, 2.194, NIRCAM, SHORT),
    "F212N": FilterStruct(2.120, 0.069, 2.226, NIRCAM, SHORT),
    "F250M": FilterStruct(2.503, 0.082, 1.302, NIRCAM, LONG),
    "F277W": FilterStruct(2.786, 0.088, 1.397, NIRCAM, LONG),
    "F300M": FilterStruct(2.996, 0.097, 1.540, NIRCAM, LONG),
    "F322W2": FilterStruct(3.247, 0.096, 1.524, NIRCAM, LONG),
    "F323N": FilterStruct(3.237, 0.106, 1.683, NIRCAM, LONG),
    "F335M": FilterStruct(3.365, 0.109, 1.730, NIRCAM, LONG),
    "F356W": FilterStruct(3.563, 0.114, 1.810, NIRCAM, LONG),
    "F360M": FilterStruct(3.621, 0.118, 1.873, NIRCAM, LONG),
    "F405N": FilterStruct(4.055, 0.132, 2.095, NIRCAM, LONG),
    "F410M": FilterStruct(4.092, 0.133, 2.111, NIRCAM, LONG),
    "F430M": FilterStruct(4.280, 0.139, 2.206, NIRCAM, LONG),
    "F444W": FilterStruct(4.421, 0.140, 2.222, NIRCAM, LONG),
    "F460M": FilterStruct(4.624, 0.151, 2.397, NIRCAM, LONG),
    "F466N": FilterStruct(4.654, 0.152, 2.413, NIRCAM, LONG),
    "F470N": FilterStruct(4.707, 0.154, 2.444, NIRCAM, LONG),
    "F480M": FilterStruct(4.834, 0.157, 2.492, NIRCAM, LONG),
    "F560W": FilterStruct(5.589, 0.207, 1.882, STAR_BUG_MIRI, NULL),
    "F770W": FilterStruct(7.528, 0.269, 2.445, STAR_BUG_MIRI, NULL),
    "F1000W": FilterStruct(9.883, 0.328, 2.982, STAR_BUG_MIRI, NULL),
    "F1130W": FilterStruct(11.298, 0.375, 3.409, STAR_BUG_MIRI, NULL),
    "F1280W": FilterStruct(12.712, 0.420, 3.818, STAR_BUG_MIRI, NULL),
    "F1500W": FilterStruct(14.932, 0.488, 4.436, STAR_BUG_MIRI, NULL),
    "F1800W": FilterStruct(17.875, 0.591, 5.373, STAR_BUG_MIRI, NULL),
    "F2100W": FilterStruct(20.563, 0.674, 6.127, STAR_BUG_MIRI, NULL),
    "F2550W": FilterStruct(25.147, 0.803, 7.300, STAR_BUG_MIRI, NULL),
}

# ZERO POINT...
ZP: Final[Dict[str, Tuple[int, int]]] = {
   "F070W"	:(3631, 0),
   "F090W"	:(3631, 0),
   "F115W"	:(3631, 0),
   "F140M"	:(3631, 0),
   "F150W"	:(3631, 0),
   "F162M"	:(3631, 0),
   "F164N"	:(3631, 0),
   "F150W2" :(3631, 0),
   "F182M"	:(3631, 0),
   "F187N"	:(3631, 0),
   "F200W"	:(3631, 0),
   "F210M"	:(3631, 0),
   "F212N"	:(3631, 0),
   "F250M"	:(3631, 0),
   "F277W"	:(3631, 0),
   "F300M"	:(3631, 0),
   "F322W2" :(3631, 0),
   "F323N"	:(3631, 0),
   "F335M"	:(3631, 0),
   "F356W"	:(3631, 0),
   "F360M"	:(3631, 0),
   "F405N"	:(3631, 0),
   "F410M"	:(3631, 0),
   "F430M"	:(3631, 0),
   "F444W"	:(3631, 0),
   "F460M"	:(3631, 0),
   "F466N"	:(3631, 0),
   "F470N"	:(3631, 0),
   "F480M"	:(3631, 0),
   "F560W"  :(3631, 0),
   "F770W"  :(3631, 0),
   "F1000W" :(3631, 0),
   "F1130W" :(3631, 0),
   "F1280W" :(3631, 0),
   "F1500W" :(3631, 0),
   "F1800W" :(3631, 0),
   "F2100W" :(3631, 0),
   "F2550W" :(3631, 0),
}
