NetSquid-FreeSpace (0.0.0)
================================

Description
-----------

This is a user contributed _snippet_ for the [NetSquid quantum network simulator](https://netsquid.org).

It provides a model for the ground-to-ground free-space channel and the downlink satellite channel, exploiting the model described in \[Vasylyev et al., PRL **108**, 220501 (2012)\].

Installation
------------

See the [INSTALL file](INSTALL.md) for instruction of how to install this snippet.

Documentation
-------------

To build and see the docs see the [docs README](docs/README.md).

Usage
-----

This snippet implements two different classes for the "quantum\_loss\_model", i.e., the FreeSpaceLossModel, for a ground-to-ground channel with no pointing error and turbulence all over the length of the channel, and the FixedSatelliteLossModel, describing a satellite-to-ground static channel, with pointing error and turbulence over the last 10 kilometers.

Contributors
------------

Matteo Schiavon - LIP6 (Sorbonne University) - matteo.schiavon@lip6.fr  
Raja Yehia - LIP6 (Sorbonne University) - raja.yehia@lip6.fr  
Tim Coopmans - TU Delft - T.J.Coopmans@tudelft.nl  

License
-------

Copyright (C) 2020- LIP6 (Sorbonne University)

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
