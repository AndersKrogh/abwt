# abwt

This reposetory contains C programs to construct a Burrows-Wheler index and FM index of a database of biological sequences. Any (small) alphabet can be used. Includes a couple of simple programs using it.

The programs all have a help option (-h).

One program is predictDNA, which calculates conditional probabilities of DNA P(base|context). It is used for this paper:
TO BE ADDED.

### License

Copyright (c) 2019 by Anders Krogh

abwt is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for more details.

You should have received a copy of the GNU General Public License along with the source code.  If not, see <http://www.gnu.org/licenses/>.

### Dependencies

My clib called "aklib" is used. Should be cloned to a directory above abwt (or  hange the Makefile).
```
git clone https://github.com/AndersKrogh/aklib.git
git clone https://github.com/AndersKrogh/abwt.git
```

