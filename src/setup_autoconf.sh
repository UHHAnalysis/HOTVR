#!/bin/bash

libtoolize -i
aclocal
automake --add-missing --copy
autoconf

