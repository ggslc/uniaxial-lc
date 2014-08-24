#  zzz.R: numerical solutions to Ericksen-Leslie equations for 1D
#  liquid crystal films: R internals
#  copyright (C) 2007  Steph Cornford
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
##

#.First.lib <- function(libname,pkgname){
#  library.dynam("EricksenLeslie",pkgname) 
#}

.onLoad <- function(libname,pkgname)
  {
    library.dynam("EricksenLeslie",pkgname,lib.loc=NULL)
  }
